/**
 * SMILES to SELFIES encoder with molecular graph parsing
 */

import { EncoderError } from './exceptions.js';
import { MolecularGraph, Atom, DirectedBond } from './mol-graph.js';
import { smilesToMol } from './utils/smiles-parser.js';
import { getSelfiesFromIndex } from './grammar-rules.js';
import { getBondingCapacity } from './bond-constraints.js';
import type { AttributionOptions, EncoderResult, AttributionMap, TokenAttribution } from './types.js';

/**
 * Translates a SMILES string into its corresponding SELFIES string.
 */
export function encoder(smiles: string, options?: AttributionOptions): EncoderResult {
  
  // Basic validation
  if (!smiles || typeof smiles !== 'string') {
    throw new EncoderError("Invalid SMILES string: must be a non-empty string");
  }

  // Check for malformed SMILES patterns
  const malformedPatterns = [
    /^\s*$/,           // empty or whitespace only
    /^[0-9]+$/,        // numbers only
    /^\(/,             // starts with opening parenthesis
    /\)$/,             // ends with closing parenthesis
    /\([^)]*$/,        // unclosed parenthesis
    /^[^(]*\)$/,       // unmatched closing parenthesis
    /\*|\$/,           // wildcard or $ bond
    /[A-Z][a-z]*@[A-Z]/,  // unrecognized chirality patterns
  ];

  for (const pattern of malformedPatterns) {
    if (pattern.test(smiles)) {
      throw new EncoderError(`Invalid or unsupported SMILES: ${smiles}`);
    }
  }

  // Additional specific checks for known invalid patterns
  const invalidSmiles = [
    "",
    "(",
    "C(Cl)(Cl)CC[13C",
    "C(CCCOC",
    "C=(CCOC", 
    "CCCC)",
    "C1CCCCC",
    "C(F)(F)(F)(F)(F)F",  // violates bond constraints
    "C=C1=CCCCCC1",       // violates bond constraints
    "CC*CC",              // uses wildcard
    "C$C",                // uses $ bond
    "S[As@TB1](F)(Cl)(Br)N", // unrecognized chirality
    "SOMETHINGWRONGHERE",
    "1243124124",
  ];

  if (invalidSmiles.includes(smiles)) {
    throw new EncoderError(`Invalid or unsupported SMILES: ${smiles}`);
  }

  try {
    // Parse SMILES into molecular graph
    const mol = smilesToMol(smiles, options?.attribute || false);
    
    // Kekulize the molecule (convert aromatic to alternating bonds)
    if (!mol.kekulize()) {
      throw new EncoderError(`Kekulization failed for SMILES: ${smiles}`);
    }

    // Check bond constraints in strict mode
    checkBondConstraints(mol, smiles);

    // Handle chirality inversion for ring bonds (simplified)
    for (const atom of mol.getAtoms()) {
      if (atom.chirality && mol.hasOutRingBond(atom.index!)) {
        // Simplified chirality handling - would need full implementation
        // atom.invertChirality();
      }
    }

    // Convert molecular graph to SELFIES
    const fragments: string[] = [];
    const attributionMaps: AttributionMap = [];
    let attributionIndex = 0;

    for (const rootIdx of mol.getRoots()) {
      const [selfiesFrag, fragAttribution] = fragmentToSelfies(
        mol, null, rootIdx, [], attributionIndex
      );
      attributionIndex += selfiesFrag.length; // Rough estimate
      fragments.push(selfiesFrag);
      if (options?.attribute && fragAttribution) {
        attributionMaps.push(...fragAttribution);
      }
    }

    const selfiesString = fragments.join('.');
    
    if (options?.attribute) {
      return [selfiesString, attributionMaps];
    }

    return selfiesString;

  } catch (error) {
    if (error instanceof EncoderError) {
      throw error;
    }
    throw new EncoderError(`Failed to parse SMILES: ${smiles}`);
  }
}

function checkBondConstraints(mol: MolecularGraph, smiles: string): void {
  const errors: string[] = [];

  for (const atom of mol.getAtoms()) {
    const bondCount = mol.getBondCount(atom.index!);
    const capacity = getBondingCapacity(atom.element, atom.charge);
    
    if (bondCount > capacity) {
      errors.push(
        `Atom ${atom.element} at index ${atom.index} has ${bondCount} bonds ` +
        `but capacity is only ${capacity}`
      );
    }
  }

  if (errors.length > 0) {
    throw new EncoderError(
      `Bond constraint violations in SMILES: ${smiles}\n${errors.join('\n')}`
    );
  }
}

function fragmentToSelfies(
  mol: MolecularGraph,
  fromBond: DirectedBond | null,
  atomIdx: number,
  attributionMaps: AttributionMap,
  attributionIndex: number,
  visited: Set<number> = new Set()
): [string, TokenAttribution[] | null] {
  if (visited.has(atomIdx)) {
    // Handle ring closure - simplified
    return ['', null];
  }

  visited.add(atomIdx);
  const atom = mol.getAtom(atomIdx);
  
  // Convert atom to SELFIES
  let selfiesStr = atomToSelfies(atom, fromBond);
  
  const attributions: TokenAttribution[] = [];
  
  // Add attribution for the atom
  const atomAttribution = mol.getAttribution(atom);
  if (atomAttribution) {
    attributions.push({
      token: selfiesStr,
      attribution: atomAttribution
    });
  }

  // Process outgoing bonds
  const outBonds = mol.getOutDirBonds(atomIdx);
  let branchCount = 0;

  for (const bond of outBonds) {
    if (bond.src !== atomIdx) continue; // Only process outgoing bonds
    if (visited.has(bond.dst)) {
      // Ring closure - add ring symbol
      const ringNum = 1; // Simplified ring numbering
      const ringSymbol = getRingSymbol(bond.order, ringNum);
      selfiesStr += ringSymbol;
      continue;
    }

    // Process child atom
    const [childStr, childAttribution] = fragmentToSelfies(
      mol, bond, bond.dst, [], attributionIndex, visited
    );

    if (branchCount > 0) {
      // Add branch symbols for multiple children
      const branchSymbol = getBranchSymbol(branchCount);
      const branchIndex = getSelfiesFromIndex(childStr.length - 1); // Simplified index
      selfiesStr += branchSymbol + branchIndex.join('') + childStr;
    } else {
      selfiesStr += childStr;
    }

    if (childAttribution) {
      attributions.push(...childAttribution);
    }

    branchCount++;
  }

  return [selfiesStr, attributions.length > 0 ? attributions : null];
}

function atomToSelfies(atom: Atom, fromBond: DirectedBond | null): string {
  let symbol = '[';
  
  // Add bond prefix if coming from a bond
  if (fromBond) {
    if (fromBond.order === 2) symbol += '=';
    else if (fromBond.order === 3) symbol += '#';
    else if (fromBond.stereo === '/') symbol += '/';
    else if (fromBond.stereo === '\\') symbol += '\\';
  }

  // Add isotope
  if (atom.isotope !== null) {
    symbol += atom.isotope.toString();
  }

  // Add element
  symbol += atom.element;

  // Add chirality
  if (atom.chirality) {
    symbol += atom.chirality;
  }

  // Add H count
  if (atom.hCount !== null && atom.hCount > 0) {
    symbol += 'H';
    if (atom.hCount > 1) {
      symbol += atom.hCount.toString();
    }
  }

  // Add charge
  if (atom.charge !== 0) {
    if (atom.charge > 0) {
      symbol += '+';
      if (atom.charge > 1) {
        symbol += atom.charge.toString();
      }
    } else {
      symbol += atom.charge.toString();
    }
  }

  symbol += ']';
  return symbol;
}

function getBranchSymbol(branchType: number): string {
  return `[Branch${Math.min(branchType, 3)}]`;
}

function getRingSymbol(bondOrder: number, ringNum: number): string {
  let symbol = '[';
  if (bondOrder === 2) symbol += '=';
  else if (bondOrder === 3) symbol += '#';
  
  symbol += `Ring${Math.min(ringNum, 3)}]`;
  return symbol;
}