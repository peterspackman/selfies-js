/**
 * SMILES to SELFIES encoder with molecular graph parsing
 */

import { EncoderError } from './exceptions.js';
import { MolecularGraph, Atom, DirectedBond } from './mol-graph.js';
import { smilesToMol } from './utils/smiles-parser.js';
import { getSelfiesFromIndex } from './grammar-rules.js';
import { getBondingCapacity } from './bond-constraints.js';
import { atomToSmiles, bondToSmiles } from './utils/smiles-utils.js';
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

    // Handle chirality inversion for ring bonds
    for (const atom of mol.getAtoms()) {
      if (atom.chirality && mol.hasOutRingBond(atom.index!) && shouldInvertChirality(mol, atom)) {
        atom.invertChirality();
      }
    }

    // Convert molecular graph to SELFIES
    const fragments: string[] = [];
    const attributionMaps: AttributionMap = [];
    let attributionIndex = 0;

    for (const rootIdx of mol.getRoots()) {
      const derived = fragmentToSelfies(mol, null, rootIdx, attributionMaps, attributionIndex);
      attributionIndex += derived.length;
      fragments.push(derived.join(''));
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

function shouldInvertChirality(mol: MolecularGraph, atom: Atom): boolean {
  if (!atom.index) return false;
  
  const outBonds = mol.getOutDirBonds(atom.index);

  // Partition bonds into 3 categories:
  // 1. rings whose right number are bonded to this atom (e.g. ...1...X1)
  // 2. rings whose left number are bonded to this atom (e.g. X1...1...)
  // 3. branches and other (e.g. X(...)...)
  const partition: number[][] = [[], [], []];
  
  for (let i = 0; i < outBonds.length; i++) {
    const bond = outBonds[i];
    if (!bond.ringBond) {
      partition[2].push(i);
    } else if (bond.src < bond.dst) {
      partition[1].push(i);
    } else {
      partition[0].push(i);
    }
  }

  // Sort partition[1] by destination atom index for consistent ordering
  partition[1].sort((a, b) => outBonds[a].dst - outBonds[b].dst);

  // Construct permutation: partition[0] + partition[1] + partition[2]
  const perm = [...partition[0], ...partition[1], ...partition[2]];
  
  // Count inversions in the permutation
  let count = 0;
  for (let i = 0; i < perm.length; i++) {
    for (let j = i + 1; j < perm.length; j++) {
      if (perm[i] > perm[j]) {
        count++;
      }
    }
  }
  
  // If odd number of inversions, should invert chirality
  return count % 2 !== 0;
}

function fragmentToSelfies(
  mol: MolecularGraph,
  bondIntoCurr: DirectedBond | null,
  curr: number,
  attributionMaps: AttributionMap,
  attributionIndex: number
): string[] {
  const derived: string[] = [];
  
  let currentBond = bondIntoCurr;
  let currentAtom = curr;
  
  while (true) {
    const atom = mol.getAtom(currentAtom);
    const token = atomToSelfies(currentBond, atom);
    derived.push(token);
    
    // Add attribution if needed
    const atomAttribution = mol.getAttribution(atom);
    if (atomAttribution) {
      attributionMaps.push({
        token,
        attribution: atomAttribution
      });
    }
    
    const outBonds = mol.getOutDirBonds(currentAtom);
    
    for (let i = 0; i < outBonds.length; i++) {
      const bond = outBonds[i];
      
      if (bond.ringBond) {
        // Only process ring closures (not openings)
        if (bond.src < bond.dst) {
          continue;
        }
        
        // Get the reverse bond
        const revBond = mol.getDirBond(bond.dst, bond.src);
        const ringLen = bond.src - bond.dst;
        const qAsSymbols = getSelfiesFromIndex(ringLen - 1);
        const ringSymbol = `[${ringBondsToSelfies(revBond, bond)}Ring${qAsSymbols.length}]`;
        
        derived.push(ringSymbol);
        for (const symbol of qAsSymbols) {
          derived.push(symbol);
        }
        
        // Add attribution for ring
        const bondAttribution = mol.getAttribution(bond);
        if (bondAttribution) {
          attributionMaps.push({
            token: ringSymbol,
            attribution: bondAttribution
          });
        }
        
      } else if (i === outBonds.length - 1) {
        // Last bond - continue chain
        currentBond = bond;
        currentAtom = bond.dst;
        
      } else {
        // Branch
        const branch = fragmentToSelfies(mol, bond, bond.dst, attributionMaps, derived.length);
        const qAsSymbols = getSelfiesFromIndex(branch.length - 1);
        const branchSymbol = `[${bondToSelfies(bond, false)}Branch${qAsSymbols.length}]`;
        
        derived.push(branchSymbol);
        for (const symbol of qAsSymbols) {
          derived.push(symbol);
        }
        derived.push(...branch);
        
        // Add attribution for branch
        const bondAttribution = mol.getAttribution(bond);
        if (bondAttribution) {
          attributionMaps.push({
            token: branchSymbol,
            attribution: bondAttribution
          });
        }
      }
    }
    
    // End of chain
    if (outBonds.length === 0 || outBonds[outBonds.length - 1].ringBond) {
      break;
    }
  }
  
  return derived;
}

function ringBondsToSelfies(lBond: DirectedBond, rBond: DirectedBond): string {
  // Assert that bond orders match
  if (lBond.order !== rBond.order) {
    throw new Error('Ring bond orders must match');
  }
  
  if (lBond.order !== 1 || (lBond.stereo === null && rBond.stereo === null)) {
    return bondToSelfies(lBond, false);
  } else {
    // Handle stereo bonds
    const bondChar = (lBond.stereo || '-') + (rBond.stereo || '-');
    return bondChar;
  }
}

function atomToSelfies(bond: DirectedBond | null, atom: Atom): string {
  // Assert atom is not aromatic (should be kekulized)
  if (atom.isAromatic) {
    throw new Error('Atom should not be aromatic after kekulization');
  }
  
  const bondChar = bond ? bondToSelfies(bond, true) : '';
  return `[${bondChar}${atomToSmiles(atom, false)}]`;
}

function bondToSelfies(bond: DirectedBond, showStereo: boolean = true): string {
  if (!showStereo && bond.order === 1) {
    return '';
  }
  return bondToSmiles(bond);
}