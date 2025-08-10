/**
 * SMILES utility functions for molecular graph conversion
 */

import { MolecularGraph, Atom, DirectedBond } from '../mol-graph.js';
import type { AttributionMap, TokenAttribution, Attribution } from '../types.js';
import { ORGANIC_SUBSET } from '../constants.js';

class RingContext {
  private nextRingNumber = 1;
  private ringBondToNumber = new Map<string, number>();
  
  getRingNumber(atomIdx1: number, atomIdx2: number): number {
    const key1 = `${atomIdx1}-${atomIdx2}`;
    const key2 = `${atomIdx2}-${atomIdx1}`;
    
    let ringNumber = this.ringBondToNumber.get(key1) ?? this.ringBondToNumber.get(key2);
    
    if (ringNumber === undefined) {
      ringNumber = this.nextRingNumber++;
      this.ringBondToNumber.set(key1, ringNumber);
      this.ringBondToNumber.set(key2, ringNumber);
    }
    
    return ringNumber;
  }
}

export function molToSmiles(
  mol: MolecularGraph, 
  attribute: boolean = false
): string | [string, AttributionMap] {
  if (mol.length === 0) {
    return attribute ? ['', []] : '';
  }

  const fragments: string[] = [];
  const attributionMaps: TokenAttribution[] = [];
  
  for (const rootIdx of mol.getRoots()) {
    const visited = new Set<number>();
    const ringContext = new RingContext();
    const [smilesFrag, fragAttribution] = moleculeAtomToSmiles(
      mol, rootIdx, visited, attribute, undefined, ringContext
    );
    fragments.push(smilesFrag);
    if (attribute && fragAttribution) {
      attributionMaps.push(...fragAttribution);
    }
  }

  const smilesString = fragments.join('.');
  
  if (attribute) {
    return [smilesString, attributionMaps];
  }
  
  return smilesString;
}

function moleculeAtomToSmiles(
  mol: MolecularGraph,
  atomIdx: number,
  visited: Set<number>,
  attribute: boolean,
  fromBond?: DirectedBond,
  ringContext?: RingContext
): [string, TokenAttribution[] | null] {
  if (visited.has(atomIdx)) {
    // This should handle ring closures, but for now just return empty
    return ['', null];
  }

  visited.add(atomIdx);
  const atom = mol.getAtom(atomIdx);
  
  // Generate atom string
  const atomStr = formatAtomForSmiles(atom);
  
  const attributions: TokenAttribution[] = [];
  if (attribute) {
    const atomAttribution = mol.getAttribution(atom);
    if (atomAttribution) {
      attributions.push({
        token: atomStr,
        attribution: atomAttribution
      });
    }
  }

  // Process outgoing bonds
  let smilesStr = atomStr;
  const outBonds = mol.getOutDirBonds(atomIdx);
  
  // Separate regular bonds from ring bonds
  const regularBonds: DirectedBond[] = [];
  const ringBonds: DirectedBond[] = [];
  
  for (const bond of outBonds) {
    if (bond.dst === atomIdx) continue; // Skip reverse bonds
    
    if (bond.ringBond) {
      ringBonds.push(bond);
    } else if (!visited.has(bond.dst)) {
      regularBonds.push(bond);
    }
  }
  
  // Add ring closure numbers for ring bonds
  for (const bond of ringBonds) {
    if (ringContext) {
      const ringNumber = ringContext.getRingNumber(atomIdx, bond.dst);
      // Ring closure numbers should not have bond order prefixes
      // Bond order is expressed through the bond path, not the closure
      smilesStr += ringNumber.toString();
    }
  }
  
  // Process regular bonds to unvisited atoms
  for (let bondIndex = 0; bondIndex < regularBonds.length; bondIndex++) {
    const bond = regularBonds[bondIndex];
    const [bondStr, bondAttribution] = moleculeBondToSmiles(bond, attribute ? mol.getAttribution(bond) : null);
    const [childStr, childAttribution] = moleculeAtomToSmiles(mol, bond.dst, visited, attribute, bond, ringContext);
    
    let bondFragment = bondStr + childStr;
    
    // Add parentheses for branching - PYTHON'S "LAST BOND WINS" STRATEGY
    // All bonds except the last one are treated as branches
    if (bondIndex < regularBonds.length - 1) {
      bondFragment = '(' + bondFragment + ')';
    }
    
    smilesStr += bondFragment;
    
    if (attribute) {
      if (bondAttribution) attributions.push(bondAttribution);
      if (childAttribution) attributions.push(...childAttribution);
    }
  }

  return [smilesStr, attribute ? attributions : null];
}

function formatAtomForSmiles(atom: Atom): string {
  // Handle organic subset shorthand
  if (ORGANIC_SUBSET.has(atom.element) && 
      atom.charge === 0 && 
      !atom.chirality && 
      !atom.isotope &&
      (atom.hCount === null || atom.hCount === 0)) {
    return atom.element;
  }

  // Build bracketed atom notation
  let result = '[';
  
  if (atom.isotope !== null) {
    result += atom.isotope.toString();
  }
  
  result += atom.element;
  
  if (atom.chirality) {
    result += atom.chirality;
  }
  
  if (atom.hCount !== null && atom.hCount > 0) {
    result += 'H';
    if (atom.hCount > 1) {
      result += atom.hCount.toString();
    }
  }
  
  if (atom.charge !== 0) {
    if (atom.charge > 0) {
      result += '+';
      result += atom.charge.toString();
    } else {
      result += atom.charge.toString();
    }
  }
  
  result += ']';
  return result;
}

function moleculeBondToSmiles(
  bond: DirectedBond,
  attribution: Attribution[] | null
): [string, TokenAttribution | null] {
  let bondStr = '';
  
  switch (bond.order) {
    case 1:
      if (bond.stereo === '/') bondStr = '/';
      else if (bond.stereo === '\\') bondStr = '\\';
      // Single bonds are usually implicit
      break;
    case 2:
      bondStr = '=';
      break;
    case 3:
      bondStr = '#';
      break;
    case 1.5:
      // Aromatic bonds are typically implicit, but shouldn't occur in kekulized form
      break;
  }

  const tokenAttribution: TokenAttribution | null = attribution ? {
    token: bondStr,
    attribution: attribution
  } : null;

  return [bondStr, tokenAttribution];
}

// Simple SMILES validation
/**
 * Convert atom to SMILES representation (for encoder)
 */
export function atomToSmiles(atom: Atom, brackets: boolean = true): string {
  let symbol = '';
  
  // Add isotope
  if (atom.isotope !== null) {
    symbol += atom.isotope.toString();
  }
  
  // Add element (uppercase for non-aromatic)
  symbol += atom.element.toUpperCase();
  
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
      symbol += atom.charge.toString();
    } else {
      symbol += atom.charge.toString();
    }
  }
  
  return brackets ? `[${symbol}]` : symbol;
}

/**
 * Convert bond to SMILES representation (for encoder)
 */
export function bondToSmiles(bond: DirectedBond): string {
  if (bond.order === 1) {
    return bond.stereo || '';
  } else if (bond.order === 2) {
    return '=';
  } else if (bond.order === 3) {
    return '#';
  } else {
    return '';
  }
}

export function isValidSmiles(smiles: string): boolean {
  if (!smiles || smiles.trim() === '') return false;
  
  // Check for basic structural issues
  const openParens = (smiles.match(/\(/g) || []).length;
  const closeParens = (smiles.match(/\)/g) || []).length;
  if (openParens !== closeParens) return false;
  
  const openBrackets = (smiles.match(/\[/g) || []).length;
  const closeBrackets = (smiles.match(/\]/g) || []).length;
  if (openBrackets !== closeBrackets) return false;
  
  // Check for invalid characters
  if (/[*$]/.test(smiles)) return false;
  
  // Check for malformed patterns
  if (/^\(/.test(smiles) || /\)$/.test(smiles)) return false;
  if (/\([^)]*$/.test(smiles) || /^[^(]*\)$/.test(smiles)) return false;
  
  return true;
}