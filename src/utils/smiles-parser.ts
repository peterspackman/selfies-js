/**
 * Simple SMILES parser for creating molecular graphs
 */

import { MolecularGraph, Atom } from '../mol-graph.js';
import { EncoderError } from '../exceptions.js';
import { ORGANIC_SUBSET, AROMATIC_SUBSET } from '../constants.js';

export function smilesToMol(smiles: string, attributable: boolean = false): MolecularGraph {
  if (!smiles || typeof smiles !== 'string') {
    throw new EncoderError('Invalid SMILES: must be a non-empty string');
  }

  const mol = new MolecularGraph(attributable);
  const fragments = smiles.split('.');
  
  for (const fragment of fragments) {
    if (fragment.trim()) {
      parseFragment(fragment.trim(), mol);
    }
  }

  return mol;
}

function parseFragment(smiles: string, mol: MolecularGraph): void {
  let i = 0;
  let currentAtom: Atom | null = null;
  let isFirstAtom = true;
  const ringClosures = new Map<string, { atom: Atom; bondOrder: number; stereo: string | null }>();
  const parenStack: Atom[] = [];

  while (i < smiles.length) {
    // Skip whitespace
    if (/\s/.test(smiles[i])) {
      i++;
      continue;
    }

    // Handle parentheses for branching
    if (smiles[i] === '(') {
      if (currentAtom) {
        parenStack.push(currentAtom);
      }
      i++;
      continue;
    }

    if (smiles[i] === ')') {
      if (parenStack.length > 0) {
        currentAtom = parenStack.pop()!;
      }
      i++;
      continue;
    }

    // Parse bond
    let bondOrder = 1;
    let stereo: string | null = null;
    
    if (i < smiles.length && /[=\-#/\\]/.test(smiles[i])) {
      switch (smiles[i]) {
        case '=':
          bondOrder = 2;
          break;
        case '#':
          bondOrder = 3;
          break;
        case '/':
          bondOrder = 1;
          stereo = '/';
          break;
        case '\\':
          bondOrder = 1;
          stereo = '\\';
          break;
        case '-':
          bondOrder = 1;
          break;
      }
      i++;
    }

    // Parse atom or ring closure
    if (i >= smiles.length) break;

    // Ring closure (number)
    if (/\d/.test(smiles[i])) {
      const ringNum = smiles[i];
      
      if (currentAtom) {
        if (ringClosures.has(ringNum)) {
          // Close the ring
          const ring = ringClosures.get(ringNum)!;
          const atom1 = ring.atom;
          const atom2 = currentAtom;
          
          // Use the bond order from either the opening or closing
          let finalBondOrder = Math.max(bondOrder, ring.bondOrder);
          
          // Handle implicit aromatic bonds for rings (e.g., c1ccccc1)
          if (finalBondOrder === 1 && atom1.isAromatic && atom2.isAromatic) {
            finalBondOrder = 1.5; // Aromatic bond
          }
          
          if (atom1.index !== null && atom2.index !== null) {
            mol.addRingBond(
              Math.min(atom1.index, atom2.index),
              Math.max(atom1.index, atom2.index),
              finalBondOrder,
              ring.stereo,
              stereo
            );
          }
          
          ringClosures.delete(ringNum);
        } else {
          // Open a ring
          ringClosures.set(ringNum, { 
            atom: currentAtom, 
            bondOrder, 
            stereo 
          });
        }
      }
      i++;
      continue;
    }

    // Parse atom
    const [atom, newI] = parseAtom(smiles, i);
    i = newI;

    const addedAtom = mol.addAtom(atom, isFirstAtom);
    
    // Connect to previous atom if not the first
    if (!isFirstAtom && currentAtom) {
      const src = currentAtom.index!;
      const dst = addedAtom.index!;
      
      // Check if both atoms are aromatic and no explicit bond order given
      let actualBondOrder = bondOrder;
      if (bondOrder === 1 && currentAtom.isAromatic && addedAtom.isAromatic) {
        actualBondOrder = 1.5; // Aromatic bond
      }
      
      mol.addBond(Math.min(src, dst), Math.max(src, dst), actualBondOrder, stereo);
    }

    currentAtom = addedAtom;
    isFirstAtom = false;
  }

  // Check for unclosed rings
  if (ringClosures.size > 0) {
    throw new EncoderError(`Unclosed rings in SMILES: ${[...ringClosures.keys()].join(', ')}`);
  }
}

function parseAtom(smiles: string, start: number): [Atom, number] {
  let i = start;

  // Bracketed atom: [...]
  if (smiles[i] === '[') {
    const endBracket = smiles.indexOf(']', i);
    if (endBracket === -1) {
      throw new EncoderError(`Unclosed bracket in SMILES at position ${i}`);
    }

    const bracketContent = smiles.slice(i + 1, endBracket);
    const atom = parseBracketedAtom(bracketContent);
    return [atom, endBracket + 1];
  }

  // Organic subset atom
  let element = '';
  if (i < smiles.length) {
    // Handle two-character elements
    if (i + 1 < smiles.length) {
      const twoChar = smiles.slice(i, i + 2);
      if (['Cl', 'Br'].includes(twoChar)) {
        element = twoChar;
        i += 2;
      } else {
        element = smiles[i];
        i++;
      }
    } else {
      element = smiles[i];
      i++;
    }
  }

  if (!element || (!ORGANIC_SUBSET.has(element.toUpperCase()) && !['Cl', 'Br'].includes(element) && !AROMATIC_SUBSET.has(element))) {
    throw new EncoderError(`Invalid element: ${element}`);
  }

  const atom = new Atom(element, element.toLowerCase() === element);
  return [atom, i];
}

function parseBracketedAtom(content: string): Atom {
  // Parse: [isotope?][element][chirality?][H[count]?][charge?]
  let i = 0;
  let isotope: number | null = null;
  let element = '';
  let chirality: string | null = null;
  let hCount: number | null = null;
  let charge = 0;

  // Parse isotope
  let isotopeStr = '';
  while (i < content.length && /\d/.test(content[i])) {
    isotopeStr += content[i];
    i++;
  }
  if (isotopeStr) {
    isotope = parseInt(isotopeStr, 10);
  }

  // Parse element
  if (i < content.length && /[A-Z]/.test(content[i])) {
    element = content[i];
    i++;
    if (i < content.length && /[a-z]/.test(content[i])) {
      element += content[i];
      i++;
    }
  } else {
    throw new EncoderError(`Invalid atom format: ${content}`);
  }

  // Parse chirality
  if (i < content.length && content[i] === '@') {
    chirality = '@';
    i++;
    if (i < content.length && content[i] === '@') {
      chirality = '@@';
      i++;
    }
  }

  // Parse H count
  if (i < content.length && content[i] === 'H') {
    i++;
    let hCountStr = '';
    while (i < content.length && /\d/.test(content[i])) {
      hCountStr += content[i];
      i++;
    }
    hCount = hCountStr ? parseInt(hCountStr, 10) : 1;
  }

  // Parse charge
  while (i < content.length && /[+-]/.test(content[i])) {
    const sign = content[i] === '+' ? 1 : -1;
    i++;
    
    let chargeStr = '';
    while (i < content.length && /\d/.test(content[i])) {
      chargeStr += content[i];
      i++;
    }
    
    const chargeVal = chargeStr ? parseInt(chargeStr, 10) : 1;
    charge += sign * chargeVal;
  }

  return new Atom(element, false, isotope, chirality, hCount, charge);
}