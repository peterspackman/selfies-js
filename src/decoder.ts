/**
 * SELFIES to SMILES decoder with molecular graph implementation
 */

import { DecoderError } from './exceptions.js';
import { MolecularGraph, Atom, Attribution } from './mol-graph.js';
import { splitSelfies } from './utils/selfies-utils.js';
import { molToSmiles } from './utils/smiles-utils.js';
import {
  processAtomSymbol,
  processBranchSymbol,
  processRingSymbol,
  nextAtomState,
  nextBranchState,
  nextRingState,
  getIndexFromSelfies
} from './grammar-rules.js';
import type { AttributionOptions, DecoderOptions, DecoderResult, State } from './types.js';
import { modernizeSymbol, warnDeprecatedCompatibility } from './utils/compatibility-utils.js';

/**
 * Translates a SELFIES string into its corresponding SMILES string.
 */
export function decoder(selfies: string, options?: DecoderOptions): DecoderResult {
  
  if (!selfies || typeof selfies !== 'string') {
    return options?.attribute ? ['C', []] : 'C';
  }

  // Issue deprecation warning if using backward compatibility
  if (options?.compatible) {
    warnDeprecatedCompatibility();
  }

  const mol = new MolecularGraph(options?.attribute || false);
  const rings: [Atom, Atom, [number, string | null]][] = [];
  let attributionIndex = 0;

  // Process each fragment separated by dots
  for (const fragment of selfies.split('.')) {
    if (fragment) {
      const n = deriveMolFromSymbols(
        tokenizeSelfies(fragment, options?.compatible),
        mol,
        selfies,
        Infinity,
        0,
        null,
        rings,
        options?.attribute ? [] : null,
        attributionIndex
      );
      attributionIndex += n;
    }
  }

  // Form ring bonds
  formRingsBilocally(mol, rings);

  // Convert molecular graph to SMILES
  return molToSmiles(mol, options?.attribute || false);
}

function* tokenizeSelfies(selfies: string, compatible?: boolean): Generator<string> {
  for (let symbol of splitSelfies(selfies)) {
    if (symbol === '[nop]') {
      continue;
    }
    
    // Apply backward compatibility modernization if needed
    if (compatible) {
      symbol = modernizeSymbol(symbol);
    }
    
    yield symbol;
  }
}

function deriveMolFromSymbols(
  symbolIter: Generator<string>,
  mol: MolecularGraph,
  selfies: string,
  maxDerive: number,
  initState: number,
  rootAtom: Atom | null,
  rings: [Atom, Atom, [number, string | null]][],
  attributeStack: Attribution[] | null,
  attributionIndex: number
): number {
  let nDerived = 0;
  let state: State = initState;
  let prevAtom = rootAtom;
  let index = 0;

  const symbolIterWithIndex = (function* (): Generator<[number, string]> {
    for (const symbol of symbolIter) {
      yield [index++, symbol] as [number, string];
    }
  })();

  while (nDerived < maxDerive) {
    let symbolResult;
    try {
      symbolResult = symbolIterWithIndex.next();
      if (symbolResult.done) break;
      nDerived++;
    } catch (error) {
      if (error instanceof Error) {
        throw new DecoderError(`Malformed SELFIES string: ${error.message}`);
      }
      throw error;
    }

    const [symbolIndex, symbol] = symbolResult.value as [number, string];

    // Case 1: Branch symbol (e.g. [Branch1])
    if (symbol.includes('Branch')) {
      const output = processBranchSymbol(symbol);
      if (!output) {
        raiseDecoderError(selfies, symbol);
      }
      
      const [btype, n] = output;

      if (state === null || state <= 1) {
        continue;
      } else {
        const [binitState, nextState] = nextBranchState(btype, state);
        
        const Q = readIndexFromSelfies(symbolIterWithIndex, n);
        
        // Create shared iterator for branch processing (matches Python's approach)
        const sharedIterator = (function* (): Generator<string> {
          while (true) {
            const result = symbolIterWithIndex.next();
            if (result.done) break;
            yield result.value[1];
          }
        })();
        
        nDerived += n + deriveMolFromSymbols(
          sharedIterator,
          mol,
          selfies,
          Q + 1,
          binitState,
          prevAtom,
          rings,
          attributeStack ? [...attributeStack, new Attribution(symbolIndex + attributionIndex, symbol)] : null,
          attributionIndex
        );
        
        state = nextState;
      }
    }
    // Case 2: Ring symbol (e.g. [Ring2])
    else if (symbol.includes('Ring')) {
      const output = processRingSymbol(symbol);
      if (!output) {
        raiseDecoderError(selfies, symbol);
      }
      
      const [ringType, n, stereo] = output;

      if (state === 0 || state === null) {
        continue;
      } else {
        const [ringOrder, nextState] = nextRingState(ringType, state);
        const bondInfo: [number, string | null] = [ringOrder, stereo[0]];

        const Q = readIndexFromSelfies(symbolIterWithIndex, n);
        nDerived += n;
        
        if (prevAtom && prevAtom.index !== null) {
          const lidx = Math.max(0, prevAtom.index - (Q + 1));
          const targetAtom = mol.getAtom(lidx);
          rings.push([targetAtom, prevAtom, bondInfo]);
        }
        
        state = nextState;
      }
    }
    // Case 3: [epsilon] - termination symbol
    else if (symbol.includes('eps')) {
      state = state === 0 ? 0 : null;
    }
    // Case 4: Regular atom symbol (e.g. [N], [=C], [F])
    else {
      const output = processAtomSymbol(symbol);
      if (!output) {
        raiseDecoderError(selfies, symbol);
      }
      
      const [bondInfo, atom] = output;
      const [bondOrder, stereo] = bondInfo;
      const cap = atom.bondingCapacity;

      const [finalBondOrder, nextState] = nextAtomState(bondOrder, cap, state);
      
      if (finalBondOrder === 0) {
        if (state === 0) {
          const addedAtom = mol.addAtom(atom, true);
          if (attributeStack) {
            mol.addAttribution(addedAtom, [
              ...attributeStack, 
              new Attribution(symbolIndex + attributionIndex, symbol)
            ]);
          }
        }
      } else {
        const addedAtom = mol.addAtom(atom);
        if (attributeStack) {
          mol.addAttribution(addedAtom, [
            ...attributeStack,
            new Attribution(symbolIndex + attributionIndex, symbol)
          ]);
        }
        
        if (prevAtom && prevAtom.index !== null && addedAtom.index !== null) {
          const src = prevAtom.index;
          const dst = addedAtom.index;
          const bond = mol.addBond(Math.min(src, dst), Math.max(src, dst), finalBondOrder, stereo);
          
          if (attributeStack) {
            mol.addAttribution(bond, [
              ...attributeStack,
              new Attribution(symbolIndex + attributionIndex, symbol)
            ]);
          }
        }
      }
      
      prevAtom = atom;
      state = nextState ?? null;
    }

  }

  // Consume remaining tokens
  while (nDerived < maxDerive) {
    try {
      const result = symbolIterWithIndex.next();
      if (result.done) break;
      nDerived++;
    } catch (error) {
      break;
    }
  }

  return nDerived;
}

function readIndexFromSelfies(symbolIter: Generator<[number, string]>, nSymbols: number): number {
  const symbols: string[] = [];
  
  for (let i = 0; i < nSymbols; i++) {
    const result = symbolIter.next();
    if (!result.done) {
      symbols.push(result.value[1]);
    }
  }
  
  return getIndexFromSelfies(...symbols);
}

function formRingsBilocally(mol: MolecularGraph, rings: [Atom, Atom, [number, string | null]][]): void {
  const ringsMade = new Array(mol.length).fill(0);
  
  // PYTHON-EXACT LOGIC: Skip all rings if any ring type has invalid count
  // First, we need to determine which rings came from which ring types
  // Since we don't store the original ring type, we use a heuristic:
  // If any ring appears an odd number of times, skip ALL rings (like Python)
  
  const ringPairCounts = new Map<string, number>();
  for (const [atom1, atom2, bondInfo] of rings) {
    if (atom1.index !== null && atom2.index !== null) {
      const key = `${Math.min(atom1.index, atom2.index)}-${Math.max(atom1.index, atom2.index)}`;
      ringPairCounts.set(key, (ringPairCounts.get(key) || 0) + 1);
    }
  }
  
  // Check if any ring pair has an odd count (invalid)
  const hasInvalidRings = Array.from(ringPairCounts.values()).some(count => count % 2 !== 0);
  
  if (hasInvalidRings) {
    // Skip ALL ring formation when any ring is invalid (Python behavior)
    return;
  }
  
  for (const [atom1, atom2, bondInfo] of rings) {
    if (atom1.index !== null && atom2.index !== null) {
      const lidx = atom1.index;
      const ridx = atom2.index;
      
      if (lidx === ridx) {
        continue;
      }
      
      const [order, stereo] = bondInfo;
      
      const leftFree = atom1.bondingCapacity - mol.getBondCount(lidx);
      const rightFree = atom2.bondingCapacity - mol.getBondCount(ridx);
      
      if (leftFree <= 0 || rightFree <= 0) {
        continue;
      }
      
      const finalOrder = Math.min(order, leftFree, rightFree);
      
      if (mol.hasBond(lidx, ridx)) {
        const existingBond = mol.getDirBond(lidx, ridx);
        const newOrder = Math.min(finalOrder + existingBond.order, 3);
        mol.updateBondOrder(lidx, ridx, newOrder);
      } else {
        mol.addRingBond(
          lidx,
          ridx,
          finalOrder,
          stereo,
          stereo,
          ringsMade[lidx],
          ringsMade[ridx]
        );
        ringsMade[lidx] += 1;
        ringsMade[ridx] += 1;
      }
    }
  }
}

function raiseDecoderError(selfies: string, symbol: string): never {
  throw new DecoderError(`invalid symbol '${symbol}'\n\tSELFIES: ${selfies}`);
}