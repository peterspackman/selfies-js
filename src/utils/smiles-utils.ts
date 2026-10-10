/**
 * SMILES utility functions for molecular graph conversion
 */

import { MolecularGraph, Atom, DirectedBond } from '../mol-graph.js';
import type { AttributionMap } from '../types.js';
import { ORGANIC_SUBSET } from '../constants.js';
import { SMILES_STEREO_BONDS } from './smiles-parser.js';

/**
 * Converts an atom into its SMILES representation.
 *
 * @param atom - the input atom
 * @param brackets - true if brackets should be added around the returned
 *   symbol (e.g. in the case of [C] or [C@@H])
 * @returns a SMILES symbol representing the input atom
 */
export function atomToSmiles(atom: Atom, brackets: boolean = true): string {
  if (atom.isAromatic) {
    throw new Error('Cannot write an aromatic atom; the molecule must be kekulized first');
  }

  if (atom.isotope === null && atom.chirality === null && atom.hCount === null && atom.charge === 0) {
    return atom.element;
  }

  let symbol = '';
  if (atom.isotope !== null) {
    symbol += atom.isotope.toString();
  }
  symbol += atom.element;
  if (atom.chirality !== null) {
    symbol += atom.chirality;
  }
  if (atom.hCount !== null && atom.hCount !== 0) {
    symbol += `H${atom.hCount}`;
  } else if (atom.hCount === 0 && atom.isotope === null && atom.chirality === null &&
             atom.charge === 0 && ORGANIC_SUBSET.has(atom.element)) {
    // Explicitly zero Hs, so that e.g. [CH0] is not read back as an organic-subset C
    symbol += 'H0';
  }
  if (atom.charge !== 0) {
    symbol += (atom.charge > 0 ? '+' : '') + atom.charge.toString();
  }

  return brackets ? `[${symbol}]` : symbol;
}

/**
 * Converts a bond into its SMILES representation.
 *
 * @param bond - the input bond
 * @returns a SMILES symbol representing the input bond
 */
export function bondToSmiles(bond: DirectedBond): string {
  if (bond.order === 1) {
    return bond.stereo !== null && SMILES_STEREO_BONDS.has(bond.stereo) ? bond.stereo : '';
  } else if (bond.order === 2) {
    return '=';
  } else if (bond.order === 3) {
    return '#';
  }
  throw new Error(`Cannot write bond of order ${bond.order}`);
}

/**
 * Converts a molecular graph into its SMILES representation, maintaining
 * the traversal order indicated by the input graph.
 *
 * @param mol - the input molecule
 * @param attribute - if an attribution should be returned
 */
export function molToSmiles(
  mol: MolecularGraph,
  attribute: boolean = false
): string | [string, AttributionMap] {
  if (!mol.isKekulized()) {
    throw new Error('Molecule must be kekulized before writing SMILES');
  }

  const fragments: string[] = [];
  let attributionMaps: AttributionMap | null = attribute ? [] : null;
  let attributionIndex = 0;
  const ringLog = new RingLog();

  for (const root of mol.getRoots()) {
    const derived: string[] = [];
    deriveSmilesFromFragment(derived, mol, root, ringLog, attributionMaps, attributionIndex);
    const fragment = derived.join('');
    attributionIndex += fragment.length + 1; // +1 for the '.' separator
    fragments.push(fragment);
  }

  const smiles = fragments.join('.');
  if (attributionMaps === null) {
    return smiles;
  }
  // Trim attribution map of empty tokens
  attributionMaps = attributionMaps.filter(a => a.token);
  return [smiles, attributionMaps];
}

/**
 * Assigns ring-closure numbers. Like the Python implementation, each ring
 * bond gets the next unused number; but once those run out at 99, numbers of
 * already-closed rings are reused, so the output stays valid SMILES (Python
 * writes e.g. %100, which SMILES parsers reject).
 */
class RingLog {
  private numbers = new Map<string, number>();
  private open = new Set<number>();
  private nextFresh = 1;

  numberFor(key: string): number {
    const existing = this.numbers.get(key);
    if (existing !== undefined) {
      // Second endpoint: the ring is now closed and its number is free again
      this.numbers.delete(key);
      this.open.delete(existing);
      return existing;
    }

    let rnum: number;
    if (this.nextFresh <= 99) {
      rnum = this.nextFresh++;
    } else {
      rnum = 1;
      while (this.open.has(rnum)) rnum++;
    }
    this.numbers.set(key, rnum);
    this.open.add(rnum);
    return rnum;
  }
}

function ringNumberToSmiles(rnum: number): string {
  if (rnum < 10) return rnum.toString();
  if (rnum < 100) return `%${rnum}`;
  return `%(${rnum})`;
}

function deriveSmilesFromFragment(
  derived: string[],
  mol: MolecularGraph,
  root: number,
  ringLog: RingLog,
  attributionMaps: AttributionMap | null,
  attributionIndex: number
): void {
  // Iterative depth-first traversal; avoids stack overflows on long chains
  let derivedLen = 0;
  // Appends a token and returns its starting character index in the full SMILES
  const push = (token: string): number => {
    const start = derivedLen + attributionIndex;
    derived.push(token);
    derivedLen += token.length;
    return start;
  };
  const record = (token: string, attribution: ReturnType<MolecularGraph['getAttribution']>, index: number) => {
    if (attributionMaps !== null) {
      attributionMaps.push({ index, token, attribution: attribution ?? [] });
    }
  };

  type Frame = { curr: number; bondIndex: number; needsClosing: boolean };
  const stack: Frame[] = [{ curr: root, bondIndex: 0, needsClosing: false }];

  while (stack.length > 0) {
    const frame = stack[stack.length - 1];
    const outBonds = mol.getOutDirBonds(frame.curr);

    if (frame.bondIndex === 0) {
      const atom = mol.getAtom(frame.curr);
      const token = atomToSmiles(atom);
      record(token, mol.getAttribution(atom), push(token));
    }

    if (frame.bondIndex < outBonds.length) {
      const bond = outBonds[frame.bondIndex];
      const isLast = frame.bondIndex === outBonds.length - 1;
      frame.bondIndex += 1;

      if (bond.ringBond) {
        const token = bondToSmiles(bond);
        record(token, mol.getAttribution(bond), push(token));
        const key = `${Math.min(bond.src, bond.dst)},${Math.max(bond.src, bond.dst)}`;
        push(ringNumberToSmiles(ringLog.numberFor(key)));
      } else {
        if (!isLast) {
          push('(');
        }
        const token = bondToSmiles(bond);
        record(token, mol.getAttribution(bond), push(token));
        stack.push({ curr: bond.dst, bondIndex: 0, needsClosing: !isLast });
      }
    } else {
      stack.pop();
      if (frame.needsClosing) {
        push(')');
      }
    }
  }
}
