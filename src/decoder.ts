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
import type { AttributionMap, DecoderOptions, State } from './types.js';
import { modernizeSymbol, warnDeprecatedCompatibility } from './utils/compatibility-utils.js';

type RingBondInfo = [number, [string | null, string | null]];
type PendingRing = [Atom, Atom, RingBondInfo];
type SymbolIterator = Iterator<[number, string]>;

/**
 * Translates a SELFIES string into its corresponding SMILES string.
 *
 * This translation is deterministic but depends on the current semantic
 * constraints. The output SMILES string is guaranteed to be syntactically
 * correct and to represent a molecule that obeys the semantic constraints.
 *
 * @param selfies - the SELFIES string to be translated
 * @param options.attribute - if true, also return an attribution map
 *   connecting SMILES tokens to the SELFIES tokens that produced them
 * @param options.compatible - if true, accept symbols from SELFIES <v2
 * @returns the SMILES string (or a [smiles, attributionMap] tuple if
 *   `attribute` is true)
 * @throws DecoderError if the input SELFIES string is malformed
 */
export function decoder(selfies: string, options: DecoderOptions & { attribute: true }): [string, AttributionMap];
export function decoder(selfies: string, options?: DecoderOptions & { attribute?: false }): string;
export function decoder(selfies: string, options?: DecoderOptions): string | [string, AttributionMap];
export function decoder(selfies: string, options?: DecoderOptions): string | [string, AttributionMap] {
  if (typeof selfies !== 'string') {
    throw new DecoderError('SELFIES input must be a string');
  }

  const attribute = options?.attribute ?? false;
  const compatible = options?.compatible ?? false;

  if (compatible) {
    warnDeprecatedCompatibility();
  }

  const mol = new MolecularGraph(attribute);
  const rings: PendingRing[] = [];
  let attributionIndex = 0;

  for (const fragment of selfies.split('.')) {
    const n = deriveMolFromSymbols(
      enumerate(tokenizeSelfies(fragment, compatible)),
      mol,
      selfies,
      Infinity,
      0,
      null,
      rings,
      attribute ? [] : null,
      attributionIndex
    );
    attributionIndex += n;
  }

  formRingsBilocally(mol, rings);
  return molToSmiles(mol, attribute);
}

function* tokenizeSelfies(selfies: string, compatible: boolean): Generator<string> {
  try {
    for (let symbol of splitSelfies(selfies)) {
      if (symbol === '[nop]') {
        continue;
      }
      if (compatible) {
        symbol = modernizeSymbol(symbol);
      }
      yield symbol;
    }
  } catch (error) {
    if (error instanceof DecoderError) throw error;
    throw new DecoderError(error instanceof Error ? error.message : String(error));
  }
}

function* enumerate<T>(iter: Iterable<T>): Generator<[number, T]> {
  let i = 0;
  for (const item of iter) {
    yield [i++, item];
  }
}

function deriveMolFromSymbols(
  symbolIter: SymbolIterator,
  mol: MolecularGraph,
  selfies: string,
  maxDerive: number,
  initState: number,
  rootAtom: Atom | null,
  rings: PendingRing[],
  attributeStack: Attribution[] | null,
  attributionIndex: number
): number {
  let nDerived = 0;
  let state: State = initState;
  let prevAtom = rootAtom;

  const attributionFor = (index: number, symbol: string): Attribution[] | null =>
    attributeStack ? [...attributeStack, new Attribution(index + attributionIndex, symbol)] : null;

  while (state !== null && nDerived < maxDerive) {
    const next = symbolIter.next();
    if (next.done) break;
    nDerived++;

    const [index, symbol] = next.value;
    let nextState: State;

    // Case 1: Branch symbol (e.g. [Branch1])
    if (symbol.slice(-4, -2) === 'ch') {
      const output = processBranchSymbol(symbol);
      if (output === null) {
        raiseDecoderError(selfies, symbol);
      }
      const [btype, n] = output;

      if (state <= 1) {
        nextState = state;
      } else {
        const [binitState, branchNextState] = nextBranchState(btype, state);
        nextState = branchNextState;

        const Q = readIndexFromSelfies(symbolIter, n);
        nDerived += n + deriveMolFromSymbols(
          symbolIter, mol, selfies, Q + 1,
          binitState, prevAtom, rings,
          attributionFor(index, symbol),
          attributionIndex
        );
      }
    }
    // Case 2: Ring symbol (e.g. [Ring2])
    else if (symbol.slice(-4, -2) === 'ng') {
      const output = processRingSymbol(symbol);
      if (output === null) {
        raiseDecoderError(selfies, symbol);
      }
      const [ringType, n, stereo] = output;

      if (state === 0) {
        nextState = state;
      } else {
        const [ringOrder, ringNextState] = nextRingState(ringType, state);
        nextState = ringNextState;

        const Q = readIndexFromSelfies(symbolIter, n);
        nDerived += n;
        const ratom = prevAtom!;
        const lidx = Math.max(0, ratom.index! - (Q + 1));
        rings.push([mol.getAtom(lidx), ratom, [ringOrder, stereo]]);
      }
    }
    // Case 3: [epsilon]
    else if (symbol.includes('eps')) {
      nextState = state === 0 ? 0 : null;
    }
    // Case 4: Regular atom symbol (e.g. [N], [=C], [F])
    else {
      const output = processAtomSymbol(symbol);
      if (output === null) {
        raiseDecoderError(selfies, symbol);
      }
      const [[bondOrder, stereo], atom] = output;

      const [finalBondOrder, atomNextState] = nextAtomState(bondOrder, atom.bondingCapacity, state);
      nextState = atomNextState;

      if (finalBondOrder === 0) {
        if (state === 0) {
          mol.addAtom(atom, true);
          mol.addAttribution(atom, attributionFor(index, symbol));
        }
      } else {
        mol.addAtom(atom);
        mol.addAttribution(atom, attributionFor(index, symbol));
        const bond = mol.addBond(prevAtom!.index!, atom.index!, finalBondOrder, stereo);
        mol.addAttribution(bond, attributionFor(index, symbol));
      }
      prevAtom = atom;
    }

    state = nextState;
  }

  // Consume remaining tokens
  while (nDerived < maxDerive) {
    if (symbolIter.next().done) break;
    nDerived++;
  }

  return nDerived;
}

function readIndexFromSelfies(symbolIter: SymbolIterator, nSymbols: number): number {
  const indexSymbols: (string | null)[] = [];
  for (let i = 0; i < nSymbols; i++) {
    const next = symbolIter.next();
    indexSymbols.push(next.done ? null : next.value[1]);
  }
  return getIndexFromSelfies(...indexSymbols);
}

function formRingsBilocally(mol: MolecularGraph, rings: PendingRing[]): void {
  const ringsMade = new Array<number>(mol.length).fill(0);

  for (const [latom, ratom, [ringOrder, [lstereo, rstereo]]] of rings) {
    const lidx = latom.index!;
    const ridx = ratom.index!;

    if (lidx === ridx) { // ring to the same atom forbidden
      continue;
    }

    const lfree = latom.bondingCapacity - mol.getBondCount(lidx);
    const rfree = ratom.bondingCapacity - mol.getBondCount(ridx);
    if (lfree <= 0 || rfree <= 0) {
      continue; // no room for ring bond
    }
    const order = Math.min(ringOrder, lfree, rfree);

    if (mol.hasBond(lidx, ridx)) {
      const bond = mol.getDirBond(lidx, ridx);
      const newOrder = Math.min(order + bond.order, 3);
      mol.updateBondOrder(lidx, ridx, newOrder);
    } else {
      mol.addRingBond(lidx, ridx, order, lstereo, rstereo, ringsMade[lidx], ringsMade[ridx]);
      ringsMade[lidx] += 1;
      ringsMade[ridx] += 1;
    }
  }
}

function raiseDecoderError(selfies: string, symbol: string): never {
  throw new DecoderError(`invalid symbol '${symbol}'\n\tSELFIES: ${selfies}`);
}
