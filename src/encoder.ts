/**
 * SMILES to SELFIES encoder with molecular graph parsing
 */

import { EncoderError, SMILESParserError } from './exceptions.js';
import { MolecularGraph, Atom, DirectedBond } from './mol-graph.js';
import { smilesToMol } from './utils/smiles-parser.js';
import { getSelfiesFromIndex } from './grammar-rules.js';
import { atomToSmiles, bondToSmiles } from './utils/smiles-utils.js';
import type { EncoderOptions, AttributionMap } from './types.js';

/**
 * Translates a SMILES string into its corresponding SELFIES string.
 *
 * This translation is deterministic and does not depend on the current
 * semantic constraints. It preserves the atom order of the input SMILES.
 * Aromatic SMILES are supported by kekulizing them before translation.
 *
 * Not supported: the wildcard `*`, the quadruple bond `$`, chirality other
 * than `@`/`@@`, ring bonds across a dot, and rings or branches whose length
 * cannot be expressed with three index symbols.
 *
 * @param smiles - the SMILES string to be translated
 * @param options.strict - if true (default), throw if the molecule violates
 *   the current semantic constraints
 * @param options.attribute - if true, also return an attribution map
 *   connecting SELFIES tokens to the SMILES tokens that produced them
 * @returns the SELFIES string (or a [selfies, attributionMap] tuple if
 *   `attribute` is true)
 * @throws EncoderError if the SMILES is invalid, cannot be kekulized, or
 *   violates the semantic constraints with `strict` set
 */
export function encoder(smiles: string, options: EncoderOptions & { attribute: true }): [string, AttributionMap];
export function encoder(smiles: string, options?: EncoderOptions & { attribute?: false }): string;
export function encoder(smiles: string, options?: EncoderOptions): string | [string, AttributionMap];
export function encoder(smiles: string, options?: EncoderOptions): string | [string, AttributionMap] {
  if (typeof smiles !== 'string') {
    throw new EncoderError('SMILES input must be a string');
  }

  const attribute = options?.attribute ?? false;
  const strict = options?.strict ?? true;

  let mol: MolecularGraph;
  try {
    mol = smilesToMol(smiles, attribute);
  } catch (error) {
    if (error instanceof SMILESParserError) {
      throw new EncoderError(`failed to parse input\n\tSMILES: ${smiles}\n${error.message}`, error);
    }
    throw error;
  }

  if (!mol.kekulize()) {
    throw new EncoderError(`kekulization failed\n\tSMILES: ${smiles}`);
  }

  if (strict) {
    checkBondConstraints(mol, smiles);
  }

  // Invert chirality of atoms where necessary,
  // such that they are restored when the SELFIES is decoded
  for (const atom of mol.getAtoms()) {
    if (atom.chirality !== null && mol.hasOutRingBond(atom.index!) && shouldInvertChirality(mol, atom)) {
      atom.invertChirality();
    }
  }

  const fragments: string[] = [];
  let attributionMaps: AttributionMap | null = attribute ? [] : null;
  let attributionIndex = 0;

  for (const root of mol.getRoots()) {
    const derived = fragmentToSelfies(mol, smiles, null, root, attributionMaps, attributionIndex);
    attributionIndex += derived.length;
    fragments.push(derived.join(''));
  }

  const selfies = fragments.join('.');
  if (attributionMaps === null) {
    return selfies;
  }
  // Trim attribution map of empty tokens
  attributionMaps = attributionMaps.filter(a => a.token);
  return [selfies, attributionMaps];
}

function checkBondConstraints(mol: MolecularGraph, smiles: string): void {
  const errors: string[] = [];

  for (const atom of mol.getAtoms()) {
    const bondCap = atom.bondingCapacity;
    const bondCount = mol.getBondCount(atom.index!);
    if (bondCount > bondCap) {
      errors.push(
        `\t[${atomToSmiles(atom)} with ${bondCount} bond(s) - a max. of ${bondCap} bond(s) was specified]`
      );
    }
  }

  if (errors.length > 0) {
    throw new EncoderError(
      `input violates the currently-set semantic constraints\n` +
      `\tSMILES: ${smiles}\n\tErrors:\n${errors.join('\n')}`
    );
  }
}

function shouldInvertChirality(mol: MolecularGraph, atom: Atom): boolean {
  const outBonds = mol.getOutDirBonds(atom.index!);

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
  partition[1].sort((a, b) => outBonds[a].dst - outBonds[b].dst);

  // Construct permutation and count its inversions
  const perm = [...partition[0], ...partition[1], ...partition[2]];
  let count = 0;
  for (let i = 0; i < perm.length; i++) {
    for (let j = i + 1; j < perm.length; j++) {
      if (perm[i] > perm[j]) {
        count++;
      }
    }
  }
  return count % 2 !== 0; // if odd permutation, should invert chirality
}

function indexSymbols(smiles: string, kind: string, index: number): string[] {
  const symbols = getSelfiesFromIndex(index);
  if (symbols.length > 3) {
    throw new EncoderError(`${kind} is too long to be represented in SELFIES\n\tSMILES: ${smiles}`);
  }
  return symbols;
}

function fragmentToSelfies(
  mol: MolecularGraph,
  smiles: string,
  bondIntoRoot: DirectedBond | null,
  root: number,
  attributionMaps: AttributionMap | null,
  attributionIndex: number
): string[] {
  const derived: string[] = [];
  const record = (token: string, o: Atom | DirectedBond): void => {
    if (attributionMaps !== null) {
      attributionMaps.push({ index: derived.length - 1 + attributionIndex, token, attribution: mol.getAttribution(o) ?? [] });
    }
  };

  let bondIntoCurr = bondIntoRoot;
  let curr = root;

  for (;;) {
    const currAtom = mol.getAtom(curr);
    const token = atomToSelfies(bondIntoCurr, currAtom);
    derived.push(token);
    record(token, currAtom);

    const outBonds = mol.getOutDirBonds(curr);
    for (let i = 0; i < outBonds.length; i++) {
      const bond = outBonds[i];

      if (bond.ringBond) {
        if (bond.src < bond.dst) {
          continue;
        }

        const revBond = mol.getDirBond(bond.dst, bond.src);
        const ringLen = bond.src - bond.dst;
        const qAsSymbols = indexSymbols(smiles, 'ring', ringLen - 1);
        const ringSymbol = `[${ringBondsToSelfies(revBond, bond)}Ring${qAsSymbols.length}]`;

        derived.push(ringSymbol);
        record(ringSymbol, bond);
        for (const symbol of qAsSymbols) {
          derived.push(symbol);
          record(symbol, bond);
        }
      } else if (i === outBonds.length - 1) {
        bondIntoCurr = bond;
        curr = bond.dst;
      } else {
        // Branch tokens are recorded before the branch symbol is inserted,
        // so shift their indices afterwards
        const start = attributionMaps?.length ?? 0;
        const branch = fragmentToSelfies(mol, smiles, bond, bond.dst, attributionMaps, derived.length + attributionIndex);
        const qAsSymbols = indexSymbols(smiles, 'branch', branch.length - 1);
        const branchSymbol = `[${bondToSelfies(bond, false)}Branch${qAsSymbols.length}]`;
        const end = attributionMaps?.length ?? 0;

        derived.push(branchSymbol);
        for (const symbol of qAsSymbols) {
          derived.push(symbol);
          record(symbol, bond);
        }
        for (let j = start; j < end; j++) {
          attributionMaps![j].index += qAsSymbols.length + 1;
        }
        if (attributionMaps !== null) {
          attributionMaps.push({
            index: derived.length - qAsSymbols.length - 1 + attributionIndex,
            token: branchSymbol,
            attribution: mol.getAttribution(bond) ?? []
          });
        }

        derived.push(...branch);
      }
    }

    // End of chain
    if (outBonds.length === 0 || outBonds[outBonds.length - 1].ringBond) {
      break;
    }
  }

  return derived;
}

function bondToSelfies(bond: DirectedBond, showStereo: boolean = true): string {
  if (!showStereo && bond.order === 1) {
    return '';
  }
  return bondToSmiles(bond);
}

function ringBondsToSelfies(lBond: DirectedBond, rBond: DirectedBond): string {
  if (lBond.order !== 1 || (lBond.stereo === null && rBond.stereo === null)) {
    return bondToSelfies(lBond, false);
  }
  return (lBond.stereo ?? '-') + (rBond.stereo ?? '-');
}

function atomToSelfies(bond: DirectedBond | null, atom: Atom): string {
  const bondChar = bond === null ? '' : bondToSelfies(bond);
  return `[${bondChar}${atomToSmiles(atom, false)}]`;
}
