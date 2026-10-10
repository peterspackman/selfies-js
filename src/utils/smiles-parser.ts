/**
 * SMILES parser for creating molecular graphs
 */

import { MolecularGraph, Atom } from '../mol-graph.js';
import { SMILESParserError } from '../exceptions.js';
import { AROMATIC_SUBSET, ELEMENTS, ORGANIC_SUBSET } from '../constants.js';
import { Attribution } from '../types.js';

const SMILES_BRACKETED_ATOM_PATTERN = new RegExp(
  "^\\[" +                       // opening square bracket [
  "(\\d*)" +                     // isotope number (optional, e.g. 123, 26)
  "([A-Za-z][a-z]?)" +           // element symbol
  "(@{0,2})" +                   // chiral tag (optional, only @ and @@ supported)
  "((?:H\\d?)?)" +               // H count (optional, e.g. H, H0, H3)
  "((?:\\++|-+|[+-]\\d+)?)" +    // charge (optional, e.g. ---, +1, ++)
  "((?::\\d+)?)" +               // atom class (optional, e.g. :12, :1)
  "\\]$"                         // closing square bracket ]
);

export const SMILES_BOND_ORDERS: Record<string, number> = {
  "-": 1, "/": 1, "\\": 1, ":": 1.5, "=": 2, "#": 3
};

export const SMILES_STEREO_BONDS = new Set(["/", "\\"]);

export enum SMILESTokenType {
  ATOM,
  BRANCH,
  RING,
  DOT
}

/**
 * A token in a SMILES string, containing a symbol (atom, branch bracket,
 * ring number, dot) and its preceding bond, if it exists (e.g. =C, %12, #N).
 */
export class SMILESToken {
  constructor(
    public bondIdx: number | null,
    public startIdx: number,
    public endIdx: number,
    public tokenType: SMILESTokenType,
    public token: string
  ) {}

  extractBondChar(smiles: string): string | null {
    return this.bondIdx === null ? null : smiles[this.bondIdx];
  }

  extractSymbol(smiles: string): string {
    return smiles.slice(this.startIdx, this.endIdx);
  }

  toString(): string {
    return this.token;
  }
}

const isAlpha = (c: string): boolean => /^[A-Za-z]$/.test(c);
const isDigit = (c: string): boolean => /^[0-9]$/.test(c);

/**
 * Splits a SMILES string into its tokens.
 *
 * @param smiles - the input SMILES string
 * @yields the tokens of the input SMILES one-by-one with order preserved
 */
export function* tokenizeSmiles(smiles: string): Generator<SMILESToken> {
  let i = 0;
  while (i < smiles.length) {
    if (smiles[i] === '.') {
      yield new SMILESToken(null, i, i + 1, SMILESTokenType.DOT, '.');
      i += 1;
      continue;
    }

    let bondIdx: number | null = null;
    if (smiles[i] in SMILES_BOND_ORDERS) {
      bondIdx = i;
      i += 1;
    }

    let token: SMILESToken;
    if (i === smiles.length) {
      throw new SMILESParserError(smiles, "hanging bond", i - 1);
    } else if (isAlpha(smiles[i])) { // organic subset elements
      const twoChar = smiles.slice(i, i + 2);
      if (twoChar === 'Br' || twoChar === 'Cl') { // two-letter elements
        token = new SMILESToken(bondIdx, i, i + 2, SMILESTokenType.ATOM, twoChar);
      } else { // one-letter elements (e.g. C, N, ...)
        token = new SMILESToken(bondIdx, i, i + 1, SMILESTokenType.ATOM, smiles[i]);
      }
    } else if (smiles[i] === '[') { // atoms encased in brackets (e.g. [NH])
      const rIdx = smiles.indexOf(']', i + 1);
      if (rIdx === -1) {
        throw new SMILESParserError(smiles, "hanging bracket [", i);
      }
      token = new SMILESToken(bondIdx, i, rIdx + 1, SMILESTokenType.ATOM, smiles.slice(i, rIdx + 1));
    } else if (smiles[i] === '(' || smiles[i] === ')') { // branch brackets
      if (bondIdx !== null) {
        throw new SMILESParserError(smiles, "hanging bond", bondIdx);
      }
      token = new SMILESToken(null, i, i + 1, SMILESTokenType.BRANCH, smiles[i]);
    } else if (isDigit(smiles[i])) { // one-digit ring number
      token = new SMILESToken(bondIdx, i, i + 1, SMILESTokenType.RING, smiles[i]);
    } else if (smiles[i] === '%') { // two-digit ring number (e.g. %12)
      const rnum = smiles.slice(i + 1, i + 3);
      if (!/^\d\d$/.test(rnum)) {
        throw new SMILESParserError(smiles, `invalid ring number '%${rnum}'`, i);
      }
      token = new SMILESToken(bondIdx, i, i + 3, SMILESTokenType.RING, smiles.slice(i, i + 3));
    } else {
      throw new SMILESParserError(smiles, `unrecognized symbol '${smiles[i]}'`, i);
    }

    yield token;
    i = token.endIdx;
  }
}

/**
 * Reads an atom from its SMILES representation.
 *
 * @param atomSymbol - a SMILES atom symbol
 * @returns the atom that the input symbol represents, or null if invalid
 */
export function smilesToAtom(atomSymbol: string): Atom | null {
  if (atomSymbol.startsWith('[') && atomSymbol.endsWith(']')) {
    // continue below
  } else if (ORGANIC_SUBSET.has(atomSymbol)) { // e.g. C, N, O, ...
    return new Atom(atomSymbol, false);
  } else if (AROMATIC_SUBSET.has(atomSymbol)) { // e.g. c, n, o, ...
    return new Atom(capitalize(atomSymbol), true);
  } else {
    return null;
  }

  // e.g. [C], [C@@H], [O-], ...
  const m = SMILES_BRACKETED_ATOM_PATTERN.exec(atomSymbol);
  if (m === null) {
    return null;
  }
  const [, isotopeStr, elementStr, chiralityStr, hCountStr, chargeStr] = m;

  const isotope = isotopeStr === '' ? null : parseInt(isotopeStr, 10);
  const isAromatic = elementStr === elementStr.toLowerCase() && AROMATIC_SUBSET.has(elementStr);
  const element = capitalize(elementStr);
  if (!ELEMENTS.has(element)) {
    return null;
  }
  const chirality = chiralityStr === '' ? null : chiralityStr;

  let hCount = 0;
  if (hCountStr !== '') {
    const s = hCountStr.slice(1); // HXXX -> XXX
    hCount = s === '' ? 1 : parseInt(s, 10);
  }

  let charge = 0;
  if (chargeStr !== '') {
    if (isDigit(chargeStr[chargeStr.length - 1])) { // (+/-)XXX
      charge = parseInt(chargeStr.slice(1), 10);
    } else { // +++... or ---...
      charge = chargeStr.length;
    }
    charge *= chargeStr[0] === '+' ? 1 : -1;
  }

  return new Atom(element, isAromatic, isotope, chirality, hCount, charge);
}

/**
 * Reads a bond from its SMILES representation.
 *
 * @param bondChar - a SMILES bond symbol
 * @returns the order and stereochemical specification of the bond
 */
export function smilesToBond(bondChar: string | null): [number, string | null] {
  const order = bondChar !== null && bondChar in SMILES_BOND_ORDERS ? SMILES_BOND_ORDERS[bondChar] : 1;
  const stereo = bondChar !== null && SMILES_STEREO_BONDS.has(bondChar) ? bondChar : null;
  return [order, stereo];
}

function capitalize(s: string): string {
  return s.charAt(0).toUpperCase() + s.slice(1).toLowerCase();
}

/**
 * Reads a molecular graph from a SMILES string.
 *
 * @param smiles - the input SMILES string
 * @param attributable - if the molecular graph needs to include attributions
 * @returns a molecular graph that the input SMILES string represents
 * @throws SMILESParserError if the input SMILES is invalid
 */
export function smilesToMol(smiles: string, attributable: boolean = false): MolecularGraph {
  if (smiles === '') {
    throw new SMILESParserError(smiles, "empty SMILES", 0);
  }

  const mol = new MolecularGraph(attributable);
  const tokens = [...tokenizeSmiles(smiles)];
  const state = { pos: 0, i: 0 };
  while (state.pos < tokens.length) {
    deriveMolFromTokens(mol, smiles, tokens, state);
  }
  return mol;
}

interface RingLogEntry {
  token: SMILESToken;
  atom: Atom;
  pos: number;
}

function deriveMolFromTokens(
  mol: MolecularGraph,
  smiles: string,
  tokens: SMILESToken[],
  state: { pos: number; i: number }
): void {
  let tok: SMILESToken | null = null;
  const prevStack: (Atom | null)[] = [null]; // previous atom on the current chain
  const branchStack: SMILESToken[] = []; // open branches
  const ringLog = new Map<string, RingLogEntry>(); // hanging ring numbers
  const nAtomsBefore = mol.length;
  let chainStart = true;

  while (state.pos < tokens.length) {
    tok = tokens[state.pos++];
    const bondChar = tok.extractBondChar(smiles);
    const symbol = tok.extractSymbol(smiles);
    const prevAtom = prevStack[prevStack.length - 1];

    if (tok.tokenType === SMILESTokenType.DOT) {
      break;
    } else if (tok.tokenType === SMILESTokenType.ATOM) {
      const curr = smilesToAtom(symbol);
      if (curr === null) {
        throw new SMILESParserError(smiles, `invalid atom symbol '${symbol}'`, tok.startIdx);
      }
      attachAtom(mol, bondChar, curr, prevAtom, state, tok);
      prevStack[prevStack.length - 1] = curr;
      chainStart = false;
    } else if (chainStart) {
      throw new SMILESParserError(smiles, "SMILES chain begins with non-atom", tok.startIdx);
    } else if (tok.tokenType === SMILESTokenType.BRANCH) {
      if (symbol === '(') {
        branchStack.push(tok);
        prevStack.push(prevAtom);
        chainStart = true;
      } else {
        if (branchStack.length === 0) {
          throw new SMILESParserError(smiles, "hanging ')' bracket", tok.startIdx);
        }
        branchStack.pop();
        prevStack.pop();
      }
    } else { // RING
      const entry = ringLog.get(symbol);
      if (entry === undefined) {
        const lpos = mol.addPlaceholderBond(prevAtom!.index!);
        ringLog.set(symbol, { token: tok, atom: prevAtom!, pos: lpos });
      } else {
        ringLog.delete(symbol);
        makeRingBonds(mol, smiles, entry, tok, prevAtom!);
      }
    }
    state.i += 1;
  }

  if (mol.length === nAtomsBefore) {
    const errIdx = (tok === null ? smiles.length : tok.startIdx) - 1;
    throw new SMILESParserError(smiles, "empty SMILES fragment", errIdx);
  }

  if (branchStack.length > 0) {
    throw new SMILESParserError(smiles, "hanging '(' bracket", branchStack[branchStack.length - 1].startIdx);
  }

  if (ringLog.size > 0) {
    const [rnum, entry] = [...ringLog.entries()][ringLog.size - 1];
    throw new SMILESParserError(smiles, `hanging ring number '${rnum}'`, entry.token.startIdx);
  }
}

function attachAtom(
  mol: MolecularGraph,
  bondChar: string | null,
  atom: Atom,
  prevAtom: Atom | null,
  state: { i: number },
  tok: SMILESToken
): void {
  const isRoot = prevAtom === null;
  if (bondChar) {
    state.i += 1;
  }
  mol.addAtom(atom, isRoot);
  mol.addAttribution(atom, [new Attribution(state.i, tok.token)]);

  if (!isRoot) {
    const [bondOrder, stereo] = smilesToBond(bondChar);
    let order = bondOrder;
    if (prevAtom.isAromatic && atom.isAromatic && bondChar === null) {
      order = 1.5; // handle implicit aromatic bonds, e.g. cc
    }
    const bond = mol.addBond(prevAtom.index!, atom.index!, order, stereo);
    mol.addAttribution(bond, [new Attribution(state.i, tok.token)]);
  }
}

function makeRingBonds(
  mol: MolecularGraph,
  smiles: string,
  left: RingLogEntry,
  rtoken: SMILESToken,
  ratom: Atom
): void {
  const latom = left.atom;
  const ltoken = left.token;

  if (latom.index === ratom.index) {
    throw new SMILESParserError(smiles, "ring bond from an atom to itself", ltoken.startIdx);
  }
  if (mol.hasBond(latom.index!, ratom.index!)) {
    throw new SMILESParserError(smiles, "ring bond specified between already-bonded atoms", ltoken.startIdx);
  }

  const lbondChar = ltoken.extractBondChar(smiles);
  const rbondChar = rtoken.extractBondChar(smiles);

  // Check that ring bonds match; swap so that if bonds[0] is null, so is bonds[1]
  const bonds = lbondChar === null ? [rbondChar, lbondChar] : [lbondChar, rbondChar];
  const bothStereo = bonds.every(b => b !== null && SMILES_STEREO_BONDS.has(b));
  if (!(bonds[0] === bonds[1] || bonds[1] === null || bothStereo)) {
    throw new SMILESParserError(smiles, "mismatched ring bonds", ltoken.startIdx);
  }

  const [lorder, lstereo] = smilesToBond(lbondChar);
  const [rorder, rstereo] = smilesToBond(rbondChar);
  let order = Math.max(lorder, rorder);
  if (latom.isAromatic && ratom.isAromatic && lbondChar === null && rbondChar === null) {
    order = 1.5; // handle implicit aromatic bonds, e.g. c1ccccc1
  }

  mol.addRingBond(
    latom.index!, ratom.index!,
    order,
    lstereo, rstereo,
    left.pos
  );
}
