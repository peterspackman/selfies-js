/**
 * Grammar rules for processing SELFIES symbols
 */

import { ELEMENTS, INDEX_ALPHABET, INDEX_CODE, ORGANIC_SUBSET } from './constants.js';
import { Atom } from './mol-graph.js';
import type { State } from './types.js';

// Cache for performance optimization
const PROCESS_ATOM_CACHE = new Map<string, [[number, string | null], () => Atom] | null>();
const PROCESS_BRANCH_CACHE = new Map<string, [number, number]>();
const PROCESS_RING_CACHE = new Map<string, [number, number, [string | null, string | null]]>();

// SELFIES atom pattern matching - updated to match Python exactly  
const SELFIES_ATOM_PATTERN = new RegExp("^\\[([=#/\\\\]?)(\\d*)([A-Z][a-z]?)([@]{0,2})((?:[H]\\d)?)((?:[+-]\\d+)?)\\]$");

export function processAtomSymbol(symbol: string): [[number, string | null], Atom] | null {
  let output = PROCESS_ATOM_CACHE.get(symbol);
  if (output === undefined) {
    output = processAtomSelfiesNoCache(symbol);
    PROCESS_ATOM_CACHE.set(symbol, output);
  }

  if (output === null) {
    return null;
  }

  const [bondInfo, atomFactory] = output;
  const atom = atomFactory();
  
  if (atom.bondingCapacity < 0) {
    return null;
  }
  
  return [bondInfo, atom];
}

export function processBranchSymbol(symbol: string): [number, number] | null {
  return PROCESS_BRANCH_CACHE.get(symbol) || null;
}

export function processRingSymbol(symbol: string): [number, number, [string | null, string | null]] | null {
  return PROCESS_RING_CACHE.get(symbol) || null;
}

export function nextAtomState(
  bondOrder: number, 
  bondCap: number, 
  state: State
): [number, State] {
  if (state === 0) {
    bondOrder = 0;
  }

  bondOrder = Math.min(bondOrder, state || 0, bondCap);
  const bondsLeft = bondCap - bondOrder;
  const nextState = bondsLeft === 0 ? null : bondsLeft;
  return [bondOrder, nextState];
}

export function nextBranchState(
  branchType: number, 
  state: State
): [number, State] {
  if (!(1 <= branchType && branchType <= 3)) {
    throw new Error(`Invalid branch type: ${branchType}`);
  }
  if (state === null || state <= 1) {
    throw new Error(`Invalid state for branching: ${state}`);
  }

  const branchInitState = Math.min(state - 1, branchType);
  const nextState = state - branchInitState;
  return [branchInitState, nextState === 0 ? null : nextState];
}

export function nextRingState(
  ringType: number, 
  state: State
): [number, State] {
  if (state === null || state <= 0) {
    throw new Error(`Invalid state for ring: ${state}`);
  }

  const bondOrder = Math.min(ringType, state);
  const bondsLeft = state - bondOrder;
  const nextState = bondsLeft === 0 ? null : bondsLeft;
  return [bondOrder, nextState];
}

export function getIndexFromSelfies(...symbols: string[]): number {
  let index = 0;
  const reversedSymbols = symbols.slice().reverse();
  
  for (let i = 0; i < reversedSymbols.length; i++) {
    const symbol = reversedSymbols[i];
    const code = INDEX_CODE[symbol] || 0;
    index += code * Math.pow(Object.keys(INDEX_CODE).length, i);
  }
  
  return index;
}

export function getSelfiesFromIndex(index: number): string[] {
  if (index < 0) {
    throw new Error('Index must be non-negative');
  } else if (index === 0) {
    return [INDEX_ALPHABET[0]];
  }

  const symbols: string[] = [];
  const base = INDEX_ALPHABET.length;
  let remaining = index;
  
  while (remaining > 0) {
    symbols.push(INDEX_ALPHABET[remaining % base]);
    remaining = Math.floor(remaining / base);
  }
  
  return symbols.reverse();
}

function smilesToBond(bondChar: string): [number, string | null] {
  switch (bondChar) {
    case '':
      return [1, null];
    case '=':
      return [2, null];
    case '#':
      return [3, null];
    case '/':
      return [1, '/'];
    case '\\':
      return [1, '\\'];
    case '-':
      return [1, null];
    default:
      throw new Error(`Unknown bond character: ${bondChar}`);
  }
}

function processAtomSelfiesNoCache(symbol: string): [[number, string | null], () => Atom] | null {
  const match = SELFIES_ATOM_PATTERN.exec(symbol);
  if (!match) {
    return null;
  }

  const [, bondChar, isotope, element, chirality, hCount, charge] = match;

  const innerSymbol = symbol.slice(1 + bondChar.length, -1);
  if (ORGANIC_SUBSET.has(innerSymbol)) {
    const atomFactory = () => new Atom(element, false);
    return [smilesToBond(bondChar), atomFactory];
  }

  if (!ELEMENTS.has(element)) {
    return null;
  }

  const isotopeNum = isotope === '' ? null : parseInt(isotope, 10);

  const chiralityStr = chirality === '' ? null : chirality;

  let hCountNum: number;
  if (hCount === '') {
    hCountNum = 0;
  } else {
    hCountNum = parseInt(hCount.slice(1), 10);
  }

  let chargeNum: number;
  if (charge === '') {
    chargeNum = 0;
  } else {
    const chargeValue = parseInt(charge.slice(1), 10);
    chargeNum = charge[0] === '+' ? chargeValue : -chargeValue;
  }

  const atomFactory = () => new Atom(
    element,
    false,
    isotopeNum,
    chiralityStr,
    hCountNum,
    chargeNum
  );

  return [smilesToBond(bondChar), atomFactory];
}

function buildAtomCache(): void {
  const commonSymbols = [
    "[#C+1]", "[#C-1]", "[#C]", "[#N+1]", "[#N]", "[#O+1]", "[#P+1]",
    "[#P-1]", "[#P]", "[#S+1]", "[#S-1]", "[#S]", "[=C+1]", "[=C-1]",
    "[=C]", "[=N+1]", "[=N-1]", "[=N]", "[=O+1]", "[=O]", "[=P+1]",
    "[=P-1]", "[=P]", "[=S+1]", "[=S-1]", "[=S]", "[Br]", "[C+1]", "[C-1]",
    "[C]", "[Cl]", "[F]", "[H]", "[I]", "[N+1]", "[N-1]", "[N]", "[O+1]",
    "[O-1]", "[O]", "[P+1]", "[P-1]", "[P]", "[S+1]", "[S-1]", "[S]"
  ];

  for (const symbol of commonSymbols) {
    PROCESS_ATOM_CACHE.set(symbol, processAtomSelfiesNoCache(symbol));
  }
}

function buildBranchCache(): void {
  for (let L = 1; L <= 3; L++) {
    for (const bondChar of ['', '=', '#']) {
      const symbol = `[${bondChar}Branch${L}]`;
      const [bondOrder] = smilesToBond(bondChar);
      PROCESS_BRANCH_CACHE.set(symbol, [bondOrder, L]);
    }
  }
}

function buildRingCache(): void {
  for (let L = 1; L <= 3; L++) {
    // [RingL], [=RingL], [#RingL]
    for (const bondChar of ['', '=', '#']) {
      const symbol = `[${bondChar}Ring${L}]`;
      const [order, stereo] = smilesToBond(bondChar);
      PROCESS_RING_CACHE.set(symbol, [order, L, [stereo, stereo]]);
    }

    // Stereo ring bonds: [\-Ring1], [/-Ring1], etc.
    const stereoChars = ['-', '/', '\\'];
    for (const lChar of stereoChars) {
      for (const rChar of stereoChars) {
        if (lChar === '-' && rChar === '-') continue;
        
        const symbol = `[${lChar}${rChar}Ring${L}]`;
        const [, lStereo] = smilesToBond(lChar);
        const [, rStereo] = smilesToBond(rChar);
        PROCESS_RING_CACHE.set(symbol, [1, L, [lStereo, rStereo]]);
      }
    }
  }
}

buildAtomCache();
buildBranchCache();
buildRingCache();

export function clearCaches(): void {
  PROCESS_ATOM_CACHE.clear();
  PROCESS_BRANCH_CACHE.clear();
  PROCESS_RING_CACHE.clear();
  buildAtomCache();
  buildBranchCache();
  buildRingCache();
}