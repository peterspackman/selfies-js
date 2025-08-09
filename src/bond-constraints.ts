/**
 * Bond constraints and alphabet management for SELFIES
 */

import { ELEMENTS, INDEX_ALPHABET } from './constants.js';
import type { SemanticConstraints } from './types.js';
import { clearCaches as clearGrammarCaches } from './grammar-rules.js';

const DEFAULT_CONSTRAINTS: SemanticConstraints = {
  "H": 1, "F": 1, "Cl": 1, "Br": 1, "I": 1,
  "B": 3, "B+1": 2, "B-1": 4,
  "O": 2, "O+1": 3, "O-1": 1,
  "N": 3, "N+1": 4, "N-1": 2,
  "C": 4, "C+1": 3, "C-1": 3,
  "P": 5, "P+1": 4, "P-1": 6,
  "S": 6, "S+1": 5, "S-1": 5,
  "?": 8
};

const PRESET_CONSTRAINTS: Record<string, SemanticConstraints> = {
  "default": { ...DEFAULT_CONSTRAINTS },
  "octet_rule": { 
    ...DEFAULT_CONSTRAINTS,
    "S": 2, "S+1": 3, "S-1": 1, 
    "P": 3, "P+1": 4, "P-1": 2 
  },
  "hypervalent": { 
    ...DEFAULT_CONSTRAINTS,
    "Cl": 7, "Br": 7, "I": 7, "N": 5 
  }
};

let currentConstraints = { ...PRESET_CONSTRAINTS.default };

// Cache for expensive operations
let semanticRobustAlphabetCache: Set<string> | null = null;
const bondingCapacityCache = new Map<string, number>();

/**
 * Returns the preset semantic constraints with the given name.
 * 
 * @param name - The preset name: 'default', 'octet_rule', or 'hypervalent'
 * @returns The preset constraints as a dictionary mapping atoms to bonding capacities
 */
export function getPresetConstraints(name: string): SemanticConstraints {
  if (!(name in PRESET_CONSTRAINTS)) {
    throw new Error(`unrecognized preset name '${name}'`);
  }
  return { ...PRESET_CONSTRAINTS[name] };
}

/**
 * Returns the semantic constraints that SELFIES is currently operating on.
 * 
 * @returns The current semantic constraints as a dictionary
 */
export function getSemanticConstraints(): SemanticConstraints {
  return { ...currentConstraints };
}

/**
 * Updates the semantic constraints that SELFIES operates on.
 * 
 * @param bondConstraints - Either a preset name or a constraints dictionary
 */
export function setSemanticConstraints(
  bondConstraints: string | SemanticConstraints = "default"
): void {
  if (typeof bondConstraints === 'string') {
    currentConstraints = getPresetConstraints(bondConstraints);
  } else if (typeof bondConstraints === 'object' && bondConstraints !== null) {
    if (!('?' in bondConstraints)) {
      throw new Error("bond_constraints missing '?' as a key");
    }

    for (const [key, value] of Object.entries(bondConstraints)) {
      const j = Math.max(key.indexOf('+'), key.indexOf('-'));
      let valid = false;
      
      if (key === '?') {
        valid = true;
      } else if (j === -1) {
        valid = ELEMENTS.has(key);
      } else {
        const element = key.slice(0, j);
        const chargeStr = key.slice(j + 1);
        valid = ELEMENTS.has(element) && /^\d+$/.test(chargeStr);
      }
      
      if (!valid) {
        throw new Error(`invalid key '${key}' in bond_constraints`);
      }

      if (!Number.isInteger(value) || value < 0) {
        throw new Error(`invalid value at bond_constraints['${key}'] = ${value}`);
      }
    }

    currentConstraints = { ...bondConstraints };
  } else {
    throw new Error("bond_constraints must be a string or object");
  }

  clearCaches();
}

function clearCaches(): void {
  semanticRobustAlphabetCache = null;
  bondingCapacityCache.clear();
  clearGrammarCaches();
}

/**
 * Returns a subset of all SELFIES symbols that are constrained
 * by SELFIES under the current semantic constraints.
 * 
 * @returns A set of SELFIES symbols that are semantically constrained
 */
export function getSemanticRobustAlphabet(): Set<string> {
  if (semanticRobustAlphabetCache !== null) {
    return new Set(semanticRobustAlphabetCache);
  }

  const alphabetSubset = new Set<string>();
  const bonds: Record<string, number> = { "": 1, "=": 2, "#": 3 };

  for (const [atom, capacity] of Object.entries(currentConstraints)) {
    for (const [bondSymbol, multiplicity] of Object.entries(bonds)) {
      if (multiplicity > capacity || atom === "?") {
        continue;
      }
      const symbol = `[${bondSymbol}${atom}]`;
      alphabetSubset.add(symbol);
    }
  }

  for (let i = 1; i <= 3; i++) {
    alphabetSubset.add(`[Ring${i}]`);
    alphabetSubset.add(`[=Ring${i}]`);
    alphabetSubset.add(`[Branch${i}]`);
    alphabetSubset.add(`[=Branch${i}]`);
    alphabetSubset.add(`[#Branch${i}]`);
  }

  for (const symbol of INDEX_ALPHABET) {
    // Parse the symbol to check if it satisfies constraints
    const match = symbol.match(/^\[([=#]?)([A-Z][a-z]?)\]$/);
    if (match) {
      const [, bondChar, atom] = match;
      const multiplicity = bondChar === '=' ? 2 : bondChar === '#' ? 3 : 1;
      const capacity = currentConstraints[atom] ?? currentConstraints["?"];
      
      if (multiplicity <= capacity) {
        alphabetSubset.add(symbol);
      }
    } else {
      alphabetSubset.add(symbol);
    }
  }

  semanticRobustAlphabetCache = new Set(alphabetSubset);
  return new Set(alphabetSubset);
}

/**
 * Returns the bonding capacity of a given atom, under the current
 * semantic constraints.
 * 
 * @param element - The element of the input atom
 * @param charge - The charge of the input atom
 * @returns The bonding capacity of the input atom
 */
export function getBondingCapacity(element: string, charge: number): number {
  let key = element;
  if (charge !== 0) {
    key += (charge > 0 ? '+' : '') + charge.toString();
  }

  const cacheKey = key;
  if (bondingCapacityCache.has(cacheKey)) {
    return bondingCapacityCache.get(cacheKey)!;
  }

  const capacity = currentConstraints[key] ?? currentConstraints["?"];
  bondingCapacityCache.set(cacheKey, capacity);
  return capacity;
}