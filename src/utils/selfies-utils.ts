/**
 * Utility functions for working with SELFIES strings
 */

/**
 * Returns the number of symbols in a given SELFIES string.
 * 
 * @param selfies - A SELFIES string
 * @returns The symbol length of the SELFIES string
 * 
 * @example
 * ```typescript
 * import { lenSelfies } from 'selfies';
 * console.log(lenSelfies("[C][=C][F].[C]")); // 5
 * ```
 */
export function lenSelfies(selfies: string): number {
  return (selfies.match(/\[/g) || []).length + (selfies.match(/\./g) || []).length;
}

/**
 * Tokenizes a SELFIES string into its individual symbols.
 * 
 * @param selfies - A SELFIES string
 * @yields The symbols of the SELFIES string one-by-one with order preserved
 * 
 * @example
 * ```typescript
 * import { splitSelfies } from 'selfies';
 * console.log([...splitSelfies("[C][=C][F].[C]")]); // ['[C]', '[=C]', '[F]', '.', '[C]']
 * ```
 */
export function* splitSelfies(selfies: string): Generator<string> {
  let leftIdx = selfies.indexOf('[');
  
  while (leftIdx >= 0 && leftIdx < selfies.length) {
    const rightIdx = selfies.indexOf(']', leftIdx + 1);
    if (rightIdx === -1) {
      throw new Error("malformed SELFIES string, hanging '[' bracket");
    }
    
    const nextSymbol = selfies.slice(leftIdx, rightIdx + 1);
    yield nextSymbol;
    
    leftIdx = rightIdx + 1;
    if (selfies.slice(leftIdx, leftIdx + 1) === '.') {
      yield '.';
      leftIdx += 1;
    }
    
    leftIdx = selfies.indexOf('[', leftIdx);
  }
}

/**
 * Constructs an alphabet from an iterable of SELFIES strings.
 * 
 * The returned alphabet is the set of all symbols that appear in the
 * SELFIES strings from the input iterable, minus the dot '.' symbol.
 * 
 * @param selfiesIter - An iterable of SELFIES strings
 * @returns An alphabet of SELFIES symbols, built from the input iterable
 * 
 * @example
 * ```typescript
 * import { getAlphabetFromSelfies } from 'selfies';
 * const selfiesList = ["[C][F][O]", "[C].[O]", "[F][F]"];
 * const alphabet = getAlphabetFromSelfies(selfiesList);
 * console.log([...alphabet].sort()); // ['[C]', '[F]', '[O]']
 * ```
 */
export function getAlphabetFromSelfies(selfiesIter: Iterable<string>): Set<string> {
  const alphabet = new Set<string>();
  
  for (const s of selfiesIter) {
    for (const symbol of splitSelfies(s)) {
      alphabet.add(symbol);
    }
  }
  
  alphabet.delete('.');
  return alphabet;
}