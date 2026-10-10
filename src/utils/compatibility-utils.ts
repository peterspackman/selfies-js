/**
 * Backward compatibility utilities for legacy (<v2) SELFIES symbols
 */

import { atomToSmiles } from './smiles-utils.js';
import { smilesToAtom } from './smiles-parser.js';

function buildUpdateTable(): Record<string, string> {
  const table: Record<string, string> = {};
  for (let L = 1; L <= 3; L++) {
    table[`[Branch${L}_1]`] = `[Branch${L}]`;
    table[`[Branch${L}_2]`] = `[=Branch${L}]`;
    table[`[Branch${L}_3]`] = `[#Branch${L}]`;
    table[`[Expl=Ring${L}]`] = `[=Ring${L}]`;
    table[`[Expl#Ring${L}]`] = `[#Ring${L}]`;
    table[`[Expl/Ring${L}]`] = `[//Ring${L}]`;
    table[`[Expl\\Ring${L}]`] = `[\\\\Ring${L}]`;
  }
  return table;
}

const SYMBOL_UPDATE_TABLE = buildUpdateTable();

/**
 * Converts a SELFIES symbol from <v2 to its latest equivalent.
 *
 * @param symbol - an old SELFIES symbol
 * @returns the latest equivalent of the input symbol, or the input symbol
 *   itself if no such equivalent exists
 */
export function modernizeSymbol(symbol: string): string {
  if (Object.prototype.hasOwnProperty.call(SYMBOL_UPDATE_TABLE, symbol)) {
    return SYMBOL_UPDATE_TABLE[symbol];
  }

  if (symbol.endsWith('expl]')) { // e.g. [XXXexpl]
    let bondChar = '';
    let atomSymbol = symbol.slice(1, -5);
    if ('=#/\\'.includes(symbol[1])) {
      bondChar = symbol[1];
      atomSymbol = symbol.slice(2, -5);
    }

    const atom = smilesToAtom(`[${atomSymbol}]`);
    if (atom !== null && !atom.isAromatic) {
      atomSymbol = atomToSmiles(atom, false); // standardize
      symbol = `[${bondChar}${atomSymbol}]`;
    }
  }

  return symbol;
}

/**
 * Warns about deprecated SELFIES usage
 */
export function warnDeprecatedCompatibility(): void {
  console.warn(
    'SELFIES: decoder() may behave differently than in previous major ' +
    'releases. We recommend using SELFIES that are up to date.'
  );
}
