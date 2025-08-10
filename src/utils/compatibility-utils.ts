/**
 * Backward compatibility utilities for legacy SELFIES symbols
 */

// Branch and Ring symbol mappings for backward compatibility
const BRANCH_RING_UPDATES: Record<string, string> = {
  // Branch symbols (levels 1, 2, 3)
  '[Branch1_1]': '[Branch1]',
  '[Branch1_2]': '[=Branch1]', 
  '[Branch1_3]': '[#Branch1]',
  '[Branch2_1]': '[Branch2]',
  '[Branch2_2]': '[=Branch2]',
  '[Branch2_3]': '[#Branch2]',
  '[Branch3_1]': '[Branch3]',
  '[Branch3_2]': '[=Branch3]',
  '[Branch3_3]': '[#Branch3]',
  
  // Ring symbols with explicit bond orders
  '[Expl=Ring1]': '[=Ring1]',
  '[Expl#Ring1]': '[#Ring1]',
  '[Expl=Ring2]': '[=Ring2]',
  '[Expl#Ring2]': '[#Ring2]',
  '[Expl=Ring3]': '[=Ring3]',
  '[Expl#Ring3]': '[#Ring3]',
  '[Expl=Ring4]': '[=Ring4]',
  '[Expl#Ring4]': '[#Ring4]',
  '[Expl=Ring5]': '[=Ring5]',
  '[Expl#Ring5]': '[#Ring5]',
  '[Expl=Ring6]': '[=Ring6]',
  '[Expl#Ring6]': '[#Ring6]',
  '[Expl=Ring7]': '[=Ring7]',
  '[Expl#Ring7]': '[#Ring7]',
  '[Expl=Ring8]': '[=Ring8]',
  '[Expl#Ring8]': '[#Ring8]',
  '[Expl=Ring9]': '[=Ring9]',
  '[Expl#Ring9]': '[#Ring9]',
};

/**
 * Modernizes legacy SELFIES symbols to current format
 */
export function modernizeSymbol(symbol: string): string {
  // Check branch/ring updates first
  if (BRANCH_RING_UPDATES[symbol]) {
    return BRANCH_RING_UPDATES[symbol];
  }
  
  // Handle atom symbols with 'expl' suffix
  if (symbol.endsWith('expl]')) {
    // Extract the part before 'expl]' and after '['
    const atomPart = symbol.slice(1, -5); // Remove '[' and 'expl]'
    return `[${atomPart}1]`;
  }
  
  // Return symbol unchanged if no modernization needed
  return symbol;
}

/**
 * Warns about deprecated SELFIES usage
 */
export function warnDeprecatedCompatibility(): void {
  console.warn(
    'SELFIES: Using backward compatibility mode is deprecated. ' +
    'Please convert your SELFIES strings to the current format. ' +
    'Behavior may differ from previous major releases.'
  );
}