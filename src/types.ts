/**
 * Type definitions for SELFIES library
 */

// Attribution tracking types
export class Attribution {
  constructor(
    public index: number,
    public token: string
  ) {}
}

export interface TokenAttribution {
  token: string;
  attribution: Attribution[];
}

export type AttributionMap = TokenAttribution[];

// Encoding types
export type EncodingType = "label" | "one_hot" | "both";

// Semantic constraints mapping
export type SemanticConstraints = Record<string, number>;

// Options for encoding/decoding with attribution
export interface AttributionOptions {
  attribute?: boolean;
}

// Options for decoding including backward compatibility
export interface DecoderOptions extends AttributionOptions {
  compatible?: boolean;
}

// Return types for functions with attribution
export type EncoderResult = string | [string, AttributionMap];
export type DecoderResult = string | [string, AttributionMap];

// Vocabulary mappings for encoding utilities
export type StringToIndex = Record<string, number>;
export type IndexToString = Record<number, string>;

// State type for SELFIES decoding - represents bonding capacity remaining
export type State = number | null; // null means no more bonds can be formed