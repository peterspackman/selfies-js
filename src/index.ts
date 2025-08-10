/**
 * SELFIES: a robust representation of semantically constrained graphs with an
 * example application in chemistry.
 *
 * SELFIES (SELF-referencIng Embedded Strings) is a general-purpose,
 * sequence-based, robust representation of semantically constrained graphs.
 * It is based on a Chomsky type-2 grammar, augmented with two self-referencing
 * functions. A main objective is to use SELFIES as direct input into machine
 * learning models, in particular in generative models, for the generation of
 * outputs with high validity.
 *
 * The code presented here is a concrete application of SELFIES in chemistry, for
 * the robust representation of molecules.
 *
 * Typical usage example:
 *     import * as selfies from 'selfies';
 *
 *     const benzene = "C1=CC=CC=C1";
 *     const benzeneSelfies = selfies.encoder(benzene);
 *     const benzeneSmiles = selfies.decoder(benzeneSelfies);
 */

export const version = "0.2.0";

// Core encoding/decoding functions
export { encoder } from './encoder.js';
export { decoder } from './decoder.js';

// Bond constraints and alphabet functions
export {
  getPresetConstraints,
  getSemanticConstraints,
  getSemanticRobustAlphabet,
  setSemanticConstraints,
  getBondingCapacity
} from './bond-constraints.js';

// Exception classes
export { EncoderError, DecoderError, SMILESParserError } from './exceptions.js';

// Utility functions
export {
  lenSelfies,
  splitSelfies,
  getAlphabetFromSelfies
} from './utils/selfies-utils.js';

export {
  selfiesToEncoding,
  batchSelfiesToFlatHot,
  encodingToSelfies,
  batchFlatHotToSelfies
} from './utils/encoding-utils.js';

// Export types for TypeScript users
export type {
  AttributionMap,
  TokenAttribution,
  Attribution,
  EncodingType,
  SemanticConstraints
} from './types.js';

// Export molecular graph classes for advanced users
export {
  Atom,
  DirectedBond,
  MolecularGraph
} from './mol-graph.js';