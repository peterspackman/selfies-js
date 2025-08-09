/**
 * Encoding utility functions for SELFIES strings
 */

import { lenSelfies, splitSelfies } from './selfies-utils.js';
import type { StringToIndex, IndexToString, EncodingType } from '../types.js';

/**
 * Converts a SELFIES string into its label (integer) and/or one-hot encoding.
 * 
 * @param selfies - The SELFIES string to be encoded
 * @param vocabStoi - Dictionary mapping SELFIES symbols to indices
 * @param padToLen - Length to pad the SELFIES string to. Defaults to -1 (no padding)
 * @param encType - Type of encoding: 'label', 'one_hot', or 'both'
 * @returns The encoded SELFIES string in the requested format
 * 
 * @example
 * ```typescript
 * import { selfiesToEncoding } from 'selfies';
 * const result = selfiesToEncoding("[C][F]", {"[C]": 0, "[F]": 1});
 * console.log(result); // [[0, 1], [[1, 0], [0, 1]]]
 * ```
 */
export function selfiesToEncoding(
  selfies: string,
  vocabStoi: StringToIndex,
  padToLen: number = -1,
  encType: EncodingType = 'both'
): number[] | number[][] | [number[], number[][]] {
  
  if (!['label', 'one_hot', 'both'].includes(encType)) {
    throw new Error("encType must be one of ('label', 'one_hot', 'both')");
  }

  // Pad with [nop]
  let paddedSelfies = selfies;
  if (padToLen > lenSelfies(selfies)) {
    const padCount = padToLen - lenSelfies(selfies);
    paddedSelfies += '[nop]'.repeat(padCount);
  }

  // Integer encode
  const integerEncoded: number[] = [];
  for (const char of splitSelfies(paddedSelfies)) {
    if (char === '.' && !('.' in vocabStoi)) {
      throw new Error(
        "The SELFIES string contains two unconnected molecules " +
        "(given by the '.' character), but vocabStoi does not " +
        "contain the '.' key. Please add it to the vocabulary " +
        "or separate the molecules."
      );
    }
    
    if (!(char in vocabStoi)) {
      throw new Error(`Symbol '${char}' not found in vocabulary`);
    }
    
    integerEncoded.push(vocabStoi[char]);
  }

  if (encType === 'label') {
    return integerEncoded;
  }

  // One-hot encode
  const vocabSize = Object.keys(vocabStoi).length;
  const oneHotEncoded: number[][] = [];
  for (const index of integerEncoded) {
    const letter = new Array(vocabSize).fill(0);
    letter[index] = 1;
    oneHotEncoded.push(letter);
  }

  if (encType === 'one_hot') {
    return oneHotEncoded;
  }

  return [integerEncoded, oneHotEncoded];
}

/**
 * Converts a label (integer) or one-hot encoding into a SELFIES string.
 * 
 * @param encoding - A label or one-hot encoding
 * @param vocabItos - Dictionary mapping indices to SELFIES symbols
 * @param encType - Type of encoding: 'label' or 'one_hot'
 * @returns The SELFIES string represented by the input encoding
 * 
 * @example
 * ```typescript
 * import { encodingToSelfies } from 'selfies';
 * const oneHot = [[0, 1, 0], [0, 0, 1], [1, 0, 0]];
 * const vocabItos = {0: "[nop]", 1: "[C]", 2: "[F]"};
 * const result = encodingToSelfies(oneHot, vocabItos, "one_hot");
 * console.log(result); // '[C][F][nop]'
 * ```
 */
export function encodingToSelfies(
  encoding: number[] | number[][],
  vocabItos: IndexToString,
  encType: 'label' | 'one_hot'
): string {
  
  if (!['label', 'one_hot'].includes(encType)) {
    throw new Error("encType must be one of ('label', 'one_hot')");
  }

  let integerEncoded: number[];
  
  if (encType === 'one_hot') {
    // Get integer encoding from one-hot
    integerEncoded = [];
    const oneHot = encoding as number[][];
    for (const row of oneHot) {
      const index = row.indexOf(1);
      if (index === -1) {
        throw new Error("Invalid one-hot encoding: no 1 found in row");
      }
      integerEncoded.push(index);
    }
  } else {
    integerEncoded = encoding as number[];
  }

  // Integer encoding -> SELFIES
  const charList = integerEncoded.map(i => {
    if (!(i in vocabItos)) {
      throw new Error(`Index ${i} not found in vocabulary`);
    }
    return vocabItos[i];
  });
  
  return charList.join('');
}

/**
 * Converts a list of SELFIES strings into its list of flattened one-hot encodings.
 * 
 * @param selfiesBatch - List of SELFIES strings to be encoded
 * @param vocabStoi - Dictionary mapping SELFIES symbols to indices
 * @param padToLen - Length to pad each SELFIES string to. Defaults to -1 (no padding)
 * @returns List of flattened one-hot encodings
 * 
 * @example
 * ```typescript
 * import { batchSelfiesToFlatHot } from 'selfies';
 * const batch = ["[C]", "[C][C]"];
 * const vocabStoi = {"[nop]": 0, "[C]": 1};
 * const result = batchSelfiesToFlatHot(batch, vocabStoi, 2);
 * console.log(result); // [[0, 1, 1, 0], [0, 1, 0, 1]]
 * ```
 */
export function batchSelfiesToFlatHot(
  selfiesBatch: string[],
  vocabStoi: StringToIndex,
  padToLen: number = -1
): number[][] {
  
  const hotList: number[][] = [];
  
  for (const selfies of selfiesBatch) {
    const oneHot = selfiesToEncoding(selfies, vocabStoi, padToLen, 'one_hot') as number[][];
    const flattened = oneHot.flat();
    hotList.push(flattened);
  }

  return hotList;
}

/**
 * Converts a list of flattened one-hot encodings into a list of SELFIES strings.
 * 
 * @param oneHotBatch - List of flattened one-hot encodings
 * @param vocabItos - Dictionary mapping indices to SELFIES symbols
 * @returns List of SELFIES strings represented by the input encodings
 * 
 * @example
 * ```typescript
 * import { batchFlatHotToSelfies } from 'selfies';
 * const batch = [[0, 1, 1, 0], [0, 1, 0, 1]];
 * const vocabItos = {0: "[nop]", 1: "[C]"};
 * const result = batchFlatHotToSelfies(batch, vocabItos);
 * console.log(result); // ['[C][nop]', '[C][C]']
 * ```
 */
export function batchFlatHotToSelfies(
  oneHotBatch: number[][],
  vocabItos: IndexToString
): string[] {
  
  const selfiesList: string[] = [];
  const vocabSize = Object.keys(vocabItos).length;

  for (const flatOneHot of oneHotBatch) {
    // Reshape to an L x M array where each column represents an alphabet
    // entry and each row is a position in the selfies
    if (flatOneHot.length % vocabSize !== 0) {
      throw new Error(
        "size of vector in oneHotBatch not divisible by the length of the vocabulary."
      );
    }

    const L = flatOneHot.length / vocabSize;
    const oneHot: number[][] = [];

    for (let i = 0; i < L; i++) {
      const start = vocabSize * i;
      const end = vocabSize * (i + 1);
      oneHot.push(flatOneHot.slice(start, end));
    }

    const selfies = encodingToSelfies(oneHot, vocabItos, 'one_hot');
    selfiesList.push(selfies);
  }

  return selfiesList;
}