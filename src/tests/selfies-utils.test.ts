import { test, describe } from 'node:test';
import assert from 'node:assert/strict';
import * as selfies from '../index.js';

interface TestEntry {
  selfies: string;
  symbols: string[];
  label: number[];
  oneHot: number[][];
}

// Dataset fixture equivalent
const DATASET: TestEntry[] = [
  {
    selfies: "",
    symbols: [],
    label: [0, 0, 0, 0],
    oneHot: [
      [1, 0, 0, 0, 0],
      [1, 0, 0, 0, 0],
      [1, 0, 0, 0, 0],
      [1, 0, 0, 0, 0]
    ]
  },
  {
    selfies: "[C][C][C]",
    symbols: ["[C]", "[C]", "[C]"],
    label: [3, 3, 3, 0],
    oneHot: [
      [0, 0, 0, 1, 0],
      [0, 0, 0, 1, 0],
      [0, 0, 0, 1, 0],
      [1, 0, 0, 0, 0]
    ]
  },
  {
    selfies: "[C].[C]",
    symbols: ["[C]", ".", "[C]"],
    label: [3, 2, 3, 0],
    oneHot: [
      [0, 0, 0, 1, 0],
      [0, 0, 1, 0, 0],
      [0, 0, 0, 1, 0],
      [1, 0, 0, 0, 0]
    ]
  },
  {
    selfies: "[C][O][C][F]",
    symbols: ["[C]", "[O]", "[C]", "[F]"],
    label: [3, 1, 3, 4],
    oneHot: [
      [0, 0, 0, 1, 0],
      [0, 1, 0, 0, 0],
      [0, 0, 0, 1, 0],
      [0, 0, 0, 0, 1]
    ]
  },
  {
    selfies: "[C][O][C]",
    symbols: ["[C]", "[O]", "[C]"],
    label: [3, 1, 3, 0],
    oneHot: [
      [0, 0, 0, 1, 0],
      [0, 1, 0, 0, 0],
      [0, 0, 0, 1, 0],
      [1, 0, 0, 0, 0]
    ]
  }
];

const VOCAB_STOI = { "[nop]": 0, "[O]": 1, ".": 2, "[C]": 3, "[F]": 4 };
const VOCAB_ITOS = Object.fromEntries(
  Object.entries(VOCAB_STOI).map(([k, v]) => [v, k])
);
const PAD_TO_LEN = 4;

// Flat hot encodings for batch testing
const DATASET_FLAT_HOTS: number[][] = DATASET.map(entry => 
  entry.oneHot.flat()
);

describe('SELFIES Utils Tests', () => {
  test('len selfies', async () => {
    for (const entry of DATASET) {
      assert.equal(selfies.lenSelfies(entry.selfies), entry.symbols.length);
    }
  });

  test('split selfies', async () => {
    for (const entry of DATASET) {
      const result = Array.from(selfies.splitSelfies(entry.selfies));
      assert.deepEqual(result, entry.symbols);
    }
  });

  test('get alphabet from selfies', async () => {
    const selfiesStrings = DATASET.map(entry => entry.selfies);
    const alphabet = selfies.getAlphabetFromSelfies(selfiesStrings);
    alphabet.add("[nop]");
    alphabet.add(".");

    const expectedAlphabet = new Set(Object.keys(VOCAB_STOI));
    assert.deepEqual(alphabet, expectedAlphabet);
  });

  test('selfies to encoding', async () => {
    for (const entry of DATASET) {
      const [label, oneHot] = selfies.selfiesToEncoding(
        entry.selfies, 
        VOCAB_STOI, 
        PAD_TO_LEN, 
        "both"
      );

      assert.deepEqual(label, entry.label);
      assert.deepEqual(oneHot, entry.oneHot);

      // Recover original selfies from label
      let recoveredFromLabel = selfies.encodingToSelfies(label, VOCAB_ITOS, "label");
      recoveredFromLabel = recoveredFromLabel.replace(/\[nop\]/g, "");
      assert.equal(recoveredFromLabel, entry.selfies);

      // Recover original selfies from one-hot
      let recoveredFromOneHot = selfies.encodingToSelfies(oneHot, VOCAB_ITOS, "one_hot");
      recoveredFromOneHot = recoveredFromOneHot.replace(/\[nop\]/g, "");
      assert.equal(recoveredFromOneHot, entry.selfies);
    }
  });

  test('selfies to flat hot', async () => {
    const batch = DATASET.map(entry => entry.selfies);
    const flatHots = selfies.batchSelfiesToFlatHot(batch, VOCAB_STOI, PAD_TO_LEN);

    assert.deepEqual(flatHots, DATASET_FLAT_HOTS);

    // Recover original selfies
    const recovered = selfies.batchFlatHotToSelfies(flatHots, VOCAB_ITOS);
    const expected = batch.map(s => s.replace(/\[nop\]/g, ""));
    assert.deepEqual(expected, recovered.map(s => s.replace(/\[nop\]/g, "")));
  });
});