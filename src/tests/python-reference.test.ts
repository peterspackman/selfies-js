import { test, describe } from 'node:test';
import assert from 'node:assert/strict';
import * as selfies from '../index.js';
import { ENCODER_CASES, ENCODER_ERROR_CASES, DECODER_CASES } from './python-reference-cases.js';

describe('Python reference equivalence', () => {
  test('encoder matches Python selfies', () => {
    selfies.setSemanticConstraints();
    for (const [smiles, expectedSelfies, expectedDecoded] of ENCODER_CASES) {
      const actual = selfies.encoder(smiles);
      assert.equal(actual, expectedSelfies, `encoder(${smiles})`);
      assert.equal(selfies.decoder(actual), expectedDecoded, `decoder(${actual})`);
    }
  });

  test('encoder rejects invalid or unsupported SMILES', () => {
    for (const smiles of ENCODER_ERROR_CASES) {
      assert.throws(() => selfies.encoder(smiles), selfies.EncoderError, `encoder(${JSON.stringify(smiles)})`);
    }
  });

  test('decoder matches Python selfies', () => {
    selfies.setSemanticConstraints();
    for (const [selfiesStr, expectedSmiles] of DECODER_CASES) {
      assert.equal(selfies.decoder(selfiesStr), expectedSmiles, `decoder(${selfiesStr})`);
    }
  });
});
