import { test, describe } from 'node:test';
import assert from 'node:assert/strict';
import * as selfies from '../index.js';

// Places where Python selfies 2.2 is wrong and this port deliberately differs.
// See "Differences from Python selfies" in the README.
describe('Fixes relative to Python selfies', () => {
  test('ring numbers stay valid past 99 rings', () => {
    const smiles = Array(105).fill('C1CC1').join('.');
    const decoded = selfies.decoder(selfies.encoder(smiles));
    // Python writes %100..%105 here, which is not valid SMILES
    assert(!/%\d{3}|%\(/.test(decoded), decoded);
    assert(decoded.startsWith('C1CC1.C2CC2.'));
    assert(decoded.endsWith('.C1CC1'));
  });

  test('multi-digit charges round-trip', () => {
    // Python's encoder writes [Fe+10] but its decoder rejects it
    assert.equal(selfies.encoder('[Fe+10]'), '[Fe+10]');
    assert.equal(selfies.decoder('[Fe+10]'), '[Fe+10]');
  });

  test('ring bond from an atom to itself is an EncoderError', () => {
    // Python raises an IndexError
    assert.throws(() => selfies.encoder('C11'), selfies.EncoderError);
  });

  test('rings or branches too long for SELFIES are rejected', () => {
    // Python silently writes an undecodable [Branch4] symbol
    const smiles = 'C(' + 'C'.repeat(5000) + ')C';
    assert.throws(() => selfies.encoder(smiles), selfies.EncoderError);
  });

  test('decoder attribution indices account for dots', () => {
    const [smiles, attribution] = selfies.decoder('[C][O].[N][C]', { attribute: true });
    assert.equal(smiles, 'CO.NC');
    for (const entry of attribution) {
      assert.equal(smiles[entry.index], entry.token);
    }
  });

  test('encoder attribution indices are correct in nested branches', () => {
    const [selfiesStr, attribution] = selfies.encoder('CC(C(F)Cl)Br', { attribute: true });
    const symbols = [...selfies.splitSelfies(selfiesStr)];
    for (const entry of attribution) {
      assert.equal(symbols[entry.index], entry.token, JSON.stringify(entry));
    }
  });
});
