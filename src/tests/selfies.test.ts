import { test, describe } from 'node:test';
import assert from 'node:assert/strict';
import * as selfies from '../index.js';

// Test configuration - these would come from command line args in Python
const TRIALS = parseInt(process.env.TRIALS || '100'); // Use 100 trials like Python test
const MAX_SELFIES_LEN = 1000;

// Large alphabet fixture equivalent
const LARGE_ALPHABET = [
  ...selfies.getSemanticRobustAlphabet(),
  "[#Br]", "[#Branch1]", "[#Branch2]", "[#Branch3]", "[#C@@H1]",
  "[#C@@]", "[#C@H1]", "[#C@]", "[#C]", "[#Cl]", "[#F]", "[#H]", "[#I]",
  "[#NH1]", "[#N]", "[#O]", "[#P]", "[#Ring1]", "[#Ring2]", "[#Ring3]",
  "[#S]", "[/Br]", "[/C@@H1]", "[/C@@]", "[/C@H1]", "[/C@]", "[/C]",
  "[/Cl]", "[/F]", "[/H]", "[/I]", "[/NH1]", "[/N]", "[/O]", "[/P]",
  "[/S]", "[=Br]", "[=Branch1]", "[=Branch2]", "[=Branch3]", "[=C@@H1]",
  "[=C@@]", "[=C@H1]", "[=C@]", "[=C]", "[=Cl]", "[=F]", "[=H]", "[=I]",
  "[=NH1]", "[=N]", "[=O]", "[=P]", "[=Ring1]", "[=Ring2]", "[=Ring3]",
  "[=S]", "[Br]", "[Branch1]", "[Branch2]", "[Branch3]", "[C@@H1]",
  "[C@@]", "[C@H1]", "[C@]", "[C]", "[Cl]", "[F]", "[H]", "[I]", "[NH1]",
  "[N]", "[O]", "[P]", "[Ring1]", "[Ring2]", "[Ring3]", "[S]", "[\\Br]",
  "[\\C@@H1]", "[\\C@@]", "[\\C@H1]", "[\\C@]", "[\\C]", "[\\Cl]",
  "[\\F]", "[\\H]", "[\\I]", "[\\NH1]", "[\\N]", "[\\O]", "[\\P]",
  "[\\S]", "[nop]"
];

// Helper function equivalent to Python's random.choices
function randomChoices<T>(population: T[], k: number): T[] {
  const result: T[] = [];
  for (let i = 0; i < k; i++) {
    const randomIndex = Math.floor(Math.random() * population.length);
    result.push(population[randomIndex]);
  }
  return result;
}

describe('SELFIES Core Tests', () => {
  test('random selfies decoder', async () => {
    // Test random SELFIES generation and decoding like Python does
    // Uses the semantic robust alphabet (same as Python test)
    const alphabet = [...selfies.getSemanticRobustAlphabet()];
    
    for (let i = 0; i < TRIALS; i++) {
      // Create random SELFIES and decode (use shorter lengths like Python)
      const randLen = Math.floor(Math.random() * 50) + 5; // 5-54 symbols like Python test
      const randSelfies = randomChoices(alphabet, randLen).join('');
      
      let smiles: string;
      try {
        smiles = selfies.decoder(randSelfies) as string;
      } catch (e) {
        assert.fail(`Decode failed for SELFIES: ${randSelfies}\nError: ${e instanceof Error ? e.message : String(e)}`);
        return;
      }
      
      // Additional validation: try to re-encode the SMILES (like Python test does)
      try {
        selfies.encoder(smiles);
        // If we get here, the SMILES was valid and could be re-encoded
      } catch (e) {
        assert.fail(`Decoded SMILES is invalid: ${smiles}\nOriginal SELFIES: ${randSelfies}\nRe-encode error: ${e instanceof Error ? e.message : String(e)}`);
      }
    }
  });

  test('nop symbol decoder', async () => {
    // Tests that the '[nop]' symbol is always skipped over.
    
    const alphabet = LARGE_ALPHABET.filter(symbol => symbol !== '[nop]');
    
    for (let i = 0; i < 100; i++) {
      // Create random SELFIES with and without [nop]
      const randLen = Math.floor(Math.random() * MAX_SELFIES_LEN) + 1;
      const randMol = randomChoices(alphabet, randLen);
      const nopCount = MAX_SELFIES_LEN - randLen;
      randMol.push(...Array(nopCount).fill('[nop]'));
      
      // Shuffle the array
      for (let j = randMol.length - 1; j > 0; j--) {
        const k = Math.floor(Math.random() * (j + 1));
        [randMol[j], randMol[k]] = [randMol[k], randMol[j]];
      }
      
      const withNops = randMol.join('');
      const withoutNops = withNops.replace(/\[nop\]/g, '');
      
      assert.equal(
        selfies.decoder(withNops) as string, 
        selfies.decoder(withoutNops) as string
      );
    }
  });

  test('get semantic constraints', async () => {
    const constraints = selfies.getSemanticConstraints();
    assert(constraints !== selfies.getSemanticConstraints()); // not alias
    assert('?' in constraints);
  });

  test('change constraints cache clear', async () => {
    const alphabet = selfies.getSemanticRobustAlphabet();
    assert.deepEqual(alphabet, selfies.getSemanticRobustAlphabet());
    assert.equal(selfies.decoder("[C][#C]"), "C#C");

    const newConstraints = selfies.getSemanticConstraints();
    newConstraints["C"] = 1;
    newConstraints["O"] = 1;
    selfies.setSemanticConstraints(newConstraints);

    // [=C] and [#C] remain, since they are index symbols; [=O] is dropped
    const newAlphabet = selfies.getSemanticRobustAlphabet();
    assert.notDeepEqual(newAlphabet, alphabet);
    assert(newAlphabet.has("[=C]"));
    assert(!newAlphabet.has("[=O]"));
    assert.equal(selfies.decoder("[C][#C]"), "CC");

    selfies.setSemanticConstraints(); // re-set alphabet
  });

  test('invalid or unsupported smiles encoder', async () => {
    const malformedSmiles = [
      "",
      "(",
      "C(Cl)(Cl)CC[13C",
      "C(CCCOC",
      "C=(CCOC",
      "CCCC)",
      "C1CCCCC",
      "C(F)(F)(F)(F)(F)F",  // violates bond constraints
      "C=C1=CCCCCC1",  // violates bond constraints
      "CC*CC",  // uses wildcard
      "C$C",  // uses $ bond
      "S[As@TB1](F)(Cl)(Br)N",  // unrecognized chirality,
      "SOMETHINGWRONGHERE",
      "1243124124",
    ];

    for (const smiles of malformedSmiles) {
      assert.throws(() => {
        selfies.encoder(smiles);
      }, selfies.EncoderError);
    }
  });

  test('malformed selfies decoder', async () => {
    assert.throws(() => {
      selfies.decoder("[O][=C][O][C][C][C][C][O][N][Branch2_3");
    }, selfies.DecoderError);
  });

  test('decoder attribution', async () => {
    // Ensure constraints are reset to default before this test
    selfies.setSemanticConstraints();
    
    const result = selfies.decoder(
      "[C][N][C][Branch1][C][P][C][C][Ring1][=Branch1]", 
      { attribute: true }
    ) as [string, selfies.AttributionMap];
    
    const [smiles, attributionMap] = result;
    assert.equal(smiles, "C1NC(P)CC1");

    // Check that P lined up
    const p = attributionMap.find(a => a.token === 'P');
    assert(p, 'Failed to find P in attribution map');
    assert.equal(smiles[p.index], 'P');
    assert.deepEqual(p.attribution.map(a => [a.index, a.token]), [[3, '[Branch1]'], [5, '[P]']]);
  });

  test('encoder attribution', async () => {
    // Ensure constraints are reset to default before this test  
    selfies.setSemanticConstraints();
    
    const smiles = "C1([O-])C=CC=C1Cl";
    const expectedIndices = [0, 3, 3, 3, 5, 7, 8, 10, null, null, 12];
    
    const result = selfies.encoder(smiles, { attribute: true }) as [string, selfies.AttributionMap];
    const [, attributionMap] = result;
    
    assert.equal(attributionMap.length, expectedIndices.length);
    for (let i = 0; i < attributionMap.length; i++) {
      const tokenAttr = attributionMap[i];
      if (expectedIndices[i] === null) {
        assert.equal(tokenAttr.attribution.length, 0);
      } else {
        assert.equal(
          expectedIndices[i], 
          tokenAttr.attribution[0].index,
          `found ${tokenAttr.attribution[0].index}; should be ${expectedIndices[i]}`
        );
      }
      if (tokenAttr.token === '[Cl]') {
        const hasClInAttribution = tokenAttr.attribution.some(
          (attr: selfies.Attribution) => attr.token === 'Cl'
        );
        assert(hasClInAttribution, 'Failed to find Cl in attribution map');
      }
    }
  });

  test('aromatic benzene encoding/decoding', async () => {
    const benzeneSelfies = selfies.encoder("c1ccccc1");
    assert.equal(benzeneSelfies, "[C][=C][C][=C][C][=C][Ring1][=Branch1]");
    assert.equal(selfies.decoder(benzeneSelfies), "C1=CC=CC=C1");
  });

  test('simple aromatic molecules', async () => {
    const testCases: [string, string][] = [
      ["c1ccoc1", "[C][C][=C][O][C][=Ring1][Branch1]"],
      ["c1ccsc1", "[C][C][=C][S][C][=Ring1][Branch1]"],
      ["c1ccncc1", "[C][=C][C][=N][C][=C][Ring1][=Branch1]"],
      ["c1cc[nH]c1", "[C][C][=C][NH1][C][=Ring1][Branch1]"],
      ["[O-][n+]1ccccc1", "[O-1][N+1][=C][C][=C][C][=C][Ring1][=Branch1]"],
    ];
    for (const [smiles, expected] of testCases) {
      assert.equal(selfies.encoder(smiles), expected, smiles);
    }
  });

  test('ring and stereo features', async () => {
    // Two-digit ring numbers, ring bond orders and chirality next to rings
    assert.equal(selfies.encoder("C%11CCCCC%11"), "[C][C][C][C][C][C][Ring1][=Branch1]");
    assert.equal(selfies.encoder("C=1CC1"), "[C][C][C][=Ring1][Ring1]");
    assert.equal(selfies.decoder("[C][C][C][=Ring1][Ring1]"), "C=1CC=1");
    assert.equal(selfies.encoder("[C@H]1(F)CCC1"), "[C@H1][Branch1][C][F][C][C][C][Ring1][Branch1]");
    assert.equal(selfies.encoder("ClC(Br)c1ccccc1"), "[Cl][C][Branch1][C][Br][C][=C][C][=C][C][=C][Ring1][=Branch1]");
    assert.equal(selfies.decoder("[C][/C][=C][\\C]"), "C/C=C\\C");
    assert.equal(selfies.decoder("[C][C][C][C][\\/Ring1][Ring2]"), "C\\1CCC/1");
  });

  test('decoder stops once bonding capacity is exhausted', async () => {
    assert.equal(selfies.decoder("[C][F][C]"), "CF");
    assert.equal(selfies.decoder("[F][F][C]"), "FF");
    assert.equal(selfies.decoder(""), "");
    assert.equal(selfies.decoder("[CH0]"), "[CH0]");
  });

  test('encoder strict mode', async () => {
    assert.throws(() => selfies.encoder("C(F)(F)(F)(F)F"), selfies.EncoderError);
    assert.equal(selfies.encoder("C(F)(F)(F)(F)F", { strict: false }), "[C][Branch1][C][F][Branch1][C][F][Branch1][C][F][Branch1][C][F][F]");
  });

  test('backward compatible decoding', async () => {
    const warn = console.warn;
    console.warn = () => {};
    try {
      assert.equal(selfies.decoder("[C][Branch1_2][C][=O][O]", { compatible: true }), "C(=O)O");
      assert.equal(selfies.decoder("[NHexpl][C]", { compatible: true }), "[NH1]C");
      assert.equal(selfies.decoder("[C][C][C][Expl/Ring1][C]", { compatible: true }), "CC=C");
    } finally {
      console.warn = warn;
    }
  });

  test('decoder python equivalence', async () => {
    // Test cases generated from Python SELFIES to ensure exact SMILES output match
    const pythonTestCases = [
      {
        selfies: '[=Branch3][#O+1][=C][=S+1][Branch3][=B+1][Ring1][O][=P][#C][N+1][#P][N]',
        expectedSmiles: '[O+1]=C=[S+1]P#C[N+1]#PN'
      },
      {
        selfies: '[=N+1][#Branch1][=S-1][#O+1][#Branch2][#B][=B-1][Branch1][O-1][N+1][O][P-1][=O]',
        expectedSmiles: '[N+1](#[O+1])B=[B-1]([N+1])O[P-1]=O'
      },
      {
        selfies: '[Ring1][=B][P+1][#Branch1][=B-1][Branch3][=Branch2][#O+1][C][Br][=Ring2]',
        expectedSmiles: 'B[P+1]Br'
      },
      {
        selfies: '[#P+1][Br][#Branch1][=S+1][=S-1][=Branch1][#C+1][#S+1][#N][=Ring3][=Ring2][=Branch2][=Ring2][S+1][S][=Ring1][=Branch3][#C][#C][#B-1][=N-1][=C+1][B+1][#Branch1][Ring3][=B][=S][B-1][=N-1][=P-1]',
        expectedSmiles: '[P+1]Br'
      },
      {
        selfies: '[I][P+1][P][S+1][#Branch3][B][=Branch2][O][C+1][#S][=O+1][#P-1][=O][=P][#C+1][=N-1][Branch1][=O][=C+1][#C+1][=P+1][O-1][=O]',
        expectedSmiles: 'I[P+1]P[S+1][C+1]=S=[O+1][P-1]=O'
      },
      {
        selfies: '[=N][=P+1][#B][#Branch2][=O][#P+1][#P][Branch3][Branch3][Branch3][=B+1][Ring3][=O+1][P+1][#S-1][Branch1][=S+1][#P-1][O+1][#B-1][=N][=C+1][=Ring1][=Ring1][B+1][Cl][#S+1]',
        expectedSmiles: 'N=[P+1]=BO[P+1]#P[S+1]#[P-1][O+1]=[B-1]=N[C+1]'
      },
      {
        selfies: '[=N-1][N+1][C+1][S+1][#C][Branch3][=Ring3][C][=S-1][=B][N][=C-1][#N+1]',
        expectedSmiles: '[N-1][N+1][C+1][S+1]#C'
      },
      {
        selfies: '[Branch1][=B-1][#B][P-1][#B][S][O+1][#N+1][=N-1][=S+1][#N+1]',
        expectedSmiles: '[B-1]#B'
      },
      {
        selfies: '[B][=Ring1][=S-1][P-1][I][N][#B][B+1][F][=C][Branch2][#O+1][#B][#S+1][=Branch1][I][=Ring2][=Branch1][H][#O+1][#C+1][P-1][C-1][#C][#S][Branch3][=B+1]',
        expectedSmiles: 'B[P-1]I'
      },
      {
        selfies: '[P-1][=Ring2][O][S-1][=Ring1][P-1][=Branch3][=S-1][I][F][#P][O-1][=Branch3][=C-1][Ring2][O+1][#P-1][B+1][Branch2][Branch3]',
        expectedSmiles: '[P-1](P)[O-1]'
      }
    ];

    for (let i = 0; i < pythonTestCases.length; i++) {
      const testCase = pythonTestCases[i];
      const actualSmiles = selfies.decoder(testCase.selfies) as string;
      
      assert.equal(
        actualSmiles,
        testCase.expectedSmiles,
        `Test ${i + 1} failed:\n` +
        `  SELFIES: ${testCase.selfies}\n` +
        `  Expected: ${testCase.expectedSmiles}\n` +
        `  Actual:   ${actualSmiles}`
      );
    }
  });
});