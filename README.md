# @peterspackman/selfies

TypeScript port of the [Python SELFIES library](https://github.com/aspuru-guzik-group/selfies) for robust molecular string representation.

## Installation

```bash
npm install @peterspackman/selfies
```

## Usage

```typescript
import * as selfies from '@peterspackman/selfies';

// SMILES to SELFIES
const benzene = "c1ccccc1";
const benzeneSelfies = selfies.encoder(benzene);  
// "[C][=C][C][=C][C][=C][Ring1][=Branch1]"

// SELFIES to SMILES
const smiles = selfies.decoder(benzeneSelfies);
// "C1=CC=CC=C1"

// Working with SELFIES
const length = selfies.lenSelfies(benzeneSelfies);  // 8
const symbols = [...selfies.splitSelfies(benzeneSelfies)];
// ['[C]', '[=C]', '[C]', '[=C]', '[C]', '[=C]', '[Ring1]', '[=Branch1]']
```

## API

- `encoder(smiles, { strict?, attribute? })` - Convert SMILES to SELFIES. Aromatic SMILES are kekulized first.
  With `strict` (default `true`) molecules violating the current semantic constraints are rejected.
  Throws `EncoderError` for invalid or unsupported SMILES (`*`, `$`, chirality other than `@`/`@@`).
- `decoder(selfies, { attribute?, compatible? })` - Convert SELFIES to SMILES. `compatible` accepts
  symbols from SELFIES < v2. Throws `DecoderError` for malformed SELFIES.
- `getSemanticConstraints()` / `setSemanticConstraints(preset | constraints)` / `getPresetConstraints(name)` -
  Configure bonding constraints (`default`, `octet_rule`, `hypervalent`)
- `getSemanticRobustAlphabet()` - Symbols constrained under the current semantic constraints
- `lenSelfies(selfies)` - Get symbol count
- `splitSelfies(selfies)` - Split into symbols
- `getAlphabetFromSelfies(selfies[])` - Extract unique symbols
- `selfiesToEncoding()` / `encodingToSelfies()` / `batchSelfiesToFlatHot()` / `batchFlatHotToSelfies()` - ML encodings

With `attribute: true`, `encoder` and `decoder` return `[string, attributionMap]`, where each entry
`{ index, token, attribution }` links an output token to the input tokens that produced it.

## Differences from Python selfies

Output matches the Python `selfies` package (v2.2); the test suite includes reference cases generated
with it. Where the Python implementation is wrong, this port fixes it instead:

- **Kekulization**: Python's perfect-matching search ignores odd cycles, so for some fused ring systems it
  returns an invalid matching or fails to kekulize. Its result is used when valid; otherwise this port
  falls back to Edmonds' blossom algorithm.
- **Ring numbers in SMILES output**: Python never reuses ring-closure numbers and writes `%100` and
  above, which is not valid SMILES. Numbering is the same up to 99, after which closed numbers are reused.
- **Multi-digit charges**: Python's encoder writes symbols such as `[Fe+10]` that its decoder rejects;
  here they decode.
- **Errors**: a ring bond from an atom to itself (`C11`) raises `EncoderError` (Python raises
  `IndexError`). Rings or branches too long to express with three index symbols raise `EncoderError`
  (Python silently writes an undecodable symbol such as `[Branch4]`).
- **Attribution indices**: `index` is the position of the output token: the symbol index in SELFIES, the
  starting character in SMILES. Python's indices are off for branch symbols, nested branches, multi-character
  SMILES tokens and fragments after a `.`.

Like Python, atom classes (e.g. `[CH4:1]`) are accepted but dropped.

## License

Apache 2.0 (same as original Python implementation)

## Citation

If you use this in research, please cite the original SELFIES paper:

```bibtex
@article{krenn2020selfies,
  title={Self-referencing embedded strings (SELFIES): A 100% robust molecular string representation},
  author={Krenn, Mario and H{\"a}se, Florian and Nigam, AkshatKumar and Friederich, Pascal and Aspuru-Guzik, Alan},
  journal={Machine Learning: Science and Technology},
  volume={1},
  number={4},
  pages={045024},
  year={2020},
  publisher={IOP Publishing}
}
```
