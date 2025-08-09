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

- `encoder(smiles)` - Convert SMILES to SELFIES
- `decoder(selfies)` - Convert SELFIES to SMILES  
- `setSemanticConstraints()` - Configure bonding constraints
- `lenSelfies(selfies)` - Get symbol count
- `splitSelfies(selfies)` - Split into symbols
- `getAlphabetFromSelfies(selfies[])` - Extract unique symbols
- `selfiesToEncoding()` / `encodingToSelfies()` - ML encodings

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
