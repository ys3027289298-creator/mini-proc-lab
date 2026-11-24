# Layer: biome

Master seed: `20261124`  |  binary: `worldgen`  |  output: `biome.u8`

## Two runs, same seed

| run | bytes | sha256 |
|---|---:|---|
| A | 65536 | `1e088bb08936bff2791e8126e6e36fe287bfdc86274ddea3bf79d4a25b1e6bfe` |
| B | 65536 | `1e088bb08936bff2791e8126e6e36fe287bfdc86274ddea3bf79d4a25b1e6bfe` |

**Verdict: bit-identical** (`cmp` clean).

## Thread counts

| threads | sha256 | identical to threads=1 |
|---|---|---|
| 1 | `1e088bb08936bff2791e8126e6e36fe287bfdc86274ddea3bf79d4a25b1e6bfe` | yes |
| 4 | `1e088bb08936bff2791e8126e6e36fe287bfdc86274ddea3bf79d4a25b1e6bfe` | yes |
| 8 | `1e088bb08936bff2791e8126e6e36fe287bfdc86274ddea3bf79d4a25b1e6bfe` | yes |

**Verdict: bit-identical at 1, 4 and 8 threads.**
