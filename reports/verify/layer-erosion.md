# Layer: erosion

Master seed: `20261124`  |  binary: `worldgen`  |  output: `erosion.f32`

## Two runs, same seed

| run | bytes | sha256 |
|---|---:|---|
| A | 262144 | `2ef80e337d568d97f92c39381e328290733102d90a6e7a4c71d1f42d0490abc7` |
| B | 262144 | `2ef80e337d568d97f92c39381e328290733102d90a6e7a4c71d1f42d0490abc7` |

**Verdict: bit-identical** (`cmp` clean).

## Thread counts

| threads | sha256 | identical to threads=1 |
|---|---|---|
| 1 | `2ef80e337d568d97f92c39381e328290733102d90a6e7a4c71d1f42d0490abc7` | yes |
| 4 | `2ef80e337d568d97f92c39381e328290733102d90a6e7a4c71d1f42d0490abc7` | yes |
| 8 | `2ef80e337d568d97f92c39381e328290733102d90a6e7a4c71d1f42d0490abc7` | yes |

**Verdict: bit-identical at 1, 4 and 8 threads.**
