# Layer: hydro

Master seed: `20261124`  |  binary: `worldgen`  |  output: `hydro.f32`

## Two runs, same seed

| run | bytes | sha256 |
|---|---:|---|
| A | 262144 | `582d4e6e4d686418d7283b286079fbe4d81aac83c34ace329a10900989090c62` |
| B | 262144 | `582d4e6e4d686418d7283b286079fbe4d81aac83c34ace329a10900989090c62` |

**Verdict: bit-identical** (`cmp` clean).

## Thread counts

| threads | sha256 | identical to threads=1 |
|---|---|---|
| 1 | `582d4e6e4d686418d7283b286079fbe4d81aac83c34ace329a10900989090c62` | yes |
| 4 | `582d4e6e4d686418d7283b286079fbe4d81aac83c34ace329a10900989090c62` | yes |
| 8 | `582d4e6e4d686418d7283b286079fbe4d81aac83c34ace329a10900989090c62` | yes |

**Verdict: bit-identical at 1, 4 and 8 threads.**
