# Layer: noise

Master seed: `20261124`  |  binary: `worldgen`  |  output: `noise.f32`

## Two runs, same seed

| run | bytes | sha256 |
|---|---:|---|
| A | 262144 | `daec11f6a9aa4465c05fa26c40b2be60fb33737c635b936616f1ec7461fd457a` |
| B | 262144 | `daec11f6a9aa4465c05fa26c40b2be60fb33737c635b936616f1ec7461fd457a` |

**Verdict: bit-identical** (`cmp` clean).

## Thread counts

| threads | sha256 | identical to threads=1 |
|---|---|---|
| 1 | `daec11f6a9aa4465c05fa26c40b2be60fb33737c635b936616f1ec7461fd457a` | yes |
| 4 | `daec11f6a9aa4465c05fa26c40b2be60fb33737c635b936616f1ec7461fd457a` | yes |
| 8 | `daec11f6a9aa4465c05fa26c40b2be60fb33737c635b936616f1ec7461fd457a` | yes |

**Verdict: bit-identical at 1, 4 and 8 threads.**
