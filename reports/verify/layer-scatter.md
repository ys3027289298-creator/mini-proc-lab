# Layer: scatter

Master seed: `20261124`  |  binary: `worldgen`  |  output: `scatter.f32`

## Two runs, same seed

| run | bytes | sha256 |
|---|---:|---|
| A | 40652 | `34e21244b33b64281bca8e8a6a5a4e015ac7ce49b3899293d1d38a81fb3e821e` |
| B | 40652 | `34e21244b33b64281bca8e8a6a5a4e015ac7ce49b3899293d1d38a81fb3e821e` |

**Verdict: bit-identical** (`cmp` clean).

## Thread counts

| threads | sha256 | identical to threads=1 |
|---|---|---|
| 1 | `34e21244b33b64281bca8e8a6a5a4e015ac7ce49b3899293d1d38a81fb3e821e` | yes |
| 4 | `34e21244b33b64281bca8e8a6a5a4e015ac7ce49b3899293d1d38a81fb3e821e` | yes |
| 8 | `34e21244b33b64281bca8e8a6a5a4e015ac7ce49b3899293d1d38a81fb3e821e` | yes |

**Verdict: bit-identical at 1, 4 and 8 threads.**
