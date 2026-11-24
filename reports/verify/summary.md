# worldgen verification report

- master seed: `20261124`
- binary: `/home/ys/feedback/uidreczz28dMgsr9Z9L/worldgen/worldgen`
- result: **PASS**

## Per-layer repeat (same seed, two runs, threads 4)

| layer | identical | sha256 |
|---|---|---|
| noise | yes | `daec11f6a9aa4465c05fa26c40b2be60fb33737c635b936616f1ec7461fd457a` |
| erosion | yes | `2ef80e337d568d97f92c39381e328290733102d90a6e7a4c71d1f42d0490abc7` |
| biome | yes | `1e088bb08936bff2791e8126e6e36fe287bfdc86274ddea3bf79d4a25b1e6bfe` |
| hydro | yes | `582d4e6e4d686418d7283b286079fbe4d81aac83c34ace329a10900989090c62` |
| scatter | yes | `34e21244b33b64281bca8e8a6a5a4e015ac7ce49b3899293d1d38a81fb3e821e` |

## Thread counts 1 / 4 / 8

| layer | threads=1 | threads=4 | threads=8 |
|---|---|---|---|
| noise | yes | yes | yes |
| erosion | yes | yes | yes |
| biome | yes | yes | yes |
| hydro | yes | yes | yes |
| scatter | yes | yes | yes |

## Seams (64x64 tiles stitched vs full world)

| layer | identical | mismatching bytes |
|---|---|---|
| noise.f32 | yes | 0 |
| erosion.f32 | yes | 0 |
| biome.u8 | yes | 0 |
| hydro.f32 | yes | 0 |
| scatter.f32 | yes | 0 |

## Master seed changes every layer

| layer | changed |
|---|---|
| noise | yes |
| erosion | yes |
| biome | yes |
| hydro | yes |
| scatter | yes |
