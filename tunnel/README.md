# tunnel

矿道中线强制对齐地形梯度：把理想中线定义为高度场梯度流流线，按推导的
逐段横向偏移预算做自适应弦长放样，并在三段难度递增的解析地形上验收。

- 推导与公式：[`DERIVATION.md`](DERIVATION.md)
  - 核心交付式：`d_allow(k) = min(e_chord, (1-cos α)/κ_k)`，α = 8°，e_chord = 0.20 m
  - 任一段 `d_max(k) > d_allow(k)` 即整段未完成
- C 核心（无第三方依赖，`-lm`）：`src/tunnel.c`
  - 三段解析地形 + 解析梯度（含有限差分对拍自检）
  - RK2 梯度流积分（2 mm 弧长步长）、宽模板 atan2 曲率
  - 旧方案（固定 4 m Euler 弦）与本站方案（自适应弦长、桩点锚定流线）
  - 逐段实测矢距、走向偏角与预算判定；`make test` 跑全部自检

```sh
make test                 # 自检（all_ok=1 才通过）
./tunnel <out_dir>        # 额外输出 CSV：流线、两种方案桩点、逐段统计、汇总
```

## Python 工具（matplotlib/numpy，独立于 C 自检）

```sh
python3 tools/plot_alignment.py <data_dir> <out.png> <stats.csv>
python3 tools/verify_acceptance.py <data_dir>   # 从 CSV 重算，不信任 C 的 pass 标志
```

## 三段地形结果（40 m 中线）

| 地形 | 落差 m | 曲率峰值 /m | 段数 | d_max max m | d_allow min m | 上限利用率 | 走向偏差 max | 判定 |
|---|---|---|---|---|---|---|---|---|
| rolling_hills | 10.5 | 0.042 | 20 | 0.020 | 0.200 | 10% | 2.3° | 达标 |
| single_ridge | 20.9 | 0.125 | 20 | 0.061 | 0.078 | 78% | 7.0° | 达标 |
| steep_bent_monocline | 49.5 | 0.228 | 21 | 0.051 | 0.043 | 90% | 7.6° | 达标 |

对照旧方案：最大横向偏差 0.57 / 2.02 / 1.77 m，超标段 7/10、7/10、9/10。
图件见 `../artifacts/mine_alignment.png`，机器可读统计见
`../artifacts/deviation_stats.csv`。
