# erosion

基于液滴的水力侵蚀，跑在 fBm 初始高度场上。`make test` 做一次短跑自检。

## 模型

单次迭代 = 在整张场上投 `drops` 个液滴。每个液滴：

1. 在随机格点落地，初始水量 1、速度 1、泥沙 0。
2. 按「惯性 × 原方向 − (1−惯性) × 梯度」下行（双线性梯度）。
3. 容量 `cap = capacity · speed · water`；沉沙 < 容量则从半径 `radius` 的刷子侵蚀，沉沙 > 容量则沉积在刷子内。
4. 水按 `evaporate` 蒸发，水耗尽或到寿命 `max_steps` 终止；剩余泥沙沉积在终止点（质量守恒）。

每次迭代的净改动先攒在 `tmp` 里、结束时一次性叠加，所以同一迭代内所有液滴看到的是同一张地表（无顺序依赖）。

**关于沟壑稳定性**：早期版本在叠加 delta 后加 0.02 的拉普拉斯平滑项，但实测它只会抹平 delta，对沟壑深度增长没有约束（深坑照样长到 −330）。现在保留 `--smooth`（默认 0）作为可选项，但默认关闭；稳定性的实际控制量是 `--capacity`——容量越大，单步侵蚀量越大，越容易在少数格点失稳。

## 测量接口

```
./erosion --n 256 --drops 65536 --iters N --seed 7 --out PREFIX [--resample]
```

- 每次迭代写 `PREFIX.<iter>.h32`（little-endian float32，行优先，边长 `n`），`0000` 是初始 fBm 面。
- stdout 输出一行 JSON：参数 + `runs[]`，每项含：
  - `iter`、`ms`（本次迭代墙钟）、`us_per_cell`
  - `residual`：与上一迭代的逐格高度差 RMS（收敛律就测这个）
  - `rms_vs_init`：与初始面的逐格高度差 RMS（看总形变）
  - `mean_abs_delta`：逐格平均 |Δh|
  - `deposit_vol` / `erode_vol`：本迭代沉积/侵蚀的总体积（格·高度单位）
- stderr 同步打印人读的每迭代行。
- 默认每次迭代复用同一组液滴落点（`sample_seed=0`），让 residual 只反映 pass 算子的收敛；`--resample` 改为每次重抽，用于对比“采样噪声”成分。

## 参数

| flag | 默认 | 含义 |
|---|---|---|
| `--max-steps` | 48 | 液滴寿命 |
| `--inertia` | 0.05 | 方向惯性 |
| `--capacity` | 0.05 | 携沙容量系数 |
| `--erode` | 0.3 | 侵蚀率 |
| `--deposit` | 0.3 | 沉积率 |
| `--evaporate` | 0.01 | 每步蒸发 |
| `--gravity` | 4.0 | 每单位落差的速度增益 |
| `--radius` | 2.0 | 侵蚀刷半径（格距） |
| `--smooth` | 0.0 | 对每步 delta 的拉普拉斯平滑（可选） |
| `--resample` | off | 每次迭代重抽液滴落点 |
