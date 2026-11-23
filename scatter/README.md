# scatter

`make test` 跑自检。

- 抖动网格散布 + 最小间距拒绝，全部由坐标哈希驱动，可复现。
- `SlopeDensity` 给出坡度密度衰减模型 `rho0 * exp(-(s/s0)^p)`，`jitter_grid` 逐候选按该概率保留，陡坡自然变稀。
- `scatter_measure` 不保存点位、只按坡度分箱累计候选数/保留数，供 `../slope` 拟合使用。
- 候选哈希按二维坐标生成（`key2`），抖动与保留判定互不相关；`tests/test_internal.c` 对该性质有回归断言。
