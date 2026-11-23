# terrain

`make test` 跑自检。

- 定容高度场（`TERRAIN_MAX` 见头文件），双线性世界坐标采样，越界钳制。
- `field_slope` 顶点中心差分坡度，`field_slope_at` 给出双线性曲面在任意世界位置的坡度（rise/run）。
- `field_generate_fbm` 用 `../noise` 的 fBm 生成高度场，供标定与演示。
