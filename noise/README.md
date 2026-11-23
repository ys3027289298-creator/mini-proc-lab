# noise

`make test` 跑自检。

- `noise_white01`：坐标哈希白噪声，[0, 1)，确定性。
- `noise_value2`：二维值噪声，smoothstep 过渡，映射到 [-1, 1]。
- `noise_fbm2`：多倍频叠加并归一化，供地形生成使用。
