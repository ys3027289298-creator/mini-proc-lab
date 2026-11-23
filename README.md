# mini-proc-lab

程序化地形与散布基础模块，C 实现，每目录独立构建可自测。

| 目录 | 内容 |
|---|---|
| `noise/` | 值噪声、fBm 叠加、确定性哈希 |
| `terrain/` | 高度场双线性采样、坡度 |
| `scatter/` | 抖动网格散布、最小间距拒绝；周期环面 dart/Bridson Poisson-disk、白噪声与晶格采样器 |
| `artifacts/analysis/` | 蓝噪声频谱—最小间距推导、三档散布测量（S(k)、最近邻分布、团簇计数）与图表，见该目录 REPORT.md |
