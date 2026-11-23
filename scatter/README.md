# scatter

`make test` 跑自检（`scatter` 原有抖动网格 + `bluenoise` 蓝噪声采样器）。

## bluenoise：周期环面上的采样器

所有距离在单位环面上周期化（wrap），坐标输出 [0,1)，首行 `# meta` 头。

```
bluenoise poisson N seed            # iid 均匀（白噪声参考）
bluenoise dart r Nmax seed [pat]    # 全域平稳 dart throwing，硬核心 r
                                    # Nmax=0 跑到饱和（连续 pat 次拒绝）
bluenoise bridson r Nmax seed [k]   # Bridson 快速 Poisson-disk
bluenoise grid n                    # n*n 方晶格（频谱参考）
bluenoise selftest                  # 内置断言
```

选型要点：

- **固定 N 的蓝噪声用 `dart`**：候选在整个环面上均匀，过程平稳，任意接受
  数下都铺满全域。
- **`bridson` 不要用于 stop-at-N**：候选只在母点 [r,2r] 环内生成，r 远小于
  域且 N 远未饱和时留下生长前沿（点集缩成一团，低频假性抬升）。保留它用于
  饱和探测与教学对照；饱和计数本身也偏保守（active list 早停）。
- 背景网格单元边 ≤ r/√2，冲突检测必须搜 5×5 邻域（3×3 会漏检，硬核心失效）。
- 硬核心密度界：方格 ρ≤1/r²，六方最密排 ρ_max=2/(√3 r²)≈1.1547/r²；
  dart throwing 饱和效率约 0.60（二维 RSA jamming）。

分析管线与蓝噪声频谱推导见 `../artifacts/analysis/`（REPORT.md）。
