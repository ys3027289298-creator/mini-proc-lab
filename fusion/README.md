# fusion — 热力梯度与高程梯度的平滑融合

## 问题

朴素融合把高程梯度项用硬掩码切入热力场：

```
F0(p) = T(h(p)) + alpha * 1[s(p) > theta] * s(p),   s = |grad h|
```

不连续来自两处：

1. 高程 `h` 在悬崖处本身跳变（幅值 `A`）⇒ 热力项 `T(h)` 同步跳变 `lambda*A`；
2. 指示函数 `1[s > theta]` 在等值线 `s = theta` 上 0/1 跳变 ⇒ 融合项跳变 `alpha*theta`。

两处跳变在热力图上表现为硬边（阶跃像素）。

## 融合式推导

目标：构造融合场 `F`，使其梯度为热力梯度与高程梯度的**连续**加权和，从构造上排除阶跃。

**第 1 步 · 高程预平滑。** 用高斯核 `G_sigma` 卷积：

```
h~ = G_sigma * h ,  G_sigma ∈ C^inf  ⇒  h~ ∈ C^inf
```

幅值 `A` 的悬崖阶跃被展宽为约 `4*sigma` 像素的连续坡，最大坡度 `A * g_sigma`，
其中 `g_sigma = 1 / (sigma * sqrt(2 pi))` 为高斯核峰值。

**第 2 步 · 坡度场。** `s = |grad h~|` 连续（在 `grad h~ != 0` 处光滑）。

**第 3 步 · 逻辑权重代替指示函数。**

```
w(s) = 1 / (1 + exp(-(s - c) / (k * r)))
c = mean(s) + gamma * std(s)        （鲁棒中心，gamma = 0.5）
r = std(s)                          （鲁棒尺度，门限不依赖绝对单位）
```

`w ∈ C^inf`，`0 < w < 1`，`w' = w(1-w)/(k r)`，`max w' = 1/(4 k r)` 有界。

**第 4 步 · 融合场。**

```
F(p) = T(h~(p)) + alpha * w(s(p)) * s(p)
```

对 `p` 求梯度（链式法则）：

```
grad F = T'(h~) * grad h~  +  alpha * [ w(s) + s * w'(s) ] * grad s
         └─ 热力梯度项 ──┘    └──── 高程梯度项（经 w 连续调制） ────┘
```

右端每个因子（`T'`、`grad h~`、`w`、`w'`、`grad s`）都是连续场，
故 `grad F` 连续、`F ∈ C^1`：**不存在 0/1 切换，阶跃在构造上被排除**。

## 无阶跃的参数条件（上界）

相邻像素最大跳变量满足

```
max|dF| <= lambda * A * g_sigma + max|grad eta~|
          + alpha * (1 + s_max / (4 k r)) * A * g_sigma / sigma
```

（第一项为热力项经平滑悬崖的坡度；第二项由
`|grad(w*s)| <= (w + s*max w') * |grad s|` 与
`|grad s| <= max|Hess h~| ≈ A * g_sigma / sigma` 得到。）

`sigma` 增大时两项分别按 `1/sigma`、`1/sigma^2` 下降，
因此对任意阶跃阈值 `tau` 都存在使 `max|dF| < tau` 的 `(sigma, k)`。
阶跃像素定义为：四邻域差分超过 `tau = 0.1 * (P99 - P1)(F0)` 的像素。

## 三档融合参数与阶跃像素统计

场景：256×256 高度场，含幅值 `A = 5` 的波浪形悬崖；`lambda = 1.2`，`alpha = 2.0`；
`tau = 0.6521`。

| 场 | sigma | k | 阶跃像素 | max\|dF\| |
|---|---:|---:|---:|---:|
| naive（硬掩码） | — | — | **1044** | 11.204 |
| tier1_light | 4.0 | 0.10 | **0** | 0.634 |
| tier2_mid | 8.0 | 0.06 | **0** | 0.296 |
| tier3_heavy | 16.0 | 0.03 | **0** | 0.153 |

三档 `max|dF|` 均小于 `tau`，阶跃像素为 0；`sigma` 越大悬崖过渡带越宽、
残余梯度越小，与上界的 `1/sigma` 衰减一致。

## 复现

```
make -C fusion run     # 渲染 + 自检，任一一档阶跃像素非零则退出码为 1
```

产物（`fusion/out/`）：`fusion_naive.png`、`fusion_tier1_light.png`、
`fusion_tier2_mid.png`、`fusion_tier3_heavy.png`（三档各一张）、
`fusion_compare.png`（2×2 热力对比图，阶跃像素以青色标出）、`stats.csv`。
