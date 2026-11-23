# thermal

点状热异常 → 有限高斯核斑块。`make test` 跑自检并输出三档面积统计与对比图 `thermal_compare.pgm`。

## 扩散推导（核半径的理论依据）

二维热传导方程，初始时刻在原点注入总能量 $E$ 的瞬时点源：

$$\frac{\partial T}{\partial t} = D\,\nabla^2 T,\qquad T(\mathbf{r},0)=E\,\delta(\mathbf{r})$$

其格林函数（基本解）为二维高斯：

$$T(r,t)=\frac{E}{4\pi D t}\exp\!\Big(-\frac{r^2}{4Dt}\Big),\qquad \sigma^2 = 2Dt$$

即扩散把点源抹成宽度 $\sigma=\sqrt{2Dt}$ 的高斯斑，峰值 $A=E/(2\pi\sigma^2)$ 随时间衰减。

**斑块区域**取可观测等温线 $T \ge \theta$（$\theta$ 为探测阈值），解出核半径：

$$\boxed{\;R=\sigma\sqrt{2\ln\frac{E}{2\pi\sigma^2\theta}}\;},\qquad \text{面积 } \mathcal{A}=\pi R^2=2\pi\sigma^2\ln\frac{E}{2\pi\sigma^2\theta}$$

可观测条件 $E > 2\pi\sigma^2\theta$（峰值必须高于阈值），否则斑块不存在——这正是"异常点太小观测不到"的定量表述。

**有限核截断**：半径 $r$ 内包含的能量份额为 $1-e^{-r^2/2\sigma^2}$。取 $R=k\sigma$ 截断，丢弃的能量份额仅 $e^{-k^2/2}$：

| k | 保留能量 | 丢弃 |
|---|---|---|
| 2 | 86.5% | 13.5% |
| 3 | 98.9% | 1.1% |
| 4 | 99.97% | 0.03% |

故核取 $R=k\sigma,\ k\ge 2$ 即有紧支撑且能量损失可控——这是核半径选取的理论依据。

**三档能量**：面积对 $E$ 只是对数增长，要让半径线性增长（$R=k\sigma,\ k=2,3,4$），能量必须指数增长：

$$E_k = 2\pi\sigma^2\theta\, e^{k^2/2}$$

本模块取 $\sigma=8$ 格、$\theta=1$，三档 $k=2,3,4$，对应 $R=16,24,32$ 格、$\mathcal{A}=256\pi,576\pi,1024\pi$。

## 自检内容

- `radius_tierN` / `area_tierN`：超采样统计的斑块半径/面积与理论值偏差 < 0.5% / 1%
- `energy_tierN`：有限核网格积分能量 = $E(1-e^{-k^2/2})$
- `sigma2_moment`：由二阶矩 $\langle r^2\rangle=2\sigma^2$ 反解 $\sigma^2$，验证扩散宽度
- `point_invisible`：未扩展的点异常仅占 1 格，面积不足一档斑块的 1%

## 对比图

`thermal_compare.pgm`：每行一档能量，左列为未扩展的点异常（单像素），右列为高斯核斑块，灰色圆环为理论等温线半径 $R$——栅格化斑块边界应与环重合。
