# postfx — 逐通道守恒校验推导

## 1. 问题定义

四个后处理 pass 顺序叠加，每个 pass 在每个通道上对全部像素作用同一个仿射映射：

| pass | 通道 $c$ 上的映射 | 说明 |
|---|---|---|
| gain | $X \mapsto g_c X$ | 逐通道增益，$b=0$ |
| lift | $X \mapsto X + l_c$ | 逐通道常量抬升 |
| overlay | $X \mapsto X + o_c(x,y)$ | 空间场叠加，$b$ 随像素变化 |
| tint | $X \mapsto t_c X + \tau_c$ | 色调仿射；**缺陷版把 G 通道写成 $X \mapsto \tau_G$（上游被丢弃）** |

"通道被吞"的含义：缺陷 tint 对 G 通道是**覆盖**而非仿射写入，$a_G s_G + b_G$ 整项丢失，只剩常量 $\tau_G$。

输入记为 $s_c(x,y)$，像素总数 $N_p=1600$（40×40）。能量定义：

$$E_2(X_c)=\sum_{p=1}^{N_p} X_c(p)^2,\qquad E_1(X_c)=\sum_{p=1}^{N_p} X_c(p)$$

## 2. 解析累加（两条独立路径）

### 2.1 路径 A：逐像素仿射复合

仿射映射对复合封闭。令单通道上 $X\mapsto aX+b$、再叠加 $X\mapsto a'X+b'$，则

$$X\mapsto a'a\,X + (a'b+b')$$

逐 pass 递推系数（每通道独立，下标 $c$ 省略）：

| pass | $A$（乘子） | $B$（常量/空间项） |
|---|---|---|
| 起点 | $1$ | $0$ |
| gain | $A\leftarrow gA$ | $B\leftarrow gB$ |
| lift | 不变 | $B\leftarrow B+l$ |
| overlay | 不变 | 逐像素 $B(x,y)\leftarrow B(x,y)+o(x,y)$ |
| tint | $A\leftarrow tA$ | $B\leftarrow tB+\tau$ |

整栈的解析输出（$o$ 的空间贡献同样被 tint 的 $t$ 缩放）：

$$\boxed{\;\hat X_c(p)=A_c\,s_c(p)+B_c+t_c\,o_c(p)\;}$$

### 2.2 路径 B：闭合式能量簿记（不碰任何像素值）

对任意 $Y=aX+b$（$b$ 可以是空间场），由内积线性性：

$$\boxed{\;E_2(Y)=a^2E_2(X)+2a\langle X,b\rangle+E_2(b)\;},\qquad E_1(Y)=aE_1(X)+\sum_p b(p)$$

簿记器只维护五个标量/矩：$E_2$、$E_1$、$\sum b$、$\sum b^2$、$\langle X,b\rangle$，
每个 pass 用上式递推一次，全程不构造输出图像。它与路径 A 用不同的运算次序回答同一问题，
因此二者一致才构成交叉校验，而不是同一算法的复读。

缺陷栈的 G 通道在最后一步违反前提（不再是 $aX+b$ 而是常量 $\tau_G$），簿记仍按正确仿射
递推，于是簿记值代表"**应当**得到的能量"，实测值与之偏离即定位缺陷。

### 2.3 吞通道的解析预测

缺陷 tint 后 G 通道恒为 $\tau_G=0.03$：

$$E_2^{bug}(G)=\sum_p \tau_G^2=N_p\,\tau_G^2=1600\times0.0009=\mathbf{1.44}$$

正确值 $409.8564$，预测相对差 $1-1.44/409.8564\approx\mathbf{99.65\%}$，与实测完全吻合。

## 3. 校验判据

逐像素（Python 端 float64 复核）：

$$r_p=\dfrac{|X^{meas}(p)-\hat X(p)|}{\max(10^{-6},\,|\hat X(p)|)},\qquad \max_p r_p \le 10^{-3}$$

逐通道能量（C 端 double 求和）：

$$r_E=\dfrac{|E_2^{meas}-E_2^{analytic}|}{\max(|E_2^{meas}|,|E_2^{analytic}|)}\le 10^{-3}$$

**任一通道任一判据超过 $10^{-3}$ 即未完成。**

## 4. 实测结果（40×40，确定性输入）

守恒栈三通道、全部五个阶段，实测 / 逐像素解析 / 闭合簿记三方一致，最差相对差
$8.6\times10^{-8}$（纯 float32 舍入量级，距 $10^{-3}$ 阈值有 4 个数量级余量）：

| 通道 | $E_2$ 输入 | $E_2$ 输出（实测） | 倍率 | 解析输出 | 相对差 | 判定 |
|---|---:|---:|---:|---:|---:|---|
| R | 927.5962 | 1735.9645 | 1.87147× | 1735.9646 | 4.62e-08 | PASS |
| G | 333.5809 | 409.8564 | 1.22866× | 409.8564 | 8.46e-08 | PASS |
| B | 497.6747 | 948.5511 | 1.90597× | 948.5512 | 1.91e-08 | PASS |

缺陷栈：R/B 与解析一致（4.6e-08、1.9e-08），**G 通道 $E_2$ 从 409.86 塌缩到 1.44**，
逐像素最大相对差 $2.75\times10^{2}$，能量相对差 $9.97\times10^{-1}$，吞通道被精确定位到 G。

## 5. 三条数字路径的独立性

| 路径 | 构造方式 | 浮点路径 |
|---|---|---|
| measured | 四个 pass 的真实像素循环，含缺陷分支 | float32 逐像素 |
| analytic | 从输入定义重算 $s_c(p)$，只走复合系数 $A,B$ | float32，运算次序不同 |
| book | 只递推矩与点积，从不构造输出图像 | double 累加 |

三者共享的只有 pass 参数表与输入定义；任一路径的实现错误都不会被另两条复制。
完整分阶段数据见 `../../artifacts/postfx/energy_table.md`，逐像素图见
`diffmap_conserving.png` 与 `diffmap_buggy.png`。

## 6. 复现

```
cd postfx && make test      # C 自检 + 二进制 dump + Python 渲染
```

产物落在 `artifacts/postfx/`：`report.txt`（控制台全量输出）、`postfx.bin`（PFX1  bundle）、
`energy_table.md`、`stats.json`、两张 PNG。
