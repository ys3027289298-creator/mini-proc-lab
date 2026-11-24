# 矿道中线强制对齐地形梯度 —— 最大允许横向偏移推导

## 1. 记号与问题定义

高度场 $z=h(x,y)$，其梯度为

$$
\mathbf g(x,y)=\nabla h=\left(\tfrac{\partial h}{\partial x},\tfrac{\partial h}{\partial y}\right),\qquad
\hat{\mathbf g}=\mathbf g/\|\mathbf g\|
$$

**理想矿道中线 = 地形梯度流流线**（沿最陡上坡方向），满足 ODE

$$
\frac{d\mathbf r}{ds}=\hat{\mathbf g}(\mathbf r(s)),\qquad \mathbf r(0)=\mathbf r_0
$$

其中 $s$ 为弧长。实际施工把中线铺成一串直线弦（折线，相邻测量桩之间为直线段），弦与理想流线之间的横向距离即**横向偏差 $d$**。

问题：给定走向偏角上限 $\alpha$ 与弦的放样/铺轨容差 $e_{\text{chord}}$，求每个施工段允许的最大横向偏移 $d_{\text{allow}}(k)$。

## 2. 走向偏角约束给出的横向偏移上界

设第 $k$ 段流线局部曲率为 $\kappa_k$（曲率半径 $R_k=1/\kappa_k$）。弦两端都落在流线上，弦切向相对流线切向的偏角记为 $\theta$（即走向偏差）。在圆弧近似下（一段内 $\kappa$ 取段内最大值，近似偏保守），弦对圆心的张角为 $2\theta$。

### 2.1 几何关系

圆弧、弦与矢距（sagitta，即弦到弧的最大横向距离 $d$）：

$$
d=R\,(1-\cos\theta)=\frac{1-\cos\theta}{\kappa}
$$

反解走向偏角：

$$
\theta(d)=\arccos(1-\kappa\,d)
$$

### 2.2 代入走向偏角上限 $\alpha$

要求 $\theta\le\alpha$，且 $1-\cos(\cdot)$ 在 $[0,\pi/2]$ 单调，于是

$$
d\;\le\;\frac{1-\cos\alpha}{\kappa_k}
$$

这就是**走向约束推出的最大允许横向偏移**。小角度展开 $1-\cos\alpha=\alpha^2/2-\alpha^4/24+\cdots$：

$$
d_{\text{head}}(k)=\frac{1-\cos\alpha}{\kappa_k}\;\approx\;\frac{\alpha^2}{2\,\kappa_k}\quad(\alpha\text{ 以弧度计})
$$

物理含义：曲率越大（高低差大、流线弯得急），允许偏移越小——这正是“矿道在高低差大处生歪”的数学原因，偏移上界与局部曲率成反比。

## 3. 放样容差约束

弦本身的测量/铺轨存在独立的直线放样容差 $e_{\text{chord}}$（本工程取 $e_{\text{chord}}=0.20\ \text{m}$，平直段 $\kappa\to0$ 时走向约束失效，全靠它兜底）。

## 4. 最大允许偏移表达式（交付式）

两条约束取更严者（min）：

$$
\boxed{\;d_{\text{allow}}(k)=\min\!\left(e_{\text{chord}},\;\dfrac{1-\cos\alpha}{\kappa_k}\right)\;}
$$

- $\alpha=8^\circ=0.13963\ \text{rad}$，$1-\cos\alpha=0.009732$；
- $\kappa_k$ = 第 $k$ 施工段内理想流线的**最大**曲率（用段内最大曲率而非均值，保证整段任何位置都达标）；
- 平直段 $\kappa_k=0$ 时第二项趋于无穷，表达式自动回退到 $e_{\text{chord}}$。

**逐段验收判据（任一段超标即未完成）**：

$$
d_{\max}(k)\le d_{\text{allow}}(k)\quad\forall k
\qquad\Longleftrightarrow\qquad
\max_k\frac{d_{\max}(k)}{d_{\text{allow}}(k)}\le 1
$$

其中 $d_{\max}(k)$ 为该段内实测横向偏差的最大值。同时直接复核走向偏角 $\theta_{\max}(k)\le\alpha$。

## 5. 曲率的离散计算

流线上点 $\mathbf r_i$ 的弧长步长为 $\Delta s$，用相邻切向 $\hat{\mathbf g}_i$ 的夹角估计：

$$
\kappa_i=\frac{\arccos\!\big(\hat{\mathbf g}_{i}\!\cdot\hat{\mathbf g}_{i+1}\big)}{\Delta s}
$$

段曲率取 $\kappa_k=\max_{i\in k}\kappa_i$。

## 6. 强制对齐的施工算法（本站方案）

旧方案（生歪来源）：等距 $L=4$ m 打桩，每桩只在该点采样一次梯度方向并直走，弦长固定 → 急弯处弦矢距超过 $d_{\text{allow}}$，且误差逐桩累积（实测三段地形分别有 7/10、7/10、9/10 段超标）。

新方案（自适应弦长，保证逐段达标）：

1. 以 $\Delta s=0.002$ m 精细积分梯度流（RK2 中点法），得到理想流线 $\{\mathbf r_i\}$ 与每点曲率 $\kappa_i$；
2. 沿流线逐节点延长当前弦，**直接实测**候选弦对窗口内所有流线节点的最大矢距（不依赖常曲率解析近似，因为一个窗口内 $\kappa$ 是变化的），当实测矢距达到段预算的安全比例 $\eta=0.9$ 时即收桩；平直段矢距始终很小，弦长封顶 $L_{\max}=2$ m：
   $$
   \text{close when}\quad \max_{j\in[i_0,i_1]}\operatorname{dist}\!\big(\mathbf r_j,\overline{\mathbf r_{i_0}\mathbf r_{i_1}}\big)\ge \eta\,d_{\text{allow}}\!\big(\max_{j\in[i_0,i_1]}\kappa_j\big)
   $$
   常曲率下的等价解析估计为 $L_k\approx 2\sqrt{2\,\eta\,d_{\text{allow}}(k)/\kappa_k}$（由 $d\approx\kappa L^2/8$ 反解），实现用实测值代替该估计只会更保守；
3. 桩点全部锚定在理想流线上（独立复核：最大锚定间距 0 m），弦为桩间直线——构造上保证 $d_{\max}(k)\le\eta\,d_{\text{allow}}(k)<d_{\text{allow}}(k)$；
4. 再独立实测每条弦的横向偏差与走向偏角做验收（`tools/verify_acceptance.py` 从 CSV 原始数据重算，不信任 C 写出的 pass 标志）。
