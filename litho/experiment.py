#!/usr/bin/env python3
# litho/experiment.py
#
# 岩性纹理辨识度实验:
#   1. 每种岩性一个噪声特征向量 f = (log2f0, octaves, gain, lacunarity, anisotropy, ridge)
#   2. 标定实验推导特征距离 d 与视觉辨识度 A 的函数关系 A(d)
#   3. 三组特征相近岩性各做一次识别测试(3 类强制选择), 阈值 0.80
#
# 输出到 litho/out/: calibration_curve.png, confusion_matrix.png,
#                    confusion_matrix.csv, collage.png, results.json
import json
import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
os.makedirs(OUT, exist_ok=True)

SIZE = 64           # 纹理块边长(有限采样 -> 类内涨落, 模拟有限视觉采样)
THRESHOLD = 0.80    # 辨识度阈值(3 类, 随机水平 1/3)
MARGIN = 0.05       # 验收余量: 实测需 >= THRESHOLD + MARGIN 才接受, 避免贴线通过
N_TRAIN = 48        # 每类训练块数
N_TEST = 96         # 每类测试块数
WEBER = 0.30        # 观察者内部噪声(韦伯分数): 每次观测 feat *= (1 + WEBER*N(0,1))

# 特征各维度的设计空间范围, 用于归一化距离
F_MIN = np.array([2.0, 1.0, 0.30, 1.5, 1.0, 0.0])
F_MAX = np.array([5.0, 6.0, 0.80, 3.0, 4.0, 1.0])
F_RANGE = F_MAX - F_MIN


def feat_dist(a, b):
    """归一化特征空间中的欧氏距离."""
    return float(np.linalg.norm((np.asarray(a) - np.asarray(b)) / F_RANGE))


# ---------------------------------------------------------------- 噪声与纹理

def _hash01(ix, iy, seed):
    h = (ix * 374761393 + iy * 668265263 + seed * 2246822519) & 0xFFFFFFFF
    h = ((h ^ (h >> 13)) * 1274126177) & 0xFFFFFFFF
    h = h ^ (h >> 16)
    return (h & 0xFFFFFF) / 16777215.0


def value_noise(x, y, seed):
    xi = np.floor(x).astype(np.int64)
    yi = np.floor(y).astype(np.int64)
    fx = x - xi
    fy = y - yi
    v00 = _hash01(xi, yi, seed)
    v10 = _hash01(xi + 1, yi, seed)
    v01 = _hash01(xi, yi + 1, seed)
    v11 = _hash01(xi + 1, yi + 1, seed)
    u = fx * fx * (3 - 2 * fx)
    v = fy * fy * (3 - 2 * fy)
    return ((v00 * (1 - u) + v10 * u) * (1 - v)
            + (v01 * (1 - u) + v11 * u) * v) * 2 - 1


def fbm(x, y, octaves, lacunarity, gain, seed):
    """octaves 可为小数: 最后一个倍频程按小数部分加权."""
    full = int(np.floor(octaves))
    frac = octaves - full
    total, norm = 0.0, 0.0
    amp, freq = 1.0, 1.0
    for i in range(full):
        total += value_noise(x * freq, y * freq, seed + i * 17) * amp
        norm += amp
        amp *= gain
        freq *= lacunarity
    if frac > 1e-6:
        total += value_noise(x * freq, y * freq, seed + full * 17) * amp * frac
        norm += amp * frac
    return total / norm if norm > 0 else np.zeros_like(x)


def texture(feat, seed):
    """按特征向量生成一幅 SIZE x SIZE 纹理, 值域 [0,1]."""
    log2f0, octaves, gain, lac, aniso, ridge = feat
    f0 = 2.0 ** log2f0
    yy, xx = np.mgrid[0:SIZE, 0:SIZE].astype(np.float64)
    x = xx * f0 / SIZE
    y = yy * f0 / SIZE * aniso          # 各向异性拉伸(层理/片理方向)
    n = fbm(x, y, octaves, lac, gain, seed)
    ridged = 1.0 - 2.0 * np.abs(n)      # 脊状变换
    t = (1 - ridge) * n + ridge * ridged
    return np.clip(t * 0.5 + 0.5, 0.0, 1.0)


# ---------------------------------------------------------------- 视觉观察者

def observe(img):
    """前注意视觉统计特征: 灰度矩 / 梯度 / 径向功率谱 / 谱方向性."""
    feats = [img.mean(), img.std()]
    gx, gy = np.gradient(img)
    gm = np.hypot(gx, gy)
    feats += [gm.mean(), gm.std()]
    p = np.abs(np.fft.fft2(img - img.mean())) ** 2
    p = np.fft.fftshift(p)
    fy, fx = np.mgrid[0:SIZE, 0:SIZE] - SIZE / 2
    r = np.hypot(fx, fy)
    total = p.sum() + 1e-12
    for lo, hi in [(2, 8), (8, 16), (16, 32), (32, 64)]:
        feats.append(p[(r >= lo) & (r < hi)].sum() / total)
    ang = np.abs(np.arctan2(fy, fx))
    horiz = p[(ang < np.pi / 8) | (ang > np.pi - np.pi / 8)].sum()
    vert = p[np.abs(ang - np.pi / 2) < np.pi / 8].sum()
    feats.append(horiz / (vert + 1e-12))
    return np.array(feats)


def make_patches(feat, n, seed0):
    rng = np.random.default_rng(seed0 ^ 0x5bd1e995)
    out = []
    for i in range(n):
        f = observe(texture(feat, seed0 + i))
        f = f * (1 + WEBER * rng.standard_normal(f.shape))  # 视觉内部噪声
        out.append(f)
    return np.array(out)


def nc_accuracy(train_x, train_y, test_x, test_y):
    """最近质心分类器(训练集标准化) -> (准确率, 逐类预测)."""
    mu, sd = train_x.mean(0), train_x.std(0) + 1e-12
    tr = (train_x - mu) / sd
    te = (test_x - mu) / sd
    cents = np.array([tr[train_y == c].mean(0) for c in np.unique(train_y)])
    d = ((te[:, None, :] - cents[None, :, :]) ** 2).sum(-1)
    pred = np.unique(train_y)[d.argmin(1)]
    return (pred == test_y).mean(), pred


# ---------------------------------------------------------------- 岩性定义
# f = (log2f0, octaves, gain, lacunarity, anisotropy, ridge)
GROUPS = {
    "A_intrusive": {  # 深成岩: 各向同性, 粒度和粗糙度相近
        "granite":  [3.20, 4.0, 0.50, 2.0, 1.0, 0.15],
        "diorite":  [3.40, 4.2, 0.53, 2.0, 1.0, 0.20],
        "gabbro":   [3.60, 4.4, 0.56, 2.0, 1.0, 0.25],
    },
    "B_sedimentary": {  # 沉积岩: 均带层理(各向异性), 层理频率相近
        "sandstone": [2.80, 3.0, 0.55, 2.2, 1.6, 0.10],
        "siltstone": [3.10, 3.4, 0.58, 2.2, 2.0, 0.12],
        "shale":     [3.40, 3.8, 0.61, 2.2, 2.4, 0.14],
    },
    "C_carbonate": {  # 碳酸盐/变质: 脊状纹理为主
        "limestone": [3.00, 4.5, 0.60, 1.9, 1.2, 0.35],
        "dolomite":  [3.20, 4.8, 0.63, 2.0, 1.3, 0.42],
        "marble":    [3.40, 5.1, 0.66, 2.1, 1.4, 0.49],
    },
}


# ---------------------------------------------------------------- 1. 标定实验

def calibration():
    """参考岩性 + 沿随机方向拉开距离 d 的变体, 测 2 类识别率, 拟合
    A(d) = 1/2 + 1/2 * (1 - exp(-(d/d0)^p))."""
    rng = np.random.default_rng(7)
    f_ref = np.array([3.5, 4.0, 0.55, 2.0, 1.0, 0.2])
    distances = [0.03, 0.06, 0.10, 0.15, 0.22, 0.32, 0.45, 0.65, 0.90, 1.20]
    n_dir = 6
    rows = []
    ref_train = make_patches(f_ref, 24, 100000)
    ref_test = make_patches(f_ref, 24, 200000)
    for d in distances:
        accs = []
        for k in range(n_dir):
            direction = rng.normal(size=6)
            direction /= np.linalg.norm(direction)
            f_var = f_ref + direction * d * F_RANGE
            f_var = np.clip(f_var, F_MIN, F_MAX)
            d_actual = feat_dist(f_ref, f_var)
            tr_x = np.vstack([ref_train, make_patches(f_var, 24, 300000 + k)])
            te_x = np.vstack([ref_test, make_patches(f_var, 24, 400000 + k)])
            y_tr = np.repeat([0, 1], 24)
            y_te = np.repeat([0, 1], 24)
            acc, _ = nc_accuracy(tr_x, y_tr, te_x, y_te)
            accs.append(acc)
            rows.append((d_actual, acc))
        print(f"  calib d={d:.2f}  acc={np.mean(accs):.3f}")
    rows = np.array(rows)
    # 网格搜索拟合 d0, p
    best = None
    for d0 in np.linspace(0.02, 2.0, 200):
        for p in np.linspace(0.3, 4.0, 200):
            pred = 1 - 0.5 * np.exp(-((rows[:, 0] / d0) ** p))
            sse = ((pred - rows[:, 1]) ** 2).sum()
            if best is None or sse < best[0]:
                best = (sse, d0, p)
    _, d0, p = best
    pred = 1 - 0.5 * np.exp(-((rows[:, 0] / d0) ** p))
    ss_res = ((rows[:, 1] - pred) ** 2).sum()
    ss_tot = ((rows[:, 1] - rows[:, 1].mean()) ** 2).sum()
    r2 = 1 - ss_res / ss_tot
    print(f"  fit: d0={d0:.3f} p={p:.3f} R2={r2:.3f}")

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.scatter(rows[:, 0], rows[:, 1], s=14, alpha=0.5, label="measured")
    dd = np.linspace(0, 1.6, 200)
    ax.plot(dd, 1 - 0.5 * np.exp(-((dd / d0) ** p)), "r-",
            label=f"A(d)=1-0.5*exp(-(d/{d0:.3f})^{p:.2f})  R2={r2:.3f}")
    ax.axhline(0.5, ls=":", c="gray")
    ax.set_xlabel("feature distance d (normalized)")
    ax.set_ylabel("2-class accuracy A")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "calibration_curve.png"), dpi=130)
    plt.close(fig)
    return d0, p, r2


def predict_group_acc(d_min, d0, p, k=3):
    """由标定关系外推 K 类强制选择辨识度."""
    return 1.0 / k + (1 - 1.0 / k) * (1 - np.exp(-((d_min / d0) ** p)))


# ---------------------------------------------------------------- 2. 三组识别测试

def run_group(name, lithos, d0, p):
    """组内 3 类识别测试; 低于阈值则把特征沿离心方向拉开后重测."""
    scale = 1.0
    history = []
    for it in range(8):
        feats = {k: np.array(v, float) for k, v in lithos.items()}
        cen = np.mean(list(feats.values()), axis=0)
        feats = {k: np.clip(cen + (f - cen) * scale, F_MIN, F_MAX)
                 for k, f in feats.items()}
        names = list(feats)
        train = {k: make_patches(feats[k], N_TRAIN, 500000 + 977 * it + i)
                 for i, k in enumerate(names)}
        test = {k: make_patches(feats[k], N_TEST, 700000 + 977 * it + i)
                for i, k in enumerate(names)}
        tr_x = np.vstack([train[k] for k in names])
        te_x = np.vstack([test[k] for k in names])
        y_tr = np.repeat(range(3), N_TRAIN)
        y_te = np.repeat(range(3), N_TEST)
        acc, pred = nc_accuracy(tr_x, y_tr, te_x, y_te)
        cm = np.zeros((3, 3), int)
        for t, q in zip(y_te, pred):
            cm[t, q] += 1
        dmin = min(feat_dist(feats[a], feats[b])
                   for i, a in enumerate(names) for b in names[i + 1:])
        history.append((it, scale, acc, dmin))
        print(f"  {name} it={it} scale={scale:.2f} d_min={dmin:.3f} acc={acc:.3f}")
        if acc >= THRESHOLD + MARGIN:
            return dict(names=names, feats={k: list(map(float, feats[k])) for k in names},
                        acc=float(acc), cm=cm.tolist(), d_min=dmin,
                        pred_acc=float(predict_group_acc(dmin, d0, p)),
                        history=history, test=test, scale=scale)
        scale *= 1.4
    return None  # 达到迭代上限仍未过阈值


def main():
    print("== calibration ==")
    d0, p, r2 = calibration()

    print("== group tests ==")
    results = {}
    for name, lithos in GROUPS.items():
        r = run_group(name, lithos, d0, p)
        if r is None:
            print(f"  {name}: FAILED to reach threshold")
            results[name] = dict(failed=True)
        else:
            results[name] = r

    if any(r.get("failed") for r in results.values()):
        print("INCOMPLETE: at least one group below threshold")
        return

    # ---------------- 全局 9 类混淆矩阵(复用各组最终特征, 重新采样) --------
    all_names, all_feats = [], []
    for g in GROUPS:
        for k in results[g]["names"]:
            all_names.append(k)
            all_feats.append(results[g]["feats"][k])
    tr_x = np.vstack([make_patches(f, N_TRAIN, 900000 + i)
                      for i, f in enumerate(all_feats)])
    te_x = np.vstack([make_patches(f, N_TEST, 950000 + i)
                      for i, f in enumerate(all_feats)])
    y_tr = np.repeat(range(9), N_TRAIN)
    y_te = np.repeat(range(9), N_TEST)
    acc9, pred9 = nc_accuracy(tr_x, y_tr, te_x, y_te)
    cm9 = np.zeros((9, 9), int)
    for t, q in zip(y_te, pred9):
        cm9[t, q] += 1
    print(f"global 9-class acc={acc9:.3f}")

    np.savetxt(os.path.join(OUT, "confusion_matrix.csv"), cm9, fmt="%d",
               delimiter=",", header=",".join(all_names))

    fig, ax = plt.subplots(figsize=(6.4, 5.6))
    im = ax.imshow(cm9 / N_TEST, cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(9), all_names, rotation=45, ha="right")
    ax.set_yticks(range(9), all_names)
    for i in range(9):
        for j in range(9):
            ax.text(j, i, f"{cm9[i, j] / N_TEST:.2f}", ha="center",
                    va="center", fontsize=7,
                    color="white" if cm9[i, j] / N_TEST < 0.5 else "black")
    ax.set_xlabel("predicted")
    ax.set_ylabel("true")
    ax.set_title(f"9-class confusion (row-normalized), acc={acc9:.3f}")
    fig.colorbar(im)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "confusion_matrix.png"), dpi=130)
    plt.close(fig)

    # ---------------- 纹理拼图: 9 岩性 x 3 样本 ----------------
    fig, axes = plt.subplots(9, 3, figsize=(6, 16))
    for i, (nm, f) in enumerate(zip(all_names, all_feats)):
        for j in range(3):
            axes[i, j].imshow(texture(np.array(f), 800000 + i * 10 + j),
                              cmap="gray", vmin=0, vmax=1)
            axes[i, j].set_xticks([])
            axes[i, j].set_yticks([])
            if j == 0:
                axes[i, j].set_ylabel(nm, rotation=0, ha="right", va="center",
                                      fontsize=9)
    fig.suptitle("lithology texture exemplars (3 per class)")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "collage.png"), dpi=130)
    plt.close(fig)

    summary = dict(d0=d0, p=p, r2=r2, threshold=THRESHOLD,
                   global_acc=float(acc9), class_names=all_names,
                   cm9=cm9.tolist(),
                   groups={g: {k: v for k, v in r.items() if k != "test"}
                           for g, r in results.items()})
    with open(os.path.join(OUT, "results.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print("done ->", OUT)


if __name__ == "__main__":
    main()
