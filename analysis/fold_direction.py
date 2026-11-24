#!/usr/bin/env python3 -I
# -*- coding: utf-8 -*-
"""fold_direction: 褶皱方向场 vs 山脊走向 / 噪声裂隙方向 对齐分析.

流程:
 1. 用 noise/ 模块同款算法 (hash2/value_noise/fbm 的 numpy 移植) 造两座不同走向的山脉;
    每座植入一条斜交谷地和一处陡崖, 作为 >20° 偏角的已知成因对照.
 2. 逐点求 Hessian 主曲率方向, 取 |特征值| 较小的特征向量作为褶皱(沿脊)方向.
 3. 在山脊线锚点播种, 按海拔从高到低向下外推方向场 (倍角向量传播).
 4. 每座山按三个海拔切片, 统计实测方向 vs 外推方向的偏角均值/标准差.
 5. 标出偏角 >20° 的格点, 按坡度/曲率符号判定是山谷还是陡崖导致.
 6. 山脚处将外推方向场与噪声纹理的裂隙方向 (噪声场结构张量沿等值线方向)
    做对齐检验, 并按海拔逐层对比判断哪个模块先偏.
"""
import numpy as np

# ---------------------------------------------------------------- noise 模块移植
# 与 noise/src/noise.c 完全相同的确定性哈希 + 值噪声 + fBm.
U32 = np.uint32


def hash2(x, y, seed):
    x = np.asarray(x, dtype=U32); y = np.asarray(y, dtype=U32)
    with np.errstate(over="ignore"):          # uint32 回绕是有意为之, 与 C 一致
        h = x * U32(374761393) + y * U32(668265263) + U32(seed) * U32(2246822519)
        h = (h ^ (h >> U32(13))) * U32(1274126177)
        return h ^ (h >> U32(16))


def rand01(h):
    return (h & U32(0xFFFFFF)).astype(np.float64) / 16777215.0


def value_noise(px, py, seed):
    xi = np.floor(px).astype(np.int64); yi = np.floor(py).astype(np.int64)
    fx = px - xi; fy = py - yi
    v00 = rand01(hash2(xi, yi, seed));     v10 = rand01(hash2(xi + 1, yi, seed))
    v01 = rand01(hash2(xi, yi + 1, seed)); v11 = rand01(hash2(xi + 1, yi + 1, seed))
    ux = fx * fx * (3 - 2 * fx); uy = fy * fy * (3 - 2 * fy)
    a = v00 + (v10 - v00) * ux
    b = v01 + (v11 - v01) * ux
    return (a + (b - a) * uy) * 2.0 - 1.0


def fbm(px, py, octaves=5, lac=2.0, gain=0.5, seed=11):
    s = np.zeros_like(px); amp = 1.0; norm = 0.0; freq = 1.0
    for i in range(octaves):
        s += value_noise(px * freq, py * freq, seed + i * 17) * amp
        norm += amp; amp *= gain; freq *= lac
    return s / norm


# ---------------------------------------------------------------- 工具
def gauss_blur(img, sigma):
    r = max(1, int(3 * sigma))
    k = np.exp(-0.5 * (np.arange(-r, r + 1) / sigma) ** 2); k /= k.sum()
    tmp = np.apply_along_axis(lambda r_: np.convolve(r_, k, mode="same"), 0, img)
    return np.apply_along_axis(lambda r_: np.convolve(r_, k, mode="same"), 1, tmp)


def ang_diff_pi(a, b):
    """无向方向夹角, 返回 [0, pi/2]."""
    d = np.abs((a - b) % np.pi)
    return np.minimum(d, np.pi - d)


def circ_std_pi(theta):
    """无向方向(周期 pi)的圆标准差, 返回角度(度)."""
    z = np.exp(2j * theta)
    R = np.abs(z.mean())
    if R >= 1.0:
        return 0.0
    return np.degrees(0.5 * np.sqrt(-2.0 * np.log(max(R, 1e-12))))


# ---------------------------------------------------------------- 造山
N = 256
H0 = 100.0
RIDGE_W = 28.0
VALLEY_DEPTH, VALLEY_W = 25.0, 6.0
CLIFF_DROP, CLIFF_W = 18.0, 3.0
NOISE_AMP = 6.0


def make_range(trend_deg, seed):
    """返回 (height, noise_only, feature_lines). 山脊过中心, 走向 trend_deg."""
    yy, xx = np.mgrid[0:N, 0:N].astype(np.float64)
    cx = cy = N / 2.0
    u = np.array([np.cos(np.radians(trend_deg)), np.sin(np.radians(trend_deg))])
    n = np.array([-u[1], u[0]])
    rx, ry = xx - cx, yy - cy
    d_perp = rx * n[0] + ry * n[1]          # 距山脊线的垂直距离
    s_along = rx * u[0] + ry * u[1]         # 沿山脊线坐标

    noise = fbm(xx * 0.03, yy * 0.03, octaves=5, seed=seed)
    mod = 0.85 + 0.30 * fbm(s_along * 0.02, np.zeros_like(s_along), octaves=3, seed=seed + 101)
    h = H0 * np.exp(-((d_perp / RIDGE_W) ** 2)) * mod + NOISE_AMP * noise

    # 斜交谷地: 走向 trend+60°, 过 (cx+20, cy-15)
    uv = np.array([np.cos(np.radians(trend_deg + 60)), np.sin(np.radians(trend_deg + 60))])
    nv = np.array([-uv[1], uv[0]])
    pvx, pvy = cx + 20 - cx, cy - 15 - cy
    dv = (rx - pvx) * nv[0] + (ry - pvy) * nv[1]
    sv = (rx - pvx) * uv[0] + (ry - pvy) * uv[1]
    valley = -VALLEY_DEPTH * np.exp(-((dv / VALLEY_W) ** 2)) * (np.abs(sv) < 70)
    h += valley

    # 陡崖: 走向 trend-40°, 过 (cx-30, cy+25), 单侧台阶
    uc = np.array([np.cos(np.radians(trend_deg - 40)), np.sin(np.radians(trend_deg - 40))])
    nc = np.array([-uc[1], uc[0]])
    pcx, pcy = -30.0, 25.0
    dc = (rx - pcx) * nc[0] + (ry - pcy) * nc[1]
    sc = (rx - pcx) * uc[0] + (ry - pcy) * uc[1]
    step = 1.0 / (1.0 + np.exp(-dc / CLIFF_W))
    cliff = -CLIFF_DROP * step * np.exp(-((sc / 45.0) ** 2))
    h += cliff

    feats = {"valley": (pvx + cx, pvy + cy, uv), "cliff": (pcx + cx, pcy + cy, uc)}
    return h, noise, feats, u


# ---------------------------------------------------------------- 主曲率方向
def fold_direction(h, sigma=1.0):
    """逐点主曲率方向: Hessian 的 |特征值| 较小者的特征向量角度 (mod pi)."""
    hs = gauss_blur(h, sigma)
    hy, hx = np.gradient(hs)
    hxx = np.gradient(hx, axis=1); hxy = np.gradient(hx, axis=0)
    hyy = np.gradient(hy, axis=0)
    # 对称 2x2 特征分解
    tr = hxx + hyy
    det = hxx * hyy - hxy * hxy
    disc = np.sqrt(np.maximum((tr / 2) ** 2 - det, 0.0))
    l1 = tr / 2 + disc; l2 = tr / 2 - disc          # l1 >= l2
    # 沿褶皱方向 = |特征值| 较小者
    small_is_l1 = np.abs(l1) < np.abs(l2)
    # l1 的特征向量: (hxy, l1-hxx); l2 的: (hxy, l2-hxx)
    vx = np.where(small_is_l1, hxy, hxy)
    vy = np.where(small_is_l1, l1 - hxx, l2 - hxx)
    theta = np.arctan2(vy, vx) % np.pi
    # |特征值| 较大者(跨褶皱方向曲率), 用于谷/脊判定
    lam_big = np.where(small_is_l1, l2, l1)
    return theta, lam_big, np.hypot(hx, hy)


# ---------------------------------------------------------------- 方向场外推
def propagate(theta, anchors, mask, h):
    """从山脊锚点按海拔降序把方向场(倍角表示)外推到整个 mask."""
    z = np.zeros((N, N, 2))
    z[..., 0] = np.cos(2 * theta); z[..., 1] = np.sin(2 * theta)
    out = np.full((N, N), np.nan)
    done = anchors.copy()
    out[anchors] = theta[anchors]
    zc = z.copy()
    order = np.argsort(-h, axis=None)
    hm = mask.ravel()
    nbr = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
    for idx in order:
        if not hm[idx]:
            continue
        y, x = divmod(idx, N)
        if done[y, x]:
            continue
        acc = np.zeros(2)
        for dy, dx in nbr:
            ny, nx = y + dy, x + dx
            if 0 <= ny < N and 0 <= nx < N and done[ny, nx]:
                acc += zc[ny, nx]
        if np.hypot(*acc) < 1e-9:
            continue
        acc /= np.hypot(*acc)
        zc[y, x] = acc
        out[y, x] = 0.5 * np.arctan2(acc[1], acc[0]) % np.pi
        done[y, x] = True
    return out


# ---------------------------------------------------------------- 裂隙方向 (噪声模块)
def fissure_direction(noise):
    """噪声场结构张量的沿等值线方向 = 裂隙走向."""
    g = gauss_blur(noise, 1.0)
    gy, gx = np.gradient(g)
    jxx = gauss_blur(gx * gx, 3.0); jyy = gauss_blur(gy * gy, 3.0); jxy = gauss_blur(gx * gy, 3.0)
    tr = jxx + jyy; det = jxx * jyy - jxy * jxy
    disc = np.sqrt(np.maximum((tr / 2) ** 2 - det, 0.0))
    lmin = tr / 2 - disc
    vx = jxy; vy = lmin - jxx
    return np.arctan2(vy, vx) % np.pi


# ---------------------------------------------------------------- 分析主流程
def analyze(name, trend_deg, seed):
    h, noise, feats, u = make_range(trend_deg, seed)
    hmax = h.max()
    trend = np.radians(trend_deg) % np.pi

    theta, lam_big, slope = fold_direction(h, sigma=2.0)
    mask = h > 0.05 * hmax
    anchors = mask & (h > 0.70 * hmax) & (lam_big < 0)
    if anchors.sum() < 20:                       # 锚点太少就放宽凸性条件
        anchors = mask & (h > 0.75 * hmax)
    ext = propagate(theta, anchors, mask, h)

    dev = np.degrees(ang_diff_pi(theta, ext))    # 实测 vs 外推 偏角

    bands = [("高海拔", 0.66, 1.01), ("中海拔", 0.33, 0.66), ("低海拔/山脚", 0.05, 0.33)]
    print(f"\n===== {name} (走向 {trend_deg}°, hmax={hmax:.1f}, 锚点 {anchors.sum()} 点) =====")
    print(f"{'切片':<10}{'格点数':>7}{'偏角均值':>9}{'偏角std':>9}"
          f"{'实测vs走向std':>13}{'外推vs走向std':>13}")
    slice_stats = []
    for label, lo, hi in bands:
        m = mask & (h >= lo * hmax) & (h < hi * hmax) & ~np.isnan(ext)
        if m.sum() == 0:
            print(f"{label:<10}  (无格点)"); slice_stats.append(None); continue
        print(f"{label:<10}{m.sum():>7}{dev[m].mean():>8.2f}°{dev[m].std():>8.2f}°"
              f"{circ_std_pi(theta[m] - trend):>12.2f}°{circ_std_pi(ext[m] - trend):>12.2f}°")
        slice_stats.append(m)

    # ---- >20° 格点标记与归因 (两个曲率尺度对照: σ=1 逐点原始, σ=2 主流程)
    yy, xx = np.mgrid[0:N, 0:N]
    feat_zone = {}
    for fname, (fx, fy, fu) in feats.items():
        d = np.abs((xx - fx) * (-fu[1]) + (yy - fy) * fu[0])
        feat_zone[fname] = mask & (d < 8.0)
    p90_slope = np.percentile(slope[mask], 90)
    print(f"\n偏角>20° 格点归因 (坡度阈值 p90={p90_slope:.2f}):")
    for sigma, tag in ((1.0, "σ=1 逐点"), (2.0, "σ=2 主流程")):
        th_s, lam_s, _ = fold_direction(h, sigma=sigma)
        dev_s = np.degrees(ang_diff_pi(th_s, ext))
        flag = mask & ~np.isnan(ext) & (dev_s > 20.0)
        n_flag = flag.sum()
        n_cliff = (flag & (slope > p90_slope)).sum()
        n_valley = (flag & (slope <= p90_slope) & (lam_s > 0)).sum()
        n_other = n_flag - n_cliff - n_valley
        rate_all = 100.0 * n_flag / mask.sum()
        print(f"  [{tag}] flagged {n_flag} 个 (占掩膜 {rate_all:.1f}%) | "
              f"陡崖型 {n_cliff}, 山谷型 {n_valley}, 噪声/其他 {n_other}")
        for fname, zone in feat_zone.items():
            r = 100.0 * (flag & zone).sum() / max(zone.sum(), 1)
            print(f"      植入{'谷地' if fname=='valley' else '陡崖'}线两侧8格内标记率 "
                  f"{r:.1f}% (全局 {rate_all:.1f}%)")
    # 主流程 flag 供后续使用
    flag = mask & ~np.isnan(ext) & (dev > 20.0)

    # ---- 山脚裂隙对齐检验
    foot = slice_stats[2]
    fiss = fissure_direction(noise)
    ext_trend = np.degrees(ang_diff_pi(ext, trend))         # 外推 vs 走向
    align = np.degrees(ang_diff_pi(ext, fiss))
    fiss_trend = np.degrees(ang_diff_pi(fiss, trend))
    print(f"\n山脚对齐检验 (外推方向场 vs 噪声裂隙方向):")
    print(f"  夹角均值 {align[foot].mean():.2f}°, std {align[foot].std():.2f}°, "
          f">20° 占比 {100.0 * (align[foot] > 20).mean():.1f}%")
    print(f"  噪声裂隙 vs 山脊走向: 均值 {fiss_trend[foot].mean():.2f}° "
          f"(各海拔: " + ", ".join(
              f"{fiss_trend[m].mean():.1f}°" if m is not None else "-"
              for m in slice_stats) + ")")
    print(f"  外推场  vs 山脊走向: 均值 " + ", ".join(
        f"{ext_trend[m].mean():.1f}°" if m is not None else "-" for m in slice_stats))

    # ---- 谁先偏: 哪个模块在更高海拔就开始偏离山脊走向
    first_dev = {}
    for mod, arr in (("terrain(外推场)", ext_trend), ("noise(裂隙)", fiss_trend)):
        first_dev[mod] = None
        for (label, lo, hi), m in zip(bands, slice_stats):
            if m is not None and arr[m].mean() > 20.0:
                first_dev[mod] = label
                break
    print(f"  谁先偏(均值首超20°的海拔层): " +
          ", ".join(f"{k}={v or '三层内未超'}" for k, v in first_dev.items()))
    return dev, flag


if __name__ == "__main__":
    np.set_printoptions(suppress=True)
    analyze("山脉A", trend_deg=25, seed=11)
    analyze("山脉B", trend_deg=100, seed=47)
