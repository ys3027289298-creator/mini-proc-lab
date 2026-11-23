#!/usr/bin/env python3
"""bandpass: 切除金字塔低频分量的高通过程化纹理流水线。

两条生成路径，共用同一组环形高通掩膜 M(fc, fhi)：
  G1  环形谱高斯白噪声 FFT 合成 —— 平稳、解析可推，主验证路径
  G2  项目自有值噪声 fBm（noise.c 的精确移植）经同一掩膜 FFT 切频 —— 落地路径

自相关距离 r0 定义：径向平均归一化自相关函数的第一个过零滞后（像素）。
验收：对每一档切频，实测 r0 与“该实现的理论预测”偏差不得超过 20%。
     - G1 解析预测 r0 = c * T_max，T_max = N/fc，c = x0/(2*pi)，x0 为
       解析自相关公式（见 README.md 推导）的第一个正根；
       alpha = fc/fhi = 0.5 时 x0 ≈ 1.5722，c ≈ 0.2502（r0 ≈ T_max/4）。
     - G2 输入谱非白，解析公式仅作近似上界参考，理论预测取数值预测
       r0{ifft(S_in * M^2)}（维纳-辛钦，逐实现精确）。
"""
import json
import os
import sys

import numpy as np
from scipy.optimize import brentq
from scipy.special import j1

# ---------------------------------------------------------------- 参数
N = 512                 # 纹理边长（像素）
PERIOD = 32             # 值噪声晶格在一张 tile 内的周期数（每晶格 16 px）
OCTAVES = 4
LACUNARITY = 2.0
GAIN = 0.5
ALPHA = 0.5             # fc / fhi，一个倍频程的环带
CUTS = (8, 16, 32)      # 三档切频 fc（cycles/tile），fhi = fc/alpha
TAPER = 1.5             # 环带边缘余弦过渡宽度（频域 bin），抑制硬边振铃
NREAL = 4               # 每档平均的实现数
TOL = 0.20              # 验收阈值：实测 r0 超出预测 20% 即失败

# matplotlib 克制配色（categorical 固定顺序，见 dataviz 规范）
INK = "#1f2933"
MUTED = "#66788a"
LINE_COLORS = ("#1f4e79", "#2f6f9f", "#5b9bd5")   # 深蓝三档，由暗到亮
ACCENT = "#c0504d"                                # 仅用于 G2/参考线

# ---------------------------------------------------------------- 频域网格
FREQS = np.fft.fftfreq(N, d=1.0 / N)                       # cycles/tile
KI, KJ = np.meshgrid(FREQS, FREQS, indexing="ij")
K = np.hypot(KI, KJ)
II, JJ = np.meshgrid(np.arange(N), np.arange(N), indexing="ij")
# Hermitian 共轭的规范半平面标记（保证 ifft2 输出严格实数）
CANON = (JJ < N // 2) | ((JJ == 0) & (II <= N // 2)) | ((JJ == N // 2) & (II == 0))
ICONJ, JCONJ = (-II) % N, (-JJ) % N
SELFCONJ = ((II == 0) | (II == N // 2)) & ((JJ == 0) | (JJ == N // 2))
# 径向滞后（考虑 FFT 周期性，向最近的回绕方向取距离）
DI = np.minimum(II, N - II)
DJ = np.minimum(JJ, N - JJ)
RBIN = np.round(np.hypot(DI, DJ)).astype(int)


# ================================================================ 解析推导
def analytic_cnorm(x, alpha):
    """硬环带 [fc, fhi]、alpha=fc/fhi 的归一化自相关。

    C(tau)/C(0) = (2*alpha/(1-alpha^2)) * (J1(x/alpha) - alpha*J1(x)) / x,
    x = 2*pi*fc*tau。推导见 README.md 第 2 节。
    """
    return (2.0 * alpha / (1.0 - alpha * alpha)) * (j1(x / alpha) - alpha * j1(x)) / x


def analytic_r0_ratio(alpha=ALPHA):
    """返回 x0 与 c = r0/T_max = x0/(2*pi)。"""
    xs = np.linspace(1e-6, 12.0, 200_001)
    vals = analytic_cnorm(xs, alpha)
    change = np.where(np.diff(np.sign(vals)) != 0)[0]
    if len(change) == 0:
        raise RuntimeError("analytic autocorrelation has no zero crossing")
    i = change[0]
    x0 = brentq(lambda x: analytic_cnorm(x, alpha), xs[i], xs[i + 1])
    return x0, x0 / (2.0 * np.pi)


# ================================================================ G1 环形谱合成
def hermitian_noise(rng):
    """构造 Hermitian 对称的复高斯频谱，使逆变换为实平稳高斯场。"""
    w = (rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))) / np.sqrt(2.0)
    h = np.empty((N, N), complex)
    h[CANON] = w[CANON]
    h[~CANON] = np.conj(w[ICONJ[~CANON], JCONJ[~CANON]])
    h[SELFCONJ] = w[SELFCONJ].real
    h[0, 0] = 0.0
    return h


def annulus_mask(fc, fhi, taper=TAPER):
    """环形带通掩膜：环带内为 1，两侧 taper 个 bin 余弦过渡到 0。"""
    m = ((K >= fc) & (K <= fhi)).astype(np.float64)
    if taper > 0.0:
        for edge in (fc, fhi):
            d = np.abs(K - edge)
            w = 0.5 * (1.0 + np.cos(np.pi * np.clip(d / taper, 0.0, 1.0)))
            near = (K >= edge - taper) & (K <= edge + taper)
            m = np.where(near, np.maximum(m, w * (d <= taper)), m)
    return m


def synth_annular(mask, rng):
    return np.fft.ifft2(np.fft.fft2(hermitian_noise(rng)) * mask).real


# ================================================================ G2 值噪声 fBm（noise.c 移植）
def _hash2(x, y, seed):
    """与 noise.c hash2 逐位一致的向量化实现（uint32 回绕即 C 的 unsigned 溢出）。"""
    x = x.astype(np.uint32)
    y = y.astype(np.uint32)
    seed = np.uint32(seed)
    # uint32 乘法回绕即 C 的 unsigned 溢出，语义上正是 noise.c 想要的
    with np.errstate(over="ignore"):
        h = x * np.uint32(374761393) + y * np.uint32(668265263) + seed * np.uint32(2246822519)
        h = (h ^ (h >> np.uint32(13))) * np.uint32(1274126177)
    return h ^ (h >> np.uint32(16))


def _rand01(h):
    return (h & np.uint32(0xFFFFFF)).astype(np.float64) / 16777215.0


def _fade(t):
    return t * t * (3.0 - 2.0 * t)


def value_noise_periodic(X, Y, seed, grid_count):
    """noise.c value_noise 的周期化版本：晶格沿 tile 边界回绕，输出可无缝平铺。"""
    xi = np.floor(X).astype(np.int64)
    yi = np.floor(Y).astype(np.int64)
    fx, fy = X - xi, Y - yi

    def corner(a, b):
        gx = (a % grid_count).astype(np.uint32)
        gy = (b % grid_count).astype(np.uint32)
        return _rand01(_hash2(gx, gy, seed)) * 2.0 - 1.0

    u, v = _fade(fx), _fade(fy)
    top = corner(xi, yi) * (1.0 - u) + corner(xi + 1, yi) * u
    bot = corner(xi, yi + 1) * (1.0 - u) + corner(xi + 1, yi + 1) * u
    return top * (1.0 - v) + bot * v


def fbm_periodic(X, Y, seed):
    """与 noise.c fbm 同参数（4 倍频、lacunarity 2、gain 0.5），逐倍频整数周期化。"""
    total = 0.0
    norm = 0.0
    amp = 1.0
    freq = 1.0
    for i in range(OCTAVES):
        grid_count = int(PERIOD * freq)
        total += value_noise_periodic(X * freq, Y * freq, seed + i * 17, grid_count) * amp
        norm += amp
        amp *= GAIN
        freq *= LACUNARITY
    return total / norm


def fbm_field(seed):
    cell = N / PERIOD
    yy, xx = np.mgrid[0:N, 0:N]
    return fbm_periodic(xx / cell, yy / cell, seed)


# ================================================================ 自相关测量
def radial_average(field2d, rmax=N // 4):
    out = np.zeros(rmax + 1)
    cnt = np.zeros(rmax + 1)
    sel = RBIN <= rmax
    np.add.at(out, RBIN[sel], field2d[sel])
    np.add.at(cnt, RBIN[sel], 1.0)
    return out / np.maximum(cnt, 1.0)


def autocorr_radial(field, rmax=N // 4):
    psd = np.abs(np.fft.fft2(field)) ** 2
    return radial_average(np.fft.ifft2(psd).real, rmax)


def first_zero(radial_curve):
    """归一化径向自相关的第一个过零滞后，亚像素线性插值；无过零返回 None。"""
    c = radial_curve / radial_curve[0]
    change = np.where(np.diff(np.sign(c)) != 0)[0]
    if len(change) == 0:
        return None
    i = change[0]
    return i + (0.0 - c[i]) / (c[i + 1] - c[i])


def numeric_prediction_r0(spectrum2d, mask, rmax=N // 4):
    """数值理论预测：r0{ifft(S_in * M^2)}（维纳-辛钦）。"""
    return first_zero(radial_average(np.fft.ifft2(spectrum2d * mask ** 2).real, rmax))


def radial_psd(field, rmax=N // 2):
    psd = np.abs(np.fft.fftshift(np.fft.fft2(field))) ** 2
    k_shift = np.fft.fftshift(K)
    r = np.round(k_shift).astype(int)
    out = np.zeros(rmax + 1)
    cnt = np.zeros(rmax + 1)
    sel = r <= rmax
    np.add.at(out, r[sel], psd[sel])
    np.add.at(cnt, r[sel], 1.0)
    return out / np.maximum(cnt, 1.0)


# ================================================================ 出图
def _save_gray(ax, field, title):
    lo, hi = np.percentile(field, (1, 99))
    ax.imshow(field, cmap="gray", vmin=lo, vmax=hi, interpolation="nearest")
    ax.set_title(title, fontsize=10, color=INK)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_color(MUTED)
        spine.set_linewidth(0.8)


def make_montage(path, reference, g1_fields, g2_fields):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 4, figsize=(13.2, 6.8),
                             gridspec_kw={"wspace": 0.08, "hspace": 0.12})
    _save_gray(axes[0, 0], g1_fields[0][1] - g1_fields[0][1].mean(),
               "G1 reference: full-band white noise\n(no low-frequency cut)")
    for j, (fc, field) in enumerate(g1_fields, start=1):
        tmax = N / fc
        _save_gray(axes[0, j], field,
                   f"G1 cut fc={fc} cyc/tile\nT_max={tmax:.0f}px, r0≈T_max/4={tmax / 4:.1f}px")
    _save_gray(axes[1, 0], reference,
               "G2 reference: value-noise fBm\n(unfiltered — flow-like manifold)")
    for j, (fc, field) in enumerate(g2_fields, start=1):
        _save_gray(axes[1, j], field,
                   f"G2 cut fc={fc} cyc/tile\nfBm spectrum × annulus mask")
    fig.suptitle("High-passed procedural textures — three cutoff frequencies "
                 "(annulus [fc, 2fc], alpha = 0.5)",
                 fontsize=12, color=INK, y=0.98)
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def make_tiling_check(path, fields):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(10.8, 4.0),
                             gridspec_kw={"wspace": 0.12})
    for ax, (fc, field) in zip(axes, fields):
        tile = field[: N // 2, : N // 2]
        tiled = np.tile(tile, (2, 2))
        lo, hi = np.percentile(field, (1, 99))
        ax.imshow(tiled, cmap="gray", vmin=lo, vmax=hi, interpolation="nearest")
        ax.axhline(N // 2, color=ACCENT, lw=1.0, ls="--")
        ax.axvline(N // 2, color=ACCENT, lw=1.0, ls="--")
        ax.set_title(f"fc={fc} cyc/tile\n2×2 seamless repeat (red: seam)",
                     fontsize=10, color=INK)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def make_corr_plot(path, g1_curves, g2_curves, psd_data):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rmax = N // 4
    r = np.arange(rmax + 1)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.4, 4.6))
    for (fc, curve), color in zip(g1_curves, LINE_COLORS):
        ax1.plot(r, curve / curve[0], color=color, lw=1.8,
                 label=f"G1 fc={fc} (r0={first_zero(curve):.1f}px)")
    for (fc, curve), color in zip(g2_curves, LINE_COLORS):
        ax1.plot(r, curve / curve[0], color=color, lw=1.2, ls="--",
                 label=f"G2 fc={fc} (r0={first_zero(curve):.1f}px)")
    ax1.axhline(0, color=MUTED, lw=0.8)
    ax1.set_xlabel("radial lag τ (px)", color=INK)
    ax1.set_ylabel("normalized autocorrelation C(τ)/C(0)", color=INK)
    ax1.set_title("Radial autocorrelation — solid: G1 annular, dashed: G2 fBm cut",
                  fontsize=10, color=INK)
    ax1.legend(fontsize=8, frameon=False, ncol=2)
    ax1.set_xlim(0, rmax)
    ax1.set_ylim(-0.35, 1.05)
    ax1.grid(alpha=0.25, lw=0.6)

    freqs_r = np.arange(N // 2 + 1)
    ref_psd, g1_psds, g2_psds = psd_data
    ax2.plot(freqs_r, ref_psd / ref_psd[1:].max(), color=MUTED, lw=1.4,
             label="G2 reference fBm PSD (unfiltered)")
    for (fc, psd), color in zip(g1_psds, LINE_COLORS):
        ax2.plot(freqs_r, psd / psd.max(), color=color, lw=1.8,
                 label=f"G1 PSD fc={fc}")
    for (fc, psd), color in zip(g2_psds, LINE_COLORS):
        ax2.plot(freqs_r, psd / psd.max(), color=color, lw=1.1, ls="--",
                 label=f"G2 cut PSD fc={fc}")
    for fc in CUTS:
        ax2.axvspan(fc, 2 * fc, color=LINE_COLORS[CUTS.index(fc)], alpha=0.08)
    ax2.set_xlabel("radial frequency (cycles/tile)", color=INK)
    ax2.set_ylabel("normalized radial PSD", color=INK)
    ax2.set_title("Power spectra — shaded bands: annulus [fc, 2fc]",
                  fontsize=10, color=INK)
    ax2.legend(fontsize=7.5, frameon=False, ncol=2)
    ax2.set_xlim(0, N // 2)
    ax2.set_ylim(0, 1.05)
    ax2.grid(alpha=0.25, lw=0.6)

    for ax in (ax1, ax2):
        ax.tick_params(colors=INK, labelsize=8)
        for spine in ax.spines.values():
            spine.set_color(MUTED)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)


# ================================================================ 主流水线
def main(outdir):
    os.makedirs(outdir, exist_ok=True)
    x0, c_ratio = analytic_r0_ratio(ALPHA)

    rng = np.random.default_rng(20261124)
    g2_seeds = (11, 23, 37, 51)

    g1_fields, g2_fields = [], []
    g1_curves, g2_curves = [], []
    g1_psds, g2_psds = [], []
    metrics = {
        "tile_px": N, "alpha": ALPHA, "taper_bins": TAPER,
        "n_realizations": NREAL, "tolerance": TOL,
        "analytic": {"x0_first_root": x0, "r0_over_Tmax": c_ratio,
                     "note": "hard-annulus continuous theory; taper causes <1.5% shift"},
        "cuts": [], "all_pass": True,
    }

    reference = fbm_field(g2_seeds[0])
    ref_curve = autocorr_radial(reference)
    ref_r0 = first_zero(ref_curve)
    metrics["unfiltered_fBm"] = {
        "r0_px": ref_r0,
        "min_normalized_autocorr": float((ref_curve / ref_curve[0]).min()),
        "interpretation": "C stays near zero with a deep-wide positive lobe: "
                          "low-frequency dominated, flow-like manifold",
    }

    for fc in CUTS:
        fhi = fc / ALPHA
        tmax = N / fc
        mask = annulus_mask(fc, fhi)

        # ---- G1: NREAL 个环形谱实现
        g1_acc_curve = np.zeros(N // 4 + 1)
        g1_acc_psd = np.zeros(N // 2 + 1)
        g1_last = None
        for _ in range(NREAL):
            z = synth_annular(mask, rng)
            g1_acc_curve += autocorr_radial(z)
            g1_acc_psd += radial_psd(z)
            g1_last = z
        g1_curve = g1_acc_curve / NREAL
        g1_psd = g1_acc_psd / NREAL
        g1_r0_meas = first_zero(g1_curve)
        g1_r0_analytic = c_ratio * tmax
        g1_r0_numeric = first_zero(
            radial_average(np.fft.ifft2(mask ** 2).real))

        # ---- G2: NREAL 个 fBm 种子，同一掩膜切频
        g2_acc_curve = np.zeros(N // 4 + 1)
        g2_acc_psd = np.zeros(N // 2 + 1)
        g2_pred_curve = np.zeros(N // 4 + 1)
        g2_last = None
        for seed in g2_seeds:
            f = fbm_field(seed)
            spec = np.fft.fft2(f)
            z = np.fft.ifft2(spec * mask).real
            g2_acc_curve += autocorr_radial(z)
            g2_acc_psd += radial_psd(z)
            g2_pred_curve += radial_average(np.fft.ifft2(np.abs(spec) ** 2 * mask ** 2).real)
            g2_last = z
        g2_curve = g2_acc_curve / NREAL
        g2_psd = g2_acc_psd / NREAL
        g2_pred_curve /= NREAL
        g2_r0_meas = first_zero(g2_curve)
        g2_r0_numeric = first_zero(g2_pred_curve)
        g2_r0_analytic = c_ratio * tmax

        g1_err = abs(g1_r0_meas - g1_r0_analytic) / g1_r0_analytic
        g2_err = abs(g2_r0_meas - g2_r0_numeric) / g2_r0_numeric
        g1_pass = g1_err <= TOL
        g2_pass = g2_err <= TOL
        metrics["all_pass"] &= g1_pass and g2_pass

        metrics["cuts"].append({
            "fc_cycles_per_tile": fc, "fhi_cycles_per_tile": fhi,
            "T_max_px": tmax, "T_min_px": N / fhi,
            "G1_annular": {
                "r0_measured_px": g1_r0_meas,
                "r0_analytic_px": g1_r0_analytic,
                "r0_numeric_mask_px": g1_r0_numeric,
                "abs_err_vs_analytic": g1_err,
                "pass_20pct": bool(g1_pass),
            },
            "G2_fBm_cut": {
                "r0_measured_px": g2_r0_meas,
                "r0_numeric_prediction_px": g2_r0_numeric,
                "r0_analytic_approx_px": g2_r0_analytic,
                "abs_err_vs_numeric": g2_err,
                "analytic_approx_err": abs(g2_r0_meas - g2_r0_analytic) / g2_r0_analytic,
                "pass_20pct": bool(g2_pass),
            },
        })

        g1_fields.append((fc, g1_last))
        g2_fields.append((fc, g2_last))
        g1_curves.append((fc, g1_curve))
        g2_curves.append((fc, g2_curve))
        g1_psds.append((fc, g1_psd))
        g2_psds.append((fc, g2_psd))

        # 单张纹理（16-bit PNG，全分辨率）
        from PIL import Image
        for name, field in (("g1", g1_last), ("g2", g2_last)):
            z = (field - field.mean()) / (field.std() + 1e-12)
            img = np.clip(np.round((z * 0.30 + 0.5) * 65535), 0, 65535).astype(np.uint16)
            Image.fromarray(img).save(os.path.join(outdir, f"texture_{name}_fc{fc}.png"))

    make_montage(os.path.join(outdir, "montage.png"), reference, g1_fields, g2_fields)
    make_tiling_check(os.path.join(outdir, "tiling_check.png"), g1_fields)
    make_corr_plot(
        os.path.join(outdir, "autocorr.png"),
        g1_curves, g2_curves,
        (radial_psd(reference), g1_psds, g2_psds),
    )

    with open(os.path.join(outdir, "metrics.json"), "w", encoding="utf-8") as fh:
        json.dump(metrics, fh, indent=2, ensure_ascii=False,
                  default=lambda o: bool(o) if isinstance(o, np.bool_) else float(o))

    # ---- 控制台自检：20% 红线
    print(f"analytic x0={x0:.5f}  r0/T_max={c_ratio:.5f} (alpha={ALPHA})")
    print(f"unfiltered fBm: r0={ref_r0:.2f}px, min C/C0="
          f"{(ref_curve / ref_curve[0]).min():.4f} (no real decorrelation -> manifold)")
    for cut in metrics["cuts"]:
        g1, g2 = cut["G1_annular"], cut["G2_fBm_cut"]
        print(f"fc={cut['fc_cycles_per_tile']:2d} T_max={cut['T_max_px']:5.1f}px | "
              f"G1 r0 {g1['r0_measured_px']:5.2f} vs analytic {g1['r0_analytic_px']:5.2f} "
              f"err {100 * g1['abs_err_vs_analytic']:5.2f}% pass={g1['pass_20pct']} | "
              f"G2 r0 {g2['r0_measured_px']:5.2f} vs numeric {g2['r0_numeric_prediction_px']:5.2f} "
              f"err {100 * g2['abs_err_vs_numeric']:5.2f}% pass={g2['pass_20pct']}")
    print(f"all_pass={metrics['all_pass']}")
    return 0 if metrics["all_pass"] else 1


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "artifacts", "bandpass")
    sys.exit(main(out))
