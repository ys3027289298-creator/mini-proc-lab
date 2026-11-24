#!/usr/bin/env python3
"""Render albedo/out/*.ppm to PNG and draw the report figures.

Run from albedo/:  python3 tools/figures.py
Inputs  : out/*.csv, out/*.ppm    (produced by `make build`)
Outputs : out/*.png
"""
import csv
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image

OUT = "out"
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
NEUTRAL = "#c9c8c1"

plt.rcParams.update({
    "font.family": ["Noto Sans CJK SC", "DejaVu Sans"],
    "font.size": 9,
    "axes.facecolor": SURFACE,
    "figure.facecolor": SURFACE,
    "axes.edgecolor": AXIS,
    "axes.labelcolor": INK2,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "text.color": INK,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "grid.linestyle": "-",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 150,
})


def read(name):
    with open(os.path.join(OUT, name + ".csv"), newline="") as f:
        return list(csv.DictReader(f))


def fget(r, k):
    return float(r[k])


def ppm(name):
    a = np.asarray(Image.open(os.path.join(OUT, name + ".ppm"))).astype(np.float64) / 255.0
    return a


# ---------------------------------------------------------------- image tiles
def save_pngs():
    for f in sorted(os.listdir(OUT)):
        if f.endswith(".ppm"):
            Image.open(os.path.join(OUT, f)).save(os.path.join(OUT, f[:-4] + ".png"))


def crop_rgb(a, x0, y0, s):
    return a[y0:y0 + s, x0:x0 + s, :]


def main():
    save_pngs()
    os.makedirs(os.path.join(OUT, "fig"), exist_ok=True)

    weather = read("weather")
    stab = {r["weather"]: r for r in read("stability")}
    lstat = read("layer_stats")
    trunc = read("truncation")
    ksweep = read("knob_independence")
    rsweep = read("resplit_check")
    scale_dec = read("scale_decision")
    scale_scan = read("scale_scan")
    cross = read("cross_terms")
    inter = read("interference")
    merged = read("merged_scheme")
    ssweep = read("structure_sweep")
    msweep = read("micro_sweep")
    corr = read("correction")

    # ============================================================= Figure 1
    # Scene, layers, clip masks, weather renders.
    alb = ppm("albedo")
    st = ppm("structure")
    mic = ppm("micro")
    clip_s = ppm("structure_clip")
    clip_m = ppm("micro_clip")
    wx = {w: ppm("weather_" + w) for w in ("clear", "overcast", "dusk")}

    # panel order and crop windows (x0, y0, size)
    run = [(10, 30, 280), (150, 120, 280), (230, 300, 280)]
    fig, axes = plt.subplots(4, 3, figsize=(9.6, 11.4))
    fig.patch.set_facecolor(SURFACE)
    for r in range(4):
        for c in range(3):
            ax = axes[r, c]
            ax.set_xticks([])
            ax.set_yticks([])
            ax.grid(False)
            for sp in ax.spines.values():
                sp.set_visible(False)
    for c, (x0, y0, s) in enumerate(run):
        axes[0, c].imshow(crop_rgb(alb, x0, y0, s), interpolation="nearest")
        axes[1, c].imshow(crop_rgb(st, x0, y0, s), interpolation="nearest")
        mm = crop_rgb(mic, x0, y0, s)
        mm = np.clip(0.5 + (mm - 0.5) * 2.2, 0, 1)
        axes[2, c].imshow(mm, interpolation="nearest")
        axes[3, c].imshow(crop_rgb(wx["clear"], x0, y0, s), interpolation="nearest")
    rows = ["合成 albedo\nS·e^m", "结构层 S\nlog 反射率 → 显示", "微噪声层 m\n±13 %（对比 ×2.2）",
            "晴天渲染\n（用于偏差测量）"]
    for r in range(4):
        axes[r, 0].set_ylabel(rows[r], fontsize=9, color=INK2, rotation=0,
                              ha="right", va="center", labelpad=18)
    axes[0, 0].set_title("样本 1", fontsize=9, color=INK2)
    axes[0, 1].set_title("样本 2", fontsize=9, color=INK2)
    axes[0, 2].set_title("样本 3", fontsize=9, color=INK2)
    fig.suptitle("图 1  三层分解：结构层承载材质与明暗，微噪声层只承载细粒（三处局部裁剪）",
                 fontsize=11.5, x=0.5, y=0.985)
    fig.tight_layout(rect=(0.02, 0, 1, 0.97))
    fig.savefig(os.path.join(OUT, "fig", "fig1_layers.png"))
    plt.close(fig)

    # ============================================================= Figure 2
    # Amplitude distribution of both layers + truncation location.
    fig = plt.figure(figsize=(9.6, 5.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], hspace=0.42, wspace=0.26)

    chcol = {"R": ORANGE, "G": AQUA, "B": BLUE}
    sch = read("scale_scan")
    ax = fig.add_subplot(gs[0, 0])
    xs = [fget(r, "sigma_px") for r in scale_scan]
    ys = [fget(r, "hp_cross_corr") for r in scale_scan]
    ax.plot(xs, ys, color=BLUE, lw=2, solid_capstyle="round")
    ax.axvline(2.75, color=INK, lw=1)
    ax.annotate("选定 σ = 2.75 px", xy=(2.75, 0.05), xytext=(4.0, 0.12),
                fontsize=8.5, color=INK2,
                arrowprops=dict(arrowstyle="-", color=AXIS, lw=0.8))
    ax.set_xlabel("高斯 σ（像素）")
    ax.set_ylabel("高通带三通道平均相关系数")
    ax.set_title("(a) 尺度扫描：小尺度各通道独立，\n大尺度三通道共同变化",
                 fontsize=9, color=INK)
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 0.75)

    # structure histogram in reflectance, with clip bounds
    ax = fig.add_subplot(gs[0, 1])
    hs = read("hist_structure")
    for name in "RGB":
        rs = [r for r in hs if r["channel"] == name]
        centers = [0.5 * (fget(r, "bin_lo") + fget(r, "bin_hi")) for r in rs]
        counts = [fget(r, "count") for r in rs]
        tot = sum(counts) * (centers[1] - centers[0])
        ax.plot(centers, [c / tot for c in counts], color=chcol[name], lw=1.6, label=name)
    ymax = ax.get_ylim()[1]
    for x, lab in ((0.030, "下限 0.030"), (0.920, "上限 0.920")):
        ax.axvline(x, color=INK, lw=1)
    ax.text(0.045, ymax * 0.86, "下限 0.030", fontsize=7.5, color=INK2, ha="left")
    ax.text(0.905, ymax * 0.86, "上限 0.920", fontsize=7.5, color=INK2, ha="right")
    ax.set_xlabel("结构层 S（反射率）")
    ax.set_ylabel("概率密度")
    ax.set_xlim(0, 1.0)
    ax.set_title("(b) 结构层幅值 0.13–0.63，上限留 0.31 余量",
                 fontsize=9, color=INK)
    ax.legend(frameon=False, fontsize=8, loc="center right")

    # micro histogram in relative percent
    ax = fig.add_subplot(gs[0, 2])
    hm = read("hist_micro")
    for name in "RGB":
        rs = [r for r in hm if r["channel"] == name]
        centers = [0.5 * (fget(r, "bin_lo") + fget(r, "bin_hi")) for r in rs]
        counts = [fget(r, "count") for r in rs]
        tot = sum(counts) * (centers[1] - centers[0])
        ax.plot(centers, [c / tot for c in counts], color=chcol[name], lw=1.6, label=name)
    lo = (np.exp(-0.130) - 1.0) * 100.0
    hi = (np.exp(0.130) - 1.0) * 100.0
    ymax = ax.get_ylim()[1]
    for x in (lo, hi):
        ax.axvline(x, color=INK, lw=1)
    ax.text(lo + 0.6, ymax * 0.9, f"截断 {lo:.1f}%", fontsize=7.5, color=INK2, ha="left")
    ax.text(hi - 0.6, ymax * 0.9, f"截断 +{hi:.1f}%", fontsize=7.5, color=INK2, ha="right")
    ax.set_xlabel("微噪声层 m（相对百分比）")
    ax.set_ylabel("概率密度")
    ax.set_title("(c) 微噪声层 σ=4.5%，约 0.5% 像素触及截断",
                 fontsize=9, color=INK)
    ax.legend(frameon=False, fontsize=8, loc="center right")

    # truncation onset
    ax = fig.add_subplot(gs[1, :])
    for layer, col, lab in (("structure", BLUE, "结构层"), ("micro", ORANGE, "微噪声层")):
        rs = [r for r in trunc if r["layer"] == layer and r["channel"] == "R"]
        g = [fget(r, "gain") for r in rs]
        cf = [fget(r, "clip_frac") * 100.0 for r in rs]
        ax.plot(g, cf, color=col, lw=2, marker="o", ms=4, label=lab + "（R 通道）")
    ax.axvline(1.0, color=INK, lw=1)
    ax.text(1.01, 9.0, "标称增益 k=1", fontsize=8.5, color=INK2, ha="left")
    ax.set_xlabel("该层缩放系数 k")
    ax.set_ylabel("被截断像素比例（%）")
    ax.set_title("(d) 截断位置：结构层到 k≈1.4 才触及上限；微噪声层在标称增益已有约 0.5% 触界",
                 fontsize=9, color=INK)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    fig.suptitle("图 2  两层各自的幅值范围与截断位置", fontsize=11.5, y=0.985)
    fig.tight_layout(rect=(0, 0, 1, 0.955))
    fig.savefig(os.path.join(OUT, "fig", "fig2_amplitudes.png"))
    plt.close(fig)

    # ============================================================= Figure 3
    # Both layers on: interference sample image + measured cross terms.
    fig = plt.figure(figsize=(9.6, 5.0))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.15], wspace=0.3)
    a = ppm("albedo")
    b = ppm("interference_sample")
    x0, y0, s = 150, 120, 280
    ax = fig.add_subplot(gs[0, 0])
    ax.imshow(crop_rgb(a, x0, y0, s), interpolation="nearest")
    ax.set_title("两层同时开启", fontsize=9, color=INK)
    ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
    ax = fig.add_subplot(gs[0, 1])
    ax.imshow(crop_rgb(b, x0, y0, s), interpolation="nearest")
    ax.set_title("交叉项（放大 ×30；中灰 = 无交叉）", fontsize=9, color=INK)
    ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)

    ax = fig.add_subplot(gs[0, 2])
    # bar chart: cross term in albedo vs display, per weather; same-channel
    labels = []
    vals_ab = []
    vals_db = []
    for w in ("clear", "overcast", "dusk"):
        labels.append({"clear": "晴天", "overcast": "阴天", "dusk": "黄昏"}[w])
        vals_ab.append(fget([r for r in cross if r["weather"] == w and r["space"] == "albedo"][0], "max_cross_pct"))
        vals_db.append(fget([r for r in cross if r["weather"] == w and r["space"] == "display"][0], "max_cross_pct"))
    x = np.arange(3)
    wdt = 0.34
    ax.bar(x - wdt / 2 - 0.01, vals_ab, wdt, color=BLUE, label="albedo 空间")
    ax.bar(x + wdt / 2 + 0.01, vals_db, wdt, color=ORANGE, label="显示空间")
    for xi, va, vb in zip(x, vals_ab, vals_db):
        ax.text(xi - wdt / 2 - 0.01, va + 0.06, f"{va:.2f}", ha="center", fontsize=8, color=INK2)
        ax.text(xi + wdt / 2 + 0.01, vb + 0.06, f"{vb:.2f}", ha="center", fontsize=8, color=INK2)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("结构层 +20% 时微噪声层变化（%）")
    ax.set_ylim(0, 2.5)
    ax.set_title("交叉项：同通道、两层同时开启", fontsize=9, color=INK)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    fig.suptitle("图 3  两层同时开启的干涉样本：结构层 +20% 只带动微噪声层 <2.2%",
                 fontsize=11.5, y=0.985)
    fig.tight_layout(rect=(0, 0, 1, 0.945))
    fig.savefig(os.path.join(OUT, "fig", "fig3_interference.png"))
    plt.close(fig)

    # ============================================================= Figure 4
    # Knob independence + sweep (what the user asked: scale each layer freely)
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.4), sharey=True)
    fig.patch.set_facecolor(SURFACE)
    names = ["clear", "overcast", "dusk"]
    cn = {"clear": "晴天", "overcast": "阴天", "dusk": "黄昏"}
    for ax, w in zip(axes, names):
        rows = [r for r in ssweep if r["weather"] == w and r["driven_channel"] == "R" and r["channel"] == "R"]
        g = [fget(r, "gain") for r in rows]
        own = [fget(r, "dev_change_pct") for r in rows]
        crossb = [fget(r, "cross_band_change_pct") for r in rows]
        ax.plot(g, own, color=BLUE, lw=2, label="红通道自身偏差变化")
        ax.plot(g, crossb, color=ORANGE, lw=2, label="微噪声层变化（串扰）")
        ax.axhline(0, color=AXIS, lw=0.8)
        ax.axvline(1.0, color=AXIS, lw=0.8)
        ax.set_title(f"{cn[w]}", fontsize=9, color=INK)
        ax.set_xlabel("结构层红通道系数 k_s")
        if ax is axes[0]:
            ax.set_ylabel("变化（%）")
    axes[0].legend(frameon=False, fontsize=7.5, loc="upper left")
    fig.suptitle("图 4  结构层系数 k_s（红通道）：红通道偏差随之变化；微噪声层几乎不动（晴天 ±0.34%，阴天 ±1.9%，黄昏 ±0.11%）",
                 fontsize=10, y=0.98)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(os.path.join(OUT, "fig", "fig4_sweep.png"))
    plt.close(fig)

    # ============================================================= Figure 5
    # Three weather conditions: channel bias (linear vs display) and stability.
    fig = plt.figure(figsize=(9.6, 5.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], width_ratios=[1.3, 1],
                          wspace=0.28, hspace=0.62)

    order = ["clear", "overcast", "dusk"]
    cnl = {"clear": "晴天", "overcast": "阴天", "dusk": "黄昏"}

    # (a) per-channel bias, linear (bars) with display overlay (diamonds)
    ax = fig.add_subplot(gs[0, 0])
    x = np.arange(3)
    for i, ch in enumerate("RGB"):
        lin = [fget([r for r in weather if r["weather"] == w and r["channel"] == ch][0], "dev_linear_pct") for w in order]
        disp = [fget([r for r in weather if r["weather"] == w and r["channel"] == ch][0], "dev_pct") for w in order]
        ax.bar(x + (i - 1) * 0.25, lin, 0.22, color=chcol[ch], label=ch)
        ax.plot(x + (i - 1) * 0.25, disp, linestyle="none", marker="D", ms=5,
                mfc=SURFACE, mec=chcol[ch], mew=1.5)
    ax.axhline(0, color=AXIS, lw=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels([cnl[w] for w in order])
    ax.set_ylabel("通道偏差 / 平均亮度（%）")
    ax.set_ylim(-22, 34)
    ax.set_title("(a) 三种天气下的通道偏差\n实心柱 = 线性偏色；空心菱形 = 色调曲线后的显示值", fontsize=9, color=INK)
    ax.legend(frameon=True, framealpha=0.0, edgecolor="none", fontsize=8, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.30))

    # (b) display micro spread across channels
    ax = fig.add_subplot(gs[0, 1])
    spread = [fget(stab[w], "display_micro_spread_pct") for w in order]
    bars = ax.bar(x, spread, 0.5, color=BLUE)
    for xi, v in zip(x, spread):
        ax.text(xi, v * 1.25 if v > 1 else v + 2, f"{v:.0f}%", ha="center", fontsize=8, color=INK2)
    ax.set_xticks(x)
    ax.set_xticklabels([cnl[w] for w in order])
    ax.set_ylabel("三通道微噪声幅值离散（%）")
    ax.set_yscale("log")
    ax.set_ylim(1, 600)
    ax.set_title("(b) 微噪声层在三种天气下的通道均匀性", fontsize=9, color=INK)

    # (c) crosstalk: layered vs merged, per-channel push
    ax = fig.add_subplot(gs[1, :])
    lay = [fget(stab[w], "layered_bias_drift_pct") for w in order]
    mer = [fget(stab[w], "merged_bias_drift_pct") for w in order]
    ax.bar(x - 0.19 - 0.01, mer, 0.36, color=ORANGE, label="合并单层（现状）")
    ax.bar(x + 0.19 + 0.01, lay, 0.36, color=BLUE, label="两层方案")
    for xi, a_, b_ in zip(x, mer, lay):
        ax.text(xi - 0.2, a_ + 0.09, f"{a_:.2f}", ha="center", fontsize=8, color=INK2)
        ax.text(xi + 0.2, b_ + 0.09, f"{b_:.2f}", ha="center", fontsize=8, color=INK2)
    ax.set_xticks(x)
    ax.set_xticklabels([cnl[w] for w in order])
    ax.set_ylabel("另一通道显示偏差变化（%）")
    ax.set_ylim(0, 7.9)
    ax.set_title("(c) 单通道 +20 % 时，另一通道偏差被带动的幅度（同口径：两边都是逐通道推）",
                 fontsize=9, color=INK)
    ax.legend(frameon=False, fontsize=8, loc="upper left", ncol=2)
    fig.suptitle("图 5  三种天气的通道偏差；阴天最稳，黄昏因色调曲线肩部最不稳",
                 fontsize=11.5, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(OUT, "fig", "fig5_weather.png"))
    plt.close(fig)
    print("figures written to", os.path.join(OUT, "fig"))


if __name__ == "__main__":
    main()
