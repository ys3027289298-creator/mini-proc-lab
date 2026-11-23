#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""读取 snowmelt C 程序输出的 float32 网格，渲染因子场、非对称统计与三条雪线对比图。

用法: python3 -I plot_snowmelt.py <grid_dir> <out_dir>
图表配色经 dataviz validate_palette 全配对校验（light/dark 均通过硬门槛）。
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow

# 显式注册本机中文字体（ttc 不在 matplotlib 默认缓存中时保证可复现）
_FONT_PATH = "/home/ys/.local/share/fonts/wqy-zenhei.ttc"
if os.path.exists(_FONT_PATH):
    font_manager.fontManager.addfont(_FONT_PATH)
_CJK = "WenQuanYi Zen Hei"

# ---- 调色板（categorical 固定顺序，勿换序）----
C_S = "#2a78d6"    # slot1 blue   —— 低雪线 / 南坡
C_N = "#eb6834"    # slot2 orange —— 中雪线 / 北坡
C_EW = "#1baf7a"   # slot3 aqua   —— 高雪线 / 东西坡
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
SURF = "#fcfcfb"
SEQ = plt.cm.Blues  # 顺序色：单色蓝 light->dark

plt.rcParams.update({
    "font.sans-serif": [_CJK, "DejaVu Sans"],
    "font.weight": 500,  # WQY ttc 仅含单一字重 500，显式指定避免回退告警
    "axes.titleweight": 500,
    "axes.labelweight": 500,
    "figure.titleweight": 500,
    "axes.unicode_minus": False,
    "font.size": 10,
    "text.color": INK,
    "axes.edgecolor": "#c3c2b7",
    "axes.labelcolor": INK2,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.8,
    "figure.facecolor": SURF,
    "axes.facecolor": SURF,
})

N = 96


def load(grid_dir, name):
    return np.fromfile(os.path.join(grid_dir, name), dtype=np.float32).reshape(N, N)


def hillshade(z, sun_h=np.deg2rad(45), sun_a=0.0):
    """与 C 端同约定：x=东 y=南，太阳方位正南，按入射余弦做晕渲。
    显示灰度做地板抬升（0.22+0.78·cosi），避免背阴截断区纯黑丢失地形。"""
    gy, gx = np.gradient(z)
    beta = np.arctan(np.hypot(gx, gy))
    phi = np.arctan2(-gx, -gy)
    cosi = (np.sin(sun_h) * np.cos(beta)
            + np.cos(sun_h) * np.sin(beta) * np.cos(sun_a - phi))
    cosi = np.clip(cosi, 0, None)
    return 0.22 + 0.78 * cosi, beta, phi


def main():
    grid_dir, out_dir = sys.argv[1], sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)
    z = load(grid_dir, "height.f32")
    F = load(grid_dir, "factor.f32")
    slope = load(grid_dir, "slope.f32")
    aspect = load(grid_dir, "aspect.f32")
    shade, _, _ = hillshade(z)
    extent = (0, N, 0, N)  # origin=lower：下方为南、上方为北

    # 坡向分组（与 C 端 aspect_class 一致）
    deg = np.rad2deg(aspect)
    flat = slope < np.deg2rad(5)
    south = (~flat) & (np.abs(np.cos(aspect)) >= np.abs(np.sin(aspect))) & (np.cos(aspect) > 0)
    north = (~flat) & (np.abs(np.cos(aspect)) >= np.abs(np.sin(aspect))) & (np.cos(aspect) < 0)
    east = (~flat) & (np.abs(np.sin(aspect)) > np.abs(np.cos(aspect))) & (np.sin(aspect) > 0)
    west = (~flat) & (np.abs(np.sin(aspect)) > np.abs(np.cos(aspect))) & (np.sin(aspect) < 0)

    groups = [("南坡", south, C_S), ("北坡", north, C_N),
              ("东坡", east, C_EW), ("西坡", west, C_EW),
              ("平坦", flat, MUTED)]
    means = {name: float(F[mask].mean()) for name, mask, _ in groups}
    asym_rel = (means["南坡"] - means["北坡"]) / (0.5 * (means["南坡"] + means["北坡"]))

    # ============ 图 1：因子场 + 非对称统计 + 理论推导曲线 ============
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle("日照角×坡向差融雪修正因子（太阳高度角 h=45°，正午正南入射）",
                 fontsize=14, color=INK)

    ax = axes[0, 0]
    ax.imshow(shade, cmap="gray", origin="lower", extent=extent, vmin=0, vmax=1)
    cs = ax.contour(z, levels=8, colors="white", linewidths=0.6, alpha=0.35,
                    origin="lower", extent=extent)
    ax.add_patch(FancyArrow(N / 2, N * 0.86, 0, -N * 0.16, width=0.25,
                            head_width=1.6, head_length=2.2, color=C_N,
                            length_includes_head=True, label="太阳入射方向"))
    ax.text(N / 2, N * 0.65, "正南日照", color=INK2, fontsize=9, ha="center",
            bbox=dict(facecolor=SURF, edgecolor="none", pad=1.2))
    ax.set_title("(a) 合成地形晕渲（白细线为等高线）", color=INK2)
    ax.set_xticks([]); ax.set_yticks([])
    ax.text(1.2, N - 4.5, "北 N", color=INK, fontsize=10,
            bbox=dict(facecolor=SURF, edgecolor="none", pad=1.5))
    ax.text(1.2, 2.5, "南 S", color=INK, fontsize=10,
            bbox=dict(facecolor=SURF, edgecolor="none", pad=1.5))

    ax = axes[0, 1]
    im = ax.imshow(F, cmap=SEQ, origin="lower", extent=extent, vmin=0, vmax=1.5)
    ax.set_title("(b) 融雪修正因子 F 场：南坡>1（融雪强），北坡≈0（背阴截断）", color=INK2)
    ax.set_xticks([]); ax.set_yticks([])
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("F（平坦基准 = 1）", color=INK2)
    cb.outline.set_linewidth(0.8)

    ax = axes[1, 0]
    names = [g[0] for g in groups]
    vals = [means[n] for n in names]
    cols = [g[2] for g in groups]
    bars = ax.bar(names, vals, color=cols, width=0.62, edgecolor="none")
    ax.axhline(1.0, color=MUTED, linewidth=1.2)
    ax.text(0.02, 1.04, "平坦基准 F=1", color=MUTED, fontsize=8.5,
            transform=ax.get_yaxis_transform())
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.03, f"{v:.3f}",
                ha="center", va="bottom", fontsize=9, color=INK)
    ax.set_ylim(0, 1.5)
    ax.set_ylabel("平均融雪修正因子 F", color=INK2)
    ax.set_title(f"(c) 各坡向平均 F：南北相对不对称量 {asym_rel*100:.1f}%（门槛 20%）",
                 color=INK2)

    ax = axes[1, 1]
    beta = np.linspace(0, np.deg2rad(60), 200)
    tanh = np.tan(np.deg2rad(45))
    f_s = np.cos(beta) + tanh * np.sin(beta)
    f_n = np.maximum(np.cos(beta) - tanh * np.sin(beta), 0)
    f_ew = np.cos(beta)
    series = [("南坡 F=cosβ+tan h·sinβ", f_s, C_S),
              ("北坡 F=max(cosβ-tan h·sinβ,0)", f_n, C_N),
              ("东西坡 F=cosβ", f_ew, C_EW)]
    for label, ys_, c in series:
        ax.plot(np.rad2deg(beta), ys_, lw=2, color=c, label=label)
        ax.text(60.6, ys_[-1], label.split()[0], color=INK2, fontsize=9, va="center")
    ax.axhline(1.0, color=MUTED, linewidth=1.2)
    ax.set_xlim(0, 70); ax.set_ylim(0, 1.5)
    ax.set_xlabel("坡度 β（度）", color=INK2)
    ax.set_ylabel("融雪修正因子 F", color=INK2)
    ax.set_title("(d) 推导式 F(β, φ) 理论曲线：南北坡随坡度单调张开", color=INK2)
    ax.legend(frameon=False, fontsize=8.5, loc="upper left")

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    p1 = os.path.join(out_dir, "factor_field.png")
    fig.savefig(p1, dpi=150)
    plt.close(fig)

    # ============ 图 2：三条雪线基准 vs 修正对比 ============
    quants = [0.45, 0.60, 0.75]
    line_colors = [C_S, C_N, C_EW]
    fig, axes = plt.subplots(3, 2, figsize=(10.5, 13.2))
    fig.suptitle("三条雪线：水平基准（左）与坡向修正后（右）对比",
                 fontsize=14, color=INK)

    for t, (q, lc) in enumerate(zip(quants, line_colors)):
        z0 = float(np.quantile(z, q))
        base = load(grid_dir, f"snow_base{t}.f32")
        adj = load(grid_dir, f"snow_adj{t}.f32")
        thr = z0 * (1 + 0.35 * (F - 1))
        z_s = float(thr[south].mean()); z_n = float(thr[north].mean())

        for col, (mask, title_prefix, contour_color) in enumerate([
            (base, f"基准雪线 {t+1}（高程分位 {q:.2f}，z0={z0:.0f} m，水平直线）", MUTED),
            (adj, f"修正雪线 {t+1}：随 F 弯曲，南坡上移/北坡下移", lc),
        ]):
            ax = axes[t, col]
            ax.imshow(shade, cmap="gray", origin="lower", extent=extent, vmin=0, vmax=1)
            snow = np.ma.masked_where(mask < 0.5, mask)
            ax.imshow(snow, cmap=matplotlib.colors.ListedColormap(["#dff0ff"]),
                      origin="lower", extent=extent, alpha=0.9, vmin=0, vmax=1)
            if col == 0:
                ax.contour(z, levels=[z0], colors=[contour_color], linewidths=1.6,
                           origin="lower", extent=extent)
            else:
                # 修正雪线边界：z(x,y) = z0·(1+K(F-1))，即 z-thr=0 的等值线
                ax.contour(z - thr, levels=[0], colors=[contour_color], linewidths=2.0,
                           origin="lower", extent=extent)
            ax.set_title(title_prefix, fontsize=10, color=INK2)
            ax.set_xticks([]); ax.set_yticks([])
            ax.text(1.2, N - 4.5, "北 N", color=INK, fontsize=9,
                    bbox=dict(facecolor=SURF, edgecolor="none", pad=1.2))
            ax.text(1.2, 2.5, "南 S", color=INK, fontsize=9,
                    bbox=dict(facecolor=SURF, edgecolor="none", pad=1.2))
            if col == 1:
                ax.text(N / 2, N - 4,
                        f"南坡雪线 {z_s:.0f} m | 北坡雪线 {z_n:.0f} m | 差距 {z_s-z_n:.0f} m",
                        ha="center", color=INK, fontsize=9,
                        bbox=dict(facecolor=SURF, edgecolor="#c3c2b7",
                                  boxstyle="round,pad=0.3", linewidth=0.8))

    fig.tight_layout(rect=(0, 0, 1, 0.97))
    p2 = os.path.join(out_dir, "snowline_compare.png")
    fig.savefig(p2, dpi=150)
    plt.close(fig)

    # ============ 表格视图（relief：aqua 浅色色块的可读性补偿 + 完整数值）============
    csv = os.path.join(out_dir, "asymmetry_table.csv")
    with open(csv, "w", encoding="utf-8") as fh:
        fh.write("坡向,样本格点数,平均F,平均高程_m\n")
        for name, mask, _ in groups:
            fh.write(f"{name},{int(mask.sum())},{F[mask].mean():.4f},{z[mask].mean():.1f}\n")
        fh.write(f"南北绝对差,_,{means['南坡']-means['北坡']:.4f},_\n")
        fh.write(f"南北相对不对称量,_,{asym_rel:.4f} ({asym_rel*100:.1f}%),_\n")
    print("wrote:", p1, p2, csv, sep="\n")
    print(f"north_south_rel={asym_rel:.4f} ({asym_rel*100:.1f}%)")


if __name__ == "__main__":
    main()
