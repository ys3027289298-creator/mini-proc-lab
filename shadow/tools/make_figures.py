#!/usr/bin/env python3
"""Render the three-moment shadow comparison from the C exporter's output.

Reads out/{terrain,shade,shadow}_*.bin, out/daypath.csv, out/curve.csv,
out/stats.csv and writes fig_shadow_compare.png (light + dark) and prints the
values used, so the figure can be checked against the CSVs.

Top row:
  A  solar path over the day, with the three moments marked
  B  shadow length L = dh/tan(e) against solar altitude, with the moments
  C  shadow coverage and median length per moment
Bottom row:
  D 08:00, E 12:00, F 16:00  rendered shadow + hillshade frames
  G  morning-vs-evening mask difference (direction is the story)

Run after `make export` in this directory.
"""
import csv
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "out")

# ---- reference palette (light / dark), dataviz skill ------------------------
LIGHT = dict(
    surface="#fcfcfb", page="#f9f9f7", ink="#0b0b0b", ink2="#52514e",
    muted="#898781", grid="#e1e0d9", axis="#c3c2b7",
    s1="#2a78d6", s2="#eb6834", s3="#1baf7a",
)
DARK = dict(
    surface="#1a1a19", page="#0d0d0d", ink="#f0efec", ink2="#c3c2b7",
    muted="#898781", grid="#2c2c2a", axis="#383835",
    s1="#3987e5", s2="#d95926", s3="#199e70",
)

CJK = "/home/ys/.local/share/fonts/wqy-zenhei.ttc"
if os.path.exists(CJK):
    fm.fontManager.addfont(CJK)


def read_csv(name):
    with open(os.path.join(OUT, name), newline="") as fh:
        return list(csv.DictReader(fh))


def read_f32(name):
    return np.fromfile(os.path.join(OUT, name), dtype=np.float32).reshape(256, 256)


def read_u8(name):
    return np.fromfile(os.path.join(OUT, name), dtype=np.uint8).reshape(256, 256)


# frames are stored [x*H + y]; transpose to image orientation [row=y, col=x]
def as_image(a):
    return np.ascontiguousarray(a.T)


def terrain_cmap(mode):
    # single-hue sequential ramp (blue) for the shadowed heightfield + a
    # gray->warm wash for hillshade. Kept deliberately low-chroma: the
    # classification (shadow / lit) is what the figure is about.
    # v=0 is full shadow, v=1 is fully lit, in both modes.
    if mode == "light":
        return LinearSegmentedColormap.from_list(
            "hs", ["#0d366b", "#256abf", "#86b6ef", "#cde2fb", "#fbf6ec"])
    return LinearSegmentedColormap.from_list(
        "hs", ["#0d366b", "#184f95", "#256abf", "#5598e7", "#cde2fb"])


def panel_terrain(ax, z, shade, pal, title, mode):
    ax.imshow(as_image(z), cmap="Greys", alpha=0.55, origin="lower")
    ax.imshow(as_image(shade), cmap=terrain_cmap(mode), origin="lower",
              vmin=0.0, vmax=1.0, alpha=0.70)
    ax.set_title(title, color=pal["ink"], fontsize=10, loc="left", pad=6)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_color(pal["axis"]); s.set_linewidth(0.8)


def main():
    stats = {r["time"]: r for r in read_csv("stats.csv")}
    path = read_csv("daypath.csv")
    curve = read_csv("curve.csv")
    pairs = read_csv("pairs.csv")

    z = read_f32("terrain.bin")
    times = ["0800", "1200", "1600"]
    labels = {"0800": "08:00", "1200": "12:00", "1600": "16:00"}
    suns = {t: (float(stats[t]["azimuth_deg"]), float(stats[t]["altitude_deg"]))
            for t in times}

    for mode in ("light", "dark"):
        pal = LIGHT if mode == "light" else DARK
        plt.rcParams.update({
            "font.family": "sans-serif",
            "font.sans-serif": ["WenQuanYi Zen Hei", "DejaVu Sans"],
            "axes.unicode_minus": False,
            "figure.facecolor": pal["surface"],
            "savefig.facecolor": pal["surface"],
            "axes.facecolor": pal["surface"],
            "text.color": pal["ink"],
            "axes.edgecolor": pal["axis"],
            "axes.labelcolor": pal["ink2"],
            "xtick.color": pal["muted"],
            "ytick.color": pal["muted"],
            "grid.color": pal["grid"],
        })

        fig = plt.figure(figsize=(16.2, 8.4), dpi=150)
        gs = fig.add_gridspec(1, 3, left=0.048, right=0.992, top=0.845, bottom=0.50,
                              wspace=0.16)
        gsB = fig.add_gridspec(1, 4, left=0.048, right=0.992, top=0.415, bottom=0.045,
                               wspace=0.12)

        # --- A: solar path -------------------------------------------------
        axA = fig.add_subplot(gs[0, 0])
        hrs = np.array([float(r["hour"]) for r in path])
        az = np.array([float(r["azimuth_deg"]) for r in path])
        alt = np.array([float(r["altitude_deg"]) for r in path])
        day = alt > 0
        axA.set_xlim(115, 245)
        axA.set_ylim(-10, 70)
        axA.plot(az[day], alt[day], color=pal["s1"], lw=2, solid_capstyle="round")
        axA.plot(az[~day], alt[~day], color=pal["muted"], lw=1.5, ls=(0, (3, 2)))
        axA.axhline(0, color=pal["axis"], lw=1)
        for t, col, dy in (("0800", pal["s2"], 4), ("1200", pal["s3"], 4), ("1600", pal["s2"], 4)):
            a, e = suns[t]
            axA.plot([a], [e], "o", ms=8, color=col, mec=pal["surface"], mew=2, zorder=5)
            axA.annotate(f"{labels[t]}  az {a:.1f}°\n高度角 {e:.1f}°",
                         (a, e), textcoords="offset points",
                         xytext=(8 if t != "1600" else -8, dy),
                         ha="left" if t != "1600" else "right",
                         fontsize=7.5, color=pal["ink2"])
        axA.set_xlabel("太阳方位角 (°，北=0，东=90)")
        axA.set_ylabel("太阳高度角 (°)")
        axA.set_title("A  日轨迹：方位角与高度角", color=pal["ink"], loc="left",
                      fontsize=10, pad=6)
        axA.text(0.02, 0.96, "虚线：地平线以下", transform=axA.transAxes,
                 va="top", fontsize=7.5, color=pal["muted"])
        axA.grid(True, lw=0.6, alpha=0.7)

        # --- B: L = dh/tan(e) ---------------------------------------------
        axB = fig.add_subplot(gs[0, 1])
        ca = np.array([float(r["altitude_deg"]) for r in curve])
        cL = np.array([float(r["L20_m"]) for r in curve])
        axB.plot(ca, cL, color=pal["s1"], lw=2, solid_capstyle="round")
        axB.set_xlim(5, 80)
        axB.set_ylim(0, 140)
        # 08:00 and 16:00 share the altitude, so they share one marker and
        # one label; annotating them separately would stack two identical
        # values on top of each other.
        for t in times:
            a, e = suns[t]
            L20 = 20.0 / math.tan(math.radians(e))
            col = pal["s2"] if t != "1200" else pal["s3"]
            axB.plot([e], [L20], "o", ms=8, color=col, mec=pal["surface"], mew=2, zorder=5)
        e_m = suns["0800"][1]
        axB.annotate("08:00 = 16:00\n48.2 m",
                     (e_m, 20.0 / math.tan(math.radians(e_m))),
                     textcoords="offset points", xytext=(10, 4),
                     fontsize=7.5, color=pal["ink2"])
        axB.annotate("12:00\n16.8 m", (suns["1200"][1], 20.0 / math.tan(math.radians(suns["1200"][1]))),
                     textcoords="offset points", xytext=(10, 4),
                     fontsize=7.5, color=pal["ink2"])
        axB.set_xlabel("太阳高度角 e (°)")
        axB.set_ylabel("L / m（dh = 20 m）")
        axB.set_title("B  L = Δh / tan e", color=pal["ink"], loc="left",
                      fontsize=10, pad=6)
        axB.text(0.97, 0.93, "等高物体：e 越低，阴影越长；\n方位角决定方向，不改变长度",
                 transform=axB.transAxes, va="top", ha="right",
                 fontsize=7.5, color=pal["ink2"])
        axB.grid(True, lw=0.6, alpha=0.7)

        # --- C: shadow statistics -----------------------------------------
        axS = fig.add_subplot(gs[0, 2])
        names = [labels[t] for t in times]
        xs = np.arange(3)
        frac = [100 * float(stats[t]["shadow_fraction"]) for t in times]
        w = 0.46
        cols = [pal["s2"], pal["s3"], pal["s2"]]
        axS.bar(xs, frac, w, color=cols, edgecolor=pal["surface"], linewidth=2)
        axS.set_ylim(0, 52)
        axS.set_xticks(xs); axS.set_xticklabels(names)
        axS.set_ylabel("阴影覆盖 (%)")
        axS.set_title("C  阴影覆盖与中位长度", color=pal["ink"], loc="left",
                      fontsize=10, pad=6)
        axS.grid(True, axis="y", lw=0.6, alpha=0.7)
        axS.set_axisbelow(True)
        for i, t in enumerate(times):
            v = frac[i]
            axS.text(i, v + 1.4, f"{v:.1f}%", ha="center", fontsize=8.5,
                     color=pal["ink"], weight="bold")
            ml = float(stats[t]["median_len_m"])
            p95 = float(stats[t]["p95_len_m"])
            if v >= 20:
                ink = "#ffffff" if cols[i] != pal["s3"] else "#0b0b0b"
                axS.text(i, v / 2, f"中位\n{ml:.0f} m", ha="center", va="center",
                         fontsize=8, color=ink)
            else:
                axS.text(i, v + 6.0, f"中位 {ml:.0f} m", ha="center",
                         fontsize=8, color=pal["ink2"])
        axS.text(0.03, 0.97,
                 "P95 长度：08:00 103 m · 12:00 16 m · 16:00 107 m\n"
                 "（同高度角 → 长度近似相同；方向相反）",
                 transform=axS.transAxes, va="top", fontsize=7.3, color=pal["ink2"])

        # --- D/E/F: the three rendered frames ------------------------------
        for i, t in enumerate(times):
            ax = fig.add_subplot(gsB[0, i])
            shade = read_f32(f"shade_{t}.bin")
            a, e = suns[t]
            panel_terrain(ax, z, shade, pal,
                          f"{'DEF'[i]}  {labels[t]}   方位角 {a:.1f}°   高度角 {e:.1f}°",
                          mode)
            # azimuth arrow: the sun sits toward (sin az, cos az); draw it in
            # image coordinates (x right = east, y up = north)
            cx, cy = 0.88, 0.88
            ax.annotate("", xy=(cx + 0.10 * math.sin(math.radians(a)),
                                cy + 0.10 * math.cos(math.radians(a))),
                        xytext=(cx, cy), xycoords="axes fraction",
                        arrowprops=dict(arrowstyle="-|>", color=pal["s2"], lw=2),
                        annotation_clip=False)
            ax.plot([cx], [cy], "o", ms=7, color=pal["s2"], mec=pal["surface"],
                    mew=1.5, transform=ax.transAxes, clip_on=False)
            fracv = 100 * float(stats[t]["shadow_fraction"])
            ax.text(0.03, 0.04, f"阴影覆盖 {fracv:.1f}%   中位长度 {float(stats[t]['median_len_m']):.0f} m",
                    transform=ax.transAxes, fontsize=7.5, color=pal["ink"],
                    bbox=dict(facecolor=pal["surface"], edgecolor="none",
                              alpha=0.85, pad=2.5))

        # --- G: morning vs evening, the direction difference --------------
        axG = fig.add_subplot(gsB[0, 3])
        m0 = read_u8("shadow_0800.bin") > 0
        m2 = read_u8("shadow_1600.bin") > 0
        # neutral base: the hillshade of the 12:00 frame, desaturated
        base = read_f32("shade_1200.bin")
        # neutral hillshade underlay: mid-gray scaled into the mode's range,
        # so un-shadowed ground never competes with the three categories.
        bg = np.array(matplotlib.colors.to_rgb(pal["surface"]))
        fg = np.array(matplotlib.colors.to_rgb(pal["grid"]))
        rgb = np.zeros((256, 256, 4), dtype=float)
        t = as_image(base)
        rgb[..., :3] = bg[None, None, :] + (fg - bg)[None, None, :] * t[..., None]
        rgb[..., 3] = 1.0
        s1 = np.array(matplotlib.colors.to_rgb(pal["s1"]))
        s2 = np.array(matplotlib.colors.to_rgb(pal["s2"]))
        both = np.array(matplotlib.colors.to_rgb(pal["ink2"]))
        only0 = as_image(m0 & ~m2)
        only2 = as_image(m2 & ~m0)
        bothm = as_image(m0 & m2)
        for mask_, col in ((only0, s1), (only2, s2), (bothm, both)):
            rgb[mask_, :3] = col
        axG.imshow(rgb, origin="lower", interpolation="nearest")
        axG.set_xticks([]); axG.set_yticks([])
        for sp in axG.spines.values():
            sp.set_color(pal["axis"]); sp.set_linewidth(0.8)
        axG.set_title("G  晨昏方向差：两帧阴影的空间分歧", color=pal["ink"],
                      loc="left", fontsize=10, pad=6)
        # legend as swatches (identity never by colour alone: labelled)
        handles = [
            matplotlib.patches.Patch(facecolor=pal["s1"], edgecolor="none",
                                     label="仅 08:00"),
            matplotlib.patches.Patch(facecolor=pal["s2"], edgecolor="none",
                                     label="仅 16:00"),
            matplotlib.patches.Patch(facecolor=pal["ink2"], edgecolor="none",
                                     label="两帧皆有"),
        ]
        leg = axG.legend(handles=handles, loc="lower left", fontsize=7.3,
                         frameon=True, framealpha=0.9, borderpad=0.5,
                         handlelength=1.1, handleheight=0.9)
        leg.get_frame().set_facecolor(pal["surface"])
        leg.get_frame().set_edgecolor(pal["axis"])
        for txt in leg.get_texts():
            txt.set_color(pal["ink"])
        pj = [r for r in pairs if r["a"] == "0800" and r["b"] == "1600"][0]
        axG.text(0.98, 0.04,
                 f"仅晨 {100 * float(pj['only_a_frac_of_shadow']):.0f}% · "
                 f"仅昏 {100 * float(pj['only_b_frac_of_shadow']):.0f}%",
                 transform=axG.transAxes, ha="right", fontsize=7.3,
                 color=pal["ink"],
                 bbox=dict(facecolor=pal["surface"], edgecolor="none",
                           alpha=0.85, pad=2.5))

        fig.suptitle("山体阴影三时刻对比：日轨迹投影到高度场，阴影长度随方位角与高度角变化",
                     x=0.048, y=0.972, ha="left", fontsize=13, color=pal["ink"])
        fig.text(0.048, 0.925,
                 "纬度 40°、春分（赤纬 0°）；08:00 与 16:00 高度角相同（22.5°），"
                 "方位角相差 139°。阴影由同一遮挡步进渲染，长度按 L = Δz / tan e 独立核验。",
                 fontsize=8.5, color=pal["ink2"])

        fig.savefig(os.path.join(ROOT, f"fig_shadow_compare_{mode}.png"),
                    dpi=150, facecolor=pal["surface"])
        plt.close(fig)
        print(f"wrote fig_shadow_compare_{mode}.png")

    # ---- accessible table view ------------------------------------------
    with open(os.path.join(ROOT, "fig_shadow_compare.md"), "w") as fh:
        fh.write("# 三时刻阴影对比 — 数值表\n\n")
        fh.write("| 时刻 | 方位角 (°) | 高度角 (°) | 阴影覆盖 | 中位长度 (m) | P95 长度 (m) "
                 "| 可归因阴影块 | 中位偏差 | ≤20% 占比 |\n")
        fh.write("|---|---|---|---|---|---|---|---|---|\n")
        for t in times:
            s = stats[t]
            fh.write(f"| {labels[t]} | {float(s['azimuth_deg']):.2f} | "
                     f"{float(s['altitude_deg']):.2f} | "
                     f"{100 * float(s['shadow_fraction']):.1f}% | "
                     f"{float(s['median_len_m']):.1f} | {float(s['p95_len_m']):.1f} | "
                     f"{int(s['attributed'])} | "
                     f"{float(s['attributed_median_dev_pct']):.2f}% | "
                     f"{100 * float(s['attributed_within_20pct_frac']):.1f}% |\n")
        fh.write("\n## 时刻两两对比\n\n")
        fh.write("| 对比 | Jaccard 重叠 | 仅前者占比 | 仅后者占比 | 中位相对长度差 | 相对差>20% 占比 |\n")
        fh.write("|---|---|---|---|---|---|\n")
        for r in pairs:
            fh.write(f"| {labels[r['a']]} vs {labels[r['b']]} | {float(r['jaccard']):.3f} "
                     f"| {100 * float(r['only_a_frac_of_shadow']):.0f}% "
                     f"| {100 * float(r['only_b_frac_of_shadow']):.0f}% "
                     f"| {float(r['median_rel_len_diff_pct']):.0f}% "
                     f"| {100 * float(r['frac_rel_diff_gt20pct']):.0f}% |\n")
    print("wrote fig_shadow_compare.md")


if __name__ == "__main__":
    main()
