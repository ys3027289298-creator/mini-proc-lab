#!/usr/bin/env python3
"""Plot the three-terrain mine-alignment results produced by ../tunnel.

Inputs  : artifacts/data/*.csv (written by the C tunnel binary)
Outputs : artifacts/mine_alignment.png  (three-terrain figure)
          artifacts/deviation_stats.csv (machine-readable acceptance table)

Run: python3 plot_alignment.py [data_dir] [out_png] [out_csv]
"""
import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

# ---- validated categorical palette (dataviz validator, light mode) ----
BLUE = "#2a78d6"      # new scheme (aligned centerline)
ORANGE = "#eb6834"    # old scheme (fixed 4 m Euler chords)
GREEN = "#0ca30c"     # pass status (text only; never adjacent to orange)
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
SURFACE = "#fcfcfb"
BUDGET_FILL = "#cde2fb"   # blue ramp step 100
RELIEF_CMAP = mcolors.LinearSegmentedColormap.from_list(
    "warm_gray", ["#e9e6dd", "#cfcabd", "#b3ada0", "#978f80", "#7c7364", "#615a4d"]
)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Noto Sans CJK SC", "DejaVu Sans"],
    "axes.edgecolor": "#c3c2b7",
    "axes.labelcolor": INK2,
    "axes.titlecolor": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "text.color": INK,
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.8,
    "font.size": 10,
})

TERRAINS = [
    ("A", "rolling_hills", "rolling_hills 缓丘"),
    ("B", "single_ridge", "single_ridge 单脊"),
    ("C", "steep_bent_monocline", "steep_bent_monocline 陡弯单斜"),
]

# analytic fields duplicated from tunnel.c (for the relief background only)
def field(name, x, y):
    if name == "rolling_hills":
        h = (0.25 * x + 0.8 * np.sin(0.14 * x) * np.cos(0.11 * y + 0.3)
             + 0.35 * np.sin(0.10 * y + 1.1))
    elif name == "single_ridge":
        q = y - (10.0 + 0.30 * x)
        h = 0.45 * x + 0.10 * y + 4.0 * np.exp(-q * q / 50.0)
    else:
        q = y - (26.0 - 0.45 * x)
        h = (1.15 * x + 6.0 * np.exp(-q * q / 24.5)
             + 1.2 * np.sin(0.18 * x + 0.7) * np.cos(0.16 * y))
    return h


def load_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def read_path(data, name, scheme):
    rows = load_csv(os.path.join(data, f"{scheme}_{name}.csv"))
    return np.array([[float(r["x"]), float(r["y"])] for r in rows])


def read_stream(data, name):
    rows = load_csv(os.path.join(data, f"stream_{name}.csv"))
    s = np.array([float(r["s"]) for r in rows])
    x = np.array([float(r["x"]) for r in rows])
    y = np.array([float(r["y"]) for r in rows])
    k = np.array([float(r["kappa"]) for r in rows])
    return s, x, y, k


def main():
    data = sys.argv[1] if len(sys.argv) > 1 else "artifacts/data"
    out_png = sys.argv[2] if len(sys.argv) > 2 else "artifacts/mine_alignment.png"
    out_csv = sys.argv[3] if len(sys.argv) > 3 else "artifacts/deviation_stats.csv"

    summary = {r["terrain"]: r for r in load_csv(os.path.join(data, "summary.csv"))}

    fig = plt.figure(figsize=(17.5, 14.8), dpi=140)
    gs = fig.add_gridspec(
        3, 3, width_ratios=[2.1, 1.0, 1.25], hspace=0.46, wspace=0.32,
        left=0.055, right=0.975, top=0.885, bottom=0.135,
    )

    table_rows = []
    for row, (tag, name, title) in enumerate(TERRAINS):
        s, sx, sy, sk = read_stream(data, name)
        new = read_path(data, name, "new")
        old = read_path(data, name, "old")
        seg = load_csv(os.path.join(data, f"segstats_{name}.csv"))
        new_seg = [r for r in seg if r["scheme"] == "new"]
        old_seg = [r for r in seg if r["scheme"] == "old"]
        sr = summary[name]

        # ---------- column 1: relief + paths ----------
        ax = fig.add_subplot(gs[row, 0])
        xg = np.linspace(0, 48, 320)
        yg = np.linspace(0, 48, 320)
        X, Y = np.meshgrid(xg, yg)
        H = field(name, X, Y)
        ax.pcolormesh(X, Y, H, cmap=RELIEF_CMAP, shading="auto", alpha=0.92)
        cs = ax.contour(X, Y, H, levels=14, colors="white", linewidths=0.7, alpha=0.55)
        ax.clabel(cs, inline=True, fontsize=6.5, fmt="%.0f", colors=INK2)

        # ideal flow line, thin dark reference
        ax.plot(sx, sy, color=INK, lw=1.1, ls=(0, (4, 2)), alpha=0.85,
                label="理想梯度流线（2 mm 步长 RK2）", zorder=4)
        # old scheme
        ax.plot(old[:, 0], old[:, 1], color=ORANGE, lw=2.0, marker="o",
                ms=5, mfc="white", mec=ORANGE, mew=1.2, zorder=5,
                label="旧方案：等距 4 m、每站单次采样")
        # new scheme
        ax.plot(new[:, 0], new[:, 1], color=BLUE, lw=2.4, marker="o",
                ms=5.5, mfc=BLUE, mec="white", mew=0.8, zorder=6,
                label="对齐中线（本站，自适应弦长）")

        # zoom window around the curvature peak: fixed 6x, kept inside domain
        ZOOM_HALF = 4.0
        ipk = int(np.argmax(sk))
        cx, cy = sx[ipk], sy[ipk]
        x0 = min(max(cx - ZOOM_HALF, 0.5), 48.0 - 2 * ZOOM_HALF - 0.5)
        y0 = min(max(cy - ZOOM_HALF, 0.5), 48.0 - 2 * ZOOM_HALF - 0.5)
        x1, y1 = x0 + 2 * ZOOM_HALF, y0 + 2 * ZOOM_HALF
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False,
                               ec=INK, lw=1.0, ls=(0, (3, 2)), zorder=7))

        ax.set_xlim(0, 48)
        ax.set_ylim(0, 48)
        ax.set_aspect("equal")
        ax.set_xlabel("x / m")
        ax.set_ylabel("y / m")
        ax.set_title(
            f"{tag}  {title}\n落差 {float(sr['drop_m']):.1f} m · "
            f"{sr['segments']} 段 · 曲率峰值 {float(sr['kappa_peak']):.3f} /m",
            loc="left", fontsize=11, pad=8)

        # ---------- column 2: zoom inset ----------
        axz = fig.add_subplot(gs[row, 1])
        axz.pcolormesh(X, Y, H, cmap=RELIEF_CMAP, shading="auto", alpha=0.92)
        axz.contour(X, Y, H, levels=14, colors="white", linewidths=0.7, alpha=0.55)
        axz.plot(sx, sy, color=INK, lw=1.0, ls=(0, (4, 2)), zorder=4)
        axz.plot(old[:, 0], old[:, 1], color=ORANGE, lw=2.0, marker="o",
                 ms=5, mfc="white", mec=ORANGE, mew=1.1, zorder=5)
        axz.plot(new[:, 0], new[:, 1], color=BLUE, lw=2.4, marker="o",
                 ms=5.5, mfc=BLUE, mec="white", mew=0.8, zorder=6)
        axz.set_xlim(x0, x1)
        axz.set_ylim(y0, y1)
        axz.set_aspect("equal")
        axz.grid(False)
        axz.set_title("局部 ×6\n曲率峰值附近", fontsize=9.5, pad=6)
        axz.tick_params(labelsize=8)

        # ---------- column 3: per-segment deviation vs budget ----------
        axd = fig.add_subplot(gs[row, 2])
        smid = np.array([(float(r["s_start"]) + float(r["s_end"])) / 2
                         for r in new_seg])
        send = np.array([float(r["s_end"]) for r in new_seg])
        budget = np.array([float(r["d_allow"]) for r in new_seg])
        dmax = np.array([float(r["d_max"]) for r in new_seg])
        # step budget edges: prepend first start so where='post' covers s=0
        s_edges = np.concatenate([[float(new_seg[0]["s_start"])], send])
        b_edges = np.concatenate([[budget[0]], budget])
        axd.fill_between(s_edges, b_edges, step="post", color=BUDGET_FILL,
                         alpha=0.75, zorder=1, label="允许偏移 d_allow(k) 区间")
        axd.step(s_edges, b_edges, where="post", color=BLUE, lw=1.4,
                 ls=":", zorder=3, label="允许偏移 d_allow(k)")
        axd.plot(smid, dmax, color=INK, lw=1.8, marker="o", ms=4.5,
                 mfc=INK, mec="white", mew=0.7, zorder=4, label="实测 d_max(k)")

        old_fail = sum(1 for r in old_seg if r["pass"] == "0")
        old_dmax = max(float(r["d_max"]) for r in old_seg)
        axd.annotate(
            f"旧方案最大横向偏差 {old_dmax:.2f} m，{old_fail}/{len(old_seg)} 段超标",
            xy=(0.04, 0.93), xycoords="axes fraction", fontsize=8.5,
            color=ORANGE, ha="left", va="top",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=ORANGE, lw=0.9, alpha=0.9))
        verdict = "全部段达标" if sr["pass"] == "1" else "存在超标段"
        axd.set_title(f"{tag}  逐段偏差 · {verdict}", fontsize=10.5,
                      color=GREEN if sr["pass"] == "1" else "#d03b3b", pad=6)
        axd.set_xlabel("沿中线里程 s / m")
        axd.set_ylabel("横向偏移 d / m")
        axd.set_xlim(0, 40)
        axd.set_ylim(0, max(0.21, float(sr["d_allow_min_m"]) * 1.6,
                            budget.max() * 1.15))
        axd.tick_params(labelsize=8.5)

        table_rows.append({
            "tag": tag, "terrain": name,
            "drop_m": float(sr["drop_m"]),
            "segments": int(sr["segments"]),
            "d_max_m": float(sr["d_max_m"]),
            "d_allow_min_m": float(sr["d_allow_min_m"]),
            "util_max_pct": 100 * float(sr["util_max"]),
            "theta_deg_max": float(sr["theta_deg_max"]),
            "old_d_max_m": float(sr["old_d_max_m"]),
            "old_failed": old_fail,
            "old_total": len(old_seg),
            "pass": sr["pass"] == "1",
        })

    # ---- shared legend (top) ----
    handles = [
        Line2D([0], [0], color=INK, lw=1.1, ls=(0, (4, 2)), label="理想梯度流线（2 mm 步长 RK2）"),
        Line2D([0], [0], color=BLUE, lw=2.4, marker="o", ms=6, mfc=BLUE,
               mec="white", label="对齐中线（本站，自适应弦长）"),
        Line2D([0], [0], color=ORANGE, lw=2.0, marker="o", ms=6, mfc="white",
               mec=ORANGE, label="旧方案：等距 4 m、每站单次采样"),
    ]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.975, 0.965),
               ncol=3, frameon=False, fontsize=10)
    fig.suptitle("矿道中线强制对齐地形梯度 · 三段地形验收",
                 x=0.055, ha="left", fontsize=17, fontweight="bold", y=0.975)
    fig.text(0.055, 0.936,
             "允许偏移 d_allow(k) = min( e_chord , (1−cos α)/κ_k )，"
             "α = 8°，e_chord = 0.20 m；实测 d_max 全程位于蓝色预算线下方即逐段达标",
             fontsize=10.5, color=INK2)

    # ---- bottom acceptance table ----
    col_labels = ["段", "地形", "落差 m", "段数", "d_max max m",
                  "d_allow min m", "上限利用率", "走向偏差 max",
                  "旧方案偏差 m（超标段）", "判定"]
    cell_text, colors = [], []
    for r in table_rows:
        cell_text.append([
            r["tag"], r["terrain"], f"{r['drop_m']:.1f}", str(r["segments"]),
            f"{r['d_max_m']:.3f}", f"{r['d_allow_min_m']:.3f}",
            f"{r['util_max_pct']:.0f}%", f"{r['theta_deg_max']:.1f}°",
            f"{r['old_d_max_m']:.2f}（{r['old_failed']}/{r['old_total']}）",
            "达标" if r["pass"] else "未达标",
        ])
        colors.append(GREEN if r["pass"] else "#d03b3b")

    axt = fig.add_axes([0.055, 0.02, 0.92, 0.078])
    axt.axis("off")
    axt.set_title("偏差统计（任一段 d_max > d_allow 即整段未完成）",
                  loc="left", fontsize=11.5, pad=8)
    tbl = axt.table(cellText=cell_text, colLabels=col_labels,
                    cellLoc="center", loc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9.5)
    tbl.scale(1, 1.55)
    for (rr, cc), cell in tbl.get_celld().items():
        cell.set_edgecolor(GRID)
        cell.set_linewidth(0.8)
        if rr == 0:
            cell.set_facecolor("#f1f0ea")
            cell.set_text_props(color=INK2, fontweight="bold")
        elif cc == 9:
            cell.set_text_props(color=colors[rr - 1], fontweight="bold")
        if rr > 0 and cc != 9:
            cell.set_text_props(color=INK)

    fig.savefig(out_png, dpi=140, facecolor=SURFACE, bbox_inches="tight")
    print(f"wrote {out_png}")

    # machine-readable acceptance table
    with open(out_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["terrain", "drop_m", "segments", "d_max_m", "d_allow_min_m",
                    "util_max_pct", "theta_deg_max", "old_d_max_m",
                    "old_failed_segs", "old_total_segs", "pass"])
        for r in table_rows:
            w.writerow([r["terrain"], r["drop_m"], r["segments"], r["d_max_m"],
                        r["d_allow_min_m"], f"{r['util_max_pct']:.2f}",
                        r["theta_deg_max"], r["old_d_max_m"], r["old_failed"],
                        r["old_total"], int(r["pass"])])
    print(f"wrote {out_csv}")

    # hard acceptance gate, mirrored from the C checks
    ok = all(r["pass"] for r in table_rows)
    print("ACCEPTANCE:", "PASS (all segments within budget on all three terrains)"
          if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
