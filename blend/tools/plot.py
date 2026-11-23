#!/usr/bin/env python3
"""Plot boundary-gradient curves and cost comparison from blend CSVs.

Inputs  : artifacts/blend/{kernel_1d,boundary_trace,kernel_metrics}.csv
Outputs : artifacts/blend/fig_gradients.png, fig_cost.png, summary.json
Run     : python3 -I tools/plot.py   (from repo root or blend/)
"""
import csv
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

# ---- CJK font (wqy ttc is not in matplotlib's default cache) ----
for p in ("/home/ys/.local/share/fonts/wqy-zenhei.ttc",):
    if os.path.exists(p):
        font_manager.fontManager.addfont(p)
plt.rcParams["font.family"] = ["WenQuanYi Zen Hei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

# validated categorical palette (dataviz reference, slots 1-3, light mode)
COLORS = {"linear": "#2a78d6", "quintic": "#eb6834", "gauss": "#1baf7a"}
LABEL = {"linear": "线性 ramp", "quintic": "五次 smootherstep",
         "gauss": "截断高斯(erf)"}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e1e0d9"
SURFACE = "#fcfcfb"

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for cand in (os.path.join(ROOT, "artifacts", "blend"),
             os.path.join(os.getcwd(), "artifacts", "blend"),
             os.path.join(ROOT, "..", "artifacts", "blend")):
    if os.path.isdir(cand):
        OUT = os.path.abspath(cand)
        break
else:
    sys.exit("artifacts/blend not found")


def load(name):
    with open(os.path.join(OUT, name), newline="") as fh:
        return list(csv.DictReader(fh))


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#c3c2b7")
    ax.grid(True, color=GRID, lw=0.8, zorder=0)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.xaxis.label.set_color(MUTED)
    ax.yaxis.label.set_color(MUTED)
    ax.title.set_color(INK)


# ================= figure 1: boundary gradient (5 panels, 2 rows) =================
rows = load("kernel_1d.csv")
u = np.array(sorted({float(r["d_over_R"]) for r in rows}))
kernels = ["linear", "quintic", "gauss"]
data = {k: {"w": np.zeros(len(u)), "g": np.zeros(len(u))} for k in kernels}
def near(value):
    return int(np.argmin(np.abs(u - value)))
for r in rows:
    k = r["kernel"]
    i = near(float(r["d_over_R"]))
    data[k]["w"][i] = float(r["weight"])
    data[k]["g"][i] = float(r["grad"])

fig = plt.figure(figsize=(13.5, 8.2), facecolor=SURFACE)
gs = fig.add_gridspec(2, 3, hspace=0.34, wspace=0.28,
                      height_ratios=[1, 1.05])
ax_w = fig.add_subplot(gs[0, 0])
ax_g = fig.add_subplot(gs[0, 1])
ax_g2 = fig.add_subplot(gs[0, 2])
ax_t = fig.add_subplot(gs[1, 0])
ax_s = fig.add_subplot(gs[1, 1])
ax_h = fig.add_subplot(gs[1, 2])

ax = ax_w
for k in kernels:
    ax.plot(u, data[k]["w"], color=COLORS[k], lw=2, label=LABEL[k])
ax.axvspan(-1, 1, color="#86b6ef", alpha=0.12, zorder=0)
ax.set_title("A. 权重曲线 w(d)，蓝带 = 混合带 |d|≤R", fontsize=10.5, loc="left")
ax.set_xlabel("归一化距离 d/R"); ax.set_ylabel("材质 A 权重 w")
ax.set_xlim(-1.2, 1.2); style(ax)

ax = ax_g
for k in kernels:
    ax.plot(u, data[k]["g"], color=COLORS[k], lw=2, label=LABEL[k])
ax.set_title("B. 核梯度 dK/dd（全量程）", fontsize=10.5, loc="left")
ax.set_xlabel("归一化距离 d/R"); ax.set_ylabel("梯度 dK/dd")
ax.set_xlim(-1.2, 1.2); style(ax)

ax = ax_g2
for k in kernels:
    ax.plot(u, data[k]["g"], color=COLORS[k], lw=2, label=LABEL[k])
    if k == "gauss":
        for e in (-1.0, 1.0):
            ax.plot([e], [data[k]["g"][near(e)]], "o", ms=7,
                    markerfacecolor="none", markeredgecolor=COLORS[k],
                    markeredgewidth=1.6, zorder=5)
ax.annotate("高斯截断缝：带外梯度 0、带内非 0",
            xy=(1.0, data["gauss"]["g"][near(1.0)]),
            xytext=(0.04, 0.78), textcoords="axes fraction", color=MUTED,
            fontsize=8.5, arrowprops=dict(arrowstyle="->", color=MUTED, lw=1))
ax.set_title("C. 带边放大：线性梯度硬跳、五次归零、高斯留缝",
             fontsize=10.5, loc="left")
ax.set_xlabel("归一化距离 d/R"); ax.set_ylabel("梯度 dK/dd")
ax.set_xlim(-1.25, -0.7); style(ax)

tr = load("boundary_trace.csv")
xs = np.array([int(r["x"]) for r in tr])
sl = np.array([float(r["slope"]) for r in tr])
h = np.array([float(r["height"]) for r in tr])
stairs = load("stairs.csv")
stair_x = [int(s["x"]) for s in stairs]

ax = ax_t
for k in kernels:
    gw = np.array([float(r[f"gw_{k}"]) for r in tr])
    ax.plot(xs, gw, color=COLORS[k], lw=2, label=LABEL[k])
for sx_ in stair_x:
    ax.axvline(sx_, color="#d03b3b", lw=1.1, ls="--", zorder=1)
ax.annotate("陡坡台阶", xy=(stair_x[0], 0.42), color="#d03b3b", fontsize=8.5)
ax.set_title("D. 实际权重梯度 ∂w/∂x（y=32，中心差分）",
             fontsize=10.5, loc="left")
ax.set_xlabel("格点 x（y=32）"); ax.set_ylabel("∂w/∂x")
style(ax)

ax = ax_s
ax.fill_between(xs, 0, sl, color="#86b6ef", alpha=0.35, zorder=1)
ax.plot(xs, sl, color="#256abf", lw=1.6)
for sx_ in stair_x:
    ax.axvline(sx_, color="#d03b3b", lw=1.1, ls="--", zorder=2)
ax.set_title("E. 同剖面坡度 |∇h|（红虚线=台阶）",
             fontsize=10.5, loc="left")
ax.set_xlabel("格点 x（y=32）"); ax.set_ylabel("坡度 |∇h|")
style(ax)

ax = ax_h
ax.plot(xs, h, color=MUTED, lw=1.6)
ax.axhline(0, color="#c3c2b7", lw=1.0)
ax.fill_between(xs, h, 0, where=(h >= 0), color="#86b6ef", alpha=0.25)
for sx_ in stair_x:
    ax.axvline(sx_, color="#d03b3b", lw=1.1, ls="--", zorder=2)
ax.set_title("F. 同剖面高度 h 与边界 h0=0", fontsize=10.5, loc="left")
ax.set_xlabel("格点 x（y=32）"); ax.set_ylabel("高度 h")
style(ax)

handles, labels = ax_w.get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False,
           fontsize=10.5, bbox_to_anchor=(0.5, 0.995))
fig.subplots_adjust(left=0.06, right=0.975, bottom=0.08, top=0.90,
                    hspace=0.42, wspace=0.30)
fig.savefig(os.path.join(OUT, "fig_gradients.png"), dpi=150)
plt.close(fig)

# ================= figure 2: cost comparison =================
mk = load("kernel_metrics.csv")
names = [r["kernel"] for r in mk]
full = np.array([float(r["full_ns_per_eval"]) for r in mk])
konly = np.array([float(r["kernel_only_ns"]) for r in mk])
xpos = np.arange(len(names)); w = 0.36

fig2, ax = plt.subplots(figsize=(8.2, 4.8), facecolor=SURFACE)
b1 = ax.bar(xpos - w / 2, konly, w, color=[COLORS[k] for k in names],
            label="核纯计算（隔离微基准）", zorder=3)
b2 = ax.bar(xpos + w / 2, full, w, color=[COLORS[k] for k in names],
            alpha=0.45, label="全场遍历（含访存，实测）", zorder=3)
for bars in (b1, b2):
    for b in bars:
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.06,
                f"{b.get_height():.2f}", ha="center", fontsize=9,
                color=INK)
ax.set_xticks(xpos)
ax.set_xticklabels([LABEL[k] for k in names])
ax.set_ylabel("每次求值成本 (ns, -O2, 7 轮取最小值)")
ax.set_title("三种核的计算成本（同一高度场，4096 格点 ×400 次/轮）",
             fontsize=11, loc="left")
ax.legend(frameon=False, fontsize=9)
style(ax)
fig2.tight_layout()
fig2.savefig(os.path.join(OUT, "fig_cost.png"), dpi=150)
plt.close(fig2)

# ================= machine-readable summary =================
base = konly[names.index("linear")]
summary = {
    "kernels": [
        {"kernel": r["kernel"],
         "full_ns": float(r["full_ns_per_eval"]),
         "kernel_only_ns": float(r["kernel_only_ns"]),
         "cost_index_vs_linear": round(
             float(r["kernel_only_ns"]) / base, 3),
         "edge_gradient_jump": max(float(r["jump_minusR"]),
                                   float(r["jump_plusR"])),
         "max_abs_second_deriv": float(r["max_abs_second_deriv"]),
         "rms_second_deriv": float(r["rms_second_deriv"])}
        for r in mk],
    "stairs": [{"coord": f"({s['x']},{s['y']})",
                "slope": float(s["slope"]),
                "height_step": float(s["height_step"]),
                "weight_delta_fixedR": float(s["weight_delta_fixedR"]),
                "weight_delta_adaptiveR": float(s["weight_delta_adaptiveR"]),
                "adaptive_R": float(s["adaptive_R"])}
               for s in stairs],
}
with open(os.path.join(OUT, "summary.json"), "w") as fh:
    json.dump(summary, fh, ensure_ascii=False, indent=2)
print("wrote fig_gradients.png fig_cost.png summary.json ->", OUT)
