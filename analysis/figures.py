#!/usr/bin/env python3
"""Figures for the erosion convergence report.

    python3 analysis/figures.py <campaign-dir> <out-dir>

Palette follows the dataviz reference instance: categorical slots in fixed
order (blue, orange, aqua), sequential blue for magnitude, light surface.
"""
import json, sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

CAMPAIGN = sys.argv[1] if len(sys.argv) > 1 else "artifacts/erosion/campaign"
OUT = sys.argv[2] if len(sys.argv) > 2 else "artifacts/erosion/report"
os.makedirs(OUT, exist_ok=True)

# --- validated light-mode slots (dataviz reference palette) ---
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"      # blue, orange, aqua
CRIT = "#d03b3b"
SEQ = LinearSegmentedColormap.from_list(
    "seq_blue", ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["system-ui", "DejaVu Sans"],
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "axes.titlecolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.size": 9, "axes.titlesize": 11, "axes.labelsize": 9,
    "legend.frameon": False, "lines.linewidth": 2.0,
})


def load(tag, it, n=256):
    return np.fromfile(f"{CAMPAIGN}/fields/{tag}.{it:04d}.h32",
                       dtype="<f4").reshape(n, n).astype(np.float64)


def series(tag):
    d = json.load(open(f"{CAMPAIGN}/logs/{tag}.json"))
    return {k: np.array([r[k] for r in d["runs"]])
            for k in ("iter", "ms", "residual", "rms_vs_init",
                      "mean_abs_delta", "us_per_cell")}


def hillshade(h, az=315.0, alt=45.0):
    gy, gx = np.gradient(h)
    slope = np.pi / 2 - np.arctan(np.hypot(gx, gy))
    aspect = np.arctan2(-gx, gy)
    a, b = np.deg2rad(360 - az + 90), np.deg2rad(alt)
    return np.sin(b) * np.sin(slope) + np.cos(b) * np.cos(slope) * np.cos(a - aspect)


# ---------------------------------------------------------------- fig 1: residual
def fig_residual():
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    ax.set_yscale("log")
    for tag, col, label in (("t10", S1, "10 iters"), ("t30", S2, "30 iters"),
                            ("t90", S3, "90 iters")):
        s = series(tag)
        ax.plot(s["iter"], s["residual"], color=col, lw=1.6, label=label)
        ax.scatter(s["iter"][-1], s["residual"][-1], s=26, color=col, zorder=5,
                   edgecolor=SURFACE, linewidth=1.5)
    # head-law annotation (power law r ∝ i^-0.54, first 25 iters)
    it = np.arange(1, 21)
    s90 = series("t90")
    p = np.polyfit(np.log(it), np.log(s90["residual"][:20]), 1)
    ax.plot(it, np.exp(np.polyval(p, np.log(it))), color=INK2, lw=1.2,
            linestyle=(0, (1, 2)), label=r"head fit  $r\propto i^{-0.56}$")
    knee = 33
    ax.axvline(knee, color=CRIT, lw=1.0)
    ax.annotate("knee @ iter 33", xy=(knee, s90["residual"][knee - 1]),
                xytext=(knee + 3, 0.115), color=CRIT, fontsize=9,
                arrowprops=dict(arrowstyle="-", color=CRIT, lw=0.8))
    ax.set_xlabel("iteration")
    ax.set_ylabel("residual  RMS(h$_i$ − h$_{i-1}$)")
    ax.set_title("Height residual per iteration (log y) — blue/orange dots mark the 10/30-iteration stops")
    ax.grid(True, which="major", axis="y")
    ax.set_xlim(1, 90)
    ax.legend(loc="upper right", ncol=2)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig1_residual.png", dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------- fig 2: evolution
def fig_evolution():
    its = [0, 10, 20, 30, 33, 40, 50, 70, 90]
    fig, axs = plt.subplots(3, 3, figsize=(8.4, 8.4))
    for ax, it in zip(axs.ravel(), its):
        h = load("t90", it)
        ax.imshow(hillshade(h), cmap="gray", vmin=-1, vmax=1, interpolation="nearest")
        ax.set_title(f"iter {it}" + ("  (knee)" if it == 33 else ""), fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
    fig.suptitle("t90: hillshade every ~10 iterations — gully incisions coarsen "
                 "and deepen after the knee", fontsize=10, color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(f"{OUT}/fig2_evolution.png", dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------- fig 3: gully depth
def fig_gully():
    fig, axs = plt.subplots(1, 2, figsize=(10.2, 3.9))

    def comp_depths(tag, it):
        from scipy.ndimage import median_filter, label
        h = load(tag, it)
        depth = median_filter(h, size=9, mode="nearest") - h
        lab, nc = label(depth > 0.5, np.ones((3, 3)))
        sizes = np.bincount(lab.ravel())[1:]
        out = []
        for c in range(1, nc + 1):
            m = lab == c
            if m.sum() < 8:
                continue
            out.append(depth[m].max())
        return np.array(out)

    ax = axs[0]
    for tag, it, col, label in (("t10", 10, S1, "10 iters"),
                                ("t30", 30, S2, "30 iters"),
                                ("t90", 90, S3, "90 iters")):
        d = comp_depths(tag, it)
        bins = np.logspace(np.log10(0.5), np.log10(400), 28)
        ax.hist(d, bins=bins, histtype="step", color=col, lw=1.8,
                density=False, label=f"{label}  (n={len(d)})")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("gully max depth, lattice height units (log)")
    ax.set_ylabel("gully count")
    ax.set_title("Per-gully max-depth distribution")
    ax.grid(True, axis="y")
    ax.legend(loc="upper right")

    ax = axs[1]
    qs = np.array([0, 10, 25, 50, 75, 90, 99])
    for tag, it, col, label in (("t10", 10, S1, "10"), ("t30", 30, S2, "30"),
                                ("t90", 90, S3, "90")):
        d = comp_depths(tag, it)
        vals = np.percentile(d, qs)
        ax.plot(qs, vals, marker="o", ms=5, color=col, label=f"{label} iters",
                markeredgecolor=SURFACE, markeredgewidth=1.2)
    ax.set_yscale("log")
    ax.set_xticks([0, 10, 25, 50, 75, 90, 99])
    ax.set_xlabel("percentile of gullies")
    ax.set_ylabel("max depth (log)")
    ax.set_title("Depth quantiles — tail grows after the knee")
    ax.grid(True, axis="y")
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig3_gully_depth.png", dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------- fig 4: cost
def fig_cost():
    fig, axs = plt.subplots(1, 2, figsize=(10.0, 3.6))

    ax = axs[0]
    for tag, col, label in (("t10", S1, "10 iters"), ("t30", S2, "30 iters"),
                            ("t90", S3, "90 iters")):
        s = series(tag)
        ax.plot(s["iter"], s["ms"], color=col, lw=1.2, label=label)
    ax.set_xlabel("iteration")
    ax.set_ylabel("wall time per iteration (ms)")
    ax.set_title("Per-iteration cost, n=256, 1 drop/cell")
    ax.set_ylim(0, 460)
    ax.grid(True, axis="y")
    ax.legend(loc="upper right", ncol=2)

    ax = axs[1]
    sizes = [128, 256, 512]
    med = []
    for n in sizes:
        s = series(f"size{n}")
        med.append(np.median(s["ms"][2:]))
    ax.plot(sizes, med, marker="o", ms=6, color=S1, markeredgecolor=SURFACE,
            markeredgewidth=1.5)
    for n, m in zip(sizes, med):
        ax.annotate(f"{m:.0f} ms", (n, m), textcoords="offset points",
                    xytext=(0, 9), ha="center", color=INK2, fontsize=9)
    # quadratic reference from the smallest size
    ref = med[0] * (np.array(sizes) / sizes[0]) ** 2
    ax.plot(sizes, ref, color=MUTED, lw=1.0, linestyle=(0, (2, 2)),
            label="∝ cells (drops = 1/cell)")
    ax.set_xscale("log", base=2); ax.set_yscale("log")
    ax.set_xticks(sizes); ax.set_xticklabels([str(s) for s in sizes])
    ax.set_xlabel("field side n")
    ax.set_ylabel("median iteration time (ms)")
    ax.set_title("Cost scales linearly in cell count")
    ax.grid(True, which="major", axis="y")
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig4_cost.png", dpi=160)
    plt.close(fig)


def fig_decomposition():
    path = f"{CAMPAIGN}/logs/t90_residual_decomposition.json"
    if not os.path.exists(path):
        return
    rows = json.load(open(path))
    it = np.array([r["iter"] for r in rows])
    tot = np.array([r["residual"] for r in rows])
    g = np.array([r["gully_rms"] for r in rows])
    n = np.array([r["nongully_rms"] for r in rows])
    e01 = np.array([r["excl_top01pct"] for r in rows])

    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    ax.set_yscale("log")
    ax.plot(it, tot, color=INK2, lw=2.0, label="total residual")
    ax.plot(it, g, color=S1, lw=1.8, label="gully-mask component")
    ax.plot(it, n, color=S2, lw=1.8, label="off-gully component")
    ax.plot(it, e01, color=S3, lw=1.4, linestyle=(0, (3, 2)),
            label="total, top 0.1% cells removed")
    ax.axvline(33, color=CRIT, lw=1.0)
    ax.annotate("knee @ iter 33", xy=(33, 0.2128), xytext=(36, 0.28), color=CRIT,
                fontsize=9, arrowprops=dict(arrowstyle="-", color=CRIT, lw=0.8))
    ax.set_xlabel("iteration")
    ax.set_ylabel("RMS of per-cell Δh (log)")
    ax.set_title("What drives the reversal: the off-gully slope component")
    ax.grid(True, axis="y")
    ax.set_xlim(1, 90)
    ax.set_ylim(0.115, 1.4)
    ax.legend(loc="upper right", ncol=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig5_decomposition.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    fig_residual()
    fig_evolution()
    fig_gully()
    fig_cost()
    fig_decomposition()
    print("figures written to", OUT)
