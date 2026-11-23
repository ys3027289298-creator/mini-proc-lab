#!/usr/bin/env python3
"""Render per-pixel conservation-audit difference maps from the postfx dump.

Usage: python3 -I render.py postfx.bin outdir

Reads the 'PFX1' bundle written by src/postfx.c (native-endian float32/double),
recomputes every difference in float64, and emits:
  diffmap_conserving.png   |measured - analytic|, one panel per channel
  diffmap_buggy.png        signed measured - analytic, one panel per channel
  stats.json               per-channel max abs/relative differences + verdict
  energy_table.md          three-channel energy change table (all stages)
"""
import json
import struct
import sys
import site
from pathlib import Path

# -I (isolated mode) suppresses the user site directory; add it back
# explicitly so numpy/matplotlib resolve while env/cwd isolation stays on.
sys.path.insert(0, site.getusersitepackages())

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# dataviz reference palette (light mode)
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_MUTED = "#52514e"
GRID = "#e1e0d9"
BLUE_600, BLUE_300, BLUE_100 = "#184f95", "#6da7ec", "#cde2fb"
RED_600, RED_300 = "#a11d1d", "#e66767"
MID_GRAY = "#f0efec"
CRITICAL = "#d03b3b"
GOOD = "#0ca30c"

CH_NAMES = ("R", "G", "B")
STAGE_NAMES = ("input", "gain", "lift", "overlay", "tint")
THRESH = 1e-3

plt.rcParams.update({
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "font.family": "sans-serif",
    "font.sans-serif": ["system-ui", "Segoe UI", "Ubuntu", "DejaVu Sans"],
    "text.color": INK,
    "axes.labelcolor": INK_MUTED,
    "xtick.color": INK_MUTED,
    "ytick.color": INK_MUTED,
    "axes.edgecolor": GRID,
})

seq_blue = LinearSegmentedColormap.from_list("seq_blue", [BLUE_100, BLUE_300, BLUE_600])
div_br = LinearSegmentedColormap.from_list("div_br", [BLUE_600, MID_GRAY, RED_600])


def load(path):
    buf = Path(path).read_bytes()
    magic, h, w, ch, stages, _ = struct.unpack_from("<6i", buf, 0)
    if magic != 0x50465831:
        raise SystemExit(f"bad magic {magic:#x}")
    n = h * w
    off = 24
    blocks = {}
    for name in ("input", "measured_ok", "measured_bug", "analytic_ok", "analytic_bug"):
        blocks[name] = np.frombuffer(buf, dtype="<f4", count=ch * n, offset=off).reshape(ch, h, w).astype(np.float64)
        off += ch * n * 4
    # stage tables: per stage 6 arrays of ch doubles (m2,a2,k2,m2b,a2b,k2b)
    tab = np.frombuffer(buf, dtype="<f8", offset=off).reshape(stages, 6, ch)
    return h, w, blocks, tab


def panel(ax, m, cmap, title):
    im = ax.imshow(m, origin="upper", interpolation="nearest",
                   cmap=cmap, norm=_lognorm(m))
    ax.set_title(title, fontsize=11, color=INK)
    ax.set_xticks([]); ax.set_yticks([])
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.ax.tick_params(labelsize=8, colors=INK_MUTED)
    cb.outline.set_edgecolor(GRID)


def _lognorm(m):
    from matplotlib.colors import LogNorm
    vmax = max(float(np.abs(m).max()), 1e-12)
    vmin = max(vmax * 1e-6, 1e-18)
    return LogNorm(vmin=vmin, vmax=vmax)


def main():
    binpath, outdir = sys.argv[1], Path(sys.argv[2])
    outdir.mkdir(parents=True, exist_ok=True)
    h, w, b, tab = load(binpath)

    stats = {"size": [h, w], "threshold": THRESH, "conserving": {}, "buggy": {}}

    # ---- figure 1: conserving stack, absolute differences ----------------
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4), constrained_layout=True)
    for c, (ax, name) in enumerate(zip(axes, CH_NAMES)):
        d = np.abs(b["measured_ok"][c] - b["analytic_ok"][c])
        denom = np.maximum(1e-6, np.abs(b["analytic_ok"][c]))
        rel = d / denom
        rmax = float(rel.max()); amax = float(d.max())
        verdict = "PASS" if rmax <= THRESH else "FAIL"
        stats["conserving"][name] = {"max_abs": amax, "max_rel": rmax, "verdict": verdict}
        panel(ax, d, seq_blue,
              f"{name}  |measured − analytic|\nmax|Δ|={amax:.2e}   max rel={rmax:.2e}   {verdict}")
        ax.text(0.02, 0.02, f"limit {THRESH:.0e} per pixel", transform=ax.transAxes,
                fontsize=8, color=GOOD if verdict == "PASS" else CRITICAL, va="bottom")
    fig.suptitle("Conserving stack — per-pixel |measured − analytic| (log scale, blue = magnitude)",
                 fontsize=13, color=INK)
    fig.savefig(outdir / "diffmap_conserving.png", dpi=150)
    plt.close(fig)

    # ---- figure 2: buggy stack, signed differences ------------------------
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4), constrained_layout=True)
    raw = [b["measured_bug"][c] - b["analytic_bug"][c] for c in range(3)]
    for c, (ax, name, d) in enumerate(zip(axes, CH_NAMES, raw)):
        denom = np.maximum(1e-6, np.abs(b["analytic_bug"][c]))
        rmax = float((np.abs(d) / denom).max()); amax = float(np.abs(d).max())
        verdict = "DETECTED" if rmax > THRESH else "clean"
        stats["buggy"][name] = {"max_abs": amax, "max_rel": rmax, "verdict": verdict}
        # per-panel scale: R/B float noise (~1e-7) must stay visible instead
        # of being flattened by G's ~1 scale on a shared colorbar
        vlim = max(amax, 1e-12)
        im = ax.imshow(d, origin="upper", interpolation="nearest",
                       cmap=div_br, vmin=-vlim, vmax=vlim)
        ax.set_title(f"{name}  measured − analytic\nmax|Δ|={amax:.2e}   max rel={rmax:.2e}   {verdict}",
                     fontsize=11, color=CRITICAL if verdict == "DETECTED" else INK)
        ax.set_xticks([]); ax.set_yticks([])
        cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cb.ax.tick_params(labelsize=8, colors=INK_MUTED)
        cb.outline.set_edgecolor(GRID)
        ax.text(0.02, 0.02, f"limit {THRESH:.0e} per pixel", transform=ax.transAxes,
                fontsize=8, color=CRITICAL if verdict == "DETECTED" else GOOD, va="bottom")
    fig.suptitle("Buggy stack — signed per-pixel difference, per-channel scale (blue = measured below analytic = swallowed)",
                 fontsize=13, color=INK)
    fig.savefig(outdir / "diffmap_buggy.png", dpi=150)
    plt.close(fig)

    # ---- energy table (markdown) ------------------------------------------
    lines = [
        "# 三通道 L2 能量变化表（E2 = Σ pixel²，N=%d 像素/通道）" % (h * w), "",
        "## 守恒栈（tint 各通道均为 a·in+b）", "",
        "| 阶段 | 通道 | 实测 E2 | 解析 E2 | 簿记 E2 | 实测vs解析 相对差 | 实测vs簿记 相对差 |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for s, sname in enumerate(STAGE_NAMES):
        for c, cn in enumerate(CH_NAMES):
            m2, a2, k2 = tab[s, 0, c], tab[s, 1, c], tab[s, 2, c]
            lines.append(f"| {sname} | {cn} | {m2:.6f} | {a2:.6f} | {k2:.6f} | "
                         f"{_rel(m2, a2):.2e} | {_rel(m2, k2):.2e} |")
    lines += ["", "## 缺陷栈（tint 的 G 通道被常量覆盖）", "",
              "| 阶段 | 通道 | 实测 E2 | 解析 E2 | 簿记 E2 | 实测vs解析 相对差 | 判定 |",
              "|---|---|---:|---:|---:|---:|---|"]
    for s, sname in enumerate(STAGE_NAMES):
        for c, cn in enumerate(CH_NAMES):
            m2, a2, k2 = tab[s, 3, c], tab[s, 4, c], tab[s, 5, c]
            r = _rel(m2, a2)
            verdict = "超限（吞通道）" if r > THRESH else "通过"
            lines.append(f"| {sname} | {cn} | {m2:.6f} | {a2:.6f} | {k2:.6f} | {r:.2e} | {verdict} |")
    lines += ["", "## 叠加前后能量变化（input → tint）", "",
              "| 栈 | 通道 | E2 输入 | E2 输出（实测） | 倍率 | 解析输出 | 相对差 |",
              "|---|---|---:|---:|---:|---:|---:|"]
    for tag, base in (("守恒", 0), ("缺陷", 3)):
        for c, cn in enumerate(CH_NAMES):
            e0, e1, ea = tab[0, base, c], tab[4, base, c], tab[4, base + 1, c]
            lines.append(f"| {tag} | {cn} | {e0:.4f} | {e1:.4f} | {e1 / e0:.5f}× | {ea:.4f} | {_rel(e1, ea):.2e} |")
    (outdir / "energy_table.md").write_text("\n".join(lines) + "\n")

    stats["all_within_threshold"] = all(v["verdict"] == "PASS" for v in stats["conserving"].values())
    stats["buggy_channel"] = [k for k, v in stats["buggy"].items() if v["verdict"] == "DETECTED"]
    (outdir / "stats.json").write_text(json.dumps(stats, indent=2))
    print(json.dumps(stats, indent=2))


def _rel(a, b):
    d = abs(float(a) - float(b))
    return d / max(abs(float(a)), abs(float(b)), 1e-300) if d > 0 else 0.0


if __name__ == "__main__":
    main()
