#!/usr/bin/env python3
"""Analyze fBm octave-boundary derivative jumps and render LOD comparisons.

Usage: python3 analyze.py <data_dir> <out_dir> [threshold_deg]

Inputs (produced by the C lod binary): meta.json, *_h.bin, *_gx.bin,
*_gy.bin, boundary_jumps.csv.
Outputs (in out_dir): stats CSVs, summary.json, and four figures.
"""
import json
import sys
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# dataviz reference palette (light surface)
INK = "#0b0b0b"
INK2 = "#52514e"
BLUE = "#2a78d6"
ORANGE = "#eb6834"
RED = "#e34948"
SURFACE = "#fcfcfb"
GRID = "#d9d8d3"

plt.rcParams.update({
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "axes.edgecolor": GRID,
    "axes.labelcolor": INK,
    "axes.titlecolor": INK,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "text.color": INK,
    "font.size": 10,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.7,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

Z = 3.0          # vertical exaggeration applied consistently to both fields
AZ = 315.0       # hillshade azimuth (deg), light from NW
ALT = 45.0       # hillshade altitude (deg)


def load(data_dir):
    with open(os.path.join(data_dir, "meta.json")) as f:
        meta = json.load(f)
    nx, ny = meta["nx"], meta["ny"]
    grids = {}
    for name in meta["arrays"]:
        grids[name] = np.fromfile(
            os.path.join(data_dir, name + ".bin"), dtype="<f4"
        ).reshape(ny, nx).astype(np.float64)
    xs = np.linspace(meta["xmin"], meta["xmax"], nx)
    ys = np.linspace(meta["ymin"], meta["ymax"], ny)
    csv = np.genfromtxt(os.path.join(data_dir, "boundary_jumps.csv"),
                        delimiter=",", names=True, dtype=None, encoding="ascii")
    return meta, grids, xs, ys, csv


def normal_angle(gx, gy):
    """Angle in degrees between normals of x-adjacent pixels.

    Normal at a pixel with analytic gradient (gx, gy): n=(-Z gx,-Z gy,1).
    Returned array has shape (ny, nx-1); column j holds the pair (j, j+1).
    """
    nx1, ny1 = -Z * gx[:, :-1], -Z * gy[:, :-1]
    nx2, ny2 = -Z * gx[:, 1:], -Z * gy[:, 1:]
    d1 = np.sqrt(nx1**2 + ny1**2 + 1.0)
    d2 = np.sqrt(nx2**2 + ny2**2 + 1.0)
    dot = (nx1 * nx2 + ny1 * ny2 + 1.0) / (d1 * d2)
    return np.degrees(np.arccos(np.clip(dot, -1.0, 1.0)))


def tier_masks(xs, meta):
    b1, b2 = meta["boundaries"]
    return {
        f"{meta['tiers'][0]} oct": xs < b1,
        f"{meta['tiers'][1]} oct": (xs >= b1) & (xs < b2),
        f"{meta['tiers'][2]} oct": xs >= b2,
    }


def tier_stats(grids, xs, meta):
    """Per-tier normal-angle statistics over all in-tier adjacent pairs."""
    masks = tier_masks(xs, meta)
    rows = []
    for field in ("hard", "smooth"):
        ang = normal_angle(grids[f"{field}_gx"], grids[f"{field}_gy"])
        # ang columns correspond to the left pixel of each x-adjacent pair
        xc = 0.5 * (xs[:-1] + xs[1:])
        for tier, m in masks.items():
            # a pair counts only when both of its pixels sit in this tier
            sel = m[:-1] & m[1:]
            a = ang[:, sel]
            kc = xc[sel]
            assert kc.size > 0
            rows.append({
                "field": field, "tier": tier,
                "pairs": int(a.size),
                "mean_deg": float(a.mean()),
                "p50_deg": float(np.percentile(a, 50)),
                "p95_deg": float(np.percentile(a, 95)),
                "p99_deg": float(np.percentile(a, 99)),
                "max_deg": float(a.max()),
                "frac_over_1deg": float((a > 1.0).mean()),
                "frac_over_5deg": float((a > 5.0).mean()),
            })
    return rows


def boundary_stats(csv, grids, xs, meta, trows):
    """Jump stats at each boundary in two metrics.

    analytic: exact one-sided gradient/height jumps at the cut (from the C
              side, CSV), both fields.
    render:   normal angles of grid-adjacent pixels. The hard field is
              evaluated on its single cut pair; the smooth field over every
              pair inside its blend band. Two render numbers are reported:
              the raw max angle, and the "excess over natural roughness":
              the p95 angle of the same octave tier measured well away from
              any boundary, i.e. the angle the terrain would have anyway.
    """
    rows = []
    dx = xs[1] - xs[0]
    ang = {f: normal_angle(grids[f"{f}_gx"], grids[f"{f}_gy"])
           for f in ("hard", "smooth")}
    # natural roughness per tier: p95 adjacent angle over ALL in-tier pairs
    # of the smooth field (tier stats already exclude cross-boundary pairs).
    # This is the angle the octave content produces without any transition.
    natural_p95 = {}
    masks = tier_masks(xs, meta)
    tmap = {(r["field"], r["tier"]): r["p95_deg"] for r in trows}
    for tier in masks:
        natural_p95[tier] = tmap[("smooth", tier)]
    for b in meta["boundaries"]:
        m = csv["boundary"] == b
        icut = int(np.argmin(np.abs(xs - b)))
        for kind, dh, dgx, dgy, an in (
            ("hard", "hard_dh", "hard_dgx", "hard_dgy", "hard_angle_deg"),
            ("smooth", "smooth_dh", "smooth_dgx", "smooth_dgy",
             "smooth_angle_deg"),
        ):
            if kind == "hard":
                band = np.array([icut - 1])      # single crossing pair
            else:
                lo = int(np.argmin(np.abs(xs - (b - meta["band_halfwidth"]))))
                hi = int(np.argmin(np.abs(xs - (b + meta["band_halfwidth"]))))
                band = np.arange(max(lo - 1, 0), min(hi, xs.size - 1))
            a = ang[kind][:, band]
            # tier above this boundary supplies the natural-roughness baseline
            above_tier = f"{b + 1} oct"
            base = natural_p95[above_tier]
            rows.append({
                "boundary": int(b), "field": kind,
                "height_step_max_abs": float(np.max(np.abs(csv[dh][m]))),
                "grad_jump_rms": float(np.sqrt(
                    csv[dgx][m] ** 2 + csv[dgy][m] ** 2).mean()),
                "grad_jump_max": float(np.max(
                    np.sqrt(csv[dgx][m] ** 2 + csv[dgy][m] ** 2))),
                "angle_mean_deg": float(csv[an][m].mean()),
                "angle_p95_deg": float(np.percentile(csv[an][m], 95)),
                "angle_max_deg": float(csv[an][m].max()),
                "render_max_adj_angle_deg": float(a.max()),
                "render_p95_adj_angle_deg": float(np.percentile(a, 95)),
                "natural_p95_adj_angle_deg": base,
                "render_excess_p95_over_natural_deg":
                    float(max(0.0, np.percentile(a, 95) - base)),
                "render_pairs": int(a.size),
                "render_dx_world": float(dx),
            })
    return rows


def hillshade(h, dx):
    gy, gx = np.gradient(h, dx)
    gx, gy = Z * gx, Z * gy
    az = np.radians(AZ)
    alt = np.radians(ALT)
    slope = np.pi / 2.0 - np.arctan(np.hypot(gx, gy))
    aspect = np.arctan2(-gy, gx)
    shade = (np.sin(alt) * np.sin(slope)
             + np.cos(alt) * np.cos(slope) * np.cos(az - aspect))
    p1, p99 = np.percentile(shade, (1, 99))
    return np.clip((shade - p1) / (p99 - p1), 0, 1)


def fig_hillshades(meta, grids, xs, ys, out):
    fig, axes = plt.subplots(2, 1, figsize=(11, 6.2), sharex=True, sharey=True)
    dx = xs[1] - xs[0]
    for ax, field, title in (
        (axes[0], "hard", "Hard octave cutoff (C0 step + derivative crease)"),
        (axes[1], "smooth", "Smooth octave blend (C1 transition, band half-width "
                            f"{meta['band_halfwidth']})"),
    ):
        ax.imshow(hillshade(grids[f"{field}_h"], dx), origin="lower",
                  cmap="gray", extent=[xs[0], xs[-1], ys[0], ys[-1]],
                  aspect="auto")
        for b in meta["boundaries"]:
            ax.axvline(b, color=RED, lw=1.4, ls="--")
        ax.set_title(title)
        ax.set_ylabel("y")
        for t, x0, x1 in ((f"{k} oct", lo, hi) for (k, lo, hi) in
                          zip(meta["tiers"],
                              [xs[0]] + list(meta["boundaries"]),
                              list(meta["boundaries"]) + [xs[-1]])):
            ax.text((x0 + x1) / 2, ys[-1] * 0.93, t, ha="center",
                    color=INK, fontsize=9,
                    bbox=dict(facecolor=SURFACE, edgecolor=GRID, alpha=0.8))
    axes[1].set_xlabel("x (octave count changes at red lines)")
    legend = [Patch(facecolor="none", edgecolor=RED, ls="--",
                    label="octave boundary")]
    axes[0].legend(handles=legend, loc="lower left", framealpha=0.9)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig1_hillshades.png"), dpi=160)
    plt.close(fig)


def fig_boundary_profiles(csv, meta, out, threshold):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, b in zip(axes, meta["boundaries"]):
        m = csv["boundary"] == b
        y = csv["y"][m]
        ax.plot(csv["hard_angle_deg"][m], y, color=RED, lw=1.6,
                label="hard cutoff")
        ax.plot(csv["smooth_angle_deg"][m], y, color=BLUE, lw=1.6,
                label="smooth blend")
        ax.axvline(threshold, color=INK2, lw=1.2, ls=":")
        ax.set_title(f"Boundary {b}→{b + 1} octaves (x = {b})")
        ax.set_xlabel("one-sided normal angle at the cut (deg)")
        ax.set_xlim(left=-1)
    axes[0].set_ylabel("y along boundary")
    axes[0].legend(framealpha=0.9)
    axes[0].annotate(f"threshold {threshold:g}°", xy=(threshold, 0.97),
                     xycoords=("data", "axes fraction"), color=INK2,
                     fontsize=8, ha="left", va="top")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig2_boundary_angle_profiles.png"), dpi=160)
    plt.close(fig)


def fig_tier_bars(trows, out):
    tiers = sorted({r["tier"] for r in trows})
    fields = ("hard", "smooth")
    x = np.arange(len(tiers))
    w = 0.36
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, key, title in (
        (axes[0], "mean_deg", "Mean adjacent-normal angle per tier"),
        (axes[1], "p99_deg", "p99 adjacent-normal angle per tier"),
    ):
        for i, fld in enumerate(fields):
            vals = [next(r[key] for r in trows
                         if r["field"] == fld and r["tier"] == t) for t in tiers]
            ax.bar(x + (i - 0.5) * w, vals, w,
                   color=RED if fld == "hard" else BLUE,
                   label=fld)
        ax.set_xticks(x)
        ax.set_xticklabels(tiers)
        ax.set_ylabel("degrees")
        ax.set_title(title)
        ax.legend(framealpha=0.9)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig3_tier_angle_stats.png"), dpi=160)
    plt.close(fig)


def fig_gradient_profile(grids, xs, ys, meta, out):
    j = np.argmin(np.abs(ys - 0.0))
    fig, axes = plt.subplots(2, 1, figsize=(11, 6.0), sharex=True)
    for ax, key, title in (
        (axes[0], "gx", "Analytic x-derivative ∂h/∂x along y = 0"),
        (axes[1], "h", "Height h along y = 0"),
    ):
        ax.plot(xs, grids[f"hard_{key}"][j], color=RED, lw=1.3,
                label="hard cutoff")
        ax.plot(xs, grids[f"smooth_{key}"][j], color=BLUE, lw=1.3,
                label="smooth blend")
        for b in meta["boundaries"]:
            ax.axvline(b, color=INK2, lw=0.9, ls="--")
        for b in meta["boundaries"]:
            ax.axvspan(b - meta["band_halfwidth"], b + meta["band_halfwidth"],
                       color=ORANGE, alpha=0.15)
        ax.set_ylabel(key)
        ax.set_title(title)
        ax.legend(framealpha=0.9, loc="upper right")
    axes[1].set_xlabel("x")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig4_derivative_profile.png"), dpi=160)
    plt.close(fig)


def write_csv(path, rows):
    keys = list(rows[0].keys())
    with open(path, "w") as f:
        f.write(",".join(keys) + "\n")
        for r in rows:
            f.write(",".join(str(r[k]) for k in keys) + "\n")


def main():
    if len(sys.argv) < 3:
        sys.exit("usage: analyze.py <data_dir> <out_dir> [threshold_deg]")
    data_dir, out_dir = sys.argv[1], sys.argv[2]
    threshold = float(sys.argv[3]) if len(sys.argv) > 3 else 5.0
    os.makedirs(out_dir, exist_ok=True)

    meta, grids, xs, ys, csv = load(data_dir)
    trows = tier_stats(grids, xs, meta)
    brows = boundary_stats(csv, grids, xs, meta, trows)
    write_csv(os.path.join(out_dir, "tier_angle_stats.csv"), trows)
    write_csv(os.path.join(out_dir, "boundary_jump_stats.csv"), brows)

    # pass/fail on the excess that the transition adds beyond the terrain's
    # natural roughness (p95 band angle minus p95 far-from-boundary angle).
    # The hard cut is also reported raw for reference; its excess is the
    # crease it genuinely introduces.
    def excess(kind):
        return max(r["render_excess_p95_over_natural_deg"]
                   for r in brows if r["field"] == kind)
    hard_excess = excess("hard")
    smooth_excess = excess("smooth")
    worst = max((r for r in brows if r["field"] == "hard"),
                key=lambda r: r["grad_jump_max"])
    verdict = {
        "threshold_deg": threshold,
        "metric": "p95 normal-angle excess in transition over natural roughness",
        "hard_excess_deg": hard_excess,
        "smooth_excess_deg": smooth_excess,
        "hard_passes": bool(hard_excess <= threshold),
        "smooth_passes": bool(smooth_excess <= threshold),
        "largest_jump_boundary": worst["boundary"],
        "largest_jump_grad_max": worst["grad_jump_max"],
        "largest_jump_analytic_angle_max_deg": worst["angle_max_deg"],
        "task_complete": bool(smooth_excess <= threshold),
    }
    with open(os.path.join(out_dir, "summary.json"), "w") as f:
        json.dump({"meta": meta, "verdict": verdict,
                   "tier_stats": trows, "boundary_stats": brows}, f, indent=2)

    fig_hillshades(meta, grids, xs, ys, out_dir)
    fig_boundary_profiles(csv, meta, out_dir, threshold)
    fig_tier_bars(trows, out_dir)
    fig_gradient_profile(grids, xs, ys, meta, out_dir)

    print(json.dumps(verdict, indent=2))
    print("\ntier stats:")
    for r in trows:
        print(f"  {r['field']:6s} {r['tier']:7s} mean={r['mean_deg']:.4f}  "
              f"p95={r['p95_deg']:.4f}  p99={r['p99_deg']:.4f}  max={r['max_deg']:.4f}")
    print("\nboundary stats (analytic one-sided + render adjacent-pixel):")
    for r in brows:
        print(f"  b={r['boundary']} {r['field']:6s} "
              f"dh_max={r['height_step_max_abs']:.5f} "
              f"|dgrad|_rms={r['grad_jump_rms']:.5f} "
              f"analytic_max={r['angle_max_deg']:.3f}deg "
              f"render_p95={r['render_p95_adj_angle_deg']:.3f}deg "
              f"natural_p95={r['natural_p95_adj_angle_deg']:.3f}deg "
              f"excess={r['render_excess_p95_over_natural_deg']:.3f}deg")


if __name__ == "__main__":
    main()
