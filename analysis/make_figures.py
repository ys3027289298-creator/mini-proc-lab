"""Figures from results.json + regenerated fields:
  fig_contours_before.png / fig_contours_after.png -- each basis on one panel,
    three iso-level ribbons, width table beneath;
  fig_gn_law.png   -- gradient RMS vs octaves with the derived law overlaid;
  fig_spectra.png  -- radial PSD power laws with derived/measured slopes.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

from noise_bases import fbm_partials

HERE = Path(__file__).resolve().parent
FIGS = HERE / "figs"
R = json.load(open(HERE / "results.json"))

BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e1e0d9"
SURFACE = "#fcfcfb"
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Noto Sans CJK SC", "Liberation Sans"],
    "font.size": 10,
    "text.color": INK, "axes.edgecolor": "#c3c2b7", "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
})

g = R["grid"]
N, DOMAIN, DELTA = g["N"], g["domain"], g["delta"]
LEVELS = g["levels"]
CELL = DOMAIN / N
rng = np.linspace(-DOMAIN / 2, DOMAIN / 2, N, dtype=np.float32)
X, Y = np.meshgrid(rng, rng, indexing="ij")
MAX_OCT = max(map(int, R["G_rms"]["value"].keys()))

# fields needed for the four panels (only render what is shown)
need = {("value", R["grid"]["ref_oct"]),
        ("gradient", R["grid"]["mismatch_oct"]),
        ("gradient", R["alignment"]["n_star_int"])}
fields = {}
for basis in ("value", "gradient"):
    partials = fbm_partials(X, Y, MAX_OCT, g["lacunarity"], g["gain"], 11, basis)
    for n, F in enumerate(partials, 1):
        if (basis, n) in need:
            fields[(basis, n)] = F

# downsampled view (2048 -> 1024) keeps PNG size sane; ribbons drawn in pixels
# of the full grid, so widths on screen are scaled by the same factor everywhere
DS = N // 1024


def hillshade(F):
    gx = np.gradient(F.astype(np.float64), CELL, axis=0)
    gy = np.gradient(F.astype(np.float64), CELL, axis=1)
    slope = np.hypot(gx, gy)
    shade = 1 - np.clip(slope / np.quantile(slope, 0.98), 0, 1) * 0.25
    return shade[::DS, ::DS]


def ribbon_masks(F):
    Fv = F[::DS, ::DS]
    return [np.abs(Fv - np.float32(lv)) <= np.float32(DELTA) for lv in LEVELS]


def contour_figure(scenario, panels, out_name, title):
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 7.2))
    extent = [-DOMAIN / 2, DOMAIN / 2, -DOMAIN / 2, DOMAIN / 2]
    for ax, (basis, n, color) in zip(axes, panels):
        F = fields[(basis, n)]
        ax.imshow(hillshade(F), cmap="gray", origin="lower",
                  extent=extent, vmin=0.7, vmax=1.0, interpolation="nearest")
        for lv, mask in zip(LEVELS, ribbon_masks(F)):
            # paint the ribbon itself (|F-level| <= delta), not its outline
            overlay = np.zeros((*mask.shape, 4))
            rgb = matplotlib.colors.to_rgb(color)
            overlay[mask] = (*rgb, 0.85)
            ax.imshow(overlay, origin="lower", extent=extent, interpolation="nearest")
        ax.set_title(f"{basis} noise · fBm N={n}", color=INK, fontsize=11)
        ax.set_xticks([]); ax.set_yticks([])
        ax.grid(False)
    fig.suptitle(title, fontsize=13, color=INK, y=0.98)
    fig.tight_layout(rect=(0, 0.27, 1, 0.95))

    # width table under the panels
    rows = [r for r in R["_width_rows"] if r["scenario"] == scenario]
    col_labels = ["basis", "N", "level", "width (px)", "width (world)",
                  "local |grad|", "2δ / (w·|grad|)"]
    cell_text = []
    for r in rows:
        # identity check w = 2 delta / |grad F|, using the gradient sampled at
        # the same contour points the width was measured at
        ratio = (2 * DELTA) / (r["width_world"] * r["grad_on_contour"])
        cell_text.append([
            r["basis"], r["octaves"], f"{r['level']:+.2f}",
            f"{r['width_px']:.2f}", f"{r['width_world']:.4f}",
            f"{r['grad_on_contour']:.3f}", f"{ratio:.3f}",
        ])
    tbl = fig.add_axes([0.05, 0.02, 0.90, 0.21])
    tbl.axis("off")
    t = tbl.table(cellText=cell_text, colLabels=col_labels, loc="upper center",
                  cellLoc="center")
    t.auto_set_font_size(False); t.set_fontsize(8.5); t.scale(1, 1.3)
    for (row, col), cell in t.get_celld().items():
        cell.set_edgecolor(GRID)
        if row == 0:
            cell.set_facecolor("#f0efec"); cell.set_text_props(color=INK)
        else:
            cell.set_facecolor(SURFACE); cell.set_text_props(color=MUTED)
    fig.savefig(FIGS / out_name, dpi=160)
    plt.close(fig)


# attach parsed widths for table use
import csv
with (HERE / "widths.csv").open() as fh:
    R["_width_rows"] = list(csv.DictReader(fh))
for r in R["_width_rows"]:
    for k in ("octaves", "n_samples"):
        r[k] = int(float(r[k]))
    for k in ("level", "delta", "width_px", "width_iqr", "width_world",
              "grad_on_contour"):
        r[k] = float(r[k])

contour_figure(
    "before",
    [("value", R["grid"]["ref_oct"], BLUE),
     ("gradient", R["grid"]["mismatch_oct"], ORANGE)],
    "fig_contours_before.png",
    "Before alignment — same fixed level band δ, different octave counts",
)
contour_figure(
    "after",
    [("value", R["grid"]["ref_oct"], BLUE),
     ("gradient", R["alignment"]["n_star_int"], ORANGE)],
    "fig_contours_after.png",
    "After alignment — gradient basis re-octaved onto the value basis' G",
)

# ---- G(N) law ----------------------------------------------------------------
ns = np.arange(1, MAX_OCT + 1)
LR = R["law_fit"]["value"]["lr"]
GAIN, LAC = g["gain"], g["lacunarity"]


def g_law(n, s):
    """Continuous form of the derived G(N) law (LR = lacunarity*gain = 1)."""
    n = np.asarray(n, dtype=float)
    if abs(LR - 1.0) < 1e-12:
        geom = np.sqrt(n)
    else:
        geom = np.sqrt(np.abs((1.0 - LR ** (2 * n)) / (1.0 - LR ** 2)))
    norm = (1.0 - GAIN ** n) / (1.0 - GAIN)
    return s * geom / norm


fig, ax = plt.subplots(figsize=(8.4, 5.0))
for basis, color, marker in (("value", BLUE, "o"), ("gradient", ORANGE, "s")):
    gs = np.array([R["G_rms"][basis][str(n)] for n in ns])
    s = R["law_fit"][basis]["scale"]
    ax.plot(ns, gs, marker, color=color, ms=7,
            label=f"{basis} measured G_rms")
    xx = np.linspace(1, MAX_OCT, 200)
    ax.plot(xx, g_law(xx, s), "-", color=color, lw=1.6, alpha=0.75,
            label=f"{basis} derived law (s={s:.3f})")
ax.axhline(R["alignment"]["target_G"], color=MUTED, ls="--", lw=1.2)
ax.annotate(f"alignment target G={R['alignment']['target_G']:.3f}",
            (1.15, R["alignment"]["target_G"]),
            textcoords="offset points", xytext=(0, 6), color=MUTED, fontsize=9)
ax.axvline(R["alignment"]["n_star_int"], color=ORANGE, ls=":", lw=1.2, alpha=0.8)
ax.annotate(f"N*={R['alignment']['n_star_real']:.2f} → {R['alignment']['n_star_int']}",
            (R["alignment"]["n_star_int"], ax.get_ylim()[0]),
            textcoords="offset points", xytext=(-46, 8), color=ORANGE, fontsize=9)
ax.set_xlabel("octave count N")
ax.set_ylabel("gradient RMS  G(N)")
ax.set_xticks(ns)
ax.set_title("Gradient magnitude vs octaves: measured and derived", fontsize=12)
ax.legend(frameon=False, fontsize=9, labelcolor=INK)
fig.tight_layout()
fig.savefig(FIGS / "fig_gn_law.png", dpi=160)
plt.close(fig)

# ---- spectra -----------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8.4, 5.0))
lo_f, hi_f = R["psd_fit_range"]
for basis, color in (("value", BLUE), ("gradient", ORANGE)):
    f, P = R["psd"][basis]
    f, P = np.array(f), np.array(P)
    m = (f >= lo_f) & (f <= hi_f)
    slope, icept = np.polyfit(np.log(f[m]), np.log(P[m]), 1)
    ax.loglog(f[1:], P[1:], color=color, lw=1.3, label=f"{basis} radial PSD")
    ax.loglog(f[m], np.exp(icept) * f[m] ** slope, "--", color=color, lw=1.4,
              label=f"{basis} fit β={-slope:.2f}")
ax.axvspan(lo_f, hi_f, color=GRID, alpha=0.5, zorder=0)
ax.annotate("fit window", ((lo_f + hi_f) / 2, ax.get_ylim()[1] * 0.6),
            ha="center", color=MUTED, fontsize=9)
ax.set_xlabel("spatial frequency f (cycles / world unit)")
ax.set_ylabel("radial power P(f)")
ax.set_title(f"Spectral slopes — derived β(2D) = 2α+2 = {R['beta_derived_2d']:.1f}",
             fontsize=12)
ax.legend(frameon=False, fontsize=9, labelcolor=INK)
fig.tight_layout()
fig.savefig(FIGS / "fig_spectra.png", dpi=160)
plt.close(fig)
print("figures written to", FIGS)
