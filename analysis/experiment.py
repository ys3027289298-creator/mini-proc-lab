"""Contour-width / gradient-magnitude / octave alignment experiment.

Pipeline
  1. render the same mountain window with two noise bases (value, gradient)
     using the shared amplitude-normalized fBm;
  2. measure the width of three iso-level ribbons rendered with a FIXED level
     band delta (marching-squares sample points + Euclidean distance transform,
     curvature-free), so width traces to the gradient w = 2*delta/|grad F|;
  3. tabulate gradient RMS vs octave count, fit the derived G(N) law and solve
     the octave count that aligns basis B to basis A;
  4. measure each basis' spectral slope two ways (per-octave gradient RMS,
     radial FFT power law) and place next to the derived value.

Outputs: results.json, widths.csv (figures are made by make_figures.py).
"""
import csv
import json
import os
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import ndimage

from noise_bases import fbm_partials, _BASES

HERE = Path(__file__).resolve().parent
OUT = Path(os.environ.get("EXP_OUT", HERE))
FIGS = OUT / "figs"
FIGS.mkdir(parents=True, exist_ok=True)

# ---- fixed experiment parameters -------------------------------------------
N = int(os.environ.get("EXP_GRID", 2048))   # grid points per side
DOMAIN = 16.0            # world units per side (same mountain window for both bases)
CELL = DOMAIN / N        # world units per pixel
LEVELS = (-0.25, 0.0, 0.25)
LAC, GAIN, SEED = 2.0, 0.5, 11
MAX_OCT = 6
FIT_OCT = 5            # octaves 1..FIT_OCT are alias-free on the N=2048 grid
REF_OCT = 6            # reference basis (value) octave count
MISMATCH_OCT = 2       # second basis (gradient) octave count before alignment
TARGET_MIN_PX = 4.0      # thinnest ribbon allowed when fixing the level band delta

rng = np.linspace(-DOMAIN / 2, DOMAIN / 2, N, dtype=np.float32)
X, Y = np.meshgrid(rng, rng, indexing="ij")


def grad_rms(F):
    gx = np.gradient(F.astype(np.float64), CELL, axis=0)
    gy = np.gradient(F.astype(np.float64), CELL, axis=1)
    g = np.hypot(gx, gy)
    m = 8                                   # interior margin (px)
    return float(np.sqrt(np.mean(g[m:-m, m:-m] ** 2))), g


def raw_octave(x, y, i, basis):
    """One un-normalized fBm octave: gain^i * basis(lacunarity^i * p)."""
    freq = np.float32(LAC ** i)
    amp = np.float32(GAIN ** i)
    return amp * _BASES[basis](x * freq, y * freq, SEED + i * 17)


# ---- 1+2: fields, gradients -------------------------------------------------
fields = {}              # (basis, n) -> F
grads = {}               # (basis, n) -> |grad|
G = defaultdict(dict)    # basis -> {n: gradient RMS}
for basis in ("value", "gradient"):
    partials = fbm_partials(X, Y, MAX_OCT, LAC, GAIN, SEED, basis)
    for n, F in enumerate(partials, start=1):
        g_rms, g = grad_rms(F)
        fields[(basis, n)] = F
        grads[(basis, n)] = g
        G[basis][n] = g_rms

# Fixed level band: thinnest typical ribbon >= TARGET_MIN_PX, capped below half
# the level spacing so adjacent iso-ribbons never merge.
scenario_keys = [("value", REF_OCT), ("gradient", MISMATCH_OCT)]
g99 = max(float(np.quantile(grads[k], 0.99)) for k in scenario_keys)
DELTA = min(TARGET_MIN_PX * CELL * g99, 0.09)
print(f"cell={CELL:.6f} world/px   g99={g99:.4f}   fixed band delta={DELTA:.5f}")

# ---- 3: marching-squares sample points + gradient-normal ribbon widths -------
def contour_points(F, level):
    """Midpoints of marching-squares segments: one sample set per iso-contour.

    Topology-free on purpose -- ribbon widths only need well-distributed points
    sitting on the level set, not chained polylines.
    """
    V = F - np.float32(level)
    inside = V > 0
    crossings = {}
    ih, jh = np.nonzero(inside[:-1, :] != inside[1:, :])
    for a, b in zip(ih.tolist(), jh.tolist()):
        t = float(V[a, b] / (V[a, b] - V[a + 1, b]))
        crossings[("h", a, b)] = np.array([a + t, float(b)])
    iv, jv = np.nonzero(inside[:, :-1] != inside[:, 1:])
    for a, b in zip(iv.tolist(), jv.tolist()):
        t = float(V[a, b] / (V[a, b] - V[a, b + 1]))
        crossings[("v", a, b)] = np.array([float(a), b + t])

    code = (inside[:-1, :-1].astype(np.uint8)
            | (inside[1:, :-1].astype(np.uint8) << 1)
            | (inside[1:, 1:].astype(np.uint8) << 2)
            | (inside[:-1, 1:].astype(np.uint8) << 3))
    ci, cj = np.nonzero(code != 0)
    ck = code[ci, cj]
    pts = []
    for a, b, c in zip(ci.tolist(), cj.tolist(), ck.tolist()):
        cell_edges = [("h", a, b), ("v", a, b), ("h", a + 1, b), ("v", a, b + 1)]
        edge_pts = [crossings[e] for e in cell_edges if e in crossings]
        if len(edge_pts) == 2:
            pts.append((edge_pts[0] + edge_pts[1]) * 0.5)
        elif len(edge_pts) == 4:  # saddle: two segments, split at centroid
            ctr = sum(edge_pts) / 4.0
            order = sorted(edge_pts,
                           key=lambda p: np.arctan2(p[1] - ctr[1], p[0] - ctr[0]))
            pts.append((order[0] + ctr) * 0.5)
            pts.append((order[2] + ctr) * 0.5)
    return np.asarray(pts) if pts else np.empty((0, 2))


def _bilinear(grid, p):
    x, y = p
    x0, y0 = int(x), int(y)
    tx, ty = x - x0, y - y0
    return (grid[x0, y0] * (1 - tx) * (1 - ty)
            + grid[x0 + 1, y0] * tx * (1 - ty)
            + grid[x0, y0 + 1] * (1 - tx) * ty
            + grid[x0 + 1, y0 + 1] * tx * ty)


def ribbon_widths(F, level, delta):
    """Ribbon width at marching-squares contour points via Euclidean distance
    transform of the |F-level|<=delta mask: w(p) = 2 d_EDT(p). Geometry-based,
    so curvature cannot bias the result. Local |grad F| sampled at the same
    points for the w = 2 delta / |grad F| identity check."""
    F64 = F.astype(np.float64)
    gx = np.gradient(F64, CELL, axis=0)
    gy = np.gradient(F64, CELL, axis=1)
    gnorm = np.hypot(gx, gy)
    pts = contour_points(F, level)
    mask = np.abs(F64 - level) <= delta
    # distance (in pixels) from every inside pixel to the band boundary
    dist_px = ndimage.distance_transform_edt(mask)
    widths, local_g = [], []
    margin = 40.0
    for k in range(0, len(pts), 2):
        p = pts[k]
        if not (margin < p[0] < N - margin and margin < p[1] < N - margin):
            continue
        x0, y0 = int(p[0]), int(p[1])
        tx, ty = p[0] - x0, p[1] - y0
        d = (dist_px[x0, y0] * (1 - tx) * (1 - ty)
             + dist_px[x0 + 1, y0] * tx * (1 - ty)
             + dist_px[x0, y0 + 1] * (1 - tx) * ty
             + dist_px[x0 + 1, y0 + 1] * tx * ty)
        if d <= 0:
            continue
        widths.append(2.0 * d)                       # full ribbon width in px
        local_g.append(_bilinear(gnorm, p))
    return np.asarray(widths), np.asarray(local_g)


def edt_width(F, level, delta):
    """Independent cross-check of ribbon_widths: same contour sample points and
    the same 2*d_EDT diameter, via a separate nearest-neighbor code path."""
    mask = np.abs(F - np.float32(level)) <= np.float32(delta)
    dist_px = ndimage.distance_transform_edt(mask)
    pts = contour_points(F, level)
    vals = []
    margin = 40.0
    for k in range(0, len(pts), 2):
        p = pts[k]
        if margin < p[0] < N - margin and margin < p[1] < N - margin:
            d = dist_px[int(round(p[0])), int(round(p[1]))]
            if d > 0:
                vals.append(2.0 * d)
    return float(np.median(vals)) if vals else float("nan")


width_rows = []
edt_check = {}
for basis, n in (("value", REF_OCT), ("gradient", MISMATCH_OCT)):
    key = (basis, n)
    for lv in LEVELS:
        w, lg = ribbon_widths(fields[key], lv, DELTA)
        if (basis, lv) == ("value", 0.0):
            edt_check[basis] = edt_width(fields[key], lv, DELTA)
        width_rows.append(dict(scenario="before", basis=basis, octaves=n, level=lv,
                               delta=DELTA, width_px=float(np.median(w)),
                               width_iqr=float(np.subtract(*np.percentile(w, [75, 25]))),
                               width_world=float(np.median(w) * CELL),
                               grad_on_contour=float(np.median(lg)),
                               n_samples=len(w)))

# ---- 4: G(N) law fit and octave alignment -----------------------------------
# The C fBm is amplitude-normalized (divided by sum_i r^i). With independent
# octaves, per-octave gradient RMS q_i = q0 (Lr)^i, hence
#   G(N) = s * sqrt(sum_i (Lr)^{2i}) / sum_i r^i
#        = s * sqrt((1-(Lr)^{2N})/(1-(Lr)^2)) * (1-r)/(1-r^N)
# where s absorbs the lattice gradient RMS and all constants.
LR = LAC * GAIN


def g_law(n, s):
    """G(N) = s * sqrt(sum_{i<N} (Lr)^{2i}) / sum_{i<N} r^i (finite sums,
    valid also at Lr = 1)."""
    scalar = np.isscalar(n)
    n = np.atleast_1d(np.asarray(n, dtype=int))
    out = np.empty(n.shape, dtype=float)
    for k in np.ndindex(n.shape):
        idx = np.arange(n[k])
        geom = np.sqrt(np.sum(LR ** (2 * idx)))
        norm = np.sum(GAIN ** idx)
        out[k] = s * geom / norm
    return out[0] if scalar else out


def g_law_cont(n, s):
    """Continuous-N version of g_law for numerical inversion."""
    n = np.asarray(n, dtype=float)
    if abs(LR - 1.0) < 1e-12:
        geom = np.sqrt(n)
    else:
        geom = np.sqrt(np.abs((1.0 - LR ** (2 * n)) / (1.0 - LR ** 2)))
    norm = (1.0 - GAIN ** n) / (1.0 - GAIN)
    return s * geom / norm


fit = {}
for basis in ("value", "gradient"):
    ns = np.arange(1, FIT_OCT + 1)
    gs = np.array([G[basis][n] for n in ns])
    s = float(np.median(gs / g_law(ns, 1.0)))
    fit[basis] = s
    pred = g_law(ns, s)
    print(f"{basis:9s} scale s={s:.5f}  max rel err of G(N) law (N<= {FIT_OCT}): "
          f"{np.max(np.abs(pred - gs) / gs):.4f}")

# Align gradient basis (B) to value basis (A) at REF_OCT.
sA, sB = fit["value"], fit["gradient"]
g_target = g_law(REF_OCT, sA)
# invert numerically (monotone in N over the relevant range)
grid_n = np.linspace(1, MAX_OCT + 2, 4001)
n_star = float(np.interp(g_target, g_law_cont(grid_n, sB), grid_n))
n_align = int(round(n_star))
print(f"alignment target G={g_target:.4f}  N*_B(real)={n_star:.3f} -> int {n_align}")

# render the aligned scenario and measure its three ribbons
for lv in LEVELS:
    for basis, n in (("value", REF_OCT), ("gradient", n_align)):
        key = (basis, n)
        w, lg = ribbon_widths(fields[key], lv, DELTA)
        width_rows.append(dict(scenario="after", basis=basis, octaves=n, level=lv,
                               delta=DELTA, width_px=float(np.median(w)),
                               width_iqr=float(np.subtract(*np.percentile(w, [75, 25]))),
                               width_world=float(np.median(w) * CELL),
                               grad_on_contour=float(np.median(lg)),
                               n_samples=len(w)))

# ---- 5: spectral slopes ------------------------------------------------------
def radial_psd(F):
    Fw = F.astype(np.float64)
    Fw = Fw - Fw.mean()
    win = np.hanning(N)[:, None] * np.hanning(N)[None, :]
    spec = np.fft.fftshift(np.fft.fft2(Fw * win))
    P = (np.abs(spec) ** 2)
    k = np.arange(N) - N // 2
    kx, ky = np.meshgrid(k, k, indexing="ij")
    rad = np.sqrt(kx ** 2 + ky ** 2).astype(int)
    bin_idx = np.arange(1, N // 2)
    Pavg = ndimage.mean(P, labels=rad, index=bin_idx)
    f_cyc = bin_idx / DOMAIN                    # cycles per world unit
    return f_cyc, np.asarray(Pavg)


psd = {}
beta_fit = {}
# interior fBm bands only: above the fundamental lattice band, below the top
# (window/alias-contaminated) octave
lo_f, hi_f = 1.25 * LAC, 0.6 * LAC ** (MAX_OCT - 1)
for basis in ("value", "gradient"):
    f, P = radial_psd(fields[(basis, MAX_OCT)])
    psd[basis] = (f, P)
    m = (f >= lo_f) & (f <= hi_f)
    slope, intercept = np.polyfit(np.log(f[m]), np.log(P[m]), 1)
    beta_fit[basis] = float(-slope)
    print(f"{basis:9s} radial-PSD fit f in [{lo_f:.2f},{hi_f:.1f}] cyc/unit  "
          f"beta={-slope:.3f}")

# per-octave gradient RMS -> empirical amplitude exponent alpha
alpha_rows = []
for basis in ("value", "gradient"):
    q, sig = [], []
    for i in range(MAX_OCT):
        c = raw_octave(X, Y, i, basis)
        gx = np.gradient(c.astype(np.float64), CELL, axis=0)
        gy = np.gradient(c.astype(np.float64), CELL, axis=1)
        q.append(np.sqrt(np.mean(gx ** 2 + gy ** 2)))
        sig.append(np.sqrt(np.mean(c.astype(np.float64) ** 2)))
    q, sig = np.array(q), np.array(sig)
    slope_q = np.polyfit(np.arange(MAX_OCT), np.log(q), 1)[0]   # = (1-alpha) ln L
    alpha_emp = 1.0 - slope_q / np.log(LAC)
    # for an isotropic 2D field the radial power law is P(f) ~ f^{-(2 alpha + 2)}
    alpha_rows.append((basis, sig.tolist(), q.tolist(), float(alpha_emp),
                       float(2 * alpha_emp + 2)))
    print(f"{basis:9s} per-octave grad RMS={np.round(q, 3)}  alpha_emp={alpha_emp:.3f}")

# estimator calibration: synthetic isotropic fields with known power law
k = np.arange(N) - N // 2
kx, ky = np.meshgrid(k, k, indexing="ij")
kr = np.maximum(np.hypot(kx, ky), 1.0)
rng_cal = np.random.default_rng(42)
calibration = {}
for beta_true in (2.0, 3.0, 4.0):
    phase = rng_cal.uniform(0, 2 * np.pi, size=(N, N))
    amp = kr ** (-beta_true / 2.0)
    Fcal = np.real(np.fft.ifft2(np.fft.ifftshift(amp * np.exp(1j * phase))))
    f, P = radial_psd(Fcal.astype(np.float32))
    m = (f >= lo_f) & (f <= hi_f)
    beta_hat = float(-np.polyfit(np.log(f[m]), np.log(P[m]), 1)[0])
    calibration[f"beta{beta_true:.0f}"] = dict(true=beta_true, measured=beta_hat,
                                               bias=beta_hat - beta_true)
    print(f"synthetic beta={beta_true:.0f} -> measured {beta_hat:.3f} "
          f"(bias {beta_hat - beta_true:+.3f})")

# ---- save --------------------------------------------------------------------
with (OUT / "widths.csv").open("w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=list(width_rows[0].keys()))
    writer.writeheader()
    writer.writerows(width_rows)

results = dict(
    grid=dict(N=N, domain=DOMAIN, cell=CELL, lacunarity=LAC, gain=GAIN,
              levels=list(LEVELS), ref_oct=REF_OCT, mismatch_oct=MISMATCH_OCT,
              delta=DELTA),
    G_rms={b: {str(n): G[b][n] for n in range(1, MAX_OCT + 1)}
           for b in ("value", "gradient")},
    law_fit={b: dict(scale=fit[b], lr=LR) for b in fit},
    alignment=dict(target_G=float(g_target), n_star_real=float(n_star),
                   n_star_int=n_align),
    edt_crosscheck=edt_check,
    beta_psd=beta_fit,
    psd_calibration=calibration,
    psd_fit_range=[lo_f, hi_f],
    alpha_per_octave=[dict(basis=b, sigma=sig, grad_rms=q, alpha_emp=a, beta_emp=be)
                      for (b, sig, q, a, be) in alpha_rows],
    # amplitude exponent alpha = -ln(r)/ln(L); isotropic 2D radial slope beta=2alpha+2
    alpha_derived=float(-np.log(GAIN) / np.log(LAC)),
    beta_derived_2d=2.0 * (-np.log(GAIN) / np.log(LAC)) + 2.0,
    psd={b: [psd[b][0].tolist(), psd[b][1].tolist()] for b in psd},
)
with (OUT / "results.json").open("w") as fh:
    json.dump(results, fh, indent=2)
print("saved results.json, widths.csv")
