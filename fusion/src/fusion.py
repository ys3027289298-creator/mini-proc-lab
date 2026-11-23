#!/usr/bin/env python3
"""fusion: smooth fusion of heat gradient and elevation gradient.

Problem: the naive composite gates the elevation-gradient term with a hard
mask 1[s > theta], so the fused heat field F0 is discontinuous along the
isoline s = theta (and across the cliff in h itself) -> hard edges.

Smooth fusion (see fusion/README.md for the full derivation):
    h~ = G_sigma * h                      (Gaussian-smoothed elevation, C^inf)
    s  = |grad h~|                        (continuous slope field)
    w  = 1 / (1 + exp(-(s - c) / (k*r)))  (logistic weight, C^inf, robust
         center c = mean(s)+0.5*std(s), scale r = std(s))
    F  = T(h~) + alpha * w(s) * s         (fused field)
    grad F = T'(h~) grad h~ + alpha (w + s w') grad s
All factors are continuous, so F is C^1: no step pixels by construction,
provided sigma is large enough that the smoothed cliff ramp fits under the
step threshold tau (bound derived in README).

A pixel is a STEP PIXEL if any 4-neighbour difference exceeds
tau = 0.1 * (P99 - P1) of the naive fused field.

Renders one heatmap per fusion tier plus a comparison figure; exits 1 if
any tier still has step pixels.
"""
import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter

N = 256          # grid size
A_CLIFF = 5.0    # cliff height (the elevation jump)
LAPSE = 1.2      # heat lapse rate dT/dh
ALPHA = 2.0      # fusion gain on the elevation-gradient term
GAMMA = 0.5      # gate center: mean(s) + GAMMA*std(s)

# three fusion tiers: (name, sigma, k)
TIERS = [
    ("tier1_light", 4.0, 0.10),
    ("tier2_mid", 8.0, 0.06),
    ("tier3_heavy", 16.0, 0.03),
]

OUT = Path(__file__).resolve().parent.parent / "out"


def fbm(n, octaves, seed, amp=1.0):
    """Deterministic bilinear value-noise fBm, matching the lab's noise module."""
    rng = np.random.default_rng(seed)
    out = np.zeros((n, n))
    total = 0.0
    a = 1.0
    t = None
    for o in range(octaves):
        g = 4 << o
        coarse = rng.random((g + 1, g + 1))
        t = np.linspace(0.0, g, n)
        i0 = np.minimum(t.astype(int), g - 1)
        f = t - i0
        i1 = i0 + 1
        fx = f[None, :]
        fy = f[:, None]
        c00 = coarse[np.ix_(i0, i0)]
        c01 = coarse[np.ix_(i0, i1)]
        c10 = coarse[np.ix_(i1, i0)]
        c11 = coarse[np.ix_(i1, i1)]
        layer = (c00 * (1 - fx) + c01 * fx) * (1 - fy) + (c10 * (1 - fx) + c11 * fx) * fy
        out += a * layer
        total += a
        a *= 0.5
    return amp * out / total


def make_terrain(n=N):
    """Rolling base terrain plus a hard cliff (true elevation jump)."""
    h = 1.2 * fbm(n, octaves=4, seed=11)
    y = np.arange(n)
    x0 = n // 2 + (18.0 * np.sin(0.06 * y) + 8.0 * np.sin(0.017 * y + 1.0)).astype(int)
    cliff = (np.arange(n)[None, :] >= x0[:, None]).astype(float)
    return h + A_CLIFF * cliff


def slope(h):
    gy, gx = np.gradient(h)
    return np.hypot(gx, gy)


def logistic_gate(s, k):
    """w = 1/(1+exp(-(s-c)/(k*r))) with robust center c and scale r."""
    c = s.mean() + GAMMA * s.std()
    r = s.std()
    return 1.0 / (1.0 + np.exp(-(s - c) / (k * r)))


def heat(h_field, eta):
    return 26.0 - LAPSE * h_field + eta


def step_mask(F, tau):
    """Pixels whose 4-neighbour max difference exceeds tau."""
    dx = np.abs(np.diff(F, axis=1))
    dy = np.abs(np.diff(F, axis=0))
    m = np.zeros(F.shape, bool)
    m[:, :-1] |= dx > tau
    m[:, 1:] |= dx > tau
    m[:-1, :] |= dy > tau
    m[1:, :] |= dy > tau
    return m


def max_step(F):
    dx = np.abs(np.diff(F, axis=1)).max()
    dy = np.abs(np.diff(F, axis=0)).max()
    return max(dx, dy)


def render(F, mask, title, path):
    fig, ax = plt.subplots(figsize=(6, 5.4), dpi=110)
    im = ax.imshow(F, cmap="inferno", origin="upper")
    if mask.any():
        ys, xs = np.nonzero(mask)
        ax.scatter(xs, ys, s=1, c="cyan", label=f"step px ({mask.sum()})")
        ax.legend(loc="lower right", fontsize=8)
    ax.set_title(title, fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])
    fig.colorbar(im, ax=ax, fraction=0.046)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    h = make_terrain()
    eta = 0.6 * fbm(N, octaves=3, seed=23)

    # --- naive fusion: hard gate on the raw (unsmoothed) slope ---
    s0 = slope(h)
    gate0 = s0 > (s0.mean() + GAMMA * s0.std())
    F0 = heat(h, eta) + ALPHA * gate0 * s0

    tau = 0.1 * (np.percentile(F0, 99) - np.percentile(F0, 1))
    m0 = step_mask(F0, tau)
    print(f"tau={tau:.4f}  (0.1 * (P99-P1) of naive field)")
    print(f"{'field':<14}{'sigma':>6}{'k':>6}{'step_px':>9}{'max|dF|':>10}")
    print(f"{'naive':<14}{'-':>6}{'-':>6}{m0.sum():>9}{max_step(F0):>10.3f}")
    render(F0, m0, f"naive hard gate — step px: {m0.sum()}", OUT / "fusion_naive.png")

    rows = [("naive", "-", "-", int(m0.sum()), f"{max_step(F0):.4f}")]
    fields = [("naive", F0, m0)]
    failures = 0

    for name, sigma, k in TIERS:
        hs = gaussian_filter(h, sigma, mode="nearest")
        s = slope(hs)
        w = logistic_gate(s, k)
        eta_s = gaussian_filter(eta, sigma / 2.0, mode="nearest")
        F = heat(hs, eta_s) + ALPHA * w * s
        m = step_mask(F, tau)
        ms = max_step(F)
        print(f"{name:<14}{sigma:>6.1f}{k:>6.2f}{m.sum():>9}{ms:>10.3f}")
        rows.append((name, f"{sigma:.1f}", f"{k:.2f}", int(m.sum()), f"{ms:.4f}"))
        render(F, m, f"{name}  sigma={sigma} k={k} — step px: {m.sum()}",
               OUT / f"fusion_{name}.png")
        fields.append((name, F, m))
        if m.sum() != 0:
            failures += 1

    # --- comparison figure: naive + three tiers ---
    fig, axes = plt.subplots(2, 2, figsize=(11, 10), dpi=110)
    for ax, (name, F, m) in zip(axes.flat, fields):
        im = ax.imshow(F, cmap="inferno", origin="upper")
        if m.any():
            ys, xs = np.nonzero(m)
            ax.scatter(xs, ys, s=1, c="cyan")
        ax.set_title(f"{name} — step px: {m.sum()}", fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(im, ax=ax, fraction=0.046)
    fig.suptitle("heat/elevation-gradient fusion: naive vs smoothed tiers", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT / "fusion_compare.png")
    plt.close(fig)

    with open(OUT / "stats.csv", "w", newline="") as fh:
        wcsv = csv.writer(fh)
        wcsv.writerow(["field", "sigma", "k", "step_pixels", "max_neighbor_diff"])
        wcsv.writerows(rows)

    print(f"all_ok={1 if failures == 0 else 0}")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
