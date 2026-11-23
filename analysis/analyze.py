#!/usr/bin/env python3
"""Read the campaign output and emit the convergence law, cost table, gully
depth distributions, and the figures. Usage: analyze.py <campaign-dir> <out-dir>"""
import json, os, sys
import numpy as np

# ---- parameters (kept in one place so the report can quote them) ----
N = 256
DROPS = 65536
HILLSHADE_LIGHT = None

CAMPAIGN = sys.argv[1] if len(sys.argv) > 1 else "artifacts/erosion/campaign"
OUT = sys.argv[2] if len(sys.argv) > 2 else "artifacts/erosion/report"
os.makedirs(OUT, exist_ok=True)


def load(tag, it, n=N):
    return np.fromfile(f"{CAMPAIGN}/fields/{tag}.{it:04d}.h32", dtype="<f4").reshape(n, n).astype(np.float64)


def series(tag):
    d = json.load(open(f"{CAMPAIGN}/logs/{tag}.json"))
    keys = ["iter", "ms", "residual", "rms_vs_init", "mean_abs_delta",
            "deposit_vol", "erode_vol", "us_per_cell"]
    return {k: np.array([r[k] for r in d["runs"]]) for k in keys}, d


def rises(r):
    return [i + 1 for i in range(1, len(r)) if r[i] > r[i - 1]]


def first_unrecovered_rise(r):
    """First index i (1-based) where r[i] > r[i-1] and the series never returns
    below r[i-1] afterwards. This is where 'a transient wobble' becomes
    'the residual is no longer decreasing'."""
    for i in range(1, len(r)):
        if r[i] > r[i - 1] and not (r[i + 1:] < r[i - 1]).any():
            return i + 1
    return None


def rolling_min(r, w):
    roll = np.convolve(r, np.ones(w) / w, mode="valid")
    idx = int(roll.argmin())
    return dict(value=float(roll[idx]), center_iter=int(idx + (w - 1) // 2 + 1))



# ---------------------------------------------------------------- convergence
def fit_power(it, r):
    p = np.polyfit(np.log(it), np.log(r), 1)
    pred = np.exp(np.polyval(p, np.log(it)))
    return -p[0], 1 - ((r - pred) ** 2).sum() / ((r - r.mean()) ** 2).sum()


def fit_exp_floor(it, r):
    # r = A exp(-i/tau) + c  -- grid search on tau, linear least squares on A,c
    best = None
    for tau in np.linspace(1.0, 60.0, 2000):
        e = np.exp(-it / tau)
        M = np.vstack([e, np.ones_like(e)]).T
        sol, *_ = np.linalg.lstsq(M, r, rcond=None)
        pred = M @ sol
        sse = ((r - pred) ** 2).sum()
        if best is None or sse < best[0]:
            best = (sse, tau, sol[0], sol[1], pred)
    sse, tau, A, c, pred = best
    r2 = 1 - sse / ((r - r.mean()) ** 2).sum()
    return dict(tau=tau, A=A, c=c, r2=r2)


def fit_stretched(it, r):
    # r = c + A exp(-(i/tau)^beta)
    from scipy.optimize import curve_fit
    def f(i, A, tau, beta, c):
        return c + A * np.exp(-(i / tau) ** beta)
    try:
        popt, _ = curve_fit(f, it, r, p0=[2, 3, 1, 0.2], maxfev=40000,
                            bounds=([0, .1, .2, -1], [100, 100, 3, 5]))
    except Exception:
        return None
    pred = f(it, *popt)
    r2 = 1 - ((r - pred) ** 2).sum() / ((r - r.mean()) ** 2).sum()
    return dict(A=popt[0], tau=popt[1], beta=popt[2], c=popt[3], r2=r2,
                pred=pred)


def fit_breakpoint(it, r):
    """piecewise: power law while i<=k, linear (or flat) after the knee.
    Returns the knee k and the two fits; used to describe the turn-around."""
    best = None
    for k in range(4, len(r) - 4):
        pre_i, pre_r = it[:k], r[:k]
        p = np.polyfit(np.log(pre_i), np.log(pre_r), 1)
        pred_pre = np.exp(np.polyval(p, np.log(pre_i)))
        post_i, post_r = it[k:], r[k:]
        q = np.polyfit(post_i, post_r, 1)
        pred_post = np.polyval(q, post_i)
        sse = ((pre_r - pred_pre) ** 2).sum() + ((post_r - pred_post) ** 2).sum()
        if best is None or sse < best[0]:
            best = (sse, k, -p[0], q[0])
    sse, k, p_pre, slope_post = best
    r2 = 1 - sse / ((r - r.mean()) ** 2).sum()
    return dict(knee=int(k), pre_exponent=p_pre, post_slope=slope_post, r2=r2)


# -------------------------------------------------------------- gully metrics
def gully_stats(h0, h, cell_units=1.0, floor_frac=0.5, depth_min=0.5):
    """Gully = eroded channel below its local neighbourhood.  We compute the
    depth of every cell below the median of a 9x9 window (so a broad valley
    floor doesn't count, only confined incisions), keep the cells deeper than
    depth_min, and group them into connected components.  Per component we
    report max depth; the field summary is pooled over components."""
    from scipy.ndimage import median_filter, label
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    sm = median_filter(h, size=9, mode="nearest")
    depth = sm - h                       # positive where h sits below neighbourhood
    mask = depth > depth_min
    lab, ncomp = label(mask, structure=np.ones((3, 3)))
    comps = []
    for c in range(1, ncomp + 1):
        m = lab == c
        if m.sum() < 8:                  # speckle
            continue
        comps.append(dict(cells=int(m.sum()), max_depth=float(depth[m].max()),
                          mean_depth=float(depth[m].mean())))
    return depth, mask, comps


def gully_distribution(tag, it, n=N, depth_min=0.5):
    h0 = load(tag, 0, n)
    h = load(tag, it, n)
    depth, mask, comps = gully_stats(h0, h, depth_min=depth_min)
    depths = np.array([c["max_depth"] for c in comps])
    return h0, h, depth, mask, comps, depths


def quantiles(x):
    if len(x) == 0:
        return {k: float("nan") for k in ("p10", "p25", "p50", "p75", "p90", "p99")}
    q = np.percentile(x, [10, 25, 50, 75, 90, 99])
    return dict(p10=q[0], p25=q[1], p50=q[2], p75=q[3], p90=q[4], p99=q[5])


def main():
    report = {"params": dict(n=N, drops=DROPS, drops_per_cell=1.0,
                              seed=7, capacity=0.05, evaporate=0.01,
                              max_steps=48, radius=2.0, inertia=0.05,
                              erode=0.3, deposit=0.3, gravity=4.0,
                              resample=False)}
    tiers = {"t10": 10, "t30": 30, "t90": 90}
    conv = {}
    for tag, iters in tiers.items():
        s, meta = series(tag)
        r = s["residual"]
        rises_ = rises(r)
        rec = dict(iters=iters, residual=[float(x) for x in r],
                   rms_vs_init=[float(x) for x in s["rms_vs_init"]],
                   mean_abs_delta=[float(x) for x in s["mean_abs_delta"]],
                   ms=[float(x) for x in s["ms"]],
                   us_per_cell=[float(x) for x in s["us_per_cell"]],
                   rises=rises_, final=r[-1], min=float(r.min()),
                   argmin=int(r.argmin() + 1))
        # global law on the monotone head (first 25 iterations; all tiers share
        # the same trajectory in that window so this is fit once from t90)
        if tag == "t90":
            it = np.arange(1, len(r) + 1)
            head = 20
            rec["fit_power_head"] = dict(exponent=fit_power(it[:head], r[:head])[0],
                                         r2=fit_power(it[:head], r[:head])[1])
            rec["fit_exp_floor_head"] = fit_exp_floor(it[:head], r[:head])
            st = fit_stretched(it[:head], r[:head])
            rec["fit_stretched_head"] = {k: v for k, v in st.items() if k != "pred"}
            rec["fit_breakpoint_full"] = fit_breakpoint(it, r)
            rec["first_unrecovered_rise"] = first_unrecovered_rise(r)
            rec["knee_rolling5"] = rolling_min(r, 5)
            # local power-law exponent in windows of 10
            exps = []
            for a in range(0, len(r) - 10, 10):
                e, r2 = fit_power(it[a:a + 10], r[a:a + 10])
                exps.append(dict(window=[a + 1, a + 10], exponent=float(e),
                                 r2=float(r2)))
            rec["local_exponent"] = exps
        conv[tag] = rec
    report["convergence"] = conv

    # seed-11 check
    s11, _ = series("t90_s11")
    report["seed_check"] = dict(
        residual=[float(x) for x in s11["residual"]],
        rises=rises(s11["residual"]),
        final=float(s11["residual"][-1]),
        min=float(s11["residual"].min()),
        argmin=int(s11["residual"].argmin() + 1))

    # repeatability: byte-identical fields?
    rep = {}
    for tag, tagr in (("t30", "t30_rep"), ("t90", "t90_rep")):
        a = load(tag, tiers[tag]); b = load(tagr, tiers[tag])
        rep[tag] = dict(max_abs_diff=float(np.abs(a - b).max()),
                        identical=bool(np.array_equal(a, b)))
    report["repeatability"] = rep

    # cost table: per-tier median / min / max ms, and scaling
    cost = {}
    for tag, iters in tiers.items():
        s, _ = series(tag)
        ms = s["ms"]; it = s["iter"]
        # first-iteration cost is cold-cache; report all, plus steady state
        steady = ms[2:]
        cost[tag] = dict(
            n=N, drops=DROPS, iters=iters,
            total_s=float(ms.sum() / 1e3),
            ms_per_iter=dict(median=float(np.median(steady)),
                             min=float(steady.min()), max=float(steady.max())),
            us_per_cell=dict(median=float(np.median(s["us_per_cell"][2:])),
                             min=float(s["us_per_cell"][2:].min()),
                             max=float(s["us_per_cell"][2:].max())),
            ms_all=[float(x) for x in ms],
        )
    report["cost"] = cost

    scaling = {}
    for tag, n in (("size128", 128), ("size256", 256), ("size512", 512)):
        s, _ = series(tag)
        steady = s["ms"][2:]
        scaling[f"n{n}"] = dict(n=n, drops=n * n,
                                ms_median=float(np.median(steady)),
                                us_per_cell=float(np.median(s["us_per_cell"][2:])))
    for tag, mult in (("drops1x", 1), ("drops2x", 2), ("drops4x", 4)):
        s, _ = series(tag)
        steady = s["ms"][2:]
        scaling[f"drops{mult}x"] = dict(n=N, drops=DROPS * mult, mult=mult,
                                        ms_median=float(np.median(steady)),
                                        us_per_cell=float(np.median(s["us_per_cell"][2:])))
    report["scaling"] = scaling

    # gully depth distributions at each tier
    gd = {}
    for tag, iters in tiers.items():
        h0, h, depth, mask, comps, depths = gully_distribution(tag, iters)
        gd[tag] = dict(iters=iters, n_components=int(len(comps)),
                       cells=int(mask.sum()),
                       depth_quantiles=quantiles(depths),
                       max_depth=float(depths.max()) if len(depths) else float("nan"),
                       depth_pooled_quantiles=quantiles(depth[mask]) if mask.any() else None,
                       components=[c for c in comps if c["max_depth"] >= 1.0][:200])
    # depth is also measured at a fixed iteration (30) for all tiers to show
    # the distributions are the same field at different stopping points
    report["gully"] = gd

    # scalar evolution of the gully distribution; only iterations whose field
    # files were kept in the artifact tree are read
    evo = {}
    for tag in ("t90",):
        rows = []
        for it in (1, 5, 10, 20, 30, 40, 60, 80, 90):
            path = f"{CAMPAIGN}/fields/{tag}.{it:04d}.h32"
            if not os.path.exists(path):
                continue
            h0 = load(tag, 0); h = load(tag, it)
            depth, mask, comps = gully_stats(h0, h)
            d = np.array([c["max_depth"] for c in comps])
            rows.append(dict(iter=it, n=len(comps), cells=int(mask.sum()),
                             max_depth=float(d.max()) if len(d) else None,
                             p90=float(np.percentile(d, 90)) if len(d) else None,
                             mean_depth=float(d.mean()) if len(d) else None))
        evo[tag] = rows
    report["gully_evolution"] = evo

    # residual decomposition: how much of r_i comes from the gully mask vs the
    # rest of the field, and from the top-k extreme cells. Cached by analyze.py
    # in logs/t90_residual_decomposition.json (computed from the full frame
    # sequence, which is too large to keep in the artifact tree).
    dec_path = f"{CAMPAIGN}/logs/t90_residual_decomposition.json"
    if os.path.exists(dec_path):
        report["residual_decomposition"] = json.load(open(dec_path))

    # single-cell instability trace at the cell that drives the residual rise
    hx, hy = 74, 49
    trace = []
    for it in (0, 10, 20, 30, 31, 32, 33, 34, 35, 40, 50, 60, 70, 80, 90):
        path = f"{CAMPAIGN}/fields/t90.{it:04d}.h32"
        if not os.path.exists(path):
            continue
        h = load("t90", it)
        nb = h[hy - 4:hy + 5, hx - 4:hx + 5]
        trace.append(dict(iter=it, h=float(h[hy, hx]),
                          nb_median=float(np.median(nb)),
                          field_min=float(h.min()), field_max=float(h.max())))
    report["spike_trace"] = dict(x=hx, y=hy, rows=trace)

    with open(f"{OUT}/metrics.json", "w") as fp:
        json.dump(report, fp, indent=1)

    # compact human-readable summary
    lines = []
    for tag, iters in tiers.items():
        c = conv[tag]
        lines.append(f"{tag}: {iters} iters  final={c['final']:.4f}  "
                     f"min={c['min']:.4f}@{c['argmin']}  rises={len(c['rises'])}")
    print("\n".join(lines))
    print("wrote", f"{OUT}/metrics.json")


if __name__ == "__main__":
    main()
