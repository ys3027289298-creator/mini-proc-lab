#!/usr/bin/env python3
"""Independent acceptance check over the tunnel CSVs.

Re-derives every quantity from the raw flow-line / station outputs rather
than trusting the pass flags written by the C binary. Exits non-zero if any
segment of any terrain violates its budget or the 8 degree heading limit.
"""
import csv
import math
import os
import sys

ALPHA_DEG = 8.0
E_CHORD = 0.20
CA = 1 - math.cos(math.radians(ALPHA_DEG))
NAMES = ["rolling_hills", "single_ridge", "steep_bent_monocline"]


def rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def main(data):
    all_ok = True
    for name in NAMES:
        stream = rows(os.path.join(data, f"stream_{name}.csv"))
        seg = rows(os.path.join(data, f"segstats_{name}.csv"))
        new_stations = rows(os.path.join(data, f"new_{name}.csv"))

        # 1. every new-scheme station must sit on the ideal flow line
        pts = [(float(r["x"]), float(r["y"])) for r in stream]
        max_anchor = 0.0
        for st in new_stations:
            sx, sy = float(st["x"]), float(st["y"])
            dmin = min(math.hypot(sx - px, sy - py) for px, py in pts)
            max_anchor = max(max_anchor, dmin)
        anchor_ok = max_anchor < 5e-3

        # 2. recompute the budget expression independently per segment
        bad = []
        for r in seg:
            if r["scheme"] != "new":
                continue
            kmax = float(r["kappa_max"])
            budget = E_CHORD if kmax <= 1e-9 else min(E_CHORD, CA / kmax)
            # CSV stores kappa/budget at finite precision; 1e-4 covers the
            # rounding, which is 3 orders below the smallest real budget (~0.04)
            if abs(budget - float(r["d_allow"])) > 1e-4:
                bad.append(f"seg {r['seg']} budget mismatch")
            if float(r["d_max"]) > budget + 1e-6:
                bad.append(f"seg {r['seg']} d_max {r['d_max']} > {budget:.4f}")
            if float(r["theta_deg"]) > ALPHA_DEG + 1e-3:
                bad.append(f"seg {r['seg']} heading {r['theta_deg']} > 8")
            # cross-check against the C-reported flag
            if r["pass"] != "1":
                bad.append(f"seg {r['seg']} C flag not pass")

        old = [r for r in seg if r["scheme"] == "old"]
        old_fail = sum(float(r["d_max"]) > float(r["d_allow"]) + 1e-6
                       or float(r["theta_deg"]) > ALPHA_DEG + 1e-3
                       for r in old)

        ok = anchor_ok and not bad
        all_ok &= ok
        print(f"{name:22s} anchors_on_streamline={anchor_ok} "
              f"(max gap {max_anchor:.2e} m)  new_segments="
              f"{sum(r['scheme']=='new' for r in seg):2d}  "
              f"violations={len(bad)}  old_scheme_failed_segs={old_fail}/{len(old)}"
              f"  -> {'PASS' if ok else 'FAIL ' + '; '.join(bad)}")

    print("INDEPENDENT ACCEPTANCE:", "PASS" if all_ok else "FAIL")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "artifacts/data"))
