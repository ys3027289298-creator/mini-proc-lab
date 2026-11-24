#!/usr/bin/env python3
"""Generate the 3 speeds x 3 directions matrix, measure every texture,
lay out the 9-row table, and check whether the two maps stayed independent.

Fixed path: spacing must be monotone in speed within every direction and
identical across directions at the same speed. Coupled path is expected to
break that (direction leaked into the spacing map); the run reports the
offending cells and exits non-zero.

Run under: python3 -I run_matrix.py [--coupled]
"""
import os
import math
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BIN = os.path.join(ROOT, "ripple")

SPEEDS = [4.0, 7.0, 10.0]
DIRS_DEG = [0.0, 45.0, 90.0]


def main():
    coupled = "--coupled" in sys.argv[1:]
    mode = "coupled" if coupled else "fixed"
    outdir = os.path.join(ROOT, "out", mode)
    os.makedirs(outdir, exist_ok=True)

    extra = ["coupled"] if coupled else []
    rows = []
    for v in SPEEDS:
        for d in DIRS_DEG:
            path = os.path.join(outdir, "s%02d_d%03d.pgm" % (int(v), int(d)))
            subprocess.run([BIN, "gen", str(v), "%.6f" % math.radians(d), path] + extra,
                           check=True)
            s, ang, s_ac = measure(path)
            rows.append((v, d, s, s_ac, ang))

    print("mode=%s" % mode)
    print("%-8s %-10s %-16s %-20s %-14s"
          % ("speed", "dir_deg", "spacing_fft_px", "spacing_autocorr_px", "dir_meas_deg"))
    for v, d, s, s_ac, ang in rows:
        print("%-8.1f %-10.0f %-16.2f %-20.2f %-14.2f" % (v, d, s, s_ac, ang))

    by = {}
    for v, d, s, _, _ in rows:
        by.setdefault(d, {})[v] = s

    print("\n-- independence checks --")
    monotone_ok = True
    for d in DIRS_DEG:
        vals = [by[d][v] for v in SPEEDS]
        ok = vals[0] < vals[1] < vals[2]
        monotone_ok &= ok
        print("spacing monotone in speed (dir=%3.0f deg): %s  %s"
              % (d, "PASS" if ok else "FAIL",
                 " < ".join("%.2f" % x for x in vals)))

    spread_ok = True
    for v in SPEEDS:
        vals = [by[d][v] for d in DIRS_DEG]
        spread = max(vals) - min(vals)
        ok = spread < 0.5
        spread_ok &= ok
        print("spacing invariant to direction (speed=%4.1f m/s): %s  spread=%.2f px %s"
              % (v, "PASS" if ok else "FAIL", spread,
                 "[" + ", ".join("%.2f" % x for x in vals) + "]"))

    print("\ndirection check: measured vs requested (circular diff, deg)")
    dir_ok = True
    for v, d, _, _, ang in rows:
        diff = abs((ang - d + 90) % 180 - 90)
        ok = diff < 3.0
        dir_ok &= ok
        if not ok:
            print("  speed=%.1f req=%.0f got=%.2f diff=%.2f FAIL" % (v, d, ang, diff))
    print("all directions within 3 deg: %s" % ("PASS" if dir_ok else "FAIL"))

    all_ok = monotone_ok and spread_ok and dir_ok
    print("\nCOUPLING_DETECTED=%s" % ("no" if all_ok else "yes"))
    print("all_ok=%s" % ("1" if all_ok else "0"))
    return 0 if all_ok else 1


def measure(path):
    sys.path.insert(0, HERE)
    import measure as m
    return m.measure(path)


if __name__ == "__main__":
    sys.exit(main())
