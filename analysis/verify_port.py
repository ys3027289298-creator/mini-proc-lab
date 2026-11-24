"""Verify the NumPy value-noise/fBm port against the C implementation."""
import subprocess
from pathlib import Path

import numpy as np

from noise_bases import value_field, fbm_field

HERE = Path(__file__).resolve().parent
BIN = HERE / "c_verify" / "port_harness"

subprocess.run(["cc", "-O2", "-std=c11", "-Wall", "-Wextra",
                str(HERE / "c_verify" / "port_harness.c"),
                "-o", str(BIN), "-lm"], check=True)
out = subprocess.run([str(BIN)], check=True, capture_output=True, text=True).stdout.splitlines()

it = iter(out)
assert next(it) == "value"
max_err_v, max_err_f = 0.0, 0.0
for _ in range(40):
    x, y, vc = map(float, next(it).split())
    vp = float(value_field(np.array(x, np.float32), np.array(y, np.float32), 7))
    max_err_v = max(max_err_v, abs(vp - vc))

assert next(it) == "fbm"
n_checked = 0
for line in it:
    n, x, y, fc = map(float, line.split())
    fp = float(fbm_field(np.array(x, np.float32), np.array(y, np.float32),
                         int(n), 2.0, 0.5, 11, "value"))
    max_err_f = max(max_err_f, abs(fp - fc))
    n_checked += 1

print(f"points checked: value=40 fbm={n_checked}")
print(f"max |py - c|  value_noise={max_err_v:.3e}  fbm={max_err_f:.3e}")
assert max_err_v < 1e-6 and max_err_f < 1e-6, "port mismatch"
print("PORT OK")
