#!/usr/bin/env python3
"""Measure average ripple spacing and dominant direction from a PGM texture.

Two independent estimators per image:
  fft   -> spectral peak gives spacing (px) and direction (deg)
  autocorr -> real-space cross-section gives a second spacing estimate
Run under: python3 -I measure.py IMAGE.pgm
"""
import sys
import site

# Isolated mode (-I) drops the user site-packages directory; re-add it
# explicitly so numpy resolves while cwd and PYTHON* paths stay excluded.
_usp = site.getusersitepackages()
if _usp not in sys.path:
    sys.path.append(_usp)

import numpy as np


def read_pgm(path):
    with open(path, "rb") as f:
        magic = f.readline().strip()
        if magic != b"P5":
            raise ValueError("only binary P5 PGM supported: %s" % magic)
        line = f.readline()
        while line.startswith(b"#"):
            line = f.readline()
        w, h = (int(v) for v in line.split())
        maxv = int(f.readline())
        if maxv != 255:
            raise ValueError("only maxval=255 supported")
        buf = f.read(w * h)
    return np.frombuffer(buf, dtype=np.uint8).reshape(h, w).astype(np.float64)


def fft_measure(img):
    n = img.shape[0]
    win = np.hanning(n)[:, None] * np.hanning(n)[None, :]
    d = (img - img.mean()) * win
    spec = np.fft.fftshift(np.fft.fft2(d))
    mag = np.abs(spec)
    cy = cx = n // 2
    mag[cy, cx] = 0.0                      # drop DC
    mag[max(cy - 2, 0):cy + 3, max(cx - 2, 0):cx + 3] = 0.0
    iy, ix = np.unravel_index(np.argmax(mag), mag.shape)

    # parabolic sub-bin interpolation along both axes
    def sub(v0, v1, v2):
        den = (v0 - 2.0 * v1 + v2)
        return 0.5 * (v0 - v2) / den if den != 0 else 0.0

    fx = (ix - cx) + sub(mag[iy, ix - 1], mag[iy, ix], mag[iy, ix + 1])
    fy = (iy - cy) + sub(mag[iy - 1, ix], mag[iy, ix], mag[iy + 1, ix])
    spacing = n / np.hypot(fx, fy)
    angle = np.degrees(np.arctan2(fy, fx)) % 180.0
    return spacing, angle


def _bilinear(a, x, y):
    x0 = int(np.floor(x))
    y0 = int(np.floor(y))
    tx, ty = x - x0, y - y0
    h, w = a.shape
    x0c, x1c = max(x0, 0), min(x0 + 1, w - 1)
    y0c, y1c = max(y0, 0), min(y0 + 1, h - 1)
    return (a[y0c, x0c] * (1 - tx) * (1 - ty) + a[y0c, x1c] * tx * (1 - ty)
            + a[y1c, x0c] * (1 - tx) * ty + a[y1c, x1c] * tx * ty)


def autocorr_spacing(img, direction_deg):
    """Cross-section of the autocorrelation along the FFT peak direction;
    first interior maximum = average real-space spacing."""
    n = img.shape[0]
    d = img - img.mean()
    f = np.fft.fft2(d)
    ac = np.fft.ifft2(f * np.conj(f)).real / (d * d).sum()
    ac = np.fft.fftshift(ac)

    a = np.radians(direction_deg)
    dx, dy = np.cos(a), np.sin(a)
    r = np.arange(1, n // 2)
    prof = np.array([_bilinear(ac, n // 2 + t * dx, n // 2 + t * dy) for t in r])

    # first interior maximum: profile turns from rising to falling
    for i in range(1, len(prof) - 1):
        if prof[i] > prof[i - 1] and prof[i] >= prof[i + 1]:
            # sub-pixel parabola fit around the peak
            v0, v1, v2 = prof[i - 1], prof[i], prof[i + 1]
            den = v0 - 2 * v1 + v2
            off = 0.5 * (v0 - v2) / den if den != 0 else 0.0
            return r[i] + off
    return float("nan")


def measure(path):
    img = read_pgm(path)
    spacing, angle = fft_measure(img)
    spacing_ac = autocorr_spacing(img, angle)
    return spacing, angle, spacing_ac


def main():
    for p in sys.argv[1:]:
        s, a, s_ac = measure(p)
        print("%s  fft_spacing_px=%.2f  autocorr_spacing_px=%.2f  direction_deg=%.2f"
              % (p, s, s_ac, a))


if __name__ == "__main__":
    main()
