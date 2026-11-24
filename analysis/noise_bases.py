"""Two 2D noise bases sharing one fBm framework.

Basis A: value noise, a faithful NumPy port of noise/src/noise.c
         (same uint32 hash, smoothstep fade 3t^2-2t^3, amplitude-normalized fBm).
Basis B: gradient (Perlin) noise, same hash family, quintic fade, identical fBm.
"""
import numpy as np

__all__ = ["value_field", "gradient_field", "fbm_field", "fbm_partials"]

_UMAX = np.uint32(0xFFFFFFFF)


def _hash2(x, y, seed):
    """Port of hash2() in noise.c, vectorized over integer lattice coords."""
    x = np.asarray(x, dtype=np.int64).astype(np.uint32)
    y = np.asarray(y, dtype=np.int64).astype(np.uint32)
    seed = np.uint32(seed)
    with np.errstate(over="ignore"):  # uint32 wraparound is the intended hash math
        h = (x * np.uint32(374761393)
             + y * np.uint32(668265263)
             + seed * np.uint32(2246822519))
        h = (h ^ (h >> np.uint32(13))) * np.uint32(1274126177)
        return h ^ (h >> np.uint32(16))


def _rand01(h):
    return (h & np.uint32(0xFFFFFF)).astype(np.float32) / np.float32(16777215.0)


def _fade3(t):
    return t * t * (np.float32(3.0) - np.float32(2.0) * t)


def _fade5(t):
    return t * t * t * (t * (t * np.float32(6.0) - np.float32(15.0)) + np.float32(10.0))


def _corners(xi, yi, fx, fy, vals, fade):
    u = fade(fx)
    v = fade(fy)
    x0 = vals[0, 0] * (1 - u) + vals[1, 0] * u
    x1 = vals[0, 1] * (1 - u) + vals[1, 1] * u
    return x0 * (1 - v) + x1 * v


def value_field(X, Y, seed):
    """Value noise in [-1, 1]; matches value_noise() in noise.c."""
    xi = np.floor(X).astype(np.int64)
    yi = np.floor(Y).astype(np.int64)
    fx = (X - xi).astype(np.float32)
    fy = (Y - yi).astype(np.float32)
    vals = np.empty((2, 2) + X.shape, dtype=np.float32)
    for j in range(2):
        for i in range(2):
            vals[i, j] = _rand01(_hash2(xi + i, yi + j, seed))
    return _corners(xi, yi, fx, fy, vals, _fade3) * np.float32(2.0) - np.float32(1.0)


_GRADS = np.array([(1, 1), (-1, 1), (1, -1), (-1, -1),
                   (1, 0), (-1, 0), (0, 1), (0, -1)], dtype=np.float32)


def gradient_field(X, Y, seed):
    """Perlin gradient noise, same hash family, quintic fade."""
    xi = np.floor(X).astype(np.int64)
    yi = np.floor(Y).astype(np.int64)
    fx = (X - xi).astype(np.float32)
    fy = (Y - yi).astype(np.float32)
    vals = np.empty((2, 2) + X.shape, dtype=np.float32)
    for j in range(2):
        for i in range(2):
            h = _hash2(xi + i, yi + j, seed) & np.uint32(7)
            g = _GRADS[h]
            vals[i, j] = g[..., 0] * (fx - np.float32(i)) + g[..., 1] * (fy - np.float32(j))
    return _corners(xi, yi, fx, fy, vals, _fade5)


_BASES = {"value": value_field, "gradient": gradient_field}


def fbm_field(X, Y, octaves, lacunarity=2.0, gain=0.5, seed=11, basis="value"):
    """Amplitude-normalized fBm, matching fbm() in noise.c (float32 accumulation)."""
    if octaves <= 0:
        return np.zeros_like(X, dtype=np.float32)
    base = _BASES[basis]
    total = np.zeros_like(X, dtype=np.float32)
    norm = np.float32(0.0)
    amp = np.float32(1.0)
    freq = np.float32(1.0)
    L = np.float32(lacunarity)
    r = np.float32(gain)
    for i in range(octaves):
        total += base(X * freq, Y * freq, seed + i * 17) * amp
        norm += amp
        amp *= r
        freq *= L
    return total / norm


def fbm_partials(X, Y, max_octaves, lacunarity=2.0, gain=0.5, seed=11, basis="value"):
    """Return normalized fBm fields for N = 1..max_octaves built incrementally.

    One noise evaluation per octave; partial sums re-used for every N.
    """
    base = _BASES[basis]
    raw_sum = np.zeros_like(X, dtype=np.float32)
    norm = np.float32(0.0)
    amp = np.float32(1.0)
    freq = np.float32(1.0)
    L = np.float32(lacunarity)
    r = np.float32(gain)
    out = []
    for i in range(max_octaves):
        raw_sum += base(X * freq, Y * freq, seed + i * 17) * amp
        norm += amp
        out.append((raw_sum / norm).copy())
        amp *= r
        freq *= L
    return out
