
/* noise: 2D value noise + fBm octave accumulation, deterministic hash. */
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

/*
 * 64-bit avalanche mixer (splitmix64 finalizer).
 *
 * The previous 32-bit hash fed rand01 from its LOW 24 bits, but a single
 * multiply by an odd constant only diffuses differences upward: sparse
 * lattice displacements (e.g. the delta 0x0c800405 reached at cell lag
 * (43,214)) kept correlation 0.19, and the final x ^= x>>16 is the identity
 * on output bits 16..23. That produced repeating tiles in the noise field.
 * A 64-bit mix with three avalanche stages and HIGH-bit output removes the
 * short-lag structure (measured residual max ~15*sigma with no seed-stable
 * lag, vs 350*sigma peaks before).
 */
static uint64_t mix64(uint64_t z) {
    z += 0x9e3779b97f4a7c15ULL;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static uint64_t hash2(int x, int y, unsigned seed) {
    uint64_t h = (uint64_t)(int64_t)x * 0xd6e8feb86659fd93ULL
               + (uint64_t)(int64_t)y * 0xa24baed4963ee407ULL
               + (uint64_t)seed * 0x9fb21c651e98df25ULL;
    return mix64(h);
}

static float rand01(uint64_t h) {
    /* high 24 bits: every output bit depends on every input bit (~0.5 flip) */
    return (float)((h >> 40) & 0xffffffULL) / 16777215.0f;
}

static float lerp(float a, float b, float t) { return a + (b - a) * t; }
static float fade(float t) { return t * t * (3.0f - 2.0f * t); }

static float value_noise(float x, float y, unsigned seed) {
    int xi = (int)floorf(x), yi = (int)floorf(y);
    float fx = x - xi, fy = y - yi;
    float v00 = rand01(hash2(xi, yi, seed));
    float v10 = rand01(hash2(xi + 1, yi, seed));
    float v01 = rand01(hash2(xi, yi + 1, seed));
    float v11 = rand01(hash2(xi + 1, yi + 1, seed));
    float u = fade(fx), v = fade(fy);
    return lerp(lerp(v00, v10, u), lerp(v01, v11, u), v) * 2.0f - 1.0f;
}

static float fbm(float x, float y, int octaves, float lacunarity, float gain, unsigned seed) {
    if (octaves <= 0) return 0.0f;
    float sum = 0.0f, amp = 1.0f, norm = 0.0f, freq = 1.0f;
    for (int i = 0; i < octaves; i++) {
        sum += value_noise(x * freq, y * freq, seed + (unsigned)i * 17u) * amp;
        norm += amp;
        amp *= gain;
        freq *= lacunarity;
    }
    return norm > 0.0f ? sum / norm : 0.0f;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* range and determinism */
    float lo = 1e9f, hi = -1e9f;
    for (int i = 0; i < 200; i++) {
        float v = value_noise(i * 0.37f, i * 0.11f, 7u);
        if (v < lo) lo = v;
        if (v > hi) hi = v;
    }
    check("range", lo >= -1.01f && hi <= 1.01f);
    check("deterministic", value_noise(1.3f, 2.7f, 7u) == value_noise(1.3f, 2.7f, 7u));
    check("seed_changes", value_noise(1.3f, 2.7f, 7u) != value_noise(1.3f, 2.7f, 8u));

    /* lattice points must vanish (value noise is 0 at integer grid after remap? check continuity instead) */
    float a = value_noise(0.5f, 0.5f, 3u);
    float b = value_noise(0.5f + 1e-4f, 0.5f, 3u);
    check("continuous", fabsf(a - b) < 0.01f);

    float f1 = fbm(0.3f, 0.4f, 4, 2.0f, 0.5f, 11u);
    float f2 = fbm(0.3f, 0.4f, 4, 2.0f, 0.5f, 11u);
    check("fbm_deterministic", f1 == f2);
    check("fbm_range", f1 >= -1.01f && f1 <= 1.01f);
    check("fbm_zero_octaves", fbm(0.3f, 0.4f, 0, 2.0f, 0.5f, 11u) == 0.0f);

    /* more octaves => higher variance of first differences */
    float d1 = 0.0f, d2 = 0.0f;
    for (int i = 0; i < 100; i++) {
        float x = i * 0.05f;
        d1 += fabsf(fbm(x, 0.5f, 1, 2.0f, 0.5f, 5u) - fbm(x + 0.01f, 0.5f, 1, 2.0f, 0.5f, 5u));
        d2 += fabsf(fbm(x, 0.5f, 5, 2.0f, 0.5f, 5u) - fbm(x + 0.01f, 0.5f, 5, 2.0f, 0.5f, 5u));
    }
    check("more_octaves_rougher", d2 > d1);

    /*
     * Anti-period regression: the old hash showed strong lattice
     * autocorrelation at short cell lags, most prominently (43,214) with
     * corr 0.17 and (91,10) with corr 0.10 (biased FFT estimator on a
     * 2048^2 field; direct estimator 0.19 / 0.11). With the 64-bit mixer
     * the direct estimator on an 800^2 lattice must stay under 0.02
     * (noise floor for n^2 pairs is ~1/n = 0.00125).
     */
    enum { LN = 800 };
    static float lattice[LN * LN];
    double mean = 0.0;
    for (int j = 0; j < LN; j++)
        for (int i = 0; i < LN; i++) {
            float r = rand01(hash2(i, j, 7u));
            lattice[j * LN + i] = r;
            mean += r;
        }
    mean /= (double)LN * LN;
    double var = 0.0;
    for (int k = 0; k < LN * LN; k++) {
        double d = lattice[k] - mean;
        var += d * d;
    }
    var /= (double)LN * LN;
    double worst = 0.0;
    const int bad_lags[][2] = {{43, 214}, {91, 10}, {430, 494}, {215, 247},
                               {86, 428}, {48, -204}, {99, -155}};
    for (unsigned q = 0; q < sizeof(bad_lags) / sizeof(bad_lags[0]); q++) {
        int dx = bad_lags[q][0], dy = bad_lags[q][1];
        double cov = 0.0;
        long cnt = 0;
        int y0 = dy > 0 ? dy : 0, y1 = dy > 0 ? LN : LN + dy;
        int x0 = dx > 0 ? dx : 0, x1 = dx > 0 ? LN : LN + dx;
        for (int j = y0; j < y1; j++)
            for (int i = x0; i < x1; i++) {
                cov += (lattice[j * LN + i] - mean)
                     * (lattice[(j - dy) * LN + (i - dx)] - mean);
                cnt++;
            }
        double corr = cov / ((double)cnt * var);
        if (fabs(corr) > worst) worst = fabs(corr);
    }
    check("no_lattice_periodicity", worst < 0.02);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
