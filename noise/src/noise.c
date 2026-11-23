
/* noise: 2D value noise + fBm octave accumulation, deterministic hash. */
#include <stdio.h>
#include <math.h>
#include <string.h>

static unsigned hash2(int x, int y, unsigned seed) {
    unsigned h = (unsigned)x * 374761393u + (unsigned)y * 668265263u + seed * 2246822519u;
    h = (h ^ (h >> 13)) * 1274126177u;
    return h ^ (h >> 16);
}

static float rand01(unsigned h) { return (float)(h & 0xffffff) / 16777215.0f; }

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

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
