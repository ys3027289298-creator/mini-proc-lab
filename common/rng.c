#include "rng.h"

#include <math.h>

#include "seed.h"

uint64_t sp_hash2d(int x, int y, uint64_t seed) {
    uint64_t h = seed;
    h ^= (uint64_t)(uint32_t)x * UINT64_C(0xD6E8FEB86659FD93);
    h ^= (uint64_t)(uint32_t)y * UINT64_C(0xA24BAED4963EE407);
    return sp_mix64(h);
}

float sp_unit(uint64_t h) {
    return (float)((h >> 40) & 0xFFFFFFu) * (1.0f / 16777216.0f);
}

float sp_lerp(float a, float b, float t) { return a + (b - a) * t; }

static float fade(float t) { return t * t * (3.0f - 2.0f * t); }

float sp_value_noise(float x, float y, uint64_t seed) {
    int xi = (int)floorf(x), yi = (int)floorf(y);
    float fx = x - (float)xi, fy = y - (float)yi;
    float v00 = sp_unit(sp_hash2d(xi, yi, seed));
    float v10 = sp_unit(sp_hash2d(xi + 1, yi, seed));
    float v01 = sp_unit(sp_hash2d(xi, yi + 1, seed));
    float v11 = sp_unit(sp_hash2d(xi + 1, yi + 1, seed));
    float u = fade(fx), v = fade(fy);
    return sp_lerp(sp_lerp(v00, v10, u), sp_lerp(v01, v11, u), v) * 2.0f - 1.0f;
}

float sp_fbm(float x, float y, int octaves, float lacunarity, float gain,
             uint64_t seed) {
    if (octaves <= 0) return 0.0f;
    float sum = 0.0f, amp = 1.0f, norm = 0.0f, freq = 1.0f;
    for (int i = 0; i < octaves; i++) {
        uint64_t oct_seed =
            sp_mix64(seed + (uint64_t)i * UINT64_C(0x9E3779B97F4A7C15));
        sum += sp_value_noise(x * freq, y * freq, oct_seed) * amp;
        norm += amp;
        amp *= gain;
        freq *= lacunarity;
    }
    return norm > 0.0f ? sum / norm : 0.0f;
}

float sp_clampf(float v, float lo, float hi) {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}
