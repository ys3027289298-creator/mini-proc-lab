#include "stones.h"

#include <math.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* shared basis                                                        */
/* ------------------------------------------------------------------ */

static unsigned mix32(unsigned h) {
    h ^= h >> 16;
    h *= 0x7feb352du;
    h ^= h >> 15;
    h *= 0x846ca68bu;
    h ^= h >> 16;
    return h;
}

static unsigned hash_xy(int x, int y, unsigned seed) {
    unsigned h = (unsigned)x * 0x9e3779b1u ^ (unsigned)y * 0x85ebca6bu ^ seed * 0xc2b2ae35u;
    return mix32(h);
}

static float u01(unsigned h) { return (float)(h >> 8) * (1.0f / 16777216.0f); }

/* Stone layout: positions chosen on a jittered grid plus a few loose
 * outliers, then a hard non-overlap pass. Every constant here is fixed. */
void stones_layout(StoneSet *set, unsigned layout_seed) {
    memset(set, 0, sizeof(*set));
    const unsigned seed = layout_seed;
    /* grid of 12 x 8 chart cells; one stone per cell, offset + permuted so the
     * arrangement reads as gravel rather than a lattice. */
    int gw = 12, gh = 8;
    for (int gy = 0; gy < gh && set->n < MAX_STONES; gy++) {
        for (int gx = 0; gx < gw && set->n < MAX_STONES; gx++) {
            unsigned h = hash_xy(gx, gy, seed);
            unsigned h2 = hash_xy(gx + 101, gy - 77, seed);

            float x = (float)gx + 0.5f + (u01(h) - 0.5f) * 0.62f;
            float y = (float)gy + 0.5f + (u01(h2) - 0.5f) * 0.62f;
            float base = 0.30f + 0.50f * u01(mix32(h + 17u));
            float rx = base * (0.86f + 0.30f * u01(mix32(h + 29u)));
            float ry = base * (0.78f + 0.32f * u01(mix32(h2 + 41u)));

            Stone *s = &set->s[set->n];
            s->x = x;
            s->y = y;
            s->rx = rx;
            s->ry = ry;
            s->rot = u01(mix32(h ^ 0x51u)) * 6.2831853f;
            for (int k = 0; k < 8; k++)
                s->e[k] = 0.86f + 0.26f * u01(mix32(h + 131u * (unsigned)(k + 1)));
            s->tone = 0.20f + 0.45f * u01(mix32(h2 ^ 0xa3u));
            s->rough = 0.55f + 0.45f * u01(mix32(h + 7u));
            s->warp = 0.20f + 0.60f * u01(mix32(h2 + 11u));
            s->variant = (int)(h % 7u);
            /* conservative bound: max radial multiplier is 1.12, radii scaled */
            s->bound = 1.16f * (rx > ry ? rx : ry) + 0.02f;
            set->n++;
        }
    }
}

static float circle_radius(const Stone *s, float ang) {
    float k = ang * (8.0f / 6.2831853f);
    int i0 = (int)floorf(k) & 7;
    int i1 = (i0 + 1) & 7;
    float t = k - floorf(k);
    float f0 = s->e[i0], f1 = s->e[i1];
    return f0 + (f1 - f0) * (t * t * (3.0f - 2.0f * t));
}

/* signed-ish coverage: >0 inside, smooth over ~0.03 units at the rim. */
static float stone_field(const Stone *s, float dx, float dy) {
    float c = cosf(s->rot), sn = sinf(s->rot);
    float lx = (dx * c + dy * sn) / s->rx;
    float ly = (-dx * sn + dy * c) / s->ry;
    float r = sqrtf(lx * lx + ly * ly);
    float ang = atan2f(ly, lx);
    float edge = circle_radius(s, ang);
    return (edge - r) / 0.045f; /* positive inside */
}

float stones_coverage(const StoneSet *set, float x, float y, int *id_out) {
    float best = -1e9f;
    int best_id = -1;
    for (int i = 0; i < set->n; i++) {
        const Stone *s = &set->s[i];
        float dx = x - s->x, dy = y - s->y;
        if (dx * dx + dy * dy > s->bound * s->bound) continue;
        float f = stone_field(s, dx, dy);
        if (f > best) { best = f; best_id = i; }
    }
    if (id_out) *id_out = best_id;
    if (best <= 0.0f) return 0.0f;
    float cov = best < 1.0f ? best : 1.0f;
    return cov * cov * (3.0f - 2.0f * cov);
}

/* ------------------------------------------------------------------ */
/* Module A: value noise + fBm                                         */
/* ------------------------------------------------------------------ */

static float fade(float t) { return t * t * (3.0f - 2.0f * t); }
static float lerpf(float a, float b, float t) { return a + (b - a) * t; }

static float value_noise_A(float x, float y, unsigned seed) {
    int xi = (int)floorf(x), yi = (int)floorf(y);
    float fx = x - (float)xi, fy = y - (float)yi;
    float v00 = u01(hash_xy(xi, yi, seed));
    float v10 = u01(hash_xy(xi + 1, yi, seed));
    float v01 = u01(hash_xy(xi, yi + 1, seed));
    float v11 = u01(hash_xy(xi + 1, yi + 1, seed));
    float u = fade(fx), w = fade(fy);
    return lerpf(lerpf(v00, v10, u), lerpf(v01, v11, u), w) * 2.0f - 1.0f;
}

static float fbm_A(float x, float y, int oct, float lac, float gain, unsigned seed) {
    float sum = 0.0f, amp = 1.0f, norm = 0.0f, freq = 1.0f;
    for (int i = 0; i < oct; i++) {
        sum += value_noise_A(x * freq, y * freq, seed + (unsigned)i * 17u) * amp;
        norm += amp;
        amp *= gain;
        freq *= lac;
    }
    return norm > 0.0f ? sum / norm : 0.0f;
}

float textureA_raw(float x, float y, unsigned seed) {
    return fbm_A(x, y, 4, 2.0f, 0.5f, seed);
}

/* ------------------------------------------------------------------ */
/* Module B: gradient noise + fBm (independent implementation)         */
/* ------------------------------------------------------------------ */

static float grad_dot(int ix, int iy, unsigned seed, float dx, float dy) {
    unsigned h = hash_xy(ix, iy, seed);
    float ang = u01(h) * 6.2831853f;
    return cosf(ang) * dx + sinf(ang) * dy;
}

static float grad_noise_B(float x, float y, unsigned seed) {
    int xi = (int)floorf(x), yi = (int)floorf(y);
    float fx = x - (float)xi, fy = y - (float)yi;
    float u = fade(fx), v = fade(fy);
    float n00 = grad_dot(xi, yi, seed, fx, fy);
    float n10 = grad_dot(xi + 1, yi, seed, fx - 1.0f, fy);
    float n01 = grad_dot(xi, yi + 1, seed, fx, fy - 1.0f);
    float n11 = grad_dot(xi + 1, yi + 1, seed, fx - 1.0f, fy - 1.0f);
    return lerpf(lerpf(n00, n10, u), lerpf(n01, n11, u), v) * 1.4142136f;
}

static float fbm_B(float x, float y, int oct, float lac, float gain, unsigned seed) {
    float sum = 0.0f, amp = 1.0f, norm = 0.0f, freq = 1.0f;
    for (int i = 0; i < oct; i++) {
        sum += grad_noise_B(x * freq, y * freq, seed ^ (0x9e3779b9u * (unsigned)(i + 1))) * amp;
        norm += amp;
        amp *= gain;
        freq *= lac;
    }
    return norm > 0.0f ? sum / norm : 0.0f;
}

float textureB_raw(float x, float y, unsigned seed) {
    /* gradient noise has a lower RMS than value noise at the same octaves;
     * this gain makes the two modules comparable in contrast, not identical. */
    return fbm_B(x, y, 4, 2.0f, 0.5f, seed) * 1.38f;
}

/* ------------------------------------------------------------------ */
/* stone surface composition                                           */
/* ------------------------------------------------------------------ */

static const float kPalette[7] = {0.34f, 0.41f, 0.48f, 0.55f, 0.62f, 0.69f, 0.76f};

static float stone_lum(const StoneSet *set, const Stone *s, float x, float y,
                       unsigned seed, int module) {
    float dx = x - s->x, dy = y - s->y;

    /* local frame for grain: rotation is part of the fixed layout */
    float c = cosf(s->rot), sn = sinf(s->rot);
    float lx = dx * c + dy * sn;
    float ly = -dx * sn + dy * c;
    float inv = 1.0f / (s->rx > 0.05f ? s->rx : 0.05f);

    /* marbling: warp the grain domain with the stone's own slow noise */
    float wseed = seed + (unsigned)(s->variant * 37 + 1);
    float wx, wy;
    if (module == 0) {
        wx = textureA_raw(lx * 1.6f + 13.0f, ly * 1.6f, wseed) * s->warp;
        wy = textureA_raw(lx * 1.6f, ly * 1.6f + 7.0f, wseed) * s->warp;
    } else {
        wx = textureB_raw(lx * 1.6f + 13.0f, ly * 1.6f, wseed) * s->warp;
        wy = textureB_raw(lx * 1.6f, ly * 1.6f + 7.0f, wseed) * s->warp;
    }

    float gx = (lx + wx * s->rx) * (5.0f * inv);
    float gy = (ly + wy * s->ry) * (5.0f * inv);
    float grain = module == 0 ? textureA_raw(gx, gy, seed + 91u)
                              : textureB_raw(gx, gy, seed + 91u);

    /* second, higher-frequency chips layer */
    float chip = module == 0 ? textureA_raw(gx * 3.1f + 4.0f, gy * 3.1f - 9.0f, seed + 211u)
                             : textureB_raw(gx * 3.1f + 4.0f, gy * 3.1f - 9.0f, seed + 211u);

    float base = kPalette[s->variant];
    float lum = base + grain * 0.20f * s->rough + chip * 0.065f * s->rough;

    /* rim darkening so pebbles read as discrete bodies */
    float r = sqrtf((dx * dx) / (s->rx * s->rx) + (dy * dy) / (s->ry * s->ry));
    float rim = (r - 0.55f) / 0.45f;
    if (rim > 0.0f) {
        float k = rim > 1.0f ? 1.0f : rim;
        lum -= 0.16f * k * k;
    }

    /* per-stone tone nudge, seeded, so stones differ beyond palette step */
    lum += (u01(hash_xy((int)(s->x * 16.0f), (int)(s->y * 16.0f), seed + 5u)) - 0.5f) * 0.09f;
    (void)set;

    if (lum < 0.0f) lum = 0.0f;
    if (lum > 1.0f) lum = 1.0f;
    return lum;
}

static float background_lum(float x, float y, unsigned seed, int module) {
    float n = module == 0 ? textureA_raw(x * 6.0f, y * 6.0f, seed + 777u)
                          : textureB_raw(x * 6.0f, y * 6.0f, seed + 777u);
    return 0.16f + 0.07f * n;
}

static float field_lum(const StoneSet *set, float x, float y, unsigned seed, int module) {
    int id = -1;
    float cov = stones_coverage(set, x, y, &id);
    float bg = background_lum(x, y, seed, module);
    if (id < 0 || cov <= 0.0f) return bg;
    float stone = stone_lum(set, &set->s[id], x, y, seed, module);
    return bg + (stone - bg) * cov;
}

float textureA_lum(const StoneSet *set, float x, float y, unsigned seed) {
    return field_lum(set, x, y, seed, 0);
}

float textureB_lum(const StoneSet *set, float x, float y, unsigned seed) {
    return field_lum(set, x, y, seed, 1);
}
