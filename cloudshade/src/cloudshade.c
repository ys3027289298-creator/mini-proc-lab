/* cloudshade: cloud thickness -> transmittance attenuation over terrain.
 *
 * Model (full derivation in README.md):
 *   clear-air extinction   beta_a(z) = BETA0 * exp(-z / HSCALE)
 *   slant optical depth    tau(d, z, mu) = (BETA0*HSCALE*(exp(-z/HSCALE) - exp(-ZTOP/HSCALE)) + KAPPA*d) / mu
 *   transmittance          T = exp(-tau)
 *   surface irradiance     E = S0 * mu * T + AMBIENT        (horizontal surface)
 * with mu = cos(zenith). Zenith angle and elevation both enter tau.
 *
 * The renderer integrates beta(s) numerically along the slant path;
 * the closed form above is the derived value. Self-test requires
 * max relative deviation of rendered vs derived irradiance <= 15%.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define GW 128
#define GH 128
#define CELL 30.0f          /* m, ground cell size */

#define S0 1361.0f          /* W/m^2, solar constant */
#define BETA0 1.75e-5f      /* 1/m, sea-level clear-air extinction */
#define HSCALE 8000.0f      /* m, density scale height */
#define KAPPA 0.004f        /* 1/m, in-cloud extinction */
#define ZBASE 2000.0f       /* m, cloud base */
#define ZTOP 30000.0f       /* m, top of integration */
#define AMBIENT 40.0f       /* W/m^2, diffuse floor */
#define ZENITH_DEG 40.0f    /* sun zenith angle for rendering */

#define TIER_N 3
static const float TIER_D[TIER_N] = { 150.0f, 600.0f, 1500.0f }; /* m, nominal thickness */

#define MAX_REL_DEV 0.15f   /* acceptance: rendered vs derived within 15% */

/* ---- noise (self-contained copy, same hash as noise/) ---- */
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

/* ---- fields ---- */
static float terrain_z(float wx, float wy) {
    float n = fbm(wx / 900.0f, wy / 900.0f, 4, 2.0f, 0.5f, 21u);
    return 1400.0f * (0.5f + 0.5f * n);   /* 0..1400 m */
}

/* cloud thickness shape in [0, 1.6], mean ~0.55, with clear gaps.
 * toff slides the cloud field over the ground (cloud drift). */
static float cloud_shape(float wx, float wy, float toff) {
    float n = fbm((wx + toff) / 1100.0f, wy / 1100.0f, 4, 2.0f, 0.5f, 77u);
    float s = 0.55f + 0.9f * n;
    if (s < 0.0f) s = 0.0f;
    if (s > 1.6f) s = 1.6f;
    return s;
}

/* ---- derived (closed form) ---- */
static float tau_closed(float d, float z, float mu) {
    float atm = BETA0 * HSCALE * (expf(-z / HSCALE) - expf(-ZTOP / HSCALE));
    return (atm + KAPPA * d) / mu;
}

static float irr_closed(float d, float z, float mu) {
    return S0 * mu * expf(-tau_closed(d, z, mu)) + AMBIENT;
}

/* ---- rendered (numerical quadrature along the slant path) ---- */
static float tau_render(float d, float z, float mu) {
    /* clear air: midpoint rule over the whole slant path (smooth integrand) */
    const int n = 256;
    float L = (ZTOP - z) / mu;      /* slant path length */
    float ds = L / n;
    float sum = 0.0f;
    for (int i = 0; i < n; i++) {
        float s = (i + 0.5f) * ds;
        sum += BETA0 * expf(-(z + s * mu) / HSCALE) * ds;
    }
    /* cloud slab: midpoint rule over its own slant interval [s1, s2] */
    float s1 = (ZBASE - z) / mu, s2 = (ZBASE + d - z) / mu;
    if (s2 > 0.0f && s1 < L) {
        if (s1 < 0.0f) s1 = 0.0f;
        if (s2 > L) s2 = L;
        const int m = 128;
        float dsc = (s2 - s1) / m;
        for (int i = 0; i < m; i++) {
            float s = s1 + (i + 0.5f) * dsc;
            (void)s;                /* in-cloud extinction is constant */
            sum += KAPPA * dsc;
        }
    }
    return sum;
}

static float irr_render(float d, float z, float mu) {
    return S0 * mu * expf(-tau_render(d, z, mu)) + AMBIENT;
}

/* ---- self-test ---- */
static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static int write_pgm(const char* path, const float* e, float eref) {
    FILE* f = fopen(path, "wb");
    if (!f) return 0;
    fprintf(f, "P5\n%d %d\n255\n", GW, GH);
    for (int y = 0; y < GH; y++)
        for (int x = 0; x < GW; x++) {
            float v = e[y * GW + x] / eref;
            if (v < 0.0f) v = 0.0f;
            if (v > 1.0f) v = 1.0f;
            fputc((int)(v * 255.0f + 0.5f), f);
        }
    fclose(f);
    return 1;
}

static void print_ascii(const float* e, float eref) {
    const char* ramp = " .:-=+*#%@";
    for (int y = 0; y < GH; y += 8) {
        for (int x = 0; x < GW; x += 2) {
            float v = e[y * GW + x] / eref;
            if (v < 0.0f) v = 0.0f;
            if (v > 1.0f) v = 1.0f;
            putchar(ramp[(int)(v * 9.0f + 0.5f)]);
        }
        putchar('\n');
    }
}

int main(void) {
    float mu = cosf(ZENITH_DEG * (float)M_PI / 180.0f);
    float eref = irr_closed(0.0f, 0.0f, mu);   /* clear sky, sea level: map white point */

    static float z[GW * GH];
    for (int y = 0; y < GH; y++)
        for (int x = 0; x < GW; x++)
            z[y * GW + x] = terrain_z(x * CELL, y * CELL);

    /* monotonicity of the derived attenuation in d, zenith, elevation */
    check("thicker_cloud_darker",
          irr_closed(150.0f, 0.0f, mu) > irr_closed(600.0f, 0.0f, mu) &&
          irr_closed(600.0f, 0.0f, mu) > irr_closed(1500.0f, 0.0f, mu));
    check("higher_zenith_darker",
          irr_closed(600.0f, 0.0f, cosf(0.0f)) > irr_closed(600.0f, 0.0f, mu) &&
          irr_closed(600.0f, 0.0f, mu) > irr_closed(600.0f, 0.0f, cosf(1.0f)));
    check("higher_elevation_brighter",
          irr_closed(600.0f, 3000.0f, mu) > irr_closed(600.0f, 0.0f, mu));

    printf("\n# three-tier surface irradiance table (closed form, W/m^2)\n");
    printf("tier  d_nom(m)  z(m)   E(theta=0)  E(theta=40)  E(theta=60)\n");
    for (int t = 0; t < TIER_N; t++)
        for (int zi = 0; zi < 3; zi++) {
            float zz = zi * 1500.0f;
            printf("%4d  %8.0f  %5.0f  %10.1f  %11.1f  %11.1f\n",
                   t + 1, (double)TIER_D[t], (double)zz,
                   (double)irr_closed(TIER_D[t], zz, cosf(0.0f)),
                   (double)irr_closed(TIER_D[t], zz, mu),
                   (double)irr_closed(TIER_D[t], zz, cosf(60.0f * (float)M_PI / 180.0f)));
        }

    printf("\n# rendered maps (theta=%.0f deg, one render per tier)\n", (double)ZENITH_DEG);
    float worst_dev = 0.0f;
    int pgm_ok = 1;
    static float emap[GW * GH];
    for (int t = 0; t < TIER_N; t++) {
        float emin = 1e9f, emax = -1e9f, esum = 0.0f, dsum = 0.0f, dev = 0.0f;
        for (int y = 0; y < GH; y++)
            for (int x = 0; x < GW; x++) {
                float zz = z[y * GW + x];
                float d = TIER_D[t] * cloud_shape(x * CELL, y * CELL, 0.0f);
                float er = irr_render(d, zz, mu);
                float ec = irr_closed(d, zz, mu);
                float rel = fabsf(er - ec) / ec;
                if (rel > dev) dev = rel;
                emap[y * GW + x] = er;
                if (er < emin) emin = er;
                if (er > emax) emax = er;
                esum += er;
                dsum += d;
            }
        int n = GW * GH;
        if (dev > worst_dev) worst_dev = dev;
        char path[128];
        snprintf(path, sizeof(path), "../artifacts/cloudshade_tier%d.pgm", t + 1);
        pgm_ok &= write_pgm(path, emap, eref);
        printf("tier %d: d_mean=%.0f m  E_mean=%.1f  E_min=%.1f  E_max=%.1f W/m^2  max_rel_dev=%.4f\n",
               t + 1, (double)(dsum / n), (double)(esum / n), (double)emin, (double)emax, (double)dev);
        print_ascii(emap, eref);
    }

    check("render_matches_derived_within_15pct", worst_dev <= MAX_REL_DEV);
    check("pgm_written", pgm_ok);

    /* the original complaint: shadow must modulate intensity, not just slide */
    check("intensity_varies_under_cloud",
          irr_render(TIER_D[0] * 1.6f, 0.0f, mu) < 0.9f * irr_render(0.0f, 0.0f, mu));

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
