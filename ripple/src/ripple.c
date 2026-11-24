
/* ripple: wind-driven sand ripple height field.
 *
 * Two independent parameter maps:
 *   wind speed -> ripple wavelength (spacing): wavelength()
 *   wind dir   -> phase offset:               phase_from_dir()
 * Direction also rotates the wave vector (pure geometry) but never touches
 * the wavelength map. The optional "coupled" mode injects a direction-
 * dependent gust factor into the spacing map on purpose, so the measurement
 * pipeline has a coupling bug to detect.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define SZ 256
#define TAU 6.283185307179586f

/* ---- deterministic value noise (same hash family as noise/) ---- */
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

/* ---- map 1: wind speed (m/s) -> ripple wavelength (px) ----
 * Bagnold-inspired: spacing grows linearly with excess wind speed above
 * the transport threshold. Pure function of speed alone. */
#define V_TH 4.0f
#define LAM0 8.0f
#define K_LAM 3.75f

static float wavelength(float speed) {
    float v = speed > V_TH ? speed : V_TH;
    return LAM0 + K_LAM * (v - V_TH);
}

/* ---- map 2: wind direction (rad) -> phase offset (rad) ----
 * Pure additive phase, function of direction alone. It shifts the pattern
 * but never enters the wave vector, so it cannot change the spacing. */
static float phase_from_dir(float dir) {
    return 0.5f * sinf(2.0f * dir);
}

/* Spacing actually used by the field. The fixed path is map 1 only.
 * Coupled path (demo bug): a gust factor leaks direction into the spacing
 * map -- this is the coupling source the measurement step must expose. */
static float ripple_lambda(float speed, float dir, int coupled) {
    float lam = wavelength(speed);
    if (coupled) {
        float gust = 1.0f + (0.10f + 0.20f * (speed - 7.0f)) * sinf(3.0f * dir);
        lam *= gust;
    }
    return lam;
}

static float ripple_height(float x, float y, float speed, float dir,
                           unsigned seed, int coupled) {
    float lam = ripple_lambda(speed, dir, coupled);
    float kx = cosf(dir), ky = sinf(dir);   /* direction rotates the wave vector */
    float along = kx * x + ky * y;
    float phi = phase_from_dir(dir);        /* direction only adds a phase offset */
    float h = 0.8f * sinf(TAU * along / lam + phi);
    h += 0.12f * fbm(x * 0.04f, y * 0.04f, 3, 2.0f, 0.5f, seed);
    return h;
}

static int write_pgm(const char *path, const unsigned char *buf) {
    FILE *f = fopen(path, "wb");
    if (!f) return -1;
    fprintf(f, "P5\n%d %d\n255\n", SZ, SZ);
    if (fwrite(buf, 1, (size_t)SZ * SZ, f) != (size_t)SZ * SZ) {
        fclose(f);
        return -1;
    }
    fclose(f);
    return 0;
}

static int gen(const char *speed_s, const char *dir_s, const char *out, int coupled) {
    float speed = strtof(speed_s, NULL);
    float dir = strtof(dir_s, NULL);
    unsigned char *buf = malloc((size_t)SZ * SZ);
    if (!buf) return 1;
    for (int y = 0; y < SZ; y++) {
        for (int x = 0; x < SZ; x++) {
            float h = ripple_height((float)x, (float)y, speed, dir, 1234u, coupled);
            float t = h / 0.95f * 0.5f + 0.5f;
            if (t < 0.0f) t = 0.0f;
            if (t > 1.0f) t = 1.0f;
            buf[y * SZ + x] = (unsigned char)(t * 255.0f + 0.5f);
        }
    }
    int rc = write_pgm(out, buf);
    free(buf);
    if (rc != 0) {
        perror("write_pgm");
        return 1;
    }
    printf("wrote %s speed=%.2f dir=%.4f coupled=%d\n", out, speed, dir, coupled);
    return 0;
}

static int failures = 0;
static void check(const char *name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(int argc, char **argv) {
    if (argc >= 2 && strcmp(argv[1], "gen") == 0) {
        if (argc < 5) {
            fprintf(stderr, "usage: %s gen SPEED DIR_RAD OUTFILE [coupled]\n", argv[0]);
            return 2;
        }
        int coupled = argc >= 6 && strcmp(argv[5], "coupled") == 0;
        return gen(argv[2], argv[3], argv[4], coupled);
    }

    /* map 1: spacing is monotone in speed and clamped at the threshold */
    check("wavelength_monotonic",
          wavelength(4.0f) < wavelength(7.0f) && wavelength(7.0f) < wavelength(10.0f));
    check("wavelength_threshold_clamp",
          fabsf(wavelength(2.0f) - wavelength(4.0f)) < 1e-6f);

    /* map 2: phase is a function of direction alone, periodic with pi */
    check("phase_varies_with_dir",
          fabsf(phase_from_dir(0.0f) - phase_from_dir(0.7854f)) > 1e-3f);
    check("phase_periodic_pi",
          fabsf(phase_from_dir(0.3f) - phase_from_dir(0.3f + 3.14159265f)) < 1e-5f);

    /* independence: fixed spacing path ignores direction entirely */
    float l0 = ripple_lambda(7.0f, 0.0f, 0);
    float l1 = ripple_lambda(7.0f, 0.7854f, 0);
    float l2 = ripple_lambda(7.0f, 1.5708f, 0);
    check("spacing_dir_invariant", l0 == l1 && l1 == l2);

    /* the deliberately coupled path must show direction leaking in */
    float c0 = ripple_lambda(7.0f, 0.0f, 1);
    float c1 = ripple_lambda(7.0f, 1.5708f, 1);
    check("coupled_leaks_dir_into_spacing", c0 != c1);

    /* geometry: dir 0 varies along x, dir pi/2 along y; same wave part */
    int orient_ok = 1;
    for (int i = 0; i < 40; i++) {
        float hx = ripple_height((float)i, 0.0f, 7.0f, 0.0f, 1u, 0);
        float hy = ripple_height(0.0f, (float)i, 7.0f, 1.5707963f, 1u, 0);
        if (fabsf(hx - hy) > 0.25f) { orient_ok = 0; break; }
    }
    check("orientation_rotates_field", orient_ok);

    /* output range and determinism */
    float lo = 1e9f, hi = -1e9f;
    for (int i = 0; i < 500; i++) {
        float h = ripple_height(i * 0.37f, i * 0.11f, 7.0f, 0.5f, 9u, 0);
        if (h < lo) lo = h;
        if (h > hi) hi = h;
    }
    check("height_range", lo > -1.0f && hi < 1.0f);
    float a = ripple_height(3.1f, 2.7f, 7.0f, 0.5f, 9u, 0);
    float b = ripple_height(3.1f, 2.7f, 7.0f, 0.5f, 9u, 0);
    check("deterministic", a == b);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
