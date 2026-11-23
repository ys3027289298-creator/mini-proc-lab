
/* lod: screen-space-error LOD policy for the heightfield terrain.
 *
 * Model (derived in README.md):
 *   rho_l(d) = K * eps_l / d                          screen-space error, px
 *   d_l      = K * eps_l / tau                        switch distance of level l
 *   f_sw     = v_r / (d_l*(sqrt(r) - 1/sqrt(r)))      switch frequency, hysteresis ratio r
 *   pop_l    = tau * (1 - eps_l / eps_{l+1})          pop amplitude at a switch, px
 *
 * Renders three tier distances, a near-quality comparison pair and a
 * far frame-to-frame flicker pair; counts LOD switches along a camera
 * path for baseline / tight / tuned policies.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define N 257
#define LMAX 3
#define CELL 2.0f
#define IMG_W 640
#define IMG_H 360
#define FOV_DEG 55.0f
#define PI 3.14159265f

#define TAU_BASE 4.0f    /* baseline error budget, px: near too coarse */
#define TAU_TUNED 1.2f   /* tuned error budget, px */
#define R_TUNED 1.35f    /* hysteresis ratio */
#define SWITCH_BUDGET 120 /* max switches over the 60 s path */

#define FPS 60
#define DUR 60

/* ---- noise (same as noise/) ---- */
static unsigned hash2(int x, int y, unsigned seed) {
    unsigned h = (unsigned)x * 374761393u + (unsigned)y * 668265263u + seed * 2246822519u;
    h = (h ^ (h >> 13)) * 1274126177u;
    return h ^ (h >> 16);
}
static float rand01(unsigned h) { return (float)(h & 0xffffff) / 16777215.0f; }
static float lerpf(float a, float b, float t) { return a + (b - a) * t; }
static float fade(float t) { return t * t * (3.0f - 2.0f * t); }

static float value_noise(float x, float y, unsigned seed) {
    int xi = (int)floorf(x), yi = (int)floorf(y);
    float fx = x - xi, fy = y - yi;
    float v00 = rand01(hash2(xi, yi, seed));
    float v10 = rand01(hash2(xi + 1, yi, seed));
    float v01 = rand01(hash2(xi, yi + 1, seed));
    float v11 = rand01(hash2(xi + 1, yi + 1, seed));
    float u = fade(fx), v = fade(fy);
    return lerpf(lerpf(v00, v10, u), lerpf(v01, v11, u), v) * 2.0f - 1.0f;
}

static float fbm(float x, float y, int octaves, float lacunarity, float gain, unsigned seed) {
    float sum = 0.0f, amp = 1.0f, norm = 0.0f, freq = 1.0f;
    for (int i = 0; i < octaves; i++) {
        sum += value_noise(x * freq, y * freq, seed + (unsigned)i * 17u) * amp;
        norm += amp;
        amp *= gain;
        freq *= lacunarity;
    }
    return norm > 0.0f ? sum / norm : 0.0f;
}

/* ---- heightfield ---- */
static float Hf[N][N];
static float eps_max[LMAX + 1], eps_rms[LMAX + 1];
static float KPROJ, D2_BASE;

static float clampf(float v, float a, float b) { return v < a ? a : v > b ? b : v; }

/* bilinear sample of the level-l mesh (vertex stride 2^l) at world (wx,wz) */
static float sample_level(int l, float wx, float wz) {
    int st = 1 << l;
    float sp = st * CELL;
    int m = (N - 1) / st;
    float gx = wx / sp, gz = wz / sp;
    int i0 = (int)floorf(gx), j0 = (int)floorf(gz);
    float tx = gx - i0, tz = gz - j0;
    if (i0 < 0) { i0 = 0; tx = 0.0f; }
    if (j0 < 0) { j0 = 0; tz = 0.0f; }
    if (i0 >= m) { i0 = m; tx = 0.0f; }
    if (j0 >= m) { j0 = m; tz = 0.0f; }
    int i1 = i0 + 1 <= m ? i0 + 1 : m;
    int j1 = j0 + 1 <= m ? j0 + 1 : m;
    float a = Hf[i0 * st][j0 * st] * (1 - tx) + Hf[i1 * st][j0 * st] * tx;
    float b = Hf[i0 * st][j1 * st] * (1 - tx) + Hf[i1 * st][j1 * st] * tx;
    return a * (1 - tz) + b * tz;
}

static void gen_terrain(void) {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            float x = i * CELL, z = j * CELL;
            float d1x = x - 340.0f, d1z = z - 170.0f;
            float d2x = x - 130.0f, d2z = z - 360.0f;
            Hf[i][j] = 10.0f
                + 24.0f * fbm(x * 0.011f, z * 0.011f, 5, 2.0f, 0.5f, 7u)
                + 46.0f * expf(-(d1x * d1x + d1z * d1z) / (2.0f * 95.0f * 95.0f))
                + 30.0f * expf(-(d2x * d2x + d2z * d2z) / (2.0f * 70.0f * 70.0f));
        }
}

static void measure_error(void) {
    for (int l = 1; l <= LMAX; l++) {
        double mx = 0.0, se = 0.0;
        long cnt = 0;
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++) {
                float d = fabsf(sample_level(l, i * CELL, j * CELL) - Hf[i][j]);
                if (d > mx) mx = d;
                se += (double)d * d;
                cnt++;
            }
        eps_max[l] = (float)mx;
        eps_rms[l] = (float)sqrt(se / cnt);
    }
}

/* ---- LOD policy ---- */
static float switch_dist(int l, float tau) { return KPROJ * eps_max[l] / tau; }

static int select_level(float d, float tau) {
    int l = 0;
    while (l < LMAX && d > switch_dist(l + 1, tau)) l++;
    return l;
}

static int select_level_hyst(float d, float tau, float r, int cur) {
    while (cur < LMAX && d > switch_dist(cur + 1, tau) * sqrtf(r)) cur++;
    while (cur > 0 && d < switch_dist(cur, tau) / sqrtf(r)) cur--;
    return cur;
}

/* ---- renderer ---- */
static void norm3(float v[3]) {
    float s = sqrtf(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (s > 0) { v[0] /= s; v[1] /= s; v[2] /= s; }
}

static void sky_col(float dy, float c[3]) {
    float t = clampf(dy * 1.6f + 0.35f, 0.0f, 1.0f);
    c[0] = lerpf(0.78f, 0.45f, t);
    c[1] = lerpf(0.82f, 0.62f, t);
    c[2] = lerpf(0.88f, 0.92f, t);
}

static void render(int level, const float o[3], const float look[3], const char* path) {
    float fwd[3] = { look[0] - o[0], look[1] - o[1], look[2] - o[2] };
    norm3(fwd);
    float right[3] = { -fwd[2], 0.0f, fwd[0] };   /* cross(fwd, world_up) */
    norm3(right);
    float up[3] = {
        right[1] * fwd[2] - right[2] * fwd[1],
        right[2] * fwd[0] - right[0] * fwd[2],
        right[0] * fwd[1] - right[1] * fwd[0]
    };
    float sun[3] = { 0.5f, 0.75f, 0.35f };
    norm3(sun);
    unsigned char* img = malloc(IMG_W * IMG_H * 3);
    if (!img) { fprintf(stderr, "oom\n"); return; }

    for (int py = 0; py < IMG_H; py++)
        for (int px = 0; px < IMG_W; px++) {
            float sx = (px + 0.5f - IMG_W / 2.0f) / KPROJ;
            float sy = (IMG_H / 2.0f - py - 0.5f) / KPROJ;
            float dir[3] = {
                fwd[0] + sx * right[0] + sy * up[0],
                fwd[1] + sx * right[1] + sy * up[1],
                fwd[2] + sx * right[2] + sy * up[2]
            };
            norm3(dir);
            float t = 1.0f, tprev = 1.0f;
            int hit = 0;
            for (int s = 0; s < 4000 && t < 1600.0f; s++) {
                float qx = o[0] + dir[0] * t, qy = o[1] + dir[1] * t, qz = o[2] + dir[2] * t;
                float f = qy - sample_level(level, qx, qz);
                if (f < 0.0f) {
                    float a = tprev, b = t;
                    for (int k = 0; k < 10; k++) {
                        float m = 0.5f * (a + b);
                        float my = o[1] + dir[1] * m;
                        if (my - sample_level(level, o[0] + dir[0] * m, o[2] + dir[2] * m) < 0.0f) b = m;
                        else a = m;
                    }
                    t = 0.5f * (a + b);
                    hit = 1;
                    break;
                }
                tprev = t;
                float step = f * 0.35f;
                float minstep = t * 0.0008f + 0.15f;
                if (step < minstep) step = minstep;
                t += step;
            }
            float c[3];
            if (!hit) {
                sky_col(dir[1], c);
            } else {
                float qx = o[0] + dir[0] * t, qy = o[1] + dir[1] * t, qz = o[2] + dir[2] * t;
                float e = (1 << level) * CELL;
                float hx = sample_level(level, qx + e, qz) - sample_level(level, qx - e, qz);
                float hz = sample_level(level, qx, qz + e) - sample_level(level, qx, qz - e);
                float n[3] = { -hx / (2 * e), 1.0f, -hz / (2 * e) };
                norm3(n);
                float dif = n[0] * sun[0] + n[1] * sun[1] + n[2] * sun[2];
                if (dif < 0.0f) dif = 0.0f;
                float slope = 1.0f - n[1];
                float w_rock = clampf(slope * 3.0f, 0.0f, 1.0f);
                float w_snow = clampf((qy - 48.0f) / 18.0f, 0.0f, 1.0f);
                float grass[3] = { 0.32f, 0.44f, 0.22f };
                float rock[3] = { 0.45f, 0.40f, 0.36f };
                float snow[3] = { 0.85f, 0.87f, 0.90f };
                for (int k = 0; k < 3; k++) {
                    c[k] = lerpf(lerpf(grass[k], rock[k], w_rock), snow[k], w_snow * (1.0f - w_rock * 0.5f));
                    c[k] *= 0.35f + 0.65f * dif;
                }
                float fog = 1.0f - expf(-t / 2400.0f);
                float sk[3];
                sky_col(dir[1], sk);
                for (int k = 0; k < 3; k++) c[k] = lerpf(c[k], sk[k], fog);
            }
            unsigned char* p = img + (py * IMG_W + px) * 3;
            for (int k = 0; k < 3; k++)
                p[k] = (unsigned char)(powf(clampf(c[k], 0.0f, 1.0f), 1.0f / 2.2f) * 255.0f);
        }

    FILE* fp = fopen(path, "wb");
    if (!fp) { fprintf(stderr, "cannot open %s\n", path); free(img); return; }
    fprintf(fp, "P6\n%d %d\n255\n", IMG_W, IMG_H);
    fwrite(img, 1, IMG_W * IMG_H * 3, fp);
    fclose(fp);
    free(img);
    printf("rendered %s (level=%d)\n", path, level);
}

static void cam_from_anchor(float dist, float az, float el, const float anchor[3], float o[3]) {
    o[0] = anchor[0] + dist * cosf(el) * cosf(az);
    o[1] = anchor[1] + dist * sinf(el);
    o[2] = anchor[2] + dist * cosf(el) * sinf(az);
}

/* ---- camera path & switch accounting ---- */
static float jitter_of(int f) { return (rand01(hash2(f, 777, 5u)) - 0.5f) * 3.0f; }

static float path_dist(float t, int f) {
    float d;
    if (t < 20.0f) d = 700.0f - 640.0f * (t / 20.0f);              /* dolly in */
    else if (t < 25.0f) d = 60.0f + (D2_BASE - 60.0f) * ((t - 20.0f) / 5.0f);
    else if (t < 50.0f) d = D2_BASE + 3.0f * sinf(2.0f * PI * 1.5f * (t - 25.0f)); /* hover */
    else d = D2_BASE + (700.0f - D2_BASE) * ((t - 50.0f) / 10.0f);  /* dolly out */
    return d + jitter_of(f);
}

static int phase_of(float t) {
    if (t < 25.0f) return 0;
    if (t < 50.0f) return 1;
    return 2;
}

int main(int argc, char** argv) {
    const char* outdir = "artifacts/lod";
    int selftest = 0;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--selftest")) selftest = 1;
        else outdir = argv[i];
    }

    gen_terrain();
    measure_error();
    KPROJ = IMG_H / (2.0f * tanf(FOV_DEG * PI / 360.0f));
    D2_BASE = switch_dist(2, TAU_BASE);

    printf("K=%.1f px\n", KPROJ);
    for (int l = 1; l <= LMAX; l++) {
        float alpha = l > 1 ? logf(eps_max[l] / eps_max[l - 1]) / logf(2.0f) : 0.0f;
        printf("eps_L%d max=%.3f m rms=%.3f m ratio=%.2f alpha=%.2f d_base=%.0f m d_tuned=%.0f m\n",
               l, eps_max[l], eps_rms[l],
               l > 1 ? eps_max[l] / eps_max[l - 1] : 0.0f, alpha,
               switch_dist(l, TAU_BASE), switch_dist(l, TAU_TUNED));
    }
    for (int l = 1; l < LMAX; l++)
        printf("pop_L%d_to_L%d base=%.2f px tuned=%.2f px\n", l, l + 1,
               TAU_BASE * (1.0f - eps_max[l] / eps_max[l + 1]),
               TAU_TUNED * (1.0f - eps_max[l] / eps_max[l + 1]));

    /* switch accounting over the camera path */
    struct { const char* name; float tau, r; } cfgs[3] = {
        { "base(tau=4,r=1)", TAU_BASE, 1.0f },
        { "tight(tau=1.2,r=1)", TAU_TUNED, 1.0f },
        { "tuned(tau=1.2,r=1.35)", TAU_TUNED, R_TUNED },
    };
    int sw_total[3] = { 0, 0, 0 };
    int sw_phase[3][3] = { {0} };
    for (int c = 0; c < 3; c++) {
        int cur = 0;
        for (int f = 0; f < DUR * FPS; f++) {
            float t = (float)f / FPS;
            float d = path_dist(t, f);
            int l = cfgs[c].r > 1.0f
                ? select_level_hyst(d, cfgs[c].tau, cfgs[c].r, cur)
                : select_level(d, cfgs[c].tau);
            if (l != cur) {
                sw_total[c]++;
                sw_phase[c][phase_of(t)]++;
                cur = l;
            }
        }
        printf("switches %-26s total=%d dolly_in=%d hover=%d dolly_out=%d budget=%d %s\n",
               cfgs[c].name, sw_total[c], sw_phase[c][0], sw_phase[c][1], sw_phase[c][2],
               SWITCH_BUDGET, sw_total[c] <= SWITCH_BUDGET ? "OK" : "OVER");
    }

    /* renders */
    float anchor[3] = { 340.0f, 0.0f, 170.0f };
    anchor[1] = sample_level(0, anchor[0], anchor[2]);
    float d1t = switch_dist(1, TAU_TUNED), d2t = switch_dist(2, TAU_TUNED), d3t = switch_dist(3, TAU_TUNED);
    float d_near = 0.55f * d1t, d_mid = sqrtf(d1t * d2t), d_far = sqrtf(d2t * d3t);
    printf("tiers near=%.0f m(L%d) mid=%.0f m(L%d) far=%.0f m(L%d)\n",
           d_near, select_level(d_near, TAU_TUNED),
           d_mid, select_level(d_mid, TAU_TUNED),
           d_far, select_level(d_far, TAU_TUNED));

    if (!selftest) {
        char path[512];
        float o[3];
        float d_nc = 0.5f * (switch_dist(1, TAU_BASE) + switch_dist(2, TAU_BASE));

        cam_from_anchor(d_near, 0.9f, 0.14f, anchor, o);
        snprintf(path, sizeof path, "%s/tier_near.ppm", outdir);
        render(select_level(d_near, TAU_TUNED), o, anchor, path);

        cam_from_anchor(d_mid, 2.2f, 0.16f, anchor, o);
        snprintf(path, sizeof path, "%s/tier_mid.ppm", outdir);
        render(select_level(d_mid, TAU_TUNED), o, anchor, path);

        cam_from_anchor(d_far, 3.6f, 0.18f, anchor, o);
        snprintf(path, sizeof path, "%s/tier_far.ppm", outdir);
        render(select_level(d_far, TAU_TUNED), o, anchor, path);

        cam_from_anchor(d_nc, 0.9f, 0.05f, anchor, o);
        snprintf(path, sizeof path, "%s/near_base.ppm", outdir);
        render(select_level(d_nc, TAU_BASE), o, anchor, path);
        snprintf(path, sizeof path, "%s/near_tuned.ppm", outdir);
        render(select_level(d_nc, TAU_TUNED), o, anchor, path);

        /* far flicker pair: identical camera, level flips as the policy
         * re-evaluates at D2_BASE -/+ 4 m -> diff isolates the LOD pop */
        cam_from_anchor(D2_BASE, 3.6f, 0.18f, anchor, o);
        snprintf(path, sizeof path, "%s/far_b0.ppm", outdir);
        render(select_level(D2_BASE - 4.0f, TAU_BASE), o, anchor, path);
        snprintf(path, sizeof path, "%s/far_t0.ppm", outdir);
        render(select_level_hyst(D2_BASE - 4.0f, TAU_TUNED, R_TUNED, 0), o, anchor, path);
        snprintf(path, sizeof path, "%s/far_b1.ppm", outdir);
        render(select_level(D2_BASE + 4.0f, TAU_BASE), o, anchor, path);
        snprintf(path, sizeof path, "%s/far_t1.ppm", outdir);
        render(select_level_hyst(D2_BASE + 4.0f, TAU_TUNED, R_TUNED, 0), o, anchor, path);
    }

    /* self-checks */
    int fails = 0;
    int ok_eps = eps_max[1] > 0.0f && eps_max[1] < eps_max[2] && eps_max[2] < eps_max[3];
    printf("check eps_increasing=%d\n", ok_eps); fails += !ok_eps;
    int ok_alpha = eps_max[2] / eps_max[1] > 1.5f && eps_max[3] / eps_max[2] < 8.0f;
    printf("check eps_ratio_sane=%d\n", ok_alpha); fails += !ok_alpha;
    int ok_exact = fabsf(sample_level(0, 34.0f, 66.0f) - Hf[17][33]) < 1e-5f;
    printf("check level0_exact=%d\n", ok_exact); fails += !ok_exact;
    int ok_less = sw_total[2] < sw_total[0];
    printf("check_tuned_fewer_switches=%d\n", ok_less); fails += !ok_less;
    int ok_budget = sw_total[2] <= SWITCH_BUDGET;
    printf("check_tuned_within_budget=%d\n", ok_budget); fails += !ok_budget;
    printf("all_ok=%d\n", fails == 0);
    return fails == 0 ? 0 : 1;
}
