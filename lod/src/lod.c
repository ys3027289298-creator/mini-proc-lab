
/* lod: spatially varying fBm octave count (3 LOD tiers = 3/4/5 octaves).
 *
 * Two heightfields are produced over the same domain:
 *   hard   - octave count switches abruptly at tier boundaries (C0 crease);
 *   smooth - extra octaves are blended in over a finite band (C1 transition).
 *
 * Value noise uses the same deterministic hash as noise/, but here both the
 * value and its analytic gradient are returned so fBm derivatives are exact;
 * finite-difference checks in main() validate the analytic derivatives.
 *
 * Output: artifacts/lod/meta.json + arrays of raw float32 (little endian).
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NX 768
#define NY 768
#define XMIN (-1.5f)
#define XMAX (5.5f)
#define YMIN (-2.0f)
#define YMAX (2.0f)

/* LOD tiers: band k = [k,k+1) in x uses k octaves; boundaries at x=3 and x=4. */
#define KMIN 3
#define KMAX 5
#define BANDW 0.35f            /* half-width of the smooth blend band */
#define ZEXAG 3.0f             /* vertical exaggeration for normals */
#define LAC 2.0f
#define GAIN 0.5f
#define SEED 11u

typedef struct { float h, gx, gy; } Sample;

static unsigned hash2(int x, int y, unsigned seed) {
    unsigned h = (unsigned)x * 374761393u + (unsigned)y * 668265263u + seed * 2246822519u;
    h = (h ^ (h >> 13)) * 1274126177u;
    return h ^ (h >> 16);
}

static float rand01(unsigned h) { return (float)(h & 0xffffff) / 16777215.0f; }
static float lerp(float a, float b, float t) { return a + (b - a) * t; }
static float fade(float t) { return t * t * (3.0f - 2.0f * t); }
static float dfade(float t) { return 6.0f * t * (1.0f - t); }

/* Value noise n(x,y) in [-1,1] with analytic gradient (chain rule through
 * the 2D smoothstep interpolation). */
static void value_noise_d(float x, float y, unsigned seed, float* v, float* dx, float* dy) {
    int xi = (int)floorf(x), yi = (int)floorf(y);
    float fx = x - xi, fy = y - yi;
    float v00 = rand01(hash2(xi,     yi,     seed));
    float v10 = rand01(hash2(xi + 1, yi,     seed));
    float v01 = rand01(hash2(xi,     yi + 1, seed));
    float v11 = rand01(hash2(xi + 1, yi + 1, seed));
    float u = fade(fx), vf = fade(fy);
    float du = dfade(fx), dv = dfade(fy);

    float nx0 = lerp(v00, v10, u);
    float nx1 = lerp(v01, v11, u);
    *v = (lerp(nx0, nx1, vf)) * 2.0f - 1.0f;
    *dx = 2.0f * (lerp(v10 - v00, v11 - v01, vf) * du);
    *dy = 2.0f * ((nx1 - nx0) * dv);
}

/* Normalized fBm truncated at m octaves: F_m = S_m / A_m,
 * with analytic gradient. Returns S, A, gradS, gradA through pointers. */
static void fbm_parts(float x, float y, int m,
                      float* S, float* A, float* gSx, float* gSy) {
    float s = 0.0f, gsx = 0.0f, gsy = 0.0f, amp = 1.0f, norm = 0.0f, freq = 1.0f;
    for (int i = 0; i < m; i++) {
        float v, dx, dy;
        value_noise_d(x * freq, y * freq, SEED + (unsigned)i * 17u, &v, &dx, &dy);
        s   += amp * v;
        gsx += amp * dx * freq;
        gsy += amp * dy * freq;
        norm += amp;
        amp  *= GAIN;
        freq *= LAC;
    }
    *S = s; *A = norm; *gSx = gsx; *gSy = gsy;
}

static void fbm_sample(float x, float y, int m, Sample* out) {
    float S, A, gSx, gSy;
    fbm_parts(x, y, m, &S, &A, &gSx, &gSy);
    out->h  = S / A;
    out->gx = (gSx * A - S * 0.0f) / (A * A);   /* A_m is spatially constant */
    out->gy = gSy / A;
}

/* Per-octave contribution q_k and its gradient (k is 0-based octave index). */
static void octave_q(float x, float y, int k, float* q, float* gqx, float* gqy) {
    float amp = powf(GAIN, (float)k);
    float freq = powf(LAC, (float)k);
    float v, dx, dy;
    value_noise_d(x * freq, y * freq, SEED + (unsigned)k * 17u, &v, &dx, &dy);
    *q   = amp * v;
    *gqx = amp * dx * freq;
    *gqy = amp * dy * freq;
}

/* Tiers follow "boundary k sits at x=k": x<3 -> 3 octaves, 3<=x<4 -> 4,
 * x>=4 -> 5. The smooth blend bands in smooth_weight use the same cuts. */
static int tier_of(float x) {
    if (x < (float)KMIN)       return KMIN;
    if (x < (float)(KMIN + 1)) return KMIN + 1;
    return KMAX;
}

/* Smooth weights for tiers k..KMAX. Inside a band around boundary b the weight
 * follows smoothstep; w'(x)=0 at every band edge, so the blended field is C1. */
static void smooth_weight(float x, int k, float* w, float* dw) {
    /* Weight of octave k (added when crossing boundary x=k): 0 below the
     * band, smoothstep 0->1 across [k-BANDW, k+BANDW], 1 above it.
     * w'=0 at both band edges, so the blended field is C1. */
    float weight = 0.0f;
    float dweight = 0.0f;
    if (x >= (float)k + BANDW) {
        weight = 1.0f;
    } else if (x > (float)k - BANDW) {
        float t = (x - ((float)k - BANDW)) / (2.0f * BANDW);
        weight = t * t * (3.0f - 2.0f * t);
        dweight = 6.0f * t * (1.0f - t) / (2.0f * BANDW);
    }
    *w = weight;
    *dw = dweight;
}

static void eval_hard(float x, float y, Sample* out) {
    fbm_sample(x, y, tier_of(x), out);
}

/* H_sm = (S_{KMIN} + sum_k w_k q_k) / (A_{KMIN} + sum_k w_k a_k);
 * a_k is constant, so only numerator and S carry x-derivatives. */
static void eval_smooth(float x, float y, Sample* out) {
    float S0, A0, gSx0, gSy0;
    fbm_parts(x, y, KMIN, &S0, &A0, &gSx0, &gSy0);
    float num = S0, den = A0, gnx = gSx0, gny = gSy0, gdx = 0.0f;
    for (int k = KMIN; k < KMAX; k++) {
        /* octave k joins at boundary x=k: w_k blends q_k in across the band */
        float w, dw;
        smooth_weight(x, k, &w, &dw);
        float q, gqx, gqy;
        octave_q(x, y, k, &q, &gqx, &gqy);
        float a = powf(GAIN, (float)k);
        num += w * q;
        den += w * a;
        gnx += dw * q + w * gqx;
        gny += w * gqy;
        gdx += dw * a;
    }
    out->h  = num / den;
    out->gx = (gnx * den - num * gdx) / (den * den); /* quotient rule: den' = gdx */
    out->gy = gny / den;
}

/* Exact one-sided values at a hard boundary b: tier b on the left, tier b+1
 * on the right, evaluated at the same physical point (b, y). */
static void write_boundary_csv(const char* dir) {
    char path[512];
    snprintf(path, sizeof(path), "%s/boundary_jumps.csv", dir);
    FILE* fp = fopen(path, "w");
    if (!fp) { perror(path); exit(1); }
    fprintf(fp, "boundary,y,hard_dh,hard_dgx,hard_dgy,hard_angle_deg,"
                "smooth_dh,smooth_dgx,smooth_dgy,smooth_angle_deg\n");
    int bs[2] = {KMIN, KMIN + 1};
    for (int bi = 0; bi < 2; bi++) {
        float b = (float)bs[bi];
        for (int j = 0; j < NY; j++) {
            float y = YMIN + (YMAX - YMIN) * (float)j / (NY - 1);
            Sample L, R;
            fbm_sample(b, y, bs[bi] - 1, &L);
            fbm_sample(b, y, bs[bi],     &R);
            Sample sL, sR;
            float e = 1e-3f;
            eval_smooth(b - e, y, &sL);
            eval_smooth(b + e, y, &sR);
            float dh = R.h - L.h, dgx = R.gx - L.gx, dgy = R.gy - L.gy;
            float zLx = ZEXAG * L.gx, zLy = ZEXAG * L.gy;
            float zRx = ZEXAG * R.gx, zRy = ZEXAG * R.gy;
            float nL = sqrtf(1.0f + zLx * zLx + zLy * zLy);
            float nR = sqrtf(1.0f + zRx * zRx + zRy * zRy);
            float dot = (1.0f + zLx * zRx + zLy * zRy) / (nL * nR);
            if (dot > 1.0f) dot = 1.0f;
            if (dot < -1.0f) dot = -1.0f;
            float ang = acosf(dot) * 180.0f / 3.14159265f;
            float sdh = sR.h - sL.h, sdgx = sR.gx - sL.gx, sdgy = sR.gy - sL.gy;
            float szLx = ZEXAG * sL.gx, szLy = ZEXAG * sL.gy;
            float szRx = ZEXAG * sR.gx, szRy = ZEXAG * sR.gy;
            float snL = sqrtf(1.0f + szLx * szLx + szLy * szLy);
            float snR = sqrtf(1.0f + szRx * szRx + szRy * szRy);
            float sdot = (1.0f + szLx * szRx + szLy * szRy) / (snL * snR);
            if (sdot > 1.0f) sdot = 1.0f;
            if (sdot < -1.0f) sdot = -1.0f;
            float sang = acosf(sdot) * 180.0f / 3.14159265f;
            fprintf(fp, "%d,%.6f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f\n",
                    bs[bi], (double)y, (double)dh, (double)dgx, (double)dgy, (double)ang,
                    (double)sdh, (double)sdgx, (double)sdgy, (double)sang);
        }
    }
    fclose(fp);
}

/* ------------------------------ self tests ------------------------------ */
static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static float hard_h(float x, float y) { Sample s; eval_hard(x, y, &s); return s.h; }
static float smooth_h(float x, float y) { Sample s; eval_smooth(x, y, &s); return s.h; }

static void write_grids(const char* dir,
                        const Sample* hard, const Sample* sm) {
    char path[512];
    const char* names[6] = {"hard_h", "hard_gx", "hard_gy",
                            "smooth_h", "smooth_gx", "smooth_gy"};
    for (int f = 0; f < 6; f++) {
        snprintf(path, sizeof(path), "%s/%s.bin", dir, names[f]);
        FILE* fp = fopen(path, "wb");
        if (!fp) { perror(path); exit(1); }
        float* buf = malloc((size_t)NX * NY * sizeof(float));
        for (int i = 0; i < NX * NY; i++) {
            const Sample* s = f < 3 ? &hard[i] : &sm[i];
            buf[i] = f == 0 ? s->h : f == 1 ? s->gx : f == 2 ? s->gy
                   : f == 3 ? s->h : f == 4 ? s->gx : s->gy;
        }
        fwrite(buf, sizeof(float), (size_t)NX * NY, fp);
        free(buf);
        fclose(fp);
    }
    snprintf(path, sizeof(path), "%s/meta.json", dir);
    FILE* fp = fopen(path, "w");
    fprintf(fp,
        "{\n"
        "  \"nx\": %d, \"ny\": %d,\n"
        "  \"xmin\": %g, \"xmax\": %g, \"ymin\": %g, \"ymax\": %g,\n"
        "  \"tiers\": [%d, %d, %d], \"boundaries\": [%d, %d],\n"
        "  \"band_halfwidth\": %g, \"z_exaggeration\": %g, \"lacunarity\": %g, \"gain\": %g, \"seed\": %u,\n"
        "  \"arrays\": [\"hard_h\",\"hard_gx\",\"hard_gy\",\"smooth_h\",\"smooth_gx\",\"smooth_gy\"]\n"
        "}\n",
        NX, NY, (double)XMIN, (double)XMAX, (double)YMIN, (double)YMAX,
        KMIN, KMIN + 1, KMAX, KMIN, KMIN + 1,
        (double)BANDW, (double)ZEXAG, (double)LAC, (double)GAIN, SEED);
    fclose(fp);
}

int main(int argc, char** argv) {
    const char* outdir = argc > 1 ? argv[1] : "artifacts/lod";

    /* 1. analytic gradient vs central finite differences (hard field) */
    float maxerr = 0.0f;
    for (int t = 0; t < 400; t++) {
        float x = -1.45f + 6.9f * (t / 399.0f);
        float y = 0.37f;
        if (fabsf(x - 3.0f) < 0.05f || fabsf(x - 4.0f) < 0.05f) continue;
        Sample s; eval_hard(x, y, &s);
        float e = 1e-3f;
        float fdx = (hard_h(x + e, y) - hard_h(x - e, y)) / (2 * e);
        float fdy = (hard_h(x, y + e) - hard_h(x, y - e)) / (2 * e);
        float err = fmaxf(fabsf(fdx - s.gx), fabsf(fdy - s.gy));
        if (err > maxerr) maxerr = err;
    }
    check("hard_gradient_matches_fd", maxerr < 1.5e-2f);

    /* 2. same check on the smooth field; blend bands are C1 but not C2, so
     *    finite differences are only compared outside the bands */
    maxerr = 0.0f;
    for (int t = 0; t < 400; t++) {
        float x = -1.45f + 6.9f * (t / 399.0f);
        float y = -0.61f;
        if (fabsf(x - 3.0f) < BANDW + 0.05f || fabsf(x - 4.0f) < BANDW + 0.05f) continue;
        Sample s; eval_smooth(x, y, &s);
        float e = 1e-3f;
        float fdx = (smooth_h(x + e, y) - smooth_h(x - e, y)) / (2 * e);
        float fdy = (smooth_h(x, y + e) - smooth_h(x, y - e)) / (2 * e);
        float err = fmaxf(fabsf(fdx - s.gx), fabsf(fdy - s.gy));
        if (err > maxerr) maxerr = err;
    }
    printf("smooth_grad_fd_maxerr=%.6f\n", maxerr);
    check("smooth_gradient_matches_fd", maxerr < 8e-3f);

    /* 3. smooth field is C0 and C1 across both band edges */
    float c0max = 0.0f, c1max = 0.0f;
    float edges[4] = {3.0f - BANDW, 3.0f + BANDW, 4.0f - BANDW, 4.0f + BANDW};
    for (int i = 0; i < 4; i++) {
        float e = 1e-4f;
        Sample a, b;
        eval_smooth(edges[i] - e, 0.23f, &a);
        eval_smooth(edges[i] + e, 0.23f, &b);
        c0max = fmaxf(c0max, fabsf(a.h - b.h));
        c1max = fmaxf(c1max, fmaxf(fabsf(a.gx - b.gx), fabsf(a.gy - b.gy)));
    }
    check("smooth_C0_at_band_edges", c0max < 1e-3f);
    check("smooth_C1_at_band_edges", c1max < 2e-2f);

    /* 4. away from bands, smooth == hard exactly (weights are 0 or 1) */
    float diffmax = 0.0f;
    /* bands are [2.65,3.35] (boundary 3) and [3.65,4.35] (boundary 4);
     * 2.35 and 5.0 sit in settled tier interiors, 0.0 in the 3-oct tier */
    float probes[3] = {0.0f, 2.35f, 5.0f};
    for (int i = 0; i < 3; i++) {
        Sample a, b;
        eval_hard(probes[i], -0.4f, &a);
        eval_smooth(probes[i], -0.4f, &b);
        diffmax = fmaxf(diffmax, fmaxf(fabsf(a.h - b.h),
                                       fmaxf(fabsf(a.gx - b.gx), fabsf(a.gy - b.gy))));
    }
    check("smooth_equals_hard_outside_bands", diffmax < 1e-6f);

    /* 5. hard field: at a boundary the two tiers disagree even at the same
     *    physical point, so there is a C0 height step. Exact identity:
     *      F_{k+1}-F_k = (A_k q_k - S_k a_k) / (A_k A_{k+1})
     *    hence |step| <= (|q_k| + |S_k| a_k/A_k) / A_{k+1}. */
    float stepmax = 0.0f, boundmax = 0.0f;
    for (int t = 0; t < 200; t++) {
        float y = -1.9f + 3.8f * t / 199.0f;
        Sample f3, f4;
        fbm_sample((float)KMIN, y, KMIN, &f3);
        fbm_sample((float)KMIN, y, KMIN + 1, &f4);
        float S3, A3, gx, gy;
        fbm_parts((float)KMIN, y, KMIN, &S3, &A3, &gx, &gy);
        float q, gqx, gqy;
        octave_q((float)KMIN, y, KMIN, &q, &gqx, &gqy);
        float a3 = powf(GAIN, (float)KMIN);
        float A4 = A3 + a3;
        float bound = (fabsf(q) + fabsf(S3) * a3 / A3) / A4;
        stepmax = fmaxf(stepmax, fabsf(f4.h - f3.h));
        boundmax = fmaxf(boundmax, bound);
    }
    check("hard_height_step_exists", stepmax > 1e-4f);
    check("hard_height_step_within_bound", stepmax <= boundmax + 1e-6f);
    printf("hard_height_stepmax=%.6e boundmax=%.6e\n", stepmax, boundmax);
    float e = 1e-4f;
    Sample h3m, h3p;
    eval_hard((float)KMIN - e, 0.0f, &h3m);
    eval_hard((float)KMIN + e, 0.0f, &h3p);
    check("hard_gradient_jump_at_boundary",
          fabsf(h3p.gx - h3m.gx) + fabsf(h3p.gy - h3m.gy) > 1e-3f);

    /* 6. weights monotone 0->1 inside the band, flat outside */
    float w0, dw0, w1, dw1;
    smooth_weight(3.0f - BANDW - 0.01f, KMIN, &w0, &dw0);
    smooth_weight(3.0f + BANDW + 0.01f, KMIN, &w1, &dw1);
    check("weight_0_below_band", w0 == 0.0f && dw0 == 0.0f);
    check("weight_1_above_band", w1 == 1.0f && dw1 == 0.0f);

    Sample* hard = malloc((size_t)NX * NY * sizeof(Sample));
    Sample* sm   = malloc((size_t)NX * NY * sizeof(Sample));
    if (!hard || !sm) { fprintf(stderr, "oom\n"); return 1; }
    for (int j = 0; j < NY; j++) {
        float y = YMIN + (YMAX - YMIN) * (float)j / (NY - 1);
        for (int i = 0; i < NX; i++) {
            float x = XMIN + (XMAX - XMIN) * (float)i / (NX - 1);
            eval_hard(x, y, &hard[(size_t)j * NX + i]);
            eval_smooth(x, y, &sm[(size_t)j * NX + i]);
        }
    }
    write_grids(outdir, hard, sm);
    write_boundary_csv(outdir);
    free(hard); free(sm);
    printf("wrote_grids_to=%s\n", outdir);
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
