
/* quant: heightfield quantization -> height error -> normal error -> shadow fracture.
 *
 * Pipeline:
 *   1. build a deterministic fBm value-noise heightfield (same noise as ../noise);
 *   2. uniform mid-tread quantize at three step sizes derived from an angular budget;
 *   3. central-difference normals, measure normal angle error vs. the unquantized field;
 *   4. ray-marched directional (sun) shadow maps, count pixels whose shadow state flips;
 *   5. locate the terrain region where the error-amplification factor is maximal;
 *   6. refine (halve delta) until shadow fracture reaches zero, or report failure.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define N 64              /* heightfield grid side */
#define CELL 1.0f         /* world units per grid cell */
#define AMP 8.0f          /* height amplitude scale */
#define OCT 5
#define SEED 11u
#define NLEVELS 3         /* mandated three step sizes */
#define NLEV 4             /* sun elevation angles probed */
#define NB 8               /* region blocks per axis */
#define MAXREFINE 10

static float href[N][N];
static int   shref[NLEV][N][N];

/* ---- value noise / fBm, identical formulation to ../noise/src/noise.c ---- */
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

static float quantize(float v, float delta) { return roundf(v / delta) * delta; }

/* central-difference gradient (mirrors terrain field_slope, full 2-cell stencil) */
static void grad(const float h[N][N], int x, int y, float *gx, float *gy) {
    int xm = x > 0 ? x - 1 : x, xp = x < N - 1 ? x + 1 : x;
    int ym = y > 0 ? y - 1 : y, yp = y < N - 1 ? y + 1 : y;
    *gx = (h[xp][y] - h[xm][y]) / ((xp - xm) * CELL);
    *gy = (h[x][yp] - h[x][ym]) / ((yp - ym) * CELL);
}

/* directional self-shadow: march toward the sun; any horizon above the view ray
 * means the point is shadowed. sun azimuth fixed along +x. */
static int shadowed(const float h[N][N], int x, int y, float elev_deg) {
    float tan_e = tanf(elev_deg * (float)M_PI / 180.0f);
    for (int t = 1; x + t < N; t++) {
        if (h[x + t][y] - h[x][y] >= tan_e * (float)t * CELL) return 1;
    }
    return 0;
}

static int cmp_float(const void *a, const void *b) {
    float fa = *(const float *)a, fb = *(const float *)b;
    return (fa > fb) - (fa < fb);
}

/* solve delta for a target worst-case normal angle e_max at reference slope s:
 * tan(e) = q / (1 + s^2 - q^2/4), q = delta/d  (worst case, 1D slope) */
static float delta_for_budget(float budget_deg, float s) {
    float t = tanf(budget_deg * (float)M_PI / 180.0f);
    float a = t / 4.0f, b = 1.0f, c = -t * (1.0f + s * s);
    float q = (-b + sqrtf(b * b - 4.0f * a * c)) / (2.0f * a);
    return q * CELL;
}

typedef struct {
    float delta, rmse_h, max_h, rms_n, max_n;
    int   max_n_x, max_n_y;
    long  frac[NLEV], lit_to_sh[NLEV], sh_to_lit[NLEV];
} Report;

static Report run_level(float delta, const float elev[NLEV]) {
    Report r;
    memset(&r, 0, sizeof(r));
    r.delta = delta;
    float hq[N][N];
    for (int x = 0; x < N; x++)
        for (int y = 0; y < N; y++)
            hq[x][y] = quantize(href[x][y], delta);

    double se = 0.0, sae = 0.0;
    long ninterior = 0;
    for (int x = 1; x < N - 1; x++)
        for (int y = 1; y < N - 1; y++) {
            float e = hq[x][y] - href[x][y];
            se += (double)e * e;
            if (fabsf(e) > r.max_h) r.max_h = fabsf(e);

            float gx0, gy0, gx1, gy1;
            grad(href, x, y, &gx0, &gy0);
            grad(hq, x, y, &gx1, &gy1);
            float dot = gx0 * gx1 + gy0 * gy1 + 1.0f;
            float n0 = sqrtf(gx0 * gx0 + gy0 * gy0 + 1.0f);
            float n1 = sqrtf(gx1 * gx1 + gy1 * gy1 + 1.0f);
            float c = dot / (n0 * n1);
            if (c > 1.0f) c = 1.0f;
            if (c < -1.0f) c = -1.0f;
            float ang = acosf(c) * 180.0f / (float)M_PI;
            sae += (double)ang * ang;
            if (ang > r.max_n) { r.max_n = ang; r.max_n_x = x; r.max_n_y = y; }
            ninterior++;
        }
    r.rmse_h = (float)sqrt(se / (double)ninterior);
    r.rms_n = (float)sqrt(sae / (double)ninterior);

    for (int e = 0; e < NLEV; e++)
        for (int x = 0; x < N; x++)
            for (int y = 0; y < N; y++) {
                int qs = shadowed(hq, x, y, elev[e]);
                if (qs != shref[e][x][y]) {
                    r.frac[e]++;
                    if (qs) r.lit_to_sh[e]++; else r.sh_to_lit[e]++;
                }
            }
    return r;
}

int main(void) {
    /* 1. reference field */
    float hmin = 1e30f, hmax = -1e30f;
    for (int x = 0; x < N; x++)
        for (int y = 0; y < N; y++) {
            href[x][y] = AMP * fbm(x * 0.07f, y * 0.07f, OCT, 2.0f, 0.5f, SEED);
            if (href[x][y] < hmin) hmin = href[x][y];
            if (href[x][y] > hmax) hmax = href[x][y];
        }
    float range = hmax - hmin;

    float gmag[N][N] = {{0}};
    float flat[(N - 2) * (N - 2)];
    int nflat = 0;
    for (int x = 1; x < N - 1; x++)
        for (int y = 1; y < N - 1; y++) {
            float gx, gy;
            grad(href, x, y, &gx, &gy);
            gmag[x][y] = sqrtf(gx * gx + gy * gy);
            flat[nflat++] = gmag[x][y];
        }
    qsort(flat, nflat, sizeof(float), cmp_float);
    float g_low = flat[(int)(0.02f * nflat)];   /* 2nd-percentile slope: near-flat ref */
    float gmax = flat[nflat - 1];

    float elev[NLEV] = {10.0f, 20.0f, 30.0f, 45.0f};
    for (int e = 0; e < NLEV; e++)
        for (int x = 0; x < N; x++)
            for (int y = 0; y < N; y++)
                shref[e][x][y] = shadowed(href, x, y, elev[e]);

    printf("grid=%dx%d cell=%.2f amp=%.1f octaves=%d seed=%u\n", N, N, CELL, AMP, OCT, SEED);
    printf("height_range=[%.4f, %.4f] span=%.4f\n", hmin, hmax, range);
    printf("gradient: g_low(p2)=%.5f g_max=%.5f\n\n", g_low, gmax);

    /* 2. region map of amplification: mean slope per NB x NB block.
     * angle amplification factor 1/sqrt(1+s^2) is largest where s is smallest. */
    float bs[NB][NB] = {{0}};
    long  bc[NB][NB] = {{0}};
    for (int x = 1; x < N - 1; x++)
        for (int y = 1; y < N - 1; y++) {
            int bx = x * NB / N, by = y * NB / N;
            bs[bx][by] += gmag[x][y];
            bc[bx][by]++;
        }
    float bmin = 1e30f; int bbx = -1, bby = -1;
    for (int i = 0; i < NB; i++)
        for (int j = 0; j < NB; j++) {
            bs[i][j] /= (float)bc[i][j];
            if (bs[i][j] < bmin) { bmin = bs[i][j]; bbx = i; bby = j; }
        }
    /* single flattest cell too */
    float amin = 1e30f; int ax = -1, ay = -1;
    for (int x = 1; x < N - 1; x++)
        for (int y = 1; y < N - 1; y++)
            if (gmag[x][y] < amin) { amin = gmag[x][y]; ax = x; ay = y; }

    printf("max-amplification region (angle factor 1/sqrt(1+s^2)):\n");
    printf("  block[%d,%d] cells x:[%d,%d) y:[%d,%d) mean_slope=%.5f factor=%.6f\n",
           bbx, bby, bbx * N / NB, (bbx + 1) * N / NB,
           bby * N / NB, (bby + 1) * N / NB,
           bmin, 1.0f / sqrtf(1.0f + bmin * bmin));
    printf("  flattest single cell @(%d,%d) slope=%.5f factor=%.6f; "
           "linear gain |de/dh|=1/d=%.4f everywhere\n\n",
           ax, ay, amin, 1.0f / sqrtf(1.0f + amin * amin), 1.0f / CELL);

    /* 3. the three mandated step sizes from angular budgets 2 / 1 / 0.5 deg */
    float budget[NLEVELS] = {2.0f, 1.0f, 0.5f};
    Report rep[NLEVELS];
    for (int k = 0; k < NLEVELS; k++) {
        float d = delta_for_budget(budget[k], g_low);
        rep[k] = run_level(d, elev);
    }

    for (int k = 0; k < NLEVELS; k++) {
        Report *r = &rep[k];
        float q = r->delta / CELL;
        float e_flat = atanf(q / (1.0f + g_low * g_low - q * q / 4.0f)) * 180.0f / (float)M_PI;
        float e_steep = atanf(q / (1.0f + gmax * gmax - q * q / 4.0f)) * 180.0f / (float)M_PI;
        printf("== level %d: delta=%.5f (budget %.1f deg, codebook ~%.0f levels) ==\n",
               k, r->delta, budget[k], range / r->delta + 1.0f);
        printf("  height_error : rms=%.5f max_abs=%.5f (theory max delta/2=%.5f)\n",
               r->rmse_h, r->max_h, r->delta / 2.0f);
        printf("  normal_error : rms=%.4f deg max=%.4f deg @(%d,%d) slope=%.4f\n",
               r->rms_n, r->max_n, r->max_n_x, r->max_n_y, gmag[r->max_n_x][r->max_n_y]);
        printf("  theory bound : %.4f deg near-flat(s=%.4f), %.4f deg steepest(s=%.4f)\n",
               e_flat, g_low, e_steep, gmax);
        for (int e = 0; e < NLEV; e++)
            printf("  fracture@%2.0fdeg: %4ld px (lit->sh %ld, sh->lit %ld)\n",
                   (double)elev[e], r->frac[e], r->lit_to_sh[e], r->sh_to_lit[e]);
        long any = 0;
        for (int e = 0; e < NLEV; e++) any += r->frac[e];
        printf("  verdict: %s\n\n", any == 0 ? "complete (zero fracture)" : "INCOMPLETE (fracture nonzero)");
    }

    /* 4. refine by halving until all four shadow maps are fracture-free */
    printf("== refinement beyond the three mandated levels ==\n");
    float d = rep[NLEVELS - 1].delta / 2.0f;
    Report done;
    int got = 0;
    for (int i = 0; i < MAXREFINE; i++) {
        Report r = run_level(d, elev);
        long any = 0;
        for (int e = 0; e < NLEV; e++) any += r.frac[e];
        printf("  delta=%.6f rms_normal=%.4f deg fracture(sun10/20/30/45)=%ld/%ld/%ld/%ld %s\n",
               d, r.rms_n, r.frac[0], r.frac[1], r.frac[2], r.frac[3],
               any == 0 ? "<- zero" : "");
        if (any == 0) { done = r; got = 1; break; }
        d /= 2.0f;
    }
    if (got)
        printf("\nzero-fracture step: delta=%.6f  height rms=%.6f max=%.6f  "
               "normal rms=%.4f deg max=%.4f deg  codebook ~%.0f levels\n",
               done.delta, done.rmse_h, done.max_h, done.rms_n, done.max_n,
               range / done.delta + 1.0f);
    else
        printf("\nrefinement cap reached without zero fracture\n");
    return 0;
}
