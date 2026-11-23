
/* blend: terrain material blend-kernel comparison on one shared heightfield.
 *
 * Three weight kernels map the signed height distance d = h - h0 over a band
 * half-width R to a material-A weight w in [0,1]:
 *   K_LINEAR  linear ramp            C0 only, gradient jumps at band edges
 *   K_QUINTIC quintic smootherstep   C2, polynomial, no transcendental ops
 *   K_GAUSS   erf ramp, truncated    smooth but needs erff; truncation seam
 *
 * Same field, same boundary h0, same R for every kernel. The program emits
 * self-test lines plus CSV traces consumed by tools/plot.py. */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define N 64            /* heightfield grid: N x N */
#define CELL 1.0f
#define R_BAND 4.0f     /* blend band half-width R, in height units */
#define H0 0.0f         /* material boundary height */
#define NSAMP 401       /* fine 1-D samples across the band */
#define BENCH_ITERS 400

typedef enum { K_LINEAR = 0, K_QUINTIC = 1, K_GAUSS = 2 } Kernel;

static const char* kname(Kernel k) {
    return k == K_LINEAR ? "linear" : k == K_QUINTIC ? "quintic" : "gauss";
}

/* ---------- deterministic value noise / fBm (mirrors noise/) ---------- */
static unsigned hash2(int x, int y, unsigned seed) {
    unsigned h = (unsigned)x * 374761393u + (unsigned)y * 668265263u
               + seed * 2246822519u;
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

static float fbm(float x, float y, int octaves, float lac, float gain,
                 unsigned seed) {
    float sum = 0.0f, amp = 1.0f, norm = 0.0f, freq = 1.0f;
    for (int i = 0; i < octaves; i++) {
        sum += value_noise(x * freq, y * freq, seed + (unsigned)i * 17u) * amp;
        norm += amp;
        amp *= gain;
        freq *= lac;
    }
    return norm > 0.0f ? sum / norm : 0.0f;
}

/* ---------- shared terrain ---------- */
typedef struct { float h[N][N]; } Field;

static void terrain_build(Field* f) {
    /* fBm broad relief plus a narrow steep ridge so the h0 boundary
     * crosses a genuine cliff (per-cell height steps up to ~5). */
    for (int y = 0; y < N; y++)
        for (int x = 0; x < N; x++) {
            float wx = x * 0.08f, wy = y * 0.08f;
            float broad = fbm(wx, wy, 4, 2.0f, 0.5f, 11u);
            float ridge = 16.0f * expf(-((x - 20.0f) * (x - 20.0f)) / 4.0f);
            f->h[x][y] = 6.0f * broad + ridge - 4.0f;
        }
}

static float field_sample(const Field* f, float gx, float gy) {
    int x0 = (int)floorf(gx), y0 = (int)floorf(gy);
    float tx = gx - x0, ty = gy - y0;
    int x1 = x0 + 1, y1 = y0 + 1;
    if (x0 < 0) x0 = 0;
    if (y0 < 0) y0 = 0;
    if (x1 >= N) x1 = N - 1;
    if (y1 >= N) y1 = N - 1;
    float a = f->h[x0][y0] * (1 - tx) + f->h[x1][y0] * tx;
    float b = f->h[x0][y1] * (1 - tx) + f->h[x1][y1] * tx;
    return a * (1 - ty) + b * ty;
}

/* central-difference height gradient magnitude at grid vertex (x,y) */
static float height_slope(const Field* f, int x, int y) {
    int xm = x > 0 ? x - 1 : x, xp = x < N - 1 ? x + 1 : x;
    int ym = y > 0 ? y - 1 : y, yp = y < N - 1 ? y + 1 : y;
    float dx = (f->h[xp][y] - f->h[xm][y]) / ((xp - xm) * CELL);
    float dy = (f->h[x][yp] - f->h[x][ym]) / ((yp - ym) * CELL);
    return sqrtf(dx * dx + dy * dy);
}

/* ---------- kernels: weight and analytic derivative dK/dd ---------- */
static float clamp01(float t) { return t < 0 ? 0 : (t > 1 ? 1 : t); }

static float kw_linear(float d, float R) {
    return clamp01(0.5f + 0.5f * d / R);
}
static float kd_linear(float d, float R) {
    return (d > -R && d < R) ? 0.5f / R : 0.0f;
}

static float kw_quintic(float d, float R) {
    float t = clamp01(0.5f + 0.5f * d / R);
    return t * t * t * (t * (t * 6.0f - 15.0f) + 10.0f);
}
/* dK/dd = (6 t^5 - 15 t^4 + 10 t^3)' wrt t, times dt/dd = 1/(2R) */
static float kd_quintic(float d, float R) {
    if (d <= -R || d >= R) return 0.0f;
    float t = 0.5f + 0.5f * d / R;
    float dkdt = 30.0f * t * t * (1.0f - t) * (1.0f - t);
    return dkdt / (2.0f * R);
}

/* Gaussian CDF ramp: w = 0.5*(1+erf(d/(sqrt2 sigma))), sigma = R/2,
 * windowed to |d| <= R (= +-2 sigma). */
#define SQRT2 1.4142135623730951f
static float kw_gauss(float d, float R) {
    if (d <= -R) return 0.0f;
    if (d >= R) return 1.0f;
    float sigma = R * 0.5f;
    return 0.5f * (1.0f + erff(d / (SQRT2 * sigma)));
}
static float kd_gauss(float d, float R) {
    if (d <= -R || d >= R) return 0.0f;
    float sigma = R * 0.5f;
    return expf(-0.5f * (d / sigma) * (d / sigma)) / (sigma * 2.5066282746310002f);
}

static float kw(Kernel k, float d, float R) {
    if (k == K_LINEAR) return kw_linear(d, R);
    if (k == K_QUINTIC) return kw_quintic(d, R);
    return kw_gauss(d, R);
}
static float kd(Kernel k, float d, float R) {
    if (k == K_LINEAR) return kd_linear(d, R);
    if (k == K_QUINTIC) return kd_quintic(d, R);
    return kd_gauss(d, R);
}

/* ---------- measurements ---------- */
typedef struct {
    float jump[2];   /* |dK/dd| just inside band edge: index0=-R, index1=+R */
    float second;    /* max |d2K/dd2| over the band (curvature proxy) */
    float rms_sec;   /* RMS second derivative */
} Smoothness;

static Smoothness smoothness_measure(Kernel k, float R) {
    Smoothness s;
    float eps = R * 1e-3f;
    s.jump[0] = kd(k, -R + eps, R) - kd(k, -R - eps, R);
    s.jump[1] = kd(k, R - eps, R) - kd(k, R + eps, R);
    s.jump[0] = fabsf(s.jump[0]);
    s.jump[1] = fabsf(s.jump[1]);
    float hi = 0.0f, acc = 0.0f;
    int n = 0;
    float h = (2.0f * R) / 1000.0f;
    for (int i = 0; i <= 1000; i++) {
        float d = -R + (2.0f * R) * i / 1000.0f;
        float d2 = (kd(k, d + h, R) - kd(k, d - h, R)) / (2.0f * h);
        float a = fabsf(d2);
        if (a > hi) hi = a;
        acc += d2 * d2;
        n++;
    }
    s.second = hi;
    s.rms_sec = sqrtf(acc / (float)n);
    return s;
}

/* one full-field blending pass, sink the result so it cannot be optimized out */
static float full_pass(Kernel k, const Field* f, float R) {
    float sink = 0.0f;
    for (int y = 0; y < N; y++)
        for (int x = 0; x < N; x++)
            sink += kw(k, f->h[x][y] - H0, R);
    return sink;
}

static double bench_ns(Kernel k, const Field* f, float R) {
    struct timespec t0, t1;
    float sink = 0.0f;
    double best = 1e300;
    for (int rep = 0; rep < 5; rep++) {
        clock_gettime(CLOCK_MONOTONIC, &t0);
        for (int i = 0; i < BENCH_ITERS; i++) sink += full_pass(k, f, R);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        double ns = ((t1.tv_sec - t0.tv_sec) * 1e9
                   + (t1.tv_nsec - t0.tv_nsec))
                  / (double)(BENCH_ITERS * N * N);
        if (ns < best) best = ns;
    }
    if (sink == -12345.0f) printf("impossible\n"); /* keep sink live */
    return best;
}

/* isolated kernel cost: arithmetic only, inputs from a small rotating
 * buffer so the field's memory traffic cannot mask the kernel work */
static double bench_kernel_only_ns(Kernel k, float R) {
    enum { M = 256 };
    float ds[M];
    for (int i = 0; i < M; i++) ds[i] = -R * 1.3f + (2.6f * R) * i / M;
    struct timespec t0, t1;
    float sink = 0.0f;
    long iters = 200000;
    double best = 1e300;
    for (int rep = 0; rep < 7; rep++) {
        clock_gettime(CLOCK_MONOTONIC, &t0);
        for (long i = 0; i < iters; i++) sink += kw(k, ds[i & (M - 1)], R);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        double ns = ((t1.tv_sec - t0.tv_sec) * 1e9
                   + (t1.tv_nsec - t0.tv_nsec)) / (double)iters;
        if (ns < best) best = ns;
    }
    if (sink == -12345.0f) printf("impossible\n");
    return best;
}

/* ---------- steep-slope stair detection ---------- */
typedef struct { int x, y; float slope, dstep, wdelta, slope_n; } Stair;

/* scan the cliff line (fixed y), report per-cell weight steps above tol */
static int stair_scan(Kernel k, const Field* f, float R, int yrow,
                      float tol, Stair* out, int maxout) {
    int n = 0;
    float wp = kw(k, f->h[0][yrow] - H0, R);
    for (int x = 1; x < N && n < maxout; x++) {
        float w = kw(k, f->h[x][yrow] - H0, R);
        float dw = fabsf(w - wp);
        float sl = height_slope(f, x, yrow);
        if (dw > tol && sl > 1.5f && n < maxout) {
            out[n].x = x; out[n].y = yrow;
            out[n].slope = sl;
            out[n].dstep = fabsf(f->h[x][yrow] - f->h[x - 1][yrow]);
            out[n].wdelta = dw;
            out[n].slope_n = sl;
            n++;
        }
        wp = w;
    }
    return n;
}

/* slope-adaptive band half-width driven by the actual per-cell height
 * step |dh| (the true cause of a weight step), not the smoothed slope.
 * Conservative bound dw <= peak*|dh|/R  =>  R >= peak*|dh|/tol.
 * peak = 0.5 (linear) or 0.9375 (quintic max dK/dd in t units).
 * Clamped to [R_BAND, 8 R_BAND]. */
static float adaptive_R(Kernel k, float hstep, float tol) {
    float peak = (k == K_LINEAR) ? 0.5f : 0.9375f;
    float R = peak * hstep / tol;
    float cap = R_BAND * 8.0f;
    return R < R_BAND ? R_BAND : (R > cap ? cap : R);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Field f;
    terrain_build(&f);

    /* ---- self-tests ---- */
    check("weight_endpoints_linear",
          kw_linear(-R_BAND, R_BAND) == 0.0f &&
          kw_linear(R_BAND, R_BAND) == 1.0f);
    check("weight_endpoints_quintic",
          kw_quintic(-R_BAND, R_BAND) == 0.0f &&
          kw_quintic(R_BAND, R_BAND) == 1.0f);
    check("weight_endpoints_gauss",
          fabsf(kw_gauss(-R_BAND, R_BAND)) < 1e-3f &&
          fabsf(kw_gauss(R_BAND, R_BAND) - 1.0f) < 1e-3f);
    check("weight_half_at_boundary",
          fabsf(kw(K_QUINTIC, 0.0f, R_BAND) - 0.5f) < 1e-6f &&
          fabsf(kw(K_GAUSS, 0.0f, R_BAND) - 0.5f) < 1e-6f &&
          fabsf(kw(K_LINEAR, 0.0f, R_BAND) - 0.5f) < 1e-6f);
    check("quintic_c2_zero_edge_deriv",
          fabsf(kd_quintic(-R_BAND + 1e-4f, R_BAND)) < 1e-4f &&
          fabsf(kd_quintic(R_BAND - 1e-4f, R_BAND)) < 1e-4f);
    check("gauss_symmetric",
          fabsf(kw_gauss(1.3f, R_BAND) + kw_gauss(-1.3f, R_BAND) - 1.0f)
          < 1e-6f);
    check("gauss_truncation_seam_present",
          kd_gauss(-R_BAND - 1e-3f, R_BAND) == 0.0f &&
          kd_gauss(-R_BAND + 1e-3f, R_BAND) > 0.0f);

    /* ---- smoothness + timing per kernel; CSV row table ---- */
    /* prefer repo-root artifacts/ when run from the module dir */
    FILE* fprobe = fopen("../artifacts/blend/.probe", "w");
    const char* dir;
    if (fprobe) { fclose(fprobe); remove("../artifacts/blend/.probe"); dir = "../artifacts/blend"; }
    else {
        if (fopen("artifacts/blend/.probe", "w")) { remove("artifacts/blend/.probe"); dir = "artifacts/blend"; }
        else dir = ".";
    }
    char path[256];
#define OPEN_CSV(var, name) \
    snprintf(path, sizeof(path), "%s/" name, dir); \
    var = fopen(path, "w")
    FILE* fc;
    FILE* ft;
    OPEN_CSV(fc, "kernel_metrics.csv");
    OPEN_CSV(ft, "kernel_1d.csv");
    if (!fc || !ft) { perror("csv"); return 2; }
    fprintf(fc, "kernel,full_ns_per_eval,kernel_only_ns,jump_minusR,"
                "jump_plusR,max_abs_second_deriv,rms_second_deriv,op_note\n");
    fprintf(ft, "d_over_R,kernel,weight,grad\n");
    double cost[3], kcost[3];
    for (int ki = 0; ki < 3; ki++) {
        Kernel k = (Kernel)ki;
        Smoothness s = smoothness_measure(k, R_BAND);
        cost[ki] = bench_ns(k, &f, R_BAND);
        kcost[ki] = bench_kernel_only_ns(k, R_BAND);
        const char* note = k == K_LINEAR ? "1 div, no transcendental"
                         : k == K_QUINTIC ? "polynomial, no transcendental"
                         : "erff+expf per eval";
        fprintf(fc, "%s,%.4f,%.4f,%.6f,%.6f,%.6f,%.6f,%s\n", kname(k),
                cost[ki], kcost[ki], s.jump[0], s.jump[1],
                s.second, s.rms_sec, note);
        for (int i = 0; i < NSAMP; i++) {
            float u = -1.2f + 2.4f * i / (NSAMP - 1); /* d/R */
            float d = u * R_BAND;
            fprintf(ft, "%.4f,%s,%.7f,%.7f\n", u, kname(k),
                    kw(k, d, R_BAND), kd(k, d, R_BAND));
        }
    }
    fclose(fc); fclose(ft);

    /* ---- on-terrain boundary gradient trace along row y=32 ---- */
    OPEN_CSV(ft, "boundary_trace.csv");
    fprintf(ft, "x,height,slope,w_linear,w_quintic,w_gauss,"
                "gw_linear,gw_quintic,gw_gauss\n");
    int yrow = 32;
    for (int x = 0; x < N; x++) {
        float h = f.h[x][yrow], d = h - H0;
        float sl = height_slope(&f, x, yrow);
        /* on-terrain gradient of w along x via central difference in x */
        float wl[3], gw[3];
        for (int ki = 0; ki < 3; ki++) {
            Kernel k = (Kernel)ki;
            wl[ki] = kw(k, d, R_BAND);
            float hm = field_sample(&f, (float)x - 1.0f, (float)yrow);
            float hp = field_sample(&f, (float)x + 1.0f, (float)yrow);
            gw[ki] = (kw(k, hp - H0, R_BAND) - kw(k, hm - H0, R_BAND))
                   / (2.0f * CELL);
        }
        fprintf(ft, "%d,%.5f,%.5f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                x, h, sl, wl[0], wl[1], wl[2], gw[0], gw[1], gw[2]);
    }
    fclose(ft);

    /* ---- steep-slope stair report (fixed R) ---- */
    Stair stairs[64];
    int nst = stair_scan(K_QUINTIC, &f, R_BAND, yrow, 0.20f, stairs, 64);
    OPEN_CSV(ft, "stairs.csv");
    fprintf(ft, "x,y,slope,height_step,weight_delta_fixedR,"
                "weight_delta_adaptiveR,adaptive_R\n");
    int adapted_violations = 0;
    for (int i = 0; i < nst; i++) {
        int x = stairs[i].x;
        float Ra = adaptive_R(K_QUINTIC, stairs[i].dstep, 0.20f);
        float wa = kw(K_QUINTIC, f.h[x][yrow] - H0, Ra);
        float wb = kw(K_QUINTIC, f.h[x - 1][yrow] - H0, Ra);
        float dwa = fabsf(wa - wb);
        if (dwa > 0.20f) adapted_violations++;
        fprintf(ft, "%d,%d,%.4f,%.4f,%.6f,%.6f,%.4f\n",
                x, yrow, stairs[i].slope, stairs[i].dstep,
                stairs[i].wdelta, dwa, Ra);
    }
    fclose(ft);
    printf("stair_cells_fixedR=%d\n", nst);
    printf("stair_cells_adaptiveR_remaining=%d\n", adapted_violations);
    check("adaptive_bandwidth_removes_steps", adapted_violations == 0);

    /* steepest vertex of the terrain, for the report */
    float smax = -1.0f; int sx = 0, sy = 0;
    for (int y = 1; y < N - 1; y++)
        for (int x = 1; x < N - 1; x++) {
            float s = height_slope(&f, x, y);
            if (s > smax) { smax = s; sx = x; sy = y; }
        }
    printf("steepest_vertex=(%d,%d) slope=%.3f height=%.3f\n",
           sx, sy, smax, f.h[sx][sy]);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
