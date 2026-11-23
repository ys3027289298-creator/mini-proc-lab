
/* collision: render-sample vs physics-sample mismatch on one heightfield.
 *
 * Two sampling paths are compared:
 *   render  : bilinear sample of the fine visual mesh (grid spacing c)
 *             at the true world position p;
 *   physics : body position p snapped to a collision lattice (spacing s),
 *             then bilinear sample of a coarse collision mesh whose
 *             vertices are a decimation (stride k, spacing d = k*c) of
 *             the same field.
 *
 * Two error upper bounds are derived and checked:
 *   B_a  a-priori analytic bound from the noise generator's Lipschitz
 *        constant (holds before the field is built);
 *   B_m  bound on the instantiated field using its exact max bilinear
 *        gradient norm (tighter, field-specific).
 * Then both paths are unified onto the shared fine sampler at the same
 * world point and penetration must be zero. */
#include <stdio.h>
#include <math.h>
#include <string.h>

/* ---------- deterministic value noise / fBm (same family as noise/) ------ */
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

#define OCTAVES 4
#define LACUNARITY 2.0f
#define GAIN 0.5f
#define NOISE_SEED 11u
#define HEIGHT_AMP 6.0f

static float fbm_height(float x, float y) {
    float sum = 0.0f, amp = 1.0f, norm = 0.0f, freq = 1.0f;
    for (int i = 0; i < OCTAVES; i++) {
        sum += value_noise(x * freq, y * freq, NOISE_SEED + (unsigned)i * 17u) * amp;
        norm += amp;
        amp *= GAIN;
        freq *= LACUNARITY;
    }
    return HEIGHT_AMP * (sum / norm);
}

/* ------------------------------ heightfields ----------------------------- */
#define N 32            /* fine grid is N x N vertices */
#define CELL 0.5f       /* fine grid spacing c (metres) */
#define STRIDE 4        /* coarse collision mesh takes every STRIDE vertex */
#define SNAP 1.0f       /* physics lattice spacing s (metres) */
#define EPS 1e-5f       /* penetration tolerance for float comparison */

typedef struct { int n; float cell; float h[N][N]; } FineField;
typedef struct {
    int n;
    float cell;
    float h[N / STRIDE + 1][N / STRIDE + 1];
} CoarseField;

static void fine_build(FineField* f) {
    f->n = N;
    f->cell = CELL;
    for (int x = 0; x < N; x++)
        for (int y = 0; y < N; y++)
            f->h[x][y] = fbm_height(x * CELL, y * CELL);
}

static void coarse_build(CoarseField* g, const FineField* f) {
    g->n = (N - 1) / STRIDE + 1;
    g->cell = CELL * STRIDE;                 /* d = k*c */
    for (int i = 0; i < g->n; i++)
        for (int j = 0; j < g->n; j++)
            g->h[i][j] = f->h[i * STRIDE][j * STRIDE];
}

static float fine_sample(const FineField* f, float wx, float wy) {
    float gx = wx / f->cell, gy = wy / f->cell;
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

static float coarse_sample(const CoarseField* g, float wx, float wy) {
    int m = g->n;
    float gx = wx / g->cell, gy = wy / g->cell;
    int x0 = (int)floorf(gx), y0 = (int)floorf(gy);
    float tx = gx - x0, ty = gy - y0;
    int x1 = x0 + 1, y1 = y0 + 1;
    if (x0 < 0) x0 = 0;
    if (y0 < 0) y0 = 0;
    if (x1 >= m) x1 = m - 1;
    if (y1 >= m) y1 = m - 1;
    float a = g->h[x0][y0] * (1 - tx) + g->h[x1][y0] * tx;
    float b = g->h[x0][y1] * (1 - tx) + g->h[x1][y1] * tx;
    return a * (1 - ty) + b * ty;
}

/* physics step: body position is snapped to the collision lattice first. */
static void physics_snap(float wx, float wy, float* sx, float* sy) {
    *sx = roundf(wx / SNAP) * SNAP;
    *sy = roundf(wy / SNAP) * SNAP;
}

/* ------------------------- Lipschitz / slope bounds ---------------------- */
/* Analytic a-priori Lipschitz constant of the generator F.
 * One value-noise octave n (lattice frequency f_n, amplitude a_n):
 *   smoothstep derivative max = 1.5, remapped corner values in [-1,1]
 *   => per-axis |dV/dx| <= 2 * 1.5 * f_n = 3 f_n,
 *   => gradient norm  <= 3*sqrt(2) * f_n.
 * fBm normalizes by norm = sum a_n and is scaled by HEIGHT_AMP:
 *   L_a = HEIGHT_AMP * 3*sqrt(2) * (sum a_n f_n) / (sum a_n). */
static float analytic_lipschitz(void) {
    float sum_af = 0.0f, amp = 1.0f, freq = 1.0f, norm = 0.0f;
    for (int i = 0; i < OCTAVES; i++) {
        sum_af += amp * freq;
        norm += amp;
        amp *= GAIN;
        freq *= LACUNARITY;
    }
    return HEIGHT_AMP * 3.0f * sqrtf(2.0f) * sum_af / norm;
}

/* Exact Euclidean Lipschitz constant of the instantiated fine bilinear
 * surface: max gradient norm over every cell. On a cell each partial is
 * affine in the other parameter, so the norm (sqrt of a convex quadratic)
 * attains its maximum at one of the four parameter corners. */
static float measured_gradient_lipschitz(const FineField* f) {
    float Lg = 0.0f;
    for (int x = 0; x + 1 < N; x++) {
        for (int y = 0; y + 1 < N; y++) {
            float sx0 = (f->h[x + 1][y]     - f->h[x][y])     / CELL;
            float sx1 = (f->h[x + 1][y + 1] - f->h[x][y + 1]) / CELL;
            float sy0 = (f->h[x][y + 1]     - f->h[x][y])     / CELL;
            float sy1 = (f->h[x + 1][y + 1] - f->h[x + 1][y]) / CELL;
            float g[4][2] = {
                { sx0, sy0 }, { sx1, sy0 }, { sx0, sy1 }, { sx1, sy1 }
            };
            for (int i = 0; i < 4; i++) {
                float n = sqrtf(g[i][0] * g[i][0] + g[i][1] * g[i][1]);
                if (n > Lg) Lg = n;
            }
        }
    }
    return Lg;
}

static float measured_edge_slope(const FineField* f) {
    float Le = 0.0f;
    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            if (x + 1 < N) {
                float r = fabsf(f->h[x + 1][y] - f->h[x][y]) / CELL;
                if (r > Le) Le = r;
            }
            if (y + 1 < N) {
                float r = fabsf(f->h[x][y + 1] - f->h[x][y]) / CELL;
                if (r > Le) Le = r;
            }
        }
    }
    return Le;
}

/* ----------------------------- self-test rig ------------------------------ */
static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

typedef struct {
    float x, y;
    float hr, hp, diff;
    const char* band;
} Probe;

#define SW 200

int main(void) {
    FineField fine;
    CoarseField coarse;
    fine_build(&fine);
    coarse_build(&coarse, &fine);

    const float c = CELL;
    const float d = coarse.cell;                  /* k*c = 2.0 m */
    const float delta = SNAP * sqrtf(0.5f);       /* |p - p'| <= s/sqrt2 */
    const float La = analytic_lipschitz();
    const float Lg = measured_gradient_lipschitz(&fine);
    const float Le = measured_edge_slope(&fine);

    /* interior sweep region: every point and its snapped image stay inside
     * both meshes, so no edge clamping participates (coarse mesh spans
     * [0, (g.n-1)*d]; keep a full cell of margin on every side). */
    const float lo = d;
    const float hi = (coarse.n - 1) * d - d;

    printf("== field ==\n");
    printf("fine_grid=%dx%d cell_c=%.3f m  coarse_stride=%d cell_d=%.3f m\n",
           N, N, c, STRIDE, d);
    printf("physics_snap_s=%.3f m  point_displacement_bound_delta=s/sqrt2=%.6f m\n",
           SNAP, delta);

    /* ---- Bound A: a-priori, generator F is La-Lipschitz.
     * Bilinear interpolation of L-Lipschitz vertex values reconstructs F
     * within L*spacing/sqrt2 (convex combination of four corners; weighted
     * mean corner distance maxes at the cell centre = spacing/sqrt2).
     *   B_a = La*(delta + (c+d)/sqrt2). */
    float Ba = La * (delta + (c + d) / sqrtf(2.0f));

    /* ---- Bound B: instantiated field, h_r is Lg-Lipschitz.
     * Coarse vertices are a subset of the fine vertices, so at a shared
     * point q, h_p(q) is a convex combination of coarse-corner values of
     * h_r, giving |h_r(q)-h_p(q)| <= Lg*d/sqrt2; plus the snap displacement
     *   |h_r(p)-h_p(p')| <= |h_r(p)-h_r(p')| + |h_r(p')-h_p(p')|
     *                     <= Lg*delta + Lg*d/sqrt2. */
    float Bm = Lg * (delta + d / sqrtf(2.0f));

    printf("\n== error upper bounds ==\n");
    printf("L_analytic(a_priori)=%.6f  L_gradient(measured)=%.6f  "
           "L_edge(measured)=%.6f\n", La, Lg, Le);
    printf("B_a=L_anal*(delta+(c+d)/sqrt2)=%.6f (a-priori)\n", Ba);
    printf("B_m=L_grad*(delta+d/sqrt2)=%.6f (instantiated field)\n", Bm);

    /* dense interior sweep: empirical worst case + height range */
    float worst = 0.0f, wx_ = 0.0f, wy_ = 0.0f;
    float hmin = 1e30f, hmax = -1e30f;
    for (int i = 0; i < SW; i++) {
        for (int j = 0; j < SW; j++) {
            float x = lo + (hi - lo) * i / (SW - 1);
            float y = lo + (hi - lo) * j / (SW - 1);
            float sx, sy;
            physics_snap(x, y, &sx, &sy);
            float hr = fine_sample(&fine, x, y);
            float e = fabsf(hr - coarse_sample(&coarse, sx, sy));
            if (e > worst) { worst = e; wx_ = x; wy_ = y; }
            if (hr < hmin) hmin = hr;
            if (hr > hmax) hmax = hr;
        }
    }
    printf("empirical_worst_diff=%.6f at (%.4f,%.4f)  render_height_range=[%.4f,%.4f]\n",
           worst, wx_, wy_, hmin, hmax);
    check("worst_inside_Bm", worst <= Bm + EPS);
    check("worst_inside_Ba", worst <= Ba + EPS);
    check("Bm_tighter_than_Ba", Bm <= Ba + EPS);

    /* ---- three probe positions: one per render-height tertile.
     * In each band pick the point with largest penetration (largest diff;
     * fall back to largest |diff| if the band has no penetration). */
    printf("\n== three probe positions (before unification) ==\n");
    const char* band_name[3] = { "low", "mid", "high" };
    Probe probes[3];
    memset(probes, 0, sizeof(probes));
    float best_abs[3] = { -1, -1, -1 }, best_pen[3] = { -1, -1, -1 };
    for (int i = 0; i < SW; i++) {
        for (int j = 0; j < SW; j++) {
            float x = lo + (hi - lo) * i / (SW - 1);
            float y = lo + (hi - lo) * j / (SW - 1);
            float sx, sy;
            physics_snap(x, y, &sx, &sy);
            float hr = fine_sample(&fine, x, y);
            float hp = coarse_sample(&coarse, sx, sy);
            float diff = hr - hp;
            int band = (int)((hr - hmin) / (hmax - hmin) * 3.0f);
            if (band > 2) band = 2;
            if (fabsf(diff) > best_abs[band]) {
                best_abs[band] = fabsf(diff);
                probes[band].x = x; probes[band].y = y;
                probes[band].hr = hr; probes[band].hp = hp; probes[band].diff = diff;
                probes[band].band = band_name[band];
            }
            /* prefer a penetrating point as the band representative */
            if (diff > EPS && diff > best_pen[band]) {
                best_pen[band] = diff;
                probes[band].x = x; probes[band].y = y;
                probes[band].hr = hr; probes[band].hp = hp; probes[band].diff = diff;
                probes[band].band = band_name[band];
            }
        }
    }
    int pen_before = 0, float_before = 0;
    for (int b = 0; b < 3; b++) {
        Probe* q = &probes[b];
        float sx, sy;
        physics_snap(q->x, q->y, &sx, &sy);
        printf("pos[%s] p=(%.3f,%.3f) snapped=(%.2f,%.2f) "
               "h_render=%.6f h_phys=%.6f height_diff=%.6f\n",
               q->band, q->x, q->y, sx, sy, q->hr, q->hp, q->diff);
        if (q->diff > EPS) pen_before++;
        else if (q->diff < -EPS) float_before++;
    }
    check("three_diff_inside_Bm",
          fabsf(probes[0].diff) <= Bm + EPS &&
          fabsf(probes[1].diff) <= Bm + EPS &&
          fabsf(probes[2].diff) <= Bm + EPS);
    check("three_bands_distinct",
          probes[0].hr < probes[1].hr && probes[1].hr < probes[2].hr);
    check("at_least_one_penetration_before", pen_before >= 1);

    /* full population: deterministic jittered bodies, penetration audit */
    printf("\n== population audit (jittered bodies) ==\n");
    int nbody = 0, pen_pop_before = 0, float_pop_before = 0;
    float max_pop_diff = 0.0f;
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            unsigned k = (unsigned)(i * 131 + j * 17) + 5u * 7919u;
            float u1 = (float)(k * 2654435761u & 0xffff) / 65535.0f;
            float u2 = (float)((k + 1) * 2654435761u & 0xffff) / 65535.0f;
            float x = lo + (i + u1) * 0.8f;
            float y = lo + (j + u2) * 0.8f;
            if (x > hi || y > hi) continue;
            float sx, sy;
            physics_snap(x, y, &sx, &sy);
            float diff = fine_sample(&fine, x, y) - coarse_sample(&coarse, sx, sy);
            nbody++;
            if (fabsf(diff) > max_pop_diff) max_pop_diff = fabsf(diff);
            if (diff > EPS) pen_pop_before++;
            else if (diff < -EPS) float_pop_before++;
        }
    }
    printf("bodies=%d penetration_before=%d floating_before=%d max_abs_diff=%.6f\n",
           nbody, pen_pop_before, float_pop_before, max_pop_diff);
    check("pop_diff_inside_Bm", max_pop_diff <= Bm + EPS);

    /* ---------------- unification: one sampler, one sample point ---------- */
    printf("\n== after unification (shared fine sampler, same world point) ==\n");
    int pen_after = 0, float_after = 0;
    for (int b = 0; b < 3; b++) {
        Probe* q = &probes[b];
        float hp_unified = fine_sample(&fine, q->x, q->y);
        float d2 = q->hr - hp_unified;
        printf("pos[%s] h_render=%.6f h_unified=%.6f height_diff=%.3e\n",
               q->band, q->hr, hp_unified, d2);
        if (d2 > EPS) pen_after++;
        else if (d2 < -EPS) float_after++;
    }
    check("three_probes_zero_penetration_after", pen_after == 0);

    int pen_pop_after = 0, float_pop_after = 0;
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            unsigned k = (unsigned)(i * 131 + j * 17) + 5u * 7919u;
            float u1 = (float)(k * 2654435761u & 0xffff) / 65535.0f;
            float u2 = (float)((k + 1) * 2654435761u & 0xffff) / 65535.0f;
            float x = lo + (i + u1) * 0.8f;
            float y = lo + (j + u2) * 0.8f;
            if (x > hi || y > hi) continue;
            float diff = fine_sample(&fine, x, y) - fine_sample(&fine, x, y);
            if (diff > EPS) pen_pop_after++;
            else if (diff < -EPS) float_pop_after++;
        }
    }
    printf("population penetration_after=%d floating_after=%d\n",
           pen_pop_after, float_pop_after);
    check("population_zero_penetration_after", pen_pop_after == 0);

    printf("\n== summary ==\n");
    printf("bounds: B_a=%.6f B_m=%.6f empirical_worst=%.6f\n", Ba, Bm, worst);
    printf("three_pos_height_diff: low=%.6f mid=%.6f high=%.6f\n",
           probes[0].diff, probes[1].diff, probes[2].diff);
    printf("penetration: three_pos before=%d after=%d | population(%d bodies) "
           "before=%d after=%d\n",
           pen_before, pen_after, nbody, pen_pop_before, pen_pop_after);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
