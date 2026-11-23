/* chunk: seams between adjacent heightfield chunks.
 *
 * Reproduces the crack between neighbouring chunks, states and checks the
 * conditions under which the boundary samples of two chunks coincide, then
 * unifies boundary sampling onto one shared lattice and measures the seams and
 * the cost before/after.  `make test` fails if any crack pixel survives.
 *
 * Crack pixel here = one seam sample (PPC per cell) whose 8-bit quantised
 * height differs between the two chunks meeting at that seam.  Four interior
 * junctions are sampled on all four arms for GRID=3.
 */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "noise.h"

#define CH   64     /* cells per chunk edge */
#define PPC  4      /* seam samples per cell (1 crack pixel = 1/PPC cell) */
#define ARM  16     /* junction arm length, in cells */
#define GRID 3      /* GRID x GRID chunks => 4 interior junctions */

/* Deliberately awkward but legal placement: at this base the two coordinate
 * construction paths round differently at some seams. */
#define BASE (1048576.0f) /* 2^20 */
#define CELL 0.30f
#define FREQ 0.35f
#define AMP  40.0f
#define OCT  5
#define SEED 1337u

typedef enum { M_NAIVE, M_DUP, M_SHARED, M_COUNT } Mode;
static const char* const MODE_NAME[M_COUNT] = { "naive", "dup", "shared" };

typedef struct {
    int n;              /* chunks per side */
    float base, cell, freq, amp;
    int oct;
    unsigned seed;
    int gn;             /* shared lattice edge = n*CH + 1 */
    float* naive;       /* n*n*CH*CH      (no boundary column) */
    float* dup;         /* n*n*(CH+1)^2   (each chunk duplicates its edges) */
    float* shared;      /* gn*gn          (one global lattice) */
} World;

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static float height(const World* w, float wx, float wy) {
    return fbm(wx * w->freq, wy * w->freq, w->oct, 2.0f, 0.5f, w->seed) * w->amp;
}

/* Single canonical path: world coordinate of global lattice index g. */
static float canon_w(const World* w, int g) { return w->base + (float)g * w->cell; }

/* Per-chunk path: origin_c + l*cell.  Each origin rounds on its own, so this
 * is not bit-identical to canon_w() at every seam. */
static float chunk_w(const World* w, int c, int l) {
    return (w->base + (float)(c * CH) * w->cell) + (float)l * w->cell;
}

static void world_init(World* w, int n, float base, float cell) {
    memset(w, 0, sizeof(*w));
    w->n = n;
    w->base = base;
    w->cell = cell;
    w->freq = FREQ;
    w->amp = AMP;
    w->oct = OCT;
    w->seed = SEED;
    w->gn = n * CH + 1;
    w->naive = malloc(sizeof(float) * (size_t)n * n * CH * CH);
    w->dup = malloc(sizeof(float) * (size_t)n * n * (CH + 1) * (CH + 1));
    w->shared = malloc(sizeof(float) * (size_t)w->gn * w->gn);
    if (!w->naive || !w->dup || !w->shared) {
        fprintf(stderr, "out of memory\n");
        exit(2);
    }
}

static void world_free(World* w) {
    free(w->naive);
    free(w->dup);
    free(w->shared);
}

static void fill_naive(World* w) {
    for (int cy = 0; cy < w->n; cy++)
        for (int cx = 0; cx < w->n; cx++)
            for (int ly = 0; ly < CH; ly++)
                for (int lx = 0; lx < CH; lx++) {
                    size_t i = (((size_t)(cy * w->n + cx) * CH) + ly) * CH + lx;
                    w->naive[i] = height(w, chunk_w(w, cx, lx), chunk_w(w, cy, ly));
                }
}

static void fill_dup(World* w) {
    for (int cy = 0; cy < w->n; cy++)
        for (int cx = 0; cx < w->n; cx++)
            for (int ly = 0; ly <= CH; ly++)
                for (int lx = 0; lx <= CH; lx++) {
                    size_t i = (((size_t)(cy * w->n + cx) * (CH + 1)) + ly) * (CH + 1) + lx;
                    w->dup[i] = height(w, chunk_w(w, cx, lx), chunk_w(w, cy, ly));
                }
}

static void fill_shared(World* w) {
    for (int gy = 0; gy < w->gn; gy++)
        for (int gx = 0; gx < w->gn; gx++)
            w->shared[(size_t)gy * w->gn + gx] = height(w, canon_w(w, gx), canon_w(w, gy));
}

static void build(World* w, Mode m) {
    if (m == M_NAIVE) fill_naive(w);
    else if (m == M_DUP) fill_dup(w);
    else fill_shared(w);
}

/* Bilinear fetch in local coordinates.  Each mode's addressing mirrors what a
 * real chunked terrain module does: per-chunk storage, local index, clamped. */
static float sample(const World* w, Mode m, int cx, int cy, float lx, float ly) {
    int x0 = (int)floorf(lx), y0 = (int)floorf(ly);
    int x1 = x0 + 1, y1 = y0 + 1;
    float tx = lx - (float)x0, ty = ly - (float)y0;
    int lim = 0, stride = 0;
    const float* f = NULL;

    if (m == M_NAIVE) {
        lim = CH - 1; stride = CH;
        f = w->naive + (size_t)(cy * w->n + cx) * CH * CH;
    } else if (m == M_DUP) {
        lim = CH; stride = CH + 1;
        f = w->dup + (size_t)(cy * w->n + cx) * (CH + 1) * (CH + 1);
    } else {
        lim = CH; stride = w->gn;
        f = w->shared + (size_t)(cy * CH) * w->gn + (size_t)(cx * CH);
    }
    if (x0 < 0) x0 = 0;
    if (x0 > lim) x0 = lim;
    if (x1 < 0) x1 = 0;
    if (x1 > lim) x1 = lim;
    if (y0 < 0) y0 = 0;
    if (y0 > lim) y0 = lim;
    if (y1 < 0) y1 = 0;
    if (y1 > lim) y1 = lim;

    float a = f[y0 * stride + x0] * (1 - tx) + f[y0 * stride + x1] * tx;
    float b = f[y1 * stride + x0] * (1 - tx) + f[y1 * stride + x1] * tx;
    return a * (1 - ty) + b * ty;
}

static int q8(float h, float amp) {
    float v = (h / (2.0f * amp) + 0.5f) * 255.0f;
    int i = (int)lrintf(v);
    return i < 0 ? 0 : (i > 255 ? 255 : i);
}

/* ---- Coincidence conditions -------------------------------------------
 *
 * A seam sample at global index G sits between the left chunk's local index CH
 * and the right chunk's local index 0, at the same tangential parameter t.  The
 * two rendered samples coincide (no crack) iff:
 *
 * C1  the boundary node is stored and addressable in both chunks: the left
 *     chunk must keep local index CH and the right chunk local index 0.  Naive
 *     chunks keep 0..CH-1 only, so a seam query is clamped a whole cell in.
 *
 * C2  both sides evaluate the height at the same world coordinate with
 *     identical bits: x_left = origin_c + CH*cell must equal x_right =
 *     origin_{c+1} + 0*cell, and both must equal the canonical single-path
 *     coordinate base + G*cell.  Equal index is not enough: per-chunk origins
 *     round independently, so the two expressions can differ by an ulp and pick
 *     different noise lattice corners / weights.
 *
 * C3  the boundary heights are then bit-identical on both sides (same height
 *     function, same argument, deterministic hash) at the same t, so the
 *     interpolation yields the same bits.
 *
 * C1 + C2 are the conditions to enforce; C3 is their observable consequence.
 * Unified sampling satisfies all three by addressing one shared lattice by
 * global index and computing each coordinate through canon_w() only.
 */

static int cond_boundary_stored(const World* w, Mode m) {
    (void)w;
    return m != M_NAIVE;
}

static int cond_same_coord(const World* w, Mode m) {
    if (m == M_SHARED) return 1; /* shared lattice is addressed only canonically */
    for (int c = 0; c < w->n - 1; c++) {
        int gi = (c + 1) * CH;
        float left = cond_boundary_stored(w, m) ? chunk_w(w, c, CH) : chunk_w(w, c, CH - 1);
        float right = chunk_w(w, c + 1, 0);
        float canon = canon_w(w, gi);
        if (left != right || left != canon) return 0;
    }
    return 1;
}

static int cond_same_value(const World* w, Mode m) {
    for (int c = 0; c < w->n - 1; c++)
        for (int l = 0; l <= CH; l++) {
            float xa, xb, ya, yb;
            if (!cond_boundary_stored(w, m)) {
                xa = sample(w, m, c, 0, CH - 1, (float)l);
                xb = sample(w, m, c + 1, 0, 0.0f, (float)l);
                ya = sample(w, m, 0, c, (float)l, CH - 1);
                yb = sample(w, m, 0, c + 1, (float)l, 0.0f);
            } else {
                xa = sample(w, m, c, 0, (float)CH, (float)l);
                xb = sample(w, m, c + 1, 0, 0.0f, (float)l);
                ya = sample(w, m, 0, c, (float)l, (float)CH);
                yb = sample(w, m, 0, c + 1, (float)l, 0.0f);
            }
            if (xa != xb || ya != yb) return 0;
        }
    return 1;
}

/* ---- Seam measurement -------------------------------------------------- */

typedef struct { long long cracks, samples; int max_lsb; float max_dh; } SeamStat;

static void seam_pair(const World* w, Mode m,
                      int ax, int ay, float alx, float aly,
                      int bx, int by, float blx, float bly,
                      SeamStat* st) {
    float ha = sample(w, m, ax, ay, alx, aly);
    float hb = sample(w, m, bx, by, blx, bly);
    int d = abs(q8(ha, w->amp) - q8(hb, w->amp));
    st->samples++;
    if (d > 0) st->cracks++;
    if (d > st->max_lsb) st->max_lsb = d;
    float dh = fabsf(ha - hb);
    if (dh > st->max_dh) st->max_dh = dh;
}

/* One junction = four seam arms: lower/upper vertical, left/right horizontal. */
static SeamStat junction_stat(const World* w, Mode m, int jx, int jy) {
    SeamStat st = { 0, 0, 0, 0.0f };
    int a = jx - 1, b = jx, c = jy - 1, d = jy;
    for (int p = 0; p <= ARM * PPC; p++) {
        float t = (float)(CH - ARM) + (float)p / PPC; /* lower / left arm */
        float u = (float)p / PPC;                     /* upper / right arm */
        seam_pair(w, m, a, c, (float)CH, t, b, c, 0.0f, t, &st);
        seam_pair(w, m, a, d, (float)CH, u, b, d, 0.0f, u, &st);
        seam_pair(w, m, a, c, t, (float)CH, a, d, t, 0.0f, &st);
        seam_pair(w, m, b, c, u, (float)CH, b, d, u, 0.0f, &st);
    }
    return st;
}

typedef struct { long long cracks; int max_lsb; float max_dh; } ModeTotals;

static ModeTotals run_mode(const World* w, Mode m, SeamStat j[2][2]) {
    ModeTotals t = { 0, 0, 0.0f };
    for (int jy = 1; jy < w->n; jy++)
        for (int jx = 1; jx < w->n; jx++) {
            j[jy - 1][jx - 1] = junction_stat(w, m, jx, jy);
            t.cracks += j[jy - 1][jx - 1].cracks;
            if (j[jy - 1][jx - 1].max_lsb > t.max_lsb) t.max_lsb = j[jy - 1][jx - 1].max_lsb;
            if (j[jy - 1][jx - 1].max_dh > t.max_dh) t.max_dh = j[jy - 1][jx - 1].max_dh;
        }
    return t;
}

static void run_preset(const char* tag, float base, float cell, int hostile) {
    World w;
    world_init(&w, GRID, base, cell);
    fill_naive(&w);
    fill_dup(&w);
    fill_shared(&w);

    const int per_junction = 4 * (ARM * PPC + 1);
    printf("preset=%s base=%g cell=%g chunks=%dx%d seam_px_per_junction=%d\n",
           tag, (double)base, (double)cell, GRID, GRID, per_junction);

    ModeTotals tot[M_COUNT];
    SeamStat js[M_COUNT][2][2];
    for (int m = 0; m < M_COUNT; m++) tot[m] = run_mode(&w, (Mode)m, js[m]);

    for (int m = 0; m < M_COUNT; m++)
        printf("mode=%s cond_boundary_stored=%d cond_same_coord=%d cond_same_value=%d "
               "crack_J11=%lld crack_J21=%lld crack_J12=%lld "
               "crack_J22=%lld crack_total=%lld max_lsb=%d max_dh=%.3f\n",
               MODE_NAME[m],
               cond_boundary_stored(&w, (Mode)m), cond_same_coord(&w, (Mode)m),
               cond_same_value(&w, (Mode)m),
               js[m][0][0].cracks, js[m][0][1].cracks,
               js[m][1][0].cracks, js[m][1][1].cracks,
               tot[m].cracks, tot[m].max_lsb, (double)tot[m].max_dh);

    check("bug_reproduced_naive_cracks", tot[M_NAIVE].cracks > 0);
    char name[64];
    for (int jy = 0; jy < 2; jy++)
        for (int jx = 0; jx < 2; jx++) {
            snprintf(name, sizeof(name), "no_crack_shared_J%d%d", jx + 1, jy + 1);
            check(name, js[M_SHARED][jy][jx].cracks == 0);
        }
    check("no_crack_shared_total", tot[M_SHARED].cracks == 0);
    check("shared_meets_all_conditions",
          cond_boundary_stored(&w, M_SHARED) && cond_same_coord(&w, M_SHARED) &&
          cond_same_value(&w, M_SHARED));

    if (hostile) {
        /* Duplicated edges alone are insufficient: per-chunk origins round
         * separately, so condition C2 fails and cracks remain.  This is what
         * makes the condition non-vacuous.  Under -ffast-math the compiler is
         * allowed to fold the two coordinate paths together, which suppresses
         * the ulp difference; skip the demonstration there rather than fail. */
#if defined(__FAST_MATH__)
        printf("note=fast_math_skips_dup_insufficiency\n");
#else
        check("dup_violates_coord_condition", !cond_same_coord(&w, M_DUP));
        check("dup_leaves_cracks", tot[M_DUP].cracks > 0);
#endif
    } else {
        check("dup_sufficient_when_coords_agree", tot[M_DUP].cracks == 0);
    }
    printf("\n");
    world_free(&w);
}

/* Timing.  All three modes run inside one process, interleaved round-robin, so
 * scheduling noise affects them equally; we report the median of a batch and
 * the median absolute deviation so a small delta can be read against its own
 * spread rather than against a single lucky run. */

typedef struct { double med, mad; } Timing;

static int cmp_double(const void* a, const void* b) {
    double x = *(const double*)a, y = *(const double*)b;
    return (x > y) - (x < y);
}

static Timing summarise(double* v, int n) {
    qsort(v, (size_t)n, sizeof(double), cmp_double);
    Timing t;
    t.med = (n & 1) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
    double* d = malloc(sizeof(double) * (size_t)n);
    for (int i = 0; i < n; i++) d[i] = fabs(v[i] - t.med);
    qsort(d, (size_t)n, sizeof(double), cmp_double);
    t.mad = (n & 1) ? d[n / 2] : 0.5 * (d[n / 2 - 1] + d[n / 2]);
    free(d);
    return t;
}

#define REPS 40
static volatile float g_sink;

/* Prevent the optimiser from hoisting seam samples out of the timing loop.
 * Applied identically to all three modes, so it does not bias the comparison. */
#if defined(__GNUC__) || defined(__clang__)
#define COMPILER_BARRIER() __asm__ volatile("" ::: "memory")
#else
#define COMPILER_BARRIER() (g_sink = g_sink)
#endif

static void time_all(World* w, Timing out[M_COUNT], int seams) {
    double* buf = malloc(sizeof(double) * REPS);
    double* accs[M_COUNT];
    for (int m = 0; m < M_COUNT; m++) accs[m] = malloc(sizeof(double) * REPS);

    for (int r = 0; r < REPS; r++) {
        for (int m = 0; m < M_COUNT; m++) {
            float acc = 0.0f;
            double t0 = now_sec();
            if (!seams) {
                build(w, (Mode)m);
            } else {
                for (int jy = 1; jy < w->n; jy++)
                    for (int jx = 1; jx < w->n; jx++) {
                        int a = jx - 1, b = jx, c = jy - 1, d = jy;
                        for (int p = 0; p <= ARM * PPC; p++) {
                            float t = (float)(CH - ARM) + (float)p / PPC;
                            float u = (float)p / PPC;
                            acc += sample(w, (Mode)m, a, c, (float)CH, t);
                            acc += sample(w, (Mode)m, b, c, 0.0f, t);
                            acc += sample(w, (Mode)m, a, d, (float)CH, u);
                            acc += sample(w, (Mode)m, b, d, 0.0f, u);
                            acc += sample(w, (Mode)m, a, c, t, (float)CH);
                            acc += sample(w, (Mode)m, a, d, t, 0.0f);
                            acc += sample(w, (Mode)m, b, c, u, (float)CH);
                            acc += sample(w, (Mode)m, b, d, u, 0.0f);
                            COMPILER_BARRIER();
                        }
                    }
                g_sink = acc;
            }
            accs[m][r] = now_sec() - t0;
        }
    }
    for (int m = 0; m < M_COUNT; m++) {
        memcpy(buf, accs[m], sizeof(double) * REPS);
        out[m] = summarise(buf, REPS);
    }
    for (int m = 0; m < M_COUNT; m++) free(accs[m]);
    free(buf);
}

int main(void) {
    run_preset("offset_2^20", BASE, CELL, 1);
    run_preset("offset_0", 0.0f, 0.35f, 0);

    World w;
    world_init(&w, GRID, BASE, CELL);
    /* Fill all three layouts before timing, so seam queries read real data. */
    fill_naive(&w);
    fill_dup(&w);
    fill_shared(&w);

    const long long nodes_naive = (long long)GRID * GRID * CH * CH;
    const long long nodes_dup = (long long)GRID * GRID * (CH + 1) * (CH + 1);
    const long long nodes_shared = (long long)w.gn * w.gn;
    const long long seam_samples = (long long)(GRID - 1) * (GRID - 1) * 4 * (ARM * PPC + 1) * 2;

    Timing tb[M_COUNT], ts[M_COUNT];
    time_all(&w, tb, 0);
    time_all(&w, ts, 1);

    printf("timing grid=%dx%d ch=%d seam_samples=%lld reps=%d (median +/- MAD, round-robin)\n",
           GRID, GRID, CH, seam_samples, REPS);
    printf("nodes_meshes naive=%lld dup=%lld shared=%lld dup_overhead_pct=%.2f "
           "shared_overhead_vs_naive_pct=%.2f\n",
           nodes_naive, nodes_dup, nodes_shared,
           (double)(nodes_dup - nodes_naive) / nodes_naive * 100.0,
           (double)(nodes_shared - nodes_naive) / nodes_naive * 100.0);
    double b_naive = tb[M_NAIVE].med, b_dup = tb[M_DUP].med, b_shared = tb[M_SHARED].med;
    double b_mad = fmax(tb[M_NAIVE].mad, tb[M_SHARED].mad);
    double bd_n = b_shared - b_naive, bd_d = b_shared - b_dup;
    printf("build_ms naive=%.3f dup=%.3f shared=%.3f (MAD<=%.3f) "
           "delta_shared_vs_naive_ms=%+.3f%s delta_shared_vs_dup_ms=%+.3f%s\n",
           b_naive * 1e3, b_dup * 1e3, b_shared * 1e3, b_mad * 1e3,
           bd_n * 1e3, fabs(bd_n) <= b_mad ? "(within noise)" : "",
           bd_d * 1e3, fabs(bd_d) <= b_mad ? "(within noise)" : "");

    double ns_naive = ts[M_NAIVE].med / (double)seam_samples * 1e9;
    double ns_dup = ts[M_DUP].med / (double)seam_samples * 1e9;
    double ns_shared = ts[M_SHARED].med / (double)seam_samples * 1e9;
    double mad_bound = fmax(ts[M_NAIVE].mad, ts[M_SHARED].mad) / (double)seam_samples * 1e9;
    double d_n = ns_shared - ns_naive, d_d = ns_shared - ns_dup;
    printf("seam_query_ns naive=%.2f dup=%.2f shared=%.2f (MAD<=%.2f) "
           "delta_shared_vs_naive_ns=%+.2f%s delta_shared_vs_dup_ns=%+.2f%s\n",
           ns_naive, ns_dup, ns_shared, mad_bound,
           d_n, fabs(d_n) <= mad_bound ? "(within noise)" : "",
           d_d, fabs(d_d) <= mad_bound ? "(within noise)" : "");
    world_free(&w);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
