#include <stdlib.h>

#include "parallel.h"
#include "rng.h"
#include "worldgen.h"

#ifdef SP_SELFTEST_NONDET
/* Defined in src/selftest_nondet.c; negative-control build only. */
float sp_selftest_nondet_mean(const float *v, size_t n, int threads);
#endif

/* Water flow: every cell drains toward its lowest 8-neighbour; flow
 * accumulation counts upstream cells.  The accumulation is computed with
 * WG_HYDRO_ITERS Jacobi passes (A_0 = 0, A_{t+1} = 1 + sum of upstream A_t),
 * so each pass reads only the previous buffer and no reduction crosses
 * threads.
 *
 * Drainage ties are broken by a small per-cell bias drawn from the hydro
 * sub-seed at global coordinates.  Without it, a flat cell would pick its
 * neighbour by iteration order and windows would disagree at seams. */

static const int DIRS_DX[8] = {1, -1, 0, 0, 1, -1, 1, -1};
static const int DIRS_DY[8] = {0, 0, 1, -1, 1, 1, -1, -1};
static const int DIRS_OPP[8] = {1, 0, 3, 2, 7, 6, 5, 4};

static float hydro_bias(int wx, int wy, uint64_t seed) {
    return 1e-3f * sp_value_noise((float)wx * 0.19f + 71.3f,
                                  (float)wy * 0.19f - 23.9f,
                                  sp_mix64(seed ^ UINT64_C(0x13198A2E03707344)));
}

typedef struct {
    const SpGrid *h;
    unsigned char *dir;
    uint64_t seed;
} DirCtx;

static void dir_rows(void *vp, int y) {
    DirCtx *c = (DirCtx *)vp;
    const SpGrid *g = c->h;
    for (int x = 0; x < g->w; x++) {
        float h0 = sp_grid_at(g, x, y) +
                   hydro_bias(g->ox + x, g->oy + y, c->seed);
        float best = h0;
        unsigned char bestdir = 255;
        for (int k = 0; k < 8; k++) {
            int nx = x + DIRS_DX[k], ny = y + DIRS_DY[k];
            if (nx < 0) nx = 0;
            if (ny < 0) ny = 0;
            if (nx >= g->w) nx = g->w - 1;
            if (ny >= g->h) ny = g->h - 1;
            if (nx == x && ny == y) continue;
            float hn = sp_grid_at(g, nx, ny) +
                       hydro_bias(g->ox + nx, g->oy + ny, c->seed);
            if (hn < best) {
                best = hn;
                bestdir = (unsigned char)k;
            }
        }
        c->dir[(size_t)y * (size_t)g->w + (size_t)x] = bestdir;
    }
}

typedef struct {
    const float *cur;
    float *next;
    const unsigned char *dir;
    int w, h;
} AccumCtx;

/* Fixed neighbour order k = 0..7; the sum is a fixed-order chain within one
 * thread, never a parallel reduction. */
static void accum_rows(void *vp, int y) {
    AccumCtx *c = (AccumCtx *)vp;
    for (int x = 0; x < c->w; x++) {
        float sum = 1.0f;
        for (int k = 0; k < 8; k++) {
            int nx = x + DIRS_DX[k], ny = y + DIRS_DY[k];
            if (nx < 0) nx = 0;
            if (ny < 0) ny = 0;
            if (nx >= c->w) nx = c->w - 1;
            if (ny >= c->h) ny = c->h - 1;
            if (nx == x && ny == y) continue;
            if (c->dir[(size_t)ny * (size_t)c->w + (size_t)nx] ==
                (unsigned char)DIRS_OPP[k]) {
                sum += c->cur[(size_t)ny * (size_t)c->w + (size_t)nx];
            }
        }
        c->next[(size_t)y * (size_t)c->w + (size_t)x] = sum;
    }
}

static float *flow_accumulation(const SpGrid *h, uint64_t seed, int threads) {
    size_t n = (size_t)h->w * (size_t)h->h;
    unsigned char *dir = (unsigned char *)malloc(n);
    float *cur = (float *)calloc(n, sizeof(float));
    float *next = (float *)malloc(n * sizeof(float));
    if (!dir || !cur || !next) abort();

    DirCtx dc = {h, dir, seed};
    sp_parallel_rows(h->h, threads, dir_rows, &dc);

    for (int it = 0; it < WG_HYDRO_ITERS; it++) {
        AccumCtx ac = {cur, next, dir, h->w, h->h};
        sp_parallel_rows(h->h, threads, accum_rows, &ac);
        float *tmp = cur;
        cur = next;
        next = tmp;
    }
    free(dir);
    free(next);
    return cur;
}

int wg_layer_hydro(SpGrid *out, const SpRun *run) {
    if (!sp_grid_alloc(out, 0, 0, run->win.w, run->win.h)) return 0;
    SpGrid region;
    if (!wg_build_height(&region, run, wg_halo_for_layer(SP_LAYER_HYDRO))) {
        sp_grid_free(out);
        return 0;
    }
    wg_erode(&region, run->sub[SP_LAYER_EROSION], run->threads);
    float *flow = flow_accumulation(&region, run->sub[SP_LAYER_HYDRO],
                                    run->threads);
#ifdef SP_SELFTEST_NONDET
    /* NEGATIVE CONTROL (see src/selftest_nondet.c): scale by a mean whose
     * value depends on thread count and completion order. */
    float norm = sp_selftest_nondet_mean(flow,
                                         (size_t)region.w * (size_t)region.h,
                                         run->threads);
    if (norm != 0.0f) {
        for (size_t i = 0; i < (size_t)region.w * (size_t)region.h; i++) {
            flow[i] /= norm;
        }
    }
#endif
    int lx0 = run->win.x0 - region.ox, ly0 = run->win.y0 - region.oy;
    for (int y = 0; y < run->win.h; y++) {
        for (int x = 0; x < run->win.w; x++) {
            out->v[(size_t)y * (size_t)run->win.w + (size_t)x] =
                flow[(size_t)(ly0 + y) * (size_t)region.w + (size_t)(lx0 + x)];
        }
    }
    free(flow);
    sp_grid_free(&region);
    return 1;
}
