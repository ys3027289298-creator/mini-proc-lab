#include <math.h>
#include <stdlib.h>

#include "parallel.h"
#include "rng.h"
#include "worldgen.h"

/* Thermal erosion.  Each cell that is more than the local talus above its
 * lowest 4-neighbour moves half of that excess downhill, to that single
 * neighbour.  A cell's new height is computed from the previous buffer only:
 *
 *     next = h - outflow + sum(inflow from neighbours that chose me)
 *
 * A cell receives material only from neighbours whose own lowest neighbour is
 * this cell, determined by comparing clamped coordinates rather than by an
 * offset table, so world-edge clamping cannot break conservation.
 *
 * Sending to a single destination is what keeps the pass stable: moving
 * half the excess in several directions at once can overshoot and oscillate.
 *
 * The talus field comes from the erosion sub-seed but is evaluated at global
 * coordinates, so it does not shift with the window. */
float wg_talus_at(int wx, int wy, uint64_t erosion_seed) {
    float n = sp_value_noise((float)wx * 0.07f + 3.1f,
                             (float)wy * 0.07f - 5.7f,
                             sp_mix64(erosion_seed ^
                                      UINT64_C(0x4F1BBCDCBFA54127)));
    return WG_TALUS * (0.85f + 0.30f * (n * 0.5f + 0.5f));
}

static const int DX[4] = {-1, 1, 0, 0};
static const int DY[4] = {0, 0, -1, 1};

typedef struct {
    int x, y;
} Cell;

static Cell clamp_cell(int x, int y, int w, int h) {
    Cell c;
    c.x = x < 0 ? 0 : (x >= w ? w - 1 : x);
    c.y = y < 0 ? 0 : (y >= h ? h - 1 : y);
    return c;
}

/* Lowest non-self 4-neighbour of (x, y) after edge clamping.  Ties are broken
 * by the fixed direction order, so the choice is reproducible. */
static Cell lowest_neighbour(const float *cur, int w, int h, int x, int y) {
    Cell best = {x, y};
    float best_h = cur[(size_t)y * (size_t)w + (size_t)x];
    for (int k = 0; k < 4; k++) {
        Cell n = clamp_cell(x + DX[k], y + DY[k], w, h);
        if (n.x == x && n.y == y) continue;
        float hn = cur[(size_t)n.y * (size_t)w + (size_t)n.x];
        if (hn < best_h) {
            best_h = hn;
            best = n;
        }
    }
    return best;
}

typedef struct {
    const float *cur;
    float *next;
    int w, h, ox, oy;
    uint64_t seed;
} ErodeCtx;

static void erode_rows(void *vp, int y) {
    ErodeCtx *c = (ErodeCtx *)vp;
    const float *cur = c->cur;
    for (int x = 0; x < c->w; x++) {
        float h0 = cur[(size_t)y * (size_t)c->w + (size_t)x];
        float talus = wg_talus_at(c->ox + x, c->oy + y, c->seed);

        float out = 0.0f;
        Cell low = lowest_neighbour(cur, c->w, c->h, x, y);
        if (low.x != x || low.y != y) {
            float d = h0 - cur[(size_t)low.y * (size_t)c->w + (size_t)low.x];
            if (d > talus) out = 0.5f * (d - talus);
        }

        /* Fixed neighbour order; the accumulator is a per-cell float chain,
         * never a cross-thread reduction. */
        float in = 0.0f;
        for (int k = 0; k < 4; k++) {
            Cell n = clamp_cell(x + DX[k], y + DY[k], c->w, c->h);
            if (n.x == x && n.y == y) continue;
            Cell nlow = lowest_neighbour(cur, c->w, c->h, n.x, n.y);
            if (nlow.x != x || nlow.y != y) continue;
            float hn = cur[(size_t)n.y * (size_t)c->w + (size_t)n.x];
            float d = hn - h0;
            if (d > 0.0f) {
                float tn = wg_talus_at(c->ox + n.x, c->oy + n.y, c->seed);
                if (d > tn) in += 0.5f * (d - tn);
            }
        }
        c->next[(size_t)y * (size_t)c->w + (size_t)x] = h0 - out + in;
    }
}

void wg_erode(SpGrid *r, uint64_t erosion_seed, int threads) {
    size_t n = (size_t)r->w * (size_t)r->h;
    float *base = r->v;
    float *cur = base;
    float *next = (float *)malloc(n * sizeof(float));
    if (!next) abort();
    for (int it = 0; it < WG_EROSION_ITERS; it++) {
        ErodeCtx c = {cur, next, r->w, r->h, r->ox, r->oy, erosion_seed};
        sp_parallel_rows(r->h, threads, erode_rows, &c);
        float *tmp = cur;
        cur = next;
        next = tmp;
    }
    /* `next` is always the buffer that is not the result, so freeing it also
     * releases the original allocation when an odd iteration count leaves the
     * result in the freshly allocated buffer. */
    r->v = cur;
    free(next);
}

float wg_slope_at(const SpGrid *r, int lx, int ly) {
    float dx = (sp_grid_at_clamped(r, lx + 1, ly) -
                sp_grid_at_clamped(r, lx - 1, ly)) * 0.5f;
    float dy = (sp_grid_at_clamped(r, lx, ly + 1) -
                sp_grid_at_clamped(r, lx, ly - 1)) * 0.5f;
    return sqrtf(dx * dx + dy * dy);
}

int wg_layer_erosion(SpGrid *out, const SpRun *run) {
    if (!sp_grid_alloc(out, 0, 0, run->win.w, run->win.h)) return 0;
    SpGrid r;
    if (!wg_build_height(&r, run, wg_halo_for_layer(SP_LAYER_EROSION))) {
        sp_grid_free(out);
        return 0;
    }
    wg_erode(&r, run->sub[SP_LAYER_EROSION], run->threads);
    int lx0 = run->win.x0 - r.ox, ly0 = run->win.y0 - r.oy;
    for (int y = 0; y < run->win.h; y++) {
        for (int x = 0; x < run->win.w; x++) {
            out->v[(size_t)y * (size_t)run->win.w + (size_t)x] =
                sp_grid_at(&r, lx0 + x, ly0 + y);
        }
    }
    sp_grid_free(&r);
    return 1;
}
