#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "parallel.h"
#include "rng.h"
#include "worldgen.h"

/* Jittered-grid scatter.
 *
 * Each candidate cell (ci, cj) gets a jittered point from the scatter
 * sub-seed.  A candidate is accepted when no strictly better-ranked valid
 * candidate lies within WG_SCATTER_MIND.  "Best in the neighbourhood wins" is
 * a local, order-independent rule: acceptance of a point depends only on the
 * cells within a fixed radius, never on the order cells are visited.  That is
 * what lets windows and thread counts agree.
 *
 * Candidate generation is split across rows; candidate storage is indexed by
 * (cj - cj0) * ncw + (ci - ci0), so no thread ever writes the same slot.  The
 * final list is gathered sequentially in (cj, ci) order. */

#define SEARCH_CELLS 2
#define MARGIN_CELLS 3

typedef struct {
    float x, y;
    uint64_t rank;
    int ci, cj;
    unsigned char valid;
    unsigned char accepted;
} Candidate;

typedef struct {
    Candidate *cand;
    int ci0, cj0, ncw, nch;
    const SpGrid *eroded;
    uint64_t seed;
    int ox, oy;
} ScatterCtx;

static Candidate *cand_at(const ScatterCtx *s, int ci, int cj) {
    return &s->cand[(size_t)(cj - s->cj0) * (size_t)s->ncw +
                    (size_t)(ci - s->ci0)];
}

static int in_range(const ScatterCtx *s, int ci, int cj) {
    return ci >= s->ci0 && ci < s->ci0 + s->ncw &&
           cj >= s->cj0 && cj < s->cj0 + s->nch;
}

static void point_of(uint64_t seed, int ci, int cj, float *x, float *y,
                     uint64_t *rank) {
    uint64_t h = sp_hash2d(ci, cj, seed);
    *x = ((float)ci + sp_unit(h)) * WG_SCATTER_CELL;
    *y = ((float)cj + sp_unit(sp_mix64(h))) * WG_SCATTER_CELL;
    *rank = sp_mix64(h ^ UINT64_C(0x243F6A8885A308D3));
}

/* Valid = above the shoreline and not on a cliff.  Both tests sample the
 * eroded height at global coordinates. */
static int valid_at(const ScatterCtx *s, float x, float y) {
    float h = sp_grid_sample(s->eroded, x, y);
    if (h < WG_SEA_LEVEL + 1.5f) return 0;
    float hx1 = sp_grid_sample(s->eroded, x + 1.0f, y);
    float hx0 = sp_grid_sample(s->eroded, x - 1.0f, y);
    float hy1 = sp_grid_sample(s->eroded, x, y + 1.0f);
    float hy0 = sp_grid_sample(s->eroded, x, y - 1.0f);
    float dx = (hx1 - hx0) * 0.5f, dy = (hy1 - hy0) * 0.5f;
    return (dx * dx + dy * dy) <= 2.5f * 2.5f;
}

static void prepare_rows(void *vp, int row) {
    ScatterCtx *s = (ScatterCtx *)vp;
    int cj = s->cj0 + row;
    for (int ci = s->ci0; ci < s->ci0 + s->ncw; ci++) {
        Candidate *c = cand_at(s, ci, cj);
        memset(c, 0, sizeof(*c));
        c->ci = ci;
        c->cj = cj;
        point_of(s->seed, ci, cj, &c->x, &c->y, &c->rank);
        c->valid = (unsigned char)valid_at(s, c->x, c->y);
    }
}

static void accept_rows(void *vp, int row) {
    ScatterCtx *s = (ScatterCtx *)vp;
    int cj = s->cj0 + row;
    for (int ci = s->ci0; ci < s->ci0 + s->ncw; ci++) {
        Candidate *c = cand_at(s, ci, cj);
        if (!c->valid) continue;
        int win = 1;
        for (int dj = -SEARCH_CELLS; dj <= SEARCH_CELLS && win; dj++) {
            for (int di = -SEARCH_CELLS; di <= SEARCH_CELLS; di++) {
                if (di == 0 && dj == 0) continue;
                int ni = ci + di, nj = cj + dj;
                if (!in_range(s, ni, nj)) continue;
                Candidate *o = cand_at(s, ni, nj);
                if (!o->valid) continue;
                float dx = o->x - c->x, dy = o->y - c->y;
                if (dx * dx + dy * dy >= WG_SCATTER_MIND * WG_SCATTER_MIND) continue;
                if (o->rank < c->rank ||
                    (o->rank == c->rank &&
                     (o->cj < cj || (o->cj == cj && ni < ci)))) {
                    win = 0;
                    break;
                }
            }
        }
        c->accepted = (unsigned char)win;
    }
}

int wg_layer_scatter(SpScatter *out, const SpRun *run) {
    int halo = wg_halo_for_layer(SP_LAYER_SCATTER);
    SpGrid region;
    if (!wg_build_height(&region, run, halo)) return 0;
    wg_erode(&region, run->sub[SP_LAYER_EROSION], run->threads);

    int ncells = WG_WORLD / (int)WG_SCATTER_CELL;
    int ci0 = (int)floorf((float)run->win.x0 / WG_SCATTER_CELL) - MARGIN_CELLS;
    int ci1 = (int)floorf((float)(run->win.x0 + run->win.w) / WG_SCATTER_CELL) +
              MARGIN_CELLS;
    int cj0 = (int)floorf((float)run->win.y0 / WG_SCATTER_CELL) - MARGIN_CELLS;
    int cj1 = (int)floorf((float)(run->win.y0 + run->win.h) / WG_SCATTER_CELL) +
              MARGIN_CELLS;
    if (ci0 < 0) ci0 = 0;
    if (cj0 < 0) cj0 = 0;
    if (ci1 > ncells - 1) ci1 = ncells - 1;
    if (cj1 > ncells - 1) cj1 = ncells - 1;

    int ncw = ci1 - ci0 + 1, nch = cj1 - cj0 + 1;
    Candidate *cand = (Candidate *)calloc((size_t)ncw * (size_t)nch,
                                          sizeof(Candidate));
    if (!cand) {
        sp_grid_free(&region);
        return 0;
    }
    ScatterCtx sc = {cand, ci0, cj0, ncw, nch, &region,
                     run->sub[SP_LAYER_SCATTER], region.ox, region.oy};
    sp_parallel_rows(nch, run->threads, prepare_rows, &sc);
    sp_parallel_rows(nch, run->threads, accept_rows, &sc);

    int cap = ncw * nch, n = 0;
    SpPoint *pts = (SpPoint *)malloc((size_t)cap * sizeof(SpPoint));
    if (!pts) {
        free(cand);
        sp_grid_free(&region);
        return 0;
    }
    /* Sequential gather in (cj, ci) order: the list itself is deterministic. */
    for (int j = 0; j < nch; j++) {
        for (int i = 0; i < ncw; i++) {
            const Candidate *c = &cand[(size_t)j * (size_t)ncw + (size_t)i];
            if (!c->accepted) continue;
            if (c->x < (float)run->win.x0 ||
                c->x >= (float)(run->win.x0 + run->win.w)) continue;
            if (c->y < (float)run->win.y0 ||
                c->y >= (float)(run->win.y0 + run->win.h)) continue;
            pts[n].x = c->x;
            pts[n].y = c->y;
            n++;
        }
    }
    out->p = pts;
    out->n = n;
    free(cand);
    sp_grid_free(&region);
    return 1;
}
