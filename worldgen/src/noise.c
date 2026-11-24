#include <math.h>

#include "parallel.h"
#include "rng.h"
#include "worldgen.h"

/* Base heightfield.  Evaluated at global world coordinates, so any window of
 * the field is a subset of the full world rather than a separate noise field.
 * Only the noise sub-seed is used; no other layer's seed leaks in here. */
static float height_at(float wx, float wy, uint64_t seed) {
    float h = 20.0f + 32.0f * sp_fbm(wx * 0.010f, wy * 0.010f, 6, 2.0f, 0.5f, seed);
    h += 5.0f * sp_fbm(wx * 0.045f + 17.3f, wy * 0.045f - 9.1f, 3, 2.0f, 0.5f,
                      sp_mix64(seed ^ UINT64_C(0x51ED270B5F31A6D1)));
    return h;
}

typedef struct {
    SpGrid *r;
    uint64_t seed;
} HeightCtx;

static void height_rows(void *vp, int row) {
    HeightCtx *c = (HeightCtx *)vp;
    SpGrid *r = c->r;
    float wy = (float)(r->oy + row);
    for (int lx = 0; lx < r->w; lx++) {
        r->v[(size_t)row * (size_t)r->w + (size_t)lx] =
            height_at((float)(r->ox + lx), wy, c->seed);
    }
}

void wg_fill_height(SpGrid *r, uint64_t noise_seed, int threads) {
    HeightCtx c = {r, noise_seed};
    sp_parallel_rows(r->h, threads, height_rows, &c);
}

int wg_build_height(SpGrid *r, const SpRun *run, int halo) {
    int ox, oy, ow, oh;
    wg_padded_bounds(&run->win, halo, &ox, &oy, &ow, &oh);
    if (!sp_grid_alloc(r, ox, oy, ow, oh)) return 0;
    wg_fill_height(r, run->sub[SP_LAYER_NOISE], run->threads);
    return 1;
}

int wg_layer_noise(SpGrid *out, const SpRun *run) {
    if (!sp_grid_alloc(out, 0, 0, run->win.w, run->win.h)) return 0;
    SpGrid r;
    if (!wg_build_height(&r, run, 0)) {
        sp_grid_free(out);
        return 0;
    }
    for (int y = 0; y < run->win.h; y++) {
        for (int x = 0; x < run->win.w; x++) {
            out->v[(size_t)y * (size_t)run->win.w + (size_t)x] =
                sp_grid_at(&r, x, y);
        }
    }
    sp_grid_free(&r);
    return 1;
}
