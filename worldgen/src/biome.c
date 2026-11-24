#include <stdlib.h>

#include "parallel.h"
#include "rng.h"
#include "worldgen.h"

/* Biome: classify each cell from the eroded height, local slope, and two
 * climate fields.  Both climate fields are fBm drawn from the biome sub-seed
 * at global coordinates, so neighbouring tiles agree at their shared edge. */
typedef struct {
    const SpGrid *eroded;
    unsigned char *out;
    uint64_t seed;
    int ox, oy;
} BiomeCtx;

static unsigned char classify(float h, float slope, float temp, float moist) {
    if (h < WG_SEA_LEVEL) return WG_BIOME_WATER;
    if (slope > 2.2f) return WG_BIOME_ROCK;
    if (temp < -0.25f) return WG_BIOME_TUNDRA;
    if (moist < -0.22f) return WG_BIOME_DESERT;
    if (moist > 0.18f) return WG_BIOME_FOREST;
    return WG_BIOME_GRASS;
}

static void biome_rows(void *vp, int row) {
    BiomeCtx *c = (BiomeCtx *)vp;
    const SpGrid *e = c->eroded;
    float wy = (float)(c->oy + row);
    for (int lx = 0; lx < e->w; lx++) {
        float wx = (float)(c->ox + lx);
        float temp = sp_fbm(wx * 0.0060f + 13.7f, wy * 0.0060f + 91.2f,
                            4, 2.0f, 0.5f, c->seed);
        float moist = sp_fbm(wx * 0.0085f - 47.3f, wy * 0.0085f + 311.9f,
                             4, 2.0f, 0.5f,
                             sp_mix64(c->seed ^ UINT64_C(0x8D58AC5B7C3A1E27)));
        float h = sp_grid_at(e, lx, row);
        float slope = wg_slope_at(e, lx, row);
        c->out[(size_t)row * (size_t)e->w + (size_t)lx] =
            classify(h, slope, temp, moist);
    }
}

int wg_layer_biome(unsigned char *out, const SpRun *run) {
    int halo = wg_halo_for_layer(SP_LAYER_BIOME);
    int ox, oy, ow, oh;
    wg_padded_bounds(&run->win, halo, &ox, &oy, &ow, &oh);

    SpGrid region;
    if (!sp_grid_alloc(&region, ox, oy, ow, oh)) return 0;
    wg_fill_height(&region, run->sub[SP_LAYER_NOISE], run->threads);
    wg_erode(&region, run->sub[SP_LAYER_EROSION], run->threads);

    unsigned char *full = (unsigned char *)malloc((size_t)ow * (size_t)oh);
    if (!full) {
        sp_grid_free(&region);
        return 0;
    }
    BiomeCtx c = {&region, full, run->sub[SP_LAYER_BIOME], ox, oy};
    sp_parallel_rows(oh, run->threads, biome_rows, &c);

    int lx0 = run->win.x0 - ox, ly0 = run->win.y0 - oy;
    for (int y = 0; y < run->win.h; y++) {
        for (int x = 0; x < run->win.w; x++) {
            out[(size_t)y * (size_t)run->win.w + (size_t)x] =
                full[(size_t)(ly0 + y) * (size_t)ow + (size_t)(lx0 + x)];
        }
    }
    free(full);
    sp_grid_free(&region);
    return 1;
}
