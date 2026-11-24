#include "worldgen.h"

/* Erosion moves material one cell per iteration, so its result at a window
 * edge is exact once the halo covers every iteration.  Biome reads slope, so
 * it needs one more cell than erosion.  Hydro runs the erosion pass and then
 * accumulates flow along paths of up to WG_HYDRO_ITERS cells, so its halo is
 * the sum plus a margin.  Scatter tests a radius of a few cells on the eroded
 * field, hence the extra margin. */
int wg_halo_for_layer(SpLayer layer) {
#ifdef SP_SELFTEST_SHORT_HALO
    /* NEGATIVE CONTROL: deliberately too small.  Window seams must then show
     * up in the outputs, proving the seam check is not vacuous. */
    switch (layer) {
        case SP_LAYER_EROSION: return 3;
        case SP_LAYER_BIOME: return 4;
        case SP_LAYER_HYDRO: return 5;
        case SP_LAYER_SCATTER: return 1;
        default: return 0;
    }
#else
    switch (layer) {
        case SP_LAYER_NOISE: return 0;
        case SP_LAYER_EROSION: return WG_EROSION_ITERS + 2;
        case SP_LAYER_BIOME: return WG_EROSION_ITERS + 3;
        case SP_LAYER_HYDRO: return WG_EROSION_ITERS + WG_HYDRO_ITERS + 4;
        case SP_LAYER_SCATTER: return WG_EROSION_ITERS + 10;
    }
    return 0;
#endif
}

void wg_padded_bounds(const SpWindow *win, int halo, int *ox, int *oy,
                      int *ow, int *oh) {
    int x0 = win->x0 - halo, y0 = win->y0 - halo;
    int x1 = win->x0 + win->w + halo, y1 = win->y0 + win->h + halo;
    if (x0 < 0) x0 = 0;
    if (y0 < 0) y0 = 0;
    if (x1 > WG_WORLD) x1 = WG_WORLD;
    if (y1 > WG_WORLD) y1 = WG_WORLD;
    *ox = x0;
    *oy = y0;
    *ow = x1 - x0;
    *oh = y1 - y0;
}
