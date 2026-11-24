#ifndef SP_WORLDGEN_H
#define SP_WORLDGEN_H

#include <stdint.h>

#include "grid.h"
#include "seed.h"

/* World is square and fixed at build time.  256 is the default; the seam tests
 * build a 512 world so a window plus its halo still leaves a margin, which is
 * what proves the halo is sufficient rather than masking the issue with world
 * edge clamping. */
#ifndef WG_WORLD
#define WG_WORLD 256
#endif

#define WG_EROSION_ITERS 64
#define WG_HYDRO_ITERS 48

#define WG_SCATTER_CELL 2.0f
#define WG_SCATTER_MIND 1.10f

#define WG_SEA_LEVEL 18.0f
#define WG_TALUS 0.45f

/* Biome ids, stored in biome.u8. */
enum {
    WG_BIOME_WATER = 0,
    WG_BIOME_DESERT = 1,
    WG_BIOME_FOREST = 2,
    WG_BIOME_GRASS = 3,
    WG_BIOME_TUNDRA = 4,
    WG_BIOME_ROCK = 5
};

typedef struct {
    int x0, y0, w, h;
} SpWindow;

typedef struct {
    uint64_t master;
    uint64_t sub[SP_LAYER_COUNT];
    int threads;
    SpWindow win;
} SpRun;

typedef struct {
    float x, y;
} SpPoint;

typedef struct {
    SpPoint *p;
    int n;
} SpScatter;

/* Halo, in world cells, each layer needs around the requested window in order
 * to be bit-identical to the full 256x256 world.  Erosion transports material
 * one cell per iteration, flow accumulation likewise. */
int wg_halo_for_layer(SpLayer layer);

/* Intersects the requested window expanded by halo with the world bounds, so
 * region edges coincide with world edges. */
void wg_padded_bounds(const SpWindow *win, int halo, int *ox, int *oy,
                      int *ow, int *oh);

/* Builds the padded noise heightfield used as the base of the other layers.
 * Returns 0 on failure. */
int wg_build_height(SpGrid *r, const SpRun *run, int halo);

/* Each layer writes win.w * win.h values and is independent of every other
 * layer; it recomputes its own inputs.  Return 0 on failure. */
int wg_layer_noise(SpGrid *out, const SpRun *run);
int wg_layer_erosion(SpGrid *out, const SpRun *run);
int wg_layer_biome(unsigned char *out, const SpRun *run);
int wg_layer_hydro(SpGrid *out, const SpRun *run);
int wg_layer_scatter(SpScatter *out, const SpRun *run);

/* Shared helpers implemented in noise.c / erosion.c. */
void wg_fill_height(SpGrid *r, uint64_t noise_seed, int threads);
void wg_erode(SpGrid *r, uint64_t erosion_seed, int threads);
float wg_talus_at(int wx, int wy, uint64_t erosion_seed);
float wg_slope_at(const SpGrid *r, int lx, int ly);

#endif
