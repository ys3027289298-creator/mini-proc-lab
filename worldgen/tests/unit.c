/* Unit tests for the worldgen layer contracts.
 *
 * These are fast structural checks; the bit-for-bit cross-run and cross-window
 * comparisons live in scripts/verify.sh.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "worldgen.h"

static int failures = 0;

static void check(const char *name, int ok) {
    printf("%-38s %s\n", name, ok ? "ok" : "FAIL");
    if (!ok) failures++;
}

static SpRun make_run(uint64_t master, int threads) {
    SpRun run;
    memset(&run, 0, sizeof(run));
    run.master = master;
    sp_derive_all(master, run.sub);
    run.threads = threads;
    run.win.x0 = 0;
    run.win.y0 = 0;
    run.win.w = WG_WORLD;
    run.win.h = WG_WORLD;
    return run;
}

static void test_seeds(void) {
    uint64_t a[SP_LAYER_COUNT], b[SP_LAYER_COUNT];
    sp_derive_all(123456789, a);
    sp_derive_all(123456789, b);
    check("sub-seeds deterministic", memcmp(a, b, sizeof(a)) == 0);

    int distinct = 1;
    for (int i = 0; i < SP_LAYER_COUNT; i++)
        for (int j = i + 1; j < SP_LAYER_COUNT; j++)
            if (a[i] == a[j]) distinct = 0;
    check("sub-seeds pairwise distinct", distinct);

    uint64_t other[SP_LAYER_COUNT];
    sp_derive_all(123456790, other);
    int moved = 1;
    for (int i = 0; i < SP_LAYER_COUNT; i++)
        if (other[i] == a[i]) moved = 0;
    check("master seed moves every layer", moved);

    check("documented noise sub-seed",
          a[SP_LAYER_NOISE] == UINT64_C(0x84E8FE37627F5BBA));
    check("documented hydro sub-seed",
          a[SP_LAYER_HYDRO] == UINT64_C(0xB69653DD63920F89));

    /* Changing the prefix/name must change the result; this guards against
     * accidental whitespace edits in the layer name table. */
    uint64_t direct = sp_subseed(123456789, "noise");
    uint64_t typo = sp_subseed(123456789, "noise ");
    check("layer name is part of the input", direct != typo);

    check("seed 0 works", sp_subseed(0, "noise") != sp_subseed(0, "biome"));
}

static void test_noise(void) {
    SpRun run = make_run(424242, 1);
    SpGrid g;
    check("noise layer builds", wg_layer_noise(&g, &run));

    int finite = 1, varied = 0;
    float lo = 1e30f, hi = -1e30f;
    for (size_t i = 0; i < (size_t)g.w * g.h; i++) {
        float v = g.v[i];
        if (!isfinite(v)) finite = 0;
        if (v < lo) lo = v;
        if (v > hi) hi = v;
    }
    varied = (hi - lo) > 1.0f;
    check("noise finite", finite);
    check("noise has relief", varied);

    /* Global coordinates: cell (0,0) of a window at (64,64) must equal cell
     * (64,64) of the full field. */
    SpRun win_run = run;
    win_run.win.x0 = 64;
    win_run.win.y0 = 64;
    win_run.win.w = 16;
    win_run.win.h = 16;
    SpGrid w;
    check("windowed noise builds", wg_layer_noise(&w, &win_run));
    int same = 1;
    for (int y = 0; y < 16; y++)
        for (int x = 0; x < 16; x++)
            if (w.v[y * 16 + x] != g.v[(64 + y) * WG_WORLD + 64 + x]) same = 0;
    check("noise is world-anchored", same);
    sp_grid_free(&w);
    sp_grid_free(&g);
}

static void test_erosion(void) {
    SpRun run = make_run(424242, 1);
    SpGrid a;
    check("erosion layer builds", wg_layer_erosion(&a, &run));

    SpGrid b;
    SpRun run4 = make_run(424242, 4);
    check("erosion builds with 4 threads", wg_layer_erosion(&b, &run4));
    check("erosion thread-invariant",
          memcmp(a.v, b.v, (size_t)a.w * a.h * sizeof(float)) == 0);
    sp_grid_free(&b);

    /* Erosion should reduce extreme local differences: compare max slope of
     * the noise layer with max slope of the eroded layer. */
    SpGrid n;
    wg_layer_noise(&n, &run);
    float slope_noise = 0.0f, slope_eroded = 0.0f;
    for (int y = 1; y < WG_WORLD - 1; y++) {
        for (int x = 1; x < WG_WORLD - 1; x++) {
            float dn = fabsf(n.v[(y + 1) * WG_WORLD + x] -
                             n.v[(y - 1) * WG_WORLD + x]);
            float de = fabsf(a.v[(y + 1) * WG_WORLD + x] -
                             a.v[(y - 1) * WG_WORLD + x]);
            if (dn > slope_noise) slope_noise = dn;
            if (de > slope_eroded) slope_eroded = de;
        }
    }
    check("erosion smooths the steepest step", slope_eroded < slope_noise);

    /* Erosion is a redistribution, not a rescale: interior totals should be
     * close, though clamping at the border leaks a little mass. */
    double sum_n = 0.0, sum_e = 0.0;
    for (int y = 8; y < WG_WORLD - 8; y++) {
        for (int x = 8; x < WG_WORLD - 8; x++) {
            sum_n += n.v[y * WG_WORLD + x];
            sum_e += a.v[y * WG_WORLD + x];
        }
    }
    double rel = fabs(sum_e - sum_n) / (fabs(sum_n) + 1e-9);
    check("erosion roughly conserves mass", rel < 0.05);
    sp_grid_free(&n);
    sp_grid_free(&a);
}

static void test_biome(void) {
    SpRun run = make_run(424242, 1);
    size_t n = (size_t)WG_WORLD * WG_WORLD;
    unsigned char *b = malloc(n);
    check("biome layer builds", b && wg_layer_biome(b, &run));

    int in_range = 1, counts[6] = {0, 0, 0, 0, 0, 0};
    for (size_t i = 0; i < n; i++) {
        if (b[i] > WG_BIOME_ROCK) in_range = 0;
        else counts[b[i]]++;
    }
    check("biome ids in range", in_range);

    int kinds = 0;
    for (int i = 0; i < 6; i++) if (counts[i] > 0) kinds++;
    check("biome is not one flat class", kinds >= 3);

    check("water exists below sea level", counts[WG_BIOME_WATER] > 0);
    check("land biomes exist", counts[WG_BIOME_GRASS] +
                                counts[WG_BIOME_FOREST] +
                                counts[WG_BIOME_DESERT] +
                                counts[WG_BIOME_TUNDRA] +
                                counts[WG_BIOME_ROCK] > 0);
    free(b);
}

static void test_hydro(void) {
    SpRun run = make_run(424242, 1);
    SpGrid a;
    check("hydro layer builds", wg_layer_hydro(&a, &run));

    SpGrid b;
    SpRun run8 = make_run(424242, 8);
    check("hydro builds with 8 threads", wg_layer_hydro(&b, &run8));
    check("hydro thread-invariant",
          memcmp(a.v, b.v, (size_t)a.w * a.h * sizeof(float)) == 0);
    sp_grid_free(&b);

    float lo = 1e30f, hi = -1e30f;
    int finite = 1, below_one = 0;
    for (size_t i = 0; i < (size_t)a.w * a.h; i++) {
        if (!isfinite(a.v[i])) finite = 0;
        if (a.v[i] < 1.0f) below_one = 1;
        if (a.v[i] < lo) lo = a.v[i];
        if (a.v[i] > hi) hi = a.v[i];
    }
    check("hydro finite", finite);
    check("hydro accumulation >= 1", !below_one);
    check("hydro has a drainage network", hi >= 8.0f);
    sp_grid_free(&a);
}

static void test_scatter(void) {
    SpRun run = make_run(424242, 1);
    SpScatter s;
    check("scatter layer builds", wg_layer_scatter(&s, &run));
    check("scatter produces points", s.n > 0);

    int min_dist_ok = 1, in_bounds = 1;
    for (int i = 0; i < s.n; i++) {
        if (s.p[i].x < 0.0f || s.p[i].x >= WG_WORLD ||
            s.p[i].y < 0.0f || s.p[i].y >= WG_WORLD) {
            in_bounds = 0;
        }
        for (int j = i + 1; j < s.n; j++) {
            float dx = s.p[i].x - s.p[j].x;
            float dy = s.p[i].y - s.p[j].y;
            if (dx * dx + dy * dy <
                WG_SCATTER_MIND * WG_SCATTER_MIND - 1e-4f) {
                min_dist_ok = 0;
            }
        }
    }
    check("scatter inside the world", in_bounds);
    check("scatter respects minimum distance", min_dist_ok);
    free(s.p);

    SpScatter s8;
    SpRun run8 = make_run(424242, 8);
    wg_layer_scatter(&s8, &run8);
    SpScatter s1;
    SpRun run1 = make_run(424242, 1);
    wg_layer_scatter(&s1, &run1);
    int same = s1.n == s8.n;
    if (same) {
        for (int i = 0; i < s1.n; i++) {
            if (s1.p[i].x != s8.p[i].x || s1.p[i].y != s8.p[i].y) same = 0;
        }
    }
    check("scatter thread-invariant", same);
    free(s1.p);
    free(s8.p);
}

static void test_window_halo(void) {
    check("halo grows with propagation",
          wg_halo_for_layer(SP_LAYER_NOISE) <
              wg_halo_for_layer(SP_LAYER_EROSION) &&
          wg_halo_for_layer(SP_LAYER_EROSION) <
              wg_halo_for_layer(SP_LAYER_HYDRO));

    SpWindow win = {8, 8, 16, 16};
    int ox, oy, ow, oh;
    wg_padded_bounds(&win, 100, &ox, &oy, &ow, &oh);
    check("padding clamps to world", ox == 0 && oy == 0 && ow == 124 && oh == 124);

    SpWindow edge = {0, 0, 4, 4};
    wg_padded_bounds(&edge, 2, &ox, &oy, &ow, &oh);
    check("window at origin stays in bounds", ox == 0 && oy == 0);
}

int main(void) {
    test_seeds();
    test_noise();
    test_erosion();
    test_biome();
    test_hydro();
    test_scatter();
    test_window_halo();
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
