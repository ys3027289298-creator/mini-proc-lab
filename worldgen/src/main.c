/*
 * worldgen CLI: one master seed, one layer, one thread count, one window.
 *
 *   worldgen --seed N --layer noise|erosion|biome|hydro|scatter|all
 *            --threads T --out DIR [--window x0,y0,w,h]
 *
 * Outputs are raw little-endian binary; see docs/output-format.md.
 */
#define _POSIX_C_SOURCE 200809L

#include <errno.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "parallel.h"
#include "seed.h"
#include "worldgen.h"

static const char *layer_name(SpLayer l) { return sp_layer_names[l]; }

static int layer_from_name(const char *s, SpLayer *out, int *all) {
    *all = 0;
    if (strcmp(s, "all") == 0) {
        *all = 1;
        return 1;
    }
    for (int i = 0; i < SP_LAYER_COUNT; i++) {
        if (strcmp(s, sp_layer_names[i]) == 0) {
            *out = (SpLayer)i;
            return 1;
        }
    }
    return 0;
}

static int parse_u64(const char *s, uint64_t *out) {
    char *end = NULL;
    errno = 0;
    unsigned long long v = strtoull(s, &end, 0);
    if (errno || end == s || *end != '\0') return 0;
    *out = (uint64_t)v;
    return 1;
}

static int parse_int(const char *s, int *out) {
    char *end = NULL;
    errno = 0;
    long v = strtol(s, &end, 0);
    if (errno || end == s || *end != '\0' || v < 1 || v > SP_MAX_THREADS) {
        return 0;
    }
    *out = (int)v;
    return 1;
}

static void usage(const char *argv0) {
    fprintf(stderr,
            "usage: %s --seed N --layer noise|erosion|biome|hydro|scatter|all\n"
            "          --threads T --out DIR [--window x0,y0,w,h]\n"
            "       %s --dump-seeds N [N ...]\n"
            "       %s --dump-config\n",
            argv0, argv0, argv0);
}

static void dump_seeds(int argc, char **argv) {
    for (int a = 0; a < argc; a++) {
        uint64_t master;
        if (!parse_u64(argv[a], &master)) {
            fprintf(stderr, "worldgen: bad master seed: %s\n", argv[a]);
            exit(2);
        }
        uint64_t sub[SP_LAYER_COUNT];
        sp_derive_all(master, sub);
        printf("master %" PRIu64 " 0x%016" PRIx64 "\n", master, master);
        for (int i = 0; i < SP_LAYER_COUNT; i++) {
            printf("  %-8s %20" PRIu64 "  0x%016" PRIx64 "\n",
                   sp_layer_names[i], sub[i], sub[i]);
        }
    }
}

static int ensure_dir(const char *path) {
    char buf[1024];
    size_t n = strlen(path);
    if (n == 0 || n >= sizeof(buf)) return 0;
    memcpy(buf, path, n + 1);
    for (size_t i = 1; i <= n; i++) {
        if (buf[i] == '/' || buf[i] == '\0') {
            char saved = buf[i];
            buf[i] = '\0';
            if (mkdir(buf, 0777) != 0 && errno != EEXIST) return 0;
            buf[i] = saved;
        }
    }
    return 1;
}

static int write_bytes(const char *dir, const char *name, const void *data,
                       size_t n) {
    char path[1024];
    snprintf(path, sizeof(path), "%s/%s", dir, name);
    FILE *f = fopen(path, "wb");
    if (!f) {
        fprintf(stderr, "worldgen: cannot write %s: %s\n", path, strerror(errno));
        return 0;
    }
    size_t got = fwrite(data, 1, n, f);
    int ok = (got == n) && (fclose(f) == 0);
    if (!ok) fprintf(stderr, "worldgen: short write to %s\n", path);
    return ok;
}

static int write_meta(const char *dir, const SpRun *run, SpLayer layer,
                      int all) {
    char path[1024];
    snprintf(path, sizeof(path), "%s/meta.txt", dir);
    FILE *f = fopen(path, "w");
    if (!f) return 0;
    fprintf(f, "format mini-proc-lab-worldgen-1\n");
    fprintf(f, "world %d\n", WG_WORLD);
    fprintf(f, "master %" PRIu64 " 0x%016" PRIx64 "\n", run->master, run->master);
    for (int i = 0; i < SP_LAYER_COUNT; i++) {
        fprintf(f, "subseed_%s %" PRIu64 " 0x%016" PRIx64 "\n",
                sp_layer_names[i], run->sub[i], run->sub[i]);
    }
    fprintf(f, "layer %s\n", all ? "all" : layer_name(layer));
    fprintf(f, "threads %d\n", run->threads);
    fprintf(f, "window %d,%d,%d,%d\n", run->win.x0, run->win.y0,
            run->win.w, run->win.h);
    fprintf(f, "erosion_iters %d\n", WG_EROSION_ITERS);
    fprintf(f, "hydro_iters %d\n", WG_HYDRO_ITERS);
    fprintf(f, "scatter_cell %.3f\n", (double)WG_SCATTER_CELL);
    fprintf(f, "scatter_mind %.3f\n", (double)WG_SCATTER_MIND);
    fprintf(f, "sea_level %.3f\n", (double)WG_SEA_LEVEL);
    fprintf(f, "talus %.3f\n", (double)WG_TALUS);
    int ok = (fclose(f) == 0);
    return ok;
}

int main(int argc, char **argv) {
    if (argc >= 2 && strcmp(argv[1], "--dump-seeds") == 0) {
        if (argc < 3) {
            usage(argv[0]);
            return 2;
        }
        dump_seeds(argc - 2, argv + 2);
        return 0;
    }
    if (argc >= 2 && strcmp(argv[1], "--dump-config") == 0) {
        printf("world %d\n", WG_WORLD);
        printf("erosion_iters %d\n", WG_EROSION_ITERS);
        printf("hydro_iters %d\n", WG_HYDRO_ITERS);
        printf("scatter_cell %.3f\n", (double)WG_SCATTER_CELL);
        printf("scatter_mind %.3f\n", (double)WG_SCATTER_MIND);
        printf("sea_level %.3f\n", (double)WG_SEA_LEVEL);
        printf("talus %.3f\n", (double)WG_TALUS);
        for (int i = 0; i < SP_LAYER_COUNT; i++) {
            printf("halo_%s %d\n", sp_layer_names[i],
                   wg_halo_for_layer((SpLayer)i));
        }
        return 0;
    }

    SpRun run;
    memset(&run, 0, sizeof(run));
    run.master = UINT64_C(0);
    run.threads = 1;
    run.win.x0 = 0;
    run.win.y0 = 0;
    run.win.w = WG_WORLD;
    run.win.h = WG_WORLD;

    SpLayer layer = SP_LAYER_NOISE;
    int all = 1;
    const char *outdir = NULL;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--seed") == 0 && i + 1 < argc) {
            if (!parse_u64(argv[++i], &run.master)) { usage(argv[0]); return 2; }
        } else if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
            if (!parse_int(argv[++i], &run.threads)) { usage(argv[0]); return 2; }
        } else if (strcmp(argv[i], "--layer") == 0 && i + 1 < argc) {
            if (!layer_from_name(argv[++i], &layer, &all)) {
                usage(argv[0]);
                return 2;
            }
        } else if (strcmp(argv[i], "--out") == 0 && i + 1 < argc) {
            outdir = argv[++i];
        } else if (strcmp(argv[i], "--window") == 0 && i + 1 < argc) {
            int v[4];
            if (sscanf(argv[++i], "%d,%d,%d,%d", &v[0], &v[1], &v[2], &v[3]) != 4) {
                usage(argv[0]);
                return 2;
            }
            run.win.x0 = v[0];
            run.win.y0 = v[1];
            run.win.w = v[2];
            run.win.h = v[3];
        } else if (strcmp(argv[i], "--help") == 0) {
            usage(argv[0]);
            return 0;
        } else {
            usage(argv[0]);
            return 2;
        }
    }

    if (!outdir) {
        usage(argv[0]);
        return 2;
    }
    if (run.win.w <= 0 || run.win.h <= 0 ||
        run.win.x0 < 0 || run.win.y0 < 0 ||
        run.win.x0 + run.win.w > WG_WORLD || run.win.y0 + run.win.h > WG_WORLD) {
        fprintf(stderr, "worldgen: window outside 0..%d\n", WG_WORLD);
        return 2;
    }
    if (!ensure_dir(outdir)) {
        fprintf(stderr, "worldgen: cannot create %s\n", outdir);
        return 1;
    }

    sp_derive_all(run.master, run.sub);
    if (!write_meta(outdir, &run, layer, all)) return 1;

    size_t cells = (size_t)run.win.w * (size_t)run.win.h;
    int ok = 1;

    SpGrid grid;
    memset(&grid, 0, sizeof(grid));
    unsigned char *biome = NULL;

    if (all || layer == SP_LAYER_NOISE) {
        ok &= wg_layer_noise(&grid, &run);
        if (ok) ok &= write_bytes(outdir, "noise.f32", grid.v,
                                  cells * sizeof(float));
        sp_grid_free(&grid);
    }
    if (all || layer == SP_LAYER_EROSION) {
        ok &= wg_layer_erosion(&grid, &run);
        if (ok) ok &= write_bytes(outdir, "erosion.f32", grid.v,
                                  cells * sizeof(float));
        sp_grid_free(&grid);
    }
    if (all || layer == SP_LAYER_BIOME) {
        biome = (unsigned char *)malloc(cells);
        ok &= biome != NULL && wg_layer_biome(biome, &run);
        if (ok) ok &= write_bytes(outdir, "biome.u8", biome, cells);
        free(biome);
    }
    if (all || layer == SP_LAYER_HYDRO) {
        ok &= wg_layer_hydro(&grid, &run);
        if (ok) ok &= write_bytes(outdir, "hydro.f32", grid.v,
                                  cells * sizeof(float));
        sp_grid_free(&grid);
    }
    if (all || layer == SP_LAYER_SCATTER) {
        SpScatter sc;
        memset(&sc, 0, sizeof(sc));
        ok &= wg_layer_scatter(&sc, &run);
        if (ok) {
            size_t bytes = sizeof(uint32_t) + (size_t)sc.n * 2 * sizeof(float);
            unsigned char *buf = (unsigned char *)malloc(bytes);
            if (!buf) {
                ok = 0;
            } else {
                uint32_t n32 = (uint32_t)sc.n;
                memcpy(buf, &n32, sizeof(n32));
                for (int i = 0; i < sc.n; i++) {
                    memcpy(buf + sizeof(n32) + (size_t)i * 2 * sizeof(float),
                           &sc.p[i].x, sizeof(float));
                    memcpy(buf + sizeof(n32) + (size_t)i * 2 * sizeof(float) +
                               sizeof(float),
                           &sc.p[i].y, sizeof(float));
                }
                ok &= write_bytes(outdir, "scatter.f32", buf, bytes);
                free(buf);
            }
            free(sc.p);
        }
    }

    if (!ok) {
        fprintf(stderr, "worldgen: failed\n");
        return 1;
    }
    printf("ok layer=%s master=%" PRIu64 " threads=%d window=%d,%d,%d,%d\n",
           all ? "all" : layer_name(layer), run.master, run.threads,
           run.win.x0, run.win.y0, run.win.w, run.win.h);
    return 0;
}
