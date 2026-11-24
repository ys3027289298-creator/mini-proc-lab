#ifndef ROCKS_ANALYSIS_H
#define ROCKS_ANALYSIS_H

#include <stdint.h>

#include "stones.h"

#define ROCKS_MAX_SEEDS 32
#define ROCKS_SLICE_N   33
#define ROCKS_MAX_COLS  12
#define ROCKS_MAX_ROWS  8
#define ROCKS_MAX_SITES (ROCKS_MAX_COLS * ROCKS_MAX_ROWS)
#define ROCKS_MAX_MARKS 32

typedef struct {
    float cx, cy, sigma, amp; /* heightfield bump, world units */
} Bump;

/* Discrete height slice shared by both texture modules. The chart is draped
 * over *this buffer*, not over the analytic bump, so "same terrain slice"
 * is a real shared input whose identity can be hashed. */
typedef struct {
    int n;
    float h[ROCKS_SLICE_N * ROCKS_SLICE_N];
    float cell_x, cell_y;
} HeightSlice;

typedef struct {
    int img_w, img_h;
    int cols, rows, cell_px;
    float world_w, world_h;
    unsigned layout_seed;
    int n_seeds;
    unsigned seeds[ROCKS_MAX_SEEDS];
    unsigned display_seed;
    float span_limit;   /* <= 0: auto (mean + 2 sigma over sites) */
    float deform_scale; /* local-stretch threshold (fraction, e.g. 0.08 = 8%) */
    float drape_k;      /* drape strength, multiplies height gradient */
    int bump_auto;      /* 1 = place the bump on the worst deviation site */
    Bump bump;
} Config;

typedef struct {
    int w, h;
    int cols, rows;
    float *span_px;     /* raw max-min luminance across seeds, per pixel */
    float *span_smooth; /* 3x3 box for display */
    float *site_span;   /* mean raw span per site */
    unsigned char *site_absurd;
    int absurd_count;
    float site_mean, site_std;
    float absurd_threshold;
    float span_domain;  /* display scale */
    float span_max;
} Deviation;

typedef struct {
    float *disp;   /* |texture displacement| per pixel, world units */
    float *deform; /* local stretch of the drape map: max(|smax-1|,|smin-1|) */
    unsigned char *over;
    float *site_deform; /* max deformation in each site */
    int cross_count;    /* sites marked (<= MAX_MARKS) */
    int cross_site[ROCKS_MAX_MARKS];
    float deform_max;
    float disp_max;
} Drape;

void config_default(Config *cfg);
void config_print(const Config *cfg);

/* Fill the configured default seed list (display seed is seeds[0]). */
void config_set_seed_count(Config *cfg, int n);

/* Luminance of one texture module (0 = value fBm, 1 = gradient fBm). */
float lum_sample(const StoneSet *set, int module, float x, float y, unsigned seed);

/* Cross-seed span analysis for one module over the whole chart. */
void deviation_analyse(const Config *cfg, const StoneSet *set, int module, Deviation *out);
void deviation_free(Deviation *d);

/* Place the bump at the centre of the site with the largest mean span.
 * Returns the site index. */
int config_place_bump_best_site(Config *cfg, const Deviation *dev);

/* Analytic bump, used only to build the slice. */
float terrain_height(const Config *cfg, float x, float y);

/* Build the shared 33x33 height slice over the chart. */
void slice_build(const Config *cfg, HeightSlice *s);
/* Bilinear height and central-difference gradient of the slice. */
float slice_height(const HeightSlice *s, float x, float y);
void slice_grad(const HeightSlice *s, float x, float y, float *gx, float *gy);

/* Drape the chart over the shared slice; computes displacement, local stretch
 * and marks. module/seed do not affect the result by construction. */
void drape_analyse(const Config *cfg, const HeightSlice *slice, Drape *out);
void drape_free(Drape *d);

/* Luminance image of the flat chart (module/seed) into out[w*h]. */
void render_lum_flat(const Config *cfg, const StoneSet *set, int module, unsigned seed,
                     float *out);
/* Luminance image draped over the shared slice, with relief shading. */
void render_lum_draped(const Config *cfg, const HeightSlice *slice, const StoneSet *set,
                       int module, unsigned seed, float *out);
/* Same, plus a warped world-space grid so local stretch is readable. */
void render_lum_draped_grid(const Config *cfg, const HeightSlice *slice, const StoneSet *set,
                            int module, unsigned seed, float *out);

uint32_t hash_floats(const float *v, int n);
float pearson_subsampled(const float *a, const float *b, int w, int h, int step);

#endif
