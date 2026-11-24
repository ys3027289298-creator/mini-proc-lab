#include "analysis.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CELL_WORLD_DEFAULT 1.0f

void config_default(Config *cfg) {
    memset(cfg, 0, sizeof(*cfg));
    cfg->cols = ROCKS_MAX_COLS;
    cfg->rows = ROCKS_MAX_ROWS;
    cfg->cell_px = 20;
    cfg->img_w = cfg->cols * cfg->cell_px;
    cfg->img_h = cfg->rows * cfg->cell_px;
    cfg->world_w = (float)cfg->cols;
    cfg->world_h = (float)cfg->rows;
    cfg->layout_seed = 1337u;
    cfg->n_seeds = 8;
    cfg->display_seed = 7u;
    cfg->span_limit = -1.0f; /* auto */
    cfg->deform_scale = 0.08f;
    cfg->drape_k = 0.25f;
    cfg->bump_auto = 1;
    cfg->bump.cx = 6.1f;
    cfg->bump.cy = 3.4f;
    cfg->bump.sigma = 1.00f;
    cfg->bump.amp = 0.55f;
    config_set_seed_count(cfg, cfg->n_seeds);
}

void config_set_seed_count(Config *cfg, int n) {
    static const unsigned defaults[ROCKS_MAX_SEEDS] = {
        7u, 19u, 23u, 31u, 47u, 61u, 83u, 97u, 113u, 131u, 149u, 167u,
        181u, 199u, 211u, 227u, 241u, 257u, 271u, 283u, 307u, 331u, 349u, 367u,
        389u, 401u, 419u, 431u, 443u, 461u, 479u, 491u,
    };
    if (n < 2) n = 2;
    if (n > ROCKS_MAX_SEEDS) n = ROCKS_MAX_SEEDS;
    cfg->n_seeds = n;
    for (int i = 0; i < n; i++) cfg->seeds[i] = defaults[i];
    cfg->display_seed = cfg->seeds[0];
}

int config_place_bump_best_site(Config *cfg, const Deviation *dev) {
    int best = -1;
    float best_v = -1.0f;
    for (int i = 0; i < cfg->cols * cfg->rows; i++) {
        if (dev->site_span[i] > best_v) { best_v = dev->site_span[i]; best = i; }
    }
    if (best < 0) return -1;
    int gx = best % ROCKS_MAX_COLS, gy = best / ROCKS_MAX_COLS;
    /* exactly on the site that deviated worst across seeds; the Gaussian may
     * run off the chart edge, which is a legitimate chart window, not an
     * error - the summary reports the centre so the clipping is explicit. */
    cfg->bump.cx = ((float)gx + 0.5f) * (cfg->world_w / (float)cfg->cols);
    cfg->bump.cy = ((float)gy + 0.5f) * (cfg->world_h / (float)cfg->rows);
    return best;
}

void config_print(const Config *cfg) {
    printf("chart=%dx%d cell=%d world=%.2fx%.2f layout_seed=%u\n", cfg->cols, cfg->rows,
           cfg->cell_px, cfg->world_w, cfg->world_h, cfg->layout_seed);
    printf("seeds=%d display=%u first=[", cfg->n_seeds, cfg->display_seed);
    for (int i = 0; i < cfg->n_seeds && i < 8; i++) printf("%s%u", i ? "," : "", cfg->seeds[i]);
    printf("%s]\n", cfg->n_seeds > 8 ? ",..." : "");
    printf("span_limit=%.3f deform_scale=%.3f drape_k=%.3f bump=(%.2f,%.2f sigma=%.2f amp=%.2f)\n",
           cfg->span_limit, cfg->deform_scale, cfg->drape_k, cfg->bump.cx, cfg->bump.cy,
           cfg->bump.sigma, cfg->bump.amp);
}

float lum_sample(const StoneSet *set, int module, float x, float y, unsigned seed) {
    return module == 0 ? textureA_lum(set, x, y, seed) : textureB_lum(set, x, y, seed);
}

/* ------------------------------------------------------------------ */
/* deviation across seeds                                              */
/* ------------------------------------------------------------------ */

void deviation_analyse(const Config *cfg, const StoneSet *set, int module, Deviation *d) {
    memset(d, 0, sizeof(*d));
    const int w = cfg->img_w, h = cfg->img_h;
    const int npix = w * h;
    d->w = w;
    d->h = h;
    d->cols = cfg->cols;
    d->rows = cfg->rows;

    float *lo = (float *)malloc(sizeof(float) * (size_t)npix);
    float *hi = (float *)malloc(sizeof(float) * (size_t)npix);
    d->span_px = (float *)calloc((size_t)npix, sizeof(float));
    d->span_smooth = (float *)calloc((size_t)npix, sizeof(float));
    d->site_span = (float *)calloc((size_t)ROCKS_MAX_SITES, sizeof(float));
    d->site_absurd = (unsigned char *)calloc(ROCKS_MAX_SITES, 1);
    if (!lo || !hi || !d->span_px || !d->span_smooth || !d->site_span || !d->site_absurd) {
        free(lo); free(hi);
        return;
    }

    for (int i = 0; i < npix; i++) { lo[i] = 1e9f; hi[i] = -1e9f; }

    for (int s = 0; s < cfg->n_seeds; s++) {
        unsigned seed = cfg->seeds[s];
        for (int y = 0; y < h; y++) {
            float v = ((float)y + 0.5f) / (float)h;
            for (int x = 0; x < w; x++) {
                float u = ((float)x + 0.5f) / (float)w;
                float lum = lum_sample(set, module, u * cfg->world_w, v * cfg->world_h, seed);
                int i = y * w + x;
                if (lum < lo[i]) lo[i] = lum;
                if (lum > hi[i]) hi[i] = lum;
            }
        }
    }

    double sum = 0.0, sum2 = 0.0;
    float maxspan = 0.0f;
    for (int i = 0; i < npix; i++) {
        float span = hi[i] - lo[i];
        d->span_px[i] = span;
        if (span > maxspan) maxspan = span;
    }
    d->span_max = maxspan;

    /* 3x3 box for display only; site statistics use the raw pixels */
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            float acc = 0.0f;
            int cnt = 0;
            for (int j = -1; j <= 1; j++) {
                for (int i = -1; i <= 1; i++) {
                    int xx = x + i, yy = y + j;
                    if (xx < 0 || yy < 0 || xx >= w || yy >= h) continue;
                    acc += d->span_px[yy * w + xx];
                    cnt++;
                }
            }
            d->span_smooth[y * w + x] = acc / (float)cnt;
        }
    }

    for (int gy = 0; gy < cfg->rows; gy++) {
        for (int gx = 0; gx < cfg->cols; gx++) {
            double ssum = 0.0;
            int cnt = 0;
            for (int y = gy * cfg->cell_px; y < (gy + 1) * cfg->cell_px; y++) {
                for (int x = gx * cfg->cell_px; x < (gx + 1) * cfg->cell_px; x++) {
                    ssum += d->span_px[y * w + x];
                    cnt++;
                }
            }
            float mean = cnt ? (float)(ssum / cnt) : 0.0f;
            d->site_span[gy * ROCKS_MAX_COLS + gx] = mean;
            sum += mean;
            sum2 += (double)mean * mean;
        }
    }
    int nsites = cfg->cols * cfg->rows;
    d->site_mean = (float)(sum / nsites);
    double var = sum2 / nsites - (double)d->site_mean * d->site_mean;
    d->site_std = (float)sqrt(var > 0.0 ? var : 0.0);

    float threshold = cfg->span_limit > 0.0f ? cfg->span_limit
                                             : d->site_mean + 2.0f * d->site_std;
    d->absurd_threshold = threshold;
    d->absurd_count = 0;
    for (int i = 0; i < nsites; i++) {
        if (d->site_span[i] > threshold) {
            d->site_absurd[i] = 1;
            d->absurd_count++;
        }
    }

    /* display domain: round the observed max up to a readable step */
    float step = 0.05f;
    float dom = ceilf(maxspan / step) * step;
    if (dom < 0.2f) dom = 0.2f;
    d->span_domain = dom;

    free(lo);
    free(hi);
}

void deviation_free(Deviation *d) {
    free(d->span_px);
    free(d->span_smooth);
    free(d->site_span);
    free(d->site_absurd);
    memset(d, 0, sizeof(*d));
}

/* ------------------------------------------------------------------ */
/* terrain bump + drape                                                */
/* ------------------------------------------------------------------ */

float terrain_height(const Config *cfg, float x, float y) {
    float dx = (x - cfg->bump.cx) / cfg->bump.sigma;
    float dy = (y - cfg->bump.cy) / cfg->bump.sigma;
    return cfg->bump.amp * expf(-0.5f * (dx * dx + dy * dy));
}

void slice_build(const Config *cfg, HeightSlice *s) {
    s->n = ROCKS_SLICE_N;
    s->cell_x = cfg->world_w / (float)(ROCKS_SLICE_N - 1);
    s->cell_y = cfg->world_h / (float)(ROCKS_SLICE_N - 1);
    for (int j = 0; j < ROCKS_SLICE_N; j++) {
        for (int i = 0; i < ROCKS_SLICE_N; i++) {
            float x = (float)i * s->cell_x;
            float y = (float)j * s->cell_y;
            s->h[j * ROCKS_SLICE_N + i] = terrain_height(cfg, x, y);
        }
    }
}

static float slice_at(const HeightSlice *s, int i, int j) {
    if (i < 0) i = 0;
    if (j < 0) j = 0;
    if (i >= s->n) i = s->n - 1;
    if (j >= s->n) j = s->n - 1;
    return s->h[j * s->n + i];
}

float slice_height(const HeightSlice *s, float x, float y) {
    float gx = x / s->cell_x, gy = y / s->cell_y;
    int i0 = (int)floorf(gx), j0 = (int)floorf(gy);
    float tx = gx - (float)i0, ty = gy - (float)j0;
    float a = slice_at(s, i0, j0), b = slice_at(s, i0 + 1, j0);
    float c = slice_at(s, i0, j0 + 1), d = slice_at(s, i0 + 1, j0 + 1);
    return (a + (b - a) * tx) + ((c + (d - c) * tx) - (a + (b - a) * tx)) * ty;
}

/* Bilinearly interpolated central-difference gradient. Using the nearest
 * node instead would make the drape map piecewise constant and paint visible
 * grid-shaped bands into the deformation field. */
void slice_grad(const HeightSlice *s, float x, float y, float *gx_out, float *gy_out) {
    float gx = x / s->cell_x, gy = y / s->cell_y;
    int i0 = (int)floorf(gx), j0 = (int)floorf(gy);
    float tx = gx - (float)i0, ty = gy - (float)j0;
    float g[2][2][2];
    for (int dj = 0; dj < 2; dj++) {
        for (int di = 0; di < 2; di++) {
            int i = i0 + di, j = j0 + dj;
            g[dj][di][0] = (slice_at(s, i + 1, j) - slice_at(s, i - 1, j)) / (2.0f * s->cell_x);
            g[dj][di][1] = (slice_at(s, i, j + 1) - slice_at(s, i, j - 1)) / (2.0f * s->cell_y);
        }
    }
    for (int c = 0; c < 2; c++) {
        float a = g[0][0][c] + (g[0][1][c] - g[0][0][c]) * tx;
        float b = g[1][0][c] + (g[1][1][c] - g[1][0][c]) * tx;
        float v = a + (b - a) * ty;
        if (c == 0) *gx_out = v; else *gy_out = v;
    }
}

/* Source lookup for a draped chart pixel.
 *
 * The chart is pulled downhill by k * grad H: the texture is anchored on the
 * heightfield, so a point that sits on a slope reads material from further
 * uphill. This is the classic first-order drape / lens model.
 *
 * The Jacobian of this map is  J = I - k * Hessian(H). On the dome top the
 * Hessian is negative-definite, so the grain is stretched there; on the
 * shoulders the two principal curvatures have opposite signs, so it is
 * compressed along one direction and stretched along the other. Both show
 * up in the deformation panel. Displacement alone (the downhill slide) does
 * not deform the grain, so it is reported separately.
 */
static void drape_src(const Config *cfg, const HeightSlice *slice, float u, float v,
                      float *su, float *sv) {
    float wx = u * cfg->world_w, wy = v * cfg->world_h;
    float gx, gy;
    slice_grad(slice, wx, wy, &gx, &gy);
    *su = (wx - cfg->drape_k * gx) / cfg->world_w;
    *sv = (wy - cfg->drape_k * gy) / cfg->world_h;
}

void render_lum_flat(const Config *cfg, const StoneSet *set, int module, unsigned seed,
                     float *out) {
    for (int y = 0; y < cfg->img_h; y++) {
        float v = ((float)y + 0.5f) / (float)cfg->img_h;
        for (int x = 0; x < cfg->img_w; x++) {
            float u = ((float)x + 0.5f) / (float)cfg->img_w;
            out[y * cfg->img_w + x] =
                lum_sample(set, module, u * cfg->world_w, v * cfg->world_h, seed);
        }
    }
}

void render_lum_draped(const Config *cfg, const HeightSlice *slice, const StoneSet *set,
                       int module, unsigned seed, float *out) {
    for (int y = 0; y < cfg->img_h; y++) {
        float v = ((float)y + 0.5f) / (float)cfg->img_h;
        for (int x = 0; x < cfg->img_w; x++) {
            float u = ((float)x + 0.5f) / (float)cfg->img_w;
            float su, sv;
            drape_src(cfg, slice, u, v, &su, &sv);
            float lum = lum_sample(set, module, su * cfg->world_w, sv * cfg->world_h, seed);
            /* relief shading from the slice itself (fixed, seed-independent) */
            float gx, gy;
            slice_grad(slice, su * cfg->world_w, sv * cfg->world_h, &gx, &gy);
            float shade = 1.0f + 0.16f * (-gx * 0.55f - gy * 0.83f);
            lum *= shade;
            out[y * cfg->img_w + x] = lum < 0.0f ? 0.0f : (lum > 1.0f ? 1.0f : lum);
        }
    }
}

void render_lum_draped_grid(const Config *cfg, const HeightSlice *slice, const StoneSet *set,
                            int module, unsigned seed, float *out) {
    for (int y = 0; y < cfg->img_h; y++) {
        float v = ((float)y + 0.5f) / (float)cfg->img_h;
        for (int x = 0; x < cfg->img_w; x++) {
            float u = ((float)x + 0.5f) / (float)cfg->img_w;
            float su, sv;
            drape_src(cfg, slice, u, v, &su, &sv);
            float wx = su * cfg->world_w, wy = sv * cfg->world_h;
            float lum = lum_sample(set, module, wx, wy, seed);
            float gx, gy;
            slice_grad(slice, wx, wy, &gx, &gy);
            float shade = 1.0f + 0.16f * (-gx * 0.55f - gy * 0.83f);
            lum *= shade;
            /* world-space grid: its spacing on screen is the local stretch.
             * 1 world unit ~ one site cell, so the deformation is legible. */
            float fx = fabsf(wx - floorf(wx + 0.5f));
            float fy = fabsf(wy - floorf(wy + 0.5f));
            if (fx < 0.012f || fy < 0.012f) lum *= 0.55f;
            if (lum < 0.0f) lum = 0.0f;
            if (lum > 1.0f) lum = 1.0f;
            out[y * cfg->img_w + x] = lum;
        }
    }
}

void drape_analyse(const Config *cfg, const HeightSlice *slice, Drape *d) {
    memset(d, 0, sizeof(*d));
    int w = cfg->img_w, h = cfg->img_h;
    int npix = w * h;
    d->disp = (float *)calloc((size_t)npix, sizeof(float));
    d->deform = (float *)calloc((size_t)npix, sizeof(float));
    d->over = (unsigned char *)calloc((size_t)npix, 1);
    d->site_deform = (float *)calloc((size_t)ROCKS_MAX_SITES, sizeof(float));
    if (!d->disp || !d->deform || !d->over || !d->site_deform) return;

    for (int y = 0; y < h; y++) {
        float v = ((float)y + 0.5f) / (float)h;
        for (int x = 0; x < w; x++) {
            float u = ((float)x + 0.5f) / (float)w;
            float su, sv;
            drape_src(cfg, slice, u, v, &su, &sv);
            float du = (su - u) * cfg->world_w;
            float dv = (sv - v) * cfg->world_h;
            float disp = sqrtf(du * du + dv * dv);

            /* local Jacobian of the drape map. stretching/compressing the
             * texture is what visibly deforms the grain; the downhill slide
             * is translation only and is reported as displacement. */
            float e = 0.5f / (float)cfg->cell_px;
            float su1, sv1, su2, sv2;
            drape_src(cfg, slice, u + e, v, &su1, &sv1);
            drape_src(cfg, slice, u, v + e, &su2, &sv2);
            float j00 = (su1 - su) / e, j10 = (sv1 - sv) / e;
            float j01 = (su2 - su) / e, j11 = (sv2 - sv) / e;
            float a = j00 * j00 + j10 * j10;
            float b = j00 * j01 + j10 * j11;
            float c = j01 * j01 + j11 * j11;
            float tr = a + c;
            float det = a * c - b * b;
            float disc = tr * tr * 0.25f - det;
            float lmax = tr * 0.5f + sqrtf(disc > 0.0f ? disc : 0.0f);
            float lmin = tr * 0.5f - sqrtf(disc > 0.0f ? disc : 0.0f);
            float smax = sqrtf(lmax > 0.0f ? lmax : 0.0f);
            float smin = sqrtf(lmin > 0.0f ? lmin : 0.0f);
            float ds = fabsf(smax - 1.0f), dn = fabsf(smin - 1.0f);
            float deform = ds > dn ? ds : dn;
            /* a vanishing/negative determinant means the drape map folded at
             * this pixel: the singular-value formula is no longer meaningful
             * and its square root blows up. That is a fold, not a gradient -
             * fall back to the displacement term alone. */
            if (!(det > 0.0f) || !(deform < 1e3f)) deform = 0.0f;

            int i = y * w + x;
            d->disp[i] = disp;
            d->deform[i] = deform;
            if (deform > cfg->deform_scale) d->over[i] = 1;
            if (deform > d->deform_max) d->deform_max = deform;
            if (disp > d->disp_max) d->disp_max = disp;
        }
    }

    /* per-site max, only sites whose *max* is over the threshold get a mark */
    for (int gy = 0; gy < cfg->rows; gy++) {
        for (int gx = 0; gx < cfg->cols; gx++) {
            float mx = 0.0f;
            for (int y = gy * cfg->cell_px; y < (gy + 1) * cfg->cell_px; y++)
                for (int x = gx * cfg->cell_px; x < (gx + 1) * cfg->cell_px; x++) {
                    float vv = d->deform[y * w + x];
                    if (vv > mx) mx = vv;
                }
            d->site_deform[gy * ROCKS_MAX_COLS + gx] = mx;
        }
    }

    /* collect over-threshold sites, strongest first, capped */
    typedef struct { float v; int idx; } Cand;
    Cand cand[ROCKS_MAX_SITES];
    int nc = 0;
    for (int i = 0; i < cfg->cols * cfg->rows; i++) {
        if (d->site_deform[i] > cfg->deform_scale) {
            cand[nc].v = d->site_deform[i];
            cand[nc].idx = i;
            nc++;
        }
    }
    for (int i = 1; i < nc; i++) {
        Cand key = cand[i];
        int j = i - 1;
        while (j >= 0 && cand[j].v < key.v) { cand[j + 1] = cand[j]; j--; }
        cand[j + 1] = key;
    }
    d->cross_count = nc < ROCKS_MAX_MARKS ? nc : ROCKS_MAX_MARKS;
    for (int i = 0; i < d->cross_count; i++) d->cross_site[i] = cand[i].idx;
}

void drape_free(Drape *d) {
    free(d->disp);
    free(d->deform);
    free(d->over);
    free(d->site_deform);
    memset(d, 0, sizeof(*d));
}

/* ------------------------------------------------------------------ */
/* reproducibility helpers                                             */
/* ------------------------------------------------------------------ */

uint32_t hash_floats(const float *v, int n) {
    uint32_t h = 2166136261u;
    const unsigned char *p = (const unsigned char *)v;
    for (size_t i = 0; i < (size_t)n * sizeof(float); i++) {
        h ^= p[i];
        h *= 16777619u;
    }
    return h;
}

float pearson_subsampled(const float *a, const float *b, int w, int h, int step) {
    if (step < 1) step = 1;
    double sa = 0, sb = 0, saa = 0, sbb = 0, sab = 0;
    long n = 0;
    for (int y = 0; y < h; y += step) {
        for (int x = 0; x < w; x += step) {
            double va = a[y * w + x], vb = b[y * w + x];
            sa += va; sb += vb; saa += va * va; sbb += vb * vb; sab += va * vb;
            n++;
        }
    }
    if (n < 2) return 0.0f;
    double cov = sab / n - (sa / n) * (sb / n);
    double va = saa / n - (sa / n) * (sa / n);
    double vb = sbb / n - (sb / n) * (sb / n);
    if (va <= 1e-12 || vb <= 1e-12) return 0.0f;
    return (float)(cov / sqrt(va * vb));
}
