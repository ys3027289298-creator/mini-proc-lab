
#include <math.h>
#include <string.h>
#include "scatter.h"

/* SplitMix-style finaliser. Candidates must be keyed by both grid coordinates:
 * the old key `i * 131 + j * 17 + seed * 7919` collided for roughly 3.5 of
 * every 4 candidates in a 512x512 grid, so jitter and acceptance draws were
 * reused and the measured density was biased. */
static unsigned h1(unsigned x) {
    x += 0x9e3779b9u;
    x ^= x >> 16;
    x *= 0x21f0aaadu;
    x ^= x >> 15;
    x *= 0x735a2d97u;
    x ^= x >> 15;
    return x;
}

static unsigned key2(int i, int j, unsigned seed) {
    return h1((unsigned)i * 0x9e3779b9u ^ h1((unsigned)j * 0x85ebca6bu ^ seed));
}

static float u01(unsigned i) { return (float)(h1(i) & 0xffff) / 65535.0f; }

void set_init(Set* s) { memset(s, 0, sizeof(*s)); }

static int set_add(Set* s, float x, float y, float mind) {
    if (s->n >= MAXP) return -1;
    for (int i = 0; i < s->n; i++) {
        float dx = s->p[i].x - x, dy = s->p[i].y - y;
        if (dx * dx + dy * dy < mind * mind) return -2;
    }
    s->p[s->n].x = x; s->p[s->n].y = y;
    return s->n++;
}

float slope_density_eval(const SlopeDensity* d, float slope) {
    if (!d) return 1.0f;
    if (d->s_cut > 0.0f && slope >= d->s_cut) return 0.0f;
    if (d->s0 <= 0.0f) return d->rho0;
    float t = powf(slope / d->s0, d->p);
    float r = d->rho0 * expf(-t);
    if (r < 0.0f) r = 0.0f;
    if (r > 1.0f) r = 1.0f;
    return r;
}

int jitter_grid(Set* s, int nx, int ny, float cell, float mind, unsigned seed,
                const SlopeDensity* dens, SlopeFn slope_fn, void* slope_ctx) {
    int rejected = 0;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            unsigned k = key2(i, j, seed);
            float x = (i + u01(k)) * cell;
            float y = (j + u01(k + 1)) * cell;
            if (dens && slope_fn) {
                float rho = slope_density_eval(dens, slope_fn(slope_ctx, x, y));
                if (u01(k + 2) >= rho) { rejected++; continue; }
            }
            int r = set_add(s, x, y, mind);
            if (r == -2) { rejected++; continue; }
            if (r == -1) return -1;
        }
    }
    return rejected;
}

int jitter_grid_flat(Set* s, int nx, int ny, float cell, float mind,
                     unsigned seed) {
    return jitter_grid(s, nx, ny, cell, mind, seed, NULL, NULL, NULL);
}

int scatter_band_counts(const Set* s, SlopeFn slope_fn, void* slope_ctx,
                        float lo, float width, int nbins, int* retained) {
    if (!s || !retained || !(width > 0.0f) || nbins <= 0) return -1;
    for (int i = 0; i < nbins; i++) retained[i] = 0;
    int outside = 0;
    for (int i = 0; i < s->n; i++) {
        float slope = slope_fn ? slope_fn(slope_ctx, s->p[i].x, s->p[i].y) : 0.0f;
        int b = (int)((slope - lo) / width);
        if (b < 0 || b >= nbins) { outside++; continue; }
        retained[b]++;
    }
    return outside;
}

int scatter_candidate_counts(int nx, int ny, float cell, unsigned seed,
                             SlopeFn slope_fn, void* slope_ctx, float lo,
                             float width, int nbins, int* exposure) {
    if (!exposure || !(width > 0.0f) || nbins <= 0) return -1;
    for (int i = 0; i < nbins; i++) exposure[i] = 0;
    int outside = 0;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            unsigned k = key2(i, j, seed);
            float x = (i + u01(k)) * cell;
            float y = (j + u01(k + 1)) * cell;
            float slope = slope_fn ? slope_fn(slope_ctx, x, y) : 0.0f;
            int b = (int)((slope - lo) / width);
            if (b < 0 || b >= nbins) { outside++; continue; }
            exposure[b]++;
        }
    }
    return outside;
}

int scatter_measure(int nx, int ny, float cell, unsigned seed, SlopeFn slope_fn,
                    void* slope_ctx, DensityFn truth, void* truth_ctx, float lo,
                    float width, int nbins, int* exposure, int* retained,
                    double* slope_sum, double* expect_sum) {
    if (!exposure || !retained || !(width > 0.0f) || nbins <= 0) return -1;
    for (int i = 0; i < nbins; i++) {
        exposure[i] = 0;
        retained[i] = 0;
        if (slope_sum) slope_sum[i] = 0.0;
        if (expect_sum) expect_sum[i] = 0.0;
    }
    int outside = 0;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            unsigned k = key2(i, j, seed);
            float x = (i + u01(k)) * cell;
            float y = (j + u01(k + 1)) * cell;
            float slope = slope_fn ? slope_fn(slope_ctx, x, y) : 0.0f;
            int b = (int)((slope - lo) / width);
            if (b < 0 || b >= nbins) { outside++; continue; }
            exposure[b]++;
            if (slope_sum) slope_sum[b] += slope;
            float rho = truth ? truth(truth_ctx, slope) : 1.0f;
            if (rho < 0.0f) rho = 0.0f;
            if (rho > 1.0f) rho = 1.0f;
            if (expect_sum) expect_sum[b] += rho;
            if (u01(k + 2) < rho) retained[b]++;
        }
    }
    return outside;
}
