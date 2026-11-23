/* scatter: deterministic jittered scatter with minimum distance rejection and
 * slope-dependent density attenuation. */
#ifndef SCATTER_H
#define SCATTER_H

#define MAXP 4096

typedef struct { float x, y; } P;

typedef struct { P p[MAXP]; int n; } Set;

/* Density-vs-slope model, `slope` in height units per world unit:
 *
 *     rho(s) = rho0 * exp(-(s / s0)^p)          (generalised Gaussian)
 *
 * `p == 1` is the exponential form, `p == 2` a Rayleigh/Gaussian-like falloff,
 * `p -> 0` approaches a step at s0. `s_cut` > 0 also hard-zeroes density at or
 * beyond that slope. */
typedef struct {
    float rho0;  /* density at zero slope, relative to the flat-grid density */
    float s0;    /* slope scale, same units as `slope` */
    float p;     /* shape exponent, > 0 */
    float s_cut; /* hard cutoff slope, <= 0 disables */
} SlopeDensity;

void set_init(Set* s);

/* Relative acceptance probability of a candidate at the given slope. */
float slope_density_eval(const SlopeDensity* d, float slope);

/* Jittered grid scatter. `slope_fn(slope_ctx, x, y)` returns the slope at a
 * world position; pass NULL for a flat field. Each candidate is retained with
 * probability `slope_density_eval`, sampled deterministically from the
 * candidate's hash, so the result is reproducible and the realised density
 * tracks the model. Returns the number of rejected candidates (spacing +
 * slope), or -1 if the set overflowed. */
typedef float (*SlopeFn)(void* ctx, float x, float y);

int jitter_grid(Set* s, int nx, int ny, float cell, float mind, unsigned seed,
                const SlopeDensity* dens, SlopeFn slope_fn, void* slope_ctx);

/* Convenience wrappers used by tests and the fit tool. */
int jitter_grid_flat(Set* s, int nx, int ny, float cell, float mind,
                     unsigned seed);

/* Apply the density model to an existing point set, bucketing retained plants
 * by their slope. Fills `retained[i]` with points whose slope falls in bucket
 * i. Returns the number of points outside the histogram range. */
int scatter_band_counts(const Set* s, SlopeFn slope_fn, void* slope_ctx,
                        float lo, float width, int nbins, int* retained);

/* Re-enumerate the candidate positions jitter_grid would generate and bucket
 * them by slope into `exposure`. Returns the count outside the range. */
int scatter_candidate_counts(int nx, int ny, float cell, unsigned seed,
                             SlopeFn slope_fn, void* slope_ctx, float lo,
                             float width, int nbins, int* exposure);

/* Retention probability as a function of slope, in [0, 1]. */
typedef float (*DensityFn)(void* ctx, float slope);

/* Simulated survey: offer every jittered candidate, keep it with probability
 * `truth(slope)`, and tally candidates/retentions per slope bucket. Unlike
 * jitter_grid nothing is stored, so the sample can exceed MAXP. If `slope_sum`
 * is non-NULL it accumulates the sum of candidate slopes per bucket so the
 * caller can use the bucket mean instead of the bucket centre; if
 * `expect_sum` is non-NULL it accumulates the sum of `truth(slope)` per bucket,
 * the exact expected count under the truth curve. Returns the number of
 * candidates outside the histogram range. */
int scatter_measure(int nx, int ny, float cell, unsigned seed, SlopeFn slope_fn,
                    void* slope_ctx, DensityFn truth, void* truth_ctx, float lo,
                    float width, int nbins, int* exposure, int* retained,
                    double* slope_sum, double* expect_sum);

#endif
