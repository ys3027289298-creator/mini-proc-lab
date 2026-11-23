/* slope_stats: fixed-threshold slope strata with observed/predicted density
 * and confidence intervals. Used by the fit report. */
#ifndef SLOPE_STATS_H
#define SLOPE_STATS_H

#include "slope_fit.h"

#define SLOPE_STRATA 3

/* Stratum boundaries in rise/run: [0, b0), [b0, b1), [b1, inf). */
extern const double slope_strata_bounds[SLOPE_STRATA - 1];

typedef struct {
    double lo, hi;      /* slope range of the layer, hi < 0 means open */
    double exposure;    /* candidate positions in the layer */
    double retained;    /* retained plants in the layer */
    double observed;    /* retained / exposure */
    double obs_lo, obs_hi; /* Wilson score interval at the report level */
    double predicted;   /* exposure-weighted model retention */
    double z;           /* observed vs predicted, normal score       */
    double pval;        /* two-sided p-value for that score          */
    double resid;       /* observed - predicted, in proportion */
} SlopeStratum;

/* Aggregate bins into the three fixed layers. `alpha` sets the Wilson level. */
void slope_strata(const SlopeBin* bins, int n, const SlopeFit* fit, double alpha,
                  SlopeStratum out[SLOPE_STRATA]);

/* Wilson score interval for k successes out of m trials. */
void slope_wilson(double k, double m, double alpha, double* lo, double* hi);

#endif
