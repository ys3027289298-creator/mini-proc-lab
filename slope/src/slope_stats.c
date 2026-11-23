#include <math.h>
#include <string.h>
#include "slope_stats.h"

/* Gentle / moderate / steep, rise per run. 0.5 and 1.0 correspond to 26.6 and
 * 45 degrees, the usual cut-offs for ground that still holds herbs. */
const double slope_strata_bounds[SLOPE_STRATA - 1] = { 0.50, 1.00 };

void slope_wilson(double k, double m, double alpha, double* lo, double* hi) {
    if (!(m > 0.0)) { *lo = 0.0; *hi = 1.0; return; }
    double z = slope_norm_quantile(1.0 - 0.5 * alpha);
    double p = k / m;
    double z2 = z * z;
    double denom = 1.0 + z2 / m;
    double centre = (p + z2 / (2.0 * m)) / denom;
    double half = z * sqrt((p * (1.0 - p) + z2 / (4.0 * m)) / m) / denom;
    double a = centre - half, b = centre + half;
    *lo = a < 0.0 ? 0.0 : a;
    *hi = b > 1.0 ? 1.0 : b;
}

static int stratum_of(double s) {
    for (int i = 0; i < SLOPE_STRATA - 1; i++)
        if (s < slope_strata_bounds[i]) return i;
    return SLOPE_STRATA - 1;
}

void slope_strata(const SlopeBin* bins, int n, const SlopeFit* fit, double alpha,
                  SlopeStratum out[SLOPE_STRATA]) {
    memset(out, 0, sizeof(SlopeStratum) * SLOPE_STRATA);
    for (int i = 0; i < SLOPE_STRATA; i++) {
        out[i].lo = i == 0 ? 0.0 : slope_strata_bounds[i - 1];
        out[i].hi = i == SLOPE_STRATA - 1 ? -1.0 : slope_strata_bounds[i];
    }
    for (int i = 0; i < n; i++) {
        if (!(bins[i].exposure > 0.0)) continue;
        SlopeStratum* s = &out[stratum_of(bins[i].slope)];
        s->exposure += bins[i].exposure;
        s->retained += bins[i].retained;
    }
    for (int i = 0; i < SLOPE_STRATA; i++) {
        SlopeStratum* s = &out[i];
        s->observed = s->exposure > 0.0 ? s->retained / s->exposure : 0.0;
        slope_wilson(s->retained, s->exposure, alpha, &s->obs_lo, &s->obs_hi);
        /* Exposure-weighted mean model density for bins in this layer. */
        double num = 0.0, den = 0.0;
        for (int j = 0; j < n; j++) {
            if (!(bins[j].exposure > 0.0)) continue;
            if (stratum_of(bins[j].slope) != i) continue;
            num += bins[j].exposure * slope_attenuation(fit->kind, fit->rho0,
                                                        fit->s0, fit->p,
                                                        bins[j].slope);
            den += bins[j].exposure;
        }
        s->predicted = den > 0.0 ? num / den : 0.0;
        s->resid = s->observed - s->predicted;
        /* Two-sided normal test of the layer's observed count against the
         * model prediction: z and p flag any layer the curve misfits. */
        double var = s->exposure * s->predicted * (1.0 - s->predicted);
        s->z = var > 1e-12 ? s->resid * s->exposure / sqrt(var) : 0.0;
        s->pval = 2.0 * (1.0 - 0.5 * (1.0 + erf(fabs(s->z) / sqrt(2.0))));
    }
}
