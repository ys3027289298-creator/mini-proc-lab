/* slope_fit: slope -> scatter-density attenuation model, binomial maximum
 * likelihood fit with confidence intervals, and deployment gating.
 *
 * Model (derivation in slope/README.md):
 *
 *     rho(s) = rho0 * exp(-(s / s0)^p)          stretched exponential
 *
 * with nested exponential (p = 1) and an alternative Hill/logistic form
 *
 *     rho(s) = rho0 / (1 + (s / s0)^p)
 *
 * for model comparison. `s` is slope magnitude (rise/run), `rho0` the retention
 * probability on flat ground relative to the candidate density, `s0` the slope
 * scale and `p` the shape exponent.
 */
#ifndef SLOPE_FIT_H
#define SLOPE_FIT_H

#include <stddef.h>

typedef enum {
    SLOPE_MODEL_EXP = 0,       /* rho0 * exp(-s / s0)              */
    SLOPE_MODEL_STRETCHED = 1, /* rho0 * exp(-(s / s0)^p)          */
    SLOPE_MODEL_HILL = 2       /* rho0 / (1 + (s / s0)^p)          */
} SlopeModelKind;

/* One slope stratum of a survey: `exposure` candidate positions were tested,
 * `retained` of them produced a plant. */
typedef struct {
    double slope;    /* representative slope of the stratum, rise/run */
    double exposure; /* number of candidate positions              */
    double retained; /* number of retained plants                  */
} SlopeBin;

typedef struct {
    SlopeModelKind kind;
    double rho0, s0, p;    /* maximum-likelihood estimates, natural scale */
    double se[3];          /* standard errors, natural scale              */
    double ci_lo[3];       /* 100*(1-alpha)% confidence interval, natural */
    double ci_hi[3];       /* scale; order is rho0, s0, p                 */
    int nparam;
    int npoints;
    int df;
    int converged;
    double nll;
    double aic;
    double aicc;
    double pearson;        /* Pearson chi-square                       */
    double pearson_per_df; /* overdispersion diagnostic                */
} SlopeFit;

/* Attenuation factor in [0, 1]. For SLOPE_MODEL_EXP, `p` is ignored. */
double slope_attenuation(SlopeModelKind kind, double rho0, double s0, double p,
                         double s);

/* Parameter count (2 for exp, 3 otherwise) and a short name. */
int slope_model_nparam(SlopeModelKind kind);
const char* slope_model_name(SlopeModelKind kind);

/* Binomial maximum-likelihood fit. `alpha` sets the two-sided CI level (use
 * 0.05 for 95%). Returns 0 on success. */
int slope_fit(SlopeModelKind kind, const SlopeBin* bins, int n, double alpha,
              SlopeFit* out);

/* Fit all three families and report the AICc-best one. */
int slope_fit_best(const SlopeBin* bins, int n, double alpha, SlopeFit* out,
                   SlopeFit all[3]);

/* Expected retention for a set of bins under a fitted model, exposure-weighted:
 * sum(m_i * rho(s_i)) / sum(m_i). */
double slope_predict_band(const SlopeFit* fit, const SlopeBin* bins, int n);

/* Check deployment parameters against the fitted CI. Returns 1 if every
 * parameter lies inside its interval, 0 otherwise (reason in `msg`, if given).
 * The comparison is inclusive on both ends. */
int slope_deploy_ok(const SlopeFit* fit, double rho0, double s0, double p,
                    char* msg, size_t msglen);

/* Quantiles used by the statistics above (also handy for band intervals). */
double slope_norm_quantile(double p);
double slope_t_quantile(double p, int df);

#endif
