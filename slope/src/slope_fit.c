#include <math.h>
#include <stdio.h>
#include <string.h>
#include "slope_fit.h"

#define TINY 1e-12

/* ------------------------------------------------------------------ */
/* special functions                                                  */
/* ------------------------------------------------------------------ */

/* Acklam's inverse normal CDF, |error| < 1.15e-9. */
double slope_norm_quantile(double p) {
    static const double a[] = { -3.969683028665376e+01, 2.209460984245205e+02,
        -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01,
        2.506628277459239e+00 };
    static const double b[] = { -5.447609879822406e+01, 1.615858368580409e+02,
        -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01 };
    static const double c[] = { -7.784894002430293e-03, -3.223964580411365e-01,
        -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00,
        2.938163982698783e+00 };
    static const double d[] = { 7.784695709041462e-03, 3.224671290700398e-01,
        2.445134137142996e+00, 3.754408661907416e+00 };
    const double plow = 0.02425, phigh = 1.0 - plow;
    double q, r;
    if (p <= 0.0) return -INFINITY;
    if (p >= 1.0) return INFINITY;
    if (p < plow) {
        q = sqrt(-2.0 * log(p));
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    }
    if (p > phigh) return -slope_norm_quantile(1.0 - p);
    q = p - 0.5;
    r = q * q;
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
           (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
}

static double betacf(double a, double b, double x) {
    const double eps = 3e-14, fpmin = 1e-300;
    double qab = a + b, qap = a + 1.0, qam = a - 1.0;
    double c = 1.0, d = 1.0 - qab * x / qap;
    if (fabs(d) < fpmin) d = fpmin;
    d = 1.0 / d;
    double h = d;
    for (int m = 1; m <= 300; m++) {
        int m2 = 2 * m;
        double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (fabs(d) < fpmin) d = fpmin;
        c = 1.0 + aa / c;
        if (fabs(c) < fpmin) c = fpmin;
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (fabs(d) < fpmin) d = fpmin;
        c = 1.0 + aa / c;
        if (fabs(c) < fpmin) c = fpmin;
        d = 1.0 / d;
        double del = d * c;
        h *= del;
        if (fabs(del - 1.0) < eps) break;
    }
    return h;
}

/* Regularised incomplete beta I_x(a, b). */
static double betai(double a, double b, double x) {
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;
    double bt = exp(lgamma(a + b) - lgamma(a) - lgamma(b) + a * log(x) +
                    b * log(1.0 - x));
    if (x < (a + 1.0) / (a + b + 2.0)) return bt * betacf(a, b, x) / a;
    return 1.0 - bt * betacf(b, a, 1.0 - x) / b;
}

/* Student-t CDF for t >= 0, integer df. */
static double t_cdf(double t, int df) {
    return 1.0 - 0.5 * betai(df / 2.0, 0.5, df / (df + t * t));
}

double slope_t_quantile(double p, int df) {
    if (df <= 0) return slope_norm_quantile(p);
    if (p <= 0.0) return -INFINITY;
    if (p >= 1.0) return INFINITY;
    if (p < 0.5) return -slope_t_quantile(1.0 - p, df);
    double lo = 0.0, hi = 1.0;
    while (t_cdf(hi, df) < p && hi < 1e12) hi *= 2.0;
    for (int i = 0; i < 200; i++) {
        double mid = 0.5 * (lo + hi);
        if (t_cdf(mid, df) < p) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

static double chi2_1_quantile(double p) {
    double z = slope_norm_quantile(0.5 + 0.5 * p);
    return z * z;
}

/* ------------------------------------------------------------------ */
/* model                                                              */
/* ------------------------------------------------------------------ */

double slope_attenuation(SlopeModelKind kind, double rho0, double s0, double p,
                         double s) {
    if (s < 0.0) s = 0.0;
    if (!(rho0 > 0.0)) return 0.0;
    if (rho0 > 1.0) rho0 = 1.0;
    if (!(s0 > 0.0)) return rho0;
    if (kind == SLOPE_MODEL_EXP) return rho0 * exp(-s / s0);
    if (!(p > 0.0)) return rho0;
    double t = pow(s / s0, p);
    if (!isfinite(t)) return 0.0;
    if (kind == SLOPE_MODEL_HILL) return rho0 / (1.0 + t);
    return rho0 * exp(-t);
}

int slope_model_nparam(SlopeModelKind kind) {
    return kind == SLOPE_MODEL_EXP ? 2 : 3;
}

const char* slope_model_name(SlopeModelKind kind) {
    switch (kind) {
    case SLOPE_MODEL_EXP: return "exp";
    case SLOPE_MODEL_HILL: return "hill";
    default: return "stretched";
    }
}

/* ------------------------------------------------------------------ */
/* binomial negative log-likelihood                                   */
/* ------------------------------------------------------------------ */

typedef struct {
    SlopeModelKind kind;
    const SlopeBin* bins;
    int n;
} NllCtx;

static double nll(const NllCtx* c, const double* prm) {
    double total = 0.0;
    for (int i = 0; i < c->n; i++) {
        double m = c->bins[i].exposure, k = c->bins[i].retained;
        if (!(m > 0.0)) continue;
        double r = slope_attenuation(c->kind, prm[0], prm[1], prm[2],
                                     c->bins[i].slope);
        if (r < TINY) r = TINY;
        if (r > 1.0 - TINY) r = 1.0 - TINY;
        if (k > 0.0) total -= k * log(r);
        if (m > k) total -= (m - k) * log1p(-r);
        if (!isfinite(total)) return INFINITY;
    }
    return total;
}

/* Optimiser works in an unconstrained space: logit(rho0), log(s0), log(p). */
static void theta_to_prm(SlopeModelKind kind, const double* th, double* prm) {
    prm[0] = 1.0 / (1.0 + exp(-th[0]));
    prm[1] = exp(th[1]);
    prm[2] = kind == SLOPE_MODEL_EXP ? 1.0 : exp(th[2]);
}

static void prm_to_theta(SlopeModelKind kind, const double* prm, double* th) {
    double r = prm[0];
    if (!(r > 1e-9)) r = 1e-9;
    if (r > 1.0 - 1e-9) r = 1.0 - 1e-9;
    th[0] = log(r / (1.0 - r));
    th[1] = log(prm[1] > 1e-9 ? prm[1] : 1e-9);
    th[2] = kind == SLOPE_MODEL_EXP ? 0.0
                                    : log(prm[2] > 1e-9 ? prm[2] : 1e-9);
}

/* Objective over free coordinates: full theta in `o->th`, free indices in
 * `o->free`. */
typedef struct {
    const NllCtx* ctx;
    double th[3];
    int free[3];
    int nfree;
} Obj;

static double obj(const double* x, void* p) {
    Obj* o = (Obj*)p;
    double th[3] = { o->th[0], o->th[1], o->th[2] };
    for (int i = 0; i < o->nfree; i++) th[o->free[i]] = x[i];
    double prm[3];
    theta_to_prm(o->ctx->kind, th, prm);
    double v = nll(o->ctx, prm);
    return isfinite(v) ? v : 1e300;
}

/* Nelder-Mead in `dim` dimensions; returns the best value, x updated. */
static double nelder_mead(Obj* o, int dim, double* x) {
    enum { MAXD = 3 };
    double s[MAXD + 1][MAXD], fs[MAXD + 1];
    const double step0 = 0.5;
    for (int i = 0; i < dim; i++) {
        memcpy(s[i], x, (size_t)dim * sizeof(double));
        s[i][i] += step0;
        fs[i] = obj(s[i], o);
    }
    memcpy(s[dim], x, (size_t)dim * sizeof(double));
    fs[dim] = obj(s[dim], o);

    for (int iter = 0; iter < 4000; iter++) {
        int ilo = 0, ihi = 0;
        for (int i = 1; i <= dim; i++) {
            if (fs[i] < fs[ilo]) ilo = i;
            if (fs[i] > fs[ihi]) ihi = i;
        }
        if (fabs(fs[ihi] - fs[ilo]) < 1e-10 * (1.0 + fabs(fs[ilo]))) break;
        double cen[MAXD] = { 0 };
        for (int i = 0; i <= dim; i++)
            if (i != ihi)
                for (int j = 0; j < dim; j++) cen[j] += s[i][j] / dim;
        double xr[MAXD], xe[MAXD], xc[MAXD], fr, fe, fc;
        for (int j = 0; j < dim; j++) xr[j] = 2.0 * cen[j] - s[ihi][j];
        fr = obj(xr, o);
        if (fr < fs[ilo]) {
            for (int j = 0; j < dim; j++) xe[j] = 3.0 * cen[j] - 2.0 * s[ihi][j];
            fe = obj(xe, o);
            if (fe < fr) { memcpy(s[ihi], xe, (size_t)dim * sizeof(double)); fs[ihi] = fe; }
            else { memcpy(s[ihi], xr, (size_t)dim * sizeof(double)); fs[ihi] = fr; }
        } else if (fr < fs[ihi]) {
            memcpy(s[ihi], xr, (size_t)dim * sizeof(double));
            fs[ihi] = fr;
        } else {
            for (int j = 0; j < dim; j++) xc[j] = 0.5 * (cen[j] + s[ihi][j]);
            fc = obj(xc, o);
            if (fc < fs[ihi]) {
                memcpy(s[ihi], xc, (size_t)dim * sizeof(double));
                fs[ihi] = fc;
            } else {
                for (int i = 0; i <= dim; i++) {
                    if (i == ilo) continue;
                    for (int j = 0; j < dim; j++)
                        s[i][j] = 0.5 * (s[i][j] + s[ilo][j]);
                    fs[i] = obj(s[i], o);
                }
            }
        }
    }
    int ib = 0;
    for (int i = 1; i <= dim; i++)
        if (fs[i] < fs[ib]) ib = i;
    memcpy(x, s[ib], (size_t)dim * sizeof(double));
    return fs[ib];
}

/* Optimise all parameters from a natural-scale start; returns nll and
 * writes the fitted natural parameters. */
static double ml_full(const NllCtx* ctx, const double* start_prm, double* out_prm) {
    int npar = slope_model_nparam(ctx->kind);
    Obj o;
    memset(&o, 0, sizeof(o));
    o.ctx = ctx;
    o.nfree = npar;
    for (int i = 0; i < npar; i++) o.free[i] = i;
    double th[3];
    prm_to_theta(ctx->kind, start_prm, th);
    for (int i = 0; i < 3; i++) o.th[i] = th[i];
    double x[3] = { th[0], th[1], th[2] };
    double f = nelder_mead(&o, npar, x);
    for (int i = 0; i < npar; i++) o.th[o.free[i]] = x[i];
    theta_to_prm(ctx->kind, o.th, out_prm);
    return f;
}

/* Profile likelihood: fix parameter `pi` at `val` (natural scale), optimise the
 * rest from `start`. Returns the constrained maximum likelihood. */
static int profile(const NllCtx* ctx, int pi, double val, const double* start_prm,
                   double* out_prm, double* out_nll) {
    int npar = slope_model_nparam(ctx->kind);
    Obj o;
    memset(&o, 0, sizeof(o));
    o.ctx = ctx;
    o.nfree = 0;
    for (int j = 0; j < npar; j++)
        if (j != pi) o.free[o.nfree++] = j;

    double start[3] = { start_prm[0], start_prm[1], start_prm[2] };
    start[pi] = val;
    double th[3];
    prm_to_theta(ctx->kind, start, th);
    for (int i = 0; i < 3; i++) o.th[i] = th[i];

    double x[3];
    for (int j = 0; j < o.nfree; j++) x[j] = th[o.free[j]];
    double f = nelder_mead(&o, o.nfree, x);
    for (int j = 0; j < o.nfree; j++) o.th[o.free[j]] = x[j];
    theta_to_prm(ctx->kind, o.th, out_prm);
    *out_nll = nll(ctx, out_prm);
    return isfinite(f) && isfinite(*out_nll);
}

/* Evaluate the profile nll when parameter `pi` is held at `val`, warm-starting
 * the nuisance optimisation from `start` (updated in place). */
static int profile_at(const NllCtx* ctx, int pi, double val, double start[3],
                      double opt[3], double* f) {
    if (!profile(ctx, pi, val, start, opt, f)) return 0;
    memcpy(start, opt, sizeof(double) * 3);
    return 1;
}

/* Find where the profile nll crosses `target`. Expands outward from the MLE
 * until the crossing is bracketed, then bisects. The profile is curved, so
 * interpolating across a single coarse step understates the interval badly. */
static double profile_bound(const NllCtx* ctx, const SlopeFit* fit, int pi,
                            double target, int sign) {
    double centre = pi == 0 ? fit->rho0 : (pi == 1 ? fit->s0 : fit->p);
    double start[3] = { fit->rho0, fit->s0, fit->p };
    double opt[3], f;

    /* Parameter-space limits. rho0 is bounded to (0, 1); the rest are open,
     * guarded only against overflow. */
    double lo_lim = pi == 0 ? 1e-6 : 1e-6;
    double hi_lim = pi == 0 ? 1.0 : 1e4;

    double a = centre, fa;
    /* At the MLE the profile nll is the minimum; if it already exceeds the
     * target the interval is degenerate (the optimiser disagreed slightly). */
    if (profile_at(ctx, pi, a, start, opt, &fa) && fa > target) return a;

    double step = pi == 0 ? 0.01 : 0.04;
    double b = a;
    int bracketed = 0;
    for (int it = 0; it < 120; it++) {
        double next = sign > 0 ? b + step * (pi == 0 ? 1.0 : (b > 0 ? b : 1.0))
                               : b - step * (pi == 0 ? 1.0 : (b > 0 ? b : 1.0));
        if (next < lo_lim) next = lo_lim;
        if (next > hi_lim) next = hi_lim;
        if (next == b) return b;
        if (!profile_at(ctx, pi, next, start, opt, &f)) return b;
        if (f >= target) { bracketed = 1; b = next; break; }
        a = next;
        b = next;
    }
    if (!bracketed) return b;

    /* Bisect [a, b] where f(a) <= target <= f(b). */
    for (int it = 0; it < 60; it++) {
        double mid = 0.5 * (a + b);
        if (!profile_at(ctx, pi, mid, start, opt, &f)) break;
        if (f >= target) b = mid;
        else a = mid;
    }
    return 0.5 * (a + b);
}

int slope_fit(SlopeModelKind kind, const SlopeBin* bins, int n, double alpha,
              SlopeFit* out) {
    if (!out || !bins || alpha <= 0.0 || alpha >= 1.0) return -1;
    if (n < slope_model_nparam(kind) + 2) return -1;
    memset(out, 0, sizeof(*out));
    out->kind = kind;
    out->nparam = slope_model_nparam(kind);

    NllCtx ctx = { kind, bins, n };

    /* Deterministic multi-start over a coarse grid. */
    const double rho_grid[] = { 0.5, 0.75, 0.92, 0.99 };
    const double s0_grid[] = { 0.1, 0.25, 0.5, 1.0, 2.0 };
    const double p_grid[] = { 0.6, 1.0, 1.6, 2.4, 3.5 };
    int nr = 4, ns = 5, np = 5;

    double best[3] = { 0.9, 1.0, 1.0 }, best_nll = INFINITY;
    for (int a = 0; a < nr; a++)
        for (int b = 0; b < ns; b++)
            for (int c = 0; c < (kind == SLOPE_MODEL_EXP ? 1 : np); c++) {
                double start[3] = { rho_grid[a], s0_grid[b],
                                    kind == SLOPE_MODEL_EXP ? 1.0 : p_grid[c] };
                double prm[3];
                double f = ml_full(&ctx, start, prm);
                if (isfinite(f) && f < best_nll) {
                    best_nll = f;
                    memcpy(best, prm, sizeof(best));
                }
            }
    if (!isfinite(best_nll)) return -1;

    out->converged = 1;
    out->nll = best_nll;
    out->rho0 = best[0];
    out->s0 = best[1];
    out->p = kind == SLOPE_MODEL_EXP ? 1.0 : best[2];
    out->aic = 2.0 * best_nll + 2.0 * out->nparam;

    /* Profile-likelihood CIs at level 1-alpha. */
    double target = best_nll + 0.5 * chi2_1_quantile(1.0 - alpha);
    for (int pi = 0; pi < out->nparam; pi++) {
        out->ci_lo[pi] = profile_bound(&ctx, out, pi, target, -1);
        out->ci_hi[pi] = profile_bound(&ctx, out, pi, target, +1);
    }
    double z = slope_norm_quantile(1.0 - 0.5 * alpha);
    for (int i = 0; i < 3; i++)
        out->se[i] = (out->ci_hi[i] - out->ci_lo[i]) / (2.0 * z);
    if (kind == SLOPE_MODEL_EXP) {
        out->ci_lo[2] = out->ci_hi[2] = 1.0;
        out->se[2] = 0.0;
    }

    /* Pearson goodness of fit and AICc over informative bins. */
    double chi = 0.0;
    int used = 0;
    for (int i = 0; i < n; i++) {
        double m = bins[i].exposure;
        if (!(m > 0.0)) continue;
        used++;
        double r = slope_attenuation(kind, out->rho0, out->s0, out->p,
                                     bins[i].slope);
        double var = m * r * (1.0 - r);
        if (var > 1e-8) {
            double d = bins[i].retained - m * r;
            chi += d * d / var;
        }
    }
    out->npoints = used;
    out->df = used - out->nparam;
    out->pearson = chi;
    out->pearson_per_df = out->df > 0 ? chi / out->df : 0.0;
    out->aicc = out->aic;
    if (used - out->nparam - 1 > 0)
        out->aicc += 2.0 * out->nparam * (out->nparam + 1.0) /
                     (double)(used - out->nparam - 1);
    return 0;
}

int slope_fit_best(const SlopeBin* bins, int n, double alpha, SlopeFit* out,
                   SlopeFit all[3]) {
    double best_aicc = INFINITY;
    int rc = -1;
    for (int k = 0; k < 3; k++) {
        SlopeFit tmp;
        int r = slope_fit((SlopeModelKind)k, bins, n, alpha, &tmp);
        if (r != 0) continue;
        if (all) all[k] = tmp;
        rc = 0;
        if (tmp.aicc < best_aicc) {
            best_aicc = tmp.aicc;
            if (out) *out = tmp;
        }
    }
    return rc;
}

double slope_predict_band(const SlopeFit* fit, const SlopeBin* bins, int n) {
    double num = 0.0, den = 0.0;
    for (int i = 0; i < n; i++) {
        if (!(bins[i].exposure > 0.0)) continue;
        num += bins[i].exposure * slope_attenuation(fit->kind, fit->rho0, fit->s0,
                                                    fit->p, bins[i].slope);
        den += bins[i].exposure;
    }
    return den > 0.0 ? num / den : 0.0;
}

int slope_deploy_ok(const SlopeFit* fit, double rho0, double s0, double p,
                    char* msg, size_t msglen) {
    const double vals[3] = { rho0, s0, p };
    const char* names[3] = { "rho0", "s0", "p" };
    int ncheck = fit->kind == SLOPE_MODEL_EXP ? 2 : 3;
    int ok = 1;
    size_t used = 0;
    if (msg && msglen) msg[0] = '\0';
    for (int i = 0; i < ncheck; i++) {
        if (!isfinite(vals[i]) || vals[i] < fit->ci_lo[i] ||
            vals[i] > fit->ci_hi[i]) {
            ok = 0;
            if (msg && used < msglen)
                used += (size_t)snprintf(msg + used, msglen - used,
                                         "%s%s=%.4g outside [%.4g, %.4g]",
                                         used ? "; " : "", names[i], vals[i],
                                         fit->ci_lo[i], fit->ci_hi[i]);
        }
    }
    return ok;
}
