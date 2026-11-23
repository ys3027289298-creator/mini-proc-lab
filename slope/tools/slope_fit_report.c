/* slope_fit_report: measure scatter density against terrain slope, fit the
 * attenuation curve, and report slope-stratified statistics.
 *
 * Data source, in order of preference:
 *   --csv PATH  measured survey bins (slope,exposure,retained per line)
 *   --fixture   built-in calibration fixture (default when no CSV is given)
 *
 * The fixture drives the real scatter + terrain code with a known ground-truth
 * curve, so the fitter and its confidence intervals can be validated end to
 * end. It is a calibration fixture, not field data: use --csv with a real
 * survey before trusting the numbers for production.
 *
 * Exit status:
 *   0  the evaluated parameter set is inside every confidence interval
 *   1  the parameter set is outside at least one interval (NOT SHIPPABLE)
 *   2  bad input
 */
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "scatter.h"
#include "slope_fit.h"
#include "slope_stats.h"
#include "terrain.h"

#define MAXBINS 512
#define MAXCSV 4096

/* ------------------------------------------------------------------ */
/* fixture                                                            */
/* ------------------------------------------------------------------ */

typedef struct {
    double rho0, s0, p;
} Truth;

static float truth_density(void* ctx, float slope) {
    const Truth* t = (const Truth*)ctx;
    return (float)slope_attenuation(SLOPE_MODEL_STRETCHED, t->rho0, t->s0, t->p,
                                    slope);
}

static float field_slope_cb(void* ctx, float x, float y) {
    return field_slope_at((const Field*)ctx, x, y);
}

/* Default fixture: a moderate slope decay that would let too much grass onto
 * steep ground if left uncorrected. */
static const Truth kTruth = { 0.92, 0.55, 1.5 };

typedef struct {
    int nx, ny;
    float cell;
    unsigned seed;
    int octaves;
    float amplitude, base_freq;
    float hist_lo, hist_width;
    int hist_bins;
} FixtureCfg;

static const FixtureCfg kFixture = {
    768, 768, 0.25f, 2024u, 6, 12.0f, 0.09f, 0.0f, 0.01f, 320
};

/* ------------------------------------------------------------------ */
/* CSV                                                                */
/* ------------------------------------------------------------------ */

static int read_csv(const char* path, SlopeBin* bins, int max, int* n) {
    FILE* f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "cannot open %s: %s\n", path, strerror(errno));
        return -1;
    }
    char line[512];
    int count = 0;
    while (fgets(line, sizeof(line), f)) {
        char* s = line;
        while (*s == ' ' || *s == '\t') s++;
        if (*s == '#' || *s == '\n' || *s == '\r' || *s == '\0') continue;
        double slope, exposure, retained;
        if (sscanf(s, "%lf %*[,;] %lf %*[,;] %lf", &slope, &exposure,
                   &retained) != 3) {
            /* also accept a header line */
            if (count == 0 && (strstr(s, "slope") || strstr(s, "Slope")))
                continue;
            fprintf(stderr, "%s: bad row: %s", path, s);
            fclose(f);
            return -1;
        }
        if (exposure <= 0.0 || retained < 0.0 || retained > exposure) {
            fprintf(stderr, "%s: invalid counts in row: %s", path, s);
            fclose(f);
            return -1;
        }
        if (count >= max) {
            fprintf(stderr, "%s: more than %d bins\n", path, max);
            fclose(f);
            return -1;
        }
        bins[count].slope = slope;
        bins[count].exposure = exposure;
        bins[count].retained = retained;
        count++;
    }
    fclose(f);
    *n = count;
    return count > 0 ? 0 : -1;
}

/* ------------------------------------------------------------------ */
/* fixture run                                                        */
/* ------------------------------------------------------------------ */

static int run_fixture(SlopeBin* bins, int* n, int* total_exposure,
                       double* expected_total) {
    static Field f;
    field_generate_fbm(&f, 256, 256, 0.5f, kFixture.amplitude, kFixture.octaves,
                       kFixture.base_freq, 77u);

    int nb = kFixture.hist_bins;
    if (nb > MAXBINS) return -1;
    static int exposure[MAXBINS], retained[MAXBINS];
    static double slope_sum[MAXBINS], expect_sum[MAXBINS];
    int outside = scatter_measure(kFixture.nx, kFixture.ny, kFixture.cell,
                                  kFixture.seed, field_slope_cb, &f,
                                  truth_density, (void*)&kTruth,
                                  kFixture.hist_lo, kFixture.hist_width, nb,
                                  exposure, retained, slope_sum, expect_sum);
    if (outside < 0) return -1;

    int count = 0, total = 0;
    double exp_total = 0.0;
    for (int i = 0; i < nb; i++) {
        if (exposure[i] <= 0) continue;
        bins[count].slope = slope_sum[i] / exposure[i];
        bins[count].exposure = exposure[i];
        bins[count].retained = retained[i];
        total += exposure[i];
        exp_total += expect_sum[i];
        count++;
    }
    *n = count;
    if (total_exposure) *total_exposure = total;
    if (expected_total) *expected_total = exp_total;
    return 0;
}

/* ------------------------------------------------------------------ */
/* printing                                                           */
/* ------------------------------------------------------------------ */

static void print_rule(const char* title) {
    printf("\n== %s ==\n", title);
}

static void print_model_table(const SlopeFit all[3]) {
    printf("%-11s %4s %12s %12s %8s %10s %9s\n", "model", "k", "nll", "aicc",
           "pearson", "pearson/df", "s0[rise/run]");
    for (int i = 0; i < 3; i++) {
        printf("%-11s %4d %12.3f %12.3f %8.2f %10.2f %9.4f\n",
               slope_model_name((SlopeModelKind)i), all[i].nparam, all[i].nll,
               all[i].aicc, all[i].pearson, all[i].pearson_per_df, all[i].s0);
    }
}

static void print_fit(const SlopeFit* f) {
    printf("formula: rho(s) = %.4f * exp(-(s / %.4f)^%.3f)\n", f->rho0, f->s0,
           f->p);
    printf("%-6s %12s %12s %12s %12s\n", "param", "estimate", "std.err",
           "95% CI low", "95% CI high");
    const char* names[3] = { "rho0", "s0", "p" };
    for (int i = 0; i < f->nparam; i++)
        printf("%-6s %12.5f %12.5f %12.5f %12.5f\n", names[i],
               i == 0 ? f->rho0 : (i == 1 ? f->s0 : f->p), f->se[i],
               f->ci_lo[i], f->ci_hi[i]);
    printf("nll=%.3f aicc=%.3f df=%d pearson=%.2f (%.2f/df)\n", f->nll, f->aicc,
           f->df, f->pearson, f->pearson_per_df);
}

static void print_strata(const SlopeStratum s[SLOPE_STRATA]) {
    printf("%-22s %9s %9s %8s %18s %9s %9s %7s\n", "slope layer", "exposure",
           "retained", "obs", "wilson CI", "model", "resid", "p");
    printf("%-22s %9s %9s %8s %18s %9s %9s %7s\n", "", "", "", "density",
           "", "density", "obs-pred", "value");
    for (int i = 0; i < SLOPE_STRATA; i++) {
        char label[32];
        if (s[i].hi < 0.0) snprintf(label, sizeof(label), "steep  >= %.2f", s[i].lo);
        else snprintf(label, sizeof(label), "layer %d  [%.2f, %.2f)", i + 1,
                      s[i].lo, s[i].hi);
        char ci[32];
        snprintf(ci, sizeof(ci), "[%.3f, %.3f]", s[i].obs_lo, s[i].obs_hi);
        printf("%-22s %9.0f %9.0f %8.3f %18s %9.3f %+9.3f %7.3g\n", label,
               s[i].exposure, s[i].retained, s[i].observed, ci, s[i].predicted,
               s[i].resid, s[i].pval);
    }
}

static void usage(const char* argv0) {
    printf("usage: %s [--csv PATH] [--fixture] [--alpha A]\n"
           "          [--rho0 X --s0 X --p X]\n"
           "  Evaluates the given parameters against the fitted confidence\n"
           "  intervals; exits 1 if any lies outside (not shippable).\n"
           "  Without --rho0/--s0/--p the shipped no-attenuation behaviour\n"
           "  (rho0=1, s0=1e9) is evaluated.\n",
           argv0);
}

int main(int argc, char** argv) {
    const char* csv = NULL;
    double alpha = 0.05;
    double rho0 = 1.0, s0 = 1e9, p = 1.0;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--csv") && i + 1 < argc) {
            csv = argv[++i];
        } else if (!strcmp(argv[i], "--fixture")) {
            csv = NULL;
        } else if (!strcmp(argv[i], "--alpha") && i + 1 < argc) {
            alpha = atof(argv[++i]);
        } else if (!strcmp(argv[i], "--rho0") && i + 1 < argc) {
            rho0 = atof(argv[++i]);
        } else if (!strcmp(argv[i], "--s0") && i + 1 < argc) {
            s0 = atof(argv[++i]);
        } else if (!strcmp(argv[i], "--p") && i + 1 < argc) {
            p = atof(argv[++i]);
        } else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            usage(argv[0]);
            return 0;
        } else {
            usage(argv[0]);
            return 2;
        }
    }
    static SlopeBin bins[MAXCSV];
    int n = 0, total = 0;
    double expected_total = 0.0;
    if (csv) {
        if (read_csv(csv, bins, MAXCSV, &n) != 0) {
            fprintf(stderr, "no usable rows in %s\n", csv);
            return 2;
        }
        for (int i = 0; i < n; i++) total += (int)bins[i].exposure;
        printf("data: %s (%d bins, %.0f candidates)\n", csv, n,
               (double)total);
    } else {
        if (run_fixture(bins, &n, &total, &expected_total) != 0) {
            fprintf(stderr, "fixture run failed\n");
            return 2;
        }
        printf("data: built-in calibration fixture (%d bins, %d candidates)\n",
               n, total);
        printf("      truth rho0=%.3f s0=%.3f p=%.3f (stretched exponential)\n",
               kTruth.rho0, kTruth.s0, kTruth.p);
        printf("      NOT field data. Supply --csv for a real survey.\n");
    }

    SlopeFit all[3];
    SlopeFit best;
    memset(all, 0, sizeof(all));
    if (slope_fit_best(bins, n, alpha, &best, all) != 0) {
        fprintf(stderr, "fit failed\n");
        return 2;
    }

    print_rule("model comparison");
    print_model_table(all);
    printf("selected by lowest AICc: %s\n", slope_model_name(best.kind));

    print_rule("fitted attenuation curve");
    print_fit(&best);

    SlopeStratum strata[SLOPE_STRATA];
    slope_strata(bins, n, &best, alpha, strata);
    print_rule("slope-stratified statistics");
    print_strata(strata);

    print_rule("deployment gate");
    printf("evaluating rho0=%.4g s0=%.4g p=%.4g against the fitted %s CIs\n",
           rho0, s0, p, slope_model_name(best.kind));
    char msg[512];
    int ok = slope_deploy_ok(&best, rho0, s0, p, msg, sizeof(msg));
    if (ok) {
        printf("PASS: parameters lie inside every confidence interval\n");
    } else {
        printf("FAIL: %s\n", msg);
        printf("      parameters outside the fitted intervals are not shippable\n");
    }
    return ok ? 0 : 1;
}
