#include <math.h>
#include <stdio.h>
#include <string.h>
#include "slope_fit.h"
#include "slope_stats.h"

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* Deterministic normal-ish sample for fixtures: sum of uniforms. */
static double gauss01(unsigned seed) {
    unsigned h = seed * 2654435761u;
    double s = 0.0;
    for (int i = 0; i < 12; i++) {
        h = h * 1103515245u + 12345u;
        s += (double)((h >> 8) & 0xffff) / 65535.0;
    }
    return (s - 6.0) / 1.0;
}

static int synth(SlopeBin* bins, int nbins, double rho0, double s0, double p,
                 double exposure, unsigned seed) {
    for (int i = 0; i < nbins; i++) {
        double slope = 0.02 + i * (0.03);
        double r = slope_attenuation(SLOPE_MODEL_STRETCHED, rho0, s0, p, slope);
        bins[i].slope = slope;
        bins[i].exposure = exposure;
        bins[i].retained = round(exposure * r + 0.35 * sqrt(exposure * r * (1 - r)) *
                                 gauss01(seed + (unsigned)i));
        if (bins[i].retained < 0) bins[i].retained = 0;
        if (bins[i].retained > exposure) bins[i].retained = exposure;
    }
    return nbins;
}

int main(void) {
    /* quantiles */
    check("norm_q_975", fabs(slope_norm_quantile(0.975) - 1.959964) < 1e-5);
    check("norm_q_025", fabs(slope_norm_quantile(0.025) + 1.959964) < 1e-5);
    check("t_q_975_df10", fabs(slope_t_quantile(0.975, 10) - 2.228139) < 2e-4);
    double z2 = slope_norm_quantile(0.975);
    check("t_converges", fabs(slope_t_quantile(0.975, 100000) - z2) < 3e-4);

    /* model behaviour */
    double r = slope_attenuation(SLOPE_MODEL_STRETCHED, 0.9, 0.5, 2.0, 0.5);
    check("atten_at_s0", fabs(r - 0.9 * exp(-1.0)) < 1e-9);
    check("atten_flat", fabs(slope_attenuation(SLOPE_MODEL_EXP, 0.9, 0.5, 1, 0) - 0.9) < 1e-9);
    check("atten_bounded", slope_attenuation(SLOPE_MODEL_HILL, 0.9, 0.5, 3, 100) < 0.1);

    /* recovery of a known stretched-exponential truth */
    SlopeBin bins[40];
    int n = synth(bins, 40, 0.90, 0.62, 1.6, 400.0, 11u);
    SlopeFit fit;
    check("fit_ok", slope_fit(SLOPE_MODEL_STRETCHED, bins, n, 0.05, &fit) == 0);
    check("fit_rho0_se", fabs(fit.rho0 - 0.90) < fit.se[0] * 4 + 0.01);
    check("fit_s0_se", fabs(fit.s0 - 0.62) < fit.se[1] * 4 + 0.01);
    check("fit_p_se", fabs(fit.p - 1.6) < fit.se[2] * 4 + 0.05);
    check("ci_order", fit.ci_lo[1] < fit.s0 && fit.s0 < fit.ci_hi[1]);
    check("ci_covered_truth", fit.ci_lo[1] <= 0.62 && 0.62 <= fit.ci_hi[1]);
    check("se_positive", fit.se[0] > 0 && fit.se[1] > 0 && fit.se[2] > 0);

    /* shrinking exposure must widen intervals */
    SlopeBin coarse[40];
    synth(coarse, 40, 0.90, 0.62, 1.6, 25.0, 11u);
    SlopeFit fitc;
    slope_fit(SLOPE_MODEL_STRETCHED, coarse, 40, 0.05, &fitc);
    check("ci_widen_low_n", (fitc.ci_hi[1] - fitc.ci_lo[1]) >
                            (fit.ci_hi[1] - fit.ci_lo[1]));

    /* model selection picks the three-parameter stretched model for
     * stretched-truth data over the nested exponential */
    SlopeFit all[3], best;
    slope_fit_best(bins, n, 0.05, &best, all);
    check("best_is_stretched", best.kind == SLOPE_MODEL_STRETCHED);
    check("aicc_nested", all[0].aicc > all[1].aicc - 1e-9);

    /* deployment gate */
    char msg[256];
    check("deploy_pass", slope_deploy_ok(&fit, fit.rho0, fit.s0, fit.p, msg, sizeof(msg)) == 1);
    check("deploy_fail", slope_deploy_ok(&fit, 1.0, 50.0, fit.p, msg, sizeof(msg)) == 0);
    check("deploy_msg_populated", strlen(msg) > 0);
    /* the shipped no-attenuation default must be rejected */
    check("deploy_reject_off", slope_deploy_ok(&fit, 1.0, 1e9, 1.0, msg, sizeof(msg)) == 0);
    check("deploy_edge_inclusive", slope_deploy_ok(&fit, fit.ci_hi[0], fit.s0, fit.p,
                                                   msg, sizeof(msg)) == 1);

    /* Regression: the profile-likelihood interval used to be taken by
     * interpolating across one 12%-wide scan step, which understated the 95%
     * interval by roughly 9x. Replicate a binomial survey and require the
     * intervals to cover the truth at close to the nominal rate. */
    {
        SlopeBin rep[40];
        int covers_rho = 0, covers_s0 = 0, covers_p = 0, reps = 10;
        unsigned rs = 424242u;
        for (int r = 0; r < reps; r++) {
            for (int i = 0; i < 40; i++) {
                double slope = 0.02 + i * 0.03;
                double pr = slope_attenuation(SLOPE_MODEL_STRETCHED, 0.92, 0.55,
                                              1.5, slope);
                int k = 0;
                for (int t = 0; t < 400; t++) {
                    rs = rs * 1664525u + 1013904223u;
                    if ((rs >> 8) / 16777216.0 < pr) k++;
                }
                rep[i].slope = slope;
                rep[i].exposure = 400;
                rep[i].retained = k;
            }
            SlopeFit rf;
            if (slope_fit(SLOPE_MODEL_STRETCHED, rep, 40, 0.05, &rf) != 0) continue;
            if (rf.ci_lo[0] <= 0.92 && 0.92 <= rf.ci_hi[0]) covers_rho++;
            if (rf.ci_lo[1] <= 0.55 && 0.55 <= rf.ci_hi[1]) covers_s0++;
            if (rf.ci_lo[2] <= 1.5 && 1.5 <= rf.ci_hi[2]) covers_p++;
            rs += 7919u;
        }
        /* 8/10 is a loose floor; the broken intervals covered 1/10. */
        check("ci_covers_rho0", covers_rho >= 8);
        check("ci_covers_s0", covers_s0 >= 8);
        check("ci_covers_p", covers_p >= 8);
    }

    /* strata */
    SlopeStratum strata[SLOPE_STRATA];
    slope_strata(bins, n, &fit, 0.05, strata);
    double sum_exp = 0;
    for (int i = 0; i < SLOPE_STRATA; i++) sum_exp += strata[i].exposure;
    check("strata_cover_all", fabs(sum_exp - (double)n * 400.0) < 1e-6);
    check("strata_model_decays",
          strata[0].predicted >= strata[1].predicted &&
          strata[1].predicted >= strata[2].predicted);
    double wlo, whi;
    slope_wilson(5, 10, 0.05, &wlo, &whi);
    check("wilson_centre", fabs(0.5 * (wlo + whi) - 0.5) < 0.02);
    check("wilson_bounds", wlo > 0.0 && whi < 1.0 && wlo < whi);
    slope_wilson(0, 10, 0.05, &wlo, &whi);
    check("wilson_zero", fabs(wlo) < 1e-12 && whi > 0.2);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
