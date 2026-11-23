#ifndef SEEDSCALE_STATS_H
#define SEEDSCALE_STATS_H
#include <math.h>
#include <stdint.h>

/* 确定性 PRNG（splitmix64），固定主种子，保证实验可复现。 */
static uint64_t ss_rng_state = 0x123456789ABCDEF0ULL;
static inline uint64_t ss_rng(void) {
    uint64_t z = (ss_rng_state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

/* 正则化不完全伽马函数（Numerical Recipes 风格），用于 χ² 检验的上侧 p。 */
#define NR_ITMAX 200
#define NR_EPS 3e-14
#define NR_FPMIN 1e-300

static double gammln(double a) { return lgamma(a); }

static double gser(double a, double x) {
    double gln = gammln(a);
    double ap = a, del = 1.0 / a, sum = del;
    for (int n = 1; n <= NR_ITMAX; n++) {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if (fabs(del) < fabs(sum) * NR_EPS) break;
    }
    return sum * exp(-x + a * log(x) - gln);
}

static double gcf(double a, double x) {
    double b = x + 1.0 - a, c = 1.0 / NR_FPMIN, d = 1.0 / b, h = d;
    for (int i = 1; i <= NR_ITMAX; i++) {
        double an = -(double)i * ((double)i - a);
        b += 2.0;
        d = an * d + b;
        if (fabs(d) < NR_FPMIN) d = NR_FPMIN;
        c = b + an / c;
        if (fabs(c) < NR_FPMIN) c = NR_FPMIN;
        d = 1.0 / d;
        double del = d * c;
        h *= del;
        if (fabs(del - 1.0) < NR_EPS) break;
    }
    return exp(-x + a * log(x) - gammln(a)) * h;
}

/* 上侧尾概率 P(χ²_df >= x)
   x<df 时级数 gser 收敛、x>df 时连分式 gcf 收敛（NR 的分支规则）。 */
static inline double chi2_p(double x, double df) {
    if (x <= 0.0) return 1.0;
    return x < df + 1.0 ? 1.0 - gser(df / 2.0, x / 2.0)
                        : gcf(df / 2.0, x / 2.0);
}

/* Pearson 相关系数的双侧 p（Fisher z 近似，n 足够大时精确）。 */
static inline double pearson_z_p(double r, long n) {
    if (r <= -0.999999) r = -0.999999;
    if (r >= 0.999999) r = 0.999999;
    double z = atanh(r) * sqrt((double)n - 3.0);
    return erfc(fabs(z) / sqrt(2.0));
}

#endif
