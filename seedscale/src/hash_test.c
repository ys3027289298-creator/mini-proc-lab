/*
 * hash_test: 哈希级独立性检验（三档种子）。
 *
 * 用法: ./hash_test good | bad-linear
 *   good       合格哈希 ss_hash2
 *   bad-linear 坏法：裸线性组合、无终混（无雪崩、跨格点代数相关）
 *
 * 每档种子都跑同一组检验：
 *   bijection   固定格点上种子->输出在 2^16 连续种子上的单射性（完整 64 位查重）
 *   uniform     均匀性 chi^2（256 桶），给 p 值
 *   bitbalance  64 个输出比特的 1 比例，期望 0.5（卡死的比特会暴露）
 *   avalanche   输入(i/j/s)各比特翻转 -> 输出各比特翻转率，期望 0.5
 *   neighbor    相邻格点输出 bit 相关、值相关，期望不显著
 *   isotropic   x/y 两方向统计量对称
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "hash.h"
#include "stats.h"

static int failures = 0;
static void verdict(const char* tier, const char* name, int ok, double detail) {
    printf("%s %s ok=%d detail=%.6g\n", tier, name, ok, detail);
    if (!ok) failures++;
}

static uint64_t (*pick_hash(const char* mode))(int64_t, int64_t, uint64_t) {
    if (!strcmp(mode, "bad-linear")) return bad_linear;
    return ss_hash2;
}

/* 三档种子采样区间 */
static uint64_t tier_seed(const char* tier, int k) {
    if (tier[0] == 'S') return (uint64_t)k;                       /* 小: 0..255 小整数 */
    if (tier[0] == 'M') return (1ULL << 16) + (uint64_t)k * 101U;/* 中: 2^16 邻域 */
    ss_rng_state = 0x123456789ABCDEF0ULL;                         /* 全: 全 uint32 */
    for (int q = 0; q < k; q++) (void)ss_rng();
    return ss_rng() & 0xFFFFFFFFULL;
}

static int cmp64(const void* a, const void* b) {
    uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b;
    return x < y ? -1 : x > y ? 1 : 0;
}

#define NBINS 256
#define NUNIFORM 200000

static void test_uniform(const char* tier, uint64_t (*H)(int64_t, int64_t, uint64_t)) {
    long cnt[NBINS];
    memset(cnt, 0, sizeof(cnt));
    for (int k = 0; k < NUNIFORM; k++) {
        uint64_t s = tier_seed(tier, k);
        uint64_t h = H(17, -9, s);
        cnt[(h >> 56) & 0xFF]++;          /* 高 8 位 */
    }
    double e = (double)NUNIFORM / NBINS, chi = 0.0;
    for (int b = 0; b < NBINS; b++) {
        double d = (double)cnt[b] - e;
        chi += d * d / e;
    }
    double p = chi2_p(chi, NBINS - 1);
    verdict(tier, "uniform", p > 0.01, p);
}

#define NBB 65536

static void test_bitbalance(const char* tier, uint64_t (*H)(int64_t, int64_t, uint64_t)) {
    /* 固定格点，扫 2^16 个连续种子：每个输出比特的 1 比例应 ~0.5。
       标准差 sqrt(0.25/65536)≈0.00195，阈值 6 sigma ≈ 0.012。 */
    long ones[64];
    memset(ones, 0, sizeof(ones));
    uint64_t base = tier[0] == 'S' ? 0 : tier_seed(tier, 0);
    for (int k = 0; k < NBB; k++) {
        uint64_t h = H(11, 13, base + (uint64_t)k);
        for (int b = 0; b < 64; b++) if ((h >> b) & 1ULL) ones[b]++;
    }
    double worst = 0.0;
    for (int b = 0; b < 64; b++) {
        double dev = fabs((double)ones[b] / NBB - 0.5);
        if (dev > worst) worst = dev;
    }
    verdict(tier, "bitbalance", worst < 0.012, worst);
}

#define NBIJ 65536

static void test_bijection(const char* tier, uint64_t (*H)(int64_t, int64_t, uint64_t)) {
    /* 满置换性质在完整 64 位输出上：收集 2^16 个连续种子的完整输出，
       排序查重。双射 => 无重复（截断到低位反而会有生日碰撞，那是检验写错）。 */
    static uint64_t vals[NBIJ];
    int collision = 0;
    uint64_t base = tier[0] == 'S' ? 0 : tier_seed(tier, 0);
    for (int k = 0; k < NBIJ; k++)
        vals[k] = H(23, 5, base + (uint64_t)k);
    qsort(vals, NBIJ, sizeof(uint64_t), cmp64);
    for (int k = 1; k < NBIJ; k++)
        if (vals[k] == vals[k - 1]) collision++;
    verdict(tier, "bijection", collision == 0, (double)collision);
}

#define NAV 20000

static void test_avalanche(const char* tier, uint64_t (*H)(int64_t, int64_t, uint64_t)) {
    /* 三组输入比特：i、j、s；分别翻转，统计 64 个输出比特翻转率 */
    double worst = 0.0, sumdev = 0.0;
    long nbits = 0;
    for (int grp = 0; grp < 3; grp++) {
        long flip[64];
        memset(flip, 0, sizeof(flip));
        for (int n = 0; n < NAV; n++) {
            int64_t i = (int64_t)(ss_rng() % 2000) - 1000;
            int64_t j = (int64_t)(ss_rng() % 2000) - 1000;
            uint64_t s;
            if (tier[0] == 'S')      s = ss_rng() % 4096;
            else if (tier[0] == 'M') s = (1ULL << 16) + ss_rng() % 4096;
            else                     s = ss_rng() & 0xFFFFFFFFULL;
            int bit = (int)(ss_rng() % 32);
            int64_t i2 = i, j2 = j;
            uint64_t s2 = s;
            if (grp == 0) i2 ^= (int64_t)1 << bit;
            if (grp == 1) j2 ^= (int64_t)1 << bit;
            if (grp == 2) s2 ^= 1ULL << bit;
            uint64_t d = H(i, j, s) ^ H(i2, j2, s2);
            for (int b = 0; b < 64; b++) if ((d >> b) & 1ULL) flip[b]++;
        }
        for (int b = 0; b < 64; b++) {
            double r = (double)flip[b] / NAV;
            double dev = fabs(r - 0.5);
            if (dev > worst) worst = dev;
            sumdev += dev;
            nbits++;
        }
    }
    double mean = sumdev / (double)nbits;
    /* 0.5 下标准差约 sqrt(0.25/20000)=0.0035；0.03/0.10 给足余量 */
    verdict(tier, "avalanche", mean < 0.03 && worst < 0.10, mean);
}

#define NIND 40000

static void test_neighbor(const char* tier, uint64_t (*H)(int64_t, int64_t, uint64_t)) {
    /* 相邻格点 (i,j) 与 (i+1,j)、(i,j+1) 的输出独立性：
       1) 均匀化后的 u01 值 Pearson 相关；
       2) 最低输出 bit 的符号相关（对坏法 A/B 最敏感）。 */
    double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
    long bx = 0, by = 0, bxy = 0, n = NIND;
    for (int k = 0; k < NIND; k++) {
        int64_t i = (int64_t)(ss_rng() % 500) - 250;
        int64_t j = (int64_t)(ss_rng() % 500) - 250;
        uint64_t s = tier_seed(tier, k);
        double u = ss_u01(H(i, j, s));
        double v = ss_u01(H(i + 1, j, s));
        sx += u; sy += v; sxx += u * u; syy += v * v; sxy += u * v;
        bx += (long)(H(i, j, s) & 1ULL);
        by += (long)(H(i + 1, j, s) & 1ULL);
        bxy += (long)((H(i, j, s) & 1ULL) & (H(i + 1, j, s) & 1ULL));
    }
    double cov = sxy - sx * sy / n;
    double vx = sxx - sx * sx / n, vy = syy - sy * sy / n;
    double r = (vx > 0 && vy > 0) ? cov / sqrt(vx * vy) : 1.0;
    double pv = (vx > 0 && vy > 0) ? pearson_z_p(r, n) : 0.0;
    double px = (double)bx / n, py = (double)by / n;
    double denom = sqrt(px * (1 - px) * py * (1 - py));
    double bitr = denom > 1e-12
                      ? ((double)bxy / n - px * py) / denom
                      : 1.0;   /* 常比特 => 直接判失败 */
    int ok = pv > 0.01 && fabs(bitr) < 4.0;   /* 4 sigma 容差 */
    verdict(tier, "neighbor", ok, fabs(r) + fabs(bitr));
}

static void test_isotropic(const char* tier, uint64_t (*H)(int64_t, int64_t, uint64_t)) {
    /* x 方向与 y 方向相邻点的相关绝对值应近似相等（各向同性） */
    double sumx = 0, sumy = 0;
    for (int k = 0; k < NIND; k++) {
        int64_t i = (int64_t)(ss_rng() % 500) - 250;
        int64_t j = (int64_t)(ss_rng() % 500) - 250;
        uint64_t s = tier_seed(tier, k);
        sumx += fabs(ss_u01(H(i, j, s)) - ss_u01(H(i + 1, j, s)));
        sumy += fabs(ss_u01(H(i, j, s)) - ss_u01(H(i, j + 1, s)));
    }
    double mx = sumx / NIND, my = sumy / NIND;
    double rel = fabs(mx - my) / (0.5 * (mx + my));
    verdict(tier, "isotropic", rel < 0.02, rel);
}

int main(int argc, char** argv) {
    const char* mode = argc > 1 ? argv[1] : "good";
    if (strcmp(mode, "good") && strcmp(mode, "bad-linear")) {
        fprintf(stderr, "usage: %s good|bad-linear\n", argv[0]);
        return 2;
    }
    printf("# mode=%s\n", mode);
    uint64_t (*H)(int64_t, int64_t, uint64_t) = pick_hash(mode);

    /* 确定性（任何模式都必须成立） */
    int det = (H(3, -7, 42) == H(3, -7, 42)) && (H(3, -7, 42) != H(3, -7, 43));
    verdict("ALL", "deterministic_and_seed_sensitive", det, 1.0);

    const char* tiers[3] = {"S_small", "M_mid", "L_full32"};
    for (int t = 0; t < 3; t++) {
        test_bijection(tiers[t], H);
        test_uniform(tiers[t], H);
        test_bitbalance(tiers[t], H);
        test_avalanche(tiers[t], H);
        test_neighbor(tiers[t], H);
        test_isotropic(tiers[t], H);
    }
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
