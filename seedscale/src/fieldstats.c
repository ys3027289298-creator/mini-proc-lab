/*
 * fieldstats: 逐种子生成值噪声高度场并输出特征统计 CSV。
 *
 * 用法:
 *   ./fieldstats <good|bad-linear> <vn|nn|sx> <out.csv> [fbm]
 *     哈希模式: good=ss_hash2, bad-linear=对照坏哈希
 *     采样模式: vn=smoothstep 值噪声(合格)
 *               nn=最近邻(采样坏法, 混叠)
 *               sx=种子耦合频率(管线坏法, 违反 C2: 格距随种子低3位变化)
 *     fbm: 可选, 单八度改为 4 八度 fBm
 *
 * 三档种子各 128 个:
 *   S: 0..127   M: 2^16 + k*101   L: 固定主种子 PRNG 抽 uint32
 *
 * 场: 32x32 晶格, 8 点/格 => 256x256 样本点, 周期环绕。
 * 自相关用 radix-2 FFT 的 Wiener-Khinchin 周期估计，横纵方向分开输出，
 * 滞后覆盖到 32 像素（一个完整晶格周期）。
 * 另输出高频能量比 hf: 波数 > 晶格 Nyquist(每周期16周波) 的频谱功率占比。
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "hash.h"
#include "stats.h"

#define CELLS 32
#define PPC 8
#define N (CELLS * PPC)          /* 256 */
#define LAGP 33                  /* 像素滞后 0..32 */
#define PERSEED 128

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static uint64_t (*g_hash)(int64_t, int64_t, uint64_t);
static int g_nearest;
static int g_seedscale;   /* sx 模式：格距随种子变化（C2 违反） */
static int g_octaves;

static double grid_val(int64_t i, int64_t j, uint64_t s, int oct) {
    uint64_t h = g_octaves > 1 ? ss_hash2_oct(i, j, s, oct)
                               : g_hash(i, j, s);
    return ss_u01(h) * 2.0 - 1.0;
}

static double value_noise(double x, double y, uint64_t s) {
    int64_t x0 = (int64_t)floor(x), y0 = (int64_t)floor(y);
    double fx = x - (double)x0, fy = y - (double)y0;
    if (g_nearest) {
        int64_t ix = fx < 0.5 ? x0 : x0 + 1;
        int64_t iy = fy < 0.5 ? y0 : y0 + 1;
        return grid_val(ix, iy, s, 0);
    }
    double u = fx * fx * (3.0 - 2.0 * fx);
    double v = fy * fy * (3.0 - 2.0 * fy);
    double a = grid_val(x0,     y0,     s, 0);
    double b = grid_val(x0 + 1, y0,     s, 0);
    double c = grid_val(x0,     y0 + 1, s, 0);
    double d = grid_val(x0 + 1, y0 + 1, s, 0);
    double ab = a + (b - a) * u, cd = c + (d - c) * u;
    return ab + (cd - ab) * v;
}

static double fbm_field(double x, double y, uint64_t s) {
    /* sx 模式（违反 C2）：种子低 3 位决定频率乘数 1..8，
       即"换种子"直接改变特征尺度。合格管线中 k 恒为 1。 */
    if (g_seedscale) {
        double k = (double)(1u + (s & 7u));
        x *= k;
        y *= k;
    }
    if (g_octaves <= 1) return value_noise(x, y, s);
    double sum = 0.0, amp = 1.0, norm = 0.0, freq = 1.0;
    for (int o = 0; o < g_octaves; o++) {
        sum += value_noise(x * freq, y * freq, s) * amp;
        norm += amp;
        amp *= 0.5;
        freq *= 2.0;
    }
    return sum / norm;
}

/* ---- 原地 radix-2 迭代 FFT，n 必须为 2 的幂 ---- */
static void fft(double* re, double* im, int n, int inv) {
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            double t = re[i]; re[i] = re[j]; re[j] = t;
            t = im[i]; im[i] = im[j]; im[j] = t;
        }
    }
    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2.0 * M_PI / len * (inv ? -1.0 : 1.0);
        double wr = cos(ang), wi = sin(ang);
        for (int i = 0; i < n; i += len) {
            double cr = 1.0, ci = 0.0;
            for (int k = 0; k < len / 2; k++) {
                double ur = re[i + k], ui = im[i + k];
                double vr = re[i + k + len / 2] * cr - im[i + k + len / 2] * ci;
                double vi = re[i + k + len / 2] * ci + im[i + k + len / 2] * cr;
                re[i + k] = ur + vr; im[i + k] = ui + vi;
                re[i + k + len / 2] = ur - vr; im[i + k + len / 2] = ui - vi;
                double ncr = cr * wr - ci * wi, nci = cr * wi + ci * wr;
                cr = ncr; ci = nci;
            }
        }
    }
    if (inv)
        for (int i = 0; i < n; i++) { re[i] /= n; im[i] /= n; }
}

/* 方向自相关：dir=0 沿 x（行内滞后），dir=1 沿 y。周期估计。 */
static void autocorr_axis(const double* f, double mean, double var,
                          double* out, int dir) {
    static double re[N], im[N];
    double acc[LAGP];
    memset(acc, 0, sizeof(acc));
    long lines = N;
    for (long line = 0; line < lines; line++) {
        for (int t = 0; t < N; t++) {
            long idx = dir == 0 ? line * N + t : (long)t * N + line;
            re[t] = f[idx] - mean;
            im[t] = 0.0;
        }
        fft(re, im, N, 0);
        for (int t = 0; t < N; t++) {
            double p = re[t] * re[t] + im[t] * im[t];
            re[t] = p; im[t] = 0.0;
        }
        fft(re, im, N, 1);
        for (int d = 0; d < LAGP; d++) acc[d] += re[d];
    }
    double scale = (double)N * (double)N * var;
    for (int d = 0; d < LAGP; d++) out[d] = acc[d] / scale;
}

/* 高频能量比：2D 周期谱中波数 > 晶格 Nyquist（16 周/场宽）的功率占比。
   最近邻采样的方波跳变把能量抛到高频，此指标专门抓采样端混叠。 */
static double spectral_hf(const double* f, double mean) {
    static double re[N * N], im[N * N], tr[N], ti[N];
    size_t tot = (size_t)N * N;
    for (size_t q = 0; q < tot; q++) { re[q] = f[q] - mean; im[q] = 0.0; }
    for (int y = 0; y < N; y++)
        fft(re + (size_t)y * N, im + (size_t)y * N, N, 0);
    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            tr[y] = re[(size_t)y * N + x];
            ti[y] = im[(size_t)y * N + x];
        }
        fft(tr, ti, N, 0);
        for (int y = 0; y < N; y++) {
            re[(size_t)y * N + x] = tr[y];
            im[(size_t)y * N + x] = ti[y];
        }
    }
    double all = 0.0, hi = 0.0;
    for (int fy = 0; fy < N; fy++)
        for (int fx = 0; fx < N; fx++) {
            int kx = fx <= N / 2 ? fx : fx - N;
            int ky = fy <= N / 2 ? fy : fy - N;
            double p = re[(size_t)fy * N + fx] * re[(size_t)fy * N + fx]
                     + im[(size_t)fy * N + fx] * im[(size_t)fy * N + fx];
            all += p;
            if (kx * kx + ky * ky > 16.0 * 16.0) hi += p;
        }
    return hi / all;
}

static uint64_t tier_seed(char tier, int k) {
    if (tier == 'S') return (uint64_t)k;
    if (tier == 'M') return (1ULL << 16) + (uint64_t)k * 101U;
    ss_rng_state = 0xD1B54A32D192ED03ULL;
    for (int q = 0; q < k; q++) (void)ss_rng();
    return ss_rng() & 0xFFFFFFFFULL;
}

int main(int argc, char** argv) {
    if (argc < 4) {
        fprintf(stderr,
            "usage: %s <good|bad-linear> <vn|nn|sx> <out.csv> [fbm]\n",
            argv[0]);
        return 2;
    }
    if (!strcmp(argv[1], "bad-linear")) g_hash = bad_linear;
    else g_hash = ss_hash2;
    g_nearest = !strcmp(argv[2], "nn");
    g_seedscale = !strcmp(argv[2], "sx");
    if (strcmp(argv[2], "vn") && !g_nearest && !g_seedscale) {
        fprintf(stderr, "sampling mode must be vn, nn or sx\n");
        return 2;
    }
    g_octaves = argc > 4 && !strcmp(argv[4], "fbm") ? 4 : 1;

    FILE* out = fopen(argv[3], "w");
    if (!out) { perror("fopen"); return 2; }
    fprintf(out, "tier,seed,octaves,mean,var,min,max,crossings,runlen,iso,hf");
    for (int d = 0; d < LAGP; d++) fprintf(out, ",acx%d", d);
    for (int d = 0; d < LAGP; d++) fprintf(out, ",acy%d", d);
    fprintf(out, ",xcorr_next\n");

    double* f = malloc((size_t)N * N * sizeof(double));
    double* g = malloc((size_t)N * N * sizeof(double));
    double acx[LAGP], acy[LAGP];
    const char tiers[3] = {'S', 'M', 'L'};

    for (int ti = 0; ti < 3; ti++) {
        char tname = tiers[ti];
        for (int k = 0; k < PERSEED; k++) {
            uint64_t s = tier_seed(tname, k);
            for (int py = 0; py < N; py++) {
                double gy = (double)py / PPC;
                for (int px = 0; px < N; px++) {
                    double gx = (double)px / PPC;
                    f[(size_t)py * N + px] = fbm_field(gx, gy, s);
                }
            }

            double sum = 0, lo = 1e300, hi = -1e300;
            for (size_t q = 0; q < (size_t)N * N; q++) {
                sum += f[q];
                if (f[q] < lo) lo = f[q];
                if (f[q] > hi) hi = f[q];
            }
            double mean = sum / (N * N);
            double ss = 0;
            int cross = 0;
            long runs = 0, run_acc = 0, cur = 0;
            for (size_t q = 0; q < (size_t)N * N; q++) {
                double d = f[q] - mean;
                ss += d * d;
                if (q > 0 && (f[q - 1] - mean) * d < 0) cross++;
                if (d >= 0) {
                    if (cur < 0) { runs++; run_acc += -cur; }
                    cur = cur >= 0 ? cur + 1 : 1;
                } else {
                    if (cur > 0) { runs++; run_acc += cur; }
                    cur = cur <= 0 ? cur - 1 : -1;
                }
            }
            double var = ss / (N * N);
            double runlen = runs ? (double)run_acc / runs : 0.0;

            autocorr_axis(f, mean, var, acx, 0);
            autocorr_axis(f, mean, var, acy, 1);
            double hf = spectral_hf(f, mean);
            double sx = 0, sy = 0;
            for (int d = 1; d < LAGP; d++) { sx += acx[d]; sy += acy[d]; }
            double iso = fabs(sx - sy) / (0.5 * (fabs(sx) + fabs(sy)) + 1e-300);

            /* 档内相邻种子场互相关（期望 ~0） */
            double xc = 0.0;
            if (k > 0) {
                uint64_t sp = tier_seed(tname, k - 1);
                for (int py = 0; py < N; py++)
                    for (int px = 0; px < N; px++)
                        g[(size_t)py * N + px] =
                            fbm_field((double)px / PPC, (double)py / PPC, sp);
                double mg = 0;
                for (size_t q = 0; q < (size_t)N * N; q++) mg += g[q];
                mg /= (N * N);
                double vg = 0, cov = 0;
                for (size_t q = 0; q < (size_t)N * N; q++) {
                    double a = f[q] - mean, b = g[q] - mg;
                    cov += a * b; vg += b * b;
                }
                xc = cov / sqrt(ss * vg);
            }

            fprintf(out, "%c,%llu,%d,%.8f,%.8f,%.6f,%.6f,%d,%.4f,%.6g,%.8f",
                    tname, (unsigned long long)s, g_octaves,
                    mean, var, lo, hi, cross, runlen, iso, hf);
            for (int d = 0; d < LAGP; d++) fprintf(out, ",%.8f", acx[d]);
            for (int d = 0; d < LAGP; d++) fprintf(out, ",%.8f", acy[d]);
            fprintf(out, ",%.8f\n", xc);
        }
        printf("tier=%c done, %d seeds, octaves=%d, hash=%s, sample=%s\n",
               tname, PERSEED, g_octaves, argv[1], argv[2]);
        fflush(stdout);
    }
    fclose(out);
    free(f); free(g);
    return 0;
}
