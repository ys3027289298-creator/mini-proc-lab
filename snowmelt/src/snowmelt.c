
/* snowmelt: 日照角 × 坡向差驱动的融雪修正因子与坡向非对称雪线。
 *
 * 模型（推导见 README.md），坐标 x=东 y=南 z=上:
 *   太阳方向 S = (sin A cos h, cos A cos h, sin h)   h=高度角, A=方位角(南=0)
 *   坡向 φ=下坡方位(南=0)，坡度 β；向上法向 n = (sinβ sinφ, sinβ cosφ, cosβ)
 *   入射余弦   cos i = S·n = sin h cosβ + cos h sinβ cos(A-φ)
 *   融雪修正因子 F = max(cos i, 0) / sin h           （平坦基准 F=1）
 *   解析推导式 F = max(cosβ + tan h sinβ cosΔφ, 0), Δφ = A-φ 为坡向差
 *   雪线修正   zline = z0 * (1 + k*(F-1))            南坡 F>1 雪线上移, 北坡 F<1 下移
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define N 96                 /* 高度场网格 N×N */
#define CELL 50.0f           /* 网格间距 (m) */

/* ---- 日照场景参数（北半球正午，南向 A=0）---- */
#define SUN_H   (45.0f * 3.14159265358979f / 180.0f)  /* 太阳高度角 45° */
#define SUN_A   0.0f                                   /* 太阳方位角：正南 */
#define K_ADJ   0.35f                                  /* 雪线对融雪因子的灵敏度 */
#define FMIN    0.0f                                   /* 背阴面入射截断 */

/* 三条基准雪线（地形高程分位数） */
static const float base_quant[3] = {0.45f, 0.60f, 0.75f};

typedef struct { float h[N][N]; } Field;

/* ---------- 确定性值噪声 / fBm（与 noise/ 同思路，本目录自包含）---------- */
static unsigned hash2(int x, int y, unsigned seed) {
    unsigned h = (unsigned)x * 374761393u + (unsigned)y * 668265263u + seed * 2246822519u;
    h = (h ^ (h >> 13)) * 1274126177u;
    return h ^ (h >> 16);
}
static float rand01(unsigned hh) { return (float)(hh & 0xffffff) / 16777215.0f; }
static float fade(float t) { return t * t * (3.0f - 2.0f * t); }

static float value_noise(float x, float y, unsigned seed) {
    int xi = (int)floorf(x), yi = (int)floorf(y);
    float fx = x - xi, fy = y - yi;
    float v00 = rand01(hash2(xi, yi, seed));
    float v10 = rand01(hash2(xi + 1, yi, seed));
    float v01 = rand01(hash2(xi, yi + 1, seed));
    float v11 = rand01(hash2(xi + 1, yi + 1, seed));
    float u = fade(fx), v = fade(fy);
    return ((v00 * (1 - u) + v10 * u) * (1 - v) +
            (v01 * (1 - u) + v11 * u) * v) * 2.0f - 1.0f;
}

static float fbm(float x, float y, int octaves, unsigned seed) {
    float sum = 0.0f, amp = 1.0f, norm = 0.0f, freq = 1.0f;
    for (int i = 0; i < octaves; i++) {
        sum += value_noise(x * freq, y * freq, seed + (unsigned)i * 17u) * amp;
        norm += amp;
        amp *= 0.5f;
        freq *= 2.0f;
    }
    return sum / norm;
}

/* 合成地形：大尺度脊线（制造稳定的南北坡）+ fBm 细节，抬升为 0..~2200 m */
static void terrain_build(Field* f) {
    memset(f, 0, sizeof(*f));
    for (int y = 0; y < N; y++) {
        for (int x = 0; x < N; x++) {
            float wx = (x - N / 2.0f) / (N / 2.0f);
            float wy = (y - N / 2.0f) / (N / 2.0f);
            float ridge = sinf(wy * 2.6f + 0.7f * fbm(wx * 1.5f, wy * 1.5f, 3, 23u));
            float detail = fbm(wx * 3.0f + 5.0f, wy * 3.0f - 2.0f, 5, 41u);
            float v = 0.62f * ridge + 0.38f * detail;          /* 约 -1..1 */
            f->h[x][y] = (v + 1.0f) * 0.5f * 2200.0f;          /* 0..2200 m */
        }
    }
}

/* 中心差分坡度 β、坡向 φ（φ=0 朝南，东正西负；β≈0 的平坦格点由 aspect_class 归类） */
static void slope_aspect(const Field* f, int x, int y, float* beta, float* phi) {
    int xm = x > 0 ? x - 1 : x, xp = x < N - 1 ? x + 1 : x;
    int ym = y > 0 ? y - 1 : y, yp = y < N - 1 ? y + 1 : y;
    float gx = (f->h[xp][y] - f->h[xm][y]) / ((xp - xm) * CELL);
    float gy = (f->h[x][yp] - f->h[x][ym]) / ((yp - ym) * CELL);
    *beta = atanf(sqrtf(gx * gx + gy * gy));
    /* 下坡方向朝东 = -gx、朝南 = -gy；坡向即下坡方位 */
    *phi = atan2f(-gx, -gy);
}

/* 数值版：直接由太阳向量与坡面法向点积求 F（作为推导式的对照基准） */
static float factor_numeric(float beta, float phi) {
    float sx = sinf(SUN_A) * cosf(SUN_H);
    float sy = cosf(SUN_A) * cosf(SUN_H);
    float sz = sinf(SUN_H);
    float nx = sinf(beta) * sinf(phi);   /* 东向分量 */
    float ny = sinf(beta) * cosf(phi);   /* 南向分量 */
    float nz = cosf(beta);
    float cosi = sx * nx + sy * ny + sz * nz;
    if (cosi < 0.0f) cosi = 0.0f;
    return cosi / sinf(SUN_H);
}

/* 推导版：F = max(cosβ + tan h sinβ cosΔφ, 0)，Δφ = A - φ */
static float factor_derived(float beta, float phi) {
    float dphi = SUN_A - phi;
    float v = cosf(beta) + tanf(SUN_H) * sinf(beta) * cosf(dphi);
    return v > FMIN ? v : FMIN;
}

static float quantile(const Field* f, float q) {
    float buf[N * N];
    int m = 0;
    for (int i = 0; i < N * N; i++) buf[m++] = ((const float*)f->h)[i];
    /* 简单插入排序，9216 个元素足够 */
    for (int i = 1; i < m; i++) {
        float key = buf[i];
        int j = i - 1;
        while (j >= 0 && buf[j] > key) { buf[j + 1] = buf[j]; j--; }
        buf[j + 1] = key;
    }
    float pos = q * (m - 1);
    int lo = (int)pos;
    float t = pos - lo;
    return buf[lo] * (1 - t) + buf[lo + 1] * t;
}

/* 坡向分组：南/北/东/西（|β| < 5° 视为平坦） */
typedef enum { ASP_FLAT = 0, ASP_S, ASP_N, ASP_E, ASP_W, ASP_NB } Asp;
static Asp aspect_class(float beta, float phi) {
    if (beta < 5.0f * 3.14159265f / 180.0f) return ASP_FLAT;
    float cp = cosf(phi); /* >0 朝南, <0 朝北 */
    float sp = sinf(phi); /* >0 朝东, <0 朝西 */
    if (fabsf(cp) >= fabsf(sp)) return cp > 0 ? ASP_S : ASP_N;
    return sp > 0 ? ASP_E : ASP_W;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(int argc, char** argv) {
    const char* out_dir = argc > 1 ? argv[1] : "out";
    char cmd[512];
    snprintf(cmd, sizeof(cmd), "mkdir -p '%s'", out_dir);
    if (system(cmd) != 0) { fprintf(stderr, "mkdir failed\n"); return 2; }

    Field f;
    terrain_build(&f);

    /* ---------------- 自检 ---------------- */
    const float DEG = 3.14159265f / 180.0f;

    /* 1. 平坦地面 F=1（两种实现一致） */
    float f_flat_n = factor_numeric(0.0f, 0.0f);
    float f_flat_d = factor_derived(0.0f, 0.0f);
    check("flat_F_is_one", fabsf(f_flat_n - 1.0f) < 1e-6f && fabsf(f_flat_d - 1.0f) < 1e-6f);

    /* 2. 同坡度 β=25° 南坡(φ=0)/北坡(φ=180°)：解析闭式值 */
    float b25 = 25.0f * DEG;
    float fS = cosf(b25) + tanf(SUN_H) * sinf(b25);
    float fN = cosf(b25) - tanf(SUN_H) * sinf(b25);
    check("south_gt_one", fabsf(factor_derived(b25, 0.0f) - fS) < 1e-6f && fS > 1.0f);
    check("north_lt_one", fabsf(factor_derived(b25, 180.0f * DEG) - fN) < 1e-6f && fN < 1.0f);

    /* 3. 东西坡 Δφ=90° → cosΔφ=0，F = cosβ（仅保留地面项） */
    float fE = factor_derived(b25, 90.0f * DEG);
    check("east_west_term", fabsf(fE - cosf(b25)) < 1e-6f);

    /* 4. 推导式 vs 数值点积，全场逐格点：绝对误差 < 1e-5，
     *    非背阴区(F>0.05)相对误差 < 0.5%（远低于两成门槛） */
    double max_abs = 0.0, max_rel = 0.0, sum_rel = 0.0;
    int rel_cnt = 0;
    for (int y = 0; y < N; y++) for (int x = 0; x < N; x++) {
        float beta, phi;
        slope_aspect(&f, x, y, &beta, &phi);
        float fn = factor_numeric(beta, phi), fd = factor_derived(beta, phi);
        float abs_err = fabsf(fn - fd);
        if (abs_err > max_abs) max_abs = abs_err;
        if (fd > 0.05f) {
            float rel = abs_err / fd;
            if (rel > max_rel) max_rel = rel;
            sum_rel += rel; rel_cnt++;
        }
    }
    printf("derivation_max_abs=%.2e max_rel=%.6f mean_rel=%.6f (n_rel=%d)\n",
           max_abs, max_rel, sum_rel / rel_cnt, rel_cnt);
    check("derivation_abs_tol", max_abs < 1e-5f);
    check("derivation_rel_within_half_pct", max_rel < 0.005f);

    /* 5. 截断：北坡且 β=60° 时 cos i<0 无日照，F=0 */
    float fShade = factor_derived(60.0f * DEG, 180.0f * DEG);
    check("clamp_shaded_to_zero", fShade == 0.0f);

    /* 6. 确定性：同一地形重建完全一致 */
    Field f2;
    terrain_build(&f2);
    check("terrain_deterministic", memcmp(&f, &f2, sizeof(Field)) == 0);

    /* ---------------- 非对称量统计 ---------------- */
    struct Bin { double fsum; int n; double zsum; } bin[ASP_NB];
    memset(bin, 0, sizeof(bin));
    float fgrid[N][N];
    float bgrid[N][N], pgrid[N][N];
    for (int y = 0; y < N; y++) for (int x = 0; x < N; x++) {
        float beta, phi;
        slope_aspect(&f, x, y, &beta, &phi);
        bgrid[x][y] = beta; pgrid[x][y] = phi;
        float fv = factor_derived(beta, phi);
        fgrid[x][y] = fv;
        Asp a = aspect_class(beta, phi);
        bin[a].fsum += fv; bin[a].n++; bin[a].zsum += f.h[x][y];
    }
    const char* bin_name[ASP_NB] = {"flat", "south", "north", "east", "west"};
    printf("\n[aspect-factor stats]\n");
    for (int i = 0; i < ASP_NB; i++)
        if (bin[i].n)
            printf("%s: n=%d mean_F=%.4f mean_z=%.1f\n",
                   bin_name[i], bin[i].n, (float)(bin[i].fsum / bin[i].n),
                   (float)(bin[i].zsum / bin[i].n));

    float mean_s = (float)(bin[ASP_S].fsum / bin[ASP_S].n);
    float mean_n = (float)(bin[ASP_N].fsum / bin[ASP_N].n);
    float asym_rel = (mean_s - mean_n) / (0.5f * (mean_s + mean_n));   /* 相对不对称量 */
    float asym_abs = mean_s - mean_n;
    printf("north_south_abs=%.4f north_south_rel=%.4f (%.1f%%)\n",
           asym_abs, asym_rel, 100.0 * asym_rel);
    check("asymmetry_over_20pct", asym_rel > 0.20f);

    /* 同坡度条件下的理论闭式不对称量（β=25°, h=45°），与数值统计对照 */
    float theo_s = cosf(b25) + tanf(SUN_H) * sinf(b25);
    float theo_n = fmaxf(cosf(b25) - tanf(SUN_H) * sinf(b25), 0.0f);
    float theo_rel = (theo_s - theo_n) / (0.5f * (theo_s + theo_n));
    printf("theoretical_beta25_rel=%.4f (%.1f%%)  F_s=%.4f F_n=%.4f\n",
           theo_rel, 100.0 * theo_rel, theo_s, theo_n);

    /* 雪线南北位移（取中间基准 z0） */
    float z0 = quantile(&f, base_quant[1]);
    float line_s = z0 * (1 + K_ADJ * (mean_s - 1));
    float line_n = z0 * (1 + K_ADJ * (mean_n - 1));
    printf("snowline_shift: south=%.1f m north=%.1f m gap=%.1f m (z0=%.1f m)\n",
           line_s, line_n, line_s - line_n, z0);

    /* ---------------- 三条雪线渲染数据输出 ---------------- */
    float zbase[3];
    for (int t = 0; t < 3; t++) zbase[t] = quantile(&f, base_quant[t]);

    char path[600];
    snprintf(path, sizeof(path), "%s/stats.txt", out_dir);
    FILE* meta = fopen(path, "w");
    fprintf(meta, "N=%d cell=%.1f sun_h_deg=45 sun_az_deg=0(south) K=%.2f\n", N, CELL, K_ADJ);
    fprintf(meta, "derivation_max_abs=%.2e max_rel=%.6f mean_rel=%.6f\n",
            max_abs, max_rel, sum_rel / rel_cnt);
    for (int i = 0; i < ASP_NB; i++)
        if (bin[i].n)
            fprintf(meta, "%s n=%d mean_F=%.4f mean_z=%.1f\n",
                    bin_name[i], bin[i].n, (float)(bin[i].fsum / bin[i].n),
                    (float)(bin[i].zsum / bin[i].n));
    fprintf(meta, "north_south_abs=%.4f north_south_rel=%.4f\n", asym_abs, asym_rel);
    fprintf(meta, "theoretical_beta25_rel=%.4f F_s=%.4f F_n=%.4f\n", theo_rel, theo_s, theo_n);
    for (int t = 0; t < 3; t++)
        fprintf(meta, "baseline%d quant=%.2f z0=%.1f\n", t, base_quant[t], zbase[t]);
    fclose(meta);

    /* 二进制网格：float32，N×N，按 y 行 x 列行主序写出 */
    snprintf(path, sizeof(path), "%s/height.f32", out_dir);
    FILE* fh = fopen(path, "wb");
    for (int y = 0; y < N; y++) for (int x = 0; x < N; x++)
        fwrite(&f.h[x][y], sizeof(float), 1, fh);
    fclose(fh);

    snprintf(path, sizeof(path), "%s/factor.f32", out_dir);
    FILE* ff = fopen(path, "wb");
    for (int y = 0; y < N; y++) for (int x = 0; x < N; x++)
        fwrite(&fgrid[x][y], sizeof(float), 1, ff);
    fclose(ff);

    snprintf(path, sizeof(path), "%s/slope.f32", out_dir);
    FILE* fs = fopen(path, "wb");
    for (int y = 0; y < N; y++) for (int x = 0; x < N; x++)
        fwrite(&bgrid[x][y], sizeof(float), 1, fs);
    fclose(fs);

    snprintf(path, sizeof(path), "%s/aspect.f32", out_dir);
    FILE* fa = fopen(path, "wb");
    for (int y = 0; y < N; y++) for (int x = 0; x < N; x++)
        fwrite(&pgrid[x][y], sizeof(float), 1, fa);
    fclose(fa);

    for (int t = 0; t < 3; t++) {
        /* 基准雪线掩膜：h >= z0 */
        snprintf(path, sizeof(path), "%s/snow_base%d.f32", out_dir, t);
        FILE* o = fopen(path, "wb");
        float one = 1.0f, zero = 0.0f;
        for (int y = 0; y < N; y++) for (int x = 0; x < N; x++)
            fwrite(f.h[x][y] >= zbase[t] ? &one : &zero, sizeof(float), 1, o);
        fclose(o);

        /* 修正雪线掩膜：h >= z0*(1+K(F-1)) */
        snprintf(path, sizeof(path), "%s/snow_adj%d.f32", out_dir, t);
        o = fopen(path, "wb");
        for (int y = 0; y < N; y++) for (int x = 0; x < N; x++) {
            float line = zbase[t] * (1.0f + K_ADJ * (fgrid[x][y] - 1.0f));
            fwrite(f.h[x][y] >= line ? &one : &zero, sizeof(float), 1, o);
        }
        fclose(o);
    }
    printf("\nwrote grids + stats to %s/\n", out_dir);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
