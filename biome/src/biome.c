
/* biome: width-adaptive biome transition band.
 *
 * 过渡带不再是定宽常量:带宽 w 同时决定
 *   - 模糊度 F(w):带内混合权重走 smoothstep,逐像素模糊指标 phi=4s(1-s)
 *   - 噪声查找增长 G(w):带内每像素 (w+1) 次 value_noise tap,带外 0 次
 * 距离场用 Felzenszwalb 精确欧氏距离变换,边界长度 L 由符号翻转边统计。
 *
 * 用法:
 *   ./biome test            自测(无参数等同 test)
 *   ./biome gen <outdir>    扫描 w=1..12,导出 csv/ppm/f32
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define N 256
#define PX (N * N)
#define WMAX 12
#define OCT 4

/* ---------------- value noise / fbm(与 noise 模块同一实现,带调用计数) ---------------- */

static long g_noise_calls = 0;

static unsigned hash2(int x, int y, unsigned seed) {
    unsigned h = (unsigned)x * 374761393u + (unsigned)y * 668265263u + seed * 2246822519u;
    h = (h ^ (h >> 13)) * 1274126177u;
    return h ^ (h >> 16);
}

static float rand01(unsigned h) { return (float)(h & 0xffffff) / 16777215.0f; }
static float lerp(float a, float b, float t) { return a + (b - a) * t; }
static float fade(float t) { return t * t * (3.0f - 2.0f * t); }

static float value_noise(float x, float y, unsigned seed) {
    g_noise_calls++;
    int xi = (int)floorf(x), yi = (int)floorf(y);
    float fx = x - xi, fy = y - yi;
    float v00 = rand01(hash2(xi, yi, seed));
    float v10 = rand01(hash2(xi + 1, yi, seed));
    float v01 = rand01(hash2(xi, yi + 1, seed));
    float v11 = rand01(hash2(xi + 1, yi + 1, seed));
    float u = fade(fx), v = fade(fy);
    return lerp(lerp(v00, v10, u), lerp(v01, v11, u), v) * 2.0f - 1.0f;
}

static float fbm(float x, float y, int octaves, float lacunarity, float gain, unsigned seed) {
    if (octaves <= 0) return 0.0f;
    float sum = 0.0f, amp = 1.0f, norm = 0.0f, freq = 1.0f;
    for (int i = 0; i < octaves; i++) {
        sum += value_noise(x * freq, y * freq, seed + (unsigned)i * 17u) * amp;
        norm += amp;
        amp *= gain;
        freq *= lacunarity;
    }
    return norm > 0.0f ? sum / norm : 0.0f;
}

/* ---------------- Felzenszwalb 1D/2D 精确欧氏距离变换 ---------------- */

static void edt_1d(const float* f, float* d, int* v, float* z, int n) {
    const float INF = 1e30f;
    int k = 0;
    v[0] = 0;
    z[0] = -INF;
    z[1] = INF;
    for (int q = 1; q < n; q++) {
        float s;
        do {
            s = ((f[q] + (float)q * q) - (f[v[k]] + (float)v[k] * v[k]))
                / (2.0f * (float)(q - v[k]));
            if (s <= z[k]) k--;
            else break;
        } while (k >= 0);
        k++;
        v[k] = q;
        z[k] = s;
        z[k + 1] = INF;
    }
    k = 0;
    for (int q = 0; q < n; q++) {
        while (z[k + 1] < (float)q) k++;
        float dv = (float)(q - v[k]);
        d[q] = dv * dv + f[v[k]];
    }
}

/* d2 输入为 0(边界点)/INF,输出到最近边界的平方欧氏距离 */
static void edt_2d(float* d2) {
    int* v = (int*)malloc(sizeof(int) * N);
    float* z = (float*)malloc(sizeof(float) * (N + 1));
    float* line = (float*)malloc(sizeof(float) * N);
    for (int y = 0; y < N; y++) {
        for (int x = 0; x < N; x++) line[x] = d2[y * N + x];
        edt_1d(line, line, v, z, N);
        for (int x = 0; x < N; x++) d2[y * N + x] = line[x];
    }
    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) line[y] = d2[y * N + x];
        edt_1d(line, line, v, z, N);
        for (int y = 0; y < N; y++) d2[y * N + x] = line[y];
    }
    free(line); free(z); free(v);
}

/* ---------------- 过渡带模型 ---------------- */

static float smoothstep(float t) { return t * t * (3.0f - 2.0f * t); }

/* 模糊指标:t=0/1 时为 0(纯群系),t=0.5 时为 1(完全混合) */
static float fuzz_phi(float s) { return 4.0f * s * (1.0f - s); }

static float clampf(float v, float lo, float hi) {
    return v < lo ? lo : (v > hi ? hi : v);
}

typedef struct {
    float g[PX];          /* 群系场,符号定群系 */
    float sd[PX];         /* 有符号距离(像素单位) */
    float ad[PX];         /* 排序后的 |d|,纯几何常量,供离散预测做二分 */
    int is_bnd[PX];       /* 边界单元格 */
    long L;               /* 边界长度(符号翻转边数) */
} BiomeField;

static int cmp_float(const void* a, const void* b) {
    float fa = *(const float*)a, fb = *(const float*)b;
    return (fa > fb) - (fa < fb);
}

static void field_build(BiomeField* bf, unsigned seed) {
    memset(bf, 0, sizeof(*bf));
    for (int y = 0; y < N; y++) {
        for (int x = 0; x < N; x++) {
            float v = fbm(x * 0.025f, y * 0.025f, OCT, 2.0f, 0.5f, seed);
            v += 0.35f * (0.5f - (float)y / N); /* 倾斜保证零等值线贯穿画面 */
            bf->g[y * N + x] = v;
        }
    }

    long L = 0;
    for (int y = 0; y < N; y++) {
        for (int x = 0; x < N; x++) {
            int i = y * N + x;
            int s0 = bf->g[i] >= 0.0f;
            if (x + 1 < N && (bf->g[i + 1] >= 0.0f) != s0) L++;
            if (y + 1 < N && (bf->g[i + N] >= 0.0f) != s0) L++;
            int nb_diff = 0;
            if (x > 0 && (bf->g[i - 1] >= 0.0f) != s0) nb_diff++;
            if (x + 1 < N && (bf->g[i + 1] >= 0.0f) != s0) nb_diff++;
            if (y > 0 && (bf->g[i - N] >= 0.0f) != s0) nb_diff++;
            if (y + 1 < N && (bf->g[i + N] >= 0.0f) != s0) nb_diff++;
            bf->is_bnd[i] = nb_diff > 0;
        }
    }
    bf->L = L;

    float* d2 = (float*)malloc(sizeof(float) * PX);
    for (long i = 0; i < PX; i++) d2[i] = bf->is_bnd[i] ? 0.0f : 1e30f;
    edt_2d(d2);
    for (long i = 0; i < PX; i++) {
        float d = sqrtf(d2[i]);
        bf->sd[i] = bf->g[i] >= 0.0f ? d : -d;
        bf->ad[i] = d;
    }
    free(d2);
    qsort(bf->ad, PX, sizeof(float), cmp_float);
}

/* 带宽 w 下的带内像素判定:|d| <= w/2 */
static int in_band(const BiomeField* bf, long i, int w) {
    return fabsf(bf->sd[i]) <= 0.5f * (float)w;
}

/* 解析连续预测: F_c(w) = (18/35) L w / P
 * 18/35 = ∫_0^1 4 s(t)(1-s(t)) dt,s=smoothstep。
 * 细带 w=1,2 时像素化使连续积分失准,误差在 CSV 中如实给出。 */
static double fuzz_pred_continuous(const BiomeField* bf, int w) {
    static const double CS = 18.0 / 35.0;
    return CS * (double)bf->L * (double)w / (double)PX;
}

/* 带内像素数:排序距离场上二分 |d| <= w/2,纯几何、不经过渲染代码 */
static long band_pixels_geometry(const BiomeField* bf, int w) {
    float lim = 0.5f * (float)w;
    long lo = 0, hi = PX;
    while (lo < hi) {
        long mid = lo + (hi - lo) / 2;
        if (bf->ad[mid] <= lim) lo = mid + 1; else hi = mid;
    }
    return lo;
}

/* 精确离散预测:对同一套几何常量直接逐点积分定义式
 * F_d(w) = (1/P) Σ_{|d|≤w/2} φ(smoothstep(1/2 + d/w))
 * 与渲染各自独立实现,浮点级别应一致。 */
static double fuzz_pred_discrete(const BiomeField* bf, int w) {
    double acc = 0.0;
    for (long i = 0; i < PX; i++) {
        float d = bf->ad[i];
        if (d > 0.5f * (float)w) break;
        double t = 0.5 + (double)d / (double)w;
        if (t > 1.0) t = 1.0;
        double s = t * t * (3.0 - 2.0 * t);
        acc += 4.0 * s * (1.0 - s);
    }
    return acc / (double)PX;
}

/* 带内每像素 (w+1) 次噪声 tap(1 个主 tap + w 个沿偏移的细节 tap) */
static int taps_per_band_pixel(int w) { return w + 1; }

static double lookup_pred_continuous(const BiomeField* bf, int w) {
    double B = (double)bf->L * (double)w; /* B(w) ≈ Lw */
    return (double)taps_per_band_pixel(w) * B;
}

static double lookup_pred_discrete(const BiomeField* bf, int w) {
    return (double)taps_per_band_pixel(w) * (double)band_pixels_geometry(bf, w);
}

/* ---------------- 渲染与统计 ---------------- */

typedef struct {
    long band_px;
    double fuzz_mean;
    long extra_lookups;
    unsigned char rgb[PX * 3];
    float phi[PX];
} BandResult;

static void render_band(const BiomeField* bf, int w, unsigned detail_seed, BandResult* out) {
    memset(out, 0, sizeof(*out));
    g_noise_calls = 0;
    /* A=草甸绿,B=沙地黄 */
    const float Ar = 74.0f, Ag = 124.0f, Ab = 68.0f;
    const float Br = 190.0f, Bg = 168.0f, Bb = 118.0f;

    for (int y = 0; y < N; y++) {
        for (int x = 0; x < N; x++) {
            long i = (long)y * N + x;
            float d = bf->sd[i];
            float s, phi = 0.0f, detail = 0.0f;
            if (in_band(bf, i, w)) {
                out->band_px++;
                float t = clampf(0.5f + d / (float)w, 0.0f, 1.0f);
                s = smoothstep(t);
                phi = fuzz_phi(s);
                for (int k = 0; k <= w; k++)
                    detail += value_noise((float)x + k * 0.7f, (float)y + k * 0.3f,
                                          detail_seed + (unsigned)k * 13u);
                detail /= (float)(w + 1);
            } else {
                s = d < 0.0f ? 0.0f : 1.0f;
            }
            out->phi[i] = phi;
            out->fuzz_mean += (double)phi;
            float bright = 1.0f + 0.10f * detail;
            float r = (Ar * (1.0f - s) + Br * s) * bright;
            float g = (Ag * (1.0f - s) + Bg * s) * bright;
            float b = (Ab * (1.0f - s) + Bb * s) * bright;
            long o = i * 3;
            out->rgb[o] = (unsigned char)clampf(r, 0.0f, 255.0f);
            out->rgb[o + 1] = (unsigned char)clampf(g, 0.0f, 255.0f);
            out->rgb[o + 2] = (unsigned char)clampf(b, 0.0f, 255.0f);
        }
    }
    out->fuzz_mean /= (double)PX;
    out->extra_lookups = g_noise_calls;
}

static void write_ppm(const char* path, const unsigned char* rgb) {
    FILE* fp = fopen(path, "wb");
    if (!fp) { perror(path); exit(1); }
    fprintf(fp, "P6\n%d %d\n255\n", N, N);
    fwrite(rgb, 1, (size_t)PX * 3, fp);
    fclose(fp);
}

static void write_f32(const char* path, const float* v) {
    FILE* fp = fopen(path, "wb");
    if (!fp) { perror(path); exit(1); }
    fwrite(v, sizeof(float), PX, fp);
    fclose(fp);
}

/* 最割裂区段:窄带渲染上,边界像素到四邻颜色距离的最大值 */
typedef struct { int x, y; double jump; double grad; } Worst;

static int worst_cmp(const void* a, const void* b) {
    double da = ((const Worst*)a)->jump, db = ((const Worst*)b)->jump;
    return (da < db) - (da > db);
}

static void find_worst(const BiomeField* bf, const BandResult* br, Worst* top, int ntop) {
    int cnt = 0;
    Worst* all = (Worst*)malloc(sizeof(Worst) * PX);
    for (int y = 1; y < N - 1; y++) {
        for (int x = 1; x < N - 1; x++) {
            long i = (long)y * N + x;
            if (!bf->is_bnd[i]) continue;
            const unsigned char* c = &br->rgb[i * 3];
            double jmax = 0.0;
            const long off[4] = { -1, 1, -N, N };
            for (int k = 0; k < 4; k++) {
                const unsigned char* q = &br->rgb[(i + off[k]) * 3];
                double dr = c[0] - q[0], dg = c[1] - q[1], db = c[2] - q[2];
                double j = sqrt(dr * dr + dg * dg + db * db);
                if (j > jmax) jmax = j;
            }
            float gx = (bf->g[i + 1] - bf->g[i - 1]) * 0.5f;
            float gy = (bf->g[i + N] - bf->g[i - N]) * 0.5f;
            all[cnt].x = x; all[cnt].y = y;
            all[cnt].jump = jmax;
            all[cnt].grad = sqrt(gx * gx + gy * gy);
            cnt++;
        }
    }
    qsort(all, cnt, sizeof(Worst), worst_cmp);
    for (int k = 0; k < ntop && k < cnt; k++) top[k] = all[k];
    free(all);
}

/* ---------------- 自测 ---------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static void run_tests(void) {
    /* EDT:单点边界,距离必须精确 */
    float* d2 = (float*)malloc(sizeof(float) * PX);
    for (long i = 0; i < PX; i++) d2[i] = 1e30f;
    d2[(N / 2) * N + N / 2] = 0.0f;
    edt_2d(d2);
    float want = sqrtf(3.0f * 3.0f + 4.0f * 4.0f);
    check("edt_single_site", fabsf(sqrtf(d2[((N / 2) + 4) * N + (N / 2) + 3]) - want) < 1e-4f);
    free(d2);

    BiomeField bf;
    field_build(&bf, 11u);
    check("field_both_signs", bf.L > 0);
    int sorted_ok = bf.ad[0] == 0.0f;
    for (long i = 1; i < PX && sorted_ok; i++)
        if (bf.ad[i] < bf.ad[i - 1]) sorted_ok = 0;
    check("sorted_distance_field", sorted_ok);

    /* 边界单元格有符号距离为 0 */
    int zero_ok = 1;
    for (long i = 0; i < PX; i++)
        if (bf.is_bnd[i] && fabsf(bf.sd[i]) > 1e-5f) { zero_ok = 0; break; }
    check("boundary_distance_zero", zero_ok);

    /* 带外 phi=0,带心 phi≈1 */
    BandResult br;
    render_band(&bf, 6, 77u, &br);
    int outside_ok = 1, mid_ok = 0;
    for (long i = 0; i < PX; i++) {
        if (!in_band(&bf, i, 6) && br.phi[i] != 0.0f) { outside_ok = 0; break; }
        if (bf.is_bnd[i] && fabsf(br.phi[i] - 1.0f) < 1e-6f) mid_ok = 1;
    }
    check("phi_zero_outside_band", outside_ok);
    check("phi_one_at_band_mid", mid_ok);

    /* 噪声计数精确等于 (w+1)*带内像素 */
    check("lookup_count_exact",
          br.extra_lookups == (long)taps_per_band_pixel(6) * br.band_px);

    /* 带宽变大 => 带内像素与平均模糊度单调不减 */
    BandResult a, b;
    render_band(&bf, 2, 77u, &a);
    render_band(&bf, 10, 77u, &b);
    check("band_grows_with_width", b.band_px > a.band_px);
    check("fuzz_grows_with_width", b.fuzz_mean > a.fuzz_mean);

    /* 门槛:离散预测对实测偏离 < 15%(浮点级独立实现,任意 w 适用);
     * 连续闭式解在 w>=4 的渐近区偏离 < 15%,细带误差留档不判失败 */
    int disc_ok = 1, cont_ok = 1, lookup_ok = 1;
    for (int w = 1; w <= WMAX; w++) {
        BandResult r;
        render_band(&bf, w, 77u, &r);
        double pd = fuzz_pred_discrete(&bf, w);
        double err = fabs(r.fuzz_mean - pd) / (r.fuzz_mean + 1e-12);
        if (err > 0.15) { disc_ok = 0; printf("  w=%d disc_err=%.4f\n", w, err); }
        double gd = lookup_pred_discrete(&bf, w);
        double lerr = fabs((double)r.extra_lookups - gd) / ((double)r.extra_lookups + 1e-12);
        if (lerr > 0.15) { lookup_ok = 0; printf("  w=%d lookup_err=%.4f\n", w, lerr); }
        if (w >= 4) {
            double pc = fuzz_pred_continuous(&bf, w);
            double cerr = fabs(r.fuzz_mean - pc) / (r.fuzz_mean + 1e-12);
            if (cerr > 0.15) { cont_ok = 0; printf("  w=%d cont_err=%.4f\n", w, cerr); }
        }
    }
    check("discrete_fuzz_within_15pct", disc_ok);
    check("continuous_fuzz_within_15pct_w_ge_4", cont_ok);
    check("discrete_lookup_within_15pct", lookup_ok);

    /* 确定性:同种子两次渲染一致 */
    BandResult c1, c2;
    render_band(&bf, 6, 77u, &c1);
    render_band(&bf, 6, 77u, &c2);
    check("render_deterministic",
          memcmp(c1.rgb, c2.rgb, sizeof(c1.rgb)) == 0 &&
          c1.fuzz_mean == c2.fuzz_mean);

    printf("all_ok=%d\n", failures == 0);
}

/* ---------------- 数据导出 ---------------- */

static void run_gen(const char* outdir) {
    BiomeField bf;
    long field_calls;
    g_noise_calls = 0;
    field_build(&bf, 11u);
    field_calls = g_noise_calls; /* 场生成成本与 w 无关,单列 */

    char path[512];
    snprintf(path, sizeof(path), "%s/meta.csv", outdir);
    FILE* fm = fopen(path, "w");
    if (!fm) { perror(path); exit(1); }
    fprintf(fm, "P,L,field_noise_calls,octaves,grid,wmax\n");
    fprintf(fm, "%d,%ld,%ld,%d,%d,%d\n", PX, bf.L, field_calls, OCT, N, WMAX);
    fclose(fm);

    snprintf(path, sizeof(path), "%s/sweep.csv", outdir);
    FILE* fs = fopen(path, "w");
    if (!fs) { perror(path); exit(1); }
    fprintf(fs, "w,band_px,fuzz_meas,fuzz_cont_pred,fuzz_disc_pred,"
                "fuzz_cont_rel_err,fuzz_disc_rel_err,"
                "extra_lookups,lookup_cont_pred,lookup_disc_pred,"
                "lookup_cont_rel_err,lookup_disc_rel_err\n");

    static const int tiers[3] = { 2, 6, 12 };
    for (int w = 1; w <= WMAX; w++) {
        BandResult r;
        render_band(&bf, w, 77u, &r);
        double fc = fuzz_pred_continuous(&bf, w);
        double fd = fuzz_pred_discrete(&bf, w);
        double gc = lookup_pred_continuous(&bf, w);
        double gd = lookup_pred_discrete(&bf, w);
        fprintf(fs, "%d,%ld,%.8f,%.8f,%.8f,%.6f,%.6f,%ld,%.1f,%.1f,%.6f,%.6f\n",
                w, r.band_px, r.fuzz_mean, fc, fd,
                fabs(r.fuzz_mean - fc) / (r.fuzz_mean + 1e-12),
                fabs(r.fuzz_mean - fd) / (r.fuzz_mean + 1e-12),
                r.extra_lookups, gc, gd,
                fabs((double)r.extra_lookups - gc) / ((double)r.extra_lookups + 1e-12),
                fabs((double)r.extra_lookups - gd) / ((double)r.extra_lookups + 1e-12));
        for (int k = 0; k < 3; k++) {
            if (w != tiers[k]) continue;
            snprintf(path, sizeof(path), "%s/render_w%d.ppm", outdir, w);
            write_ppm(path, r.rgb);
            snprintf(path, sizeof(path), "%s/phi_w%d.f32", outdir, w);
            write_f32(path, r.phi);
            if (w == tiers[0]) {
                Worst top[10];
                find_worst(&bf, &r, top, 10);
                char wp[512];
                snprintf(wp, sizeof(wp), "%s/worst.csv", outdir);
                FILE* fw = fopen(wp, "w");
                fprintf(fw, "rank,x,y,color_jump,g_grad,scan_width\n");
                for (int q = 0; q < 10; q++)
                    fprintf(fw, "%d,%d,%d,%.4f,%.6f,%d\n",
                            q + 1, top[q].x, top[q].y, top[q].jump, top[q].grad, w);
                fclose(fw);
            }
        }
    }
    fclose(fs);
    printf("gen_ok outdir=%s L=%ld field_calls=%ld\n", outdir, bf.L, field_calls);
}

int main(int argc, char** argv) {
    if (argc >= 3 && strcmp(argv[1], "gen") == 0) {
        run_gen(argv[2]);
        return 0;
    }
    run_tests();
    return failures == 0 ? 0 : 1;
}
