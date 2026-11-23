/* river: D8 flow direction from slope-aspect, confluence merge, reverse-slope audit.
 *
 * 管道: value-noise/fBm 高度场 -> priority-flood 填洼 -> D8 最陡下降流向
 *       -> 流量累积 -> 汇流点(入度>=2) -> 阈值提取河道 -> SVG 河道图
 * 自检: make test ; 三档出图: make generate
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define W 80
#define H 80
#define NCELL (W * H)
#define EPS 1e-4f            /* 填洼抬升步长: 保证下游严格低于上游 */
#define RIVER_THRESH 100     /* 汇水面积(格点数)达到该值即为河道 */

/* ---------- value noise / fBm (与 noise/ 模块同一套确定性哈希) ---------- */

static unsigned hash2(int x, int y, unsigned seed) {
    unsigned h = (unsigned)x * 374761393u + (unsigned)y * 668265263u + seed * 2246822519u;
    h = (h ^ (h >> 13)) * 1274126177u;
    return h ^ (h >> 16);
}

static float rand01(unsigned h) { return (float)(h & 0xffffff) / 16777215.0f; }
static float lerp(float a, float b, float t) { return a + (b - a) * t; }
static float fade(float t) { return t * t * (3.0f - 2.0f * t); }

static float value_noise(float x, float y, unsigned seed) {
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

/* ---------- 栅格与 D8 罗盘 ---------- */

/* 固定扫描顺序即并列时的裁决顺序, 保证可复现: E, SE, S, SW, W, NW, N, NE */
static const int DX[8] = {1, 1, 0, -1, -1, -1, 0, 1};
static const int DY[8] = {0, 1, 1, 1, 0, -1, -1, -1};
static const char *DNAME[8] = {"E", "SE", "S", "SW", "W", "NW", "N", "NE"};
static const char *dname(int d) { return d >= 0 && d < 8 ? DNAME[d] : "OUT"; }

typedef struct {
    float h[NCELL];      /* 填洼后的高程 */
    float raw[NCELL];    /* 原始噪声高程 */
    int   dir[NCELL];    /* D8 流向, 0..7; 边界出口为 -1 */
    int   acc[NCELL];    /* 汇水面积(含自身) */
    int   indeg[NCELL];  /* 流入邻居数 */
} Grid;

static int idx(int x, int y) { return y * W + x; }
static int inside(int x, int y) { return x >= 0 && y >= 0 && x < W && y < H; }
static int is_border(int x, int y) { return x == 0 || y == 0 || x == W - 1 || y == H - 1; }

/* ---------- 三档坡度高度场 ---------- */

typedef struct {
    const char *key, *name;
    int octaves;
    float amp, freq, tilt;
} Tier;

/* 同一基础倾斜 + 不同 fBm 层数/振幅/频率: 噪声越重, 坡度方差越大。
 * 频率均取到 3.5 以上, 保证局部起伏能持续拨转流向, 避免长距离大直道。 */
static const Tier TIERS[3] = {
    {"gentle", "缓坡(3 层 fBm, 弱起伏)", 3,  5.0f, 3.5f, 0.12f},
    {"mid",    "中坡(4 层 fBm, 中起伏)", 4,  8.0f, 4.5f, 0.12f},
    {"steep",  "陡坡(6 层 fBm, 强起伏)", 6, 14.0f, 5.5f, 0.12f},
};

static void terrain_build(Grid* g, const Tier* t, unsigned seed) {
    memset(g, 0, sizeof(*g));
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            float n = fbm((float)x * t->freq / (float)W,
                          (float)y * t->freq / (float)H,
                          t->octaves, 2.0f, 0.5f, seed);
            float h = -t->tilt * (float)y + t->amp * n; /* -y 倾斜: 水朝 +y 流 */
            g->raw[idx(x, y)] = h;
            g->h[idx(x, y)] = h;
        }
    }
}

/* ---------- 最小二叉堆 (priority-flood 用, 键为高程) ---------- */

typedef struct { int c[NCELL]; int n; } Heap;

static void heap_push(Heap* q, const float* key, int c) {
    int i = q->n++;
    q->c[i] = c;
    while (i > 0) {
        int p = (i - 1) / 2;
        if (key[q->c[p]] <= key[q->c[i]]) break;
        int tmp = q->c[p]; q->c[p] = q->c[i]; q->c[i] = tmp;
        i = p;
    }
}

static int heap_pop(Heap* q, const float* key) {
    int top = q->c[0];
    q->c[0] = q->c[--q->n];
    int i = 0;
    for (;;) {
        int l = i * 2 + 1, r = l + 1, s = i;
        if (l < q->n && key[q->c[l]] < key[q->c[s]]) s = l;
        if (r < q->n && key[q->c[r]] < key[q->c[s]]) s = r;
        if (s == i) break;
        int tmp = q->c[s]; q->c[s] = q->c[i]; q->c[i] = tmp;
        i = s;
    }
    return top;
}

/* priority-flood (Barnes et al. 2014, epsilon 版):
 * 从边界向内泛洪, 低洼处被抬升到"入口高程 + EPS",
 * 于是每个内部格点都存在严格更低的下游, 且不制造新的分水岭。 */
static int priority_flood(Grid* g) {
    Heap q; q.n = 0;
    unsigned char* done = calloc(NCELL, 1);
    int filled = 0;
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++)
            if (is_border(x, y)) {
                int c = idx(x, y);
                done[c] = 1;
                heap_push(&q, g->h, c);
            }
    while (q.n > 0) {
        int c = heap_pop(&q, g->h);
        int cx = c % W, cy = c / W;
        for (int d = 0; d < 8; d++) {
            int nx = cx + DX[d], ny = cy + DY[d];
            if (!inside(nx, ny)) continue;
            int nc = idx(nx, ny);
            if (done[nc]) continue;
            done[nc] = 1;
            if (g->h[nc] <= g->h[c]) {        /* 洼地/平地: 抬到刚好低于?不,高于 c 一个 EPS */
                g->h[nc] = g->h[c] + EPS;
                filled++;
            }
            heap_push(&q, g->h, nc);
        }
    }
    free(done);
    return filled;
}

/* ---------- 坡度与流向 ---------- */

/* 中心差分坡度(无量纲, cell=1): |∇h|; 用于三档坡度统计 */
static void slope_stats(const Grid* g, float* mean_out, float* max_out) {
    double sum = 0.0; float mx = 0.0f;
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            int xm = x > 0 ? x - 1 : x, xp = x < W - 1 ? x + 1 : x;
            int ym = y > 0 ? y - 1 : y, yp = y < H - 1 ? y + 1 : y;
            float dx = (g->h[idx(xp, y)] - g->h[idx(xm, y)]) / (float)(xp - xm);
            float dy = (g->h[idx(x, yp)] - g->h[idx(x, ym)]) / (float)(yp - ym);
            float s = sqrtf(dx * dx + dy * dy);
            sum += s;
            if (s > mx) mx = s;
        }
    }
    *mean_out = (float)(sum / NCELL);
    *max_out = mx;
}

/* D8: 流向 = 坡降方向上最陡的 8 邻居。
 * 推导: 连续场中下降方向为 -∇h; 离散格网用 8 方向中 (h_c - h_n)/L_n 最大者逼近,
 * 对角邻居 L=sqrt(2), 正交邻居 L=1。只接受严格下降; 边界无更低邻居时为出口(-1)。
 * 并列时按罗盘固定顺序取先出现者。 */
static void flow_d8(Grid* g) {
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            int c = idx(x, y);
            g->dir[c] = -1;
            float best = 0.0f;
            for (int d = 0; d < 8; d++) {
                int nx = x + DX[d], ny = y + DY[d];
                if (!inside(nx, ny)) continue;
                int nc = idx(nx, ny);
                float L = (DX[d] == 0 || DY[d] == 0) ? 1.0f : 1.41421356f;
                float drop = (g->h[c] - g->h[nc]) / L;
                if (drop > best) { best = drop; g->dir[c] = d; }
            }
        }
    }
}

typedef struct { int from, to; } BadEdge;

/* 逆坡审计: 沿每条流向边核验 h_to < h_from; 内部格点无下游视为汇水死角。
 * 返回问题总数, 前 cap 条写入 out 并区分 kind: 1=逆坡边 2=内部死角。
 * 每条记录的坐标/高程/流向均可直接用于定位报告。 */
static int reverse_audit(const Grid* g, BadEdge* out, int* kind, int cap) {
    int bad = 0;
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            int c = idx(x, y);
            if (g->dir[c] < 0) {
                if (!is_border(x, y)) {
                    if (bad < cap) { out[bad].from = c; out[bad].to = -1; kind[bad] = 2; }
                    bad++;
                }
                continue;
            }
            int nc = idx(x + DX[g->dir[c]], y + DY[g->dir[c]]);
            if (!(g->h[nc] < g->h[c])) {
                if (bad < cap) { out[bad].from = c; out[bad].to = nc; kind[bad] = 1; }
                bad++;
            }
        }
    }
    return bad;
}

/* 故障注入(演示用): 把内部格点 (fx,fy) 的流向强行改为强制罗盘方向,
 * 用于验证逆坡审计能否判定"未完成"并给出位置。返回 0=注入成功。 */
static int inject_fault(Grid* g, int fx, int fy, int force_dir) {
    if (!inside(fx, fy) || is_border(fx, fy) || force_dir < 0 || force_dir >= 8) return -1;
    int nx = fx + DX[force_dir], ny = fy + DY[force_dir];
    if (!inside(nx, ny)) return -1;
    g->dir[idx(fx, fy)] = force_dir;
    return 0;
}

/* 在 (x,y) 的 8 邻居中找高程最高者作为上坡方向; 该邻居不高于自身(局部峰顶)时返回 -1 */
static int pick_uphill_dir(const Grid* g, int x, int y) {
    int c = idx(x, y), best = -1;
    float best_gain = 0.0f;
    for (int d = 0; d < 8; d++) {
        int nx = x + DX[d], ny = y + DY[d];
        if (!inside(nx, ny)) continue;
        float gain = g->h[idx(nx, ny)] - g->h[c];
        if (gain > best_gain) { best_gain = gain; best = d; }
    }
    return best;
}

/* ---------- 汇水面积与汇流点 ---------- */

static int cmp_float_desc(const void* a, const void* b, void* arg) {
    const float* h = arg;
    float ha = h[*(const int*)a], hb = h[*(const int*)b];
    return (ha < hb) - (ha > hb);
}

static void accumulate(Grid* g) {
    int* order = malloc(NCELL * sizeof(int));
    for (int i = 0; i < NCELL; i++) { order[i] = i; g->acc[i] = 1; g->indeg[i] = 0; }
    qsort_r(order, NCELL, sizeof(int), cmp_float_desc, g->h);
    for (int i = 0; i < NCELL; i++) {
        int c = order[i];
        if (g->dir[c] >= 0) {
            int to = idx(c % W + DX[g->dir[c]], c / W + DY[g->dir[c]]);
            g->acc[to] += g->acc[c];
            g->indeg[to]++;
        }
    }
    free(order);
}

/* 汇流合并条件: 格点入度 >= 2, 即两条及以上独立支流在该点并入同一下游链路。
 * 入度 1 只是河道延伸, 不计汇流。 */
static int count_confluences(const Grid* g, int* out_cells, int cap) {
    int n = 0;
    for (int i = 0; i < NCELL; i++)
        if (g->indeg[i] >= 2 && n < cap) out_cells[n++] = i;
    return n;
}

/* ---------- 河道形态统计 (防"大直道": 方向多样性与最长直线段) ---------- */

static void river_shape(const Grid* g, int* dirs_used, int* max_run, int* river_cells) {
    int used = 0, cells = 0;
    for (int i = 0; i < NCELL; i++) if (g->acc[i] >= RIVER_THRESH) cells++;
    for (int d = 0; d < 8; d++) {
        int hit = 0;
        for (int i = 0; i < NCELL; i++)
            if (g->acc[i] >= RIVER_THRESH && g->dir[i] == d) { hit = 1; break; }
        if (hit) used |= 1 << d;
    }
    int worst = 0;
    for (int s = 0; s < NCELL; s++) {
        if (g->acc[s] < RIVER_THRESH || g->dir[s] < 0) continue;
        int run = 0, c = s, d = g->dir[s];
        while (g->acc[c] >= RIVER_THRESH && g->dir[c] == d) {
            run++;
            c = idx(c % W + DX[d], c / W + DY[d]);
            if (!inside(c % W, c / W)) break;
        }
        if (run > worst) worst = run;
    }
    *dirs_used = 0;
    for (int d = 0; d < 8; d++) if (used & (1 << d)) (*dirs_used)++;
    *max_run = worst;
    *river_cells = cells;
}

/* ---------- ASCII 河道图 (终端速览; 与 SVG 同口径, 只标河道网络上的汇流点) ---------- */

static void ascii_map(const Grid* g, const int* conf, int nconf, FILE* fp) {
    unsigned char* isconf = calloc(NCELL, 1);
    for (int i = 0; i < nconf; i++) isconf[conf[i]] = 1;
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            int c = idx(x, y);
            char ch = '.';
            if (g->acc[c] >= RIVER_THRESH) ch = '#';
            if (isconf[c]) ch = 'o';
            if (is_border(x, y) && g->dir[c] < 0) ch = '~';
            fputc(ch, fp);
        }
        fputc('\n', fp);
    }
    free(isconf);
}

/* ---------- SVG 河道图 ---------- */

static int hex_shade(const Grid* g, int c, float lo, float hi) {
    float t = (g->h[c] - lo) / (hi - lo + 1e-12f);
    if (t < 0.0f) t = 0.0f;
    if (t > 1.0f) t = 1.0f;
    /* 中性灰阶: 低 #d9d7d1 -> 高 #8f8c84 */
    int lo_c[3] = {217, 215, 209}, hi_c[3] = {143, 140, 132};
    int v[3];
    for (int k = 0; k < 3; k++) v[k] = lo_c[k] + (int)((hi_c[k] - lo_c[k]) * t);
    return (v[0] << 16) | (v[1] << 8) | v[2];
}

static void write_svg(const Grid* g, const Tier* t, const char* path,
                      int filled, float mean_s, float max_s, int nconf) {
    FILE* fp = fopen(path, "w");
    if (!fp) { perror(path); exit(3); }
    const int cell = 8, pad = 24, title_h = 64, legend_h = 30;
    int mw = W * cell, mh = H * cell;
    int fw = mw + pad * 2, fh = mh + pad * 2 + title_h + legend_h;
    float lo = g->h[0], hi = g->h[0];
    for (int i = 1; i < NCELL; i++) { if (g->h[i] < lo) lo = g->h[i]; if (g->h[i] > hi) hi = g->h[i]; }

    fprintf(fp, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\" viewBox=\"0 0 %d %d\" font-family=\"system-ui,-apple-system,'Segoe UI',sans-serif\">\n", fw, fh, fw, fh);
    fprintf(fp, "<rect width=\"100%%\" height=\"100%%\" fill=\"#fcfcfb\"/>\n");
    fprintf(fp, "<text x=\"%d\" y=\"30\" font-size=\"18\" font-weight=\"600\" fill=\"#0b0b0b\">河道网络 · %s</text>\n", pad, t->name);
    fprintf(fp, "<text x=\"%d\" y=\"50\" font-size=\"12\" fill=\"#52514e\">D8 最陡下降流向 · 汇水阈值 %d 格 · 平均坡度 %.3f · 最大坡度 %.3f · 填洼 %d 格 · 汇流点 %d 个</text>\n",
            pad, RIVER_THRESH, mean_s, max_s, filled, nconf);

    int top = pad + title_h;
    /* 高程底图 */
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++) {
            int rgb = hex_shade(g, idx(x, y), lo, hi);
            fprintf(fp, "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" fill=\"#%06x\"/>\n",
                    pad + x * cell, top + y * cell, cell, cell, rgb);
        }
    /* 河道: 先描表面色隔离环, 再画蓝色线 (2px 表面间隙规则) */
    for (int pass = 0; pass < 2; pass++) {
        const char* col = pass == 0 ? "#fcfcfb" : "#3987e5";
        float wd = pass == 0 ? 3.6f : 2.2f;
        for (int y = 0; y < H; y++)
            for (int x = 0; x < W; x++) {
                int c = idx(x, y);
                if (g->acc[c] < RIVER_THRESH || g->dir[c] < 0) continue;
                int nx = x + DX[g->dir[c]], ny = y + DY[g->dir[c]];
                fprintf(fp, "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"%s\" stroke-width=\"%.1f\" stroke-linecap=\"round\"/>\n",
                        pad + (x + 0.5f) * cell, top + (y + 0.5f) * cell,
                        pad + (nx + 0.5f) * cell, top + (ny + 0.5f) * cell, col, wd);
            }
    }
    /* 汇流点(仅实质性汇流): 橙色描边圆环, 不填充, 蓝色河道保持为视觉主体 */
    for (int i = 0; i < NCELL; i++) {
        if (g->indeg[i] < 2 || g->acc[i] < RIVER_THRESH) continue;
        int strong = 0;
        for (int d = 0; d < 8; d++) {
            int nx = i % W + DX[d], ny = i / W + DY[d];
            if (!inside(nx, ny)) continue;
            /* 邻居 n 沿 d 方向抵达 c, 即 n 的流向为 d 的反方向 (d+4)%8 */
            if (g->dir[idx(nx, ny)] == (d + 4) % 8 && g->acc[idx(nx, ny)] >= RIVER_THRESH) strong++;
        }
        if (strong < 2) continue;
        int x = i % W, y = i / W;
        fprintf(fp, "<circle cx=\"%d\" cy=\"%d\" r=\"4.2\" fill=\"#fcfcfb\" stroke=\"#d95926\" stroke-width=\"2\"/>\n",
                pad + x * cell + cell / 2, top + y * cell + cell / 2);
    }
    /* 图例 */
    int ly = top + mh + pad + 14;
    fprintf(fp, "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"#3987e5\" stroke-width=\"2.2\" stroke-linecap=\"round\"/>\n", pad, ly, pad + 26, ly);
    fprintf(fp, "<text x=\"%d\" y=\"%d\" font-size=\"12\" fill=\"#52514e\">河道(汇水面积≥%d)</text>\n", pad + 32, ly + 4, RIVER_THRESH);
    fprintf(fp, "<circle cx=\"%d\" cy=\"%d\" r=\"4.2\" fill=\"#fcfcfb\" stroke=\"#d95926\" stroke-width=\"2\"/>\n", pad + 170, ly);
    fprintf(fp, "<text x=\"%d\" y=\"%d\" font-size=\"12\" fill=\"#52514e\">实质汇流点(两条入流均成河)</text>\n", pad + 180, ly + 4);
    fprintf(fp, "<rect x=\"%d\" y=\"%d\" width=\"12\" height=\"12\" fill=\"#d9d7d1\" stroke=\"#c3c2b7\"/>\n", pad + 320, ly - 9);
    fprintf(fp, "<rect x=\"%d\" y=\"%d\" width=\"12\" height=\"12\" fill=\"#8f8c84\" stroke=\"#c3c2b7\"/>\n", pad + 334, ly - 9);
    fprintf(fp, "<text x=\"%d\" y=\"%d\" font-size=\"12\" fill=\"#52514e\">高程 低→高</text>\n", pad + 352, ly + 4);
    fprintf(fp, "</svg>\n");
    fclose(fp);
}

/* ---------- 一档完整流程 ---------- */

typedef struct {
    int filled, nconf, nconf_draw, river_cells, dirs_used, max_run, reverse_bad, interior_sinks;
    float mean_s, max_s;
    BadEdge bad[16]; int bad_kind[16], bad_total;
} TierStat;

static TierStat run_tier(const Tier* t, unsigned seed, const char* outdir, int write_files,
                         int fx, int fy, int force_dir) {
    Grid g;
    terrain_build(&g, t, seed);
    TierStat st;
    st.filled = priority_flood(&g);
    slope_stats(&g, &st.mean_s, &st.max_s);
    flow_d8(&g);
    if (fx >= 0) inject_fault(&g, fx, fy, force_dir);

    st.bad_total = reverse_audit(&g, st.bad, st.bad_kind, 16);
    st.reverse_bad = 0; st.interior_sinks = 0;
    for (int i = 0; i < st.bad_total && i < 16; i++) {
        if (st.bad_kind[i] == 1) st.reverse_bad++;
        if (st.bad_kind[i] == 2) st.interior_sinks++;
    }
    if (st.bad_total > 16) st.reverse_bad += st.bad_total - 16; /* 超出记录上限的兜底计数 */

    accumulate(&g);
    int all_conf[NCELL];
    st.nconf = count_confluences(&g, all_conf, NCELL);
    /* 制图只标实质性汇流: 入度>=2 且至少两条入流各自的汇水面积都达到河道阈值
     * (排除只有细沟并入的小汇合); 先反查每个格点的入流邻居。 */
    int conf[NCELL], n_draw = 0;
    for (int k = 0; k < st.nconf; k++) {
        int c = all_conf[k], strong = 0;
        for (int d = 0; d < 8; d++) {
            int nx = c % W + DX[d], ny = c / W + DY[d];
            if (!inside(nx, ny)) continue;
            int nc = idx(nx, ny);
            if (g.dir[nc] == (d + 4) % 8 && g.acc[nc] >= RIVER_THRESH) strong++;
        }
        if (strong >= 2 && g.acc[c] >= RIVER_THRESH) conf[n_draw++] = c;
    }
    st.nconf_draw = n_draw;
    river_shape(&g, &st.dirs_used, &st.max_run, &st.river_cells);

    if (write_files) {
        char path[512];
        snprintf(path, sizeof(path), "%s/river_%s.svg", outdir, t->key);
        write_svg(&g, t, path, st.filled, st.mean_s, st.max_s, st.nconf_draw);
        snprintf(path, sizeof(path), "%s/river_%s.txt", outdir, t->key);
        FILE* fp = fopen(path, "w");
        if (!fp) { perror(path); exit(3); }
        fprintf(fp, "# %s  '#'=河道 'o'=实质汇流点(两条入流均成河) '~'=边界出口 '.'=其他 (全图入度>=2 共 %d 个)\n",
                t->name, st.nconf);
        ascii_map(&g, conf, n_draw, fp);
        fclose(fp);
    }
    return st;
}

/* ---------- 自检 ---------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* 测试用可配置尺寸网格, 与主网格同算法的独立小实现 */
#define TW 7
#define TH 7
#define TN (TW * TH)

static int tidx(int x, int y) { return y * TW + x; }

static void tiny_heap_push(int* q, int* n, const float* key, int c) {
    int i = (*n)++; q[i] = c;
    while (i > 0) {
        int p = (i - 1) / 2;
        if (key[q[p]] <= key[q[i]]) break;
        int tmp = q[p]; q[p] = q[i]; q[i] = tmp; i = p;
    }
}
static int tiny_heap_pop(int* q, int* n, const float* key) {
    int top = q[0]; q[0] = q[--(*n)];
    int i = 0;
    for (;;) {
        int l = i * 2 + 1, r = l + 1, s = i;
        if (l < *n && key[q[l]] < key[q[s]]) s = l;
        if (r < *n && key[q[r]] < key[q[s]]) s = r;
        if (s == i) break;
        int tmp = q[s]; q[s] = q[i]; q[i] = tmp; i = s;
    }
    return top;
}

static int tiny_flood(float* h) {
    int q[TN], qn = 0, filled = 0;
    unsigned char done[TN] = {0};
    for (int y = 0; y < TH; y++)
        for (int x = 0; x < TW; x++)
            if (x == 0 || y == 0 || x == TW - 1 || y == TH - 1) {
                done[tidx(x, y)] = 1; tiny_heap_push(q, &qn, h, tidx(x, y));
            }
    while (qn > 0) {
        int c = tiny_heap_pop(q, &qn, h), cx = c % TW, cy = c / TW;
        for (int d = 0; d < 8; d++) {
            int nx = cx + DX[d], ny = cy + DY[d];
            if (nx < 0 || ny < 0 || nx >= TW || ny >= TH) continue;
            int nc = tidx(nx, ny);
            if (done[nc]) continue;
            done[nc] = 1;
            if (h[nc] <= h[c]) { h[nc] = h[c] + EPS; filled++; }
            tiny_heap_push(q, &qn, h, nc);
        }
    }
    return filled;
}

static void tiny_d8(const float* h, int* dir) {
    for (int y = 0; y < TH; y++)
        for (int x = 0; x < TW; x++) {
            int c = tidx(x, y); dir[c] = -1; float best = 0.0f;
            for (int d = 0; d < 8; d++) {
                int nx = x + DX[d], ny = y + DY[d];
                if (nx < 0 || ny < 0 || nx >= TW || ny >= TH) continue;
                float L = (DX[d] == 0 || DY[d] == 0) ? 1.0f : 1.41421356f;
                float drop = (h[c] - h[tidx(nx, ny)]) / L;
                if (drop > best) { best = drop; dir[c] = d; }
            }
        }
}

static void tests(void) {
    /* 1. 含坑地形: 填洼后无内部死角, 所有流向边严格下降 */
    float h[TN];
    for (int y = 0; y < TH; y++)
        for (int x = 0; x < TW; x++)
            h[tidx(x, y)] = 5.0f - 0.1f * y + 0.05f * (float)(x - 3);
    h[tidx(3, 3)] = -9.0f; /* 深坑 */
    int filled = tiny_flood(h);
    check("pit_gets_filled", filled >= 1);
    int dir[TN];
    tiny_d8(h, dir);
    int sinks = 0, uphill = 0;
    for (int y = 0; y < TH; y++)
        for (int x = 0; x < TW; x++) {
            int c = tidx(x, y);
            if (dir[c] < 0) {
                if (!(x == 0 || y == 0 || x == TW - 1 || y == TH - 1)) sinks++;
                continue;
            }
            int nc = tidx(x + DX[dir[c]], y + DY[dir[c]]);
            if (!(h[nc] < h[c])) uphill++;
        }
    check("no_interior_sink_after_flood", sinks == 0);
    check("no_uphill_edge_after_flood", uphill == 0);

    /* 2. D8 最陡方向: 平面 h = -0.1*y (朝 +y 单向倾斜) -> 全部流向 S */
    float plane[TN]; int pd[TN];
    for (int i = 0; i < TN; i++) plane[i] = 0.0f;
    for (int y = 0; y < TH; y++)
        for (int x = 0; x < TW; x++)
            plane[tidx(x, y)] = -0.1f * y;
    tiny_d8(plane, pd);
    int all_south = 1;
    for (int y = 0; y < TH - 1; y++)
        for (int x = 0; x < TW; x++)
            if (pd[tidx(x, y)] != 2) all_south = 0;
    check("d8_plane_points_south", all_south);

    /* 3. 对角最陡: h 沿 SE 对角下降, 中心格点应选 SE(1) */
    float dia[TN]; int dd[TN];
    for (int y = 0; y < TH; y++)
        for (int x = 0; x < TW; x++)
            dia[tidx(x, y)] = -0.2f * (x + y);
    tiny_d8(dia, dd);
    check("d8_diagonal_pick", dd[tidx(3, 3)] == 1);

    /* 4. 逆坡审计能定位人为制造的逆坡边 */
    float bad_h[TN]; int bad_d[TN];
    for (int i = 0; i < TN; i++) { bad_h[i] = plane[i]; bad_d[i] = -1; }
    tiny_d8(bad_h, bad_d);
    int bx = 3, by = 2;
    bad_d[tidx(bx, by)] = 6; /* 强行指向 N: 那里高程更高 -> 逆坡 */
    int found = 0, located = -1;
    for (int y = 0; y < TH; y++)
        for (int x = 0; x < TW; x++) {
            int c = tidx(x, y);
            if (bad_d[c] < 0) continue;
            int nc = tidx(x + DX[bad_d[c]], y + DY[bad_d[c]]);
            if (!(bad_h[nc] < bad_h[c])) { found++; located = c; }
        }
    check("reverse_slope_detected", found == 1 && located == tidx(bx, by));
    if (located >= 0)
        printf("  reverse located at cell (%d,%d) forced dir %s (upstream h=%.2f >= downstream h=%.2f)\n",
               located % TW, located / TW, dname(bad_d[located]),
               bad_h[located], bad_h[tidx(bx + DX[bad_d[located]], by + DY[bad_d[located]])]);

    /* 5. 汇流条件: Y 型双支流地形(天然无坑, 不需填洼)。
     *    在向南倾斜的基底上刻两条严格单调降低的沟槽, 于 (3,3) 汇合:
     *    左支 (2,1)->(2,2)->(3,3), 右支 (4,1)->(4,2)->(3,3),
     *    合流后 (3,3)->(3,4)->(3,5)-> 边界出口 (3,6)。沟槽全程最低且单调,
     *    priority-flood 不应改动任何格点。 */
    float y[TN]; int yd[TN];
    for (int yy = 0; yy < TH; yy++)
        for (int xx = 0; xx < TW; xx++)
            y[tidx(xx, yy)] = -(float)yy;
    const int carve[][2] = {{2,1},{2,2},{4,1},{4,2},{3,3},{3,4},{3,5},{3,6}};
    for (unsigned i = 0; i < sizeof(carve) / sizeof(carve[0]); i++)
        y[tidx(carve[i][0], carve[i][1])] -= 1.5f;
    int y_filled = tiny_flood(y);
    tiny_d8(y, yd);
    int y_indeg[TN] = {0};
    for (int yy = 0; yy < TH; yy++)
        for (int xx = 0; xx < TW; xx++) {
            int c = tidx(xx, yy);
            if (yd[c] >= 0) y_indeg[tidx(xx + DX[yd[c]], yy + DY[yd[c]])]++;
        }
    check("y_shape_needs_no_fill", y_filled == 0);
    if (getenv("RIVER_DEBUG")) {
        for (int yy = 0; yy < TH; yy++) {
            for (int xx = 0; xx < TW; xx++)
                printf("%6.2f/%-3s", y[tidx(xx, yy)], dname(yd[tidx(xx, yy)]));
            printf("\n");
        }
        printf("indeg(3,3)=%d\n", y_indeg[tidx(3, 3)]);
    }
    check("confluence_indeg_two", y_indeg[tidx(3, 3)] >= 2);
    int y_acc[TN];
    for (int i = 0; i < TN; i++) y_acc[i] = 1;
    int yorder[TN];
    for (int i = 0; i < TN; i++) yorder[i] = i;
    for (int a = 0; a < TN; a++)
        for (int b = a + 1; b < TN; b++)
            if (y[yorder[b]] > y[yorder[a]]) { int t = yorder[a]; yorder[a] = yorder[b]; yorder[b] = t; }
    for (int i = 0; i < TN; i++) {
        int c = yorder[i];
        if (yd[c] >= 0) y_acc[tidx(c % TW + DX[yd[c]], c / TW + DY[yd[c]])] += y_acc[c];
    }
    int y_total = 0;
    for (int yy = 0; yy < TH; yy++)
        for (int xx = 0; xx < TW; xx++)
            if (yd[tidx(xx, yy)] < 0)
                y_total += y_acc[tidx(xx, yy)];   /* 所有边界出口累积之和 = 全部格点 */
    check("y_shape_drains_all_cells", y_total == TN);

    /* 5b. 实质汇流过滤(强入流 = 入流邻居汇水面积 >= 2, 即不止单点细沟):
     *    (3,3) 两条沟槽支流各自汇集多格, 强入流 >= 2;
     *    (3,4) 虽有侧向细沟汇入, 但细沟面积仅 1, 强入流 < 2。 */
    int strong_at_join = 0;
    for (int d = 0; d < 8; d++) {
        int nx = 3 + DX[d], ny = 3 + DY[d];
        if (nx < 0 || ny < 0 || nx >= TW || ny >= TH) continue;
        if (yd[tidx(nx, ny)] == (d + 4) % 8 && y_acc[tidx(nx, ny)] >= 2) strong_at_join++;
    }
    check("strong_filter_join", strong_at_join >= 2);
    int strong_at_pipe = 0;
    for (int d = 0; d < 8; d++) {
        int nx = 3 + DX[d], ny = 4 + DY[d];
        if (nx < 0 || ny < 0 || nx >= TW || ny >= TH) continue;
        if (yd[tidx(nx, ny)] == (d + 4) % 8 && y_acc[tidx(nx, ny)] >= 2) strong_at_pipe++;
    }
    check("strong_filter_pipe_excluded", strong_at_pipe < 2);

    /* 6. 三档真实地形: 零逆坡、零内部死角、有河有汇流, 且河道不是单一方向大直道 */
    int all_clean = 1, all_rivers = 1, shape_ok = 1;
    TierStat ss[3];
    for (int k = 0; k < 3; k++) {
        ss[k] = run_tier(&TIERS[k], 42u, NULL, 0, -1, -1, -1);
        if (ss[k].reverse_bad != 0 || ss[k].interior_sinks != 0) all_clean = 0;
        if (ss[k].river_cells == 0 || ss[k].nconf < 1) all_rivers = 0;
        if (ss[k].dirs_used < 3 || ss[k].max_run >= W) shape_ok = 0;
    }
    check("three_tiers_no_reverse", all_clean);
    check("three_tiers_have_rivers", all_rivers);
    check("rivers_not_single_straight_line", shape_ok);

    /* 7. 坡度档位单调: 平均坡度 缓 < 中 < 陡 */
    check("slope_tiers_ordered", ss[0].mean_s < ss[1].mean_s && ss[1].mean_s < ss[2].mean_s);

    /* 8. 确定性: 同种子两次结果一致 */
    TierStat a2 = run_tier(&TIERS[1], 7u, NULL, 0, -1, -1, -1);
    TierStat b2 = run_tier(&TIERS[1], 7u, NULL, 0, -1, -1, -1);
    check("deterministic", a2.river_cells == b2.river_cells && a2.nconf == b2.nconf &&
                           a2.filled == b2.filled && a2.max_run == b2.max_run);

    /* 9. 故障注入: 强制内部格点 (40,30) 向北(逆整体倾斜), 审计必须报逆坡并定位到该格点 */
    TierStat f = run_tier(&TIERS[1], 42u, NULL, 0, 40, 30, 6);
    int fault_located = 0;
    for (int i = 0; i < f.bad_total && i < 16; i++)
        if (f.bad_kind[i] == 1 && f.bad[i].from == idx(40, 30)) fault_located = 1;
    check("fault_injection_flagged_reverse", f.reverse_bad >= 1 && fault_located);

    /* 10. 每档真实地形都存在带"上坡邻居"的内部格点(自动注入演示的前提) */
    int auto_ok = 1;
    for (int k = 0; k < 3; k++) {
        Grid ag; terrain_build(&ag, &TIERS[k], 42u);
        priority_flood(&ag); flow_d8(&ag);
        int any = 0;
        for (int y = 1; y < H - 1 && !any; y++)
            for (int x = 1; x < W - 1 && !any; x++)
                if (pick_uphill_dir(&ag, x, y) >= 0) any = 1;
        if (!any) auto_ok = 0;
    }
    check("uphill_neighbor_exists_each_tier", auto_ok);

    printf("all_ok=%d\n", failures == 0);
}

/* ---------- 逆坡演示模式: 人为注入一条逆流, 验证"判未完成 + 指出位置" ---------- */

static int demo_fault(int tier_k, int fx, int fy, int force_dir) {
    const Tier* t = &TIERS[tier_k];
    Grid g;
    terrain_build(&g, t, 42u);
    int filled = priority_flood(&g);
    float mean_s, max_s;
    slope_stats(&g, &mean_s, &max_s);
    flow_d8(&g);
    if (!inside(fx, fy) || is_border(fx, fy)) {
        printf("fault inject failed: (%d,%d) 不是内部格点\n", fx, fy);
        return 3;
    }
    int orig_dir = g.dir[idx(fx, fy)];
    int use_dir = force_dir;
    int auto_relocated = 0;

    /* 自动模式或指定方向不构成上坡时, 取该格点最陡上坡邻居;
     * 若该格点是局部峰顶(无更高邻居), 在全图找一个有上坡邻居的内部格点。 */
    int chosen = pick_uphill_dir(&g, fx, fy);
    if (use_dir < 0) {
        if (chosen < 0) {
            for (int y = 1; y < H - 1 && chosen < 0; y++)
                for (int x = 1; x < W - 1 && chosen < 0; x++) {
                    int d = pick_uphill_dir(&g, x, y);
                    if (d >= 0) { fx = x; fy = y; chosen = d; auto_relocated = 1; }
                }
        }
        use_dir = chosen;
    } else {
        int nx = fx + DX[use_dir], ny = fy + DY[use_dir];
        if (!inside(nx, ny) || g.h[idx(nx, ny)] <= g.h[idx(fx, fy)]) {
            printf("注意: 指定方向 %s 不构成上坡, 改用 (%d,%d) 的最陡上坡方向 %s\n",
                   dname(use_dir), fx, fy, chosen >= 0 ? dname(chosen) : "无");
            if (chosen < 0) {
                for (int y = 1; y < H - 1 && chosen < 0; y++)
                    for (int x = 1; x < W - 1 && chosen < 0; x++) {
                        int d = pick_uphill_dir(&g, x, y);
                        if (d >= 0) { fx = x; fy = y; chosen = d; auto_relocated = 1; }
                    }
            }
            use_dir = chosen;
        }
    }
    if (use_dir < 0) {
        printf("fault inject failed: 全图无可用上坡邻居(退化地形)\n");
        return 3;
    }
    int rc = inject_fault(&g, fx, fy, use_dir);
    if (rc != 0) {
        printf("fault inject failed: cell/方向非法\n");
        return 3;
    }
    printf("fault demo on [%s] seed=42\n", t->key);
    if (auto_relocated) printf("原目标为局部峰顶, 已自动改取内部格点 (%d,%d)\n", fx, fy);
    printf("注入: 内部格点 (%d,%d) 原流向 %s -> 强制改为 %s\n",
           fx, fy, dname(orig_dir), dname(use_dir));
    BadEdge bad[16]; int kind[16];
    int total = reverse_audit(&g, bad, kind, 16);
    for (int i = 0; i < total && i < 16; i++) {
        int c = bad[i].from, ux = c % W, uy = c / W;
        if (kind[i] == 1) {
            int nc = bad[i].to;
            printf("逆坡位置 #%d: 格点 (%d,%d) 沿 %s 流向 (%d,%d), 上游高程 %.4f >= 下游高程 %.4f (高差 %+.4f)\n",
                   i + 1, ux, uy, dname(g.dir[c]), nc % W, nc / W,
                   g.h[c], g.h[nc], g.h[nc] - g.h[c]);
        } else {
            printf("死角位置 #%d: 内部格点 (%d,%d) 无下游\n", i + 1, ux, uy);
        }
    }
    if (total > 16) printf("...另有 %d 条问题未展开\n", total - 16);
    printf("审计结论: %s (问题边 %d 条, 填洼 %d 格)\n",
           total > 0 ? "未完成: 存在逆坡流动" : "完成: 无逆坡", total, filled);
    return total > 0 ? 2 : 0;
}

/* ---------- 出图模式 ---------- */

static int generate(const char* outdir) {
    char report[512];
    snprintf(report, sizeof(report), "%s/report.txt", outdir);
    FILE* rp = fopen(report, "w");
    if (!rp) { perror(report); return 3; }
    int fail = 0;
    printf("outdir=%s\n", outdir);
    fprintf(rp, "mini-proc-lab river 生成报告 (种子 42, 网格 %dx%d, 汇水阈值 %d)\n\n", W, H, RIVER_THRESH);
    for (int k = 0; k < 3; k++) {
        TierStat st = run_tier(&TIERS[k], 42u, outdir, 1, -1, -1, -1);
        printf("[%s] mean_slope=%.4f max_slope=%.4f filled=%d reverse=%d sinks=%d rivers=%d confluences=%d(on-network %d) dirs=%d max_run=%d\n",
               TIERS[k].key, st.mean_s, st.max_s, st.filled, st.reverse_bad, st.interior_sinks,
               st.river_cells, st.nconf, st.nconf_draw, st.dirs_used, st.max_run);
        fprintf(rp, "== %s ==\n", TIERS[k].name);
        fprintf(rp, "平均坡度 %.4f, 最大坡度 %.4f, 填洼格点 %d\n", st.mean_s, st.max_s, st.filled);
        fprintf(rp, "逆坡边 %d, 内部死角 %d -> %s\n", st.reverse_bad, st.interior_sinks,
                (st.reverse_bad == 0 && st.interior_sinks == 0) ? "通过(无逆坡流动)" : "未完成");
        fprintf(rp, "河道格点 %d, 入度>=2 格点 %d, 其中实质汇流点(两条入流均达河道阈值) %d, 河道使用流向数 %d, 最长同向直段 %d 格\n\n",
                st.river_cells, st.nconf, st.nconf_draw, st.dirs_used, st.max_run);
        if (st.reverse_bad != 0 || st.interior_sinks != 0) fail = 1;
    }
    fprintf(rp, "图例: river_<tier>.svg 为河道图, river_<tier>.txt 为 ASCII 图(#=河道 o=汇流点 ~=出口)。\n");
    fclose(rp);
    if (fail) {
        printf("RESULT=未完成: 存在逆坡流动或内部死角, 详见 %s\n", report);
        return 2;
    }
    printf("RESULT=完成: 三档均无逆坡流动, 报告 %s\n", report);
    return 0;
}

int main(int argc, char** argv) {
    if (argc >= 2 && strcmp(argv[1], "generate") == 0) {
        const char* outdir = argc >= 3 ? argv[2] : "river_out";
        return generate(outdir);
    }
    if (argc >= 2 && strcmp(argv[1], "demo-fault") == 0) {
        /* 用法: demo-fault [tier 0..2] [x] [y] [dir 0..7, -1=自动取最陡上坡] */
        int tier = argc >= 3 ? atoi(argv[2]) : 1;
        int fx = argc >= 4 ? atoi(argv[3]) : 40;
        int fy = argc >= 5 ? atoi(argv[4]) : 30;
        int fd = argc >= 6 ? atoi(argv[5]) : -1; /* -1: 自动选最陡上坡邻居 */
        if (tier < 0 || tier > 2) tier = 1;
        return demo_fault(tier, fx, fy, fd);
    }
    tests();
    return failures == 0 ? 0 : 1;
}
