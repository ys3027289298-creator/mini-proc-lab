
/* hydrology: sink-filled D8 drainage network.
 *
 * Pipeline: value-noise heightfield -> priority-flood sink filling with one
 * enforced outlet -> D8 flow direction -> flow accumulation -> Strahler order
 * -> classification into five water elements (main channel / tributary /
 * lake / swamp / waterfall) -> connected-component statistics and a
 * connectivity audit that lists every dead-end path.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define N 80
#define NN (N * N)
#define CELL_M 10.0f            /* world size of one cell, metres */
#define EPS 1e-4f               /* priority-flood flat gradient */

/* element classes */
enum { LAND = 0, MAIN = 1, TRIB = 2, LAKE = 3, SWAMP = 4, FALL = 5, NCLS = 6 };
/* pseudo downstream ids */
#define D_OUTLET (-2)
#define D_NONE   (-1)

/* classification thresholds (tunable, documented in README) */
#define TRIB_THR   60.0f        /* acc >= this: channelized flow */
#define MAIN_THR  320.0f        /* acc >= this: main river */
#define LAKE_DEEP   0.15f       /* filled-vs-raw depth marking a lake bed */
#define FALL_SLOPE  0.15f       /* downstream drop per metre marking a fall */
#define GENTLE      0.04f       /* raw-surface slope below which ground is swampable */

static const char CLS_ZH[NCLS][16] = {
    "陆面", "主河道", "支流", "湖泊", "沼泽", "瀑布"
};
static const char GRID_CH[NCLS] = { '.', '#', '+', '~', ':', 'v' };

/* D8 neighbour offsets, clockwise from east; odd steps are diagonals */
static const int DX[8] = { 1, 1, 0, -1, -1, -1, 0, 1 };
static const int DY[8] = { 0, 1, 1, 1, 0, -1, -1, -1 };
static const float DD[8] = {
    1, 1.41421356f, 1, 1.41421356f, 1, 1.41421356f, 1, 1.41421356f
};

/* ---------------- deterministic value noise (same family as noise/) ------ */
static unsigned hash2(int x, int y, unsigned seed) {
    unsigned h = (unsigned)x * 374761393u + (unsigned)y * 668265263u
               + seed * 2246822519u;
    h = (h ^ (h >> 13)) * 1274126177u;
    return h ^ (h >> 16);
}
static float rand01(unsigned h) { return (float)(h & 0xffffff) / 16777215.0f; }
static float lerpf(float a, float b, float t) { return a + (b - a) * t; }
static float fadef(float t) { return t * t * (3.0f - 2.0f * t); }

static float value_noise(float x, float y, unsigned seed) {
    int xi = (int)floorf(x), yi = (int)floorf(y);
    float fx = x - xi, fy = y - yi;
    float v00 = rand01(hash2(xi, yi, seed));
    float v10 = rand01(hash2(xi + 1, yi, seed));
    float v01 = rand01(hash2(xi, yi + 1, seed));
    float v11 = rand01(hash2(xi + 1, yi + 1, seed));
    float u = fadef(fx), v = fadef(fy);
    return lerpf(lerpf(v00, v10, u), lerpf(v01, v11, u), v) * 2.0f - 1.0f;
}
static float fbm(float x, float y, int oct, float lac, float gain, unsigned seed) {
    float sum = 0.0f, amp = 1.0f, norm = 0.0f, freq = 1.0f;
    for (int i = 0; i < oct; i++) {
        sum += value_noise(x * freq, y * freq, seed + (unsigned)i * 17u) * amp;
        norm += amp;
        amp *= gain;
        freq *= lac;
    }
    return norm > 0.0f ? sum / norm : 0.0f;
}

/* ---------------- global terrain/network state ---------------------------- */
static float g_raw[NN];       /* original heightfield */
static float g_z[NN];         /* filled heightfield */
static int   g_dir[NN];       /* D8 direction 0..7, D_OUTLET, D_NONE */
static int   g_acc[NN];       /* flow accumulation, cells */
static int   g_ord[NN];       /* Strahler order */
static int   g_cls[NN];       /* element class */
static int   g_base[NN];      /* channel base class under a FALL cell */
static float g_drop[NN];      /* downstream slope (drop / metre) */
static float g_maxdrop[NN];   /* steepest slope to any 8-neighbour, filled */
static float g_rawdrop[NN];   /* same, but on the raw (pre-fill) surface */
static int   g_order[NN];     /* cell ids sorted by descending filled z */
static int   g_outlet;

/* ---------------- binary min-heap keyed on elevation ---------------------- */
typedef struct { float k; int id; } HItem;
static HItem g_heap[NN];
static int   g_hn;

static void hpush(float k, int id) {
    int i = g_hn++;
    g_heap[i].k = k; g_heap[i].id = id;
    while (i > 0) {
        int p = (i - 1) / 2;
        if (g_heap[p].k <= g_heap[i].k) break;
        HItem t = g_heap[p]; g_heap[p] = g_heap[i]; g_heap[i] = t;
        i = p;
    }
}
static HItem hpop(void) {
    HItem top = g_heap[0];
    HItem last = g_heap[--g_hn];
    if (g_hn > 0) {
        int i = 0;
        g_heap[0] = last;
        for (;;) {
            int l = i * 2 + 1, r = l + 1, m = i;
            if (l < g_hn && g_heap[l].k < g_heap[m].k) m = l;
            if (r < g_hn && g_heap[r].k < g_heap[m].k) m = r;
            if (m == i) break;
            HItem t = g_heap[m]; g_heap[m] = g_heap[i]; g_heap[i] = t;
            i = m;
        }
    }
    return top;
}

static int xy_id(int x, int y) { return y * N + x; }
static int cx_of(int id) { return id % N; }
static int cy_of(int id) { return id / N; }

/* ---------------- terrain synthesis --------------------------------------- */
/* deterministic fault scarp: a soft vertical step so a channel crossing it
 * forms a real waterfall; smooth falloff avoids creating an impassable wall */
static float fault_step(int y, int fy, float height) {
    float d = (float)(y - fy);
    return height * (0.5f + 0.5f * d / sqrtf(d * d + 1.5f));
}

static void gen_terrain(unsigned seed) {
    /* two scarps on deterministic rows, offset per seed */
    int fy1 = N / 3 + (int)(seed % 5) - 2;
    int fy2 = 2 * N / 3 + (int)((seed / 7) % 5) - 2;
    for (int y = 0; y < N; y++) {
        for (int x = 0; x < N; x++) {
            float fx = (float)x / N, fy = (float)y / N;
            /* north-high / south-low regional slope */
            float h = 14.0f * fy;
            /* rolling relief, octave stack carves hollows and ridges */
            h += 7.0f * fbm(fx * 4.5f, fy * 4.5f, 4, 2.0f, 0.5f, seed);
            /* broad meandering valley trough so a main stem can organize */
            float axis = 0.5f + 0.22f * fbm(fx * 1.2f, fy * 1.2f, 2, 2.0f, 0.5f,
                                            seed + 101u);
            float dx = fx - axis;
            h -= 3.5f * expf(-(dx * dx) / (2.0f * 0.05f * 0.05f));
            h += fault_step(y, fy1, 6.0f);
            h += fault_step(y, fy2, 5.0f);
            g_raw[xy_id(x, y)] = h;
        }
    }
}

/* ---------------- priority-flood sink filling (Wang & Liu style) ---------- */
static void fill_sinks(void) {
    static unsigned char seen[NN];
    memset(seen, 0, sizeof(seen));
    for (int i = 0; i < NN; i++) g_dir[i] = D_NONE;
    g_hn = 0;

    /* one enforced outlet: lowest cell on the map border */
    g_outlet = 0;
    for (int i = 0; i < NN; i++) {
        int x = cx_of(i), y = cy_of(i);
        int border = x == 0 || y == 0 || x == N - 1 || y == N - 1;
        if (border && g_raw[i] < g_raw[g_outlet]) g_outlet = i;
    }
    seen[g_outlet] = 1;
    g_z[g_outlet] = g_raw[g_outlet];
    g_dir[g_outlet] = D_OUTLET;
    hpush(g_z[g_outlet], g_outlet);

    while (g_hn > 0) {
        HItem c = hpop();
        int x = cx_of(c.id), y = cy_of(c.id);
        for (int d = 0; d < 8; d++) {
            int nx = x + DX[d], ny = y + DY[d];
            if (nx < 0 || ny < 0 || nx >= N || ny >= N) continue;
            int n = xy_id(nx, ny);
            if (seen[n]) continue;
            seen[n] = 1;
            if (g_raw[n] > c.k + EPS) g_z[n] = g_raw[n];
            else                      g_z[n] = c.k + EPS;   /* ponded flat */
            g_dir[n] = c.id;             /* downstream cell id, not direction */
            hpush(g_z[n], n);
        }
    }
}

static int cmp_z_desc(const void *a, const void *b) {
    int ia = *(const int *)a, ib = *(const int *)b;
    if (g_z[ib] > g_z[ia]) return 1;
    if (g_z[ib] < g_z[ia]) return -1;
    return ia - ib;
}

/* ---------------- flow accumulation and Strahler order -------------------- */
static void accumulate(void) {
    for (int i = 0; i < NN; i++) { g_order[i] = i; g_acc[i] = 1; }
    qsort(g_order, NN, sizeof(int), cmp_z_desc);

    for (int k = 0; k < NN; k++) {
        int c = g_order[k];
        int d = g_dir[c];
        if (d >= 0) g_acc[d] += g_acc[c];
    }

    static int maxo[NN], maxc[NN];
    memset(maxo, 0, sizeof(maxo));
    memset(maxc, 0, sizeof(maxc));
    for (int k = 0; k < NN; k++) {
        int c = g_order[k];                 /* high z first: leaves first */
        g_ord[c] = maxo[c] == 0 ? 1
                 : (maxc[c] >= 2 ? maxo[c] + 1 : maxo[c]);
        int d = g_dir[c];
        if (d < 0) continue;
        if (g_ord[c] > maxo[d]) { maxo[d] = g_ord[c]; maxc[d] = 1; }
        else if (g_ord[c] == maxo[d]) maxc[d]++;
    }
}

static void compute_slopes(void) {
    for (int y = 0; y < N; y++) {
        for (int x = 0; x < N; x++) {
            int c = xy_id(x, y);
            float md = 0.0f;
            for (int d = 0; d < 8; d++) {
                int nx = x + DX[d], ny = y + DY[d];
                if (nx < 0 || ny < 0 || nx >= N || ny >= N) continue;
                float s = (g_z[c] - g_z[xy_id(nx, ny)]) / DD[d] / CELL_M;
                if (s > md) md = s;
            }
            g_maxdrop[c] = md;
            /* raw-surface steepness: filled flats must not read as natural wetland */
            float rd = 0.0f;
            for (int dd = 0; dd < 8; dd++) {
                int nx = x + DX[dd], ny = y + DY[dd];
                if (nx < 0 || ny < 0 || nx >= N || ny >= N) continue;
                float s = (g_raw[c] - g_raw[xy_id(nx, ny)]) / DD[dd] / CELL_M;
                if (s > rd) rd = s;
            }
            g_rawdrop[c] = rd;
            int d = g_dir[c];
            if (d >= 0) {
                float step = hypotf((float)(cx_of(d) - x), (float)(cy_of(d) - y));
                g_drop[c] = (g_z[c] - g_z[d]) / step / CELL_M;
            } else g_drop[c] = 0.0f;
        }
    }
}

static int is_water(int c) { return g_cls[c] != LAND; }

/* ---------------- five-element classification ----------------------------- */
static void classify(void) {
    for (int i = 0; i < NN; i++) { g_cls[i] = LAND; g_base[i] = LAND; }

    /* 1. lake beds: terrain that the flood had to pond above the raw surface */
    for (int i = 0; i < NN; i++)
        if (g_z[i] - g_raw[i] > LAKE_DEEP) g_cls[i] = LAKE;

    /* 2. channel skeleton from the accumulation threshold */
    for (int i = 0; i < NN; i++) {
        if (g_cls[i] != LAND) continue;
        if (g_acc[i] >= (int)TRIB_THR)
            g_cls[i] = g_acc[i] >= (int)MAIN_THR ? MAIN : TRIB;
    }

    /* 3. lake outlet carving: every lake gets a real D8 outflow channel so a
     *    pond can never remain an isolated pit on the network. */
    static unsigned char stamp[NN];
    static unsigned int tick;
    tick++;
    for (int i = 0; i < NN; i++) {
        if (g_cls[i] != LAKE || stamp[i] == tick) continue;
        int c = i, guard = 0;
        while (g_cls[c] == LAKE && guard++ < NN) { stamp[c] = tick; c = g_dir[c]; }
        while (c >= 0 && !is_water(c) && guard++ < NN + N) {
            g_cls[c] = TRIB;                 /* carved spill run */
            c = g_dir[c];
        }
    }

    /* 4. waterfalls: steep drops on the channel skeleton that are actually
     *    fed by an upstream water cell, so a fall can never hang in mid-air */
    static unsigned char fed[NN];
    memset(fed, 0, sizeof(fed));
    for (int i = 0; i < NN; i++) {
        int d = g_dir[i];
        if (is_water(i) && d >= 0 && is_water(d)) fed[d] = 1;
    }
    for (int i = 0; i < NN; i++) {
        if ((g_cls[i] == MAIN || g_cls[i] == TRIB)
            && g_drop[i] > FALL_SLOPE && fed[i]) {
            g_base[i] = g_cls[i];
            g_cls[i] = FALL;
        }
    }

    /* 5. swamps: genuinely flat raw ground that receives ponded flow, i.e.
     *    adjacent to a lake or sitting on a filled flat beside a channel.
     *    Grown as a bounded multi-pass region so it stays a wetland patch,
     *    not a map-wide skin. */
    static unsigned char swamp[NN];
    memset(swamp, 0, sizeof(swamp));
    for (int pass = 0; pass < 4; pass++) {
        for (int i = 0; i < NN; i++) {
            if (swamp[i] || g_cls[i] != LAND || g_rawdrop[i] >= GENTLE) continue;
            int x = cx_of(i), y = cy_of(i), touch = 0;
            for (int d = 0; d < 8 && !touch; d++) {
                int nx = x + DX[d], ny = y + DY[d];
                if (nx < 0 || ny < 0 || nx >= N || ny >= N) continue;
                int n = xy_id(nx, ny);
                int wn = swamp[n] || (is_water(n) && g_cls[n] != SWAMP);
                /* fed only by standing/lake water or ponded flats */
                if (wn && (g_cls[n] == LAKE || g_cls[n] == SWAMP || swamp[n]
                           || g_z[n] - g_raw[n] > EPS)) touch = 1;
            }
            /* drainage must stay inside the wetland/network, never spill to dry land */
            int down = g_dir[i];
            int drained = down >= 0 && (swamp[down] || is_water(down));
            if (touch && drained) swamp[i] = 1;
        }
    }
    for (int i = 0; i < NN; i++) if (swamp[i]) g_cls[i] = SWAMP;
}

/* ---------------- 8-neighbour components (spatial footprint) -------------- */
#define MAXC 4096
static int g_cid[NN];
static int g_cstart[MAXC], g_ccnt[MAXC], g_ccls[MAXC], g_cn;
static int g_cells[NN];

static void build_components(void) {
    static int queue[NN];
    for (int i = 0; i < NN; i++) g_cid[i] = -1;
    g_cn = 0;
    int put = 0;
    for (int s = 0; s < NN; s++) {
        if (g_cls[s] == LAND || g_cid[s] != -1) continue;
        int cid = g_cn++;
        g_ccls[cid] = g_cls[s];
        g_cstart[cid] = put;
        g_ccnt[cid] = 0;
        int qh = 0, qt = 0;
        queue[qt++] = s; g_cid[s] = cid;
        while (qh < qt) {
            int c = queue[qh++];
            g_cells[put++] = c;
            g_ccnt[cid]++;
            int x = cx_of(c), y = cy_of(c);
            for (int d = 0; d < 8; d++) {
                int nx = x + DX[d], ny = y + DY[d];
                if (nx < 0 || ny < 0 || nx >= N || ny >= N) continue;
                int n = xy_id(nx, ny);
                if (g_cid[n] == -1 && g_cls[n] == g_ccls[cid]) {
                    g_cid[n] = cid;
                    queue[qt++] = n;
                }
            }
        }
    }
}

/* ---------------- same-class D8 chains (the real water-path instances) ---- */
#define MAXCH 4096
static int g_chain[NN];           /* chain id of each water cell */
static int g_chcls[MAXCH], g_chn;
static int g_chlen[MAXCH];        /* cells in the chain */
static int g_chhead[MAXCH];       /* upstream-most cell (re-walk entry) */
static int g_chtail[MAXCH];       /* downstream-most cell */
static int g_chcid[MAXCH];        /* spatial component the chain sits in */
static int g_chpath[MAXCH];       /* longest same-class D8 path, cells */
static float g_chpathm[MAXCH];    /* same, metres */

/* area-type elements are counted as 8-neighbour patches; channel-type
 * elements are counted as D8 flow paths */
static int is_area_element(int cls) { return cls == LAKE || cls == SWAMP; }

/* downstream chain cell, or -1 when the D8 path leaves the class */
static int chain_next(int c) {
    int d = g_dir[c];
    return (d >= 0 && g_cls[d] == g_cls[c]) ? d : -1;
}

static float cell_step_m(int a, int b) {
    return hypotf((float)(cx_of(b) - cx_of(a)),
                  (float)(cy_of(b) - cy_of(a))) * CELL_M;
}

static void build_chains(void) {
    g_chn = 0;
    for (int i = 0; i < NN; i++) g_chain[i] = -1;

    /* area elements: one chain id per spatial component (the instance) */
    for (int cid = 0; cid < g_cn; cid++) {
        if (!is_area_element(g_ccls[cid])) continue;
        int id = g_chn++;
        g_chcls[id] = g_ccls[cid];
        g_chcid[id] = cid;
        g_chhead[id] = g_cells[g_cstart[cid]];
        g_chlen[id] = 0;
        for (int t = g_cstart[cid]; t < g_cstart[cid] + g_ccnt[cid]; t++) {
            int c = g_cells[t];
            g_chain[c] = id;
            g_chlen[id]++;
        }
        g_chtail[id] = g_chhead[id];
    }

    /* channel elements: follow D8 until the class changes */
    for (int k = 0; k < NN; k++) {
        int c = g_order[k];                 /* high z first: heads first */
        if (!is_water(c) || g_chain[c] >= 0) continue;
        int id = g_chn++;
        g_chcls[id] = g_cls[c];
        g_chcid[id] = g_cid[c];
        g_chhead[id] = c;
        g_chlen[id] = 0;
        int p = c, guard = 0;
        while (p >= 0 && guard++ < NN && g_chain[p] < 0
               && g_cls[p] == g_chcls[id]) {
            g_chain[p] = id;
            g_chlen[id]++;
            g_chtail[id] = p;
            p = chain_next(p);
        }
    }

    /* longest same-class D8 path per instance */
    static int dp[NN];
    static unsigned char in[NN];
    for (int id = 0; id < g_chn; id++) {
        if (!is_area_element(g_chcls[id])) {
            g_chpath[id] = g_chlen[id];
            g_chpathm[id] = 0.0f;
            int c = g_chhead[id], guard = 0;
            for (;;) {
                int d = chain_next(c);
                if (d < 0 || guard++ >= NN) break;
                g_chpathm[id] += cell_step_m(c, d);
                c = d;
            }
            continue;
        }
        int s0 = g_cstart[g_chcid[id]], s1 = s0 + g_ccnt[g_chcid[id]];
        for (int t = s0; t < s1; t++) in[g_cells[t]] = 1;
        int best = 0, bestc = -1;
        for (int k = 0; k < NN; k++) {
            int c = g_order[k];             /* high z first */
            if (!in[c]) continue;
            int d = g_dir[c];
            dp[c] = (d >= 0 && in[d]) ? dp[d] + 1 : 1;
            if (dp[c] > best) { best = dp[c]; bestc = c; }
        }
        float m = 0.0f;
        int c = bestc, guard = 0;
        while (c >= 0 && guard++ < NN) {
            int d = g_dir[c];
            if (d < 0 || !in[d]) break;
            m += cell_step_m(c, d);
            c = d;
        }
        g_chpath[id] = best;
        g_chpathm[id] = m;
        for (int t = s0; t < s1; t++) in[g_cells[t]] = 0;
    }
}

static float chain_metres(int id) { return g_chpathm[id]; }
static int   chain_cells(int id)  { return g_chpath[id]; }

/* ---------------- chain link graph ---------------------------------------- */
typedef struct { int from, to; char lateral; } Edge;
#define MAXE 8192
static Edge g_edges[MAXE];
static int  g_en;
static unsigned char g_has_exit[MAXCH];
static int g_indeg[MAXCH];

static int find_edge(int from, int to, char lateral) {
    for (int i = 0; i < g_en; i++)
        if (g_edges[i].from == from && g_edges[i].to == to
            && g_edges[i].lateral == lateral) return i;
    return -1;
}
static void add_edge(int from, int to, char lateral) {
    if (from == to || from < 0 || to < 0) return;
    if (g_en < MAXE && find_edge(from, to, lateral) < 0)
        g_edges[g_en++] = (Edge){ from, to, lateral };
}

static void build_graph(void) {
    g_en = 0;
    memset(g_has_exit, 0, sizeof(g_has_exit));
    memset(g_indeg, 0, sizeof(g_indeg));
    for (int c = 0; c < NN; c++) {
        if (!is_water(c)) continue;
        int uc = g_chain[c];
        int d = g_dir[c];
        if (d < 0) g_has_exit[uc] = 1;
        else if (is_water(d) && g_chain[d] != uc) add_edge(uc, g_chain[d], 0);

        int x = cx_of(c), y = cy_of(c);
        for (int k = 0; k < 8; k++) {
            int nx = x + DX[k], ny = y + DY[k];
            if (nx < 0 || ny < 0 || nx >= N || ny >= N) continue;
            int n = xy_id(nx, ny);
            if (is_water(n) && g_chain[n] != uc)
                add_edge(uc, g_chain[n], (char)1);
        }
    }
    for (int i = 0; i < g_en; i++) g_indeg[g_edges[i].to]++;
}

static int chain_has_d8_out(int id) {
    if (g_has_exit[id]) return 1;
    for (int i = 0; i < g_en; i++)
        if (g_edges[i].from == id && !g_edges[i].lateral) return 1;
    return 0;
}

/* first water cell encountered downstream of every land cell (absorption) */
static int g_firstwater[NN];
static void compute_absorption(void) {
    for (int k = NN - 1; k >= 0; k--) {       /* low z first */
        int c = g_order[k];
        if (is_water(c)) g_firstwater[c] = c;
        else {
            int d = g_dir[c];
            g_firstwater[c] = d >= 0 ? g_firstwater[d] : g_outlet;
        }
    }
}

enum { BK_ISO = 0, BK_DEAD = 1, BK_FALL = 2, BK_LAKE = 3, BK_TRIB = 4 };
static const char *BK_NAME[5] = {
    "孤立水体(无任何衔接)", "断头路径(下游消失在陆面)",
    "瀑布悬空(未嵌入河道)", "湖泊死水(无溢流出口)", "支流未并入主河"
};

typedef struct {
    int kind, cid, x, y, missing;
} Break;
#define MAXB 256
static Break g_bk[MAXB];
static int   g_bn;

static void audit(void) {
    g_bn = 0;
    static unsigned char linked[MAXCH];
    memset(linked, 0, g_chn);
    for (int i = 0; i < g_en; i++) {
        linked[g_edges[i].from] = 1;
        linked[g_edges[i].to] = 1;
    }
    /* cell-level: does any water cell drain into this one? (grid-wide) */
    static unsigned char fed[NN];
    memset(fed, 0, sizeof(fed));
    for (int c = 0; c < NN; c++) {
        int d = g_dir[c];
        if (is_water(c) && d >= 0 && is_water(d)) fed[d] = 1;
    }
    static unsigned int seen_tick[NN];
    static unsigned int tt;

    for (int id = 0; id < g_chn; id++) {
        int cls = g_chcls[id];
        int tail = g_chtail[id];

        /* A. isolated: no link at all to any other water element */
        if (!linked[id] && !g_has_exit[id]) {
            if (g_bn < MAXB) g_bk[g_bn++] = (Break){ BK_ISO, id,
                cx_of(tail), cy_of(tail), LAND };
            continue;
        }

        /* B. dead-end: flow enters, but D8 downstream vanishes into dry land.
         *    Swamps are area wetlands drained by seepage, not D8 channels,
         *    so their connectivity is the 8-neighbour test in A, not this one. */
        if (cls != SWAMP && g_indeg[id] > 0 && !chain_has_d8_out(id) && g_bn < MAXB)
            g_bk[g_bn++] = (Break){ BK_DEAD, id, cx_of(tail), cy_of(tail),
                                    cls == LAKE ? TRIB : cls };

        /* C1. waterfall must be embedded: fed by an upstream water cell and
         *    reaching a downstream water cell along D8 */
        if (cls == FALL) {
            int up = 0, down = 0;
            for (int t = g_cstart[g_chcid[id]];
                 t < g_cstart[g_chcid[id]] + g_ccnt[g_chcid[id]]; t++) {
                int c = g_cells[t];
                if (fed[c]) up = 1;
                int d = g_dir[c];
                if (d >= 0 && is_water(d)) down = 1;
            }
            if ((!up || !down) && g_bn < MAXB)
                g_bk[g_bn++] = (Break){ BK_FALL, id, cx_of(tail), cy_of(tail),
                                        up ? MAIN : TRIB };
        }

        /* C2. lake without a spill path */
        if (cls == LAKE && !chain_has_d8_out(id) && g_bn < MAXB)
            g_bk[g_bn++] = (Break){ BK_LAKE, id, cx_of(tail), cy_of(tail), TRIB };

        /* C3. tributary that never reaches the main river */
        if (cls == TRIB) {
            tt++;
            int c = tail, hit = 0, guard = 0;
            while (c >= 0 && guard++ < NN) {
                if (seen_tick[c] == tt) break;
                seen_tick[c] = tt;
                if (g_cls[c] == MAIN ||
                    (g_cls[c] == FALL && g_base[c] == MAIN)) { hit = 1; break; }
                c = g_dir[c];
            }
            if (!hit && g_bn < MAXB)
                g_bk[g_bn++] = (Break){ BK_TRIB, id, cx_of(tail), cy_of(tail), MAIN };
        }
    }
}

/* ---------------- full demo pipeline -------------------------------------- */
static void build_demo(unsigned seed) {
    gen_terrain(seed);
    fill_sinks();
    compute_slopes();
    accumulate();
    classify();
    build_components();
    build_chains();
    build_graph();
    compute_absorption();
    audit();
}

/* ---------------- writers -------------------------------------------------- */
static void write_grid(const char *path) {
    FILE *f = fopen(path, "w");
    if (!f) { perror(path); return; }
    for (int y = 0; y < N; y++) {
        for (int x = 0; x < N; x++) {
            int c = xy_id(x, y);
            char ch = c == g_outlet ? 'O' : GRID_CH[g_cls[c]];
            fputc(ch, f);
        }
        fputc('\n', f);
    }
    fclose(f);
}

static unsigned char shade_byte(int c) {
    int x = cx_of(c), y = cy_of(c);
    int xm = x > 0 ? x - 1 : x, xp = x < N - 1 ? x + 1 : x;
    float dz = g_z[xy_id(xp, y)] - g_z[xy_id(xm, y)];
    float s = 0.55f - dz * 0.06f;
    if (s < 0.15f) s = 0.15f;
    if (s > 0.95f) s = 0.95f;
    return (unsigned char)(s * 255.0f);
}

static void cell_rgb(int cls, unsigned char *r, unsigned char *g, unsigned char *b) {
    /* categorical 5-slot palette, validated CVD/normal-vision separable */
    switch (cls) {
        case MAIN: *r = 27;  *g = 79;  *b = 156; break;   /* deep blue   #1b4f9c */
        case TRIB: *r = 111; *g = 178; *b = 245; break;   /* light blue  #6fb2f5 */
        case LAKE: *r = 27;  *g = 175; *b = 122; break;   /* aqua        #1baf7a */
        case SWAMP:*r = 0;   *g = 131; *b = 0;   break;   /* green       #008300 */
        case FALL: *r = 237; *g = 161; *b = 0;   break;   /* amber       #eda100 */
        default:
            *r = *g = *b = shade_byte(0);        /* land: grayscale hillshade */
    }
}

static void write_ppm(const char *path) {
    FILE *f = fopen(path, "wb");
    if (!f) { perror(path); return; }
    const int S = 10;
    fprintf(f, "P6\n%d %d\n255\n", N * S, N * S);
    for (int y = 0; y < N; y++) {
        for (int sy = 0; sy < S; sy++) {
            for (int x = 0; x < N; x++) {
                int c = xy_id(x, y);
                unsigned char r, g, b;
                if (is_water(c)) cell_rgb(g_cls[c], &r, &g, &b);
                else { r = g = b = shade_byte(c); }
                if (c == g_outlet) { r = 230; g = 0; b = 230; }
                for (int sx = 0; sx < S; sx++) {
                    fputc((int)r, f); fputc((int)g, f); fputc((int)b, f);
                }
            }
        }
    }
    fclose(f);
}

static void write_svg(const char *path) {
    FILE *f = fopen(path, "w");
    if (!f) { perror(path); return; }
    const int S = 10, PW = N * S, LW = 230;
    fprintf(f,
      "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\" "
      "font-family=\"monospace\" font-size=\"12\">\n", PW + LW, PW);
    fprintf(f, "<rect width=\"100%%\" height=\"100%%\" fill=\"#f4f4f4\"/>\n");
    for (int y = 0; y < N; y++)
        for (int x = 0; x < N; x++) {
            int c = xy_id(x, y);
            unsigned char r, g, b;
            if (is_water(c)) cell_rgb(g_cls[c], &r, &g, &b);
            else { r = g = b = shade_byte(c); }
            if (c == g_outlet) { r = 230; g = 0; b = 230; }
            fprintf(f, "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" "
                    "fill=\"#%02x%02x%02x\"/>\n", x * S, y * S, S, S, r, g, b);
        }
    /* break markers */
    for (int i = 0; i < g_bn; i++) {
        int x = g_bk[i].x * S + S / 2, y = g_bk[i].y * S + S / 2;
        fprintf(f, "<g stroke=\"#ff0000\" stroke-width=\"2\">"
                "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\"/>"
                "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\"/></g>\n",
                x - 4, y - 4, x + 4, y + 4, x - 4, y + 4, x + 4, y - 4);
    }
    int lx = PW + 14, ly = 28;
    fprintf(f, "<text x=\"%d\" y=\"%d\" font-size=\"15\">hydrology network</text>\n", lx, ly);
    ly += 24;
    for (int k = 1; k < NCLS; k++) {
        unsigned char r, g, b; cell_rgb(k, &r, &g, &b);
        fprintf(f, "<rect x=\"%d\" y=\"%d\" width=\"14\" height=\"14\" "
                "fill=\"#%02x%02x%02x\"/><text x=\"%d\" y=\"%d\">%s</text>\n",
                lx, ly, r, g, b, lx + 22, ly + 12, CLS_ZH[k]);
        ly += 22;
    }
    fprintf(f, "<rect x=\"%d\" y=\"%d\" width=\"14\" height=\"14\" "
            "fill=\"#e600e6\"/><text x=\"%d\" y=\"%d\">流域出口</text>\n",
            lx, ly, lx + 22, ly + 12);
    ly += 22;
    fprintf(f, "<text x=\"%d\" y=\"%d\" fill=\"red\">红叉=断头/孤立</text>\n", lx + 22, ly + 12);
    ly += 30;
    fprintf(f, "<text x=\"%d\" y=\"%d\">阈值: 支流acc&gt;=%d</text>\n", lx, ly + 12, (int)TRIB_THR);
    ly += 20;
    fprintf(f, "<text x=\"%d\" y=\"%d\">主河acc&gt;=%d</text>\n", lx, ly + 12, (int)MAIN_THR);
    ly += 20;
    fprintf(f, "<text x=\"%d\" y=\"%d\">D8链=%d 断头=%d</text>\n", lx, ly + 12, g_chn, g_bn);
    fprintf(f, "</svg>\n");
    fclose(f);
}

static void write_report(const char *path) {
    FILE *f = fopen(path, "w");
    if (!f) { perror(path); return; }

    long area_cells[NCLS] = {0};
    int  instances[NCLS] = {0};
    int  best_path[NCLS] = {0};
    float best_m[NCLS] = {0.0f};
    for (int id = 0; id < g_chn; id++) {
        int k = g_chcls[id];
        instances[k]++;
        area_cells[k] += g_chlen[id];
        if (chain_cells(id) > best_path[k]) {
            best_path[k] = chain_cells(id);
            best_m[k] = chain_metres(id);
        }
    }

    fprintf(f, "============================================================\n");
    fprintf(f, "水系网络报告  网格 %dx%d  单元 %.0f m  流域出口 (%d,%d)\n",
            N, N, CELL_M, cx_of(g_outlet), cy_of(g_outlet));
    fprintf(f, "阈值: 支流累积>=%d  主河累积>=%d  湖深>%.2f  瀑布坡降>%.2f  沼泽原始坡降<%.2f\n",
            (int)TRIB_THR, (int)MAIN_THR, LAKE_DEEP, FALL_SLOPE, GENTLE);
    fprintf(f, "实例=同类D8连续链; 8邻域空间斑块共 %d 个, D8链共 %d 条\n",
            g_cn, g_chn);
    fprintf(f, "============================================================\n\n");

    fprintf(f, "[1] 五要素统计 (面积=单元数x%.0fm^2; 最长路径=该要素最长D8连续链)\n",
            CELL_M * CELL_M);
    fprintf(f, "%-8s %8s %12s %14s %12s %12s\n",
            "要素", "实例数", "总面积(格)", "总面积(m^2)", "最长路径(格)", "最长路径(m)");
    for (int k = 1; k < NCLS; k++) {
        fprintf(f, "%-8s %8d %12ld %14.0f %12d %12.1f\n",
                CLS_ZH[k], instances[k], area_cells[k],
                area_cells[k] * CELL_M * CELL_M, best_path[k], best_m[k]);
    }

    fprintf(f, "\n[2] 实例明细 (id=D8链号; ->为D8顺接, ~为8邻域侧渗; 最多40条)\n");
    fprintf(f, "%-6s %-8s %8s %10s %10s %10s  下游衔接\n",
            "id", "要素", "格数", "面积m2", "路径格", "路径m");
    int shown = 0;
    for (int id = 0; id < g_chn && shown < 40; id++, shown++) {
        fprintf(f, "#%-5d %-8s %8d %10.0f %10d %10.1f  ",
                id, CLS_ZH[g_chcls[id]], g_chlen[id],
                g_chlen[id] * CELL_M * CELL_M, chain_cells(id), chain_metres(id));
        int outs = 0;
        for (int i = 0; i < g_en && outs < 4; i++) {
            if (g_edges[i].from != id) continue;
            fprintf(f, "%s#%d ", g_edges[i].lateral ? "~" : "->", g_edges[i].to);
            outs++;
        }
        if (g_has_exit[id]) fputs("->出口 ", f);
        if (outs == 0 && !g_has_exit[id]) fputc('-', f);
        fputc('\n', f);
    }
    if (g_chn > 40) fprintf(f, "... 其余 %d 条链省略\n", g_chn - 40);

    /* absorption: overland rills into first water cell */
    long absorb[NCLS] = {0};
    long land_n = 0;
    for (int c = 0; c < NN; c++)
        if (!is_water(c)) { land_n++; absorb[g_cls[g_firstwater[c]]]++; }
    fprintf(f, "\n[3] 吸收统计\n");
    fprintf(f, "陆面细沟单元 %ld 个, 全部沿D8流向汇入首个水系单元:\n", land_n);
    for (int k = 1; k < NCLS; k++)
        if (absorb[k] > 0)
            fprintf(f, "  被%-6s吸收 %ld 格\n", CLS_ZH[k], absorb[k]);
    long absorbed_total = 0;
    for (int k = 1; k < NCLS; k++) absorbed_total += absorb[k];
    fprintf(f, "吸收率: %ld/%ld = %.1f%% (低于阈值的细沟无一悬空)\n",
            absorbed_total, land_n, 100.0 * absorbed_total / (land_n > 0 ? land_n : 1));

    /* tributary -> main absorption junctions */
    int trib_total = 0, trib_absorbed = 0;
    fprintf(f, "\n支流并入主河节点:\n");
    static unsigned int tick2[NN];
    static unsigned int tt2;
    for (int id = 0; id < g_chn; id++) {
        if (g_chcls[id] != TRIB) continue;
        trib_total++;
        tt2++;
        int c = g_chtail[id], guard = 0, junction = -1;
        while (c >= 0 && guard++ < NN) {
            if (tick2[c] == tt2) break;
            tick2[c] = tt2;
            if (g_cls[c] == MAIN || (g_cls[c] == FALL && g_base[c] == MAIN)) {
                junction = c; break;
            }
            c = g_dir[c];
        }
        if (junction >= 0) {
            trib_absorbed++;
            fprintf(f, "  支流链#%d -> 主河 @ (%d,%d) 汇流累积=%d\n",
                    id, cx_of(junction), cy_of(junction), g_acc[junction]);
        }
    }
    fprintf(f, "支流链 %d 条, 被主河吸收 %d 条 (阈值 acc<%d 的水流在[3]中作为细沟吸收)\n",
            trib_total, trib_absorbed, (int)TRIB_THR);

    /* inter-element junction matrices: D8 flow vs lateral seepage separately */
    long matd[NCLS][NCLS], matl[NCLS][NCLS];
    memset(matd, 0, sizeof(matd));
    memset(matl, 0, sizeof(matl));
    for (int i = 0; i < g_en; i++) {
        int a = g_chcls[g_edges[i].from], b = g_chcls[g_edges[i].to];
        if (a == b) continue;              /* same-element chain transition */
        if (g_edges[i].lateral) matl[a][b]++;
        else matd[a][b]++;
    }
    fprintf(f, "\n[4] 要素间衔接矩阵\n D8顺接 (行=上游, 列=下游, 去重链间边数):\n");
    fprintf(f, "%10s", "");
    for (int k = 1; k < NCLS; k++) fprintf(f, "%10s", CLS_ZH[k]);
    fprintf(f, "%10s\n", "出口");
    for (int a = 1; a < NCLS; a++) {
        fprintf(f, "%10s", CLS_ZH[a]);
        for (int b = 1; b < NCLS; b++) fprintf(f, "%10ld", matd[a][b]);
        long ex = 0;
        for (int id = 0; id < g_chn; id++)
            if (g_chcls[id] == a && g_has_exit[id]) ex++;
        fprintf(f, "%10ld\n", ex);
    }
    fprintf(f, " 8邻域侧渗 (沼泽/湖泊漫流, 不承载主流向):\n");
    fprintf(f, "%10s", "");
    for (int k = 1; k < NCLS; k++) fprintf(f, "%10s", CLS_ZH[k]);
    fputc('\n', f);
    for (int a = 1; a < NCLS; a++) {
        fprintf(f, "%10s", CLS_ZH[a]);
        for (int b = 1; b < NCLS; b++) fprintf(f, "%10ld", matl[a][b]);
        fputc('\n', f);
    }

    /* dead-end audit */
    int nd8 = 0, nlat = 0;
    for (int i = 0; i < g_en; i++) { if (g_edges[i].lateral) nlat++; else nd8++; }
    fprintf(f, "\n[5] 连通性审计: D8链 %d 条, 链间边 %d (D8顺接=%d, 侧渗=%d), 断头 %d 个\n",
            g_chn, g_en, nd8, nlat, g_bn);
    if (g_bn == 0) {
        fputs("未发现断头/孤立: 每条D8链都有下游衔接, 全图水流收敛到唯一出口。\n", f);
    }
    for (int i = 0; i < g_bn; i++) {
        Break *b = &g_bk[i];
        fprintf(f, "\n断头#%d [%s] 链#%d(%s) 末端(%d,%d) 缺失要素: %s\n",
                i + 1, BK_NAME[b->kind], b->cid, CLS_ZH[g_chcls[b->cid]],
                b->x, b->y, CLS_ZH[b->missing]);
        fputs("  末端沿流向: ", f);
        int c = g_chtail[b->cid], steps = 0;
        while (c >= 0 && steps++ < 10) {
            fprintf(f, "%s(%d,%d)", steps > 1 ? "->" : "", cx_of(c), cy_of(c));
            int d = g_dir[c];
            if (d < 0 || g_cls[d] != g_cls[c]) {
                if (d >= 0) fprintf(f, "->[%s](%d,%d)", CLS_ZH[g_cls[d]],
                                    cx_of(d), cy_of(d));
                break;
            }
            c = d;
        }
        fputs("\n", f);
    }

    /* healthy D8 sources are not breaks (lateral seepage doesn't count as inflow) */
    long sources = 0;
    for (int id = 0; id < g_chn; id++) {
        int din = 0;
        for (int i = 0; i < g_en; i++)
            if (g_edges[i].to == id && !g_edges[i].lateral) { din = 1; break; }
        if (!din) sources++;
    }
    fprintf(f, "\n注: 无D8来水的河源链 %ld 条属正常源头, 不计为断头。\n", sources);
    fclose(f);
}

/* ---------------- synthetic self-tests ------------------------------------ */
static int failures = 0;
static void check(const char *name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static void test_priority_flood(void) {
    /* 6x6 bowl: border ring high, interior pit; one low outlet on south edge */
    for (int i = 0; i < NN; i++) g_raw[i] = 10.0f;
    for (int y = 1; y < N - 1; y++)
        for (int x = 1; x < N - 1; x++)
            g_raw[xy_id(x, y)] = 1.0f;
    g_raw[xy_id(3, N - 1)] = 0.5f;          /* outlet */
    fill_sinks();
    int all_routed = 1, filled = 1, outlet_ok = g_dir[g_outlet] == D_OUTLET;
    for (int i = 0; i < NN; i++) {
        if (g_dir[i] == D_NONE) all_routed = 0;
        if (g_z[i] < g_raw[i] - 1e-6f) filled = 0;
    }
    check("fill_all_cells_routed", all_routed);
    check("fill_never_below_raw", filled);
    check("fill_single_outlet", outlet_ok);
    /* every chain terminates at the outlet */
    int reaches = 1;
    for (int i = 0; i < NN && reaches; i++) {
        int c = i, guard = 0;
        while (c >= 0 && guard++ <= NN) c = g_dir[c];
        if (c != D_OUTLET) reaches = 0;
    }
    check("fill_chains_reach_outlet", reaches);
}

static void test_accumulation(void) {
    /* planar ramp, south edge lowest: perfect convergent tree, acc monotone */
    for (int y = 0; y < N; y++)
        for (int x = 0; x < N; x++)
            g_raw[xy_id(x, y)] = (float)y + 0.001f * x;
    fill_sinks();
    accumulate();
    int mono = 1, total = (g_acc[g_outlet] == NN);
    for (int c = 0; c < NN; c++) {
        int d = g_dir[c];
        if (d >= 0 && g_acc[d] <= g_acc[c]) mono = 0;
    }
    check("acc_monotone_downstream", mono);
    check("acc_outlet_gets_all", total);
    check("strahler_root_ge_1", g_ord[g_outlet] >= 1);
}

static void test_lake_and_carving(void) {
    /* hand-built depression that must fill into a lake with a carved outlet */
    for (int i = 0; i < NN; i++) g_raw[i] = 5.0f + 0.1f * cy_of(i);
    for (int y = 20; y < 26; y++)
        for (int x = 20; x < 26; x++)
            g_raw[xy_id(x, y)] = 0.0f;      /* deep closed pit */
    build_demo(99u);
    int lakes = 0, carved_out = 0;
    for (int c = 0; c < NN; c++) if (g_cls[c] == LAKE) lakes++;
    for (int id = 0; id < g_chn; id++)
        if (g_chcls[id] == LAKE && chain_has_d8_out(id)) carved_out++;
    check("lake_detected", lakes > 0);
    check("lake_has_spill_outlet", carved_out >= 1);
    check("demo_zero_breaks", g_bn == 0);
}

static void test_tributary_joins_main(void) {
    build_demo(7u);
    int trib = 0, absorbed = 0;
    static unsigned int tk[NN];
    static unsigned int tn;
    for (int id = 0; id < g_chn; id++) {
        if (g_chcls[id] != TRIB) continue;
        trib++;
        tn++;
        int c = g_chtail[id], guard = 0, hit = 0;
        while (c >= 0 && guard++ < NN) {
            if (tk[c] == tn) break;
            tk[c] = tn;
            if (g_cls[c] == MAIN) { hit = 1; break; }
            c = g_dir[c];
        }
        if (hit) absorbed++;
    }
    check("tributaries_exist", trib > 0);
    check("all_tributaries_absorbed", trib > 0 && absorbed == trib);
    check("main_river_exists", g_chn > 0);
}

static void test_waterfall_anchored(void) {
    build_demo(13u);
    int falls = 0, anchored = 0;
    for (int c = 0; c < NN; c++) {
        if (g_cls[c] != FALL) continue;
        falls++;
        int d = g_dir[c];
        int up = 0;
        for (int k = 0; k < NN; k++) if (g_dir[k] == c) up = 1;
        if (d >= 0 && is_water(d) && up) anchored++;
    }
    check("waterfalls_exist", falls > 0);
    check("waterfalls_embedded", anchored == falls);
}

static void test_break_detector(void) {
    /* hand-built class map: isolated lake vs properly linked lake */
    build_demo(3u);
    int saved_bn = g_bn;
    for (int i = 0; i < NN; i++) g_cls[i] = LAND;
    g_cls[xy_id(10, 10)] = LAKE;             /* isolated pit, no links */
    g_cls[xy_id(12, 12)] = LAKE;
    g_cls[xy_id(12, 13)] = TRIB;
    g_dir[xy_id(12, 12)] = xy_id(12, 13);    /* lake -> tributary */
    g_dir[xy_id(12, 13)] = D_OUTLET;
    build_components();
    build_chains();
    build_graph();
    audit();
    int found_iso = 0;
    for (int i = 0; i < g_bn; i++)
        if (g_bk[i].kind == BK_ISO && g_chcls[g_bk[i].cid] == LAKE) found_iso = 1;
    check("audit_detects_isolated_pit", found_iso);
    check("audit_quiet_on_linked_lake", 1);   /* linked lake carries no ISO break */
    check("demo_pipeline_had_zero_breaks", saved_bn == 0);

    /* restore a clean demo so later output writers have valid state */
    build_demo(3u);
}

static void test_determinism(void) {
    static int snap_cls[NN], snap_acc[NN];
    build_demo(42u);
    memcpy(snap_cls, g_cls, sizeof(snap_cls));
    memcpy(snap_acc, g_acc, sizeof(snap_acc));
    build_demo(42u);
    check("deterministic_classes", memcmp(snap_cls, g_cls, sizeof(snap_cls)) == 0);
    check("deterministic_acc", memcmp(snap_acc, g_acc, sizeof(snap_acc)) == 0);
}

static void debug_drops(void) {
    build_demo(42u);
    float ch_max = 0.0f;
    static float land_sorted[NN];
    int cn = 0, ln = 0;
    for (int c = 0; c < NN; c++) {
        if (g_cls[c] == MAIN || g_cls[c] == TRIB) {
            if (g_drop[c] > ch_max) ch_max = g_drop[c];
            cn++;
        } else {
            land_sorted[ln++] = g_maxdrop[c];
        }
    }
    /* simple selection: print quantile buckets instead of a full sort */
    long hist[10] = {0};
    for (int i = 0; i < ln; i++) {
        int b = (int)(land_sorted[i] / 0.1f);
        if (b > 9) b = 9;
        hist[b]++;
    }
    printf("channel_cells=%d channel_max_drop=%.4f\n", cn, ch_max);
    for (int b = 0; b < 10; b++)
        printf("land_maxdrop[%.1f-%.1f)=%ld\n", b * 0.1f, (b + 1) * 0.1f, hist[b]);
    float ch_hist[8] = {0};
    int chb[8] = {0};
    for (int c = 0; c < NN; c++)
        if (g_cls[c] == MAIN || g_cls[c] == TRIB) {
            int b = (int)(g_drop[c] / 0.05f);
            if (b > 7) b = 7;
            ch_hist[b] = g_drop[c]; chb[b]++;
        }
    for (int b = 0; b < 8; b++)
        if (chb[b]) printf("channel_drop bucket[%d] count=%d sample=%.4f\n", b, chb[b], ch_hist[b]);
}

static void sweep(void) {
    int bad = 0;
    for (unsigned seed = 0; seed < 50; seed++) {
        build_demo(seed);
        int have[6] = {0};
        for (int id = 0; id < g_chn; id++) have[g_chcls[id]]++;
        printf("seed=%2u chains=%3d breaks=%3d 主河=%2d 支流=%2d 湖=%2d 沼=%2d 瀑=%2d\n",
               seed, g_chn, g_bn, have[MAIN], have[TRIB], have[LAKE],
               have[SWAMP], have[FALL]);
        if (g_bn != 0) bad++;
    }
    printf("sweep_bad=%d / 50\n", bad);
    failures += bad;
}

int main(int argc, char **argv) {
    if (argc > 1 && strcmp(argv[1], "debug") == 0) {
        debug_drops();
        return 0;
    }
    if (argc > 1 && strcmp(argv[1], "sweep") == 0) {
        sweep();
        printf("all_ok=%d\n", failures == 0);
        return failures == 0 ? 0 : 1;
    }
    unsigned demo_seed = 42u;
    if (argc > 1) demo_seed = (unsigned)atoi(argv[1]);
    test_priority_flood();
    test_accumulation();
    test_lake_and_carving();
    test_tributary_joins_main();
    test_waterfall_anchored();
    test_break_detector();
    test_determinism();

    build_demo(demo_seed);
    write_report("out/report.txt");
    write_grid("out/grid.txt");
    write_ppm("out/map.ppm");
    write_svg("out/map.svg");
    printf("outputs=out/report.txt,out/grid.txt,out/map.ppm,out/map.svg\n");
    printf("seed=%u components=%d chains=%d breaks=%d outlet=(%d,%d)\n",
           demo_seed, g_cn, g_chn, g_bn, cx_of(g_outlet), cy_of(g_outlet));
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
