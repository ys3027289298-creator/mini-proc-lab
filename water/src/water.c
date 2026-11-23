
/* water: subsidence pools — local minima -> catchments -> coupled depth solve.
 *
 * pipeline: fBm base terrain -> subsidence bowls at 3 amplitudes ->
 * local minima (plateau-merged) -> steepest-descent catchment tracing ->
 * pool level from the mutual constraint between catchment area A and
 * subsidence-shaped capacity:
 *
 *     V(L)      = sum_i max(0, L - h_i)        (depression volume below L)
 *     V(L_pool) = min(RUNOFF * A, V(L_spill))  (demand vs capacity)
 *
 * catchment area A pushes the level up (demand term), subsidence depth caps
 * it (capacity term up to spill level); the solved L gives pool shape and
 * per-cell depth. Per tier the report intersects three views: catchment
 * derivation table, pool shape map (catchment vs water overlay), and depth
 * distribution. A catchment with demand but no wet cell is "missed" -> the
 * run is incomplete and the missed cell coordinates are listed.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define W 24
#define H 24
#define MAXREG 128
#define NCELLS (W * H)
#define EPS 1e-6f
#define WET_EPS 1e-4f
#define RUNOFF 0.6f /* water volume demanded per unit catchment area */

static float base_h[W][H];
static float h[W][H];
static float depth[W][H];

/* --- deterministic value noise, same hash family as noise/ --- */
static unsigned hash2(int x, int y, unsigned seed) {
    unsigned u = (unsigned)x * 374761393u + (unsigned)y * 668265263u + seed * 2246822519u;
    u = (u ^ (u >> 13)) * 1274126177u;
    return u ^ (u >> 16);
}
static float rand01(unsigned u) { return (float)(u & 0xffffff) / 16777215.0f; }
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

static float fbm(float x, float y, int oct, unsigned seed) {
    float sum = 0.0f, amp = 1.0f, norm = 0.0f, freq = 1.0f;
    for (int i = 0; i < oct; i++) {
        sum += value_noise(x * freq, y * freq, seed + (unsigned)i * 17u) * amp;
        norm += amp;
        amp *= 0.5f;
        freq *= 2.0f;
    }
    return norm > 0.0f ? sum / norm : 0.0f;
}

/* --- terrain + subsidence --- */
static void build_base(void) {
    for (int x = 0; x < W; x++)
        for (int y = 0; y < H; y++)
            base_h[x][y] = fbm(x * 0.18f, y * 0.18f, 4, 42u) * 1.2f;
}

static void apply_subsidence(float amp) {
    const float cx[2] = { 7.0f, 17.0f };
    const float cy[2] = { 9.0f, 15.0f };
    const float rr[2] = { 4.0f, 5.0f };
    for (int x = 0; x < W; x++)
        for (int y = 0; y < H; y++) {
            float s = 0.0f;
            for (int k = 0; k < 2; k++) {
                float dx = (float)x - cx[k], dy = (float)y - cy[k];
                s += amp * expf(-(dx * dx + dy * dy) / (rr[k] * rr[k]));
            }
            h[x][y] = base_h[x][y] - s;
            depth[x][y] = 0.0f;
        }
}

/* --- regions: local minima, plateau-merged via union-find --- */
static int region_of[W][H];  /* region id or -1 */
static int catch_id[W][H];   /* region the cell drains to; -1 off-map; -2 unknown */
static int nreg;
static int reg_ax[MAXREG], reg_ay[MAXREG]; /* anchor = lowest cell of region */
static float reg_hmin[MAXREG];
static int catch_n[MAXREG];
static int catch_cells[MAXREG][NCELLS]; /* linear index y*W+x */
static float reg_avail[MAXREG], reg_cap[MAXREG], reg_level[MAXREG];
static float reg_vol[MAXREG], reg_spill[MAXREG], reg_dmax[MAXREG];

static int uf_parent[NCELLS];
static int uf_find(int a) {
    while (uf_parent[a] != a) { uf_parent[a] = uf_parent[uf_parent[a]]; a = uf_parent[a]; }
    return a;
}
static void uf_union(int a, int b) {
    a = uf_find(a); b = uf_find(b);
    if (a != b) uf_parent[b] = a;
}

static const int DX[8] = { -1, 0, 1, -1, 1, -1, 0, 1 };
static const int DY[8] = { -1, -1, -1, 0, 0, 1, 1, 1 };

static int in_bounds(int x, int y) { return x >= 0 && y >= 0 && x < W && y < H; }

static int no_lower(int x, int y) {
    for (int k = 0; k < 8; k++) {
        int nx = x + DX[k], ny = y + DY[k];
        if (!in_bounds(nx, ny)) continue;
        if (h[nx][ny] < h[x][y] - EPS) return 0;
    }
    return 1;
}

static void find_regions(void) {
    static int root_id[NCELLS];
    for (int i = 0; i < NCELLS; i++) { uf_parent[i] = i; root_id[i] = -1; }
    memset(region_of, 0xff, sizeof(region_of));
    /* border minima drain off the map edge: no closed depression, no region */
    for (int x = 1; x < W - 1; x++)
        for (int y = 1; y < H - 1; y++) {
            if (!no_lower(x, y)) continue;
            for (int k = 0; k < 8; k++) {
                int nx = x + DX[k], ny = y + DY[k];
                if (nx < 1 || ny < 1 || nx >= W - 1 || ny >= H - 1) continue;
                if (no_lower(nx, ny) && fabsf(h[nx][ny] - h[x][y]) <= EPS)
                    uf_union(y * W + x, ny * W + nx);
            }
        }
    nreg = 0;
    for (int x = 1; x < W - 1; x++)
        for (int y = 1; y < H - 1; y++) {
            if (!no_lower(x, y)) continue;
            int root = uf_find(y * W + x);
            if (root_id[root] < 0) {
                if (nreg >= MAXREG) break;
                root_id[root] = nreg;
                reg_hmin[nreg] = h[x][y];
                reg_ax[nreg] = x;
                reg_ay[nreg] = y;
                nreg++;
            }
            int id = root_id[root];
            region_of[x][y] = id;
            if (h[x][y] < reg_hmin[id]) {
                reg_hmin[id] = h[x][y];
                reg_ax[id] = x;
                reg_ay[id] = y;
            }
        }
}

/* --- steepest-descent flow tracing with memoization --- */
static int flow_step(int x, int y, int* nx, int* ny) {
    float best = h[x][y] - EPS;
    int bx = -1, by = -1;
    for (int k = 0; k < 8; k++) {
        int cx = x + DX[k], cy = y + DY[k];
        if (!in_bounds(cx, cy)) continue;
        if (h[cx][cy] < best) { best = h[cx][cy]; bx = cx; by = cy; }
    }
    if (bx < 0) return 0;
    *nx = bx; *ny = by;
    return 1;
}

static int trace_cell(int sx, int sy) {
    static int path[NCELLS];
    int n = 0, x = sx, y = sy, result;
    for (;;) {
        if (region_of[x][y] >= 0) { result = region_of[x][y]; break; }
        if (catch_id[x][y] != -2) { result = catch_id[x][y]; break; }
        int nx, ny;
        if (!flow_step(x, y, &nx, &ny)) { result = -1; break; } /* border minimum: off-map */
        path[n++] = y * W + x;
        x = nx; y = ny;
    }
    for (int i = 0; i < n; i++) catch_id[path[i] % W][path[i] / W] = result;
    return result;
}

static void build_catchments(void) {
    for (int x = 0; x < W; x++)
        for (int y = 0; y < H; y++)
            catch_id[x][y] = region_of[x][y] >= 0 ? region_of[x][y] : -2;
    for (int r = 0; r < nreg; r++) catch_n[r] = 0;
    for (int x = 0; x < W; x++)
        for (int y = 0; y < H; y++) {
            int r = trace_cell(x, y);
            if (r >= 0 && catch_n[r] < NCELLS)
                catch_cells[r][catch_n[r]++] = y * W + x;
        }
}

/* --- coupled solve: V(L) = min(RUNOFF * A, V(spill)) --- */
static int cmp_float(const void* a, const void* b) {
    float d = *(const float*)a - *(const float*)b;
    return d < 0 ? -1 : d > 0;
}

static void solve_region(int r) {
    static float hs[NCELLS];
    int n = catch_n[r];
    for (int i = 0; i < n; i++) {
        int idx = catch_cells[r][i];
        hs[i] = h[idx % W][idx / W];
    }
    qsort(hs, (size_t)n, sizeof(float), cmp_float);

    /* spill level: lowest rim crossing, max(inside, outside) per edge */
    float spill = 1e30f;
    for (int i = 0; i < n; i++) {
        int idx = catch_cells[r][i], x = idx % W, y = idx / W;
        if (x == 0 || y == 0 || x == W - 1 || y == H - 1) {
            if (h[x][y] < spill) spill = h[x][y];
        }
        for (int k = 0; k < 8; k++) {
            int nx = x + DX[k], ny = y + DY[k];
            if (!in_bounds(nx, ny) || catch_id[nx][ny] == r) continue;
            float s = fmaxf(h[x][y], h[nx][ny]);
            if (s < spill) spill = s;
        }
    }

    float cap = 0.0f;
    for (int i = 0; i < n; i++)
        if (hs[i] < spill) cap += spill - hs[i];
    float avail = RUNOFF * (float)n;       /* demand term: catchment area A */
    float target = fminf(avail, cap);      /* mutual constraint: demand vs capacity */

    /* walk sorted heights: in segment i, V(L) = L*(i+1) - prefix */
    float prefix = 0.0f, L = hs[0];
    int solved = 0;
    for (int i = 0; i < n; i++) {
        float top = (i + 1 < n && hs[i + 1] < spill) ? hs[i + 1] : spill;
        prefix += hs[i];
        if (top * (i + 1) - prefix >= target - 1e-9f) {
            L = (target + prefix) / (float)(i + 1);
            solved = 1;
            break;
        }
    }
    if (!solved || L > spill) L = spill;

    float vol = 0.0f, dmax = 0.0f;
    for (int i = 0; i < n; i++) {
        int idx = catch_cells[r][i], x = idx % W, y = idx / W;
        float d = L - h[x][y];
        if (d < 0.0f) d = 0.0f;
        depth[x][y] = d;
        vol += d;
        if (d > dmax) dmax = d;
    }
    reg_avail[r] = avail;
    reg_cap[r] = cap;
    reg_level[r] = L;
    reg_vol[r] = vol;
    reg_spill[r] = spill;
    reg_dmax[r] = dmax;
}

/* --- report: derivation table, shape map, depth distribution --- */
static void print_derivation(void) {
    printf("catchment derivation (A -> demand, capacity from subsidence, coupled level):\n");
    for (int r = 0; r < nreg; r++)
        printf("  region %2d anchor=(%2d,%2d) hmin=%+.3f catch=%3d area=%6.1f"
               " | demand=%7.2f cap=%7.2f spill=%+.3f"
               " | level=%+.3f vol=%7.2f dmax=%.3f\n",
               r, reg_ax[r], reg_ay[r], reg_hmin[r], catch_n[r], (float)catch_n[r],
               reg_avail[r], reg_cap[r], reg_spill[r],
               reg_level[r], reg_vol[r], reg_dmax[r]);
}

static void print_map(float tier_dmax) {
    printf("pool shape map: '.' outside catchment, ',' dry catchment, '1'-'9' water depth tier\n");
    printf("    ");
    for (int x = 0; x < W; x++) printf("%d", x % 10);
    printf("\n");
    for (int y = 0; y < H; y++) {
        printf("%3d ", y);
        for (int x = 0; x < W; x++) {
            if (depth[x][y] > WET_EPS && tier_dmax > 0.0f) {
                int t = (int)ceilf(depth[x][y] / tier_dmax * 9.0f);
                if (t < 1) t = 1;
                if (t > 9) t = 9;
                putchar('0' + t);
            } else if (catch_id[x][y] >= 0) {
                putchar(',');
            } else {
                putchar('.');
            }
        }
        printf("\n");
    }
}

static void print_depth_dist(float tier_dmax) {
    int wet = 0, bins[6] = { 0 };
    float sum = 0.0f;
    for (int x = 0; x < W; x++)
        for (int y = 0; y < H; y++) {
            if (depth[x][y] <= WET_EPS) continue;
            wet++;
            sum += depth[x][y];
            int b = tier_dmax > 0.0f ? (int)(depth[x][y] / tier_dmax * 6.0f) : 0;
            if (b > 5) b = 5;
            bins[b]++;
        }
    printf("depth distribution: wet=%d mean=%.3f max=%.3f\n", wet,
           wet ? sum / wet : 0.0f, tier_dmax);
    for (int b = 0; b < 6; b++) {
        float lo = tier_dmax * b / 6.0f, hi = tier_dmax * (b + 1) / 6.0f;
        printf("  [%5.2f,%5.2f) %4d ", lo, hi, bins[b]);
        for (int i = 0; i < bins[b] && i < 40; i++) putchar('#');
        printf("\n");
    }
}

/* --- missed-catchment detector: demand > 0 but no wet cell -> incomplete --- */
static int find_missed(void) {
    int missed = 0;
    for (int r = 0; r < nreg; r++) {
        if (reg_avail[r] <= WET_EPS) continue;
        int wet = 0;
        for (int i = 0; i < catch_n[r]; i++) {
            int idx = catch_cells[r][i];
            if (depth[idx % W][idx / W] > WET_EPS) wet++;
        }
        if (wet == 0) {
            missed++;
            printf("missed catchment region=%d anchor=(%d,%d) cells=%d coords:",
                   r, reg_ax[r], reg_ay[r], catch_n[r]);
            for (int i = 0; i < catch_n[r]; i++) {
                if (i == 40) { printf(" ..."); break; }
                int idx = catch_cells[r][i];
                printf(" (%d,%d)", idx % W, idx / W);
            }
            printf("\n");
        }
    }
    return missed;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    build_base();
    const float amps[3] = { 0.5f, 1.5f, 3.0f };
    const char* names[3] = { "light", "mid", "heavy" };
    int total_missed = 0;
    float tier_dmax[3] = { 0 }, tier_vol[3] = { 0 };
    static float saved_level[MAXREG];
    int saved_nreg = 0;

    for (int t = 0; t < 3; t++) {
        apply_subsidence(amps[t]);
        find_regions();
        build_catchments();
        for (int r = 0; r < nreg; r++) solve_region(r);

        for (int r = 0; r < nreg; r++) {
            if (reg_dmax[r] > tier_dmax[t]) tier_dmax[t] = reg_dmax[r];
            tier_vol[t] += reg_vol[r];
        }
        printf("=== tier %s amp=%.1f: regions=%d total_vol=%.2f max_depth=%.3f\n",
               names[t], amps[t], nreg, tier_vol[t], tier_dmax[t]);
        print_derivation();
        print_map(tier_dmax[t]);
        print_depth_dist(tier_dmax[t]);
        int m = find_missed();
        if (m > 0) printf("tier %s INCOMPLETE: %d missed catchment(s)\n", names[t], m);
        total_missed += m;
        if (t == 0) {
            saved_nreg = nreg;
            for (int r = 0; r < nreg; r++) saved_level[r] = reg_level[r];
        }
    }

    /* determinism: rerun tier light, levels must match */
    apply_subsidence(amps[0]);
    find_regions();
    build_catchments();
    for (int r = 0; r < nreg; r++) solve_region(r);
    int det = nreg == saved_nreg;
    for (int r = 0; r < nreg && det; r++)
        if (fabsf(reg_level[r] - saved_level[r]) > 1e-5f) det = 0;
    check("deterministic", det);

    /* mutual constraint residual: pooled volume == min(demand, capacity) */
    int res_ok = 1;
    for (int r = 0; r < nreg; r++) {
        float want = fminf(reg_avail[r], reg_cap[r]);
        if (fabsf(reg_vol[r] - want) > 1e-3f * (1.0f + want)) res_ok = 0;
    }
    check("constraint_volume", res_ok);

    /* detector self-test: wipe region 0's water, expect it flagged with coords */
    printf("self-test: expect one missed report below\n");
    if (nreg > 0) {
        for (int i = 0; i < catch_n[0]; i++) {
            int idx = catch_cells[0][i];
            depth[idx % W][idx / W] = 0.0f;
        }
        check("detector_flags_missed", find_missed() >= 1);
    }

    check("no_missed_catchments", total_missed == 0);
    check("pools_nonempty", tier_vol[0] > 0 && tier_vol[1] > 0 && tier_vol[2] > 0);
    check("deeper_subsidence_deeper_pool", tier_dmax[2] > tier_dmax[0]);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
