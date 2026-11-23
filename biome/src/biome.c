/* biome: climate-driven biome map with a controlled boundary transition width.
 *
 * A cell has climate (T, H); biome i sits at (T_i, H_i).  The score is the
 * weighted squared climate distance
 *
 *     s_i = -( w_t (T - T_i)^2 + w_h (H - H_i)^2 ),      w_t + w_h = 1
 *
 * and the map is the softmax of s / tau.  tau alone sets the transition
 * width; the weights decide which climate axis controls it.
 *
 * Along the boundary between biomes a and b the margin is m = s_a - s_b and
 * the two dominant classes give P_a = sigmoid(m / tau).  P_a runs 0.9 -> 0.1
 * while m runs over 2 ln(9) tau, so the perpendicular width is
 *
 *     W = 2 ln(9) tau / |grad m|,   grad m = 2 w_t dT grad T + 2 w_h dH grad H.   (1)
 *
 * That is the weight <=> width relation: the temperature part of |grad m| is
 * 2 w_t dT |grad T|, so a temperature-controlled boundary has W ~ 1/w_t;
 * raising an axis's weight narrows the boundaries it controls.
 *
 * Derivation of the weights: over the boundary set, require the two RMS
 * contributions to |grad m| to be equal,
 *     w_t rms(2 dT |grad T|) = w_h rms(2 dH |grad H|)  =>  w_t : w_h = B : A.
 * The boundary set moves with the weights, so this is iterated to a fixed
 * point.  With the weights fixed, tau = W |grad m| / (2 ln 9) is calibrated
 * until the measured median width equals the target.
 *
 * Runs: hard baseline (equal weights), hard derived, W = 3, 10, 24, and an
 * equal-weight W = 10 control.  Every boundary cell whose transition closes
 * inside 1.5 cells is listed with coordinates.
 */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>

#define N 256
#define NB 6
#define MAXCROSS (N * (N - 1) * 2)

static const double LN9 = 2.1972245773362196;
static const double K_WIDTH = 2.0 * 2.1972245773362196;   /* 2 ln 9 */

typedef struct { const char* name; float t, h; } Anchor;

/* anchors in normalized climate space: temperature 0..1, humidity 0..1 */
static const Anchor BIOMES[NB] = {
    { "ice",            0.10f, 0.22f },
    { "taiga",          0.26f, 0.58f },
    { "grassland",      0.48f, 0.18f },
    { "temperate_rain", 0.58f, 0.82f },
    { "desert",         0.90f, 0.10f },
    { "tropical_rain",  0.90f, 0.72f },
};

/* ---------- climate fields (same value-noise / fBm as ../noise) ---------- */

static unsigned hash2(int x, int y, unsigned seed) {
    unsigned h = (unsigned)x * 374761393u + (unsigned)y * 668265263u + seed * 2246822519u;
    h = (h ^ (h >> 13)) * 1274126177u;
    return h ^ (h >> 16);
}

static float rand01(unsigned h) { return (float)(h & 0xffffff) / 16777215.0f; }
static float lerpf(float a, float b, float t) { return a + (b - a) * t; }
static float fade(float t) { return t * t * (3.0f - 2.0f * t); }

static float value_noise(float x, float y, unsigned seed) {
    int xi = (int)floorf(x), yi = (int)floorf(y);
    float fx = x - xi, fy = y - yi;
    float v00 = rand01(hash2(xi, yi, seed));
    float v10 = rand01(hash2(xi + 1, yi, seed));
    float v01 = rand01(hash2(xi, yi + 1, seed));
    float v11 = rand01(hash2(xi + 1, yi + 1, seed));
    float u = fade(fx), v = fade(fy);
    return lerpf(lerpf(v00, v10, u), lerpf(v01, v11, u), v) * 2.0f - 1.0f;
}

static float fbm(float x, float y, int octaves, float lacunarity, float gain, unsigned seed) {
    float sum = 0.0f, amp = 1.0f, norm = 0.0f, freq = 1.0f;
    for (int i = 0; i < octaves; i++) {
        sum += value_noise(x * freq, y * freq, seed + (unsigned)i * 17u) * amp;
        norm += amp;
        amp *= gain;
        freq *= lacunarity;
    }
    return norm > 0.0f ? sum / norm : 0.0f;
}

/* ---------- state ---------- */

static float T[N][N], H[N][N];
static float gTx[N][N], gTy[N][N], gHx[N][N], gHy[N][N];
static float S[NB][N][N];          /* scores */
static float P[NB][N][N];          /* softmax probabilities */
static float G[N][N];              /* |grad P| */
static float g_pair_q50[NB][NB];      /* median |grad m| by family, W-independent */

static float clampf(float v, float lo, float hi) { return v < lo ? lo : v > hi ? hi : v; }

static void normalize01(float f[N][N]) {
    float lo = 1e9f, hi = -1e9f;
    for (int y = 0; y < N; y++) for (int x = 0; x < N; x++) {
        if (f[y][x] < lo) lo = f[y][x];
        if (f[y][x] > hi) hi = f[y][x];
    }
    float k = (hi > lo) ? 0.94f / (hi - lo) : 0.0f;
    for (int y = 0; y < N; y++) for (int x = 0; x < N; x++)
        f[y][x] = 0.03f + (f[y][x] - lo) * k;
}

static void gen_fields(void) {
    /* temperature is the slow axis, humidity varies about 2.4x faster */
    for (int y = 0; y < N; y++) for (int x = 0; x < N; x++) {
        T[y][x] = fbm(x * 0.020f, y * 0.020f, 4, 2.0f, 0.5f, 101u);
        H[y][x] = fbm(x * 0.048f + 31.7f, y * 0.048f + 17.3f, 4, 2.0f, 0.5f, 227u);
    }
    normalize01(T);
    normalize01(H);
}

static void compute_gradients(void) {
    for (int y = 0; y < N; y++) for (int x = 0; x < N; x++) {
        int xm = x > 0 ? x - 1 : x, xp = x < N - 1 ? x + 1 : x;
        int ym = y > 0 ? y - 1 : y, yp = y < N - 1 ? y + 1 : y;
        gTx[y][x] = (T[y][xp] - T[y][xm]) / (float)(xp - xm);
        gTy[y][x] = (T[yp][x] - T[ym][x]) / (float)(yp - ym);
        gHx[y][x] = (H[y][xp] - H[y][xm]) / (float)(xp - xm);
        gHy[y][x] = (H[yp][x] - H[ym][x]) / (float)(yp - ym);
    }
}

static void compute_scores(float wt, float wh) {
    for (int y = 0; y < N; y++) for (int x = 0; x < N; x++) {
        float t = T[y][x], h = H[y][x];
        for (int i = 0; i < NB; i++) {
            float dt = t - BIOMES[i].t, dh = h - BIOMES[i].h;
            S[i][y][x] = -(wt * dt * dt + wh * dh * dh);
        }
    }
}

static int label_at(int x, int y) {
    int best = 0;
    for (int i = 1; i < NB; i++) if (S[i][y][x] > S[best][y][x]) best = i;
    return best;
}

static void hard_labels(int lab[N][N]) {
    for (int y = 0; y < N; y++) for (int x = 0; x < N; x++) lab[y][x] = label_at(x, y);
}

static void onehot_from_labels(const int lab[N][N]) {
    for (int y = 0; y < N; y++) for (int x = 0; x < N; x++)
        for (int i = 0; i < NB; i++) P[i][y][x] = (i == lab[y][x]) ? 1.0f : 0.0f;
}

static void softmax(float tau) {
    for (int y = 0; y < N; y++) for (int x = 0; x < N; x++) {
        float mx = S[0][y][x];
        for (int i = 1; i < NB; i++) if (S[i][y][x] > mx) mx = S[i][y][x];
        float sum = 0.0f;
        for (int i = 0; i < NB; i++) {
            P[i][y][x] = expf((S[i][y][x] - mx) / tau);
            sum += P[i][y][x];
        }
        for (int i = 0; i < NB; i++) P[i][y][x] /= sum;
    }
}

/* ---------- transition width along the boundary normal ---------------- */

/* ---------- small stats helpers ---------- */

static int cmp_f(const void* a, const void* b) {
    float x = *(const float*)a, y = *(const float*)b;
    return (x > y) - (x < y);
}

static float quant(const float* v, int n, float p) {
    if (n <= 0) return 0.0f;
    int k = (int)floorf(p * (n - 1) + 0.5f);
    if (k < 0) k = 0;
    if (k >= n) k = n - 1;
    return v[k];
}

static int boundary_neighbor(const int lab[N][N], int x, int y) {
    int a = lab[y][x];
    if (x + 1 < N && lab[y][x + 1] != a) return lab[y][x + 1];
    if (x > 0 && lab[y][x - 1] != a) return lab[y][x - 1];
    if (y + 1 < N && lab[y + 1][x] != a) return lab[y + 1][x];
    if (y > 0 && lab[y - 1][x] != a) return lab[y - 1][x];
    return -1;
}

static float margin_grad(int a, int b, int x, int y) {
    int xm = x > 0 ? x - 1 : x, xp = x < N - 1 ? x + 1 : x;
    int ym = y > 0 ? y - 1 : y, yp = y < N - 1 ? y + 1 : y;
    float dx = ((S[a][y][xp] - S[b][y][xp]) - (S[a][y][xm] - S[b][y][xm])) / (float)(xp - xm);
    float dy = ((S[a][yp][x] - S[b][yp][x]) - (S[a][ym][x] - S[b][ym][x])) / (float)(yp - ym);
    return sqrtf(dx * dx + dy * dy);
}

/* ---------- weight derivation ---------------------------------------- */

/* RMS of the two axis contributions to |grad m| over boundary cells:
 *   a = 2 |dT| |grad T|,  b = 2 |dH| |grad H|. */
static void boundary_terms(float wt, float wh, double* A, double* B, long* nout) {
    int lab[N][N];
    compute_scores(wt, wh);
    hard_labels(lab);
    double sa = 0.0, sb = 0.0;
    long n = 0;
    for (int y = 0; y < N; y++) for (int x = 0; x < N; x++) {
        int b = boundary_neighbor(lab, x, y);
        if (b < 0) continue;
        int a = lab[y][x];
        double dT = fabs((double)BIOMES[a].t - BIOMES[b].t);
        double dH = fabs((double)BIOMES[a].h - BIOMES[b].h);
        double gT = hypot(gTx[y][x], gTy[y][x]);
        double gH = hypot(gHx[y][x], gHy[y][x]);
        sa += 4.0 * dT * dT * gT * gT;
        sb += 4.0 * dH * dH * gH * gH;
        n++;
    }
    *A = n ? sqrt(sa / (double)n) : 0.0;
    *B = n ? sqrt(sb / (double)n) : 0.0;
    *nout = n;
}

/* ---------- transition width along the boundary normal ---------------- */

static float m_bilin(int a, int b, float fx, float fy);

/* signed margin at the anchor position shifted by t0 along the normal */
static float margin_at0(int a, int b, float x, float y, float nx, float ny, float t0) {
    return m_bilin(a, b, x + t0 * nx, y + t0 * ny);
}

/* bilinear sample of the margin m = s_a - s_b at a fractional position */
static float m_bilin(int a, int b, float fx, float fy) {
    fx = clampf(fx, 0.0f, (float)(N - 2));
    fy = clampf(fy, 0.0f, (float)(N - 2));
    int x0 = (int)floorf(fx), y0 = (int)floorf(fy);
    float tx = fx - x0, ty = fy - y0;
    int x1 = x0 + 1 < N ? x0 + 1 : x0;
    int y1 = y0 + 1 < N ? y0 + 1 : y0;
    float m00 = S[a][y0][x0] - S[b][y0][x0];
    float m10 = S[a][y0][x1] - S[b][y0][x1];
    float m01 = S[a][y1][x0] - S[b][y1][x0];
    float m11 = S[a][y1][x1] - S[b][y1][x1];
    return (m00 * (1 - tx) + m10 * tx) * (1 - ty) + (m01 * (1 - tx) + m11 * tx) * ty;
}

/* the bilinear field is only reliable one cell in from the border */
static int safe(float fx, float fy) {
    return fx >= 1.0f && fy >= 1.0f && fx <= (float)(N - 2) && fy <= (float)(N - 2);
}

/* march from the anchor to the first |m| > M, then bisect that last step for
 * the exact crossing.  Returns the distance in cells; on reaching the map edge
 * (or the cap) first, returns the distance travelled and sets *ok = 0.  The
 * returned value is then a lower bound on the true band edge. */
static float band_edge(int a, int b, float x, float y, float nx, float ny,
                       float t0, int sg, float M, float cap, int* ok) {
    const float step = 0.25f;
    float prev_t = 0.0f, prev_m = margin_at0(a, b, x, y, nx, ny, t0);
    *ok = 0;
    for (float t = step; t <= cap; t += step) {
        float fx = x + (t0 + sg * t) * nx, fy = y + (t0 + sg * t) * ny;
        if (!safe(fx, fy)) return t - step;
        float m = m_bilin(a, b, fx, fy);
        if (fabsf(m) > M) {
            float lo = prev_t, hi = t;
            for (int it = 0; it < 20; it++) {
                float mid = 0.5f * (lo + hi);
                float fm = m_bilin(a, b, x + (t0 + sg * mid) * nx, y + (t0 + sg * mid) * ny);
                if (fabsf(fm) > M) hi = mid; else lo = mid;
            }
            *ok = 1;
            return 0.5f * (lo + hi);
        }
        prev_t = t; prev_m = m;
    }
    (void)prev_m;
    return cap;
}

/* perpendicular 10-90 width of the a/b transition, in cells.
 *
 * The band is |m| <= M with M = ln(9) tau: inside it the two-class share
 * P_a/(P_a+P_b) lies between 0.1 and 0.9.  The margin's zero on the normal
 * through the anchor is t0 = -m0 / |grad m| (linear anchor); each side is then
 * marched out to |m| = M.  For a locally linear margin this is exactly eq. (1):
 * 2 ln9 tau / |grad m|.
 *
 * *closed = 2 when both sides reached |m| = M inside the map.  Otherwise the
 * returned width stops at the map edge and is a lower bound.  *capped always
 * returns that lower-bound width, which is what the calibration bisects on:
 * it is monotone in tau even when the band runs off the map. */
static float measure_width(int a, int b, int x, int y, float tau, int* closed, float* capped) {
    float q = margin_grad(a, b, x, y);
    *closed = 0;
    *capped = 0.0f;
    if (q < 1e-9f) return 0.0f;
    int xm = x > 0 ? x - 1 : x, xp = x < N - 1 ? x + 1 : x;
    int ym = y > 0 ? y - 1 : y, yp = y < N - 1 ? y + 1 : y;
    float gx = ((S[a][y][xp] - S[b][y][xp]) - (S[a][y][xm] - S[b][y][xm])) / (float)(xp - xm);
    float gy = ((S[a][yp][x] - S[b][yp][x]) - (S[a][ym][x] - S[b][ym][x])) / (float)(yp - ym);
    float nx = gx / q, ny = gy / q;
    float t0 = -(S[a][y][x] - S[b][y][x]) / q;
    /* linear anchoring is only meaningful if the zero is within a few cells;
     * a tangential sliver has |grad m| ~ 0 and t0 blows up */
    if (fabsf(t0) > 6.0f) { *closed = 3; return 0.0f; }
    float M = (float)(LN9 * tau);
    float cap = 12.0f * M / q + 8.0f;
    int ok0 = 0, ok1 = 0;
    float d0 = band_edge(a, b, (float)x, (float)y, nx, ny, t0, +1, M, cap, &ok0);
    float d1 = band_edge(a, b, (float)x, (float)y, nx, ny, t0, -1, M, cap, &ok1);
    *capped = d0 + d1;
    *closed = (ok0 && ok1) ? 2 : (ok0 || ok1) ? 1 : 0;
    return *closed == 2 ? d0 + d1 : 0.0f;
}

/* ---------- per-run bookkeeping -------------------------------------- */

typedef struct {
    unsigned short x, y;
    unsigned char axis, a, b, kind;   /* 0 = cut, 1 = slope step */
    float width, dp;
} StepRec;

typedef struct {
    const char* tag;
    int run_index;
    float wt, wh, W, tau;
    int hard;
    long crossings, censored, steps, near_target, slivers, slope_steps;
    float w_mean, w_p05, w_p50, w_p95, w_cv;
    float q50, pred_w50, frac_hit, frac_closed, corr;
    long pair_diluted;
    long pair_n[NB][NB], pair_compressed[NB][NB], pair_sliver[NB][NB];
    double pair_wmed[NB][NB];
    float dp50, dp95, dp_max;
    float gb_mean, gb_p95, gb_max, g_theory;
    int gx_worst, gy_worst;
    long boundary_cells, band_cells;
    long hard_area[NB];
    double soft_area[NB];
    StepRec* step;
    long step_count, step_cap;
} Run;

static void step_push(Run* r, int x, int y, int axis, int a, int b, float w, float dp, int kind) {
    if (r->step_count >= r->step_cap) {
        r->step_cap = r->step_cap ? r->step_cap * 2 : 256;
        r->step = realloc(r->step, (size_t)r->step_cap * sizeof(StepRec));
        if (!r->step) { perror("realloc"); exit(1); }
    }
    StepRec* s = &r->step[r->step_count++];
    s->x = (unsigned short)x; s->y = (unsigned short)y;
    s->axis = (unsigned char)axis; s->a = (unsigned char)a; s->b = (unsigned char)b;
    s->kind = (unsigned char)kind;
    s->width = w; s->dp = dp;
}

static void gradient_field(void) {
    for (int y = 0; y < N; y++) for (int x = 0; x < N; x++) {
        int xm = x > 0 ? x - 1 : x, xp = x < N - 1 ? x + 1 : x;
        int ym = y > 0 ? y - 1 : y, yp = y < N - 1 ? y + 1 : y;
        float g2 = 0.0f;
        for (int i = 0; i < NB; i++) {
            float dx = (P[i][y][xp] - P[i][y][xm]) / (float)(xp - xm);
            float dy = (P[i][yp][x] - P[i][ym][x]) / (float)(yp - ym);
            g2 += dx * dx + dy * dy;
        }
        G[y][x] = sqrtf(g2);
    }
}

static void boundary_gradient_stats(const int lab[N][N], Run* r) {
    static float v[N * N];
    int n = 0;
    double sum = 0.0;
    for (int y = 0; y < N; y++) for (int x = 0; x < N; x++) {
        if (boundary_neighbor(lab, x, y) < 0) continue;
        v[n++] = G[y][x];
        sum += G[y][x];
    }
    qsort(v, (size_t)n, sizeof(float), cmp_f);
    r->gb_mean = n ? (float)(sum / n) : 0.0f;
    r->gb_p95 = quant(v, n, 0.95f);
    r->gb_max = quant(v, n, 1.0f);
    r->gx_worst = r->gy_worst = 0;
    for (int y = 0; y < N; y++) for (int x = 0; x < N; x++) {
        if (boundary_neighbor(lab, x, y) < 0) continue;
        if (G[y][x] == r->gb_max) { r->gx_worst = x; r->gy_worst = y; }
    }
    /* peak |grad P| of a smooth two-class ramp of perpendicular width W:
     * dP/ds = P(1-P)/tau and ds/dx = |grad m|, with W = 2 ln9 tau/|grad m|;
     * at P = 0.5 the slope is |grad m|/(4 tau) = ln9/(2 W) ~ 1.1/W. */
    r->g_theory = r->W > 0.0f ? (float)(LN9 / (2.0 * r->W)) : 0.0f;
}

static void area_counts(Run* r, const int lab[N][N]) {
    memset(r->hard_area, 0, sizeof(r->hard_area));
    memset(r->soft_area, 0, sizeof(r->soft_area));
    r->boundary_cells = 0;
    r->band_cells = 0;
    for (int y = 0; y < N; y++) for (int x = 0; x < N; x++) {
        r->hard_area[lab[y][x]]++;
        for (int i = 0; i < NB; i++) r->soft_area[i] += P[i][y][x];
        float pmax = P[0][y][x];
        for (int i = 1; i < NB; i++) if (P[i][y][x] > pmax) pmax = P[i][y][x];
        if (pmax < 0.9f) r->band_cells++;
        if (boundary_neighbor(lab, x, y) >= 0) r->boundary_cells++;
    }
}

static float pearson(const float* x, const float* y, int n) {
    if (n < 2) return 0.0f;
    double mx = 0.0, my = 0.0;
    for (int i = 0; i < n; i++) { mx += x[i]; my += y[i]; }
    mx /= n; my /= n;
    double sxy = 0.0, sxx = 0.0, syy = 0.0;
    for (int i = 0; i < n; i++) {
        double dx = x[i] - mx, dy = y[i] - my;
        sxy += dx * dy; sxx += dx * dx; syy += dy * dy;
    }
    return (sxx > 1e-12 && syy > 1e-12) ? (float)(sxy / sqrt(sxx * syy)) : 0.0f;
}

static void pair_compressed_add(Run* r, int a, int b) {
    if (a > b) { int t = a; a = b; b = t; }
    r->pair_compressed[a][b]++;
}

typedef struct {
    unsigned short x, y;
    unsigned char axis, a, b, closed;
    float w, q, dp;
} Cand;

static void collect(Run* r, const int lab[N][N]) {
    static Cand cand[MAXCROSS];
    static float qall[N * N];
    int nc = 0, nq = 0;

    r->crossings = r->censored = r->slivers = 0;
    r->dp_max = 0.0f;
    r->step_count = 0;
    r->pair_diluted = 0;

    /* median |grad m| over every boundary pair, used by the step test */
    if (!r->hard) {
        for (int y = 0; y < N; y++) for (int x = 0; x < N - 1; x++)
            if (lab[y][x] != lab[y][x + 1]) qall[nq++] = margin_grad(lab[y][x], lab[y][x + 1], x, y);
        for (int x = 0; x < N; x++) for (int y = 0; y < N - 1; y++)
            if (lab[y][x] != lab[y + 1][x]) qall[nq++] = margin_grad(lab[y][x], lab[y + 1][x], x, y);
        qsort(qall, (size_t)nq, sizeof(float), cmp_f);
        r->q50 = quant(qall, nq, 0.5f);
    }

    /* pass 1: measure every crossing, keep the pair's widths for its median */
    for (int axis = 0; axis < 2; axis++) {
        for (int j = 0; j < N; j++) {
            for (int i = 0; i < N - 1; i++) {
                int x = axis == 0 ? i : j;
                int y = axis == 0 ? j : i;
                int xn = axis == 0 ? i + 1 : x;
                int yn = axis == 0 ? y : i + 1;
                int a = lab[y][x], b = lab[yn][xn];
                if (a == b) continue;
                if (r->crossings >= MAXCROSS) { fprintf(stderr, "crossing buffer full\n"); exit(2); }
                r->crossings++;

                float dp = 0.0f;
                for (int k = 0; k < NB; k++) dp += fabsf(P[k][yn][xn] - P[k][y][x]);
                if (dp > r->dp_max) r->dp_max = dp;

                Cand c = { 0 };
                c.x = (unsigned short)x; c.y = (unsigned short)y;
                c.axis = (unsigned char)axis;
                c.a = (unsigned char)a; c.b = (unsigned char)b;
                c.dp = dp; c.closed = 2;
                c.w = 1.0f;                       /* hard map: a one-cell cliff */
                c.q = margin_grad(a, b, x, y);
                if (!r->hard) {
                    if (P[a][y][x] + P[b][y][x] < 0.5f) { r->pair_diluted++; continue; }
                    int cl;
                    float capped;
                    c.w = measure_width(a, b, x, y, r->tau, &cl, &capped);
                    c.closed = (unsigned char)cl;
                    if (c.closed == 3) {
                        r->slivers++;
                        int a2 = a < b ? a : b, b2 = a < b ? b : a;
                        r->pair_sliver[a2][b2]++;
                        continue;
                    }
                    if (c.closed != 2) { r->censored++; continue; }

                }
                cand[nc++] = c;
            }
        }
    }

    /* pass 2: statistics and classification.
     *   cut   - band still closes inside 1.5 cells: a hard edge in the map.
     *   step  - width collapsed to under half the family median while |grad m|
     *           jumps at least 2x along the boundary tangent: a crease.
     * Crossings that are censored or diluted are already excluded. */
    static float widths[MAXCROSS], dps[MAXCROSS], sv[MAXCROSS];
    static float famw[NB][NB][4096];
    static long famn[NB][NB];
    memset(famn, 0, sizeof(famn));
    int nw = 0;
    for (int k = 0; k < nc; k++) {
        Cand* c = &cand[k];
        widths[nw] = c->w;
        dps[nw] = c->dp;
        nw++;
        /* cut: the transition closes inside 1.5 cells regardless of W.
         * gradient step: the local slope predicts a wide ramp (pred >= 3
         * cells) but the actual transition is less than half of it - the
         * margin's slope changes inside the band, so |grad P| has a terrace. */
        int is_cut = c->w <= 1.5f;
        float pred = c->q > 1e-9f ? (float)(K_WIDTH * r->tau / c->q) : 1e6f;
        int is_step = !r->hard && !is_cut && pred >= 3.0f && c->w < 0.5f * pred;
        if (is_cut || is_step) {
            step_push(r, c->x, c->y, c->axis, c->a, c->b, c->w, c->dp,
                      is_step ? 1 : 0);
            if (!r->hard) {
                int a2 = c->a < c->b ? c->a : c->b, b2 = c->a < c->b ? c->b : c->a;
                pair_compressed_add(r, a2, b2);
            }
            if (is_step) r->slope_steps++;
        }
        if (!r->hard) {
            int a2 = c->a < c->b ? c->a : c->b, b2 = c->a < c->b ? c->b : c->a;
            r->pair_n[a2][b2]++;
            if (famn[a2][b2] < 4096) famw[a2][b2][famn[a2][b2]++] = c->w;
        }
    }

    if (!r->hard) {
        for (int a = 0; a < NB; a++) for (int b = a + 1; b < NB; b++) {
            if (!famn[a][b]) continue;
            qsort(famw[a][b], (size_t)famn[a][b], sizeof(float), cmp_f);
            r->pair_wmed[a][b] = quant(famw[a][b], (int)famn[a][b], 0.5f);
        }
    }
    if (nw <= 0) return;

    memcpy(sv, widths, (size_t)nw * sizeof(float));
    qsort(sv, (size_t)nw, sizeof(float), cmp_f);
    r->w_p05 = quant(sv, nw, 0.05f);
    r->w_p50 = quant(sv, nw, 0.50f);
    r->w_p95 = quant(sv, nw, 0.95f);
    double sum = 0.0;
    for (int i = 0; i < nw; i++) sum += widths[i];
    r->w_mean = (float)(sum / nw);
    double var = 0.0;
    for (int i = 0; i < nw; i++) { double d = widths[i] - r->w_mean; var += d * d; }
    r->w_cv = r->w_mean > 0.0f ? (float)(sqrt(var / nw) / r->w_mean) : 0.0f;

    memcpy(sv, dps, (size_t)nw * sizeof(float));
    qsort(sv, (size_t)nw, sizeof(float), cmp_f);
    r->dp50 = quant(sv, nw, 0.50f);
    r->dp95 = quant(sv, nw, 0.95f);

    if (!r->hard) {
        static float preds[MAXCROSS];
        for (int k = 0; k < nc; k++)
            preds[k] = cand[k].q > 1e-9f ? (float)(K_WIDTH * r->tau / cand[k].q) : 1e6f;
        static float pv[MAXCROSS], wv2[MAXCROSS];
        memcpy(pv, preds, (size_t)nc * sizeof(float));
        qsort(pv, (size_t)nc, sizeof(float), cmp_f);
        r->pred_w50 = quant(pv, nc, 0.5f);
        for (int k = 0; k < nc; k++) wv2[k] = cand[k].w;
        r->corr = pearson(preds, wv2, nc);
        long hit = 0;
        for (int k = 0; k < nc; k++) if (fabsf(cand[k].w - r->W) <= 0.5f * r->W) hit++;
        r->frac_hit = nc ? (float)hit / (float)nc : 0.0f;
        r->frac_closed = (float)nc / (float)(nc + r->censored);
    } else {
        r->frac_closed = 1.0f;
    }
}

/* config is read from r on entry; measurements are reset */
static void run_case(Run* r, int hard) {
    int lab[N][N];
    const char* tag = r->tag;
    int idx = r->run_index;
    float wt = r->wt, wh = r->wh, W = r->W, tau = r->tau;
    StepRec* step = r->step;
    long cap = r->step_cap;
    memset(r, 0, sizeof(*r));
    r->tag = tag; r->run_index = idx; r->wt = wt; r->wh = wh; r->W = W; r->hard = hard;
    r->step = step; r->step_cap = cap;

    compute_scores(wt, wh);
    hard_labels(lab);
    if (hard) { r->tau = 0.0f; onehot_from_labels(lab); }
    else      { r->tau = tau; softmax(tau); }

    gradient_field();
    boundary_gradient_stats(lab, r);
    area_counts(r, lab);
    collect(r, lab);
}

/* tau such that the median measured width equals W.  The anchor set is fixed
 * (the hard-map boundary crossings); the objective is the median of the
 * capped width, which counts a band that runs off the map as its lower bound
 * and is therefore monotone in tau.  Bisection on log tau. */
static float calibrate_tau(float wt, float wh, float W, int verbose) {
    int lab[N][N];
    compute_scores(wt, wh);
    hard_labels(lab);

    static int ax[4096], ay[4096], aa[4096], ab[4096];
    int na = 0;
    for (int y = 1; y < N - 1 && na < 4096; y += 3)
        for (int x = 1; x < N - 2 && na < 4096; x += 2) {
            if (lab[y][x] == lab[y][x + 1]) continue;
            ax[na] = x; ay[na] = y; aa[na] = lab[y][x]; ab[na] = lab[y][x + 1];
            na++;
        }

    static float qv[4096];
    for (int i = 0; i < na; i++) qv[i] = margin_grad(aa[i], ab[i], ax[i], ay[i]);
    qsort(qv, (size_t)na, sizeof(float), cmp_f);
    float q50 = quant(qv, na, 0.5f);

    float lo = 1e-5f, hi = 1.0f;
    float best = (float)(W * q50 / K_WIDTH);
    for (int it = 0; it < 30; it++) {
        float mid = sqrtf(lo * hi);
        static float wv[4096];
        int n = 0;
        for (int i = 0; i < na; i++) {
            int cl;
            float capped;
            measure_width(aa[i], ab[i], ax[i], ay[i], mid, &cl, &capped);
            wv[n++] = capped;
        }
        qsort(wv, (size_t)n, sizeof(float), cmp_f);
        float mw = quant(wv, n, 0.5f);
        if (verbose && it % 6 == 0)
            printf("    bisect tau=%.6f median capped width=%.2f (anchors=%d)\n", mid, mw, n);
        if (mw < W) lo = mid; else hi = mid;
        if (fabsf(mw - W) <= 0.02f * W) { best = mid; break; }
        best = mid;
    }
    if (verbose) printf("    anchor count=%d q50=%.5f -> tau=%.6f\n", na, q50, best);
    return best;
}

/* ---------- output ---------------------------------------------------- */

static FILE* open_file(const char* name) {
    char path[256];
    snprintf(path, sizeof(path), "reports/%s", name);
    FILE* f = fopen(path, "w");
    if (!f) { perror(path); exit(1); }
    return f;
}

static const unsigned char PAL[NB][3] = {
    { 210, 235, 255 },  /* ice */
    {  40, 100,  70 },  /* taiga */
    { 170, 200,  90 },  /* grassland */
    {  30, 130, 110 },  /* temperate_rain */
    { 230, 200, 120 },  /* desert */
    {  20, 140,  40 },  /* tropical_rain */
};

static void write_ppm(const Run* r) {
    static unsigned char img[N][N][3];
    for (int y = 0; y < N; y++) for (int x = 0; x < N; x++) {
        float acc[3] = { 0, 0, 0 };
        for (int i = 0; i < NB; i++) for (int k = 0; k < 3; k++)
            acc[k] += P[i][y][x] * PAL[i][k];
        for (int k = 0; k < 3; k++)
            img[y][x][k] = (unsigned char)clampf(acc[k], 0.0f, 255.0f);
    }
    char path[256];
    snprintf(path, sizeof(path), "reports/field_%s.ppm", r->tag);
    FILE* f = fopen(path, "wb");
    if (!f) { perror(path); return; }
    fprintf(f, "P6\n%d %d\n255\n", N, N);
    for (int y = 0; y < N; y++) fwrite(img[y], 1, (size_t)N * 3, f);
    fclose(f);
}

static void write_steps(const Run* r) {
    char name[96];
    snprintf(name, sizeof(name), "steps_%s.csv", r->tag);
    FILE* f = open_file(name);
    fprintf(f, "x,y,axis,biome_a,biome_b,width_cells,dp,kind\n");
    for (long i = 0; i < r->step_count; i++) {
        const StepRec* s = &r->step[i];
        fprintf(f, "%d,%d,%c,%s,%s,%.2f,%.4f,%s\n", s->x, s->y, s->axis ? 'v' : 'h',
                BIOMES[s->a].name, BIOMES[s->b].name, s->width, s->dp,
                s->kind ? "gstep" : "cut");
    }
    fclose(f);
}

int main(void) {
    float eq1_err = 0.0f;
    mkdir("reports", 0755);
    gen_fields();
    compute_gradients();
    printf("biome transition experiment: %dx%d grid, %d biomes\n", N, N, NB);

    /* ---- 1) weights: fixed point of w_t : w_h = B : A ---- */
    struct { float wt, wh; double A, B; long n; } it[24];
    int n_it = 0;
    it[0].wt = 0.5f; it[0].wh = 0.5f;
    boundary_terms(it[0].wt, it[0].wh, &it[0].A, &it[0].B, &it[0].n);
    for (;;) {
        double A = it[n_it].A, B = it[n_it].B;
        double wt_new = B / (A + B);
        double prev = it[n_it].wt;
        n_it++;
        it[n_it].wt = (float)wt_new;
        it[n_it].wh = (float)(1.0 - wt_new);
        boundary_terms(it[n_it].wt, it[n_it].wh, &it[n_it].A, &it[n_it].B, &it[n_it].n);
        if (fabs(wt_new - prev) < 1e-4 || n_it >= 23) break;
    }
    float wt = it[n_it].wt, wh = it[n_it].wh;

    printf("\n[weight derivation]  w_t : w_h = B : A  (equal RMS axis contributions)\n");
    printf("  it  w_t       w_h       A (temp)   B (humid)  boundary_cells\n");
    for (int i = 0; i <= n_it; i++)
        printf("  %2d  %.6f  %.6f  %.8f  %.8f  %8ld\n",
               i, it[i].wt, it[i].wh, it[i].A, it[i].B, it[i].n);
    printf("  fixed point after %d iterations: w_t=%.6f w_h=%.6f  (|dw|=%.2e)\n",
           n_it, wt, wh, fabs(it[n_it].wt - it[n_it - 1].wt));
    printf("  width relation (1): W = 2 ln9 * tau / |grad m|,\n"
           "                      |grad m| = 2 w_t dT gradT + 2 w_h dH gradH\n");

    /* ---- 2) tau per target width ---- */
    printf("\n[tau calibration]\n");
    float targets[3] = { 3.0f, 10.0f, 24.0f };
    float taus[3];
    for (int i = 0; i < 3; i++) {
        printf("  W=%.0f\n", targets[i]);
        taus[i] = calibrate_tau(wt, wh, targets[i], 1);
        printf("  -> tau=%.6f\n", taus[i]);
    }

    /* ---- 3) runs ---- */
    Run runs[6];
    memset(runs, 0, sizeof(runs));
    runs[0] = (Run){ .tag = "hard_baseline", .run_index = 0, .wt = 0.5f, .wh = 0.5f, .hard = 1 };
    runs[1] = (Run){ .tag = "hard_derived",  .run_index = 1, .wt = wt, .wh = wh, .hard = 1 };
    runs[2] = (Run){ .tag = "w3",  .run_index = 2, .wt = wt, .wh = wh, .W = 3.0f,  .tau = taus[0] };
    runs[3] = (Run){ .tag = "w10", .run_index = 3, .wt = wt, .wh = wh, .W = 10.0f, .tau = taus[1] };
    runs[4] = (Run){ .tag = "w24", .run_index = 4, .wt = wt, .wh = wh, .W = 24.0f, .tau = taus[2] };
    runs[5] = (Run){ .tag = "eq_w10", .run_index = 5, .wt = 0.5f, .wh = 0.5f, .W = 10.0f,
                     .tau = calibrate_tau(0.5f, 0.5f, 10.0f, 0) };

    /* median |grad m| per boundary family on the derived-weight hard map */
    {
        int lab[N][N];
        static float qv[NB][NB][N * N / 4];
        static long qn[NB][NB];
        compute_scores(wt, wh);
        hard_labels(lab);
        memset(qn, 0, sizeof(qn));
        for (int y = 0; y < N; y++) for (int x = 0; x < N - 1; x++) {
            int a = lab[y][x], b = lab[y][x + 1];
            if (a == b) continue;
            int a2 = a < b ? a : b, b2 = a < b ? b : a;
            if (qn[a2][b2] < N * N / 4)
                qv[a2][b2][qn[a2][b2]++] = margin_grad(a, b, x, y);
        }
        for (int a = 0; a < NB; a++) for (int b = a + 1; b < NB; b++) {
            if (!qn[a][b]) continue;
            qsort(qv[a][b], (size_t)qn[a][b], sizeof(float), cmp_f);
            g_pair_q50[a][b] = quant(qv[a][b], (int)qn[a][b], 0.5f);
        }
    }
    for (int i = 0; i < 6; i++) {
        run_case(&runs[i], runs[i].hard);
        write_ppm(&runs[i]);
        write_steps(&runs[i]);
    }

    printf("\n[boundary transition width / gradient]\n");
    printf("  %-14s %4s %9s %7s %6s %6s %6s %6s %7s %7s %8s %8s %6s %6s %6s\n",
           "run", "W", "tau", "cross", "w05", "w50", "w95", "cv",
           "dp50", "dp95", "|gP|mean", "|gP|p95", "cuts", "gstep", "sliver");
    for (int i = 0; i < 6; i++) {
        const Run* r = &runs[i];
        printf("  %-14s %4.1f %9.6f %7ld %6.2f %6.2f %6.2f %6.2f %7.3f %7.3f %8.4f %8.4f %6ld %6ld %6ld\n",
               r->tag, r->W, r->tau, r->crossings,
               r->w_p05, r->w_p50, r->w_p95, r->w_cv,
               r->dp50, r->dp95, r->gb_mean, r->gb_p95,
               r->step_count - r->slope_steps, r->slope_steps, r->slivers);
    }
    printf("  w05/w50/w95  perpendicular width of the a/b margin band |m| <= ln9*tau, in cells\n"
           "  cv           coefficient of variation of that width over closed crossings\n"
           "  dp           adjacent-cell L1 probability jump\n"
           "  |gP|         |grad P| over all biome probabilities, boundary cells only\n"
           "  cuts         boundary cells whose transition closes inside 1.5 cells;\n"
           "  gstep        boundary cells where the local slope predicts a ramp >= 3 cells but\n"
           "               the measured band is under half of that: |grad P| terraces there;\n"
           "  sliver       boundary pairs with |grad m| ~ 0 (linear anchor > 6 cells away), excluded\n");
    printf("  worst |gP| cell per run: ");
    for (int i = 0; i < 6; i++)
        printf("%s=(%d,%d,%.3f) ", runs[i].tag, runs[i].gx_worst, runs[i].gy_worst,
               runs[i].gb_max);
    printf("\n");
    printf("  censored (band ran off the map): ");
    for (int i = 0; i < 6; i++) printf("%s=%ld ", runs[i].tag, runs[i].censored);
    printf("\n");

    printf("\n[measured vs eq.(1); share of crossings within 50%% of target]\n");
    printf("  %-14s %6s %9s %9s %9s %9s %7s %7s\n",
           "run", "W", "pred50", "w50", "within50", "closed", "r", "diluted");
    for (int i = 2; i <= 5; i++) {
        const Run* r = &runs[i];
        printf("  %-14s %6.1f %9.2f %9.2f %8.1f%% %8.1f%% %7.2f %7ld\n",
               r->tag, r->W, r->pred_w50, r->w_p50, 100.0f * r->frac_hit,
               100.0f * r->frac_closed, r->corr, r->pair_diluted);
    }
    printf("  pred50 = median over crossings of 2 ln9 tau / |grad m|, with |grad m| the\n"
           "  crossing cell's own central difference (biased high: the band averages a\n"
           "  gentler slope than the single cell at the crossing, and censoring removes the\n"
           "  widest survivors). It brackets the trend, not the absolute width; tau is\n"
           "  calibrated on the measured median, not on this number.\n"
           "  pred50 = median over crossings of 2 ln9 tau / |grad m|;\n"
           "  closed = crossings whose 10-90 band closed inside the map;\n"
           "  r = Pearson(pred, measured) over those crossings;\n"
           "  diluted = crossings dropped because P[a]+P[b] < 0.5 (a third biome rules).\n");

    printf("\n[median width per boundary family, cells]\n");
    printf("  %-30s %6s %8s %7s %7s %7s %8s %6s %6s %6s\n",
           "pair (dT,dH)", "n", "q50", "w3", "w10", "w24", "w24/w3", "pred24",
           "cut", "sliver");
    for (int a = 0; a < NB; a++) for (int b = a + 1; b < NB; b++) {
        long n = runs[2].pair_n[a][b];
        long n24 = runs[4].pair_n[a][b];
        if (!n && !n24) continue;
        char label[64];
        snprintf(label, sizeof(label), "%s/%s (%.2f,%.2f)", BIOMES[a].name, BIOMES[b].name,
                 fabsf(BIOMES[a].t - BIOMES[b].t), fabsf(BIOMES[a].h - BIOMES[b].h));
        float q = g_pair_q50[a][b];
        float m3 = (float)runs[2].pair_wmed[a][b];
        float m10 = (float)runs[3].pair_wmed[a][b];
        float m24 = (float)runs[4].pair_wmed[a][b];
        float pred24 = q > 1e-9f ? (float)(K_WIDTH * runs[4].tau / q) : 0.0f;
        printf("  %-30s %6ld %8.5f %7.2f %7.2f %7.2f %8.2f %6.1f %6ld %6ld\n",
               label, n, q, m3, m10, m24, m3 > 0.0f ? m24 / m3 : 0.0f, pred24,
               runs[4].pair_compressed[a][b], runs[4].pair_sliver[a][b]);
    }
    printf("  q50 = median |grad m|; pred24 = 2 ln9 tau24 / q50; cut/sliver counted at W=24.\n"
           "  Every family's median should scale with the requested W (w24/w3 ~ 5.7 for\n"
           "  the 24/3 target ratio); a family that stays flat is one the knob cannot reach.\n");

    printf("\n[biome area: hard cells / soft mass]\n");
    printf("  %-16s %10s %10s", "biome", "hard_der", "hard_base");
    for (int i = 2; i <= 4; i++) printf(" %10s", runs[i].tag);
    printf(" %10s\n", runs[5].tag);
    for (int b = 0; b < NB; b++) {
        printf("  %-16s %10ld %10ld", BIOMES[b].name, runs[1].hard_area[b], runs[0].hard_area[b]);
        for (int i = 2; i <= 5; i++) printf(" %10.1f", runs[i].soft_area[b]);
        printf("\n");
    }
    printf("  %-16s %10ld %10ld", "total", (long)N * N, (long)N * N);
    for (int i = 2; i <= 5; i++) {
        double s = 0.0;
        for (int b = 0; b < NB; b++) s += runs[i].soft_area[b];
        printf(" %10.1f", s);
    }
    printf("\n  hard = argmax cells; soft = sum of P_i over the map (expected count).\n");

    /* eq. (1) is exact for a locally linear margin; check it on synthetic ones */
    {
        float save_S[NB][N][N];
        memcpy(save_S, S, sizeof(S));
        int a = 0, b = 1;
        for (int y = 0; y < N; y++) for (int xx = 0; xx < N; xx++) {
            S[a][y][xx] = 0.01f * ((float)(xx - 128) * 0.866f + (float)(y - 128) * 0.5f);
            S[b][y][xx] = 0.0f;
        }
        float tau = 0.02f, expected = (float)(K_WIDTH * tau / 0.01f);
        int cl;
        float capped;
        float got = measure_width(a, b, 128, 128, tau, &cl, &capped);
        printf("\n[eq.(1) on a synthetic 30-degree linear margin]\n");
        printf("  tau=%.4f |grad m|=0.0100 expected width=%.4f measured=%.4f (closed=%d)\n",
               tau, expected, got, cl);
        memcpy(S, save_S, sizeof(S));
        eq1_err = fabsf(got - expected) / expected;
    }

    printf("\n[staircase coordinates: cuts (width <= 1.5 cells) and gradient steps]\n");
    for (int i = 0; i < 6; i++) {
        const Run* r = &runs[i];
        printf("  %-14s listed=%ld (cuts=%ld, gstep=%ld) / %ld crossings\n",
               r->tag, r->step_count, r->step_count - r->slope_steps,
               r->slope_steps, r->crossings);
        int shown = 0, maxshow = (i <= 1) ? 4 : 8;
        for (long k = 0; k < r->step_count && shown < maxshow; k++) {
            const StepRec* s = &r->step[k];
            printf("      (%3d,%3d) %c %-5s %s/%s width=%.2f dp=%.3f\n",
                   s->x, s->y, s->axis ? 'v' : 'h', s->kind ? "gstep" : "cut",
                   BIOMES[s->a].name, BIOMES[s->b].name, s->width, s->dp);
            shown++;
        }
        if (r->step_count > maxshow)
            printf("      ... %ld more in reports/steps_%s.csv\n",
                   r->step_count - maxshow, r->tag);
    }

    /* ---- CSVs ---- */
    FILE* f = open_file("weights.csv");
    fprintf(f, "quantity,value\n");
    for (int i = 0; i <= n_it; i++)
        fprintf(f, "iter%d_wt,%.8f\niter%d_wh,%.8f\niter%d_A,%.10f\niter%d_B,%.10f\niter%d_n,%ld\n",
                i, it[i].wt, i, it[i].wh, i, it[i].A, i, it[i].B, i, it[i].n);
    fprintf(f, "final_wt,%.8f\nfinal_wh,%.8f\n", wt, wh);
    fprintf(f, "width_relation,\"W = 2*ln(9)*tau/|grad m|\"\n");
    fprintf(f, "weight_relation,\"w_t:w_h = B:A\"\n");
    fclose(f);

    f = open_file("runs.csv");
    fprintf(f, "tag,hard,wt,wh,W,tau,crossings,censored,listed,width_p05,width_p50,width_p95,"
               "width_mean,width_cv,q50,pred_width_50,frac_within_50pct,dp_median,dp_p95,dp_max,"
               "gb_mean,gb_p95,gb_max,gb_theory,frac_closed,corr_pred_width,pair_diluted,"
               "slivers,gradient_steps,cuts,boundary_cells,band_cells_lt09\n");
    for (int i = 0; i < 6; i++) {
        const Run* r = &runs[i];
        fprintf(f, "%s,%d,%.6f,%.6f,%.2f,%.8f,%ld,%ld,%ld,%.3f,%.3f,%.3f,%.3f,%.4f,"
                   "%.6f,%.3f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,"
                   "%ld,%ld,%ld,%ld,%ld,%ld\n",
                r->tag, r->hard, r->wt, r->wh, r->W, r->tau, r->crossings, r->censored,
                r->step_count, r->w_p05, r->w_p50, r->w_p95, r->w_mean, r->w_cv,
                r->q50, r->hard ? 0.0f : r->pred_w50, r->hard ? 0.0f : r->frac_hit,
                r->dp50, r->dp95, r->dp_max, r->gb_mean, r->gb_p95, r->gb_max,
                r->g_theory, r->frac_closed, r->corr, r->pair_diluted,
                r->slivers, r->slope_steps, r->step_count - r->slope_steps,
                r->boundary_cells, r->band_cells);
    }
    fclose(f);

    f = open_file("area.csv");
    fprintf(f, "run,biome,hard_cells,soft_cells,soft_frac\n");
    for (int i = 0; i < 6; i++) for (int b = 0; b < NB; b++)
        fprintf(f, "%s,%s,%ld,%.3f,%.6f\n", runs[i].tag, BIOMES[b].name,
                runs[i].hard_area[b], runs[i].soft_area[b],
                runs[i].soft_area[b] / (double)(N * N));
    fclose(f);

    /* ---- self checks ---- */
    int failures = 0;
#define CHECK(name, ok) do { printf("%s=%d\n", name, (ok) ? 1 : 0); if (!(ok)) failures++; } while (0)
    int in_range = 1;
    for (int y = 0; y < N; y += 7) for (int x = 0; x < N; x += 7)
        if (T[y][x] < 0.0f || T[y][x] > 1.0f || H[y][x] < 0.0f || H[y][x] > 1.0f) in_range = 0;
    CHECK("fields_normalized", in_range);
    CHECK("weights_sum_one", fabs(wt + wh - 1.0f) < 1e-6f);
    CHECK("weights_positive", wt > 0.0f && wh > 0.0f);
    CHECK("weight_fixed_point", fabs(it[n_it].wt - it[n_it - 1].wt) < 1e-3f);
    CHECK("faster_axis_downweighted", wh < wt);

    compute_scores(wt, wh);
    softmax(runs[3].tau);
    int argmax_ok = 1, prob_ok = 1;
    for (int y = 0; y < N; y += 11) for (int x = 0; x < N; x += 11) {
        int a = label_at(x, y), best = 0;
        for (int i = 1; i < NB; i++) if (P[i][y][x] > P[best][y][x]) best = i;
        if (a != best) argmax_ok = 0;
        float s = 0.0f;
        for (int i = 0; i < NB; i++) s += P[i][y][x];
        if (fabsf(s - 1.0f) > 1e-5f) prob_ok = 0;
    }
    CHECK("softmax_keeps_argmax", argmax_ok);
    CHECK("softmax_sums_one", prob_ok);
    CHECK("tau_increases_with_width",
          runs[2].tau < runs[3].tau && runs[3].tau < runs[4].tau);
    CHECK("width_monotonic", runs[2].w_p50 < runs[3].w_p50 && runs[3].w_p50 < runs[4].w_p50);
    CHECK("width_hits_target", fabsf(runs[2].w_p50 - 3.0f) < 1.0f &&
                               fabsf(runs[3].w_p50 - 10.0f) < 2.0f &&
                               fabsf(runs[4].w_p50 - 24.0f) < 5.0f);
    CHECK("cuts_vanish_with_width", runs[4].step_count < runs[2].step_count);
    CHECK("eq1_ratio_stable",
          runs[3].w_p50 / runs[3].pred_w50 > 1.0f && runs[3].w_p50 / runs[3].pred_w50 < 3.0f);
    CHECK("hard_has_steps_everywhere", runs[0].step_count == runs[0].crossings);
    CHECK("soft_band_grows", runs[2].band_cells < runs[3].band_cells &&
                             runs[3].band_cells < runs[4].band_cells);
    CHECK("gradient_falls_with_width", runs[4].gb_p95 < runs[2].gb_p95);
    CHECK("eq1_predicts_median", fabsf(runs[3].pred_w50 - runs[3].w_p50) < 0.5f * runs[3].w_p50);
    CHECK("most_bands_close", runs[3].frac_closed > 0.8f);
    CHECK("eq1_exact_on_linear_margin", eq1_err < 1e-3f);
    CHECK("gradient_near_ramp_theory", runs[3].gb_p95 <= 4.0f * runs[3].g_theory);
    long tot = 0;
    for (int b = 0; b < NB; b++) tot += runs[1].hard_area[b];
    CHECK("hard_area_total", tot == (long)N * N);
    double stot = 0.0;
    for (int b = 0; b < NB; b++) stot += runs[4].soft_area[b];
    CHECK("soft_area_total", fabs(stot - (double)N * N) < 1e-2);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
