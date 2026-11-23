
/* warp: domain warping of a heightfield, amplitude capped by derived bounds.
 *
 * Warp map: w(p) = p + A * q(p), q = (g1, g2), g_i fBm in [-1, 1].
 *
 * 1) amplitude vs distortion: displacement d(p) = A q(p), |d| <= sqrt(2) A;
 *    strain eps(A) = A * Lq with Lq = max ||J_q||_2; folds (non-injective w,
 *    i.e. mush) become possible only when A * Lq >= 1.
 * 2) gradient bound: grad F = J_w^T grad H  =>  |grad F| <= G_H (1 + A Lq);
 *    a grid with cell h resolves gradients up to G_allow = 1/(2h).
 * 3) cap: A_max = min(EPS_MAX/Lq, D_MAX*lambda_min/sqrt(2),
 *                     (G_allow/G_H - 1)/Lq), all verified numerically below.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define N 256
#define DOM 16.0f
#define CELL (DOM / (float)N)

/* terrain: ramp + fBm detail, so the coastline starts nearly straight */
#define T_FREQ 0.2f
#define T_OCT 2
#define T_AMP 0.04f
#define RAMP 0.055f
#define T_SEED 11u

/* warp vector field */
#define W_FREQ 0.2f
#define W_OCT 2
#define W_SEED1 101u
#define W_SEED2 202u

#define EPS_MAX 0.5f   /* max allowed strain (relative feature distortion) */
#define D_MAX 0.5f     /* max displacement, in units of lambda_min */
#define SEA 0.5f       /* coastline level */

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

static float H[N][N], Q1[N][N], Q2[N][N], F[N][N];

static void build_terrain(void) {
    for (int y = 0; y < N; y++)
        for (int x = 0; x < N; x++) {
            float wx = x * CELL, wy = y * CELL;
            H[y][x] = 0.5f + RAMP * (wy - DOM * 0.5f)
                    + T_AMP * fbm(wx * T_FREQ, wy * T_FREQ, T_OCT, 2.0f, 0.5f, T_SEED);
        }
}

static void build_warp_field(void) {
    for (int y = 0; y < N; y++)
        for (int x = 0; x < N; x++) {
            float wx = x * CELL, wy = y * CELL;
            Q1[y][x] = fbm(wx * W_FREQ, wy * W_FREQ, W_OCT, 2.0f, 0.5f, W_SEED1);
            Q2[y][x] = fbm(wx * W_FREQ, wy * W_FREQ, W_OCT, 2.0f, 0.5f, W_SEED2);
        }
}

static float sample_h(float wx, float wy) {
    float gx = wx / CELL, gy = wy / CELL;
    if (gx < 0) gx = 0;
    if (gy < 0) gy = 0;
    if (gx > N - 1) gx = N - 1;
    if (gy > N - 1) gy = N - 1;
    int x0 = (int)floorf(gx), y0 = (int)floorf(gy);
    float tx = gx - x0, ty = gy - y0;
    int x1 = x0 + 1, y1 = y0 + 1;
    if (x1 >= N) x1 = N - 1;
    if (y1 >= N) y1 = N - 1;
    float a = H[y0][x0] * (1 - tx) + H[y0][x1] * tx;
    float b = H[y1][x0] * (1 - tx) + H[y1][x1] * tx;
    return a * (1 - ty) + b * ty;
}

static void build_warped(float A) {
    for (int y = 0; y < N; y++)
        for (int x = 0; x < N; x++)
            F[y][x] = sample_h(x * CELL + A * Q1[y][x], y * CELL + A * Q2[y][x]);
}

/* max singular value of 2x2 [[a,b],[c,d]] */
static float svmax2(float a, float b, float c, float d) {
    float t = a * a + b * b + c * c + d * d;
    float det = a * d - b * c;
    float disc = t * t - 4.0f * det * det;
    return sqrtf(0.5f * (t + sqrtf(disc > 0.0f ? disc : 0.0f)));
}

/* Lq = max ||J_q||_2 over the grid (central differences) */
static float measure_lq(void) {
    float m = 0.0f;
    for (int y = 1; y < N - 1; y++)
        for (int x = 1; x < N - 1; x++) {
            float a = (Q1[y][x + 1] - Q1[y][x - 1]) / (2 * CELL);
            float b = (Q1[y + 1][x] - Q1[y - 1][x]) / (2 * CELL);
            float c = (Q2[y][x + 1] - Q2[y][x - 1]) / (2 * CELL);
            float d = (Q2[y + 1][x] - Q2[y - 1][x]) / (2 * CELL);
            float s = svmax2(a, b, c, d);
            if (s > m) m = s;
        }
    return m;
}

/* max strain of J_w = I + A J_q, minus 1 */
static float measure_strain(float A) {
    float m = 0.0f;
    for (int y = 1; y < N - 1; y++)
        for (int x = 1; x < N - 1; x++) {
            float a = 1.0f + A * (Q1[y][x + 1] - Q1[y][x - 1]) / (2 * CELL);
            float b = A * (Q1[y + 1][x] - Q1[y - 1][x]) / (2 * CELL);
            float c = A * (Q2[y][x + 1] - Q2[y][x - 1]) / (2 * CELL);
            float d = 1.0f + A * (Q2[y + 1][x] - Q2[y - 1][x]) / (2 * CELL);
            float s = svmax2(a, b, c, d) - 1.0f;
            if (s > m) m = s;
        }
    return m;
}

static float grad_max(float X[N][N]) {
    float m = 0.0f;
    for (int y = 1; y < N - 1; y++)
        for (int x = 1; x < N - 1; x++) {
            float gx = (X[y][x + 1] - X[y][x - 1]) / (2 * CELL);
            float gy = (X[y + 1][x] - X[y - 1][x]) / (2 * CELL);
            float g = sqrtf(gx * gx + gy * gy);
            if (g > m) m = g;
        }
    return m;
}

/* max |grad| of the bilinear interpolant of X (per-cell edge differences) */
static float grad_max_interp(float X[N][N]) {
    float m = 0.0f;
    for (int y = 0; y < N - 1; y++)
        for (int x = 0; x < N - 1; x++) {
            float gx = fmaxf(fabsf(X[y][x + 1] - X[y][x]),
                             fabsf(X[y + 1][x + 1] - X[y + 1][x])) / CELL;
            float gy = fmaxf(fabsf(X[y + 1][x] - X[y][x]),
                             fabsf(X[y + 1][x + 1] - X[y][x + 1])) / CELL;
            float g = sqrtf(gx * gx + gy * gy);
            if (g > m) m = g;
        }
    return m;
}

/* ---- coastline extraction: marching squares at SEA level ---- */

#define NSEGMAX (2 * (N - 1) * (N - 1))
#define NEDGE (2 * N * (N - 1))

enum { EB, ER, ET, EL };

static float sx1[NSEGMAX], sy1[NSEGMAX], sx2[NSEGMAX], sy2[NSEGMAX];
static int se1[NSEGMAX], se2[NSEGMAX], nseg;
static int e2s[NEDGE][2];

static int edge_id(int kind, int x, int y) {
    switch (kind) {
    case EB: return y * (N - 1) + x;                 /* (x,y)-(x+1,y) */
    case ET: return (y + 1) * (N - 1) + x;
    case EL: return (N - 1) * N + y * N + x;         /* (x,y)-(x,y+1) */
    default:  return (N - 1) * N + y * N + (x + 1);
    }
}

static void edge_point(int kind, int x, int y, float* px, float* py) {
    float v0, v1, t;
    if (kind == EB)      { v0 = F[y][x];     v1 = F[y][x + 1]; }
    else if (kind == ET) { v0 = F[y + 1][x]; v1 = F[y + 1][x + 1]; }
    else if (kind == EL) { v0 = F[y][x];     v1 = F[y + 1][x]; }
    else                 { v0 = F[y][x + 1]; v1 = F[y + 1][x + 1]; }
    t = (SEA - v0) / (v1 - v0);
    if (kind == EB)      { *px = (x + t) * CELL; *py = y * CELL; }
    else if (kind == ET) { *px = (x + t) * CELL; *py = (y + 1) * CELL; }
    else if (kind == EL) { *px = x * CELL;       *py = (y + t) * CELL; }
    else                 { *px = (x + 1) * CELL; *py = (y + t) * CELL; }
}

static void add_segment(int e1, float x1, float y1, int e2, float x2, float y2) {
    int s = nseg++;
    sx1[s] = x1; sy1[s] = y1; sx2[s] = x2; sy2[s] = y2;
    se1[s] = e1; se2[s] = e2;
    if (e2s[e1][0] < 0) e2s[e1][0] = s; else e2s[e1][1] = s;
    if (e2s[e2][0] < 0) e2s[e2][0] = s; else e2s[e2][1] = s;
}

/* segment endpoint pairs per marching-squares case (bit0=v00,1=v10,2=v11,3=v01) */
static const signed char seg_tab[16][4] = {
    {-1, -1, -1, -1}, {EL, EB, -1, -1}, {EB, ER, -1, -1}, {EL, ER, -1, -1},
    {ER, ET, -1, -1}, {EB, EL, ER, ET}, {EB, ET, -1, -1}, {EL, ET, -1, -1},
    {ET, EL, -1, -1}, {EB, ET, -1, -1}, {EB, ER, EL, ET}, {ER, ET, -1, -1},
    {EL, ER, -1, -1}, {EB, ER, -1, -1}, {EB, EL, -1, -1}, {-1, -1, -1, -1}
};

static void build_segments(void) {
    nseg = 0;
    for (int i = 0; i < NEDGE; i++) e2s[i][0] = e2s[i][1] = -1;
    for (int y = 0; y < N - 1; y++)
        for (int x = 0; x < N - 1; x++) {
            int c = (F[y][x] > SEA)
                  | ((F[y][x + 1] > SEA) << 1)
                  | ((F[y + 1][x + 1] > SEA) << 2)
                  | ((F[y + 1][x] > SEA) << 3);
            for (int k = 0; k < 4 && seg_tab[c][k] >= 0; k += 2) {
                int k1 = seg_tab[c][k], k2 = seg_tab[c][k + 1];
                float x1, y1, x2, y2;
                edge_point(k1, x, y, &x1, &y1);
                edge_point(k2, x, y, &x2, &y2);
                add_segment(edge_id(k1, x, y), x1, y1, edge_id(k2, x, y), x2, y2);
            }
        }
}

/* chain segments; returns longest polyline, *npts points (malloc'd) */
static void longest_polyline(float** outx, float** outy, int* npts) {
    static char vis[NSEGMAX];
    float* bx = malloc(sizeof(float) * (size_t)(nseg + 2));
    float* by = malloc(sizeof(float) * (size_t)(nseg + 2));
    float* tx = malloc(sizeof(float) * (size_t)(2 * nseg + 4));
    float* ty = malloc(sizeof(float) * (size_t)(2 * nseg + 4));
    int best = 0;
    memset(vis, 0, (size_t)nseg);
    for (int s0 = 0; s0 < nseg; s0++) {
        if (vis[s0]) continue;
        int lo = nseg, hi = nseg + 1;   /* polyline occupies tx[lo..hi) */
        vis[s0] = 1;
        tx[lo] = sx1[s0]; ty[lo] = sy1[s0];
        tx[hi] = sx2[s0]; ty[hi] = sy2[s0];
        /* extend head via se1, tail via se2 */
        for (int dir = 0; dir < 2; dir++) {
            int cur = s0, e = dir == 0 ? se1[s0] : se2[s0];
            for (;;) {
                int s2 = (e2s[e][0] == cur) ? e2s[e][1] : e2s[e][0];
                if (s2 < 0 || vis[s2]) break;
                vis[s2] = 1;
                float px, py;
                if (se1[s2] == e) { px = sx2[s2]; py = sy2[s2]; e = se2[s2]; }
                else              { px = sx1[s2]; py = sy1[s2]; e = se1[s2]; }
                if (dir == 0) { tx[--lo] = px; ty[lo] = py; }
                else          { tx[++hi] = px; ty[hi] = py; }
                cur = s2;
            }
        }
        int n = hi - lo + 1;
        if (n > best) {
            best = n;
            memcpy(bx, tx + lo, sizeof(float) * (size_t)n);
            memcpy(by, ty + lo, sizeof(float) * (size_t)n);
        }
    }
    free(tx); free(ty);
    *outx = bx; *outy = by; *npts = best;
}

/* resample to uniform spacing, then curvature via circumcircle of triples */
static float poly_metrics(float* px, float* py, int n, float ds,
                          float* out_len, float* out_bend, float* out_cmax) {
    float* ux = malloc(sizeof(float) * (size_t)n);
    float* uy = malloc(sizeof(float) * (size_t)n);
    float total = 0.0f;
    for (int i = 1; i < n; i++)
        total += hypotf(px[i] - px[i - 1], py[i] - py[i - 1]);
    int m = 0;
    float target = 0.0f, acc = 0.0f;
    ux[m] = px[0]; uy[m] = py[0]; m++;
    for (int i = 1; i < n; i++) {
        float seg = hypotf(px[i] - px[i - 1], py[i] - py[i - 1]);
        while (acc + seg >= target + ds && target + ds <= total) {
            target += ds;
            float t = (target - acc) / seg;
            ux[m] = px[i - 1] + (px[i] - px[i - 1]) * t;
            uy[m] = py[i - 1] + (py[i] - py[i - 1]) * t;
            m++;
        }
        acc += seg;
    }
    int k = 3;
    float sum = 0.0f, cmax = 0.0f;
    int cnt = 0;
    for (int i = k; i < m - k; i++) {
        float d1x = ux[i] - ux[i - k], d1y = uy[i] - uy[i - k];
        float d2x = ux[i + k] - ux[i], d2y = uy[i + k] - uy[i];
        float cross = d1x * d2y - d1y * d2x;
        float l1 = hypotf(d1x, d1y), l2 = hypotf(d2x, d2y);
        float l3 = hypotf(ux[i + k] - ux[i - k], uy[i + k] - uy[i - k]);
        if (l1 < 1e-9f || l2 < 1e-9f || l3 < 1e-9f) continue;
        float curv = fabsf(2.0f * cross / (l1 * l2 * l3));
        sum += curv; cnt++;
        if (curv > cmax) cmax = curv;
    }
    float chord = hypotf(px[n - 1] - px[0], py[n - 1] - py[0]);
    *out_len = total;
    *out_bend = chord > 1e-3f ? total / chord : -1.0f;
    *out_cmax = cmax;
    free(ux); free(uy);
    return cnt > 0 ? sum / cnt : 0.0f;
}

static int write_bmp(const char* path) {
    FILE* f = fopen(path, "wb");
    if (!f) return 0;
    int row = (N * 3 + 3) & ~3;
    unsigned char hdr[54] = {0};
    unsigned long imgsz = (unsigned long)row * N, filesz = 54 + imgsz;
    unsigned long off = 54, hs = 40;
    int w = N, h = N;
    unsigned short planes = 1, bpp = 24;
    hdr[0] = 'B'; hdr[1] = 'M';
    memcpy(hdr + 2, &filesz, 4); memcpy(hdr + 10, &off, 4);
    memcpy(hdr + 14, &hs, 4); memcpy(hdr + 18, &w, 4); memcpy(hdr + 22, &h, 4);
    memcpy(hdr + 26, &planes, 2); memcpy(hdr + 28, &bpp, 2);
    fwrite(hdr, 1, 54, f);
    unsigned char* buf = calloc((size_t)row, 1);
    for (int y = N - 1; y >= 0; y--) {
        for (int x = 0; x < N; x++) {
            int land = F[y][x] > SEA;
            int edge = (x + 1 < N && (F[y][x + 1] > SEA) != land)
                    || (y + 1 < N && (F[y + 1][x] > SEA) != land);
            unsigned char r, g, b;
            if (edge) { r = 255; g = 0; b = 0; }
            else {
                float v = F[y][x];
                if (v < 0) v = 0;
                if (v > 1) v = 1;
                unsigned char gray = (unsigned char)(v * 255.0f);
                if (land) { r = gray; g = gray; b = gray; }
                else      { r = gray / 3; g = gray / 3; b = gray / 2; }
            }
            buf[x * 3 + 0] = b; buf[x * 3 + 1] = g; buf[x * 3 + 2] = r;
        }
        fwrite(buf, 1, (size_t)row, f);
    }
    free(buf);
    fclose(f);
    return 1;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(int argc, char** argv) {
    float frac = argc > 1 ? (float)atof(argv[1]) : 0.0f;
    const char* out = argc > 2 ? argv[2] : "warp.bmp";

    build_terrain();
    build_warp_field();

    float lq = measure_lq();
    /* analytic: per-octave |dn/dx| <= 1.5*2*f = 3f, so per component
     * L_g <= 3*sqrt(2)*O*f; vector field: Lq <= sqrt(2)*L_g = 6*O*f */
    float lq_analytic = 6.0f * W_OCT * W_FREQ;
    float gh = grad_max_interp(H);
    float lambda_min = 1.0f / (T_FREQ * (float)(1 << (T_OCT - 1)));
    float g_allow = 1.0f / (2.0f * CELL);

    float cap_strain = EPS_MAX / lq;
    float cap_feat = D_MAX * lambda_min / sqrtf(2.0f);
    float cap_grad = (g_allow / gh - 1.0f) / lq;
    float a_max = fminf(cap_strain, fminf(cap_feat, cap_grad));
    float A = frac * a_max;

    build_warped(A);
    float strain = measure_strain(A);
    float gf = grad_max(F);
    float grad_bound = gh * (1.0f + A * lq);
    float feat_d = sqrtf(2.0f) * A / lambda_min;

    build_segments();
    float* px = NULL; float* py = NULL; int npts = 0;
    longest_polyline(&px, &py, &npts);
    float clen = 0, bend = 0, cmax = 0;
    float cmean = npts > 8 ? poly_metrics(px, py, npts, 2.0f * CELL, &clen, &bend, &cmax) : 0.0f;

    printf("A=%.4f A_max=%.4f frac=%.3f\n", A, a_max, frac);
    printf("lq_meas=%.4f lq_analytic=%.4f\n", lq, lq_analytic);
    printf("cap_strain=%.4f cap_feature=%.4f cap_gradient=%.4f\n", cap_strain, cap_feat, cap_grad);
    printf("strain=%.4f eps_max=%.2f\n", strain, EPS_MAX);
    printf("grad_meas=%.4f grad_bound=%.4f g_allow=%.3f\n", gf, grad_bound, g_allow);
    printf("feature_D=%.4f d_max=%.2f\n", feat_d, D_MAX);
    printf("contour_pts=%d contour_len=%.3f bendiness=%.4f curv_mean=%.5f curv_max=%.5f\n",
           npts, clen, bend, cmean, cmax);

    check("lq_within_analytic", lq <= lq_analytic + 1e-4f);
    check("strain_ok", strain <= EPS_MAX + 1e-4f);
    check("grad_ok", gf <= grad_bound * 1.05f + 1e-6f);
    check("feature_ok", feat_d <= D_MAX + 1e-6f);
    check("amplitude_within_cap", A <= a_max + 1e-9f);
    check("bmp_written", write_bmp(out));
    printf("all_ok=%d\n", failures == 0);
    free(px); free(py);
    return failures == 0 ? 0 : 1;
}
