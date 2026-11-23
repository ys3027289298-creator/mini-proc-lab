/* stitch: seamless stitching of terrain tiles generated with different seeds.
 *
 * Boundary conditions for an invisible seam (full derivation in ../README.md):
 *   C0: H_A(p) = H_B(p)    for p on the shared edge  (no height crack)
 *   C1: dH_A/dn = dH_B/dn  normal derivative match   (no lighting seam)
 *
 * Construction: every tile samples one global stitched function
 *   H(p) = sum_k w_k(p) * N_k(p),   sum_k w_k(p) = 1,  w_k in C^1,
 * where N_k is tile k's own seeded fBm and w_k a smoothstep partition of
 * unity overlapping MU of a tile width. A shared edge vertex evaluates the
 * identical float expression no matter which tile asks -> C0 bit-exact.
 * H is C^1 as a function, so cross-seam (ghost-layer) normals agree
 * bit-exactly as well. A clamped per-tile renderer keeps only an
 * O(kappa*cell) finite-difference truncation residual, bounded below.
 */
#include <stdio.h>
#include <math.h>

#define TW 16                          /* vertices per tile side, shared edges */
#define NT 2                           /* 2x2 tiles: vertical + horizontal + corner junctions */
#define GW (NT * (TW - 1) + 1)         /* global vertex grid side */

#define CELL 1.0f
#define MU 0.25f                       /* blend overlap, fraction of a tile width */
#define OCT 4
#define LAC 2.0f
#define GAIN 0.5f
#define AMP 3.0f
#define FREQ 0.045f

#define PI 3.14159265358979323846f

/* noise: same algorithm as ../../noise/src/noise.c, copied so terrain/ builds standalone */
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

typedef struct { float h[TW][TW]; } Tile;
typedef struct { float x, y, z; } Vec3;

/* weight of the tile owning local coordinate u (u may reach MU past either end):
 * 1 in the interior, smoothstep ramp across the MU overlap, 0 outside.
 * w_i(u) + w_{i+1}(u) = 1 in the overlap because smooth01(t)+smooth01(1-t)=1. */
static float smooth01(float t) { return t * t * (3.0f - 2.0f * t); }

static float axis_weight(float u) {
    if (u < -MU || u > 1.0f + MU) return 0.0f;
    if (u < MU) return smooth01((u + MU) / (2.0f * MU));
    if (u > 1.0f - MU) return smooth01((1.0f + MU - u) / (2.0f * MU));
    return 1.0f;
}

/* Global stitched height: sum over the (<=3x3) tiles whose weight support
 * covers p. Loop order is fixed, so every tile evaluating a shared vertex
 * runs the identical float expression -> bit-identical results. */
static float stitched_height(unsigned base, float wx, float wy) {
    const float T = (TW - 1) * CELL;
    float gx = wx / T, gy = wy / T;
    int ix = (int)floorf(gx), iy = (int)floorf(gy);
    float sum = 0.0f;
    for (int k = ix - 1; k <= ix + 1; k++) {
        float wk = axis_weight(gx - (float)k);
        if (wk == 0.0f) continue;
        for (int l = iy - 1; l <= iy + 1; l++) {
            float wl = axis_weight(gy - (float)l);
            if (wl == 0.0f) continue;
            sum += wk * wl * fbm(wx * FREQ, wy * FREQ, OCT, LAC, GAIN, hash2(k, l, base));
        }
    }
    return AMP * sum;
}

static void tile_gen(Tile* t, unsigned base, int tx, int ty) {
    for (int x = 0; x < TW; x++)
        for (int y = 0; y < TW; y++)
            t->h[x][y] = stitched_height(base,
                (float)(tx * (TW - 1) + x) * CELL,
                (float)(ty * (TW - 1) + y) * CELL);
}

static Vec3 mknormal(float dx, float dy) {
    float inv = 1.0f / sqrtf(dx * dx + dy * dy + 1.0f);
    Vec3 n = { -dx * inv, -dy * inv, inv };
    return n;
}

static float angle_deg(Vec3 a, Vec3 b) {
    /* atan2(|a x b|, a . b): accurate near 0, unlike acos(dot) */
    float cx = a.y * b.z - a.z * b.y;
    float cy = a.z * b.x - a.x * b.z;
    float cz = a.x * b.y - a.y * b.x;
    float cross = sqrtf(cx * cx + cy * cy + cz * cz);
    float dot = a.x * b.x + a.y * b.y + a.z * b.z;
    return atan2f(cross, dot) * 180.0f / PI;
}

/* per-tile clamped central difference, same convention as field_slope() */
static Vec3 tile_normal(const Tile* t, int x, int y) {
    int xm = x > 0 ? x - 1 : x, xp = x < TW - 1 ? x + 1 : x;
    int ym = y > 0 ? y - 1 : y, yp = y < TW - 1 ? y + 1 : y;
    float dx = (t->h[xp][y] - t->h[xm][y]) / ((float)(xp - xm) * CELL);
    float dy = (t->h[x][yp] - t->h[x][ym]) / ((float)(yp - ym) * CELL);
    return mknormal(dx, dy);
}

typedef struct { float dh, clamp_deg, ghost_deg; } SeamStats;

/* vertical seam: L's right column (x=TW-1) against R's left column (x=0) */
static SeamStats seam_vert(const Tile* L, const Tile* R) {
    SeamStats s = { 0.0f, 0.0f, 0.0f };
    for (int y = 0; y < TW; y++) {
        float dh = fabsf(L->h[TW - 1][y] - R->h[0][y]);
        if (dh > s.dh) s.dh = dh;

        float a = angle_deg(tile_normal(L, TW - 1, y), tile_normal(R, 0, y));
        if (a > s.clamp_deg) s.clamp_deg = a;

        /* ghost-layer normals: cross-seam central difference, each tile
         * pairing its own interior column with the neighbor's first column */
        int ym = y > 0 ? y - 1 : y, yp = y < TW - 1 ? y + 1 : y;
        float dx = (R->h[1][y] - L->h[TW - 2][y]) / (2.0f * CELL);
        Vec3 gL = mknormal(dx, (L->h[TW - 1][yp] - L->h[TW - 1][ym]) / ((float)(yp - ym) * CELL));
        Vec3 gR = mknormal(dx, (R->h[0][yp] - R->h[0][ym]) / ((float)(yp - ym) * CELL));
        float g = angle_deg(gL, gR);
        if (g > s.ghost_deg) s.ghost_deg = g;
    }
    return s;
}

/* horizontal seam: B's bottom row (y=TW-1) against U's top row (y=0) */
static SeamStats seam_horz(const Tile* B, const Tile* U) {
    SeamStats s = { 0.0f, 0.0f, 0.0f };
    for (int x = 0; x < TW; x++) {
        float dh = fabsf(B->h[x][TW - 1] - U->h[x][0]);
        if (dh > s.dh) s.dh = dh;

        float a = angle_deg(tile_normal(B, x, TW - 1), tile_normal(U, x, 0));
        if (a > s.clamp_deg) s.clamp_deg = a;

        int xm = x > 0 ? x - 1 : x, xp = x < TW - 1 ? x + 1 : x;
        float dy = (U->h[x][1] - B->h[x][TW - 2]) / (2.0f * CELL);
        Vec3 gB = mknormal((B->h[xp][TW - 1] - B->h[xm][TW - 1]) / ((float)(xp - xm) * CELL), dy);
        Vec3 gU = mknormal((U->h[xp][0] - U->h[xm][0]) / ((float)(xp - xm) * CELL), dy);
        float g = angle_deg(gB, gU);
        if (g > s.ghost_deg) s.ghost_deg = g;
    }
    return s;
}

/* corner where 4 tiles meet: t00(SE) t10(SW) t01(NE) t11(NW) */
static SeamStats seam_corner(const Tile* t00, const Tile* t10,
                             const Tile* t01, const Tile* t11) {
    SeamStats s = { 0.0f, 0.0f, 0.0f };
    float hv[4] = { t00->h[TW - 1][TW - 1], t10->h[0][TW - 1],
                    t01->h[TW - 1][0],      t11->h[0][0] };
    Vec3 clamp_n[4] = { tile_normal(t00, TW - 1, TW - 1), tile_normal(t10, 0, TW - 1),
                        tile_normal(t01, TW - 1, 0),      tile_normal(t11, 0, 0) };
    /* ghost normals from each tile's perspective, reading cross-tile neighbors */
    float c2 = 2.0f * CELL;
    Vec3 ghost_n[4] = {
        mknormal((t10->h[1][TW - 1] - t00->h[TW - 2][TW - 1]) / c2,
                 (t01->h[TW - 1][1] - t00->h[TW - 1][TW - 2]) / c2),
        mknormal((t10->h[1][TW - 1] - t00->h[TW - 2][TW - 1]) / c2,
                 (t11->h[0][1] - t10->h[0][TW - 2]) / c2),
        mknormal((t11->h[1][0] - t01->h[TW - 2][0]) / c2,
                 (t01->h[TW - 1][1] - t00->h[TW - 1][TW - 2]) / c2),
        mknormal((t11->h[1][0] - t01->h[TW - 2][0]) / c2,
                 (t11->h[0][1] - t10->h[0][TW - 2]) / c2),
    };
    for (int i = 0; i < 4; i++)
        for (int j = i + 1; j < 4; j++) {
            float dh = fabsf(hv[i] - hv[j]);
            if (dh > s.dh) s.dh = dh;
            float a = angle_deg(clamp_n[i], clamp_n[j]);
            if (a > s.clamp_deg) s.clamp_deg = a;
            float g = angle_deg(ghost_n[i], ghost_n[j]);
            if (g > s.ghost_deg) s.ghost_deg = g;
        }
    return s;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    const unsigned base = 20261124u;
    static Tile t[NT][NT];
    for (int ty = 0; ty < NT; ty++)
        for (int tx = 0; tx < NT; tx++)
            tile_gen(&t[tx][ty], base, tx, ty);

    /* baseline: independent per-tile seeds, no blending -> the reported fault */
    float base_dh = 0.0f;
    for (int y = 0; y < TW; y++) {
        float wx = (TW - 1) * CELL, wy = (float)y * CELL;
        float ha = AMP * fbm(wx * FREQ, wy * FREQ, OCT, LAC, GAIN, hash2(0, 0, base));
        float hb = AMP * fbm(wx * FREQ, wy * FREQ, OCT, LAC, GAIN, hash2(1, 0, base));
        float d = fabsf(ha - hb);
        if (d > base_dh) base_dh = d;
    }
    printf("baseline_s1_max_dh=%.6g\n", base_dh);

    /* three junctions: S1 vertical, S2 horizontal, S3 four-tile corner */
    SeamStats s1 = seam_vert(&t[0][0], &t[1][0]);
    SeamStats s2 = seam_horz(&t[0][0], &t[0][1]);
    SeamStats s3 = seam_corner(&t[0][0], &t[1][0], &t[0][1], &t[1][1]);
    printf("s1_max_dh=%.6g s1_clamped_deg=%.6g s1_ghost_deg=%.6g\n", s1.dh, s1.clamp_deg, s1.ghost_deg);
    printf("s2_max_dh=%.6g s2_clamped_deg=%.6g s2_ghost_deg=%.6g\n", s2.dh, s2.clamp_deg, s2.ghost_deg);
    printf("s3_max_dh=%.6g s3_clamped_deg=%.6g s3_ghost_deg=%.6g\n", s3.dh, s3.clamp_deg, s3.ghost_deg);

    /* kappa: max |second derivative| over the assembled global grid,
     * estimated by central second differences */
    float G[GW][GW];
    for (int ty = 0; ty < NT; ty++)
        for (int tx = 0; tx < NT; tx++)
            for (int x = 0; x < TW; x++)
                for (int y = 0; y < TW; y++)
                    G[tx * (TW - 1) + x][ty * (TW - 1) + y] = t[tx][ty].h[x][y];
    float kappa = 0.0f;
    for (int x = 1; x < GW - 1; x++)
        for (int y = 1; y < GW - 1; y++) {
            float kx = fabsf(G[x + 1][y] - 2.0f * G[x][y] + G[x - 1][y]) / (CELL * CELL);
            float ky = fabsf(G[x][y + 1] - 2.0f * G[x][y] + G[x][y - 1]) / (CELL * CELL);
            if (kx > kappa) kappa = kx;
            if (ky > kappa) kappa = ky;
        }
    printf("kappa=%.6g\n", kappa);

    /* thresholds (derivations in ../README.md):
     * dh:      construction is bit-exact; tolerance is float-rounding headroom
     *          eps(1.2e-7) * |h|max(~AMP) rounded up -> 1e-6*AMP
     * ghost:   Weber 2% luminance JND, |dn| ~= dtheta -> 0.02 rad = 1.1459 deg
     * clamped: forward/backward difference err <= (cell/2)*kappa each,
     *          gap <= kappa*cell; 1.5x safety for sampled-kappa underestimate */
    const float DH_THRESH = 1e-6f * AMP;
    const float GHOST_THRESH_DEG = 0.02f * 180.0f / PI;
    const float CLAMP_THRESH_DEG = 1.5f * kappa * CELL * 180.0f / PI;
    printf("dh_thresh=%.6g ghost_thresh_deg=%.6g clamp_thresh_deg=%.6g\n",
           DH_THRESH, GHOST_THRESH_DEG, CLAMP_THRESH_DEG);

    check("s1_c0_height", s1.dh <= DH_THRESH);
    check("s1_c1_ghost_normal", s1.ghost_deg <= GHOST_THRESH_DEG);
    check("s1_clamped_is_truncation_only", s1.clamp_deg <= CLAMP_THRESH_DEG);
    check("s2_c0_height", s2.dh <= DH_THRESH);
    check("s2_c1_ghost_normal", s2.ghost_deg <= GHOST_THRESH_DEG);
    check("s2_clamped_is_truncation_only", s2.clamp_deg <= CLAMP_THRESH_DEG);
    check("s3_c0_height", s3.dh <= DH_THRESH);
    check("s3_c1_ghost_normal", s3.ghost_deg <= GHOST_THRESH_DEG);
    check("s3_clamped_is_truncation_only", s3.clamp_deg <= CLAMP_THRESH_DEG);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
