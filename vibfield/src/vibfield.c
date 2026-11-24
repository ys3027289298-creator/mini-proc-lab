
/* vibfield: heightmap -> predicted vibration states, density-based frequency
 * extrapolation into a border buffer, out-of-bounds cell reporting and a
 * rockfall check. Low-altitude fine texture creates pseudo-vibration, which is
 * isolated by 3x3 smoothing and suppressed on low, gentle, rough cells.
 *
 * Usage:
 *   vibfield                 self-test (check lines, all_ok)
 *   vibfield runs FILTERED 10 machine-readable RUN/OOB records, 10 texture seeds
 *   vibfield runs RAW 10      same pipeline without pseudo-vibration suppression
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define N 40                /* core heightmap size */
#define B 4                 /* buffer width around the core */
#define S (N + 2 * B)       /* working grid size */
#define CELL 1.0f

/* states */
#define CALM 0
#define GENTLE 1
#define ACTIVE 2
#define CLIFF 3

/* thresholds (tuned against the built-in terrain) */
#define S_CLIFF   1.50f     /* smoothed slope marking a cliff face */
#define S_ACTIVE  0.45f     /* smoothed slope marking active vibration */
#define S_GENTLE  0.12f     /* smoothed slope marking gentle vibration */
#define C_ACTIVE  0.80f     /* smoothed |laplacian| marking active vibration */
#define Z_LOW     2.20f     /* below this elevation texture dominates easily */
#define R_TEX     0.05f     /* raw-minus-smoothed height residual = texture */
#define RHO_WIN   5         /* density window (square, odd) */
#define RHO_MIN   0.06f     /* minimum extrapolated density that carries a tone */
#define F_BASE    1.00f     /* Hz floor for a vibrating cell */
#define K_RHO     12.0f     /* Hz per unit density */
#define F_ALARM   4.50f     /* buffer frequency above which a cell is OOB */
#define S_BLOCK   0.12f     /* flatter than this arrests a rolling rock */
#define MAX_OOB   1024

typedef struct {
    float base[S][S];      /* texture-free macro terrain (ramp + cliff band) */
    float raw[S][S];       /* heightmap as generated, texture included */
    float sm[S][S];        /* 3x3 box-smoothed heightmap */
    int   st[S][S];        /* final predicted vibration state */
    int   pseudo[S][S];    /* 1 where raw vibration was suppressed as texture */
    float rho[S][S];       /* vibrating-cell density, extrapolated into buffer */
    float freq[S][S];      /* predicted frequency, buffer included */
} Grid;

/* ---- deterministic value noise (same family as the noise module) ---- */

static unsigned hash2(int x, int y, unsigned seed) {
    unsigned h = (unsigned)x * 374761393u + (unsigned)y * 668265263u
               + seed * 2246822519u;
    h = (h ^ (h >> 13)) * 1274126177u;
    return h ^ (h >> 16);
}

static float rand01(unsigned h) { return (float)(h & 0xffffff) / 16777215.0f; }

static float lerpf(float a, float b, float t) { return a + (b - a) * t; }
static float fade(float t) { return t * t * (3.0f - 2.0f * t); }
static float smoothstepf(float e0, float e1, float x) {
    float t = (x - e0) / (e1 - e0);
    if (t < 0.0f) t = 0.0f;
    if (t > 1.0f) t = 1.0f;
    return t * t * (3.0f - 2.0f * t);
}

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

static float fbm2(float x, float y, unsigned seed) {
    float a = value_noise(x, y, seed) * 0.65f;
    float b = value_noise(x * 2.0f, y * 2.0f, seed + 17u) * 0.35f;
    return a + b;
}

/* ---- heightmap: fixed gentle ramp x diagonal cliff intersection, per-seed
 *      texture. The cliff band crosses the west edge in the low zone, so its
 *      steep face and the gentle textured slope intersect right at the buffer. */

static void height_build(Grid* g, unsigned seed) {
    memset(g, 0, sizeof(*g));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            /* gentle ramp, slope ~0.12 */
            float zbase = 0.10f * (float)i + 0.06f * (float)j;
            /* diagonal cliff band: a 7 m drop spread over ~2.5 cells */
            float d = (float)i + 0.6f * (float)j - 20.0f;
            zbase -= 7.0f * smoothstepf(0.0f, 2.5f, d);
            g->base[B + i][B + j] = zbase;
            /* low-altitude fine texture, fading out above the low zone */
            float amp = 0.06f + 0.45f * (1.0f - smoothstepf(0.8f, Z_LOW,
                                                             zbase - g->base[B][B]));
            float tex = fbm2((float)i * 0.45f, (float)j * 0.45f, seed) * amp;
            g->raw[B + i][B + j] = zbase + tex;
        }
    }
}

static int clamp_core(int k) {
    if (k < B) return B;
    if (k >= B + N) return B + N - 1;
    return k;
}

/* Low zone = at most Z_LOW above the terrain floor (the lowest base cell).
 * Absolute elevation is unusable because the cliff drop shifts the whole map
 * below any fixed cutoff. */
static float terrain_floor(const Grid* g) {
    float z = g->base[B][B];
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            float b = g->base[B + i][B + j];
            if (b < z) z = b;
        }
    return z;
}

static int in_low_zone(const Grid* g, int r, int c) {
    return g->sm[r][c] < terrain_floor(g) + Z_LOW;
}

static void smooth3(Grid* g) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            float sum = 0.0f;
            for (int di = -1; di <= 1; di++)
                for (int dj = -1; dj <= 1; dj++)
                    sum += g->raw[clamp_core(B + i + di)][clamp_core(B + j + dj)];
            g->sm[B + i][B + j] = sum / 9.0f;
        }
    }
}

static float slope_at(const float h[S][S], int r, int c) {
    int rm = r > 0 ? r - 1 : r, rp = r < S - 1 ? r + 1 : r;
    int cm = c > 0 ? c - 1 : c, cp = c < S - 1 ? c + 1 : c;
    float dr = (h[rp][c] - h[rm][c]) / (float)(rp - rm);
    float dc = (h[r][cp] - h[r][cm]) / (float)(cp - cm);
    return sqrtf(dr * dr + dc * dc);
}

static float lap_at(const float h[S][S], int r, int c) {
    int rm = r > 0 ? r - 1 : r, rp = r < S - 1 ? r + 1 : r;
    int cm = c > 0 ? c - 1 : c, cp = c < S - 1 ? c + 1 : c;
    return h[rm][c] + h[rp][c] + h[r][cm] + h[r][cp] - 4.0f * h[r][c];
}

static int classify(float slope, float curv) {
    if (slope >= S_CLIFF) return CLIFF;
    if (slope >= S_ACTIVE || fabsf(curv) >= C_ACTIVE) return ACTIVE;
    if (slope >= S_GENTLE) return GENTLE;
    return CALM;
}

/* Convert the heightmap to predicted vibration states. FILTERED classifies the
 * smoothed field and suppresses texture-only vibration; RAW trusts the raw map.
 * A cell is pseudo-vibration when the raw map looks ACTIVE+ but that energy
 * vanishes under smoothing, at low elevation: it is texture, not structure. */
static void states_build(Grid* g, int filtered) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int r = B + i, c = B + j;
            int raw_st = classify(slope_at(g->raw, r, c), lap_at(g->raw, r, c));
            if (!filtered) { g->st[r][c] = raw_st; continue; }

            int st = classify(slope_at(g->sm, r, c), lap_at(g->sm, r, c));
            float residual = fabsf(g->raw[r][c] - g->sm[r][c]);
            if (raw_st >= ACTIVE && st <= GENTLE &&
                in_low_zone(g, r, c) && residual >= R_TEX) {
                st = CALM;
                g->pseudo[r][c] = 1;
            }
            g->st[r][c] = st;
        }
    }
}

static void density_build(Grid* g) {
    int half = RHO_WIN / 2;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int hit = 0;
            for (int di = -half; di <= half; di++)
                for (int dj = -half; dj <= half; dj++) {
                    int rr = clamp_core(B + i + di), cc = clamp_core(B + j + dj);
                    if (g->st[rr][cc] >= ACTIVE) hit++;
                }
            g->rho[B + i][B + j] = (float)hit / (float)(RHO_WIN * RHO_WIN);
        }
    }
}

/* Linear gradient extrapolation from the two outermost core layers out to the
 * buffer boundary. Vibration density cannot grow away from the source, so the
 * outward gradient is clipped to <= 0 and density stays inside [0, 1]. The
 * first extrapolated cell equals the boundary value (no jump). */
static void density_extrapolate(Grid* g) {
    float (*rho)[S] = g->rho;
    for (int c = B; c < B + N; c++) {
        float gt = fminf(0.0f, rho[B][c] - rho[B + 1][c]);
        float gb = fminf(0.0f, rho[B + N - 1][c] - rho[B + N - 2][c]);
        for (int d = 1; d <= B; d++) {
            rho[B - d][c] = fmaxf(0.0f, rho[B][c] + gt * (d - 1));
            rho[B + N - 1 + d][c] = fmaxf(0.0f, rho[B + N - 1][c] + gb * (d - 1));
        }
    }
    for (int r = B; r < B + N; r++) {
        float gl = fminf(0.0f, rho[r][B] - rho[r][B + 1]);
        float gr = fminf(0.0f, rho[r][B + N - 1] - rho[r][B + N - 2]);
        for (int d = 1; d <= B; d++) {
            rho[r][B - d] = fmaxf(0.0f, rho[r][B] + gl * (d - 1));
            rho[r][B + N - 1 + d] = fmaxf(0.0f, rho[r][B + N - 1] + gr * (d - 1));
        }
    }
    for (int r = 0; r < S; r++) {
        for (int c = 0; c < S; c++) {
            int in_r = r >= B && r < B + N, in_c = c >= B && c < B + N;
            if (in_r && in_c) continue;
            if (!in_r && !in_c) {
                int dr = r < B ? B - r : r - (B + N - 1);
                int dc = c < B ? B - c : c - (B + N - 1);
                rho[r][c] = dr >= dc ? rho[r][clamp_core(c)]
                                     : rho[clamp_core(r)][c];
            } else if (!in_r) {
                rho[r][c] = rho[r][clamp_core(c)];
            } else {
                rho[r][c] = rho[clamp_core(r)][c];
            }
            if (rho[r][c] > 1.0f) rho[r][c] = 1.0f;
        }
    }
}

static void frequency_build(Grid* g) {
    for (int r = 0; r < S; r++) {
        for (int c = 0; c < S; c++) {
            int in_core = r >= B && r < B + N && c >= B && c < B + N;
            if (in_core) {
                g->freq[r][c] = g->st[r][c] == CALM
                                    ? 0.0f
                                    : F_BASE + K_RHO * g->rho[r][c];
            } else {
                g->freq[r][c] = g->rho[r][c] >= RHO_MIN
                                    ? F_BASE + K_RHO * g->rho[r][c]
                                    : 0.0f;
            }
        }
    }
}

typedef struct {
    int r, c;               /* working-grid coordinates */
    char side;              /* N / S / W / E */
    int dist;               /* cells outside the core boundary */
    float rho, freq;
    int src_state;
    int s1, s2, s3, s4, risk; /* rockfall check steps + verdict */
} Oob;

/* Rockfall check for one buffer cell:
 * 1. source cell (nearest core cell) is a CLIFF face;
 * Rockfall is a structural process, so gates 1-3 read the texture-free base
 * terrain rather than the textured/smoothed fields:
 * 1. the source cell (nearest core cell) is a cliff face on the base terrain;
 * 2. base elevation drops toward the buffer (the rock rolls outward);
 * 3. a 3-cell inward trace stays steep and monotonic, no flat arrest zone;
 * 4. the extrapolated frequency clears F_ALARM (implicit: the cell is OOB).
 * A rockfall verdict requires all four gates. */
static void rockfall_check(const Grid* g, const Oob* probe, int out[4]) {
    int sr = clamp_core(probe->r), sc = clamp_core(probe->c);
    out[0] = slope_at(g->base, sr, sc) >= S_CLIFF;

    int dr = 0, dc = 0;
    if (probe->side == 'N') dr = -1;
    else if (probe->side == 'S') dr = 1;
    else if (probe->side == 'W') dc = -1;
    else dc = 1;
    /* outward downhill: the cell just inside is higher than the boundary */
    int ir = clamp_core(sr - dr), ic = clamp_core(sc - dc);
    out[1] = g->base[ir][ic] > g->base[sr][sc];

    int open = 1;
    int pr = sr, pc = sc;
    for (int k = 0; k < 3; k++) {
        int nr = clamp_core(pr - dr), nc = clamp_core(pc - dc);
        if (nr == pr && nc == pc) break;
        if (!(g->base[nr][nc] > g->base[pr][pc])) { open = 0; break; }
        if (slope_at(g->base, nr, nc) < S_BLOCK) { open = 0; break; }
        pr = nr; pc = nc;
    }
    out[2] = open;
    out[3] = probe->freq > F_ALARM;
}

static int oob_collect(const Grid* g, Oob* out, int maxout) {
    int n = 0;
    for (int r = 0; r < S; r++) {
        for (int c = 0; c < S; c++) {
            int in_core = r >= B && r < B + N && c >= B && c < B + N;
            if (in_core || g->freq[r][c] <= F_ALARM) continue;
            if (n >= maxout) return n;
            Oob* o = &out[n++];
            o->r = r; o->c = c;
            if (r < B) { o->side = 'N'; o->dist = B - r; }
            else if (r >= B + N) { o->side = 'S'; o->dist = r - (B + N - 1); }
            else if (c < B) { o->side = 'W'; o->dist = B - c; }
            else { o->side = 'E'; o->dist = c - (B + N - 1); }
            o->rho = g->rho[r][c];
            o->freq = g->freq[r][c];
            o->src_state = g->st[clamp_core(r)][clamp_core(c)];
            int steps[4];
            rockfall_check(g, o, steps);
            o->s1 = steps[0]; o->s2 = steps[1]; o->s3 = steps[2]; o->s4 = steps[3];
            o->risk = steps[0] && steps[1] && steps[2] && steps[3];
        }
    }
    return n;
}

static void grid_run(Grid* g, unsigned seed, int filtered, Oob* oob, int maxout,
                     int* n_oob, int* n_pseudo) {
    height_build(g, seed);
    smooth3(g);
    states_build(g, filtered);
    density_build(g);
    density_extrapolate(g);
    frequency_build(g);
    *n_oob = oob_collect(g, oob, maxout);
    int p = 0;
    for (int r = 0; r < S; r++)
        for (int c = 0; c < S; c++)
            p += g->pseudo[r][c];
    *n_pseudo = p;
}

static const char* state_name(int st) {
    static const char* nm[4] = {"CALM", "GENTLE", "ACTIVE", "CLIFF"};
    return nm[st];
}

static void print_map(const Grid* g) {
    const char ch[4] = {'.', '~', 'o', '#'};
    printf("map(.:calm ~:gentle o:active #:cliff P:pseudo-suppressed)\n");
    for (int i = 0; i < N; i++) {
        char line[N + 1];
        for (int j = 0; j < N; j++) {
            int r = B + i, c = B + j;
            line[j] = g->pseudo[r][c] ? 'P' : ch[g->st[r][c]];
        }
        line[N] = 0;
        printf("%s\n", line);
    }
}

/* ---- self-test ---- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(int argc, char** argv) {
    if (argc == 4 && strcmp(argv[1], "runs") == 0) {
        int filtered = strcmp(argv[2], "RAW") != 0;
        int runs = atoi(argv[3]);
        for (int k = 0; k < runs; k++) {
            Grid g;
            Oob oob[MAX_OOB];
            int n_oob, n_pseudo;
            unsigned seed = 100u + (unsigned)k;
            grid_run(&g, seed, filtered, oob, MAX_OOB, &n_oob, &n_pseudo);
            printf("RUN,%d,%u,%d,%d\n", k, seed, n_oob, n_pseudo);
            for (int i = 0; i < n_oob; i++) {
                Oob* o = &oob[i];
                printf("OOB,%d,%d,%d,%c,%d,%.4f,%.4f,%s,%d,%d,%d,%d,%d\n",
                       k, o->r, o->c, o->side, o->dist, o->rho, o->freq,
                       state_name(o->src_state), o->s1, o->s2, o->s3, o->s4,
                       o->risk);
            }
        }
        return 0;
    }

    Grid g, g2;
    Oob oob[MAX_OOB], oob2[MAX_OOB];
    int n_oob, n_pseudo, n2, p2;

    /* texture removal: in the low zone, raw-minus-smoothed height residuals
     * are substantial while the smoothed gentle-ramp slope stays small */
    grid_run(&g, 100u, 1, oob, MAX_OOB, &n_oob, &n_pseudo);
    float res_energy = 0.0f, sm_slope = 0.0f; int nlow = 0;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            int r = B + i, c = B + j;
            if (!in_low_zone(&g, r, c)) continue;
            res_energy += fabsf(g.raw[r][c] - g.sm[r][c]);
            sm_slope += slope_at(g.sm, r, c);
            nlow++;
        }
    check("smoothing_strips_texture",
          nlow > 0 && (res_energy / nlow) > R_TEX &&
          (sm_slope / nlow) < S_ACTIVE);
    check("pseudo_detected", n_pseudo >= 8);
    int pseudo_low_only = 1;
    for (int r = 0; r < S && pseudo_low_only; r++)
        for (int c = 0; c < S; c++)
            if (g.pseudo[r][c] && !in_low_zone(&g, r, c)) pseudo_low_only = 0;
    check("pseudo_only_in_low_zone", pseudo_low_only);
    print_map(&g);

    /* cliff band intersects the gentle ramp and reaches the low, textured zone */
    int cliffs = 0, low_cliffs = 0;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            int r = B + i, c = B + j;
            if (g.st[r][c] == CLIFF) {
                cliffs++;
                if (in_low_zone(&g, r, c)) low_cliffs++;
            }
        }
    check("cliff_present", cliffs >= 4);
    check("cliff_reaches_low_zone", low_cliffs >= 1);

    /* raw mode trusts texture: every cell the filtered run suppresses as
     * pseudo-vibration is classified ACTIVE+ by the raw run, and no
     * suppression happens there at all. */
    grid_run(&g2, 100u, 0, oob2, MAX_OOB, &n2, &p2);
    int suppressed_but_raw_active = 0, low_n = 0;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            int r = B + i, c = B + j;
            if (!in_low_zone(&g, r, c)) continue;
            low_n++;
            if (g.pseudo[r][c] && g2.st[r][c] >= ACTIVE)
                suppressed_but_raw_active++;
        }
    check("filter_reduces_texture_vibration",
          low_n > 0 && suppressed_but_raw_active == n_pseudo && n_pseudo > 0);
    check("raw_mode_has_no_suppression", p2 == 0);

    /* extrapolation: bounded density, finite values, continuous at boundary */
    int nonneg = 1, bounded = 1, finite = 1;
    float boundary_gap = 0.0f;
    for (int r = 0; r < S; r++)
        for (int c = 0; c < S; c++) {
            if (g.rho[r][c] < 0.0f) nonneg = 0;
            if (g.rho[r][c] > 1.0f) bounded = 0;
            if (!isfinite(g.rho[r][c]) || !isfinite(g.freq[r][c])) finite = 0;
        }
    for (int c = B; c < B + N; c++) {
        boundary_gap += fabsf(g.rho[B - 1][c] - g.rho[B][c]);
        boundary_gap += fabsf(g.rho[B + N][c] - g.rho[B + N - 1][c]);
    }
    for (int r = B; r < B + N; r++) {
        boundary_gap += fabsf(g.rho[r][B - 1] - g.rho[r][B]);
        boundary_gap += fabsf(g.rho[r][B + N] - g.rho[r][B + N - 1]);
    }
    check("extrap_bounded_finite", nonneg && bounded && finite);
    check("extrap_boundary_continuous", boundary_gap < 1e-4f * (4 * N));

    /* rockfall verdicts exist and every step is a hard gate */
    int risks = 0, gated = 1;
    for (int i = 0; i < n_oob; i++) {
        risks += oob[i].risk;
        if (oob[i].risk && !(oob[i].s1 && oob[i].s2 && oob[i].s3 && oob[i].s4))
            gated = 0;
    }
    check("rockfall_cells_found", risks > 0);
    check("rockfall_gates_and", gated);

    /* determinism: same seed reproduces the exact OOB set */
    grid_run(&g2, 100u, 1, oob2, MAX_OOB, &n2, &p2);
    int same = n2 == n_oob && p2 == n_pseudo;
    for (int i = 0; same && i < n_oob; i++)
        if (oob[i].r != oob2[i].r || oob[i].c != oob2[i].c ||
            oob[i].risk != oob2[i].risk)
            same = 0;
    check("deterministic", same);

    /* ten texture seeds over the fixed terrain: rockfall OOB cells are
     * structural -> 10/10 repeat; texture-only OOB cells move with the seed */
    int seen[S][S] = {{0}}, risk_seen[S][S] = {{0}};
    int totals[10];
    for (int k = 0; k < 10; k++) {
        Grid gk;
        Oob ok[MAX_OOB];
        int nk, pk;
        grid_run(&gk, 100u + (unsigned)k, 1, ok, MAX_OOB, &nk, &pk);
        totals[k] = nk;
        for (int i = 0; i < nk; i++) {
            seen[ok[i].r][ok[i].c]++;
            if (ok[i].risk) risk_seen[ok[i].r][ok[i].c]++;
        }
    }
    int structural_stable = 1, risk_total = 0;
    for (int r = 0; r < S; r++)
        for (int c = 0; c < S; c++)
            if (risk_seen[r][c] > 0) {
                risk_total++;
                if (risk_seen[r][c] != 10) structural_stable = 0;
            }
    check("rockfall_cells_repeat_10_of_10", structural_stable && risk_total > 0);
    int varying = 0;
    for (int k = 1; k < 10; k++)
        if (totals[k] != totals[0]) varying = 1;
    check("texture_seeds_move_oob_count", varying);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
