
/* veins: low-resolution vein density field + density-driven Poisson sampling.
 *
 * The old scatter pipeline placed one rock per grid cell independently, so
 * each rock type ended up uniformly spread with no clustered ore veins.
 * Here the pipeline is:
 *   1. throw sparse anchors on a coarse (low-resolution) lattice,
 *   2. splat Gaussian kernels of width sigma onto the coarse grid, building
 *      a smooth vein-density field (bilinear-sampled at fine positions),
 *   3. Poisson-disk sampling: accept a candidate only where density >= TAU
 *      and it keeps distance R0 from every accepted point, hard-capped at
 *      NMAX points.
 *
 * Smoothing scale sigma has two failure regimes:
 *   - too coarse: neighbouring splats merge into one featureless band, the
 *     vein degenerates into a straight line;
 *   - too fine: the thresholded region shrinks to isolated specks and the
 *     sample set degenerates back into independent scatter.
 *
 * The scale is fixed jointly by two curves:
 *   length vs scale : L(s) = cL * s
 *   count  vs scale : N(s) = min(NMAX, cN * s^2)
 * cL and cN are calibrated once on an uncapped pilot run at s = SIGMA0.
 * A target length L* asks for s* = L* / cL; the second curve says how many
 * samples that demands. If cN*s*^2 > NMAX the sample cap binds: the field
 * keeps growing but the points cannot fill it, continuity (measured on the
 * point graph, link radius RLINK) fragments and the longest vein shrinks.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define DOMAIN    64.0f   /* world is [0, DOMAIN) x [0, DOMAIN) */
#define CCELL     4.0f    /* coarse (low-res) field cell size */
#define CGW       16      /* coarse grid is CGW x CGW nodes */
#define ANCH      3       /* ANCH x ANCH jittered anchor lattice */
#define TAU       0.5f    /* density threshold for vein ground */
#define R0        2.5f    /* Poisson minimum distance */
#define RLINK     (1.35f * R0) /* points within this link into one vein */
#define NMAX      160     /* hard cap on accepted samples */
#define CSTEP     1.5f    /* candidate stream grid step */
#define NCELL     ((int)(DOMAIN / CSTEP))
#define SIGMA0    4.0f    /* pilot scale used to calibrate cL, cN */

typedef struct { float x, y; } P;

typedef struct {
    P p[NMAX];
    int n;            /* accepted samples */
    int veins;        /* graph components with >= 2 points */
    int singles;      /* components with 1 point */
    int capped;       /* NMAX was hit */
    float mean_len;   /* mean component diameter over veins */
    float longest;    /* largest component diameter */
} Run;

static unsigned h2(int x, int y, unsigned salt) {
    unsigned h = (unsigned)x * 374761393u + (unsigned)y * 668265263u
               + salt * 2246822519u;
    h = (h ^ (h >> 13)) * 1274126177u;
    return h ^ (h >> 16);
}

static float u01(int x, int y, unsigned salt) {
    return (float)(h2(x, y, salt) & 0xffffff) / 16777215.0f;
}

/* ANCH x ANCH jittered anchors, spread over the whole domain. */
static void make_anchors(P* a, unsigned seed) {
    float step = DOMAIN / ANCH;
    int k = 0;
    for (int i = 0; i < ANCH; i++)
        for (int j = 0; j < ANCH; j++) {
            float jx = (u01(i, j, seed) - 0.5f) * step * 0.8f;
            float jy = (u01(j, i, seed + 1u) - 0.5f) * step * 0.8f;
            a[k].x = (i + 0.5f) * step + jx;
            a[k].y = (j + 0.5f) * step + jy;
            k++;
        }
}

/* Low-resolution density field: Gaussian splats accumulated on coarse nodes. */
static void make_field(float field[CGW][CGW], const P* a, int na, float sigma) {
    memset(field, 0, sizeof(float) * CGW * CGW);
    float s2 = 2.0f * sigma * sigma;
    for (int k = 0; k < na; k++) {
        int gx0 = (int)(a[k].x / CCELL), gy0 = (int)(a[k].y / CCELL);
        int reach = (int)ceilf(3.0f * sigma / CCELL);
        for (int gx = gx0 - reach; gx <= gx0 + reach; gx++)
            for (int gy = gy0 - reach; gy <= gy0 + reach; gy++) {
                if (gx < 0 || gy < 0 || gx >= CGW || gy >= CGW) continue;
                float wx = (gx + 0.5f) * CCELL;
                float wy = (gy + 0.5f) * CCELL;
                float dx = wx - a[k].x, dy = wy - a[k].y;
                field[gx][gy] += expf(-(dx * dx + dy * dy) / s2);
            }
    }
}

static float field_sample(const float field[CGW][CGW], float x, float y) {
    float gx = x / CCELL - 0.5f, gy = y / CCELL - 0.5f;
    int x0 = (int)floorf(gx), y0 = (int)floorf(gy);
    float tx = gx - x0, ty = gy - y0;
    int x1 = x0 + 1, y1 = y0 + 1;
    if (x0 < 0) x0 = 0;
    if (y0 < 0) y0 = 0;
    if (x1 >= CGW) x1 = CGW - 1;
    if (y1 >= CGW) y1 = CGW - 1;
    float a = field[x0][y0] * (1 - tx) + field[x1][y0] * tx;
    float b = field[x0][y1] * (1 - tx) + field[x1][y1] * tx;
    return a * (1 - ty) + b * ty;
}

static int far_enough(const Run* r, float x, float y) {
    for (int i = 0; i < r->n; i++) {
        float dx = r->p[i].x - x, dy = r->p[i].y - y;
        if (dx * dx + dy * dy < R0 * R0) return 0;
    }
    return 1;
}

/* Union-find over the accepted points. */
static int uf_find(int* parent, int i) {
    while (parent[i] != i) { parent[i] = parent[parent[i]]; i = parent[i]; }
    return i;
}

static void uf_union(int* parent, int a, int b) {
    a = uf_find(parent, a); b = uf_find(parent, b);
    if (a != b) parent[b] = a;
}

static void analyze(Run* r) {
    int parent[NMAX];
    for (int i = 0; i < r->n; i++) parent[i] = i;
    float link2 = RLINK * RLINK;
    for (int i = 0; i < r->n; i++)
        for (int j = i + 1; j < r->n; j++) {
            float dx = r->p[i].x - r->p[j].x;
            float dy = r->p[i].y - r->p[j].y;
            if (dx * dx + dy * dy <= link2) uf_union(parent, i, j);
        }

    int root[NMAX], size[NMAX];
    memset(root, -1, sizeof(root));
    memset(size, 0, sizeof(size));
    for (int i = 0; i < r->n; i++) {
        int rep = uf_find(parent, i);
        if (root[rep] < 0) root[rep] = rep;
        size[rep]++;
    }

    r->veins = 0; r->singles = 0; r->mean_len = 0.0f; r->longest = 0.0f;
    for (int rep = 0; rep < r->n; rep++) {
        if (root[rep] != rep) continue;
        if (size[rep] == 1) { r->singles++; continue; }
        float diam = 0.0f;
        for (int i = 0; i < r->n; i++) {
            if (uf_find(parent, i) != rep) continue;
            for (int j = i + 1; j < r->n; j++) {
                if (uf_find(parent, j) != rep) continue;
                float dx = r->p[i].x - r->p[j].x;
                float dy = r->p[i].y - r->p[j].y;
                float d = sqrtf(dx * dx + dy * dy);
                if (d > diam) diam = d;
            }
        }
        r->veins++;
        r->mean_len += diam;
        if (diam > r->longest) r->longest = diam;
    }
    if (r->veins > 0) r->mean_len /= r->veins;
}

static void run_once(float sigma, unsigned seed, Run* r) {
    P anchors[ANCH * ANCH];
    make_anchors(anchors, seed);
    float field[CGW][CGW];
    make_field(field, anchors, ANCH * ANCH, sigma);

    memset(r, 0, sizeof(*r));
    int ncell = NCELL;
    int ncand = ncell * ncell;

    /* deterministic shuffled visitation: when the cap binds, thinning is
     * spread over every vein instead of filling the scan's first half. */
    int order[NCELL * NCELL];
    for (int k = 0; k < ncand; k++) order[k] = k;
    for (int k = ncand - 1; k > 0; k--) {
        int j = (int)(h2(k, (int)seed, 99u) % (unsigned)(k + 1));
        int tmp = order[k]; order[k] = order[j]; order[j] = tmp;
    }

    for (int q = 0; q < ncand; q++) {
        int k = order[q];
        int gx = k / ncell, gy = k % ncell;
        float x = (gx + u01(gx, gy, seed + 7u)) * CSTEP;
        float y = (gy + u01(gy, gx, seed + 13u)) * CSTEP;
        if (field_sample(field, x, y) < TAU) continue;
        if (!far_enough(r, x, y)) continue;
        if (r->n >= NMAX) { r->capped = 1; break; }
        r->p[r->n].x = x; r->p[r->n].y = y;
        r->n++;
    }
    analyze(r);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* --- self tests ---------------------------------------------------- */
    Run a, b;
    run_once(5.0f, 11u, &a);
    run_once(5.0f, 11u, &b);
    check("deterministic",
          a.n == b.n && a.veins == b.veins &&
          fabsf(a.p[0].x - b.p[0].x) < 1e-9f &&
          fabsf(a.longest - b.longest) < 1e-6f);

    Run small, big;
    run_once(3.0f, 11u, &small);
    run_once(6.0f, 11u, &big);
    check("count_grows_uncapped", big.n > small.n);
    check("length_grows_uncapped", big.longest > small.longest);
    check("veins_exist", small.veins > 0 && big.veins > 0);

    Run huge;
    run_once(40.0f, 11u, &huge);
    check("cap_enforced", huge.n == NMAX && huge.capped);

    int viol = 0;
    for (int i = 0; i < big.n; i++)
        for (int j = i + 1; j < big.n; j++) {
            float dx = big.p[i].x - big.p[j].x;
            float dy = big.p[i].y - big.p[j].y;
            if (dx * dx + dy * dy < R0 * R0 - 1e-5f) viol++;
        }
    check("min_distance_held", viol == 0);

    /* every accepted point must sit on thresholded vein ground */
    P chk_anchors[ANCH * ANCH];
    make_anchors(chk_anchors, 11u);
    float chk_field[CGW][CGW];
    make_field(chk_field, chk_anchors, ANCH * ANCH, 6.0f);
    int off_ground = 0;
    for (int i = 0; i < big.n; i++)
        if (field_sample(chk_field, big.p[i].x, big.p[i].y) < TAU - 1e-6f)
            off_ground++;
    check("points_on_vein_ground", off_ground == 0);

    /* --- calibrate the two curves on the uncapped pilot ---------------- */
    Run pilot;
    run_once(SIGMA0, 11u, &pilot);
    float cL = pilot.mean_len / SIGMA0;          /* L(s) = cL * s        */
    float cN = (float)pilot.n / (SIGMA0 * SIGMA0);/* N(s) = cN * s^2     */
    float s_cap = sqrtf((float)NMAX / cN);        /* scale where cap binds */
    float L_cap = cL * s_cap;

    printf("\nformula_length: L(s) = %.4f * s\n", cL);
    printf("formula_count : N(s) = min(%d, %.4f * s^2)\n", NMAX, cN);
    printf("cap_point     : s_cap=%.3f  L_cap=%.3f\n\n", s_cap, L_cap);

    /* --- three target lengths ------------------------------------------ */
    const float targets[3] = { 6.0f, 13.0f, 30.0f };
    Run tier[3];
    float sigma_star[3], n_req[3];

    printf(" tier  target_L  sigma_star  N_required  N_actual |"
           " vein_count  mean_length  longest_vein\n");
    for (int t = 0; t < 3; t++) {
        sigma_star[t] = targets[t] / cL;
        n_req[t] = cN * sigma_star[t] * sigma_star[t];
        run_once(sigma_star[t], 11u, &tier[t]);
        printf("  %d     %7.2f  %9.3f  %10.2f  %8d |  %9d  %11.3f  %12.3f\n",
               t + 1, targets[t], sigma_star[t], n_req[t], tier[t].n,
               tier[t].veins, tier[t].mean_len, tier[t].longest);
    }

    printf("\nnine_numbers (count, mean_length, longest per tier):\n");
    for (int t = 0; t < 3; t++)
        printf("  tier%d: count=%d mean_length=%.3f longest=%.3f\n",
               t + 1, tier[t].veins, tier[t].mean_len, tier[t].longest);

    int inverted = tier[2].longest < tier[0].longest;
    printf("\nlongest_tier_inverted=%d (tier3=%.3f vs tier1=%.3f)\n",
           inverted, tier[2].longest, tier[0].longest);
    if (inverted) {
        /* cap diagnosis: actual smoothing scale next to actual sample count */
        printf("cap_diagnosis_tier3: actual_sigma=%.3f actual_samples=%d"
               " (requested_N=%.1f, s_cap=%.3f, singles=%d, capped=%d)\n",
               sigma_star[2], tier[2].n, n_req[2], s_cap,
               tier[2].singles, tier[2].capped);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
