
/* bluenoise: periodic 2D point samplers on the unit torus.
 *
 *   poisson N seed            uniform iid (white noise reference)
 *   dart r Nmax seed [pat]    stationary dart throwing, hard minimum r;
 *                             Nmax=0 runs to saturation (pat rejections)
 *   bridson r Nmax seed [k]   Bridson fast Poisson-disk (saturation probe)
 *   grid n                    n*n square lattice (spectral reference)
 *   selftest                  built-in checks
 *
 * Output: one '# meta' header then 'x,y' rows, coordinates in [0,1).
 * All distances are periodic (wrap at the unit square), so downstream FFT
 * and nearest-neighbour measurements need no windowing or edge correction.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef struct { double x, y; } Pt;

/* deterministic RNG: splitmix64 */
static unsigned long long sm_state;
static unsigned long long splitmix64(void) {
    unsigned long long z = (sm_state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
static double rand01(void) { return (splitmix64() >> 11) * (1.0 / 9007199254740992.0); }

static double wrap1(double v) { return v - floor(v); }

/* periodic squared distance on the unit torus, wrapped to [-0.5, 0.5) */
static double pdist2(double x1, double y1, double x2, double y2) {
    double dx = x1 - x2, dy = y1 - y2;
    dx -= floor(dx + 0.5);
    dy -= floor(dy + 0.5);
    return dx * dx + dy * dy;
}

typedef struct {
    Pt *p;          /* accepted points */
    int n, cap;
    int *active;    /* active-list indices, swap-removal */
    int na;
    int *cell;      /* spatial hash: point index or -1; <=1 point per cell */
    int m;          /* cells per axis, cell side = 1/m <= r/sqrt(2) */
    double r, r2;
} Disk;

static int disk_init(Disk *d, double r) {
    d->r = r; d->r2 = r * r;
    d->m = (int)ceil(sqrt(2.0) / r);
    if (d->m < 2) d->m = 2;
    d->cap = 4096;
    d->p = malloc((size_t)d->cap * sizeof(Pt));
    d->active = malloc((size_t)d->cap * sizeof(int));
    d->cell = malloc((size_t)d->m * d->m * sizeof(int));
    if (!d->p || !d->active || !d->cell) return -1;
    for (int i = 0; i < d->m * d->m; i++) d->cell[i] = -1;
    d->n = 0; d->na = 0;
    return 0;
}

static void disk_free(Disk *d) {
    free(d->p); free(d->active); free(d->cell);
}

static int disk_occupied(Disk *d, double x, double y) {
    int gx = (int)(x * d->m) % d->m;
    int gy = (int)(y * d->m) % d->m;
    /* cell side s=1/m <= r/sqrt(2); a neighbour within r can then sit at
     * most 2 cells away on either axis, so a 5x5 window is sufficient. */
    for (int ix = -2; ix <= 2; ix++) {
        int cx = (gx + ix + d->m) % d->m;
        for (int iy = -2; iy <= 2; iy++) {
            int cy = (gy + iy + d->m) % d->m;
            int q = d->cell[cy * d->m + cx];
            if (q >= 0 && pdist2(x, y, d->p[q].x, d->p[q].y) < d->r2)
                return 1;
        }
    }
    return 0;
}

static int disk_add(Disk *d, double x, double y) {
    if (d->n >= d->cap) {
        d->cap *= 2;
        d->p = realloc(d->p, (size_t)d->cap * sizeof(Pt));
        d->active = realloc(d->active, (size_t)d->cap * sizeof(int));
        if (!d->p || !d->active) return -1;
    }
    int gx = (int)(x * d->m) % d->m;
    int gy = (int)(y * d->m) % d->m;
    int idx = d->n;
    d->p[idx].x = x; d->p[idx].y = y;
    d->cell[gy * d->m + gx] = idx;
    d->active[d->na++] = idx;
    d->n++;
    return 0;
}

/* Periodic dart throwing: candidates are uniform on the whole torus, so the
 * process is stationary at every accepted point (no growing front). Stops at
 * nmax accepted points, or after `patience` consecutive rejections when
 * nmax == 0 (saturation). This is the sampler used for fixed-N blue noise. */
static int dart_throw(double r, int nmax, unsigned long long seed,
                      long long patience, int emit) {
    sm_state = seed;
    Disk d;
    if (disk_init(&d, r)) return 1;
    long long fails = 0, darts = 0;
    while ((nmax == 0 || d.n < nmax) && fails < patience) {
        double x = rand01(), y = rand01();
        darts++;
        if (disk_occupied(&d, x, y)) { fails++; continue; }
        if (disk_add(&d, x, y)) { disk_free(&d); return 1; }
        fails = 0;
    }
    if (emit) {
        printf("# sampler=dart r=%.10f nmax=%d n=%d seed=%llu darts=%lld\n",
               r, nmax, d.n, seed, darts);
        for (int i = 0; i < d.n; i++)
            printf("%.10f,%.10f\n", d.p[i].x, d.p[i].y);
    }
    disk_free(&d);
    return 0;
}

/* run dart throw silently into caller-owned arrays for selftest */
static int run_dart(double r, int nmax, unsigned long long seed,
                    long long patience, Pt *out, int *out_n, int cap) {
    sm_state = seed;
    Disk d;
    if (disk_init(&d, r)) return -1;
    long long fails = 0;
    while ((nmax == 0 || d.n < nmax) && fails < patience) {
        double x = rand01(), y = rand01();
        if (disk_occupied(&d, x, y)) { fails++; continue; }
        if (disk_add(&d, x, y)) { disk_free(&d); return -1; }
        fails = 0;
    }
    if (d.n > cap) { disk_free(&d); return -1; }
    memcpy(out, d.p, (size_t)d.n * sizeof(Pt));
    *out_n = d.n;
    disk_free(&d);
    return 0;
}

/* Bridson dart throwing. Stops at Nmax accepted points, or at saturation
 * (empty active list) when Nmax == 0. k = darts per active point.
 * Fast near jamming, but stop-at-N output inherits a growing front, so fixed
 * N blue-noise sets should use dart_throw instead. */
static int bridson(double r, int nmax, unsigned long long seed, int k) {
    sm_state = seed;
    Disk d;
    if (disk_init(&d, r)) return 1;
    if (disk_add(&d, rand01(), rand01())) { disk_free(&d); return 1; }

    while (d.na > 0 && (nmax == 0 || d.n < nmax)) {
        int ai = (int)(splitmix64() % (unsigned long long)d.na);
        Pt base = d.p[d.active[ai]];
        int accepted = 0;
        for (int t = 0; t < k && (nmax == 0 || d.n < nmax); t++) {
            double ang = 6.283185307179586 * rand01();
            double rad = r * sqrt(1.0 + 3.0 * rand01()); /* uniform in [r, 2r] */
            double x = wrap1(base.x + rad * cos(ang));
            double y = wrap1(base.y + rad * sin(ang));
            if (!disk_occupied(&d, x, y)) {
                if (disk_add(&d, x, y)) { disk_free(&d); return 1; }
                accepted = 1;
                break;
            }
        }
        if (!accepted) {
            d.active[ai] = d.active[--d.na]; /* swap-remove */
        }
    }

    printf("# sampler=bridson r=%.10f nmax=%d n=%d seed=%llu k=%d m=%d\n",
           r, nmax, d.n, seed, k, d.m);
    for (int i = 0; i < d.n; i++)
        printf("%.10f,%.10f\n", d.p[i].x, d.p[i].y);
    disk_free(&d);
    return 0;
}

static int poisson(int n, unsigned long long seed) {
    sm_state = seed;
    printf("# sampler=poisson r=0 nmax=%d n=%d seed=%llu\n", n, n, seed);
    for (int i = 0; i < n; i++)
        printf("%.10f,%.10f\n", rand01(), rand01());
    return 0;
}

static int square_grid(int n) {
    printf("# sampler=grid r=%.10f nmax=%d n=%d seed=0\n",
           1.0 / n, n * n, n * n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            printf("%.10f,%.10f\n", (i + 0.5) / n, (j + 0.5) / n);
    return 0;
}

static int min_violations(const Pt *p, int n, double r) {
    int bad = 0;
    double r2 = r * r;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (pdist2(p[i].x, p[i].y, p[j].x, p[j].y) < r2 - 1e-12) bad++;
    return bad;
}

static int failures = 0;
static void check(const char *name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* run bridson silently into caller-owned arrays for selftest */
static int run_bridson(double r, int nmax, unsigned long long seed, int k,
                       Pt *out, int *out_n, int cap) {
    sm_state = seed;
    Disk d;
    if (disk_init(&d, r)) return -1;
    if (disk_add(&d, rand01(), rand01())) { disk_free(&d); return -1; }
    while (d.na > 0 && (nmax == 0 || d.n < nmax)) {
        int ai = (int)(splitmix64() % (unsigned long long)d.na);
        Pt base = d.p[d.active[ai]];
        int accepted = 0;
        for (int t = 0; t < k && (nmax == 0 || d.n < nmax); t++) {
            double ang = 6.283185307179586 * rand01();
            double rad = r * sqrt(1.0 + 3.0 * rand01());
            double x = wrap1(base.x + rad * cos(ang));
            double y = wrap1(base.y + rad * sin(ang));
            if (!disk_occupied(&d, x, y)) {
                if (disk_add(&d, x, y)) { disk_free(&d); return -1; }
                accepted = 1;
                break;
            }
        }
        if (!accepted) d.active[ai] = d.active[--d.na];
    }
    if (d.n > cap) { disk_free(&d); return -1; }
    memcpy(out, d.p, (size_t)d.n * sizeof(Pt));
    *out_n = d.n;
    disk_free(&d);
    return 0;
}

static int selftest(void) {
    /* periodic distance wraps both axes */
    check("pdist_wrap_x", fabs(sqrt(pdist2(0.05, 0.5, 0.95, 0.5)) - 0.1) < 1e-9);
    check("pdist_wrap_y", fabs(sqrt(pdist2(0.5, 0.02, 0.5, 0.98)) - 0.04) < 1e-9);
    check("pdist_direct", fabs(sqrt(pdist2(0.2, 0.2, 0.5, 0.6)) - 0.5) < 1e-9);

    enum { CAP = 20000 };
    Pt *buf = malloc((size_t)CAP * sizeof(Pt));
    int n1, n2, n3;

    /* saturated disk: hard minimum distance, deterministic, packing bounds.
     * Centre density bounds for exclusion diameter r on the unit torus:
     * trivial square bound rho < 1/r^2; densest (hex lattice) bound
     * rho_max = 2/(sqrt(3) r^2) ~= 1.1547/r^2 ~= 461.9 for r=0.05. */
    run_bridson(0.05, 0, 42u, 30, buf, &n1, CAP);
    check("sat_count_range", n1 > 200 && n1 < 462); /* >1/(2r^2)=200 loose lower bound, hex ~462 upper */
    check("sat_no_violation", min_violations(buf, n1, 0.05) == 0);

    Pt *buf2 = malloc((size_t)CAP * sizeof(Pt));
    run_bridson(0.05, 0, 42u, 30, buf2, &n2, CAP);
    check("deterministic", n1 == n2 &&
          fabs(buf[0].x - buf2[0].x) < 1e-12 && fabs(buf[1].y - buf2[1].y) < 1e-12);
    free(buf2);

    /* capped run also enforces the disk and stops at the cap */
    run_bridson(0.05, 120, 7u, 30, buf, &n3, CAP);
    check("cap_stops_at_nmax", n3 == 120);
    check("cap_no_violation", min_violations(buf, n3, 0.05) == 0);

    /* smaller disk admits strictly more points at saturation */
    int n_small;
    run_bridson(0.03, 0, 42u, 30, buf, &n_small, CAP);
    check("smaller_r_more_points", n_small > n1);
    check("small_no_violation", min_violations(buf, n_small, 0.03) == 0);

    /* loose disk stop-at-N below the cap works (dart thrower, stationary) */
    int n_loose;
    run_dart(0.0078125, 300, 11u, 1000000, buf, &n_loose, CAP);
    check("loose_cap_count", n_loose == 300);
    check("loose_no_violation", min_violations(buf, n_loose, 0.0078125) == 0);

    /* dart throwing with a loose disk must fill the whole torus: the
     * accepted points hit every 4x4 quadrant (no Bridson-style front blob) */
    run_dart(0.0078125, 300, 11u, 1000000, buf, &n_loose, CAP);
    int quad[4][4] = {{0}};
    for (int i = 0; i < n_loose; i++)
        quad[(int)(buf[i].y * 4)][(int)(buf[i].x * 4)]++;
    int minquad = quad[0][0];
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            if (quad[i][j] < minquad) minquad = quad[i][j];
    check("dart_covers_domain", minquad > 0);
    check("dart_deterministic", ({
        Pt b2[CAP]; int nn2;
        run_dart(0.0078125, 300, 11u, 1000000, b2, &nn2, CAP);
        nn2 == n_loose &&
        fabs(buf[0].x - b2[0].x) < 1e-12 && fabs(buf[5].y - b2[5].y) < 1e-12;
    }));

    /* square lattice: n^2 points with periodic spacing 1/n */
    Pt g[4];
    int gn = 0;
    int n = 2;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            g[gn].x = (i + 0.5) / n; g[gn].y = (j + 0.5) / n; gn++;
        }
    check("grid_count", gn == 4);
    check("grid_spacing", min_violations(g, gn, 0.5 - 1e-9) == 0 &&
          min_violations(g, gn, 0.5 + 1e-9) > 0);

    free(buf);
    printf("all_ok=%d\n", failures == 0);
    return failures ? 1 : 0;
}

int main(int argc, char **argv) {
    if (argc >= 2 && strcmp(argv[1], "selftest") == 0)
        return selftest();

    if (argc >= 3 && strcmp(argv[1], "poisson") == 0) {
        int n = atoi(argv[2]);
        unsigned long long seed = argc > 3 ? strtoull(argv[3], NULL, 10) : 1;
        if (n <= 0) return 1;
        return poisson(n, seed);
    }

    if (argc >= 4 && strcmp(argv[1], "dart") == 0) {
        double r = atof(argv[2]);
        int nmax = atoi(argv[3]);
        unsigned long long seed = argc > 4 ? strtoull(argv[4], NULL, 10) : 1;
        long long patience = argc > 5 ? atoll(argv[5]) : 1000000;
        if (!(r > 0 && r < 0.5) || nmax < 0 || patience <= 0) return 1;
        return dart_throw(r, nmax, seed, patience, 1);
    }

    if (argc >= 4 && strcmp(argv[1], "bridson") == 0) {
        double r = atof(argv[2]);
        int nmax = atoi(argv[3]);
        unsigned long long seed = argc > 4 ? strtoull(argv[4], NULL, 10) : 1;
        int k = argc > 5 ? atoi(argv[5]) : 30;
        if (!(r > 0 && r < 0.5) || nmax < 0 || k <= 0) return 1;
        return bridson(r, nmax, seed, k);
    }

    if (argc >= 3 && strcmp(argv[1], "grid") == 0) {
        int n = atoi(argv[2]);
        if (n <= 0) return 1;
        return square_grid(n);
    }

    fprintf(stderr,
        "usage: %s poisson N seed | dart r Nmax seed [patience] |\n"
        "       bridson r Nmax seed [k] | grid n | selftest\n",
        argv[0]);
    return 1;
}
