/* slope_grid: grid-convergence study for bilinear sampling + central-
 * difference slope estimation across a steep cliff.
 *
 * Terrain: z(x,y) = A*tanh((x-X0)/W) on [0,L]^2 (y-independent, so the
 * cliff is a vertical escarpment along x). Analytic slope magnitude is
 * s(x) = (A/W) sech^2((x-X0)/W).
 *
 * Pipeline mirrored from terrain.c: nodal heights -> central-difference
 * slope (field_slope) -> compared against the analytic slope. A "stair
 * pixel" is a slope-map pixel whose absolute error exceeds TOL*S_MAX,
 * S_MAX = A/W being the peak analytic slope.
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

static const double L   = 16.0; /* domain size                      */
static const double X0  = 8.0;  /* cliff center                     */
static const double A   = 4.0;  /* cliff half-height (total drop 8) */
static const double W   = 0.8;  /* cliff transition width           */
static const double TOL = 0.10; /* acceptable relative slope error  */

static double z_exact(double x) { return A * tanh((x - X0) / W); }
static double s_exact(double x) {
    double c = cosh((x - X0) / W);
    return (A / W) / (c * c);
}

/* central-difference slope on a row of nodal heights, clamped borders,
 * identical to field_slope() in terrain.c (dy term is identically 0). */
static double cd_slope(const double* z, int i, int n, double h) {
    int im = i > 0 ? i - 1 : i, ip = i < n - 1 ? i + 1 : i;
    return (z[ip] - z[im]) / ((ip - im) * h);
}

typedef struct {
    int n; double h;
    double slope_max;      /* max computed slope            */
    double band_mean;      /* mean computed slope, |x-X0|<2W */
    double err_max;        /* max |s_computed - s_exact|    */
    double err_rms;        /* rms error over all nodes      */
    double peak_rel;       /* |slope_max - S_MAX| / S_MAX   */
    int stair_cols;        /* bad columns in the row        */
    long stair_pixels;     /* bad pixels in the N x N map   */
} Stats;

/* stair pixel: local relative error above TOL where the true slope is
 * non-negligible (> 1% of S_MAX) -- i.e. a visibly wrong slope texel. */
static int is_stair(double err, double s_exact, double smax) {
    return s_exact > 0.01 * smax && err > TOL * s_exact;
}

static Stats run_grid(int n) {
    double h = L / (n - 1);
    double smax = A / W;
    double* z = malloc((size_t)n * sizeof *z);
    for (int i = 0; i < n; i++) z[i] = z_exact(i * h);

    Stats st = { .n = n, .h = h, .slope_max = 0, .band_mean = 0,
                 .err_max = 0, .err_rms = 0, .peak_rel = 0,
                 .stair_cols = 0, .stair_pixels = 0 };
    int band = 0;
    for (int i = 0; i < n; i++) {
        double x = i * h;
        double sc = fabs(cd_slope(z, i, n, h));
        double se = s_exact(x);
        double e = fabs(sc - se);
        if (sc > st.slope_max) st.slope_max = sc;
        if (e > st.err_max) st.err_max = e;
        st.err_rms += e * e;
        if (fabs(x - X0) < 2 * W) { st.band_mean += sc; band++; }
        if (is_stair(e, se, smax)) st.stair_cols++;
    }
    st.band_mean /= band > 0 ? band : 1;
    st.err_rms = sqrt(st.err_rms / n);
    st.peak_rel = fabs(st.slope_max - smax) / smax;
    st.stair_pixels = (long)st.stair_cols * n; /* stripes run full height */
    free(z);
    return st;
}

int main(void) {
    double smax = A / W;
    /* crest truncation error of the CD+bilinear pipeline:
     * eps_rel ~= h^2 |z'''| / (6 s_max) = h^2 / (3 W^2)            */
    double h_min = W * sqrt(3.0 * TOL);
    int n_min = (int)ceil(L / h_min) + 1;

    printf("== theory ==\n");
    printf("S_MAX=A/W=%.4f  z'''_crest=-2A/W^3=%.4f\n", smax, -2 * A / (W * W * W));
    printf("eps_rel(h) ~= h^2/(3 W^2) = %.4f h^2\n", 1.0 / (3 * W * W));
    printf("tol=%.2f -> h_min = W*sqrt(3*tol) = %.4f -> N_min = ceil(L/h_min)+1 = %d\n",
           TOL, h_min, n_min);

    printf("\n== three tiers ==\n");
    printf("%4s %8s %9s %9s %9s %9s %9s %7s %9s\n",
           "N", "h", "smax_cmp", "band_avg", "err_max", "err_rms",
           "peak_rel", "stair", "below_min");
    int tiers[] = {16, 32, 64};
    for (int t = 0; t < 3; t++) {
        Stats st = run_grid(tiers[t]);
        printf("%4d %8.4f %9.4f %9.4f %9.4f %9.4f %9.4f %7ld %9s\n",
               st.n, st.h, st.slope_max, st.band_mean, st.err_max,
               st.err_rms, st.peak_rel, st.stair_pixels,
               st.h > h_min ? "YES" : "no");
    }

    printf("\n== stair pixels below the minimum grid (N < %d) ==\n", n_min);
    for (int t = 0; t < 3; t++) {
        Stats st = run_grid(tiers[t]);
        if (st.h > h_min)
            printf("N=%d (h=%.4f): %d bad columns x %d rows = %ld stair pixels\n",
                   st.n, st.h, st.stair_cols, st.n, st.stair_pixels);
    }

    printf("\n== convergence (max error vs h) ==\n");
    printf("%4s %8s %11s %7s\n", "N", "h", "err_max", "order");
    double ph = 0, pe = 0;
    for (int n = 16; n <= 256; n *= 2) {
        Stats st = run_grid(n);
        double order = ph > 0 ? log(pe / st.err_max) / log(ph / st.h) : 0;
        printf("%4d %8.4f %11.5f %7.2f\n", n, st.h, st.err_max, order);
        ph = st.h; pe = st.err_max;
    }

    printf("\n== crest slope profile, N=16 (stair plateaus) ==\n");
    {
        int n = 16; double h = L / (n - 1);
        double* z = malloc((size_t)n * sizeof *z);
        for (int i = 0; i < n; i++) z[i] = z_exact(i * h);
        for (int i = 0; i < n; i++)
            printf("x=%5.2f  computed=%7.4f  exact=%7.4f\n",
                   i * h, fabs(cd_slope(z, i, n, h)), s_exact(i * h));
        free(z);
    }
    return 0;
}
