
/* thermal: point heat anomaly expanded into a finite Gaussian kernel.
 *
 * Diffusion basis (2D heat equation, instantaneous point source of energy E):
 *   T(r,t) = E/(4 pi D t) * exp(-r^2/(4 D t))   =>  Gaussian, sigma^2 = 2 D t.
 * Patch region = isotherm T >= theta (detection threshold):
 *   R = sigma * sqrt(2 ln(E/(2 pi sigma^2 theta))),   area = pi R^2.
 * The kernel is truncated at R (finite support); truncation drops the fraction
 * exp(-R^2/(2 sigma^2)) of the total energy, which is why R = k*sigma, k >= 2.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define N 160   /* display grid side, cell = 1 */
#define SS 16   /* subsamples per cell axis for area statistics */
#define PI 3.14159265358979323846

#define SIGMA 8.0   /* sqrt(2 D t), diffusion length in cells */
#define THETA 1.0   /* detection threshold (isotherm level) */
#define NTIER 3

typedef struct { int k; double R, peak, E; } Tier;

/* finite Gaussian kernel: compact support r <= R */
static double kern(double r2, double peak, double R) {
    if (r2 > R * R) return 0.0;
    return peak * exp(-r2 / (2.0 * SIGMA * SIGMA));
}

/* supersampled count of the T >= theta region */
static double measure_area(const Tier* t) {
    double cx = N / 2.0, cy = N / 2.0, area = 0.0;
    int ext = (int)ceil(t->R) + 1;
    for (int iy = (int)cy - ext; iy <= (int)cy + ext; iy++)
        for (int ix = (int)cx - ext; ix <= (int)cx + ext; ix++) {
            int hit = 0;
            for (int sy = 0; sy < SS; sy++)
                for (int sx = 0; sx < SS; sx++) {
                    double x = ix + (sx + 0.5) / SS - cx;
                    double y = iy + (sy + 0.5) / SS - cy;
                    if (kern(x * x + y * y, t->peak, t->R) >= THETA) hit++;
                }
            area += (double)hit / (SS * SS);
        }
    return area;
}

/* total energy of the finite kernel on the unit-cell grid */
static double grid_energy(const Tier* t) {
    double cx = N / 2.0, cy = N / 2.0, s = 0.0;
    for (int iy = 0; iy < N; iy++)
        for (int ix = 0; ix < N; ix++) {
            double dx = ix + 0.5 - cx, dy = iy + 0.5 - cy;
            s += kern(dx * dx + dy * dy, t->peak, t->R);
        }
    return s; /* cell area = 1 */
}

/* recover sigma^2 from the second moment of an untruncated Gaussian:
 * <r^2> = 2 sigma^2, so sigma^2 = M2 / (2 M0) */
static double sigma2_from_moment(void) {
    double cx = N / 2.0, cy = N / 2.0, m0 = 0.0, m2 = 0.0;
    for (int iy = 0; iy < N; iy++)
        for (int ix = 0; ix < N; ix++) {
            double dx = ix + 0.5 - cx, dy = iy + 0.5 - cy;
            double g = exp(-(dx * dx + dy * dy) / (2.0 * SIGMA * SIGMA));
            m0 += g;
            m2 += (dx * dx + dy * dy) * g;
        }
    return m2 / (2.0 * m0);
}

/* comparison figure: rows = energy tiers, left = raw point, right = kernel
 * with the theoretical radius R drawn as a mid-gray ring */
static int write_figure(const Tier* tiers, const char* path) {
    const int P = N, G = 8, W = 2 * P + 3 * G, H = NTIER * P + (NTIER + 1) * G;
    static unsigned char buf[(2 * 160 + 3 * 8) * (3 * 160 + 4 * 8)];
    memset(buf, 24, (size_t)W * H);
    for (int i = 0; i < NTIER; i++) {
        const Tier* t = &tiers[i];
        int oy = G + i * (P + G);
        /* left panel: the unexpanded anomaly is a single pixel */
        buf[(oy + P / 2) * W + (G + P / 2)] = 255;
        /* right panel: kernel normalized to its own peak */
        int ox = 2 * G + P;
        for (int y = 0; y < P; y++)
            for (int x = 0; x < P; x++) {
                double dx = x + 0.5 - P / 2.0, dy = y + 0.5 - P / 2.0;
                double r2 = dx * dx + dy * dy, r = sqrt(r2);
                double v = kern(r2, t->peak, t->R) / t->peak;
                unsigned char c = (unsigned char)(v * 255.0 + 0.5);
                if (fabs(r - t->R) < 0.8) c = 160; /* theoretical isotherm */
                buf[(oy + y) * W + (ox + x)] = c;
            }
    }
    FILE* f = fopen(path, "wb");
    if (!f) return -1;
    fprintf(f, "P5\n%d %d\n255\n", W, H);
    fwrite(buf, 1, (size_t)W * H, f);
    fclose(f);
    return 0;
}

static void print_ascii(const Tier* t) {
    static const char ramp[] = " .:-=+*#%@";
    for (int row = 0; row < 32; row++) {
        for (int col = 0; col < 64; col++) {
            double x = (col + 0.5) * N / 64.0 - N / 2.0;
            double y = (row + 0.5) * N / 32.0 - N / 2.0;
            double v = kern(x * x + y * y, t->peak, t->R) / t->peak;
            putchar(ramp[(int)(v * 9.999)]);
        }
        putchar('\n');
    }
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static double rel(double got, double want) { return fabs(got - want) / want; }

int main(void) {
    static const int ks[NTIER] = {2, 3, 4};
    Tier tiers[NTIER];

    printf("diffusion: sigma=sqrt(2Dt)=%.1f cells, theta=%.1f, R=k*sigma\n",
           SIGMA, THETA);
    printf("tier k   E            peak      R_theory R_meas   A_theory  A_meas    err%%\n");
    for (int i = 0; i < NTIER; i++) {
        Tier* t = &tiers[i];
        t->k = ks[i];
        t->R = t->k * SIGMA;
        t->peak = THETA * exp(0.5 * t->k * t->k);       /* peak = E/(2 pi sigma^2) */
        t->E = 2.0 * PI * SIGMA * SIGMA * t->peak;      /* total energy */
        double a_meas = measure_area(t);
        double a_theory = PI * t->R * t->R;
        double r_meas = sqrt(a_meas / PI);
        printf("%d    %d   %.4e  %.4e  %.3f    %.3f   %.2f   %.2f   %.3f\n",
               i + 1, t->k, t->E, t->peak, t->R, r_meas, a_theory, a_meas,
               100.0 * rel(a_meas, a_theory));
        char nm[32];
        snprintf(nm, sizeof nm, "radius_tier%d", i + 1);
        check(nm, rel(r_meas, t->R) < 0.005);
        snprintf(nm, sizeof nm, "area_tier%d", i + 1);
        check(nm, rel(a_meas, a_theory) < 0.01);
        /* energy kept inside the finite kernel: E * (1 - exp(-k^2/2)) */
        snprintf(nm, sizeof nm, "energy_tier%d", i + 1);
        check(nm, rel(grid_energy(t), t->E * (1.0 - exp(-0.5 * t->k * t->k))) < 0.01);
    }

    /* sigma^2 = 2Dt recovered from the second moment of the kernel */
    check("sigma2_moment", rel(sigma2_from_moment(), SIGMA * SIGMA) < 0.005);

    /* motivation: the raw point anomaly covers one cell, far below any patch */
    check("point_invisible", 1.0 < 0.01 * PI * tiers[0].R * tiers[0].R);

    /* area grows with energy tier */
    check("area_monotonic", tiers[0].R < tiers[1].R && tiers[1].R < tiers[2].R);

    if (write_figure(tiers, "thermal_compare.pgm") == 0)
        printf("figure=thermal_compare.pgm\n");
    else
        check("figure_written", 0);

    printf("tier 2 (k=%d, R=%.0f) patch, peak-normalized:\n", tiers[1].k, tiers[1].R);
    print_ascii(&tiers[1]);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
