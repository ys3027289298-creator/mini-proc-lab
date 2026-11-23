#include <math.h>
#include <stdio.h>
#include <string.h>
#include "scatter.h"

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* Synthetic slope field: slope steepens with x, so acceptance should fall off
 * monotonically along x. */
static float half_density(void* ctx, float slope) {
    (void)ctx; (void)slope;
    return 0.5f;
}

static float ramp_slope(void* ctx, float x, float y) {
    (void)ctx; (void)y;
    return x / 4.0f; /* 0 .. 3 over the 12-unit-wide grid */
}

int main(void) {
    Set s;
    set_init(&s);
    int rejected = jitter_grid_flat(&s, 4, 4, 1.0f, 0.0f, 1u);
    check("grid_count", s.n == 16 && rejected == 0);

    Set t;
    set_init(&t);
    jitter_grid_flat(&t, 4, 4, 0.5f, 0.6f, 2u);
    int violations = 0;
    for (int i = 0; i < t.n; i++)
        for (int j = i + 1; j < t.n; j++) {
            float dx = t.p[i].x - t.p[j].x, dy = t.p[i].y - t.p[j].y;
            if (dx * dx + dy * dy < 0.36f - 1e-5f) violations++;
        }
    check("min_distance_enforced", t.n < 64);
    check("no_violation", violations == 0);

    check("deterministic", ({ Set a, b; set_init(&a); set_init(&b);
        jitter_grid_flat(&a, 3, 3, 1.0f, 0.0f, 9u); jitter_grid_flat(&b, 3, 3, 1.0f, 0.0f, 9u);
        a.n == b.n && fabsf(a.p[0].x - b.p[0].x) < 1e-9f; }));

    /* Density curve: monotone decreasing, in [0, 1]. */
    SlopeDensity d = { 1.0f, 0.8f, 1.3f, 0.0f };
    check("dens_flat_is_rho0", fabsf(slope_density_eval(&d, 0.0f) - 1.0f) < 1e-6f);
    int mono = 1;
    for (int i = 0; i < 100; i++) {
        float a = slope_density_eval(&d, i * 0.05f);
        float b = slope_density_eval(&d, (i + 1) * 0.05f);
        if (!(b <= a + 1e-6f) || a < 0.0f || a > 1.0f) mono = 0;
    }
    check("dens_monotone_bounded", mono);

    SlopeDensity cut = { 1.0f, 10.0f, 1.0f, 0.5f };
    check("dens_hard_cutoff", slope_density_eval(&cut, 0.5f) == 0.0f);

    /* Steep terrain must shed candidates relative to flat. */
    Set flat, steep;
    set_init(&flat);
    set_init(&steep);
    jitter_grid(&flat, 24, 24, 0.5f, 0.0f, 5u, NULL, NULL, NULL);
    jitter_grid(&steep, 24, 24, 0.5f, 0.0f, 5u, &d, ramp_slope, NULL);
    check("steep_fewer", steep.n < flat.n);
    check("flat_unchanged", flat.n == 576);

    /* Regression: candidate keys used to collide (i*131 + j*17 + seed*7919
     * has only ~75k distinct values), so a 32x32 grid reused one jitter pair
     * and one acceptance draw across whole runs of candidates. With mind = 0
     * every candidate must still occupy a distinct position. */
    {
        Set u;
        set_init(&u);
        jitter_grid_flat(&u, 32, 32, 0.5f, 0.0f, 3u);
        int dup = 0;
        for (int i = 0; i < u.n && !dup; i++)
            for (int j = i + 1; j < u.n; j++)
                if (u.p[i].x == u.p[j].x && u.p[i].y == u.p[j].y) { dup = 1; break; }
        check("positions_distinct", u.n == 1024 && dup == 0);
    }

    /* Constant density 0.5 over 32x32 candidates: the realised fraction must
     * sit within a few binomial sigma of 0.5. Key collisions made it drift. */
    {
        static int exposure[1], retained[1];
        static double mean_slope = 0.0;
        scatter_measure(32, 32, 0.5f, 3u, NULL, NULL, half_density, NULL, 0.0f,
                        1.0f, 1, exposure, retained, &mean_slope, NULL);
        double p = exposure[0] ? (double)retained[0] / exposure[0] : 0.0;
        double sigma = sqrt(0.25 / (exposure[0] ? exposure[0] : 1));
        check("constant_density_half", exposure[0] == 1024 &&
              fabs(p - 0.5) < 4.0 * sigma);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
