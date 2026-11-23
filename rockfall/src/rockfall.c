
/* rockfall: slope-aware stone scatter.
 *
 * A stone of mass m = 1 kg is released on a slope of angle theta with an
 * initial kinetic energy E0 drawn from the collapse model E0 ~ Gamma(2, beta).
 * With effective friction mu (rolling included) and mu > tan(theta) the stone
 * decelerates at a = g (mu cos theta - sin theta) and stops after
 *
 *     S = E0 / (m g (mu cos theta - sin theta))
 *
 * hence S ~ Gamma(2, lambda) with lambda = beta / (m g (mu cos theta - sin theta)),
 * distribution peak (mode) s* = lambda, mean 2 lambda.
 *
 * The program simulates three (theta, mu, beta) groups, reports endpoint
 * statistics, writes SVG top views, and fails if the histogram peak deviates
 * from the derived mode by more than 30%.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979f
#define GRAV 9.8f
#define MASS 1.0f

#define GRID_N 32               /* start zone: GRID_N x GRID_N jittered grid */
#define NSTONE (GRID_N * GRID_N)
#define CELL 0.5f               /* 32 x 32 cells -> same 16 m start zone */
#define PHI_MAX (15.0f * PI / 180.0f) /* lateral spread around fall line */

#define MAXBINS 256

typedef struct { float x, y; } P;

typedef struct {
    const char* name;
    float theta_deg;  /* slope angle */
    float mu;         /* effective friction coefficient */
    float beta;       /* energy scale of Gamma(2, beta), joules */
    unsigned seed;
    const char* color;
} Group;

typedef struct {
    P start[NSTONE];
    P end[NSTONE];
    float s[NSTONE];  /* roll distance per stone */
    int n;
    float mean, std;
    float mode_hist;  /* smoothed-histogram peak (reference) */
    float mode_kde;   /* Gaussian KDE peak: the reported distribution peak */
    float mode_derived;
    float max_s;
} Sim;

static unsigned h1(unsigned i) {
    unsigned h = i * 2654435761u;
    h ^= h >> 15;
    h *= 2246822519u;
    h ^= h >> 13;
    return h;
}

/* uniform in (0, 1] so logf() is safe */
static float u01(unsigned i) {
    float u = (float)(h1(i) & 0xffff) / 65535.0f;
    return u < 1e-6f ? 1e-6f : u;
}

static float mode_derived(const Group* g) {
    float th = g->theta_deg * PI / 180.0f;
    float decel = GRAV * (g->mu * cosf(th) - sinf(th)); /* m = 1 kg */
    return g->beta / (MASS * decel);
}

static void sim_run(const Group* g, Sim* r) {
    float th = g->theta_deg * PI / 180.0f;
    float decel = GRAV * (g->mu * cosf(th) - sinf(th));
    r->n = 0;
    r->max_s = 0.0f;
    for (int i = 0; i < GRID_N; i++) {
        for (int j = 0; j < GRID_N; j++) {
            int k = r->n;
            unsigned base = g->seed * 7919u + (unsigned)(i * 131 + j * 17);
            float sx = (i + u01(base + 0)) * CELL;
            float sy = (j + u01(base + 1)) * CELL;
            /* E0 ~ Gamma(2, beta): sum of two exponentials */
            float e0 = -g->beta * logf(u01(base + 2) * u01(base + 3));
            float s = e0 / (MASS * decel);
            float phi = (u01(base + 4) - 0.5f) * 2.0f * PHI_MAX;
            r->start[k].x = sx; r->start[k].y = sy;
            r->end[k].x = sx + s * cosf(phi);
            r->end[k].y = sy + s * sinf(phi);
            r->s[k] = s;
            if (s > r->max_s) r->max_s = s;
            r->n = ++k;
        }
    }
    float mean = 0.0f;
    for (int i = 0; i < r->n; i++) mean += r->s[i];
    mean /= r->n;
    float var = 0.0f;
    for (int i = 0; i < r->n; i++) {
        float d = r->s[i] - mean;
        var += d * d;
    }
    r->mean = mean;
    r->std = sqrtf(var / r->n);
    r->mode_derived = mode_derived(g);

    /* histogram with bin width = derived mode / 3, 3-bin moving average */
    float w = r->mode_derived / 3.0f;
    int nbins = (int)(r->max_s / w) + 2;
    if (nbins > MAXBINS) nbins = MAXBINS;
    int cnt[MAXBINS];
    memset(cnt, 0, sizeof(cnt));
    for (int i = 0; i < r->n; i++) {
        int b = (int)(r->s[i] / w);
        if (b >= nbins) b = nbins - 1;
        cnt[b]++;
    }
    int best = 0;
    float bestv = -1.0f;
    for (int b = 0; b < nbins; b++) {
        float c = (float)cnt[b] * 2.0f;
        if (b > 0) c += (float)cnt[b - 1];
        if (b < nbins - 1) c += (float)cnt[b + 1];
        if (c > bestv) { bestv = c; best = b; }
    }
    r->mode_hist = (best + 0.5f) * w;

    /* distribution peak: argmax of Gaussian KDE, Silverman bandwidth */
    float h = 0.9f * r->std * powf((float)r->n, -0.2f);
    float bestx = 0.0f;
    bestv = -1.0f;
    for (int k = 0; k <= 400; k++) {
        float x = 0.05f * r->mode_derived + 4.0f * r->mode_derived * k / 400.0f;
        float v = 0.0f;
        for (int i = 0; i < r->n; i++) {
            float z = (x - r->s[i]) / h;
            if (z > 6.0f || z < -6.0f) continue;
            v += expf(-0.5f * z * z);
        }
        if (v > bestv) { bestv = v; bestx = x; }
    }
    r->mode_kde = bestx;
}

static int write_svg(const char* path, const Group* g, const Sim* r) {
    float x0 = 1e30f, y0 = 1e30f, x1 = -1e30f, y1 = -1e30f;
    for (int i = 0; i < r->n; i++) {
        if (r->start[i].x < x0) x0 = r->start[i].x;
        if (r->start[i].y < y0) y0 = r->start[i].y;
        if (r->end[i].x > x1) x1 = r->end[i].x;
        if (r->end[i].y > y1) y1 = r->end[i].y;
        if (r->end[i].y < y0) y0 = r->end[i].y;
        if (r->end[i].x > x1) x1 = r->end[i].x;
    }
    x0 -= 1.0f; y0 -= 1.0f; x1 += 1.0f; y1 += 1.0f;
    FILE* f = fopen(path, "w");
    if (!f) return -1;
    fprintf(f, "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"%.2f %.2f %.2f %.2f\""
               " width=\"900\" style=\"background:#fafaf7\">\n", x0, y0, x1 - x0, y1 - y0);
    fprintf(f, "<rect x=\"0\" y=\"0\" width=\"%.1f\" height=\"%.1f\" fill=\"#e8e2d4\""
               " stroke=\"#999\" stroke-width=\"0.06\"/>\n", GRID_N * CELL, GRID_N * CELL);
    /* predicted peak line: mean start x + derived mode */
    float mx = 8.0f + r->mode_derived;
    fprintf(f, "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"#d00\""
               " stroke-width=\"0.08\" stroke-dasharray=\"0.4 0.3\"/>\n", mx, y0, mx, y1);
    fprintf(f, "<text x=\"%.2f\" y=\"%.2f\" font-size=\"0.9\" fill=\"#d00\">mode=%.2fm</text>\n",
            mx + 0.2f, y0 + 1.0f, r->mode_derived);
    for (int i = 0; i < r->n; i++)
        fprintf(f, "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"%s\""
                   " stroke-width=\"0.05\" opacity=\"0.25\"/>\n",
                r->start[i].x, r->start[i].y, r->end[i].x, r->end[i].y, g->color);
    for (int i = 0; i < r->n; i++)
        fprintf(f, "<circle cx=\"%.2f\" cy=\"%.2f\" r=\"0.14\" fill=\"#555\"/>\n",
                r->start[i].x, r->start[i].y);
    for (int i = 0; i < r->n; i++)
        fprintf(f, "<circle cx=\"%.2f\" cy=\"%.2f\" r=\"0.18\" fill=\"%s\"/>\n",
                r->end[i].x, r->end[i].y, g->color);
    fprintf(f, "<text x=\"%.2f\" y=\"%.2f\" font-size=\"1.0\" fill=\"#333\">"
               "theta=%.0f deg, mu=%.2f, beta=%.1f J, n=%d</text>\n",
            x0 + 0.3f, y1 - 0.6f, g->theta_deg, g->mu, g->beta, r->n);
    fprintf(f, "</svg>\n");
    fclose(f);
    return 0;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Group groups[3] = {
        { "A", 20.0f, 0.50f,  8.0f, 11u, "#1f6fb2" },
        { "B", 30.0f, 0.65f,  8.0f, 23u, "#c2551f" },
        { "C", 15.0f, 0.35f, 12.0f, 37u, "#3d8a3d" },
    };
    Sim r[3];
    FILE* stats = fopen("out/stats.txt", "w");

    for (int gi = 0; gi < 3; gi++) {
        Group* g = &groups[gi];
        float th = g->theta_deg * PI / 180.0f;
        check(g->name, g->mu > tanf(th)); /* must decelerate on the slope */
        sim_run(g, &r[gi]);

        float dev = fabsf(r[gi].mode_kde - r[gi].mode_derived) / r[gi].mode_derived;
        float dev2 = fabsf(r[gi].mean / 2.0f - r[gi].mode_derived) / r[gi].mode_derived;
        char key[64];
        printf("group %s: theta=%.0f mu=%.2f beta=%.1f n=%d mean=%.3f std=%.3f "
               "mode_kde=%.3f mode_hist=%.3f mode_mean/2=%.3f mode_derived=%.3f "
               "dev_kde=%.1f%% dev_mean/2=%.1f%%\n",
               g->name, g->theta_deg, g->mu, g->beta, r[gi].n, r[gi].mean, r[gi].std,
               r[gi].mode_kde, r[gi].mode_hist, r[gi].mean / 2.0f, r[gi].mode_derived,
               dev * 100.0f, dev2 * 100.0f);
        if (stats)
            fprintf(stats, "group %s: theta=%.0f mu=%.2f beta=%.1f n=%d mean=%.3f std=%.3f "
                           "mode_kde=%.3f mode_hist=%.3f mode_mean/2=%.3f mode_derived=%.3f "
                           "dev_kde=%.1f%% dev_mean/2=%.1f%%\n",
                    g->name, g->theta_deg, g->mu, g->beta, r[gi].n, r[gi].mean, r[gi].std,
                    r[gi].mode_kde, r[gi].mode_hist, r[gi].mean / 2.0f, r[gi].mode_derived,
                    dev * 100.0f, dev2 * 100.0f);

        snprintf(key, sizeof(key), "peak_within_30pct_%s", g->name);
        check(key, dev <= 0.30f);
        snprintf(key, sizeof(key), "peak_mean2_within_30pct_%s", g->name);
        check(key, dev2 <= 0.30f);
        snprintf(key, sizeof(key), "mean_near_2x_mode_%s", g->name);
        check(key, fabsf(r[gi].mean - 2.0f * r[gi].mode_derived) /
                   (2.0f * r[gi].mode_derived) <= 0.15f);

        char path[64];
        snprintf(path, sizeof(path), "out/rockfall_%s.svg", g->name);
        snprintf(key, sizeof(key), "svg_written_%s", g->name);
        check(key, write_svg(path, g, &r[gi]) == 0);
    }
    if (stats) fclose(stats);

    /* determinism: same seed, same endpoints */
    Sim a, b;
    sim_run(&groups[0], &a);
    sim_run(&groups[0], &b);
    check("deterministic", a.n == b.n &&
          fabsf(a.end[0].x - b.end[0].x) < 1e-9f &&
          fabsf(a.end[123].y - b.end[123].y) < 1e-9f);

    /* physics direction: steeper effective deceleration -> shorter runout */
    check("mode_ordering", r[0].mode_derived < r[1].mode_derived &&
                           r[1].mode_derived < r[2].mode_derived);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
