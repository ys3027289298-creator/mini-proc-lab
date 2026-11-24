
/*
 * puddle: sub-pixel puddle coverage for distance LOD color blending.
 *
 * The same grassland region is rendered at three distances. The
 * existence-based (majority) path keeps a coarse cell as water only when
 * >= 50% of it is water, so puddles smaller than a cell are discarded and
 * the average color dries out with distance. The coverage path measures
 * the water fraction inside each coarse cell and folds it into a
 * grass/water lerp, conserving the high-resolution coverage (and therefore
 * the region's average color) at every distance.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define WORLD    16.0f   /* side length of the captured region, world units */
#define REF_RES  512     /* high-resolution reference grid */
#define MAXP     64
#define SUB      8       /* sub-samples per coarse cell axis */
#define THRESH   0.5f    /* majority threshold of the existence-based path */
#define TOL      0.008f  /* coverage / color conservation tolerance */
#define NTIERS  3

static const float GRASS[3] = {0.36f, 0.55f, 0.18f};
static const float WATER[3] = {0.10f, 0.25f, 0.45f};

typedef struct { float x, y, r; } Puddle;
typedef struct { Puddle p[MAXP]; int n; } Set;

static unsigned h1(unsigned i) {
    unsigned h = i * 2654435761u;
    h ^= h >> 15;
    return h;
}

static float u01(unsigned i) { return (float)(h1(i) & 0xffff) / 65535.0f; }

static void set_init(Set* s) { memset(s, 0, sizeof(*s)); }

/* min distance 0.9 > 2 * rmax keeps puddle disks from overlapping. */
static int set_add(Set* s, float x, float y, float r, float mind) {
    if (s->n >= MAXP) return -1;
    for (int i = 0; i < s->n; i++) {
        float dx = s->p[i].x - x, dy = s->p[i].y - y;
        if (dx * dx + dy * dy < mind * mind) return -2;
    }
    s->p[s->n].x = x; s->p[s->n].y = y; s->p[s->n].r = r;
    return s->n++;
}

static int make_puddles(Set* s, unsigned seed) {
    int g = 8, rejected = 0;
    float cell = WORLD / g;
    for (int i = 0; i < g; i++) {
        for (int j = 0; j < g; j++) {
            unsigned k = (unsigned)(i * 131 + j * 17) + seed * 7919u;
            float x = (i + u01(k)) * cell;
            float y = (j + u01(k + 1)) * cell;
            float r = 0.12f + 0.16f * u01(k + 2);
            if (set_add(s, x, y, r, 0.9f) == -2) rejected++;
        }
    }
    return rejected;
}

static int water_at(const Set* s, float x, float y) {
    for (int i = 0; i < s->n; i++) {
        float dx = x - s->p[i].x, dy = y - s->p[i].y;
        if (dx * dx + dy * dy <= s->p[i].r * s->p[i].r) return 1;
    }
    return 0;
}

static float reference_coverage(const Set* s) {
    float step = WORLD / REF_RES;
    long hits = 0;
    for (int iy = 0; iy < REF_RES; iy++)
        for (int ix = 0; ix < REF_RES; ix++)
            hits += water_at(s, (ix + 0.5f) * step, (iy + 0.5f) * step);
    return (float)hits / (REF_RES * REF_RES);
}

/* disks never overlap, so the exact coverage is the sum of disk areas. */
static float analytic_coverage(const Set* s) {
    float area = 0.0f;
    for (int i = 0; i < s->n; i++)
        area += 3.14159265f * s->p[i].r * s->p[i].r;
    return area / (WORLD * WORLD);
}

typedef struct {
    const char* name;
    int res;
    int visible;       /* puddles owning at least one majority-water cell */
    float cov;         /* mean sub-pixel water coverage over all cells */
    float hard_frac;   /* fraction of cells classified water by majority */
    float rgb_hard[3]; /* existence-based average color */
    float rgb_blend[3];/* coverage-blended average color */
    float mix_hard;    /* effective water mix projected back from rgb_hard */
    float mix_blend;
} Tier;

static unsigned char hard_grid[256 * 256];
static float cov_grid[256 * 256];

static void write_ppm(const char* path, int res, int blend) {
    FILE* f = fopen(path, "wb");
    if (!f) { perror(path); return; }
    fprintf(f, "P6\n%d %d\n255\n", res, res);
    for (int i = 0; i < res * res; i++) {
        float c = blend ? cov_grid[i] : (hard_grid[i] ? 1.0f : 0.0f);
        unsigned char px[3];
        for (int k = 0; k < 3; k++) {
            float v = GRASS[k] + (WATER[k] - GRASS[k]) * c;
            px[k] = (unsigned char)(fminf(1.0f, fmaxf(0.0f, v)) * 255.0f + 0.5f);
        }
        fwrite(px, 1, 3, f);
    }
    fclose(f);
}

static void write_reference_ppm(const Set* s, const char* path) {
    FILE* f = fopen(path, "wb");
    if (!f) { perror(path); return; }
    fprintf(f, "P6\n%d %d\n255\n", REF_RES, REF_RES);
    float step = WORLD / REF_RES;
    for (int iy = 0; iy < REF_RES; iy++) {
        for (int ix = 0; ix < REF_RES; ix++) {
            int w = water_at(s, (ix + 0.5f) * step, (iy + 0.5f) * step);
            const float* c = w ? WATER : GRASS;
            unsigned char px[3];
            for (int k = 0; k < 3; k++)
                px[k] = (unsigned char)(c[k] * 255.0f + 0.5f);
            fwrite(px, 1, 3, f);
        }
    }
    fclose(f);
}

/* project an rgb color back onto the grass->water axis */
static float mix_of(const float rgb[3]) {
    float num = 0.0f, den = 0.0f;
    for (int k = 0; k < 3; k++) {
        float d = WATER[k] - GRASS[k];
        num += (rgb[k] - GRASS[k]) * d;
        den += d * d;
    }
    return num / den;
}

static int imin(int a, int b) { return a < b ? a : b; }
static int imax(int a, int b) { return a > b ? a : b; }

static void render(const Set* s, const char* name, int res, Tier* t) {
    float cell = WORLD / res;
    t->name = name;
    t->res = res;
    float cov_sum = 0.0f, hard_sum = 0.0f;

    for (int iy = 0; iy < res; iy++) {
        for (int ix = 0; ix < res; ix++) {
            int hits = 0;
            for (int a = 0; a < SUB; a++)
                for (int b = 0; b < SUB; b++) {
                    float x = (ix + (a + 0.5f) / SUB) * cell;
                    float y = (iy + (b + 0.5f) / SUB) * cell;
                    hits += water_at(s, x, y);
                }
            float c = (float)hits / (SUB * SUB);
            int is_water = c >= THRESH;
            hard_grid[iy * res + ix] = (unsigned char)is_water;
            cov_grid[iy * res + ix] = c;
            cov_sum += c;
            hard_sum += is_water ? 1.0f : 0.0f;
        }
    }

    int total = res * res;
    t->cov = cov_sum / total;
    t->hard_frac = hard_sum / total;

    int vis = 0;
    for (int i = 0; i < s->n; i++) {
        int x0 = imax(0, (int)floorf((s->p[i].x - s->p[i].r) / cell));
        int x1 = imin(res - 1, (int)floorf((s->p[i].x + s->p[i].r) / cell));
        int y0 = imax(0, (int)floorf((s->p[i].y - s->p[i].r) / cell));
        int y1 = imin(res - 1, (int)floorf((s->p[i].y + s->p[i].r) / cell));
        for (int iy = y0; iy <= y1; iy++)
            for (int ix = x0; ix <= x1; ix++)
                if (hard_grid[iy * res + ix]) { vis++; goto found; }
found:;
    }
    t->visible = vis;

    for (int k = 0; k < 3; k++) {
        t->rgb_hard[k]  = GRASS[k] + (WATER[k] - GRASS[k]) * t->hard_frac;
        t->rgb_blend[k] = GRASS[k] + (WATER[k] - GRASS[k]) * t->cov;
    }
    t->mix_hard = mix_of(t->rgb_hard);
    t->mix_blend = mix_of(t->rgb_blend);
}

static int first_bad_tier(const Tier* tiers, float ref, int use_blend) {
    for (int i = 0; i < NTIERS; i++) {
        float m = use_blend ? tiers[i].mix_blend : tiers[i].mix_hard;
        if (fabsf(m - ref) > TOL) return i;
    }
    return -1;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(int argc, char** argv) {
    int dump = argc > 1;
    Set s;
    set_init(&s);
    int rejected = make_puddles(&s, 3u);

    float cov_ref = reference_coverage(&s);
    float cov_ana = analytic_coverage(&s);

    static const int res[NTIERS] = {256, 48, 12};
    static const char* names[NTIERS] = {"near", "mid", "far"};
    Tier tiers[NTIERS];
    for (int i = 0; i < NTIERS; i++) {
        render(&s, names[i], res[i], &tiers[i]);
        if (dump) {
            char path[256];
            snprintf(path, sizeof path, "%s/%s_hard.ppm", argv[1], names[i]);
            write_ppm(path, res[i], 0);
            snprintf(path, sizeof path, "%s/%s_blend.ppm", argv[1], names[i]);
            write_ppm(path, res[i], 1);
        }
    }
    if (dump) {
        char ref_path[256];
        snprintf(ref_path, sizeof ref_path, "%s/reference.ppm", argv[1]);
        write_reference_ppm(&s, ref_path);
    }

    printf("region=%.2fx%.2f puddles=%d rejected=%d\n",
           WORLD, WORLD, s.n, rejected);
    printf("coverage reference=%.4f analytic=%.4f\n", cov_ref, cov_ana);
    printf("reference rgb=(%.3f,%.3f,%.3f)\n",
           GRASS[0] + (WATER[0] - GRASS[0]) * cov_ref,
           GRASS[1] + (WATER[1] - GRASS[1]) * cov_ref,
           GRASS[2] + (WATER[2] - GRASS[2]) * cov_ref);
    printf("%-4s %4s %7s %8s %9s %9s %9s %8s %8s  %-19s %-19s\n",
           "tier", "res", "cell", "visible", "coverage", "hard_mix",
           "blend_mix", "d_hard", "d_blend", "hard_rgb", "blend_rgb");
    for (int i = 0; i < NTIERS; i++) {
        printf("%-4s %4d %7.4f %8d %9.4f %9.4f %9.4f %8.4f %8.4f  (%.3f,%.3f,%.3f)  (%.3f,%.3f,%.3f)\n",
               tiers[i].name, tiers[i].res, WORLD / tiers[i].res,
               tiers[i].visible, tiers[i].cov, tiers[i].mix_hard,
               tiers[i].mix_blend,
               fabsf(tiers[i].mix_hard - cov_ref),
               fabsf(tiers[i].mix_blend - cov_ref),
               tiers[i].rgb_hard[0], tiers[i].rgb_hard[1], tiers[i].rgb_hard[2],
               tiers[i].rgb_blend[0], tiers[i].rgb_blend[1], tiers[i].rgb_blend[2]);
    }

    int hard_bad = first_bad_tier(tiers, cov_ref, 0);
    int blend_bad = first_bad_tier(tiers, cov_ref, 1);
    printf("diagnosis:\n");
    printf("  majority path first biased tier = %s\n",
           hard_bad < 0 ? "none" : names[hard_bad]);
    if (hard_bad >= 0)
        printf("  there hard_mix=%.4f vs coverage=%.4f (delta=%.4f); far hard_mix=%.4f\n",
               tiers[hard_bad].mix_hard, cov_ref,
               fabsf(tiers[hard_bad].mix_hard - cov_ref),
               tiers[NTIERS - 1].mix_hard);
    printf("  coverage blend first biased tier = %s (tol=%.4f)\n",
           blend_bad < 0 ? "none" : names[blend_bad], TOL);

    /* self tests */
    check("analytic_matches_reference", fabsf(cov_ana - cov_ref) < 0.01f);

    check("near_all_puddles_visible", tiers[0].visible == s.n);
    check("far_puddles_filtered", tiers[2].visible == 0);
    check("visibility_nonincreasing",
          tiers[0].visible >= tiers[1].visible &&
          tiers[1].visible >= tiers[2].visible);
    check("hard_path_dries",
          tiers[2].mix_hard + TOL < cov_ref);

    for (int i = 0; i < NTIERS; i++) {
        char label[64];
        snprintf(label, sizeof label, "coverage_conserved_%s", tiers[i].name);
        check(label, fabsf(tiers[i].cov - cov_ref) <= TOL);
        snprintf(label, sizeof label, "blend_color_stable_%s", tiers[i].name);
        int ok = fabsf(tiers[i].mix_blend - cov_ref) <= TOL;
        for (int k = 0; k < 3 && ok; k++) {
            float refk = GRASS[k] + (WATER[k] - GRASS[k]) * cov_ref;
            ok = fabsf(tiers[i].rgb_blend[k] - refk) <= TOL;
        }
        check(label, ok);
    }

    check("diagnostic_flags_mid_for_majority", hard_bad == 1);
    check("diagnostic_no_bias_for_blend", blend_bad == -1);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
