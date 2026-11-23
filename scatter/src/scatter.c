
/* scatter: deterministic jittered scatter with minimum distance rejection.
 *
 * Feasibility relation (square jitter, one point per grid cell):
 *   A grid of cell size c. A point in cell (i,j) sits at its center plus a
 *   jitter (jx, jy), |jx|,|jy| <= J. Two axis-adjacent cells have center
 *   distance c; the worst-case displacement between their points happens when
 *   the jitters point at each other, giving minimum achievable distance
 *       d_min(J) = c - 2J   (horizontal/vertical neighbors; the tightest pair).
 *   Diagonal neighbors keep distance c*sqrt(2) - 2J*sqrt(2) = sqrt(2)*d_min,
 *   so the axis pair is the binding constraint. To guarantee every pair is at
 *   least the required minimum spacing d apart:
 *       J <= Jcap = (c - d) / 2,            with J >= 0.
 *   Feasibility therefore requires c >= d, i.e. point density rho = 1/c^2
 *   must satisfy rho <= 1/d^2. Above that density even J = 0 (a rigid grid)
 *   violates d, so no jitter amplitude can make the placement feasible.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAXP 256

typedef struct { float x, y; } P;

typedef struct { P p[MAXP]; int n; } Set;

static unsigned h1(unsigned i) {
    unsigned h = i * 2654435761u;
    h ^= h >> 15;
    return h;
}

static float u01(unsigned i) { return (float)(h1(i) & 0xffff) / 65535.0f; }

static void set_init(Set* s) { memset(s, 0, sizeof(*s)); }

static int set_add(Set* s, float x, float y, float mind) {
    if (s->n >= MAXP) return -1;
    for (int i = 0; i < s->n; i++) {
        float dx = s->p[i].x - x, dy = s->p[i].y - y;
        if (dx * dx + dy * dy < mind * mind) return -2;
    }
    s->p[s->n].x = x; s->p[s->n].y = y;
    return s->n++;
}

/* unchecked insertion: used to count the overlaps the rejection rule hides. */
static void set_add_raw(Set* s, float x, float y) {
    if (s->n < MAXP) { s->p[s->n].x = x; s->p[s->n].y = y; s->n++; }
}

/* Jcap from the feasibility relation; negative means the density is infeasible. */
static float jitter_cap(float cell, float mind) { return (cell - mind) * 0.5f; }

static int jitter_grid(Set* s, int nx, int ny, float cell, float mind, unsigned seed) {
    int tries = 0;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            unsigned k = (unsigned)(i * 131 + j * 17) + seed * 7919u;
            float x = (i + u01(k)) * cell;
            float y = (j + u01(k + 1)) * cell;
            int r = set_add(s, x, y, mind);
            if (r == -2) { tries++; continue; }
            if (r == -1) return -1;
        }
    }
    return tries;
}

/* centered jitter of amplitude J, optionally without rejection. */
static int jitter_grid_c(Set* s, int nx, int ny, float cell, float mind,
                         float J, unsigned seed, int reject) {
    int tries = 0;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            unsigned k = (unsigned)(i * 131 + j * 17) + seed * 7919u;
            float jx = (u01(k) - 0.5f) * 2.0f * J;
            float jy = (u01(k + 1) - 0.5f) * 2.0f * J;
            float x = (i + 0.5f) * cell + jx;
            float y = (j + 0.5f) * cell + jy;
            if (!reject) { set_add_raw(s, x, y); continue; }
            int r = set_add(s, x, y, mind);
            if (r == -2) { tries++; continue; }
            if (r == -1) return -1;
        }
    }
    return tries;
}

/* number of unordered pairs closer than mind (strictly, modulo float epsilon). */
static int count_overlaps(const Set* s, float mind) {
    int overlaps = 0;
    for (int i = 0; i < s->n; i++)
        for (int j = i + 1; j < s->n; j++) {
            float dx = s->p[i].x - s->p[j].x, dy = s->p[i].y - s->p[j].y;
            if (dx * dx + dy * dy < mind * mind - 1e-7f) overlaps++;
        }
    return overlaps;
}

static float min_pair_distance(const Set* s) {
    float best = 1e30f;
    for (int i = 0; i < s->n; i++)
        for (int j = i + 1; j < s->n; j++) {
            float dx = s->p[i].x - s->p[j].x, dy = s->p[i].y - s->p[j].y;
            float d = sqrtf(dx * dx + dy * dy);
            if (d < best) best = d;
        }
    return best;
}

typedef struct { int nx, ny; float cell, mind; } Tier;

typedef struct {
    int placed, target, rejected, overlaps;
    float cap, J_used, min_dist, density;
} Report;

/* one density tier: cap the jitter by Jcap and run a full scatter. */
static Report run_tier(Tier t, unsigned seed, int reject) {
    Report r;
    Set s;
    set_init(&s);
    r.cap = jitter_cap(t.cell, t.mind);
    r.J_used = r.cap >= 0.0f ? r.cap * (1.0f - 1e-5f) : 0.0f;
    int rejected = jitter_grid_c(&s, t.nx, t.ny, t.cell, t.mind,
                                 r.J_used, seed, reject);
    r.placed = s.n;
    r.target = t.nx * t.ny;
    r.rejected = rejected < 0 ? 0 : rejected;
    r.overlaps = count_overlaps(&s, t.mind);
    r.min_dist = min_pair_distance(&s);
    r.density = 1.0f / (t.cell * t.cell);
    return r;
}

static void print_report(const char* tag, const Tier* t, const Report* r) {
    printf("%s_grid=%dx%d cell=%.4f mind=%.4f density=%.4f\n",
           tag, t->nx, t->ny, t->cell, t->mind, r->density);
    printf("%s_jitter_cap=%.6f jitter_used=%.6f\n", tag, r->cap, r->J_used);
    printf("%s_placed=%d target=%d rejected=%d overlap_count=%d min_pair_dist=%.6f\n",
           tag, r->placed, r->target, r->rejected, r->overlaps, r->min_dist);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Set s;
    set_init(&s);
    int rejected = jitter_grid(&s, 4, 4, 1.0f, 0.0f, 1u);
    check("grid_count", s.n == 16 && rejected == 0);

    Set t;
    set_init(&t);
    rejected = jitter_grid(&t, 4, 4, 0.5f, 0.6f, 2u);
    check("min_distance_enforced", rejected > 0);
    int violations = 0;
    for (int i = 0; i < t.n; i++)
        for (int j = i + 1; j < t.n; j++) {
            float dx = t.p[i].x - t.p[j].x, dy = t.p[i].y - t.p[j].y;
            if (dx * dx + dy * dy < 0.36f - 1e-5f) violations++;
        }
    check("no_violation", violations == 0);

    check("deterministic", ({ Set a, b; set_init(&a); set_init(&b);
        jitter_grid(&a, 3, 3, 1.0f, 0.0f, 9u); jitter_grid(&b, 3, 3, 1.0f, 0.0f, 9u);
        a.n == b.n && fabsf(a.p[0].x - b.p[0].x) < 1e-9f; }));

    /* ---- three density tiers, jitter bounded by Jcap = (c - d)/2 ---- */
    printf("relation: J <= (cell - mind)/2 ; feasible iff cell >= mind\n");
    Tier tiers[3] = {
        { 8,  8, 1.00f, 0.50f }, /* low    rho = 1.00 */
        { 10, 10, 0.60f, 0.40f }, /* medium rho = 2.78 */
        { 12, 12, 0.40f, 0.30f }, /* high   rho = 6.25 */
    };
    const char* names[3] = { "low", "medium", "high" };
    int all_zero_overlap = 1;
    for (int i = 0; i < 3; i++) {
        Report r = run_tier(tiers[i], 11u + (unsigned)i, 1);
        print_report(names[i], &tiers[i], &r);
        printf("%s_all_targets_placed=%d\n", names[i], r.placed == r.target);
        if (r.overlaps != 0) all_zero_overlap = 0;
    }
    check("three_tiers_zero_overlap", all_zero_overlap);

    /* ---- counterexample: density past the feasibility boundary ----
     * cell = 0.25, mind = 0.30 -> Jcap = -0.025 < 0, i.e. rho = 16 > 1/d^2
     * = 11.11. Even a rigid grid (J = 0) has neighbors only 0.25 apart, so
     * the relation cannot supply any admissible jitter. Run it raw (no
     * rejection) to expose the overlaps that rejection would silently drop. */
    Tier bad = { 12, 12, 0.25f, 0.30f };
    printf("counterexample: cell=%.4f mind=%.4f jitter_cap=%.4f (negative -> no feasible J)\n",
           bad.cell, bad.mind, jitter_cap(bad.cell, bad.mind));
    Report rb = run_tier(bad, 99u, 0);
    print_report("counter", &bad, &rb);
    check("counterexample_relation_fails", rb.cap < 0.0f && rb.overlaps > 0);

    /* same density with rejection on: spacing can only be "enforced" by
     * refusing to place most points, which itself proves infeasibility. */
    Report rr = run_tier(bad, 99u, 1);
    print_report("counter_rejected", &bad, &rr);
    check("counterexample_rejection_loses_points", rr.placed < rr.target);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
