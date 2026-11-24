
/* tunnel: force-align a mine centerline to the terrain gradient.
 *
 * Pieces:
 *   - three analytic heightfields with exact gradients (no finite-diff drift);
 *   - ideal centerline = gradient flow line, integrated by RK2 at fine arc step;
 *   - per-segment curvature and the derived lateral-offset budget
 *       d_allow(k) = min(e_chord, (1-cos alpha)/kappa_k);
 *   - old scheme: fixed 4 m Euler chords, one gradient sample per station;
 *   - new scheme: adaptive chord length, every station anchored on the flow line;
 *   - independent per-segment measurement of lateral deviation and heading error.
 *
 * `make test` runs the embedded self-tests. Run with an output dir as argv[1]
 * to additionally dump CSVs for plotting.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846f

/* ---- engineering constants (see DERIVATION.md) ---- */
#define DS          0.002f   /* flow-line integration arc step, m          */
#define DOMAIN      48.0f    /* square heightfield domain [0,DOMAIN]^2 m   */
#define MARGIN      1.5f     /* stop the flow before leaving the domain    */
#define PATH_LEN    40.0f    /* target centerline arc length, m            */
#define ALPHA_DEG   8.0f     /* max allowed heading error                  */
#define E_CHORD     0.20f    /* straight-chord construction tolerance, m   */
#define L_OLD       4.0f     /* old scheme: fixed station spacing, m       */
#define L_MAX       2.0f     /* new scheme: max chord length, m            */
#define ETA         0.9f     /* budget utilisation target when closing     */

#define NMAX 32000          /* 64 m / 0.002 m */
#define MAXSEG 256

typedef struct { float x, y; } V2;

static float alpha_rad(void) { return ALPHA_DEG * PI / 180.0f; }

/* ---------------- analytic terrain ---------------- */

typedef enum { T_ROLLING = 0, T_RIDGE = 1, T_MONOCLINE = 2 } TerrainId;

typedef struct {
    TerrainId id;
    const char* name;
    V2 start;
} Terrain;

static float terrain_h(const Terrain* t, float x, float y) {
    switch (t->id) {
    case T_ROLLING:
        return 0.25f * x
             + 0.8f * sinf(0.14f * x) * cosf(0.11f * y + 0.3f)
             + 0.35f * sinf(0.10f * y + 1.1f);
    case T_RIDGE:
        return 0.45f * x + 0.10f * y
             + 4.0f * expf(-((y - (10.0f + 0.30f * x)) *
                             (y - (10.0f + 0.30f * x))) / (2.0f * 25.0f));
    case T_MONOCLINE: {
        float q = y - (26.0f - 0.45f * x);
        return 1.15f * x
             + 6.0f * expf(-(q * q) / (2.0f * 12.25f))
             + 1.2f * sinf(0.18f * x + 0.7f) * cosf(0.16f * y);
    }}
    return 0.0f;
}

/* exact analytic gradient */
static V2 terrain_grad(const Terrain* t, float x, float y) {
    V2 g;
    switch (t->id) {
    case T_ROLLING:
        g.x = 0.25f
            + 0.8f * 0.14f * cosf(0.14f * x) * cosf(0.11f * y + 0.3f);
        g.y = 0.8f * (-0.11f) * sinf(0.14f * x) * sinf(0.11f * y + 0.3f)
            + 0.35f * 0.10f * cosf(0.10f * y + 1.1f);
        return g;
    case T_RIDGE: {
        float q = y - (10.0f + 0.30f * x);
        float e = expf(-(q * q) / 50.0f);
        g.x = 0.45f + 4.0f * e * (2.0f * q * 0.30f / 50.0f);
        g.y = 0.10f + 4.0f * e * (-2.0f * q / 50.0f);
        return g;
    }
    case T_MONOCLINE: {
        float q = y - (26.0f - 0.45f * x);
        float e = expf(-(q * q) / 24.5f);
        /* dq/dx = +0.45, dq/dy = 1 */
        g.x = 1.15f
            + 6.0f * e * (-2.0f * q * 0.45f / 24.5f)
            + 1.2f * 0.18f * cosf(0.18f * x + 0.7f) * cosf(0.16f * y);
        g.y = 6.0f * e * (-2.0f * q / 24.5f)
            + 1.2f * (-0.16f) * sinf(0.18f * x + 0.7f) * sinf(0.16f * y);
        return g;
    }}
    g.x = g.y = 0.0f;
    return g;
}

/* ---------------- vector helpers ---------------- */

static V2 v2(float x, float y) { V2 r = { x, y }; return r; }
static V2 vadd(V2 a, V2 b) { return v2(a.x + b.x, a.y + b.y); }
static V2 vscale(V2 a, float k) { return v2(a.x * k, a.y * k); }
static float vdot(V2 a, V2 b) { return a.x * b.x + a.y * b.y; }
static float vlen(V2 a) { return sqrtf(vdot(a, a)); }
static V2 vnorm(V2 a) {
    float l = vlen(a);
    return l > 1e-12f ? v2(a.x / l, a.y / l) : v2(1.0f, 0.0f);
}
static float clampf(float v, float lo, float hi) {
    return v < lo ? lo : (v > hi ? hi : v);
}
static float angle_between(V2 a, V2 b) {
    return acosf(clampf(vdot(vnorm(a), vnorm(b)), -1.0f, 1.0f));
}
/* perpendicular distance from p to the infinite line through a-b */
static float point_line_dist(V2 p, V2 a, V2 b) {
    V2 d = vadd(b, vscale(a, -1.0f));
    float l = vlen(d);
    if (l < 1e-12f) return vlen(vadd(p, vscale(a, -1.0f)));
    return fabsf(d.x * (a.y - p.y) - (a.x - p.x) * d.y) / l;
}

/* ---------------- ideal gradient flow line ---------------- */

typedef struct {
    V2 p[NMAX];      /* positions */
    V2 t[NMAX];      /* unit tangents (= normalized gradient) */
    float s[NMAX];   /* arc length */
    int n;
} Flow;

static V2 grad_hat(const Terrain* tr, V2 p) {
    return vnorm(terrain_grad(tr, p.x, p.y));
}

/* RK2 (midpoint) integration along the steepest-up direction */
static int flow_build(const Terrain* tr, Flow* f) {
    f->n = 0;
    V2 p = tr->start;
    float s = 0.0f;
    f->p[0] = p;
    f->t[0] = grad_hat(tr, p);
    f->s[0] = 0.0f;
    f->n = 1;
    while (s < PATH_LEN && f->n < NMAX - 1) {
        V2 k1 = grad_hat(tr, p);
        V2 mid = vadd(p, vscale(k1, 0.5f * DS));
        V2 k2 = grad_hat(tr, mid);
        p = vadd(p, vscale(k2, DS));
        s += DS;
        if (p.x < MARGIN || p.y < MARGIN ||
            p.x > DOMAIN - MARGIN || p.y > DOMAIN - MARGIN) break;
        f->p[f->n] = p;
        f->t[f->n] = k2;
        f->s[f->n] = s;
        f->n++;
    }
    return f->n;
}

/* curvature at node i. A three-point Menger formula on 2 mm-spaced float32
 * positions is swamped by position quantisation, so instead differentiate the
 * analytic unit tangents: turning angle over a 2*STENCIL*DS arc window, with
 * atan2 (accurate where the angle is small) instead of acos. */
#define KAPPA_STENCIL 50  /* half-window in nodes; full window = 0.2 m */
static float flow_curvature(const Flow* f, int i) {
    if (i <= 0 || i >= f->n - 1) return 0.0f;
    int m = KAPPA_STENCIL;
    if (i - m < 0) m = i;
    if (i + m > f->n - 1) m = (f->n - 1) - i;
    if (m <= 0) return 0.0f;
    V2 ta = f->t[i - m], tb = f->t[i + m];
    float cross = ta.x * tb.y - ta.y * tb.x;
    float dot = ta.x * tb.x + ta.y * tb.y;
    return fabsf(atan2f(cross, dot)) / (2.0f * (float)m * DS);
}

/* ---------------- offset budget: the delivered expression ---------------- */

static float d_allow(float kappa_max) {
    float head = (1.0f - cosf(alpha_rad())) / kappa_max;
    return kappa_max <= 1e-9f ? E_CHORD : fminf(E_CHORD, head);
}

/* ---------------- measured stats over an arc window [i0,i1] ---------------- */

typedef struct {
    float kappa_max;
    float d_max;       /* max lateral deviation of flow points from chord */
    float theta_max;   /* max heading error, radians                      */
    float budget;      /* d_allow for the window                          */
    int pass;
} WindowStat;

static WindowStat measure_window(const Flow* f, int i0, int i1, V2 a, V2 b) {
    WindowStat w;
    memset(&w, 0, sizeof(w));
    V2 chord = vadd(b, vscale(a, -1.0f));
    for (int i = i0; i <= i1 && i < f->n; i++) {
        float k = flow_curvature(f, i);
        if (k > w.kappa_max) w.kappa_max = k;
        float d = point_line_dist(f->p[i], a, b);
        if (d > w.d_max) w.d_max = d;
        float th = angle_between(chord, f->t[i]);
        if (th > w.theta_max) w.theta_max = th;
    }
    w.budget = d_allow(w.kappa_max);
    w.pass = w.d_max <= w.budget &&
             w.theta_max <= alpha_rad() + 1e-6f;
    return w;
}

/* ---------------- OLD scheme: fixed 4 m Euler chords ---------------- */

typedef struct {
    V2 p[MAXSEG + 1];
    int i0[MAXSEG], i1[MAXSEG]; /* flow-node window of each chord */
    WindowStat seg[MAXSEG];
    int n;                      /* number of chords */
    int all_pass;
} Scheme;

static void scheme_old(const Terrain* tr, const Flow* f, Scheme* sc) {
    sc->n = 0;
    V2 p = f->p[0];
    int steps_per = (int)(L_OLD / DS + 0.5f);
    int nchords = (f->n - 1) / steps_per;
    for (int k = 0; k < nchords && sc->n < MAXSEG; k++) {
        V2 dir = grad_hat(tr, p);              /* one sample per station */
        V2 q = vadd(p, vscale(dir, L_OLD));    /* straight 4 m Euler step */
        int i0 = k * steps_per;
        int i1 = (k + 1) * steps_per;
        if (i1 >= f->n) i1 = f->n - 1;
        sc->p[k] = p;
        sc->i0[k] = i0;
        sc->i1[k] = i1;
        sc->seg[k] = measure_window(f, i0, i1, p, q);
        sc->n++;
        p = q;
    }
    sc->p[sc->n] = p;
    sc->all_pass = 1;
    for (int k = 0; k < sc->n; k++)
        if (!sc->seg[k].pass) sc->all_pass = 0;
}

/* ---------------- NEW scheme: adaptive chords on the flow line ---------------- */

/* largest sagitta of flow nodes [i0..i1] away from the candidate chord */
static float window_worst_sagitta(const Flow* f, int i0, int i1) {
    float worst = 0.0f;
    for (int j = i0; j <= i1; j++) {
        float dj = point_line_dist(f->p[j], f->p[i0], f->p[i1]);
        if (dj > worst) worst = dj;
    }
    return worst;
}

static float window_kappa_max(const Flow* f, int i0, int i1) {
    float kmax = 0.0f;
    for (int j = i0; j <= i1 && j < f->n; j++) {
        float k = flow_curvature(f, j);
        if (k > kmax) kmax = k;
    }
    return kmax;
}

static void scheme_new(const Flow* f, Scheme* sc) {
    sc->n = 0;
    int i0 = 0;
    while (i0 < f->n - 1 && sc->n < MAXSEG) {
        int i1 = i0 + 1;
        /* grow the chord until the measured sagitta hits its budget share
           or the chord reaches L_MAX; the budget uses the window's own
           max curvature, so every closing decision is self-contained */
        while (i1 < f->n - 1) {
            float arc = f->s[i1] - f->s[i0];
            if (arc >= L_MAX) break;
            float budget = d_allow(window_kappa_max(f, i0, i1 + 1));
            if (window_worst_sagitta(f, i0, i1 + 1) >= ETA * budget) break;
            i1++;
        }
        sc->p[sc->n] = f->p[i0];
        sc->i0[sc->n] = i0;
        sc->i1[sc->n] = i1;
        sc->seg[sc->n] = measure_window(f, i0, i1, f->p[i0], f->p[i1]);
        sc->n++;
        i0 = i1;
    }
    sc->p[sc->n] = f->p[i0];
    sc->all_pass = 1;
    for (int k = 0; k < sc->n; k++)
        if (!sc->seg[k].pass) sc->all_pass = 0;
}

/* ---------------- terrain-level summary ---------------- */

typedef struct {
    float drop;
    float kappa_peak;
    float d_max;
    float d_allow_min;
    float util_max;       /* max d_max / budget over segments */
    float theta_deg_max;
} Summary;

static Summary summarize(const Flow* f, const Scheme* sc) {
    Summary su;
    memset(&su, 0, sizeof(su));
    su.d_allow_min = 1e9f;
    for (int i = 0; i < f->n; i++) {
        float k = flow_curvature(f, i);
        if (k > su.kappa_peak) su.kappa_peak = k;
    }
    su.drop = 0.0f; /* overwritten by the caller with the endpoint height drop */
    for (int k = 0; k < sc->n; k++) {
        if (sc->seg[k].d_max > su.d_max) su.d_max = sc->seg[k].d_max;
        if (sc->seg[k].budget < su.d_allow_min) su.d_allow_min = sc->seg[k].budget;
        float u = sc->seg[k].d_max / sc->seg[k].budget;
        if (u > su.util_max) su.util_max = u;
        if (sc->seg[k].theta_max > su.theta_deg_max)
            su.theta_deg_max = sc->seg[k].theta_max;
    }
    su.theta_deg_max = su.theta_deg_max * 180.0f / PI;
    return su;
}

/* ---------------- CSV output ---------------- */

static void write_csv(const char* dir, const Terrain* tr,
                      const Flow* f, const Scheme* old, const Scheme* nw) {
    char path[512];
    FILE* fp;

    snprintf(path, sizeof(path), "%s/stream_%s.csv", dir, tr->name);
    fp = fopen(path, "w");
    fprintf(fp, "s,x,y,h,kappa\n");
    for (int i = 0; i < f->n; i++)
        fprintf(fp, "%.4f,%.5f,%.5f,%.5f,%.6f\n",
                f->s[i], f->p[i].x, f->p[i].y,
                terrain_h(tr, f->p[i].x, f->p[i].y), flow_curvature(f, i));
    fclose(fp);

    snprintf(path, sizeof(path), "%s/new_%s.csv", dir, tr->name);
    fp = fopen(path, "w");
    fprintf(fp, "x,y\n");
    for (int i = 0; i <= nw->n; i++)
        fprintf(fp, "%.5f,%.5f\n", nw->p[i].x, nw->p[i].y);
    fclose(fp);

    snprintf(path, sizeof(path), "%s/old_%s.csv", dir, tr->name);
    fp = fopen(path, "w");
    fprintf(fp, "x,y\n");
    for (int i = 0; i <= old->n; i++)
        fprintf(fp, "%.5f,%.5f\n", old->p[i].x, old->p[i].y);
    fclose(fp);

    /* per-segment measurement table, one row per chord of each scheme */
    snprintf(path, sizeof(path), "%s/segstats_%s.csv", dir, tr->name);
    fp = fopen(path, "w");
    fprintf(fp, "scheme,seg,s_start,s_end,arc_len,kappa_max,d_allow,"
                "d_max,theta_deg,pass\n");
    for (int k = 0; k < nw->n; k++)
        fprintf(fp, "new,%d,%.4f,%.4f,%.4f,%.7f,%.6f,%.6f,%.4f,%d\n",
                k, f->s[nw->i0[k]], f->s[nw->i1[k]],
                f->s[nw->i1[k]] - f->s[nw->i0[k]],
                nw->seg[k].kappa_max, nw->seg[k].budget,
                nw->seg[k].d_max, nw->seg[k].theta_max * 180.0f / PI,
                nw->seg[k].pass);
    for (int k = 0; k < old->n; k++)
        fprintf(fp, "old,%d,%.4f,%.4f,%.4f,%.7f,%.6f,%.6f,%.4f,%d\n",
                k, f->s[old->i0[k]], f->s[old->i1[k]],
                f->s[old->i1[k]] - f->s[old->i0[k]],
                old->seg[k].kappa_max, old->seg[k].budget,
                old->seg[k].d_max, old->seg[k].theta_max * 180.0f / PI,
                old->seg[k].pass);
    fclose(fp);
}

/* ---------------- self tests ---------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static int nearf(float a, float b, float tol) { return fabsf(a - b) < tol; }

int main(int argc, char** argv) {
    const char* outdir = argc > 1 ? argv[1] : NULL;
    const Terrain terrains[3] = {
        { T_ROLLING,  "rolling_hills",       v2(6.0f, 24.0f) },
        { T_RIDGE,    "single_ridge",        v2(3.0f, 24.0f) },
        { T_MONOCLINE,"steep_bent_monocline",v2(3.0f, 12.0f) },
    };

    /* ---- test 1: analytic gradients match central finite differences ---- */
    {
        int ok = 1;
        float pts[][2] = { {6, 7}, {20, 20}, {33, 11}, {40, 35} };
        for (int ti = 0; ti < 3; ti++) {
            for (unsigned j = 0; j < sizeof(pts) / sizeof(pts[0]); j++) {
                /* e=1e-2 keeps the subtraction well above float32 round-off */
                float x = pts[j][0], y = pts[j][1], e = 1e-2f;
                V2 g = terrain_grad(&terrains[ti], x, y);
                float gx = (terrain_h(&terrains[ti], x + e, y) -
                            terrain_h(&terrains[ti], x - e, y)) / (2 * e);
                float gy = (terrain_h(&terrains[ti], x, y + e) -
                            terrain_h(&terrains[ti], x, y - e)) / (2 * e);
                if (!nearf(g.x, gx, 1e-3f) || !nearf(g.y, gy, 1e-3f)) ok = 0;
            }
        }
        check("grad_matches_finite_diff", ok);
    }

    /* ---- test 2: budget expression, flat limit and identity ---- */
    check("budget_flat_floor", nearf(d_allow(0.0f), E_CHORD, 1e-7f));
    {
        float kappa = 0.12f;
        float expect = (1.0f - cosf(alpha_rad())) / kappa;
        check("budget_heading_branch", nearf(d_allow(kappa), expect, 1e-7f));
        check("budget_min_of_two",
              d_allow(0.001f) == E_CHORD &&
              nearf(d_allow(0.5f), (1 - cosf(alpha_rad())) / 0.5f, 1e-7f));
        /* sagitta identity: arc of curvature kappa, chord half-angle alpha */
        float k = 0.1f, d = (1 - cosf(alpha_rad())) / k;
        float theta = acosf(1 - k * d);
        check("sagitta_identity", nearf(theta, alpha_rad(), 1e-6f));
    }

    /* ---- test 3: curvature estimator on a manufactured circle ---- */
    {
        Flow f; memset(&f, 0, sizeof(f));
        float R = 7.0f, n = 2000;
        f.n = (int)n;
        for (int i = 0; i < f.n; i++) {
            float a = (float)i * DS / R;
            f.p[i] = v2(R * cosf(a), R * sinf(a));
            f.t[i] = v2(-sinf(a), cosf(a));
            f.s[i] = (float)i * DS;
        }
        float kmeas = flow_curvature(&f, 1000);
        check("curvature_on_circle",
              fabsf(kmeas - 1.0f / R) / (1.0f / R) < 0.01f);
    }

    /* ---- test 4: straight manufactured path has zero curvature ---- */
    {
        Flow f; memset(&f, 0, sizeof(f));
        V2 p = v2(5, 5), dir = vnorm(v2(1, 1));
        f.n = 0;
        for (int i = 0; i < 1000; i++) {
            f.p[f.n] = p; f.t[f.n] = dir; f.s[f.n] = i * DS; f.n++;
            p = vadd(p, vscale(dir, DS));
        }
        float kmax = 0;
        for (int i = 1; i < f.n - 1; i++) {
            float kk = flow_curvature(&f, i);
            if (kk > kmax) kmax = kk;
        }
        check("straight_line_zero_curvature", kmax < 1e-9f);
    }

    /* ---- end to end on the three terrains ---- */
    FILE* summary = NULL;
    if (outdir) {
        char path[512];
        snprintf(path, sizeof(path), "%s/summary.csv", outdir);
        summary = fopen(path, "w");
        fprintf(summary, "id,terrain,drop_m,kappa_peak,segments,"
                "d_max_m,d_allow_min_m,util_max,theta_deg_max,old_d_max_m,"
                "old_theta_deg_max,old_segments,old_failed_segs,pass\n");
    }
    int all_terrains_pass = 1;
    for (int ti = 0; ti < 3; ti++) {
        const Terrain* tr = &terrains[ti];
        static Flow f;
        static Scheme old, nw;
        int n = flow_build(tr, &f);
        scheme_old(tr, &f, &old);
        scheme_new(&f, &nw);

        char buf[96];
        snprintf(buf, sizeof(buf), "flow_length_%s", tr->name);
        check(buf, f.s[n - 1] >= PATH_LEN - 0.01f);

        Summary su = summarize(&f, &nw);
        su.drop = terrain_h(tr, f.p[f.n - 1].x, f.p[f.n - 1].y) -
                  terrain_h(tr, f.p[0].x, f.p[0].y);

        int old_fail = 0;
        float old_dmax = 0, old_th = 0;
        for (int k = 0; k < old.n; k++) {
            if (!old.seg[k].pass) old_fail++;
            if (old.seg[k].d_max > old_dmax) old_dmax = old.seg[k].d_max;
            if (old.seg[k].theta_max > old_th) old_th = old.seg[k].theta_max;
        }
        old_th *= 180.0f / PI;

        snprintf(buf, sizeof(buf), "new_scheme_all_segments_pass_%s", tr->name);
        check(buf, nw.all_pass);
        snprintf(buf, sizeof(buf), "new_heading_within_alpha_%s", tr->name);
        int head_ok = 1;
        for (int k = 0; k < nw.n; k++)
            if (nw.seg[k].theta_max > alpha_rad() + 1e-6f) head_ok = 0;
        check(buf, head_ok);
        if (!nw.all_pass) all_terrains_pass = 0;

        printf("[%s] drop=%.2fm kappa_peak=%.4f/m segs=%d d_max=%.4fm "
               "d_allow_min=%.4fm util=%.2f%% theta_max=%.2fdeg "
               "old_dmax=%.3fm old_fail_segs=%d/%d pass=%d\n",
               tr->name, su.drop, su.kappa_peak, nw.n, su.d_max,
               su.d_allow_min, 100 * su.util_max, su.theta_deg_max,
               old_dmax, old_fail, old.n, nw.all_pass);

        if (summary)
            fprintf(summary, "%d,%s,%.3f,%.5f,%d,%.4f,%.4f,%.4f,%.3f,%.4f,%.3f,%d,%d,%d\n",
                    ti, tr->name, su.drop, su.kappa_peak, nw.n,
                    su.d_max, su.d_allow_min, su.util_max, su.theta_deg_max,
                    old_dmax, old_th, old.n, old_fail, nw.all_pass);
        if (outdir) write_csv(outdir, tr, &f, &old, &nw);
    }
    /* determinism: rebuild terrain 0 and compare */
    {
        static Flow f1, f2;
        flow_build(&terrains[0], &f1);
        flow_build(&terrains[0], &f2);
        int same = f1.n == f2.n;
        for (int i = 0; same && i < f1.n; i++)
            if (!nearf(f1.p[i].x, f2.p[i].x, 1e-9f) ||
                !nearf(f1.p[i].y, f2.p[i].y, 1e-9f)) same = 0;
        check("deterministic", same);
    }
    check("all_three_terrains_pass", all_terrains_pass);
    if (summary) fclose(summary);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
