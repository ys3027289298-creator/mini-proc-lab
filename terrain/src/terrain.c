
/* terrain: heightfield sampling with slope and bilinear interpolation.
 *
 * Gradient/normal estimation: interior nodes use the 2nd-order centered
 * difference; boundary nodes use a 3-point one-sided 2nd-order stencil so the
 * O(h^2) truncation order (and hence the normal-angle bound, see
 * BOUNDARY_NORMALS.md) is retained at all four grid edges. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define W 16
#define H 16
#define PI 3.14159265358979323846

typedef struct { float h[W][H]; float cell; } Field;

/* boundary one-sided stencil needs two interior nodes per axis */
_Static_assert(W >= 3 && H >= 3, "boundary gradient stencil requires W,H >= 3");

static void field_init(Field* f, float cell) {
    memset(f, 0, sizeof(*f));
    f->cell = cell;
}

static void field_set(Field* f, int x, int y, float v) {
    if (x < 0 || y < 0 || x >= W || y >= H) return;
    f->h[x][y] = v;
}

static float field_sample(const Field* f, float wx, float wy) {
    float gx = wx / f->cell, gy = wy / f->cell;
    int x0 = (int)floorf(gx), y0 = (int)floorf(gy);
    float tx = gx - x0, ty = gy - y0;
    int x1 = x0 + 1, y1 = y0 + 1;
    if (x0 < 0) x0 = 0; if (y0 < 0) y0 = 0;
    if (x1 >= W) x1 = W - 1; if (y1 >= H) y1 = H - 1;
    float a = f->h[x0][y0] * (1 - tx) + f->h[x1][y0] * tx;
    float b = f->h[x0][y1] * (1 - tx) + f->h[x1][y1] * tx;
    return a * (1 - ty) + b * ty;
}

/* Boundary-aware gradient.
 * interior : D0 = (f+ - f-)/(2h),          |e| <= M3 h^2/6
 * edge     : Df = (-3f0 + 4f1 - f2)/(2h),  |e| <= M3 h^2/3
 *            Db = ( 3f0 - 4f-1 + f-2)/(2h) (same order/constant) */
static void field_grad(const Field* f, int x, int y, float* gx, float* gy) {
    float hc = f->cell;
    if (x > 0 && x < W - 1) {
        *gx = (f->h[x + 1][y] - f->h[x - 1][y]) / (2.0f * hc);
    } else if (x == 0) {
        *gx = (-3.0f * f->h[0][y] + 4.0f * f->h[1][y] - f->h[2][y]) / (2.0f * hc);
    } else {
        *gx = ( 3.0f * f->h[W - 1][y] - 4.0f * f->h[W - 2][y]
                + f->h[W - 3][y]) / (2.0f * hc);
    }
    if (y > 0 && y < H - 1) {
        *gy = (f->h[x][y + 1] - f->h[x][y - 1]) / (2.0f * hc);
    } else if (y == 0) {
        *gy = (-3.0f * f->h[x][0] + 4.0f * f->h[x][1] - f->h[x][2]) / (2.0f * hc);
    } else {
        *gy = ( 3.0f * f->h[x][H - 1] - 4.0f * f->h[x][H - 2]
                + f->h[x][H - 3]) / (2.0f * hc);
    }
}

/* Legacy boundary handling: index clamping degrades the centered stencil to a
 * 1st-order one-sided difference at the edges (|e| = M2 h/2 + O(h^2)). */
static void field_grad_old(const Field* f, int x, int y, float* gx, float* gy) {
    int xm = x > 0 ? x - 1 : x, xp = x < W - 1 ? x + 1 : x;
    int ym = y > 0 ? y - 1 : y, yp = y < H - 1 ? y + 1 : y;
    *gx = (f->h[xp][y] - f->h[xm][y]) / ((xp - xm) * f->cell);
    *gy = (f->h[x][yp] - f->h[x][ym]) / ((yp - ym) * f->cell);
}

static float field_slope(const Field* f, int x, int y) {
    float gx, gy;
    field_grad(f, x, y, &gx, &gy);
    return sqrtf(gx * gx + gy * gy);
}

/* Unit normal of the heightfield surface z = f(x,y): n = (-fx,-fy,1)/|.| */
static void grad_normal(float gx, float gy, float n[3]) {
    float inv = 1.0f / sqrtf(gx * gx + gy * gy + 1.0f);
    n[0] = -gx * inv;
    n[1] = -gy * inv;
    n[2] = inv;
}

static float angle_deg(const float a[3], const float b[3]) {
    float d = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    if (d > 1.0f) d = 1.0f;
    if (d < -1.0f) d = -1.0f;
    return acosf(d) * (float)(180.0 / PI);
}

/* ---- analytic cubic reference surface ------------------------------------
 * f = .01 x^3 + .02 y^3 + .3 x^2 + .4 y^2 + .1 xy + 1
 * constant third derivatives make the O(h^2) stencil bound tight. */
static float surf_f(double x, double y) {
    return (float)(0.01 * x * x * x + 0.02 * y * y * y + 0.3 * x * x
                   + 0.4 * y * y + 0.1 * x * y + 1.0);
}
static void surf_grad(double x, double y, double* gx, double* gy) {
    *gx = 0.03 * x * x + 0.6 * x + 0.1 * y;   /* fxxx = 0.06 */
    *gy = 0.06 * y * y + 0.8 * y + 0.1 * x;   /* fyyy = 0.12 */
}
#define M3X 0.06
#define M3Y 0.12

static void build_surface(Field* f, double h) {
    field_init(f, (float)h);
    for (int x = 0; x < W; x++)
        for (int y = 0; y < H; y++)
            f->h[x][y] = surf_f(x * h, y * h);
}

/* normal-angle bound: alpha = asin(eps/(1 - s eps)), eps = |grad error| */
static double bound_alpha_deg(double eps, double s) {
    double v = eps / (1.0 - s * eps);
    if (v >= 1.0) return 180.0;
    return asin(v) * (180.0 / PI);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Field f;
    field_init(&f, 1.0f);
    field_set(&f, 0, 0, 0.0f);
    field_set(&f, 1, 0, 1.0f);
    field_set(&f, 0, 1, 2.0f);
    field_set(&f, 1, 1, 3.0f);

    check("vertex_exact", fabsf(field_sample(&f, 0.0f, 0.0f)) < 1e-6f);
    check("bilinear_mid", fabsf(field_sample(&f, 0.5f, 0.5f) - 1.5f) < 1e-5f);
    check("bilinear_x", fabsf(field_sample(&f, 0.5f, 0.0f) - 0.5f) < 1e-5f);
    check("clamp_outside", field_sample(&f, -5.0f, -5.0f) == f.h[0][0]);

    check("slope_zero_flat", field_slope(&f, 8, 8) < 1e-6f);

    Field g;
    field_init(&g, 1.0f);
    for (int x = 0; x < W; x++) for (int y = 0; y < H; y++) g.h[x][y] = (float)x;
    check("slope_x_unit", fabsf(field_slope(&g, 5, 5) - 1.0f) < 1e-5f);

    /* ---- boundary normal validation on the analytic surface -------------- */
    const double h = 0.5;
    Field s;
    build_surface(&s, h);

    /* one-sided stencil component bounds: |e_x| <= h^2/3 M3x, etc. */
    double ebx = h * h / 3.0 * M3X;
    double eby = h * h / 3.0 * M3Y;
    double eps_bnd = sqrt(ebx * ebx + eby * eby);
    /* centered stencil component bounds: h^2/6 */
    double icx = h * h / 6.0 * M3X, icy = h * h / 6.0 * M3Y;
    double eps_in = sqrt(icx * icx + icy * icy);

    struct { const char* name; int x, y; } edges[4] = {
        {"west",  0,     7},
        {"east",  W - 1, 7},
        {"south", 7,     0},
        {"north", 7,     H - 1},
    };

    int all_new_inside = 1, all_old_outside = 1;
    printf("\nedge normal angles vs analytic normal (h=%.2f, eps_bnd=%.6f):\n",
           h, eps_bnd);
    for (int i = 0; i < 4; i++) {
        int x = edges[i].x, y = edges[i].y;
        double px = x * h, py = y * h, agx, agy;
        surf_grad(px, py, &agx, &agy);
        double smag = sqrt(agx * agx + agy * agy);
        double alpha = bound_alpha_deg(eps_bnd, smag);

        float na[3];
        grad_normal((float)agx, (float)agy, na);

        float gxo, gyo, no[3];
        field_grad_old(&s, x, y, &gxo, &gyo);
        grad_normal(gxo, gyo, no);
        float ang_old = angle_deg(na, no);

        float gxn, gyn, nn[3];
        field_grad(&s, x, y, &gxn, &gyn);
        grad_normal(gxn, gyn, nn);
        float ang_new = angle_deg(na, nn);

        double ex = (double)gxn - agx, ey = (double)gyn - agy;
        printf("  %-5s bound=%.4f deg  old_angle=%.4f  new_angle=%.4f"
               "  |grad_err|=%.2e\n",
               edges[i].name, alpha, ang_old, ang_new, sqrt(ex * ex + ey * ey));

        /* angles above the bound fail the task; slack absorbs float rounding */
        if (!(ang_new <= alpha + 5e-4)) all_new_inside = 0;
        if (!(ang_old > alpha)) all_old_outside = 0;
    }
    check("edge_new_angle_inside_bound_all4", all_new_inside);
    check("edge_old_angle_exceeds_bound_all4", all_old_outside);

    /* interior node must satisfy the tighter centered-difference bound */
    {
        int x = 7, y = 7;
        double agx, agy;
        surf_grad(x * h, y * h, &agx, &agy);
        float gx, gy;
        field_grad(&s, x, y, &gx, &gy);
        double ex = (double)gx - agx, ey = (double)gy - agy;
        double err = sqrt(ex * ex + ey * ey);
        double alpha_in = bound_alpha_deg(eps_in, sqrt(agx * agx + agy * agy));
        float na[3], nn[3];
        grad_normal((float)agx, (float)agy, na);
        grad_normal(gx, gy, nn);
        printf("\ninterior (7,7): |grad_err|=%.3e bound_eps=%.3e,"
               " angle=%.4f bound=%.4f deg\n",
               err, eps_in, angle_deg(na, nn), alpha_in);
        check("interior_inside_centered_bound", err <= eps_in + 1e-6);
    }

    /* convergence order: new stencil O(h^2), old one O(h), at the west edge */
    printf("\nconvergence at west edge (0,7):\n");
    double hs[3] = {1.0, 0.5, 0.25};
    double err_new[3], err_old[3];
    for (int k = 0; k < 3; k++) {
        Field t;
        build_surface(&t, hs[k]);
        double agx, agy;
        surf_grad(0.0, 7 * hs[k], &agx, &agy);
        float gxn, gyn, gxo, gyo;
        field_grad(&t, 0, 7, &gxn, &gyn);
        field_grad_old(&t, 0, 7, &gxo, &gyo);
        err_new[k] = hypot((double)gxn - agx, (double)gyn - agy);
        err_old[k] = hypot((double)gxo - agx, (double)gyo - agy);
        printf("  h=%.2f  old_err=%.6e  new_err=%.6e\n",
               hs[k], err_old[k], err_new[k]);
    }
    double ratio_new = err_new[0] / err_new[1];
    double ratio_new2 = err_new[1] / err_new[2];
    double ratio_old = err_old[0] / err_old[1];
    double ratio_old2 = err_old[1] / err_old[2];
    printf("  ratios new(h/2): %.3f %.3f (expect ~4), old: %.3f %.3f (expect ~2)\n",
           ratio_new, ratio_new2, ratio_old, ratio_old2);
    check("new_stencil_order_h2",
          ratio_new > 3.6 && ratio_new < 4.4 &&
          ratio_new2 > 3.6 && ratio_new2 < 4.4);
    check("old_stencil_order_h1",
          ratio_old > 1.85 && ratio_old < 2.15 &&
          ratio_old2 > 1.85 && ratio_old2 < 2.15);

    printf("\nall_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
