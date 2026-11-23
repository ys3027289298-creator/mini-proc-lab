
/* terrain: heightfield sampling with slope and bilinear interpolation. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define W 16
#define H 16

typedef struct { float h[W][H]; float cell; } Field;

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

static float field_slope(const Field* f, int x, int y) {
    int xm = x > 0 ? x - 1 : x, xp = x < W - 1 ? x + 1 : x;
    int ym = y > 0 ? y - 1 : y, yp = y < H - 1 ? y + 1 : y;
    float dx = (f->h[xp][y] - f->h[xm][y]) / ((xp - xm) * f->cell);
    float dy = (f->h[x][yp] - f->h[x][ym]) / ((yp - ym) * f->cell);
    return sqrtf(dx * dx + dy * dy);
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

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
