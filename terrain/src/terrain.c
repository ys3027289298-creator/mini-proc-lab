
#include <math.h>
#include <string.h>
#include "noise.h"
#include "terrain.h"

void field_init(Field* f, int w, int h, float cell) {
    memset(f, 0, sizeof(*f));
    f->wc = w;
    f->hc = h;
    f->cell = cell;
}

void field_set(Field* f, int x, int y, float v) {
    if (x < 0 || y < 0 || x >= f->wc || y >= f->hc) return;
    f->h[x][y] = v;
}

float field_sample(const Field* f, float wx, float wy) {
    float gx = wx / f->cell, gy = wy / f->cell;
    int x0 = (int)floorf(gx), y0 = (int)floorf(gy);
    float tx = gx - x0, ty = gy - y0;
    int x1 = x0 + 1, y1 = y0 + 1;
    if (x0 < 0) x0 = 0;
    if (y0 < 0) y0 = 0;
    if (x1 >= f->wc) x1 = f->wc - 1;
    if (y1 >= f->hc) y1 = f->hc - 1;
    float a = f->h[x0][y0] * (1 - tx) + f->h[x1][y0] * tx;
    float b = f->h[x0][y1] * (1 - tx) + f->h[x1][y1] * tx;
    return a * (1 - ty) + b * ty;
}

float field_slope(const Field* f, int x, int y) {
    int xm = x > 0 ? x - 1 : x, xp = x < f->wc - 1 ? x + 1 : x;
    int ym = y > 0 ? y - 1 : y, yp = y < f->hc - 1 ? y + 1 : y;
    float dx = (f->h[xp][y] - f->h[xm][y]) / ((xp - xm) * f->cell);
    float dy = (f->h[x][yp] - f->h[x][ym]) / ((yp - ym) * f->cell);
    return sqrtf(dx * dx + dy * dy);
}

float field_slope_at(const Field* f, float wx, float wy) {
    float gx = wx / f->cell, gy = wy / f->cell;
    int x0 = (int)floorf(gx), y0 = (int)floorf(gy);
    if (x0 < 0) x0 = 0;
    if (y0 < 0) y0 = 0;
    if (x0 > f->wc - 2) x0 = f->wc - 2;
    if (y0 > f->hc - 2) y0 = f->hc - 2;
    int x1 = x0 + 1, y1 = y0 + 1;
    /* Mean gradient over the four cell edges. */
    float dx = ((f->h[x1][y0] - f->h[x0][y0]) + (f->h[x1][y1] - f->h[x0][y1])) * 0.5f;
    float dy = ((f->h[x0][y1] - f->h[x0][y0]) + (f->h[x1][y1] - f->h[x1][y0])) * 0.5f;
    dx /= f->cell;
    dy /= f->cell;
    return sqrtf(dx * dx + dy * dy);
}

void field_generate_fbm(Field* f, int w, int h, float cell, float amplitude,
                        int octaves, float base_freq, unsigned seed) {
    field_init(f, w, h, cell);
    for (int x = 0; x < w; x++)
        for (int y = 0; y < h; y++) {
            float n = noise_fbm2(x * cell * base_freq, y * cell * base_freq,
                                 octaves, 2.0f, 0.5f, seed);
            f->h[x][y] = n * amplitude;
        }
}
