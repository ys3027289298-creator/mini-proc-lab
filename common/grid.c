#include "grid.h"

#include <math.h>
#include <stdlib.h>

int sp_grid_alloc(SpGrid *g, int ox, int oy, int w, int h) {
    g->ox = ox;
    g->oy = oy;
    g->w = w;
    g->h = h;
    g->v = (float *)calloc((size_t)w * (size_t)h, sizeof(float));
    return g->v != NULL;
}

void sp_grid_free(SpGrid *g) {
    free(g->v);
    g->v = NULL;
}

void sp_grid_fill(SpGrid *g, float value) {
    size_t n = (size_t)g->w * (size_t)g->h;
    for (size_t i = 0; i < n; i++) g->v[i] = value;
}

float sp_grid_at(const SpGrid *g, int lx, int ly) {
    return g->v[(size_t)ly * (size_t)g->w + (size_t)lx];
}

float sp_grid_at_clamped(const SpGrid *g, int lx, int ly) {
    if (lx < 0) lx = 0;
    if (ly < 0) ly = 0;
    if (lx >= g->w) lx = g->w - 1;
    if (ly >= g->h) ly = g->h - 1;
    return sp_grid_at(g, lx, ly);
}

float sp_grid_sample(const SpGrid *g, float wx, float wy) {
    float gx = wx - (float)g->ox;
    float gy = wy - (float)g->oy;
    int x0 = (int)floorf(gx), y0 = (int)floorf(gy);
    float tx = gx - (float)x0, ty = gy - (float)y0;
    float a = sp_grid_at_clamped(g, x0, y0) * (1.0f - tx) +
              sp_grid_at_clamped(g, x0 + 1, y0) * tx;
    float b = sp_grid_at_clamped(g, x0, y0 + 1) * (1.0f - tx) +
              sp_grid_at_clamped(g, x0 + 1, y0 + 1) * tx;
    return a * (1.0f - ty) + b * ty;
}
