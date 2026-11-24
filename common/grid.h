#ifndef SP_GRID_H
#define SP_GRID_H

#include <stddef.h>

/* Float grid positioned in world space: (ox, oy) is the world coordinate of
 * v[0].  Full-world layers use ox = oy = 0; a windowed run keeps the world
 * position so noise and slope stay global. */
typedef struct {
    int ox, oy;
    int w, h;
    float *v;
} SpGrid;

int sp_grid_alloc(SpGrid *g, int ox, int oy, int w, int h);
void sp_grid_free(SpGrid *g);
void sp_grid_fill(SpGrid *g, float value);

float sp_grid_at(const SpGrid *g, int lx, int ly);
float sp_grid_at_clamped(const SpGrid *g, int lx, int ly);
float sp_grid_sample(const SpGrid *g, float wx, float wy);

#endif
