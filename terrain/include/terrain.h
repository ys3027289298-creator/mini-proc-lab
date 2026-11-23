/* terrain: heightfield sampling with slope and bilinear interpolation. */
#ifndef TERRAIN_H
#define TERRAIN_H

#define TERRAIN_MAX 256

/* Fixed-capacity heightfield with bilinear sampling and central-difference
 * slope. `wc`/`hc` give the live extent (<= TERRAIN_MAX). */
typedef struct {
    float h[TERRAIN_MAX][TERRAIN_MAX];
    int wc, hc;
    float cell; /* world units per grid spacing */
} Field;

void field_init(Field* f, int w, int h, float cell);

/* Set a vertex; out-of-range coordinates are ignored. */
void field_set(Field* f, int x, int y, float v);

/* Bilinear world-space lookup, clamped to the field bounds. */
float field_sample(const Field* f, float wx, float wy);

/* Slope magnitude at a vertex, in height units per world unit. One-sided at
 * the border. */
float field_slope(const Field* f, int x, int y);

/* Slope magnitude of the bilinear surface at a world position, in height units
 * per world unit. */
float field_slope_at(const Field* f, float wx, float wy);

/* Fill an w x h field with fBm scaled to `amplitude` height units, with the
 * base octave having wavelength 1/base_freq world units. Deterministic for a
 * given seed. */
void field_generate_fbm(Field* f, int w, int h, float cell, float amplitude,
                        int octaves, float base_freq, unsigned seed);

#endif
