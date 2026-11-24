#ifndef SP_RNG_H
#define SP_RNG_H

#include <stdint.h>

/* Hash of two signed lattice coordinates and a 64-bit seed.  Every layer
 * draws from world coordinates through this function, never from a mutable
 * stream or a tile-local counter. */
uint64_t sp_hash2d(int x, int y, uint64_t seed);

/* Uniform in [0, 1), 24 significant bits so the value is exact in float. */
float sp_unit(uint64_t h);

/* Signed value noise in [-1, 1], C1 continuous. */
float sp_value_noise(float x, float y, uint64_t seed);

/* Fractal sum of value noise; octave seeds are derived from `seed`. */
float sp_fbm(float x, float y, int octaves, float lacunarity, float gain,
             uint64_t seed);

float sp_lerp(float a, float b, float t);
float sp_clampf(float v, float lo, float hi);

#endif
