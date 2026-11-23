/* noise: 2D value noise + fBm octave accumulation, deterministic hash. */
#ifndef NOISE_H
#define NOISE_H

/* Hash-based white noise in [0, 1). Deterministic for a given (x, y, seed). */
float noise_white01(int x, int y, unsigned seed);

/* 2D value noise with smoothstep fade, remapped to [-1, 1]. */
float noise_value2(float x, float y, unsigned seed);

/* fBm: octaves summed with amplitude gain and frequency lacunarity, normalised
 * to [-1, 1]. Returns 0 for octaves <= 0. */
float noise_fbm2(float x, float y, int octaves, float lacunarity, float gain,
                 unsigned seed);

#endif
