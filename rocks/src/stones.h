#ifndef ROCKS_STONES_H
#define ROCKS_STONES_H

/* Two independent texture modules over one shared stone-layout description.
 *
 * The layout (position, radius, rotation, vertex jitter, palette role) is
 * generated from a FIXED layout seed and is identical for every texture seed
 * and for both modules. Only the surface grain is seeded, so a luminance
 * difference between two runs at the same coordinate is a pure texture
 * difference - not a differently-placed pebble.
 */

#define MAX_STONES 512

typedef struct {
    float x, y;     /* centre, texture-space units */
    float rx, ry;   /* radii */
    float rot;      /* radians */
    float e[8];     /* unit-circle radial jitter per vertex */
    float tone;     /* 0..1 pebble fill level */
    float rough;    /* 0..1 grain contrast */
    float warp;     /* 0..1 internal marbling amount */
    int variant;    /* palette step index */
    float bound;    /* conservative max extent from centre */
} Stone;

typedef struct {
    Stone s[MAX_STONES];
    int n;
} StoneSet;

/* Deterministic stone layout. layout_seed is fixed by the caller. */
void stones_layout(StoneSet *set, unsigned layout_seed);

/* Sample stone coverage (0..1) and stone id at a texture-space point. */
float stones_coverage(const StoneSet *set, float x, float y, int *id_out);

/* ------------------------------------------------------------------ */
/* Module A: value-noise fBm grain                                     */
/* Module B: gradient(Perlin-style) noise fBm grain                    */
/* Both take (x, y, seed) and return roughly [-1, 1].                  */
/* ------------------------------------------------------------------ */

float textureA_raw(float x, float y, unsigned seed);
float textureB_raw(float x, float y, unsigned seed);

/* Full stone-surface luminance in [0,1] for a given module at (x,y). */
float textureA_lum(const StoneSet *set, float x, float y, unsigned seed);
float textureB_lum(const StoneSet *set, float x, float y, unsigned seed);

#endif
