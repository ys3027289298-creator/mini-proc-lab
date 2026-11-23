
/* noise: 2D value noise + fBm octave accumulation, deterministic hash. */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "noise.h"

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* range and determinism */
    float lo = 1e9f, hi = -1e9f;
    for (int i = 0; i < 200; i++) {
        float v = value_noise(i * 0.37f, i * 0.11f, 7u);
        if (v < lo) lo = v;
        if (v > hi) hi = v;
    }
    check("range", lo >= -1.01f && hi <= 1.01f);
    check("deterministic", value_noise(1.3f, 2.7f, 7u) == value_noise(1.3f, 2.7f, 7u));
    check("seed_changes", value_noise(1.3f, 2.7f, 7u) != value_noise(1.3f, 2.7f, 8u));

    /* lattice points must vanish (value noise is 0 at integer grid after remap? check continuity instead) */
    float a = value_noise(0.5f, 0.5f, 3u);
    float b = value_noise(0.5f + 1e-4f, 0.5f, 3u);
    check("continuous", fabsf(a - b) < 0.01f);

    float f1 = fbm(0.3f, 0.4f, 4, 2.0f, 0.5f, 11u);
    float f2 = fbm(0.3f, 0.4f, 4, 2.0f, 0.5f, 11u);
    check("fbm_deterministic", f1 == f2);
    check("fbm_range", f1 >= -1.01f && f1 <= 1.01f);
    check("fbm_zero_octaves", fbm(0.3f, 0.4f, 0, 2.0f, 0.5f, 11u) == 0.0f);

    /* more octaves => higher variance of first differences */
    float d1 = 0.0f, d2 = 0.0f;
    for (int i = 0; i < 100; i++) {
        float x = i * 0.05f;
        d1 += fabsf(fbm(x, 0.5f, 1, 2.0f, 0.5f, 5u) - fbm(x + 0.01f, 0.5f, 1, 2.0f, 0.5f, 5u));
        d2 += fabsf(fbm(x, 0.5f, 5, 2.0f, 0.5f, 5u) - fbm(x + 0.01f, 0.5f, 5, 2.0f, 0.5f, 5u));
    }
    check("more_octaves_rougher", d2 > d1);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
