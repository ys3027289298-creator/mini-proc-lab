#include <math.h>
#include <stdio.h>
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
        float v = noise_value2(i * 0.37f, i * 0.11f, 7u);
        if (v < lo) lo = v;
        if (v > hi) hi = v;
    }
    check("range", lo >= -1.01f && hi <= 1.01f);
    check("deterministic", noise_value2(1.3f, 2.7f, 7u) == noise_value2(1.3f, 2.7f, 7u));
    check("seed_changes", noise_value2(1.3f, 2.7f, 7u) != noise_value2(1.3f, 2.7f, 8u));

    float a = noise_value2(0.5f, 0.5f, 3u);
    float b = noise_value2(0.5f + 1e-4f, 0.5f, 3u);
    check("continuous", fabsf(a - b) < 0.01f);

    /* white noise must stay in [0, 1) and decorrelate across lattice points */
    int white_ok = 1;
    for (int y = 0; y < 16; y++)
        for (int x = 0; x < 16; x++) {
            float v = noise_white01(x, y, 4u);
            if (!(v >= 0.0f && v < 1.0f)) white_ok = 0;
        }
    check("white_range", white_ok);
    check("white_distinct", noise_white01(1, 2, 4u) != noise_white01(2, 1, 4u));

    float f1 = noise_fbm2(0.3f, 0.4f, 4, 2.0f, 0.5f, 11u);
    float f2 = noise_fbm2(0.3f, 0.4f, 4, 2.0f, 0.5f, 11u);
    check("fbm_deterministic", f1 == f2);
    check("fbm_range", f1 >= -1.01f && f1 <= 1.01f);
    check("fbm_zero_octaves", noise_fbm2(0.3f, 0.4f, 0, 2.0f, 0.5f, 11u) == 0.0f);

    /* more octaves => higher variance of first differences */
    float d1 = 0.0f, d2 = 0.0f;
    for (int i = 0; i < 100; i++) {
        float x = i * 0.05f;
        d1 += fabsf(noise_fbm2(x, 0.5f, 1, 2.0f, 0.5f, 5u) - noise_fbm2(x + 0.01f, 0.5f, 1, 2.0f, 0.5f, 5u));
        d2 += fabsf(noise_fbm2(x, 0.5f, 5, 2.0f, 0.5f, 5u) - noise_fbm2(x + 0.01f, 0.5f, 5, 2.0f, 0.5f, 5u));
    }
    check("more_octaves_rougher", d2 > d1);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
