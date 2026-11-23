#include <math.h>
#include <stdio.h>
#include "terrain.h"

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Field f;
    field_init(&f, 16, 16, 1.0f);
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
    field_init(&g, 16, 16, 1.0f);
    for (int x = 0; x < 16; x++) for (int y = 0; y < 16; y++) g.h[x][y] = (float)x;
    check("slope_x_unit", fabsf(field_slope(&g, 5, 5) - 1.0f) < 1e-5f);
    check("slope_at_matches", fabsf(field_slope_at(&g, 5.4f, 6.3f) - field_slope(&g, 5, 6)) < 1e-5f);

    Field gen;
    field_generate_fbm(&gen, 64, 64, 1.0f, 10.0f, 4, 0.25f, 3u);
    float lo = 1e9f, hi = -1e9f;
    for (int x = 0; x < 64; x++)
        for (int y = 0; y < 64; y++) {
            if (gen.h[x][y] < lo) lo = gen.h[x][y];
            if (gen.h[x][y] > hi) hi = gen.h[x][y];
        }
    check("gen_range", lo >= -10.01f && hi <= 10.01f);
    check("gen_varies", hi - lo > 0.5f);
    check("gen_deterministic", gen.h[3][4] == ({
        Field other; field_generate_fbm(&other, 64, 64, 1.0f, 10.0f, 4, 0.25f, 3u);
        other.h[3][4]; }));

    Field wide;
    field_generate_fbm(&wide, 32, 48, 0.5f, 10.0f, 4, 0.25f, 3u);
    check("gen_dims", wide.wc == 32 && wide.hc == 48 && wide.cell == 0.5f);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
