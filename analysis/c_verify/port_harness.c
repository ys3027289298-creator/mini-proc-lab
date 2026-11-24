/* Expose value_noise()/fbm() of noise/src/noise.c logic for port verification.
   The functions there are static, so the original source is #included directly. */
#include <stdio.h>
#include <math.h>
#define main noise_c_main_unused
#include "../../noise/src/noise.c"
#undef main

int main(void) {
    printf("value\n");
    for (int k = 0; k < 40; k++) {
        float x = (k * 0.137f) - 2.5f;
        float y = (k * 0.211f) - 3.1f;
        printf("%.8f %.8f %.8f\n", x, y, value_noise(x, y, 7u));
    }
    printf("fbm\n");
    int octs[] = {1, 2, 3, 5, 8};
    for (unsigned k = 0; k < sizeof(octs) / sizeof(octs[0]); k++) {
        int n = octs[k];
        for (int j = 0; j < 10; j++) {
            float x = j * 0.31f - 1.2f;
            float y = j * 0.17f + 0.4f;
            printf("%d %.8f %.8f %.8f\n", n, x, y,
                   fbm(x, y, n, 2.0f, 0.5f, 11u));
        }
    }
    return 0;
}
