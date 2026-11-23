/* Internal invariants of the scatter RNG. Includes the source to reach the
 * static hash, since the property under test is not observable through the
 * public API. */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "../src/scatter.c"

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static int cmp_u32(const void* a, const void* b) {
    unsigned x = *(const unsigned*)a, y = *(const unsigned*)b;
    return (x > y) - (x < y);
}

int main(void) {
    enum { N = 512 };
    static unsigned keys[512 * 512];
    int n = 0;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) keys[n++] = key2(i, j, 2024u);

    qsort(keys, (size_t)n, sizeof(*keys), cmp_u32);
    int unique = 0;
    for (int i = 0; i < n; i++)
        if (i == 0 || keys[i] != keys[i - 1]) unique++;
    /* Only birthday-paradox duplicates are acceptable: the old key
     * i*131 + j*17 + seed*7919 yielded ~75k distinct values for 262k cells. */
    check("keys_essentially_unique", unique > n - n / 1000);

    /* The jitter and acceptance draws must be uncorrelated. */
    double sxy = 0, sx = 0, sy = 0, sxx = 0, syy = 0, m = 0;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            unsigned k = key2(i, j, 2024u);
            double a = u01(k), b = u01(k + 2);
            sxy += a * b; sx += a; sy += b; sxx += a * a; syy += b * b; m += 1;
        }
    double cov = sxy / m - (sx / m) * (sy / m);
    double va = sxx / m - (sx / m) * (sx / m);
    double vb = syy / m - (sy / m) * (sy / m);
    check("draws_uncorrelated", fabs(cov / sqrt(va * vb)) < 0.01);

    /* Uniformity of the acceptance draw, 100 equal bins. */
    {
        enum { B = 100 };
        static int hist[B];
        memset(hist, 0, sizeof(hist));
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++) {
                unsigned k = key2(i, j, 2024u);
                int b = (int)(u01(k + 2) * B);
                if (b < 0) b = 0;
                if (b >= B) b = B - 1;
                hist[b]++;
            }
        double e = (double)n / B, chi = 0;
        for (int b = 0; b < B; b++) chi += (hist[b] - e) * (hist[b] - e) / e;
        /* chi2(99) 99.9th percentile is ~148 */
        check("draw_uniform", chi < 150.0);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
