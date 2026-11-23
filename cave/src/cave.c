/* cave: 3D fBm density carving + connected-component (island) analysis.
 *
 * Usage:
 *   ./cave            run self tests
 *   ./cave <t>        carve with threshold t, print connectivity stats
 *   ./cave sweep      scan thresholds, print volume/connectivity per step
 */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define NX 64
#define NY 64
#define NZ 64
#define NCELLS (NX * NY * NZ)

#define IDX(x, y, z) (((z) * NY + (y)) * NX + (x))

static unsigned hash3(int x, int y, int z, unsigned seed) {
    unsigned h = (unsigned)x * 374761393u + (unsigned)y * 668265263u
               + (unsigned)z * 1442695041u + seed * 2246822519u;
    h = (h ^ (h >> 13)) * 1274126177u;
    return h ^ (h >> 16);
}

static float rand01(unsigned h) { return (float)(h & 0xffffff) / 16777215.0f; }

static float lerp(float a, float b, float t) { return a + (b - a) * t; }
static float fade(float t) { return t * t * (3.0f - 2.0f * t); }

static float value_noise3(float x, float y, float z, unsigned seed) {
    int xi = (int)floorf(x), yi = (int)floorf(y), zi = (int)floorf(z);
    float fx = x - xi, fy = y - yi, fz = z - zi;
    float c[2][2][2];
    for (int dz = 0; dz < 2; dz++)
        for (int dy = 0; dy < 2; dy++)
            for (int dx = 0; dx < 2; dx++)
                c[dz][dy][dx] = rand01(hash3(xi + dx, yi + dy, zi + dz, seed));
    float u = fade(fx), v = fade(fy), w = fade(fz);
    float x00 = lerp(c[0][0][0], c[0][0][1], u);
    float x10 = lerp(c[0][1][0], c[0][1][1], u);
    float x01 = lerp(c[1][0][0], c[1][0][1], u);
    float x11 = lerp(c[1][1][0], c[1][1][1], u);
    float y0 = lerp(x00, x10, v), y1 = lerp(x01, x11, v);
    return lerp(y0, y1, w) * 2.0f - 1.0f;
}

static float fbm3(float x, float y, float z, int octaves, float lacunarity, float gain, unsigned seed) {
    if (octaves <= 0) return 0.0f;
    float sum = 0.0f, amp = 1.0f, norm = 0.0f, freq = 1.0f;
    for (int i = 0; i < octaves; i++) {
        sum += value_noise3(x * freq, y * freq, z * freq, seed + (unsigned)i * 17u) * amp;
        norm += amp;
        amp *= gain;
        freq *= lacunarity;
    }
    return norm > 0.0f ? sum / norm : 0.0f;
}

/* world state */
static float density[NCELLS];
static unsigned char cave[NCELLS];
static int label[NCELLS];
static int queue_[NCELLS];
static int comp_size[NCELLS]; /* sizes indexed by component id - 1 */
static int ncomp = 0;

#define NOISE_SCALE 0.08f
#define OCTAVES 4
#define SEED 42u

static void density_init(void) {
    for (int z = 0; z < NZ; z++)
        for (int y = 0; y < NY; y++)
            for (int x = 0; x < NX; x++)
                density[IDX(x, y, z)] = fbm3(x * NOISE_SCALE, y * NOISE_SCALE, z * NOISE_SCALE,
                                             OCTAVES, 2.0f, 0.5f, SEED);
}

static int carve(float t) {
    int n = 0;
    for (int i = 0; i < NCELLS; i++) {
        cave[i] = density[i] > t;
        n += cave[i];
    }
    return n;
}

/* 6-connectivity BFS labeling. O(NCELLS) time, O(NCELLS) memory. */
static int label_components(void) {
    static const int d[6][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
    memset(label, 0, sizeof(label));
    ncomp = 0;
    for (int start = 0; start < NCELLS; start++) {
        if (!cave[start] || label[start]) continue;
        ncomp++;
        int head = 0, tail = 0, size = 0;
        queue_[tail++] = start;
        label[start] = ncomp;
        while (head < tail) {
            int cur = queue_[head++];
            size++;
            int cx = cur % NX, cy = (cur / NX) % NY, cz = cur / (NX * NY);
            for (int k = 0; k < 6; k++) {
                int nx = cx + d[k][0], ny = cy + d[k][1], nz = cz + d[k][2];
                if (nx < 0 || ny < 0 || nz < 0 || nx >= NX || ny >= NY || nz >= NZ) continue;
                int ni = IDX(nx, ny, nz);
                if (cave[ni] && !label[ni]) {
                    label[ni] = ncomp;
                    queue_[tail++] = ni;
                }
            }
        }
        comp_size[ncomp - 1] = size;
    }
    return ncomp;
}

static double now_ms(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1000.0 + ts.tv_nsec / 1e6;
}

static int cmp_desc(const void* a, const void* b) { return *(const int*)b - *(const int*)a; }

static void run_once(float t) {
    density_init();
    int carved = carve(t);
    double t0 = now_ms();
    int c = label_components();
    double label_ms = now_ms() - t0;

    int largest = 0;
    for (int i = 0; i < c; i++) if (comp_size[i] > largest) largest = comp_size[i];
    int islands = c > 0 ? c - 1 : 0;
    int island_cells = carved - largest;

    qsort(comp_size, c, sizeof(int), cmp_desc);
    printf("t=%.3f carved=%d p=%.4f components=%d largest=%d R=%.4f islands=%d island_cells=%d label_ms=%.2f\n",
           t, carved, NCELLS ? (float)carved / NCELLS : 0.0f, c, largest,
           carved ? (float)largest / carved : 0.0f, islands, island_cells, label_ms);
    printf("island_sizes_top5=");
    for (int i = 1; i < c && i <= 5; i++) printf("%s%d", i > 1 ? "," : "", comp_size[i]);
    printf("%s", c <= 1 ? "(none)" : "");
    printf("\n");
}

static void sweep(void) {
    density_init();
    printf("t carved p components largest R islands\n");
    for (float t = -0.2f; t <= 0.8001f; t += 0.05f) {
        int carved = carve(t);
        int c = label_components();
        int largest = 0;
        for (int i = 0; i < c; i++) if (comp_size[i] > largest) largest = comp_size[i];
        printf("%.2f %d %.4f %d %d %.4f %d\n",
               t, carved, (float)carved / NCELLS, c, largest,
               carved ? (float)largest / carved : 0.0f, c > 0 ? c - 1 : 0);
    }
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(int argc, char** argv) {
    if (argc > 1 && strcmp(argv[1], "sweep") == 0) {
        sweep();
        return 0;
    }
    if (argc > 1) {
        run_once((float)atof(argv[1]));
        return 0;
    }

    /* self tests */
    check("fbm3_deterministic", fbm3(0.3f, 0.4f, 0.5f, 4, 2.0f, 0.5f, 11u) == fbm3(0.3f, 0.4f, 0.5f, 4, 2.0f, 0.5f, 11u));
    check("fbm3_range", fabsf(fbm3(0.3f, 0.4f, 0.5f, 4, 2.0f, 0.5f, 11u)) <= 1.01f);

    density_init();

    carve(-2.0f); /* everything carved -> single component */
    check("full_carve_one_component", label_components() == 1 && comp_size[0] == NCELLS);

    carve(2.0f); /* nothing carved -> no components */
    check("empty_carve_zero_components", label_components() == 0);

    /* two isolated single-cell caves -> 2 components, 1 island */
    memset(cave, 0, sizeof(cave));
    cave[IDX(1, 1, 1)] = 1;
    cave[IDX(NX - 2, NY - 2, NZ - 2)] = 1;
    check("two_cells_two_components", label_components() == 2);
    check("two_cells_unit_sizes", comp_size[0] == 1 && comp_size[1] == 1);

    /* adjacent cells merge into one component */
    memset(cave, 0, sizeof(cave));
    cave[IDX(5, 5, 5)] = 1;
    cave[IDX(6, 5, 5)] = 1;
    cave[IDX(7, 6, 5)] = 1; /* diagonal only: not 6-connected */
    check("adjacent_merges", label_components() == 2 && (comp_size[0] == 2 || comp_size[1] == 2));

    /* labeling cost scales linearly: full grid must label in well under 1s */
    carve(-2.0f);
    double t0 = now_ms();
    label_components();
    check("label_full_fast", now_ms() - t0 < 1000.0);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
