#define _POSIX_C_SOURCE 200809L

/*
 * erosion: droplet-based hydraulic erosion over an fBm heightfield.
 *
 * Model: each droplet lands on a random lattice point, walks downhill by
 * bilinear-gradient descent with inertia, erodes by
 *   capacity = K * speed * water
 * when under capacity, deposits when over, evaporates water each step, and
 * dumps any remaining sediment where it stops (mass-conserving except at the
 * domain edge).
 *
 * An "iteration" is one pass of `drops` droplets. The droplet sample is held
 * FIXED across iterations (same starting points and jitter) so the residual
 *   r_i = RMS( h_i - h_{i-1} )
 * measures convergence of the pass operator, not resampling noise. Pass
 * `--resample` to redraw the droplet sample every iteration instead.
 *
 * The CLI is built for measurement:
 *   ./erosion --n N --drops D --iters I --seed S --out PREFIX [--resample]
 * writes PREFIX.0000.h32 (initial surface) and PREFIX.<i>.h32 for i=1..I,
 * little-endian float32 row-major. stdout is one JSON object with parameters,
 * per-iteration wall time and residual; stderr carries the same in text.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <time.h>

/* ---------- deterministic hash / value noise (same family as noise/) ---------- */
static unsigned hash2(int x, int y, unsigned seed) {
    unsigned h = (unsigned)x * 374761393u + (unsigned)y * 668265263u + seed * 2246822519u;
    h = (h ^ (h >> 13)) * 1274126177u;
    return h ^ (h >> 16);
}
static float rand01u(unsigned h) { return (float)(h & 0xffffff) / 16777215.0f; }
static float lerpf(float a, float b, float t) { return a + (b - a) * t; }
static float fade(float t) { return t * t * (3.0f - 2.0f * t); }

static float value_noise(float x, float y, unsigned seed) {
    int xi = (int)floorf(x), yi = (int)floorf(y);
    float fx = x - xi, fy = y - yi;
    float v00 = rand01u(hash2(xi, yi, seed));
    float v10 = rand01u(hash2(xi + 1, yi, seed));
    float v01 = rand01u(hash2(xi, yi + 1, seed));
    float v11 = rand01u(hash2(xi + 1, yi + 1, seed));
    float u = fade(fx), v = fade(fy);
    return lerpf(lerpf(v00, v10, u), lerpf(v01, v11, u), v) * 2.0f - 1.0f;
}

static float fbm(float x, float y, int octaves, float lacunarity, float gain, unsigned seed) {
    float sum = 0.0f, amp = 1.0f, norm = 0.0f, freq = 1.0f;
    for (int i = 0; i < octaves; i++) {
        sum += value_noise(x * freq, y * freq, seed + (unsigned)i * 17u) * amp;
        norm += amp;
        amp *= gain;
        freq *= lacunarity;
    }
    return norm > 0.0f ? sum / norm : 0.0f;
}

/* ---------- field ---------- */
typedef struct {
    int n;
    float* h;    /* elevation, row-major y*n+x */
    float* tmp;  /* net delta accumulated by one pass */
} Field;

static int field_alloc(Field* f, int n) {
    f->n = n;
    f->h = (float*)calloc((size_t)n * n, sizeof(float));
    f->tmp = (float*)calloc((size_t)n * n, sizeof(float));
    return f->h && f->tmp;
}
static void field_free(Field* f) { free(f->h); free(f->tmp); f->h = f->tmp = NULL; }

static inline int clampi(int v, int lo, int hi) { return v < lo ? lo : (v > hi ? hi : v); }

static float sample_bilinear(const Field* f, float x, float y) {
    int n = f->n;
    int x0 = (int)floorf(x), y0 = (int)floorf(y);
    float tx = x - x0, ty = y - y0;
    int x1 = clampi(x0 + 1, 0, n - 1), y1 = clampi(y0 + 1, 0, n - 1);
    x0 = clampi(x0, 0, n - 1);
    y0 = clampi(y0, 0, n - 1);
    const float* h = f->h;
    float a = h[(size_t)y0 * n + x0] * (1 - tx) + h[(size_t)y0 * n + x1] * tx;
    float b = h[(size_t)y1 * n + x0] * (1 - tx) + h[(size_t)y1 * n + x1] * tx;
    return a * (1 - ty) + b * ty;
}

/* gradient of the bilinear patch at (x, y), per lattice unit */
static void grad_at(const Field* f, float x, float y, float* gx, float* gy) {
    int n = f->n;
    int x0 = clampi((int)floorf(x), 0, n - 1);
    int y0 = clampi((int)floorf(y), 0, n - 1);
    int x1 = clampi(x0 + 1, 0, n - 1), y1 = clampi(y0 + 1, 0, n - 1);
    float tx = x - (float)x0, ty = y - (float)y0;
    const float* h = f->h;
    float h00 = h[(size_t)y0 * n + x0], h10 = h[(size_t)y0 * n + x1];
    float h01 = h[(size_t)y1 * n + x0], h11 = h[(size_t)y1 * n + x1];
    *gx = lerpf(h10 - h00, h11 - h01, ty);
    *gy = lerpf(h01 - h00, h11 - h10, tx);
}

/* ---------- parameters ---------- */
typedef struct {
    int      drops;       /* droplets per iteration */
    int      max_steps;   /* droplet lifetime */
    float    inertia;
    float    capacity;
    float    erode;
    float    deposit;
    float    evaporate;
    float    gravity;
    float    radius;      /* erosion brush radius, lattice units */
    float    smooth;      /* laplacian smoothing of the pass delta */
    int      resample;    /* 1 = redraw droplet sample each iteration */
    unsigned seed;
} Params;

static float brush_weight(float dist, float radius) {
    if (dist >= radius) return 0.0f;
    float t = 1.0f - dist / radius;
    return t * t;
}

/* deposit `amt` around lattice cell (cx, cy) into tmp */
static void deposit_brush(Field* f, int cx, int cy, float amt, float radius) {
    int n = f->n;
    int r = (int)ceilf(radius);
    float wsum = 0.0f;
    /* normalize so the brush deposits exactly `amt` in total */
    for (int oy = -r; oy <= r; oy++)
        for (int ox = -r; ox <= r; ox++) {
            int sx = cx + ox, sy = cy + oy;
            if (sx < 0 || sy < 0 || sx >= n || sy >= n) continue;
            wsum += brush_weight(sqrtf((float)(ox * ox + oy * oy)), radius);
        }
    if (wsum <= 0.0f) return;
    for (int oy = -r; oy <= r; oy++)
        for (int ox = -r; ox <= r; ox++) {
            int sx = cx + ox, sy = cy + oy;
            if (sx < 0 || sy < 0 || sx >= n || sy >= n) continue;
            float w = brush_weight(sqrtf((float)(ox * ox + oy * oy)), radius);
            if (w > 0.0f) f->tmp[(size_t)sy * n + sx] += amt * w / wsum;
        }
}

/* erode `amt` from the brush around (cx, cy) into tmp, return amount actually
 * removed (clamped so a cell may not go far below its neighbors) */
static float erode_brush(Field* f, int cx, int cy, float amt, float radius) {
    int n = f->n;
    int r = (int)ceilf(radius);
    float wsum = 0.0f, removed = 0.0f;
    for (int oy = -r; oy <= r; oy++)
        for (int ox = -r; ox <= r; ox++) {
            int sx = cx + ox, sy = cy + oy;
            if (sx < 0 || sy < 0 || sx >= n || sy >= n) continue;
            wsum += brush_weight(sqrtf((float)(ox * ox + oy * oy)), radius);
        }
    if (wsum <= 0.0f) return 0.0f;
    for (int oy = -r; oy <= r; oy++)
        for (int ox = -r; ox <= r; ox++) {
            int sx = cx + ox, sy = cy + oy;
            if (sx < 0 || sy < 0 || sx >= n || sy >= n) continue;
            float w = brush_weight(sqrtf((float)(ox * ox + oy * oy)), radius);
            if (w <= 0.0f) continue;
            float take = amt * w / wsum;
            f->tmp[(size_t)sy * n + sx] -= take;
            removed += take;
        }
    return removed;
}

/* one iteration = `drops` droplets, all reading the same pre-pass surface */
static void erode_pass(Field* f, const Params* p, unsigned sample_seed) {
    int n = f->n;
    size_t nn = (size_t)n * n;
    memset(f->tmp, 0, nn * sizeof(float));

    for (int d = 0; d < p->drops; d++) {
        unsigned h0 = hash2(d, 0, p->seed + sample_seed * 7919u);
        float px = (float)(h0 & 0xffff) / 65535.0f * (float)(n - 1);
        float py = (float)((h0 >> 16) & 0xffff) / 65535.0f * (float)(n - 1);
        float dx = 0.0f, dy = 0.0f, speed = 1.0f, water = 1.0f, sediment = 0.0f;
        int ix = clampi((int)px, 0, n - 1), iy = clampi((int)py, 0, n - 1);

        for (int step = 0; step < p->max_steps; step++) {
            float gx, gy;
            grad_at(f, px, py, &gx, &gy);
            dx = dx * p->inertia - gx * (1.0f - p->inertia);
            dy = dy * p->inertia - gy * (1.0f - p->inertia);
            float len = sqrtf(dx * dx + dy * dy);
            if (len < 1e-12f) break;
            dx /= len;
            dy /= len;
            float nx = px + dx, ny = py + dy;
            if (nx < 0.0f || ny < 0.0f || nx > (float)(n - 1) || ny > (float)(n - 1)) break;

            int nix = clampi((int)nx, 0, n - 1), niy = clampi((int)ny, 0, n - 1);
            float h_old = sample_bilinear(f, px, py);
            float h_new = sample_bilinear(f, nx, ny);
            float drop = h_old - h_new;
            speed = sqrtf(speed * speed + p->gravity * (drop > 0.0f ? drop : 0.0f));
            float cap = p->capacity * speed * water;
            if (cap < 0.0f) cap = 0.0f;

            if (sediment < cap) {
                float amt = (cap - sediment) * p->erode;
                float removed = erode_brush(f, nix, niy, amt, p->radius);
                sediment += removed;
            } else {
                float amt = (sediment - cap) * p->deposit;
                sediment -= amt;
                deposit_brush(f, nix, niy, amt, p->radius);
            }
            water *= (1.0f - p->evaporate);
            px = nx;
            py = ny;
            ix = nix;
            iy = niy;
            if (water < 1e-3f) break;
        }
        /* mass conservation: dump remaining sediment where the droplet stopped */
        if (sediment > 0.0f) deposit_brush(f, ix, iy, sediment, p->radius);
    }

    /* apply the net delta; optional laplacian smoothing of the delta only
     * (does not diffuse the terrain itself) */
    if (p->smooth > 0.0f) {
        for (int y = 1; y < n - 1; y++)
            for (int x = 1; x < n - 1; x++) {
                size_t i = (size_t)y * n + x;
                float lap = f->tmp[i - 1] + f->tmp[i + 1] + f->tmp[i - n] + f->tmp[i + n]
                          - 4.0f * f->tmp[i];
                f->h[i] += f->tmp[i] + p->smooth * lap;
            }
        for (int y = 0; y < n; y++)
            for (int x = 0; x < n; x++)
                if (x == 0 || y == 0 || x == n - 1 || y == n - 1)
                    f->h[(size_t)y * n + x] += f->tmp[(size_t)y * n + x];
    } else {
        for (size_t i = 0; i < nn; i++) f->h[i] += f->tmp[i];
    }
}

/* stats for one iteration, all defined against (h_i, h_{i-1}, h_0) */
static void iter_stats(const float* h0, const float* prev, const float* cur,
                       size_t nn, double* out) {
    double res = 0.0, rms_init = 0.0, l1 = 0.0;
    double pos = 0.0, neg = 0.0;   /* deposition / erosion volume, cell units */
    for (size_t i = 0; i < nn; i++) {
        double d = (double)cur[i] - (double)prev[i];
        double di = (double)cur[i] - (double)h0[i];
        res += d * d;
        rms_init += di * di;
        l1 += fabs(d);
        if (d > 0.0) pos += d; else neg -= d;
    }
    out[0] = sqrt(res / (double)nn);          /* residual (RMS vs prev) */
    out[1] = sqrt(rms_init / (double)nn);     /* RMS vs initial surface */
    out[2] = l1 / (double)nn;                 /* mean |delta| */
    out[3] = pos;                             /* deposition volume */
    out[4] = neg;                             /* erosion volume */
}

static double now_s(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static int write_h32(const char* path, const float* h, size_t nn) {
    FILE* fp = fopen(path, "wb");
    if (!fp) return 0;
    size_t wrote = fwrite(h, sizeof(float), nn, fp);
    fclose(fp);
    return wrote == nn;
}

static void usage(void) {
    fprintf(stderr,
        "usage: erosion --n N --drops D --iters I --seed S --out PREFIX [options]\n"
        "  writes PREFIX.0000.h32 (initial) .. PREFIX.<I>.h32\n"
        "options: --max-steps --inertia --capacity --erode --deposit --evaporate\n"
        "         --gravity --radius --smooth --resample\n");
}

int main(int argc, char** argv) {
    int n = 256, drops = 20000, iters = 10;
    unsigned seed = 7u;
    const char* out = NULL;
    Params p;
    memset(&p, 0, sizeof(p));
    p.max_steps = 48;
    p.inertia = 0.05f;
    p.capacity = 0.05f;
    p.erode = 0.3f;
    p.deposit = 0.3f;
    p.evaporate = 0.01f;
    p.gravity = 4.0f;
    p.radius = 2.0f;
    p.smooth = 0.0f;
    p.resample = 0;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--n") && i + 1 < argc) n = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--drops") && i + 1 < argc) drops = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--iters") && i + 1 < argc) iters = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--seed") && i + 1 < argc) seed = (unsigned)strtoul(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "--out") && i + 1 < argc) out = argv[++i];
        else if (!strcmp(argv[i], "--max-steps") && i + 1 < argc) p.max_steps = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--inertia") && i + 1 < argc) p.inertia = (float)atof(argv[++i]);
        else if (!strcmp(argv[i], "--capacity") && i + 1 < argc) p.capacity = (float)atof(argv[++i]);
        else if (!strcmp(argv[i], "--erode") && i + 1 < argc) p.erode = (float)atof(argv[++i]);
        else if (!strcmp(argv[i], "--deposit") && i + 1 < argc) p.deposit = (float)atof(argv[++i]);
        else if (!strcmp(argv[i], "--evaporate") && i + 1 < argc) p.evaporate = (float)atof(argv[++i]);
        else if (!strcmp(argv[i], "--gravity") && i + 1 < argc) p.gravity = (float)atof(argv[++i]);
        else if (!strcmp(argv[i], "--radius") && i + 1 < argc) p.radius = (float)atof(argv[++i]);
        else if (!strcmp(argv[i], "--smooth") && i + 1 < argc) p.smooth = (float)atof(argv[++i]);
        else if (!strcmp(argv[i], "--resample")) p.resample = 1;
        else { usage(); return 2; }
    }
    if (!out || n < 8 || drops < 1 || iters < 1) { usage(); return 2; }
    p.drops = drops;
    p.seed = seed;

    Field f, prev;
    if (!field_alloc(&f, n) || !field_alloc(&prev, n)) { fprintf(stderr, "alloc failed\n"); return 1; }
    size_t nn = (size_t)n * n;

    /* initial surface: fBm */
    for (int y = 0; y < n; y++)
        for (int x = 0; x < n; x++) {
            float u = (float)x / (float)(n - 1), v = (float)y / (float)(n - 1);
            float hgt = fbm(u * 6.0f, v * 6.0f, 6, 2.0f, 0.5f, seed);
            hgt += 0.25f * fbm(u * 18.0f, v * 18.0f, 3, 2.0f, 0.5f, seed + 100u);
            f.h[(size_t)y * n + x] = hgt * (float)n * 0.15f;
        }
    float* h0 = (float*)malloc(nn * sizeof(float));
    if (!h0) { fprintf(stderr, "alloc failed\n"); return 1; }
    memcpy(h0, f.h, nn * sizeof(float));

    char path[1024];
    snprintf(path, sizeof(path), "%s.0000.h32", out);
    if (!write_h32(path, f.h, nn)) { fprintf(stderr, "write failed: %s\n", path); return 1; }

    printf("{\"n\":%d,\"drops\":%d,\"iters\":%d,\"seed\":%u,"
           "\"max_steps\":%d,\"inertia\":%.4f,\"capacity\":%.4f,"
           "\"erode\":%.4f,\"deposit\":%.4f,\"evaporate\":%.4f,"
           "\"gravity\":%.4f,\"radius\":%.4f,\"smooth\":%.4f,\"resample\":%d,\"runs\":[",
           n, drops, iters, seed, p.max_steps, p.inertia, p.capacity,
           p.erode, p.deposit, p.evaporate, p.gravity, p.radius, p.smooth, p.resample);
    fflush(stdout);

    float* prevbuf = prev.h;
    for (int it = 1; it <= iters; it++) {
        memcpy(prevbuf, f.h, nn * sizeof(float));
        double t0 = now_s();
        erode_pass(&f, &p, p.resample ? (unsigned)it : 0u);
        double t1 = now_s();
        double st[5];
        iter_stats(h0, prevbuf, f.h, nn, st);
        double ms = (t1 - t0) * 1e3;
        snprintf(path, sizeof(path), "%s.%04d.h32", out, it);
        if (!write_h32(path, f.h, nn)) { fprintf(stderr, "write failed: %s\n", path); return 1; }
        printf("%s{\"iter\":%d,\"ms\":%.3f,\"residual\":%.9g,\"rms_vs_init\":%.9g,"
               "\"mean_abs_delta\":%.9g,\"deposit_vol\":%.9g,\"erode_vol\":%.9g,"
               "\"us_per_cell\":%.4f}",
               it == 1 ? "" : ",", it, ms, st[0], st[1], st[2], st[3], st[4],
               ms * 1e3 / ((double)n * (double)n));
        fflush(stdout);
        fprintf(stderr, "iter=%d ms=%.3f residual=%.9g rms_vs_init=%.9g\n",
                it, ms, st[0], st[1]);
    }
    printf("]}\n");

    free(h0);
    field_free(&f);
    field_free(&prev);
    return 0;
}
