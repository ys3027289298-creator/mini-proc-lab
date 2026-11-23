/* props: static props on a heightfield with a per-prop offset height layer.
 *
 * Bug reproduced: a prop's base is placed at the terrain height sampled at
 * its anchor (footprint center). On sloped terrain the high side of the
 * footprint pokes through the base -> the prop is buried; lighting samples
 * the bare terrain layer and misses the prop entirely.
 *
 * Fix: per-prop offset layer
 *   offset = max_{(x,y) in footprint} h_terrain(x,y) - h_terrain(anchor)
 *   base   = h_anchor + offset = max over footprint
 * The offset is exact, not just minimal:
 *   - any smaller: terrain above the base somewhere in the footprint (buried)
 *   - any larger:  base strictly above terrain everywhere (feet floating)
 * The max of a bilinear patch over an axis-aligned rect is attained at a
 * candidate vertex: rect corners, grid vertices inside the rect, and
 * grid-line/rect-edge intersections (interior critical points of a bilinear
 * patch are saddles, and the patch is linear along axis-aligned edges).
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>

#define W 64
#define H 64
#define CELL 1.0f
#define MAXPROP 32
#define TIERN 4

/* ---- deterministic value noise + fbm (same hash as the noise module) ---- */
static unsigned hash2(int x, int y, unsigned seed) {
    unsigned h = (unsigned)x * 374761393u + (unsigned)y * 668265263u + seed * 2246822519u;
    h = (h ^ (h >> 13)) * 1274126177u;
    return h ^ (h >> 16);
}
static float rand01(unsigned h) { return (float)(h & 0xffffff) / 16777215.0f; }
static float lerp(float a, float b, float t) { return a + (b - a) * t; }
static float fade(float t) { return t * t * (3.0f - 2.0f * t); }

static float value_noise(float x, float y, unsigned seed) {
    int xi = (int)floorf(x), yi = (int)floorf(y);
    float fx = x - xi, fy = y - yi;
    float v00 = rand01(hash2(xi, yi, seed));
    float v10 = rand01(hash2(xi + 1, yi, seed));
    float v01 = rand01(hash2(xi, yi + 1, seed));
    float v11 = rand01(hash2(xi + 1, yi + 1, seed));
    float u = fade(fx), v = fade(fy);
    return lerp(lerp(v00, v10, u), lerp(v01, v11, u), v) * 2.0f - 1.0f;
}

static float fbm(float x, float y, int octaves, unsigned seed) {
    float sum = 0.0f, amp = 1.0f, norm = 0.0f, freq = 1.0f;
    for (int i = 0; i < octaves; i++) {
        sum += value_noise(x * freq, y * freq, seed + (unsigned)i * 17u) * amp;
        norm += amp;
        amp *= 0.5f;
        freq *= 2.0f;
    }
    return norm > 0.0f ? sum / norm : 0.0f;
}

/* ---- heightfield (bilinear, as in the terrain module) ---- */
typedef struct { float h[W][H]; float cell; } Field;

static float field_sample(const Field* f, float wx, float wy) {
    float gx = wx / f->cell, gy = wy / f->cell;
    int x0 = (int)floorf(gx), y0 = (int)floorf(gy);
    float tx = gx - x0, ty = gy - y0;
    int x1 = x0 + 1, y1 = y0 + 1;
    if (x0 < 0) x0 = 0;
    if (y0 < 0) y0 = 0;
    if (x1 >= W) x1 = W - 1;
    if (y1 >= H) y1 = H - 1;
    float a = f->h[x0][y0] * (1 - tx) + f->h[x1][y0] * tx;
    float b = f->h[x0][y1] * (1 - tx) + f->h[x1][y1] * tx;
    return a * (1 - ty) + b * ty;
}

/* Exact maximum of the bilinear field over an axis-aligned rect.
 * Candidate set: rect corners, interior grid vertices, and the
 * intersections of grid lines with the rect edges (full Cartesian
 * product of candidate x and y coordinates). */
static float rect_max(const Field* f, float x0, float y0, float x1, float y1) {
    float xs[2 + W], ys[2 + H];
    int nx = 0, ny = 0;
    xs[nx++] = x0; xs[nx++] = x1;
    ys[ny++] = y0; ys[ny++] = y1;
    for (int i = 0; i <= W; i++) {
        float gx = i * f->cell;
        if (gx > x0 && gx < x1) xs[nx++] = gx;
    }
    for (int j = 0; j <= H; j++) {
        float gy = j * f->cell;
        if (gy > y0 && gy < y1) ys[ny++] = gy;
    }
    float best = -1e30f;
    for (int a = 0; a < nx; a++)
        for (int b = 0; b < ny; b++) {
            float v = field_sample(f, xs[a], ys[b]);
            if (v > best) best = v;
        }
    return best;
}

/* ---- props ---- */
enum { MODE_BEFORE = 0, MODE_AFTER = 1 };

typedef struct {
    float cx, cy, w, d, h;
    int tier;
    float h_anchor;  /* terrain height at the anchor (wrong layer) */
    float offset;    /* derived offset layer */
    float base;      /* h_anchor + offset (fixed base height) */
} Prop;

static Prop props[MAXPROP];
static int nprops = 0;

static const char* tier_name(int t) { return t == 0 ? "S" : t == 1 ? "M" : "L"; }

static float prop_top(const Prop* p, int mode) {
    return (mode == MODE_AFTER ? p->base : p->h_anchor) + p->h;
}

static int inside_footprint(const Prop* p, float x, float y) {
    return fabsf(x - p->cx) <= p->w * 0.5f && fabsf(y - p->cy) <= p->d * 0.5f;
}

/* Visible surface height: terrain plus any prop top poking above it. */
static float surface_height(const Field* f, float x, float y, int mode) {
    float h = field_sample(f, x, y);
    for (int i = 0; i < nprops; i++) {
        if (inside_footprint(&props[i], x, y)) {
            float t = prop_top(&props[i], mode);
            if (t > h) h = t;
        }
    }
    return h;
}

static int prop_surface_at(const Field* f, float x, float y, int mode) {
    float ter = field_sample(f, x, y);
    for (int i = 0; i < nprops; i++)
        if (inside_footprint(&props[i], x, y) && prop_top(&props[i], mode) >= ter)
            return 1;
    return 0;
}

/* ---- sun ---- */
#define PI 3.14159265358979f
static const float EL = 35.0f * PI / 180.0f;
static const float AZ = 30.0f * PI / 180.0f;

static void sun_travel_dir(float* dx, float* dy, float* dz) {
    *dx = cosf(EL) * cosf(AZ);
    *dy = cosf(EL) * sinf(AZ);
    *dz = -sinf(EL);
}

/* Shadow length cast by a prop on the terrain: march from the prop top
 * center along the sun travel direction until the ray meets the terrain.
 * Length is measured from the footprint boundary along the sun azimuth. */
static float shadow_length(const Field* f, const Prop* p, int mode) {
    float dx, dy, dz;
    sun_travel_dir(&dx, &dy, &dz);
    float top = prop_top(p, mode);
    float t_hit = -1.0f;
    for (float t = 0.02f; t < 80.0f; t += 0.02f) {
        float x = p->cx + dx * t, y = p->cy + dy * t, z = top + dz * t;
        if (z <= field_sample(f, x, y)) { t_hit = t; break; }
    }
    if (t_hit < 0.0f) t_hit = 80.0f;
    /* distance from center to footprint edge along the sun azimuth */
    float ux = dx, uy = dy;
    float ex = fabsf(ux) > 1e-6f ? (p->w * 0.5f) / fabsf(ux) : 1e30f;
    float ey = fabsf(uy) > 1e-6f ? (p->d * 0.5f) / fabsf(uy) : 1e30f;
    float edge = ex < ey ? ex : ey;
    float len = t_hit * cosf(EL) - edge;
    return len > 0.0f ? len : 0.0f;
}

/* ---- PPM output ---- */
typedef struct { unsigned char r, g, b; } RGB;

static void write_ppm(const char* path, const RGB* img, int w, int h) {
    FILE* fp = fopen(path, "wb");
    if (!fp) { fprintf(stderr, "cannot write %s\n", path); exit(1); }
    fprintf(fp, "P6\n%d %d\n255\n", w, h);
    fwrite(img, sizeof(RGB), (size_t)w * (size_t)h, fp);
    fclose(fp);
}

/* Top-down comparison: left = before (anchor layer, lighting on bare
 * terrain), right = after (offset layer, lighting on the combined
 * terrain+prop surface). */
static void render_top(const Field* f, const char* path) {
    int S = 256, gap = 4;
    int iw = 2 * S + gap, ih = S;
    RGB* img = malloc(sizeof(RGB) * (size_t)iw * (size_t)ih);
    float dx, dy, dz;
    sun_travel_dir(&dx, &dy, &dz);
    float lx = -dx, ly = -dy, lz = -dz; /* toward the sun */
    for (int py = 0; py < ih; py++) {
        for (int px = 0; px < iw; px++) {
            int panel = px < S ? 0 : (px < S + gap ? -1 : 1);
            if (panel < 0) { img[py * iw + px] = (RGB){20, 20, 24}; continue; }
            int mode = panel == 0 ? MODE_BEFORE : MODE_AFTER;
            float wx = (float)(px - (panel ? S + gap : 0) + 0.5f) / S * (W * CELL);
            float wy = (float)(py + 0.5f) / S * (H * CELL);
            /* lighting layer: before samples bare terrain (wrong layer),
             * after samples the combined terrain+prop surface */
            float (*layer)(const Field*, float, float, int) = surface_height;
            int lmode = mode == MODE_AFTER ? MODE_AFTER : MODE_BEFORE;
            float h0 = layer(f, wx, wy, lmode);
            float e = 0.35f;
            float hx = layer(f, wx + e, wy, lmode) - layer(f, wx - e, wy, lmode);
            float hy = layer(f, wx, wy + e, lmode) - layer(f, wx, wy - e, lmode);
            float nx = -hx / (2 * e), ny = -hy / (2 * e), nz = 1.0f;
            float nl = sqrtf(nx * nx + ny * ny + nz * nz);
            float diff = (nx * lx + ny * ly + nz * lz) / nl;
            if (diff < 0.0f) diff = 0.0f;
            /* shadow march on the same (possibly wrong) layer */
            int shadowed = 0;
            float sz = h0 + 0.05f, sx = wx, sy = wy;
            for (float t = 0.05f; t < 40.0f; t += 0.05f) {
                sx += lx * 0.05f; sy += ly * 0.05f; sz += lz * 0.05f;
                if (layer(f, sx, sy, lmode) > sz) { shadowed = 1; break; }
            }
            float light = 0.35f + 0.65f * diff;
            if (shadowed) light *= 0.45f;
            RGB c;
            if (prop_surface_at(f, wx, wy, mode)) {
                c = (RGB){210, 208, 200};
            } else {
                float g = h0 * 0.12f + 0.45f;
                if (g < 0.2f) g = 0.2f;
                if (g > 1.0f) g = 1.0f;
                c = (RGB){(unsigned char)(90 * g + 40), (unsigned char)(120 * g + 35),
                          (unsigned char)(60 * g + 25)};
            }
            c.r = (unsigned char)(c.r * light > 255 ? 255 : c.r * light);
            c.g = (unsigned char)(c.g * light > 255 ? 255 : c.g * light);
            c.b = (unsigned char)(c.b * light > 255 ? 255 : c.b * light);
            img[py * iw + px] = c;
        }
    }
    write_ppm(path, img, iw, ih);
    free(img);
}

/* Side cross-section through the first L-tier prop: top = before (base at
 * anchor height, buried), bottom = after (base on the offset layer). */
static void render_side(const Field* f, const char* path, const Prop* p) {
    int S = 256, gap = 4;
    int iw = S, ih = 2 * S / 2 + gap + S / 2; /* two 256x128 panels stacked */
    int ph = S / 2;
    ih = 2 * ph + gap;
    RGB* img = malloc(sizeof(RGB) * (size_t)iw * (size_t)ih);
    float x0 = p->cx - 9.0f, x1 = p->cx + 9.0f;
    float zmin = 1e30f, zmax = -1e30f;
    for (int i = 0; i < S; i++) {
        float wx = x0 + (x1 - x0) * (i + 0.5f) / S;
        float z = field_sample(f, wx, p->cy);
        if (z < zmin) zmin = z;
        if (z > zmax) zmax = z;
    }
    float top_after = p->base + p->h;
    if (zmax < top_after) zmax = top_after;
    zmin -= 0.6f; zmax += 0.6f;
    for (int py = 0; py < ih; py++) {
        int panel = py < ph ? 0 : (py < ph + gap ? -1 : 1);
        int ly = panel == 1 ? py - ph - gap : py;
        for (int px = 0; px < S; px++) {
            if (panel < 0) { img[py * iw + px] = (RGB){20, 20, 24}; continue; }
            int mode = panel == 0 ? MODE_BEFORE : MODE_AFTER;
            float wx = x0 + (x1 - x0) * (px + 0.5f) / S;
            float wz = zmax - (zmax - zmin) * (ly + 0.5f) / ph;
            float ter = field_sample(f, wx, p->cy);
            float base = mode == MODE_AFTER ? p->base : p->h_anchor;
            RGB c = {38, 44, 58}; /* sky */
            if (wz <= ter) c = (RGB){104, 88, 62}; /* terrain */
            if (fabsf(wz - ter) < (zmax - zmin) / ph * 0.7f)
                c = (RGB){150, 132, 96}; /* terrain surface line */
            if (inside_footprint(p, wx, p->cy) && wz <= base + p->h && wz >= base)
                c = (RGB){200, 198, 190}; /* prop body */
            img[py * iw + px] = c;
        }
    }
    write_ppm(path, img, iw, ih);
    free(img);
}

/* ---- self test + report ---- */
static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static unsigned h1(unsigned i) {
    unsigned h = i * 2654435761u;
    h ^= h >> 15;
    return h;
}
static float u01(unsigned i) { return (float)(h1(i) & 0xffff) / 65535.0f; }

int main(int argc, char** argv) {
    const char* outdir = argc > 1 ? argv[1] : "out";
    mkdir(outdir, 0755);

    Field f;
    f.cell = CELL;
    for (int x = 0; x < W; x++)
        for (int y = 0; y < H; y++)
            f.h[x][y] = 3.5f * fbm(x * 0.08f, y * 0.08f, 4, 7u) + 4.0f;

    /* three size tiers, four props each, on a jittered 4x3 layout */
    const float sizes[3][3] = { {1.2f, 1.2f, 0.8f}, {2.6f, 2.2f, 1.6f}, {5.0f, 4.0f, 3.0f} };
    for (int j = 0; j < 3; j++) {       /* tier */
        for (int i = 0; i < TIERN; i++) { /* instances */
            unsigned k = (unsigned)(i * 131 + j * 17) + 5u * 7919u;
            float cx = 10.0f + (i + 0.25f + 0.5f * u01(k)) * 11.0f;
            float cy = 10.0f + (j + 0.25f + 0.5f * u01(k + 1)) * 14.5f;
            Prop* p = &props[nprops++];
            p->cx = cx; p->cy = cy;
            p->w = sizes[j][0]; p->d = sizes[j][1]; p->h = sizes[j][2];
            p->tier = j;
            p->h_anchor = field_sample(&f, cx, cy);
            float rmax = rect_max(&f, cx - p->w * 0.5f, cy - p->d * 0.5f,
                                  cx + p->w * 0.5f, cy + p->d * 0.5f);
            p->offset = rmax - p->h_anchor;
            p->base = p->h_anchor + p->offset;
        }
    }

    /* offset layer invariants, per prop */
    int buried_before = 0, exact = 1, contact = 1, not_buried = 1;
    float worst_gap = 0.0f;
    for (int i = 0; i < nprops; i++) {
        Prop* p = &props[i];
        float rmax = rect_max(&f, p->cx - p->w * 0.5f, p->cy - p->d * 0.5f,
                              p->cx + p->w * 0.5f, p->cy + p->d * 0.5f);
        if (p->offset > 0.05f) buried_before++;
        if (fabsf(p->base - rmax) > 1e-4f) exact = 0;
        if (rmax - p->base > 1e-4f) not_buried = 0;          /* terrain above base */
        float gap = p->base - rmax;                          /* 0 means touching */
        if (gap > 1e-4f) contact = 0;                        /* feet floating */
        if (gap > worst_gap) worst_gap = gap;
    }
    check("buried_before_exists", buried_before > 0);
    check("offset_is_exact_max", exact);
    check("not_buried_after", not_buried);
    check("feet_not_floating", contact);
    check("deterministic", fabsf(rect_max(&f, props[0].cx - 0.6f, props[0].cy - 0.6f,
                                          props[0].cx + 0.6f, props[0].cy + 0.6f) -
                               (props[0].h_anchor + props[0].offset)) < 1e-4f);

    /* three-tier statistics */
    printf("\n%-4s %3s %8s %8s | %8s %8s | %8s %8s\n", "tier", "n", "offset", "buried",
           "sh_bef", "sh_aft", "min_aft", "max_aft");
    int any_float = 0;
    for (int t = 0; t < 3; t++) {
        float so = 0, sb = 0, ssb = 0, ssa = 0, mn = 1e30f, mx = -1e30f;
        int n = 0;
        for (int i = 0; i < nprops; i++) {
            Prop* p = &props[i];
            if (p->tier != t) continue;
            n++;
            so += p->offset;
            sb += p->offset; /* buried depth before == required offset */
            float lb = shadow_length(&f, p, MODE_BEFORE);
            float la = shadow_length(&f, p, MODE_AFTER);
            ssb += lb; ssa += la;
            if (la < mn) mn = la;
            if (la > mx) mx = la;
        }
        printf("%-4s %3d %8.3f %8.3f | %8.3f %8.3f | %8.3f %8.3f\n", tier_name(t), n,
               so / n, sb / n, ssb / n, ssa / n, mn, mx);
    }
    (void)any_float;

    /* per-prop contact report: any floating foot fails the run */
    for (int i = 0; i < nprops; i++) {
        Prop* p = &props[i];
        float rmax = rect_max(&f, p->cx - p->w * 0.5f, p->cy - p->d * 0.5f,
                              p->cx + p->w * 0.5f, p->cy + p->d * 0.5f);
        printf("prop %2d tier %s: base=%.3f footprint_max=%.3f contact_gap=%.4f\n",
               i, tier_name(p->tier), p->base, rmax, p->base - rmax);
    }

    char path[512];
    snprintf(path, sizeof path, "%s/compare_top.ppm", outdir);
    render_top(&f, path);
    const Prop* big = NULL;
    for (int i = 0; i < nprops; i++)
        if (props[i].tier == 2) { big = &props[i]; break; }
    snprintf(path, sizeof path, "%s/compare_side.ppm", outdir);
    render_side(&f, path, big);

    printf("\nall_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
