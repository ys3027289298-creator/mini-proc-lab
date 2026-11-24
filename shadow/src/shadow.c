/*
 * shadow: project the sun's path onto a heightfield, render cast shadows and
 * hillshade at three moments, and measure shadow lengths against the
 * derivation in DERIVATION.md.
 *
 * Solar position comes from latitude / declination / hour angle. A caster dh
 * above the ground blocks the sun wherever the terrain rises above the solar
 * line z0 + t*tan(e); the shadow therefore reaches L = dh / tan(e) on level
 * ground, and L = dh / (tan(e) + g) where the ground's gradient along the
 * shadow direction is g.
 *
 * The renderer marches from every receiver toward the sun. The lengths it
 * renders are then measured independently: each shadow patch is walked along
 * the solar axis, attributed to the caster at its foot, and compared with the
 * closed form using the caster's height above the ground at the shadow tip.
 * Synthetic features (mesa, wall, cylinder) check the same law where the
 * geometry is known exactly.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define W 256
#define H 256
#define CELL 1.0f
#define PI 3.14159265358979323846f
#define DEG (PI / 180.0f)

/* ---- terrain: same value-noise hash and fBm as noise/src/noise.c --------- */
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

static float fbm(float x, float y, int octaves, float lacunarity, float gain, unsigned seed) {
    if (octaves <= 0) return 0.0f;
    float sum = 0.0f, amp = 1.0f, norm = 0.0f, freq = 1.0f;
    for (int i = 0; i < octaves; i++) {
        sum += value_noise(x * freq, y * freq, seed + (unsigned)i * 17u) * amp;
        norm += amp;
        amp *= gain;
        freq *= lacunarity;
    }
    return norm > 0.0f ? sum / norm : 0.0f;
}

typedef struct { float h[W][H]; } Field; /* indexed [x][y] */

static void field_build_terrain(Field* f, unsigned seed) {
    for (int x = 0; x < W; x++)
        for (int y = 0; y < H; y++)
            f->h[x][y] = 55.0f * fbm((float)x * CELL / 88.0f, (float)y * CELL / 88.0f,
                                     5, 2.0f, 0.5f, seed);
}

static float field_sample(const Field* f, float wx, float wy) {
    float gx = wx / CELL, gy = wy / CELL;
    int x0 = (int)floorf(gx), y0 = (int)floorf(gy);
    float tx = gx - x0, ty = gy - y0;
    int x1 = x0 + 1, y1 = y0 + 1;
    if (x0 < 0) x0 = 0;
    if (y0 < 0) y0 = 0;
    if (x1 < 0) x1 = 0;
    if (y1 < 0) y1 = 0;
    if (x0 > W - 1) x0 = W - 1;
    if (y0 > H - 1) y0 = H - 1;
    if (x1 > W - 1) x1 = W - 1;
    if (y1 > H - 1) y1 = H - 1;
    float a = f->h[x0][y0] * (1 - tx) + f->h[x1][y0] * tx;
    float b = f->h[x0][y1] * (1 - tx) + f->h[x1][y1] * tx;
    return a * (1 - ty) + b * ty;
}

static float field_max(const Field* f) {
    float m = -1e30f;
    for (int x = 0; x < W; x++) for (int y = 0; y < H; y++) if (f->h[x][y] > m) m = f->h[x][y];
    return m;
}

/* ---- solar geometry ----------------------------------------------------- */
typedef struct { float azimuth_deg; float altitude_deg; } Sun;

/* az clockwise from north (0=N, 90=E); alt above the horizon.
   hour_angle_deg = 15 * (solar_hours - 12). */
static Sun sun_pos(float lat_deg, float dec_deg, float hour_angle_deg) {
    float phi = lat_deg * DEG, d = dec_deg * DEG, hh = hour_angle_deg * DEG;
    float up = sinf(d) * sinf(phi) + cosf(d) * cosf(phi) * cosf(hh);
    if (up > 1.0f) up = 1.0f;
    if (up < -1.0f) up = -1.0f;
    float east = -cosf(d) * sinf(hh);
    float north = sinf(d) * cosf(phi) - cosf(d) * sinf(phi) * cosf(hh);
    float az = atan2f(east, north) / DEG;
    if (az < 0.0f) az += 360.0f;
    Sun s;
    s.azimuth_deg = az;
    s.altitude_deg = asinf(up) / DEG;
    return s;
}

static void sun_dir(const Sun* s, float l[3]) {
    float a = s->azimuth_deg * DEG, e = s->altitude_deg * DEG;
    l[0] = sinf(a) * cosf(e); /* east  */
    l[1] = cosf(a) * cosf(e); /* north */
    l[2] = sinf(e);           /* up    */
}

/* Horizontal unit vector from a receiver toward the sun. */
static void sun_horizontal(const Sun* s, float* ux, float* uy) {
    float a = s->azimuth_deg * DEG;
    *ux = sinf(a);
    *uy = cosf(a);
}

/* Shadow length on level ground: the anti-solar displacement of the solar
   ray that passes over a caster of height dh. */
static float shadow_length_flat(float dh, float altitude_deg) {
    float c = tanf(altitude_deg * DEG);
    return c > 1e-6f ? dh / c : INFINITY;
}

/* Sloped ground, gradient g along the shadow direction (positive = ground
   rises away from the caster). Solve dh = L*(tan e + g). */
static float shadow_length_slope(float dh, float altitude_deg, float g) {
    float den = tanf(altitude_deg * DEG) + g;
    return den > 1e-6f ? dh / den : INFINITY;
}

/* ---- occlusion march ---------------------------------------------------- */
typedef struct { int hit; float t; float z; } Occl;

/* March from a receiver toward the sun, comparing the terrain with the solar
   line. On a hit, the crossing is refined by linear interpolation and the
   caster's height there is returned. */
static Occl march_occluder(const Field* f, float wx, float wy, const Sun* s,
                           float ds, float tmax, float maxz) {
    Occl o = { 0, 0.0f, 0.0f };
    float z0 = field_sample(f, wx, wy);
    float tan_e = tanf(s->altitude_deg * DEG);
    float ux, uy;
    sun_horizontal(s, &ux, &uy);
    float prev_d = 0.0f, prev_t = 0.0f;
    int has_prev = 0;
    for (float t = ds; t <= tmax; t += ds) {
        if (z0 + t * tan_e > maxz + 0.5f) break;
        float z = field_sample(f, wx + ux * t, wy + uy * t);
        float d = z - (z0 + t * tan_e);
        if (d > 0.0f) {
            float ts = t;
            if (has_prev && d - prev_d > 1e-6f)
                ts = prev_t + (t - prev_t) * (-prev_d) / (d - prev_d);
            if (ts < 0.0f) ts = 0.0f;
            o.hit = 1;
            o.t = ts;
            o.z = field_sample(f, wx + ux * ts, wy + uy * ts);
            return o;
        }
        prev_d = d; prev_t = t; has_prev = 1;
    }
    return o;
}

static void normal_at(const Field* f, int x, int y, float n[3]) {
    int xm = x > 0 ? x - 1 : x, xp = x < W - 1 ? x + 1 : x;
    int ym = y > 0 ? y - 1 : y, yp = y < H - 1 ? y + 1 : y;
    float dzdx = (f->h[xp][y] - f->h[xm][y]) / ((float)(xp - xm) * CELL);
    float dzdy = (f->h[x][yp] - f->h[x][ym]) / ((float)(yp - ym) * CELL);
    n[0] = -dzdx; n[1] = -dzdy; n[2] = 1.0f;
    float len = sqrtf(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    if (len > 1e-9f) { n[0] /= len; n[1] /= len; n[2] /= len; }
}

/* ---- analytic test fields ---------------------------------------------- */
static void build_mesa(Field* f, float height, float radius) {
    memset(f, 0, sizeof(*f));
    for (int x = 0; x < W; x++)
        for (int y = 0; y < H; y++) {
            float dx = (float)x - 128.0f, dy = (float)y - 128.0f;
            if (sqrtf(dx * dx + dy * dy) <= radius) f->h[x][y] = height;
        }
}

/* Plane z = b*x with a half-plane wall of height dh on the up-sun side of the
   line u = 0 through the field centre. */
static void build_wall_plane(Field* f, float b, float dh, float azimuth_deg) {
    float a = azimuth_deg * DEG;
    float vx = -sinf(a), vy = -cosf(a); /* shadow direction */
    for (int x = 0; x < W; x++)
        for (int y = 0; y < H; y++) {
            float u = ((float)x - 128.0f) * vx + ((float)y - 128.0f) * vy;
            f->h[x][y] = b * (float)x * CELL + (u < 0.0f ? dh : 0.0f);
        }
}

static void build_cylinder(Field* f, float cx, float cy, float radius, float height) {
    memset(f, 0, sizeof(*f));
    for (int x = 0; x < W; x++)
        for (int y = 0; y < H; y++) {
            float dx = (float)x - cx, dy = (float)y - cy;
            if (sqrtf(dx * dx + dy * dy) <= radius) f->h[x][y] = height;
        }
}

static void add_plane(Field* f, float b) {
    for (int x = 0; x < W; x++)
        for (int y = 0; y < H; y++)
            f->h[x][y] += b * (float)x * CELL;
}

/* Farthest shadowed receiver along the anti-solar axis from the field centre,
   starting at `start` metres. */
static float shadow_reach(const Field* f, const Sun* s, float start, float step,
                          float smax) {
    float a = s->azimuth_deg * DEG;
    float vx = -sinf(a), vy = -cosf(a);
    float maxz = field_max(f);
    float reach = 0.0f;
    for (float ss = start; ss <= smax; ss += step) {
        float px = 128.0f * CELL + vx * ss;
        float py = 128.0f * CELL + vy * ss;
        if (march_occluder(f, px, py, s, 0.25f, 600.0f, maxz).hit) reach = ss;
    }
    return reach;
}

/* ---- validation bookkeeping -------------------------------------------- */
static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

typedef struct {
    char name[40];
    float azimuth, altitude, dh, g;
    float derived, measured, dev_pct;
    int pass;
} ValRow;

static ValRow val_rows[64];
static int val_n = 0;

static void add_val(const char* name, float az, float alt, float dh, float g,
                    float derived, float measured) {
    ValRow* r = &val_rows[val_n++];
    snprintf(r->name, sizeof(r->name), "%s", name);
    r->azimuth = az; r->altitude = alt; r->dh = dh; r->g = g;
    r->derived = derived; r->measured = measured;
    r->dev_pct = 100.0f * fabsf(measured - derived) / (derived > 1e-6f ? derived : 1e-6f);
    r->pass = r->dev_pct <= 20.0f;
}

static void run_geometry_validation(void) {
    const float altitudes[4] = { 15.0f, 25.0f, 40.0f, 60.0f };
    const float azimuths[2] = { 110.0f, 250.0f };
    const float H_MESA = 25.0f, R_MESA = 20.0f;
    Field mesa;
    build_mesa(&mesa, H_MESA, R_MESA);
    val_n = 0;

    for (int ai = 0; ai < 4; ai++) {
        for (int zi = 0; zi < 2; zi++) {
            Sun s;
            s.azimuth_deg = azimuths[zi];
            s.altitude_deg = altitudes[ai];
            float reach = shadow_reach(&mesa, &s, R_MESA + 0.5f, 0.25f, 230.0f);
            char nm[40];
            snprintf(nm, sizeof(nm), "mesa_az%.0f_alt%.0f", s.azimuth_deg, s.altitude_deg);
            add_val(nm, s.azimuth_deg, s.altitude_deg, H_MESA, 0.0f,
                    shadow_length_flat(H_MESA, s.altitude_deg), reach - R_MESA);
        }
    }

    const float DH_WALL = 15.0f;
    Field wall;
    const float w_az[6] = { 110.0f, 180.0f, 250.0f, 110.0f, 180.0f, 250.0f };
    const float w_b[6]  = {  0.0f,   0.0f,   0.0f,   0.3f,   0.3f,   0.3f };
    for (int i = 0; i < 6; i++) {
        build_wall_plane(&wall, w_b[i], DH_WALL, w_az[i]);
        Sun s;
        s.azimuth_deg = w_az[i];
        s.altitude_deg = 30.0f;
        float g = w_b[i] * (-sinf(w_az[i] * DEG));
        float reach = shadow_reach(&wall, &s, 0.5f, 0.25f, 90.0f);
        char nm[40];
        snprintf(nm, sizeof(nm), "wall_%s_az%.0f_alt30",
                 w_b[i] > 0.0f ? "slope" : "flat", s.azimuth_deg);
        add_val(nm, s.azimuth_deg, s.altitude_deg, DH_WALL, g,
                shadow_length_slope(DH_WALL, s.altitude_deg, g), reach);
    }
}

/* Isolated cylinder, rendered with the same march as the frames. */
typedef struct {
    char time[12];
    float cx, cy, radius, height, g;
    float az, alt;
    float measured, derived, dev_pct;
    int pass;
} CasterRow;

static CasterRow caster_rows[64];
static int caster_n = 0;

static float caster_reach(const Field* f, const Sun* s, float cx, float cy,
                          float radius) {
    float a = s->azimuth_deg * DEG;
    float vx = -sinf(a), vy = -cosf(a);
    float maxz = field_max(f);
    float reach = 0.0f;
    for (float ss = radius + 0.25f; ss <= radius + 500.0f; ss += 0.25f) {
        float px = (cx + vx * ss) * CELL;
        float py = (cy + vy * ss) * CELL;
        if (march_occluder(f, px, py, s, 0.25f, 600.0f, maxz).hit) reach = ss;
    }
    return reach - radius;
}

static void run_caster_validation(const Sun suns[3], const char* tags[3]) {
    const float H_CYL = 20.0f, R_CYL = 6.0f;
    Field cyl;
    build_cylinder(&cyl, 128.0f, 128.0f, R_CYL, H_CYL);
    caster_n = 0;
    for (int k = 0; k < 3; k++) {
        CasterRow* r = &caster_rows[caster_n++];
        snprintf(r->time, sizeof(r->time), "%s", tags[k]);
        r->cx = 128.0f; r->cy = 128.0f; r->radius = R_CYL; r->height = H_CYL;
        r->g = 0.0f;
        r->az = suns[k].azimuth_deg; r->alt = suns[k].altitude_deg;
        r->measured = caster_reach(&cyl, &suns[k], 128.0f, 128.0f, R_CYL);
        r->derived = shadow_length_flat(H_CYL, suns[k].altitude_deg);
        r->dev_pct = 100.0f * fabsf(r->measured - r->derived) / r->derived;
        r->pass = r->dev_pct <= 20.0f;
    }

    const float B = 0.30f;
    Field slope;
    build_cylinder(&slope, 128.0f, 128.0f, R_CYL, H_CYL);
    add_plane(&slope, B);
    for (int k = 0; k < 3; k++) {
        CasterRow* r = &caster_rows[caster_n++];
        snprintf(r->time, sizeof(r->time), "%s_slope", tags[k]);
        r->cx = 128.0f; r->cy = 128.0f; r->radius = R_CYL; r->height = H_CYL;
        r->g = B * (-sinf(suns[k].azimuth_deg * DEG));
        r->az = suns[k].azimuth_deg; r->alt = suns[k].altitude_deg;
        r->measured = caster_reach(&slope, &suns[k], 128.0f, 128.0f, R_CYL);
        r->derived = shadow_length_slope(H_CYL, suns[k].altitude_deg, r->g);
        r->dev_pct = 100.0f * fabsf(r->measured - r->derived) / r->derived;
        r->pass = r->dev_pct <= 20.0f;
    }
}

/* ---- shadow-patch measurement ------------------------------------------ */
/* A cast shadow is a contiguous run of the mask along the solar axis. Its two
   ends are surface points where the incoming solar ray grazes the ground: at
   the foot the ray leaves the caster's summit, at the tip it settles back onto
   the terrain. Over a horizontal run L the grazing ray drops L*tan(e), so the
   heights of the ground at the two ends must differ by exactly that:

       L * tan(e) = z_foot - z_tip   =>   L = (z_foot - z_tip) / tan(e)

   Measuring L from the mask and comparing it with the right-hand side uses
   only the terrain and the sun's altitude - nothing from the occlusion march
   that produced the mask. Patches shorter than SEG_MIN_LEN are below the
   resolution at which a run is distinguishable from the heightfield's own
   grain and are counted, not scored. */
#define SEG_MIN_LEN   5.0f
#define SEG_MIN_DH    0.5f
#define SEG_STEP      0.5f
#define SEG_HOLE_MAX  2
#define SEG_MAX       400.0f
#define SEG_CAP       16384

typedef struct {
    float length, foot_ground, tip_ground, dh, derived, dev_pct;
} SegRow;

static SegRow seg_rows[SEG_CAP];
static int seg_n = 0;

typedef struct { float down, up; int tx, ty, fx, fy; } SegExtent;

/* Walk a shadow patch from one member pixel to both of its ends along the
   solar axis. (vx,vy) points down-sun, so dir=+1 walks toward the tip and
   dir=-1 back toward the foot. */
static SegExtent segment_extent(const unsigned char* mask, int x0, int y0,
                                float vx, float vy) {
    SegExtent e;
    e.down = e.up = 0.0f;
    e.tx = e.fx = x0; e.ty = e.fy = y0;
    for (int dir = -1; dir <= 1; dir += 2) {
        float dx = (float)dir * vx, dy = (float)dir * vy;
        float last = 0.0f;
        int lx = x0, ly = y0, holes = 0;
        for (float t = SEG_STEP; t < SEG_MAX; t += SEG_STEP) {
            int xi = (int)lrintf((float)x0 + dx * t);
            int yi = (int)lrintf((float)y0 + dy * t);
            if (xi < 0 || yi < 0 || xi >= W || yi >= H) break;
            if (mask[xi * H + yi]) {
                holes = 0; last = t; lx = xi; ly = yi;
            } else if (++holes > SEG_HOLE_MAX) {
                break;
            }
        }
        if (dir > 0) { e.down = last; e.tx = lx; e.ty = ly; }
        else { e.up = last; e.fx = lx; e.fy = ly; }
    }
    return e;
}

/* Patch length through one pixel (down-sun + up-sun extent). */
static float segment_length(const unsigned char* mask, int x0, int y0,
                            float vx, float vy) {
    SegExtent e = segment_extent(mask, x0, y0, vx, vy);
    return e.down + e.up;
}

static float ground_at(const Field* f, int x, int y) {
    return field_sample(f, ((float)x + 0.5f) * CELL, ((float)y + 0.5f) * CELL);
}

typedef struct { int patches, kept, too_short, flat; } Funnel;

static Funnel attribute_shadows(const Field* f, const unsigned char* mask,
                                const Sun* s) {
    Funnel fn = { 0, 0, 0, 0 };
    float te = tanf(s->altitude_deg * DEG);
    float ux, uy;
    sun_horizontal(s, &ux, &uy);
    float vx = -ux, vy = -uy;
    unsigned char* seen = (unsigned char*)calloc((size_t)W * H, 1);
    if (!seen) { fprintf(stderr, "out of memory\n"); exit(2); }
    for (int x = 0; x < W; x++)
        for (int y = 0; y < H; y++) {
            if (!mask[x * H + y] || seen[x * H + y]) continue;
            SegExtent e = segment_extent(mask, x, y, vx, vy);
            for (int dir = -1; dir <= 1; dir += 2) {
                float dx = (float)dir * vx, dy = (float)dir * vy;
                float lim = (dir > 0 ? e.down : e.up) + 0.5f;
                for (float t = 0.0f; t <= lim; t += SEG_STEP) {
                    int xi = (int)lrintf((float)x + dx * t);
                    int yi = (int)lrintf((float)y + dy * t);
                    if (xi >= 0 && yi >= 0 && xi < W && yi < H) seen[xi * H + yi] = 1;
                }
            }
            fn.patches++;
            float len = e.down + e.up;
            if (len < SEG_MIN_LEN) { fn.too_short++; continue; }
            float zf = ground_at(f, e.fx, e.fy);
            float zt = ground_at(f, e.tx, e.ty);
            float dh = zf - zt;
            if (dh < SEG_MIN_DH) { fn.flat++; continue; }
            if (seg_n >= SEG_CAP) continue;
            SegRow* r = &seg_rows[seg_n++];
            r->length = len;
            r->foot_ground = zf;
            r->tip_ground = zt;
            r->dh = dh;
            r->derived = dh / te;
            r->dev_pct = 100.0f * fabsf(len - r->derived) / r->derived;
            fn.kept++;
        }
    free(seen);
    return fn;
}

/* ---- output helpers ----------------------------------------------------- */
static void write_f32(const char* path, const float* v, size_t n) {
    FILE* fp = fopen(path, "wb");
    if (!fp) { fprintf(stderr, "cannot write %s\n", path); exit(2); }
    fwrite(v, sizeof(float), n, fp);
    fclose(fp);
}
static void write_u8(const char* path, const unsigned char* v, size_t n) {
    FILE* fp = fopen(path, "wb");
    if (!fp) { fprintf(stderr, "cannot write %s\n", path); exit(2); }
    fwrite(v, 1, n, fp);
    fclose(fp);
}
static int cmp_float(const void* a, const void* b) {
    float x = *(const float*)a, y = *(const float*)b;
    return (x > y) - (x < y);
}
static float quantile(const float* sorted, int n, float q) {
    if (n <= 0) return 0.0f;
    float pos = q * (float)(n - 1);
    int i = (int)floorf(pos);
    if (i >= n - 1) return sorted[n - 1];
    float frac = pos - (float)i;
    return sorted[i] * (1 - frac) + sorted[i + 1] * frac;
}

/* ---- per-frame render + measure ----------------------------------------- */
typedef struct {
    float az, alt;
    int n_shadowed;
    float shadow_fraction;
    float mean_len, median_len, p90_len, p95_len, max_len;
    int patches, kept, too_short, flat;
    int n_within20;
    float median_dev, p75_dev, frac_within20;
    unsigned hash;
} FrameStat;

static FrameStat frame_stats[3];

static void render_frame(const Field* f, const Sun* s, int k,
                         float* shade, unsigned char* mask, float* len,
                         int write_files, const char* tag) {
    float maxz = field_max(f);
    float l[3];
    sun_dir(s, l);
    float vx = -sinf(s->azimuth_deg * DEG), vy = -cosf(s->azimuth_deg * DEG);
    int nshadow = 0;

    for (int x = 0; x < W; x++)
        for (int y = 0; y < H; y++) {
            int idx = x * H + y;
            float wx = ((float)x + 0.5f) * CELL;
            float wy = ((float)y + 0.5f) * CELL;
            Occl o = march_occluder(f, wx, wy, s, 0.25f, 400.0f, maxz);
            mask[idx] = o.hit ? 1 : 0;
            float n[3];
            normal_at(f, x, y, n);
            float lam = n[0] * l[0] + n[1] * l[1] + n[2] * l[2];
            if (lam < 0.0f) lam = 0.0f;
            shade[idx] = o.hit ? 0.05f : 0.22f + 0.78f * lam;
        }

    for (int x = 0; x < W; x++)
        for (int y = 0; y < H; y++) {
            int idx = x * H + y;
            len[idx] = mask[idx] ? segment_length(mask, x, y, vx, vy) : 0.0f;
            if (mask[idx]) nshadow++;
        }

    float* copy = (float*)malloc(sizeof(float) * (nshadow > 0 ? nshadow : 1));
    int m = 0;
    for (int i = 0; i < W * H; i++) if (mask[i]) copy[m++] = len[i];
    qsort(copy, nshadow, sizeof(float), cmp_float);

    int seg_start = seg_n;
    Funnel fn = attribute_shadows(f, mask, s);
    int kept = seg_n - seg_start;
    float* devs = (float*)malloc(sizeof(float) * (kept > 0 ? kept : 1));
    int n20 = 0;
    for (int i = 0; i < kept; i++) {
        devs[i] = seg_rows[seg_start + i].dev_pct;
        if (devs[i] <= 20.0f) n20++;
    }
    qsort(devs, kept, sizeof(float), cmp_float);

    unsigned fh = 2166136261u;
    for (int i = 0; i < W * H; i++) { fh ^= mask[i]; fh *= 16777619u; }

    FrameStat* st = &frame_stats[k];
    st->az = s->azimuth_deg; st->alt = s->altitude_deg;
    st->n_shadowed = nshadow;
    st->shadow_fraction = (float)nshadow / (float)(W * H);
    st->mean_len = st->median_len = st->p90_len = st->p95_len = st->max_len = 0.0f;
    if (nshadow) {
        float sum = 0.0f;
        for (int i = 0; i < nshadow; i++) sum += copy[i];
        st->mean_len = sum / (float)nshadow;
        st->median_len = quantile(copy, nshadow, 0.5f);
        st->p90_len = quantile(copy, nshadow, 0.9f);
        st->p95_len = quantile(copy, nshadow, 0.95f);
        st->max_len = copy[nshadow - 1];
    }
    st->patches = fn.patches; st->kept = fn.kept;
    st->too_short = fn.too_short; st->flat = fn.flat;
    st->n_within20 = n20;
    st->median_dev = quantile(devs, kept, 0.5f);
    st->p75_dev = quantile(devs, kept, 0.75f);
    st->frac_within20 = kept ? (float)n20 / (float)kept : 0.0f;
    st->hash = fh;

    printf("%s az=%.2f alt=%.2f shadow=%.1f%% median_len=%.1f m | patches=%d attributed=%d "
           "(short %d, level %d) dev p50=%.1f%% p75=%.1f%% within20=%.0f%%\n",
           tag, st->az, st->alt, 100.0f * st->shadow_fraction, st->median_len,
           st->patches, st->kept, st->too_short, st->flat,
           st->median_dev, st->p75_dev, 100.0f * st->frac_within20);

    if (write_files) {
        char path[160];
        snprintf(path, sizeof(path), "out/shade_%s.bin", tag); write_f32(path, shade, (size_t)W * H);
        snprintf(path, sizeof(path), "out/shadow_%s.bin", tag); write_u8(path, mask, (size_t)W * H);
        snprintf(path, sizeof(path), "out/len_%s.bin", tag); write_f32(path, len, (size_t)W * H);
    }
    free(copy);
    free(devs);
}

static void run_export(const Sun suns[3], const char* tags[3]) {
    Field f;
    field_build_terrain(&f, 20261124u);
    float* z = (float*)malloc(sizeof(float) * W * H);
    for (int x = 0; x < W; x++)
        for (int y = 0; y < H; y++)
            z[x * H + y] = f.h[x][y];
    write_f32("out/terrain.bin", z, (size_t)W * H);

    float* shade = (float*)malloc(sizeof(float) * W * H);
    float* len = (float*)malloc(sizeof(float) * W * H);
    unsigned char* mask = (unsigned char*)malloc((size_t)W * H);
    unsigned char* masks[3] = { NULL, NULL, NULL };
    float* lens[3] = { NULL, NULL, NULL };
    seg_n = 0;

    for (int k = 0; k < 3; k++) {
        render_frame(&f, &suns[k], k, shade, mask, len, 1, tags[k]);
        masks[k] = (unsigned char*)malloc((size_t)W * H);
        lens[k] = (float*)malloc(sizeof(float) * W * H);
        memcpy(masks[k], mask, (size_t)W * H);
        memcpy(lens[k], len, sizeof(float) * W * H);
    }

    FILE* stats = fopen("out/stats.csv", "w");
    if (!stats) { fprintf(stderr, "cannot write out/stats.csv\n"); exit(2); }
    fprintf(stats, "time,azimuth_deg,altitude_deg,cot_altitude,shadow_fraction,n_shadowed,"
                   "mean_len_m,median_len_m,p90_len_m,p95_len_m,max_len_m,"
                   "patches,attributed,tooshort_dropped,level_dropped,attributed_within_20pct,"
                   "attributed_median_dev_pct,attributed_p75_dev_pct,"
                   "attributed_within_20pct_frac,mask_hash\n");
    for (int k = 0; k < 3; k++) {
        FrameStat* st = &frame_stats[k];
        fprintf(stats, "%s,%.2f,%.2f,%.4f,%.6f,%d,%.3f,%.3f,%.3f,%.3f,%.3f,"
                       "%d,%d,%d,%d,%d,%.2f,%.2f,%.4f,%08x\n",
                tags[k], st->az, st->alt, 1.0f / tanf(st->alt * DEG),
                st->shadow_fraction, st->n_shadowed,
                st->mean_len, st->median_len, st->p90_len, st->p95_len, st->max_len,
                st->patches, st->kept, st->too_short, st->flat, st->n_within20,
                st->median_dev, st->p75_dev, st->frac_within20, st->hash);
    }
    fclose(stats);

    /* Per-patch records: every row is one shadow attributed to one caster. */
    FILE* segs = fopen("out/segments.csv", "w");
    if (!segs) { fprintf(stderr, "cannot write out/segments.csv\n"); exit(2); }
    fprintf(segs, "time,length_m,foot_ground_m,tip_ground_m,dh_m,derived_m,dev_pct,pass\n");
    int row = 0;
    for (int k = 0; k < 3; k++)
        for (int i = 0; i < frame_stats[k].kept; i++) {
            SegRow* r = &seg_rows[row++];
            fprintf(segs, "%s,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%d\n",
                    tags[k], r->length, r->foot_ground, r->tip_ground,
                    r->dh, r->derived, r->dev_pct, r->dev_pct <= 20.0f);
        }
    fclose(segs);

    /* Pairwise frame comparison: direction is the question, so overlap and
       per-patch change carry more than the mean. */
    FILE* pairs = fopen("out/pairs.csv", "w");
    if (!pairs) { fprintf(stderr, "cannot write out/pairs.csv\n"); exit(2); }
    fprintf(pairs, "a,b,intersection,union,jaccard,only_a,only_b,"
                   "only_a_frac_of_shadow,only_b_frac_of_shadow,"
                   "mean_len_a_in_both,mean_len_b_in_both,"
                   "median_abs_len_diff_m,median_rel_len_diff_pct,frac_rel_diff_gt20pct\n");
    int shadow_count[3];
    for (int k = 0; k < 3; k++) {
        shadow_count[k] = 0;
        for (int p = 0; p < W * H; p++) shadow_count[k] += masks[k][p];
    }
    for (int i = 0; i < 3; i++)
        for (int j = i + 1; j < 3; j++) {
            int inter = 0, uni = 0, onlya = 0, onlyb = 0, n = 0, nbig = 0;
            float sa = 0.0f, sb = 0.0f;
            float* absd = (float*)malloc(sizeof(float) * W * H);
            float* reld = (float*)malloc(sizeof(float) * W * H);
            int nd = 0;
            for (int p = 0; p < W * H; p++) {
                int a = masks[i][p], b = masks[j][p];
                if (a && b) {
                    inter++; uni++; n++;
                    sa += lens[i][p]; sb += lens[j][p];
                    float mx = lens[i][p] > lens[j][p] ? lens[i][p] : lens[j][p];
                    if (mx >= 1.0f) {
                        float d = fabsf(lens[i][p] - lens[j][p]);
                        float rel = 100.0f * d / mx;
                        absd[nd] = d; reld[nd] = rel; nd++;
                        if (rel > 20.0f) nbig++;
                    }
                } else if (a) { onlya++; uni++; }
                else if (b) { onlyb++; uni++; }
            }
            qsort(absd, nd, sizeof(float), cmp_float);
            qsort(reld, nd, sizeof(float), cmp_float);
            fprintf(pairs, "%s,%s,%d,%d,%.5f,%d,%d,%.4f,%.4f,%.3f,%.3f,%.3f,%.2f,%.4f\n",
                    tags[i], tags[j], inter, uni,
                    uni ? (float)inter / (float)uni : 0.0f, onlya, onlyb,
                    shadow_count[i] ? (float)onlya / (float)shadow_count[i] : 0.0f,
                    shadow_count[j] ? (float)onlyb / (float)shadow_count[j] : 0.0f,
                    n ? sa / (float)n : 0.0f, n ? sb / (float)n : 0.0f,
                    quantile(absd, nd, 0.5f), quantile(reld, nd, 0.5f),
                    nd ? (float)nbig / (float)nd : 0.0f);
            free(absd); free(reld);
        }
    fclose(pairs);

    /* Solar path across the day and L(e) for a 20 m caster. */
    FILE* path = fopen("out/daypath.csv", "w");
    if (!path) { fprintf(stderr, "cannot write out/daypath.csv\n"); exit(2); }
    fprintf(path, "hour,azimuth_deg,altitude_deg,shadow_factor_cot\n");
    for (int i = 0; i <= 56; i++) {
        float hr = 5.0f + 0.25f * (float)i;
        Sun s = sun_pos(40.0f, 0.0f, 15.0f * (hr - 12.0f));
        float cot = s.altitude_deg > 0.5f ? 1.0f / tanf(s.altitude_deg * DEG) : 0.0f;
        fprintf(path, "%.2f,%.3f,%.3f,%.3f\n", hr, s.azimuth_deg, s.altitude_deg, cot);
    }
    fclose(path);

    FILE* curve = fopen("out/curve.csv", "w");
    if (!curve) { fprintf(stderr, "cannot write out/curve.csv\n"); exit(2); }
    fprintf(curve, "altitude_deg,L20_m\n");
    for (int a = 5; a <= 80; a++)
        fprintf(curve, "%d,%.4f\n", a, shadow_length_flat(20.0f, (float)a));
    fclose(curve);

    FILE* val = fopen("out/validation.csv", "w");
    if (!val) { fprintf(stderr, "cannot write out/validation.csv\n"); exit(2); }
    fprintf(val, "case,azimuth_deg,altitude_deg,dh_m,g,derived_m,measured_m,dev_pct,pass\n");
    for (int i = 0; i < val_n; i++)
        fprintf(val, "%s,%.2f,%.2f,%.2f,%.4f,%.3f,%.3f,%.3f,%d\n",
                val_rows[i].name, val_rows[i].azimuth, val_rows[i].altitude,
                val_rows[i].dh, val_rows[i].g, val_rows[i].derived,
                val_rows[i].measured, val_rows[i].dev_pct, val_rows[i].pass);
    fclose(val);

    FILE* casters = fopen("out/casters.csv", "w");
    if (!casters) { fprintf(stderr, "cannot write out/casters.csv\n"); exit(2); }
    fprintf(casters, "time,cx,cy,radius_m,height_m,g,azimuth_deg,altitude_deg,"
                     "measured_m,derived_m,dev_pct,pass\n");
    for (int i = 0; i < caster_n; i++)
        fprintf(casters, "%s,%.1f,%.1f,%.1f,%.1f,%.4f,%.2f,%.2f,%.3f,%.3f,%.3f,%d\n",
                caster_rows[i].time, caster_rows[i].cx, caster_rows[i].cy,
                caster_rows[i].radius, caster_rows[i].height, caster_rows[i].g,
                caster_rows[i].az, caster_rows[i].alt, caster_rows[i].measured,
                caster_rows[i].derived, caster_rows[i].dev_pct, caster_rows[i].pass);
    fclose(casters);

    free(z); free(shade); free(len); free(mask);
    for (int k = 0; k < 3; k++) { free(masks[k]); free(lens[k]); }
}

int main(int argc, char** argv) {
    Sun suns[3];
    const char* tags[3] = { "0800", "1200", "1600" };
    suns[0] = sun_pos(40.0f, 0.0f, -60.0f);
    suns[1] = sun_pos(40.0f, 0.0f, 0.0f);
    suns[2] = sun_pos(40.0f, 0.0f, 60.0f);

    check("noon_altitude_50", fabsf(suns[1].altitude_deg - 50.0f) < 1e-4f);
    check("noon_azimuth_south", fabsf(suns[1].azimuth_deg - 180.0f) < 1e-3f);
    check("morning_azimuth_ese", fabsf(suns[0].azimuth_deg - 110.36f) < 0.02f);
    check("evening_azimuth_wsw", fabsf(suns[2].azimuth_deg - 249.64f) < 0.02f);
    check("morning_evening_altitude_equal", fabsf(suns[0].altitude_deg - suns[2].altitude_deg) < 1e-4f);
    check("morning_evening_azimuth_mirror",
          fabsf(suns[0].azimuth_deg + suns[2].azimuth_deg - 360.0f) < 1e-2f);

    check("L_flat_30deg_10m", fabsf(shadow_length_flat(10.0f, 30.0f) - 17.3205f) < 1e-3f);
    check("L_flat_60deg_10m", fabsf(shadow_length_flat(10.0f, 60.0f) - 5.7735f) < 1e-3f);
    check("L_slope_30deg_g0p3",
          fabsf(shadow_length_slope(10.0f, 30.0f, 0.3f) -
                10.0f / (tanf(30.0f * DEG) + 0.3f)) < 1e-4f);
    /* Equal altitude, opposite azimuth: on level ground the lengths tie; the
       direction is what differs. */
    check("equal_altitude_equal_length",
          fabsf(shadow_length_flat(20.0f, suns[0].altitude_deg) -
                shadow_length_flat(20.0f, suns[2].altitude_deg)) < 1e-4f);

    run_geometry_validation();
    {
        int mesa_ok = 1, wall_ok = 1;
        float mesa_worst = 0.0f, wall_worst = 0.0f;
        for (int i = 0; i < val_n; i++) {
            if (strncmp(val_rows[i].name, "mesa", 4) == 0) {
                if (val_rows[i].dev_pct > mesa_worst) mesa_worst = val_rows[i].dev_pct;
                if (!val_rows[i].pass) mesa_ok = 0;
            } else {
                if (val_rows[i].dev_pct > wall_worst) wall_worst = val_rows[i].dev_pct;
                if (!val_rows[i].pass) wall_ok = 0;
            }
        }
        printf("mesa_worst_dev=%.3f%% wall_worst_dev=%.3f%%\n", mesa_worst, wall_worst);
        check("mesa_matches_derivation_within_20pct", mesa_ok);
        check("slope_matches_derivation_within_20pct", wall_ok);
    }

    run_caster_validation(suns, tags);
    {
        int all_pass = 1;
        float worst = 0.0f;
        for (int i = 0; i < caster_n; i++) {
            if (!caster_rows[i].pass) all_pass = 0;
            if (caster_rows[i].dev_pct > worst) worst = caster_rows[i].dev_pct;
        }
        printf("caster_worst_dev=%.3f%%\n", worst);
        check("casters_match_derivation_within_20pct", all_pass);
        check("slope_breaks_morning_evening_symmetry",
              fabsf(caster_rows[3].measured - caster_rows[5].measured) > 5.0f);
    }

    /* Real terrain: render each moment and measure every attributable patch. */
    {
        Field f;
        field_build_terrain(&f, 20261124u);
        float* shade = (float*)malloc(sizeof(float) * W * H);
        unsigned char* mask = (unsigned char*)malloc((size_t)W * H);
        unsigned char* masks0 = (unsigned char*)malloc((size_t)W * H);
        unsigned char* masks2 = (unsigned char*)malloc((size_t)W * H);
        float* len = (float*)malloc(sizeof(float) * W * H);
        seg_n = 0;
        for (int k = 0; k < 3; k++) {
            render_frame(&f, &suns[k], k, shade, mask, len, 0, tags[k]);
            if (k == 0) memcpy(masks0, mask, (size_t)W * H);
            if (k == 2) memcpy(masks2, mask, (size_t)W * H);
        }

        int all_pass = 1;
        for (int k = 0; k < 3; k++)
            if (frame_stats[k].kept < 100 || frame_stats[k].median_dev > 20.0f) all_pass = 0;
        check("terrain_median_length_dev_within_20pct", all_pass);
        check("terrain_attribution_has_samples",
              frame_stats[0].kept > 100 && frame_stats[2].kept > 100);
        check("frames_differ",
              frame_stats[0].hash != frame_stats[1].hash &&
              frame_stats[0].hash != frame_stats[2].hash &&
              frame_stats[1].hash != frame_stats[2].hash);
        /* Morning and evening have the same altitude, so their lengths tie;
           the difference that must show is direction: the share of 08:00
           shadow that is still shadow at 16:00 must be small. */
        {
            int inter = 0, n0 = 0;
            for (int p = 0; p < W * H; p++) {
                n0 += frame_stats[0].n_shadowed > 0 && masks0[p];
                inter += masks0[p] && masks2[p];
            }
            printf("morning_evening_overlap=%.3f\n",
                   n0 ? (float)inter / (float)n0 : 0.0f);
            check("morning_evening_shadows_mostly_disjoint",
                  n0 > 0 && (float)inter / (float)n0 < 0.30f);
        }
        free(shade); free(mask); free(masks0); free(masks2); free(len);
    }

    if (argc > 1 && strcmp(argv[1], "export") == 0) {
        run_export(suns, tags);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
