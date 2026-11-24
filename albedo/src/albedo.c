
/* albedo: two-layer albedo texture with decoupled per-channel faders.
 *
 *   L(x,y,c) = log A(x,y,c)                      log reflectance (optical density)
 *   S(x,y,c) = gauss(L, SIGMA_SPLIT)             structure layer (low-pass)
 *   m(x,y,c) = L - S                             micro layer (high-pass, zero mean)
 *
 *   L' = Sbar_c + ks_c * (S_c - Sbar_c) + km_c * m_c
 *   A  = clamp( exp(L'), ALB_LO, ALB_HI )
 *
 * The two layers are the two halves of one filter split, so in the log domain
 * they are additive and orthogonal by construction: a gain on S cannot move m,
 * and a gain on m cannot move S.  The log domain also makes m *relative* - a
 * +0.05 gain on m is +5 % reflectance everywhere, whatever the material
 * brightness - which is what "micro noise" means.
 *
 * Truncation is part of each layer, not an accident of an 8-bit buffer:
 *   structure  clamped to [log 0.030, log 0.920] reflectance
 *   micro      clamped to [-M_CLIP, +M_CLIP] log units = [-12.2 %, +13.9 %]
 *
 *   ./albedo test    self-checks: orthogonality, clipping, weather sanity
 *   ./albedo build   writes out/ CSV reports and out/ PPM images (512x512)
 *
 * Weather: directional sun + sky ambient, ray-marched hard shadows, per-pixel
 * albedo*irradiance, luminance-preserving filmic tonemap, one fixed white
 * balance calibrated once on a neutral card under neutral light.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define N 512
#define SIGMA_SPLIT 2.75f   /* px; largest sigma keeping >=80 % of the 28 px
                               fine texture in the structure layer (out/scale_decision.csv) */
#define M_SIGMA 0.045f      /* micro layer design sigma, log units (4.5 %)     */
#define M_CLIP 0.130f       /* micro truncation, log units                     */
#define S_REFL_LO 0.030f    /* structure truncation, reflectance               */
#define S_REFL_HI 0.920f
#define LS_LO (-3.5065579f) /* logf(0.030) */
#define LS_HI (-0.0833816f) /* logf(0.920) */
#define ALB_LO 0.015f       /* final compose clamp; below both layer clamps    */
#define ALB_HI 0.980f

/* ------------------------------------------------------------------ noise */
static unsigned hash2(int x, int y, unsigned seed) {
    unsigned h = (unsigned)x * 374761393u + (unsigned)y * 668265263u + seed * 2246822519u;
    h = (h ^ (h >> 13)) * 1274126177u;
    return h ^ (h >> 16);
}
static float rand01(unsigned h) { return (float)(h & 0xffffff) / 16777215.0f; }
static float lerpf(float a, float b, float t) { return a + (b - a) * t; }
static float fade(float t) { return t * t * (3.0f - 2.0f * t); }

static float value_noise(float x, float y, unsigned seed) {
    int xi = (int)floorf(x), yi = (int)floorf(y);
    float fx = x - xi, fy = y - yi;
    float v00 = rand01(hash2(xi, yi, seed));
    float v10 = rand01(hash2(xi + 1, yi, seed));
    float v01 = rand01(hash2(xi, yi + 1, seed));
    float v11 = rand01(hash2(xi + 1, yi + 1, seed));
    float u = fade(fx), v = fade(fy);
    return lerpf(lerpf(v00, v10, u), lerpf(v01, v11, u), v) * 2.0f - 1.0f;
}
static float fbm(float x, float y, int oct, float lac, float gain, unsigned seed) {
    float sum = 0.0f, amp = 1.0f, norm = 0.0f, freq = 1.0f;
    for (int i = 0; i < oct; i++) {
        sum += value_noise(x * freq, y * freq, seed + (unsigned)i * 17u) * amp;
        norm += amp; amp *= gain; freq *= lac;
    }
    return sum / norm;
}
static float smoothstepf(float e0, float e1, float x) {
    float t = (x - e0) / (e1 - e0);
    t = t < 0 ? 0 : (t > 1 ? 1 : t);
    return t * t * (3.0f - 2.0f * t);
}
static float clampf(float x, float lo, float hi) { return x < lo ? lo : (x > hi ? hi : x); }

/* ------------------------------------------------------------- gaussian */
static void gauss_kernel(float* k, int* r, float sigma) {
    *r = (int)ceilf(3.0f * sigma);
    float sum = 0.0f;
    for (int i = -(*r); i <= (*r); i++) {
        float v = expf(-0.5f * (i / sigma) * (i / sigma));
        k[i + *r] = v; sum += v;
    }
    for (int i = 0; i <= 2 * (*r); i++) k[i] /= sum;
}
static int reflect_i(int i, int n) {
    if (i < 0) i = -i - 1;
    if (i >= n) i = 2 * n - i - 1;
    return i < 0 ? 0 : (i >= n ? n - 1 : i);
}
static void blur_plane(const float* in, float* out, float* tmp, float sigma) {
    float k[128]; int r;
    gauss_kernel(k, &r, sigma);
    for (int y = 0; y < N; y++)
        for (int x = 0; x < N; x++) {
            float s = 0.0f;
            for (int i = -r; i <= r; i++) s += in[y * N + reflect_i(x + i, N)] * k[i + r];
            tmp[y * N + x] = s;
        }
    for (int y = 0; y < N; y++)
        for (int x = 0; x < N; x++) {
            float s = 0.0f;
            for (int i = -r; i <= r; i++) s += tmp[reflect_i(y + i, N) * N + x] * k[i + r];
            out[y * N + x] = s;
        }
}
/* blur 3 interleaved planes; tmp must hold 2*N*N floats */
static void blur_rgb(const float* in, float* out, float* tmp, float sigma) {
    float* plane = tmp;
    float* blurred = tmp + N * N;
    for (int c = 0; c < 3; c++) {
        for (int i = 0; i < N * N; i++) plane[i] = in[i * 3 + c];
        blur_plane(plane, blurred, plane, sigma);
        for (int i = 0; i < N * N; i++) out[i * 3 + c] = blurred[i];
    }
}

/* ------------------------------------------------------------------ scene */
typedef struct {
    float* LS;   /* structure layer, log reflectance, N*N*3  */
    float* lm;   /* micro layer, log units (relative), N*N*3 */
    float* A0;   /* nominal composed albedo, N*N*3           */
    float* H;    /* height field, N*N                        */
} Scene;

static void scene_build(Scene* sc, unsigned seed) {
    int n = N * N;
    sc->LS = calloc((size_t)n * 3, sizeof(float));
    sc->lm = calloc((size_t)n * 3, sizeof(float));
    sc->A0 = calloc((size_t)n * 3, sizeof(float));
    sc->H = calloc((size_t)n, sizeof(float));
    float* A = calloc((size_t)n * 3, sizeof(float));
    float* L = calloc((size_t)n * 3, sizeof(float));
    float* elev = calloc((size_t)n, sizeof(float));
    float* moist = calloc((size_t)n, sizeof(float));
    float* strata = calloc((size_t)n, sizeof(float));
    float* sn = calloc((size_t)n, sizeof(float));
    float* mraw = calloc((size_t)n * 3, sizeof(float));

    for (int y = 0; y < N; y++)
        for (int x = 0; x < N; x++) {
            int i = y * N + x;
            float fx = (float)x, fy = (float)y;
            float s1 = fbm(fx / 128.0f, fy / 128.0f, 4, 2.0f, 0.55f, seed);
            float s2 = fbm(fx / 64.0f + 11.3f, fy / 64.0f + 4.7f, 4, 2.0f, 0.55f, seed + 101);
            elev[i] = clampf((0.62f * s1 + 0.38f * s2 + 0.55f) / 1.17f, 0.0f, 1.0f);
            moist[i] = clampf((fbm(fx / 150.0f + 3.1f, fy / 150.0f - 8.2f, 3, 2.0f, 0.55f, seed + 211) + 0.6f) / 1.2f, 0.0f, 1.0f);
            strata[i] = fbm(fx / 40.0f, fy / 14.0f, 2, 2.0f, 0.5f, seed + 307);
            sn[i] = fbm(fx / 40.0f + 5.5f, fy / 40.0f + 9.1f, 2, 2.0f, 0.5f, seed + 401);
        }

    for (int i = 0; i < n; i++) {
        float e = elev[i], mo = moist[i], st = strata[i];
        float m_soil = smoothstepf(0.35f, 0.62f, mo) * (1.0f - smoothstepf(0.55f, 0.78f, e));
        float m_moss = smoothstepf(0.52f, 0.72f, mo) * smoothstepf(0.30f, 0.55f, e)
                     * (1.0f - smoothstepf(0.72f, 0.88f, e));
        float m_rock = fmaxf(0.0f, 1.0f - m_soil - m_moss);
        float wsum = m_rock + m_soil + m_moss;
        m_rock /= wsum; m_soil /= wsum; m_moss /= wsum;
        float m_snow = smoothstepf(0.70f, 0.86f, e + 0.07f * sn[i]);
        float cr = (0.34f + 0.05f * st) * m_rock + (0.27f + 0.03f * st) * m_soil
                 + (0.18f + 0.02f * st) * m_moss;
        float cg = (0.325f + 0.04f * st) * m_rock + (0.215f + 0.02f * st) * m_soil
                 + (0.27f + 0.04f * st) * m_moss;
        float cb = (0.30f + 0.03f * st) * m_rock + (0.15f + 0.02f * st) * m_soil
                 + (0.15f + 0.02f * st) * m_moss;
        A[i * 3 + 0] = lerpf(cr, 0.56f + 0.03f * st, m_snow);
        A[i * 3 + 1] = lerpf(cg, 0.575f + 0.03f * st, m_snow);
        A[i * 3 + 2] = lerpf(cb, 0.600f + 0.03f * st, m_snow);
        float lum = 0.2126f * A[i * 3] + 0.7152f * A[i * 3 + 1] + 0.0722f * A[i * 3 + 2];
        sc->H[i] = 70.0f * e + 8.0f * lum;
    }

    /* surface texture belongs to the structure layer (>= 14 px) */
    for (int c = 0; c < 3; c++)
        for (int y = 0; y < N; y++)
            for (int x = 0; x < N; x++) {
                int i = y * N + x;
                float f = fbm(x / 28.0f, y / 28.0f, 2, 2.0f, 0.50f, seed + 451u + (unsigned)c);
                A[i * 3 + c] *= (1.0f + 0.08f * f);
            }

    /* micro: fine grain, mostly per-channel independent, 1-3 px band.
       The shared surface texture lives in the structure layer instead, so a
       high-pass at a small scale is dominated by this grain and a high-pass
       at a larger scale picks up the shared texture - that crossover is what
       out/scale_scan.csv measures. */
    for (int c = 0; c < 3; c++) {
        const float inv = 1.0f / sqrtf(0.28f * 0.28f + 0.72f * 0.72f);
        for (int y = 0; y < N; y++)
            for (int x = 0; x < N; x++) {
                int i = y * N + x;
                float common = fbm(x / 2.2f, y / 2.2f, 1, 2.0f, 0.5f, seed + 700u);
                unsigned h = hash2(x, y, seed + 601u + (unsigned)c);
                float own = 0.75f * (rand01(h) * 2.0f - 1.0f)
                          + 0.25f * fbm(x / 1.8f, y / 1.8f, 2, 2.0f, 0.6f, seed + 800u + (unsigned)c);
                mraw[i * 3 + c] = (0.28f * common + 0.72f * own) * inv;
            }
    }
    /* per-channel zero-mean, unit-std, modulated by local brightness; then
       scale to the design sigma in log units and compose */
    for (int c = 0; c < 3; c++) {
        double sum = 0, sq = 0;
        for (int i = 0; i < n; i++) { sum += mraw[i * 3 + c]; sq += (double)mraw[i * 3 + c] * mraw[i * 3 + c]; }
        double mean = sum / n, var = sq / n - mean * mean, sd = var > 0 ? sqrt(var) : 1;
        for (int i = 0; i < n; i++) {
            float lum = 0.2126f * A[i * 3] + 0.7152f * A[i * 3 + 1] + 0.0722f * A[i * 3 + 2];
            mraw[i * 3 + c] = (float)((mraw[i * 3 + c] - mean) / sd) * (0.65f + 0.7f * (1.0f - lum));
        }
        sum = sq = 0;
        for (int i = 0; i < n; i++) { sum += mraw[i * 3 + c]; sq += (double)mraw[i * 3 + c] * mraw[i * 3 + c]; }
        mean = sum / n; var = sq / n - mean * mean; sd = var > 0 ? sqrt(var) : 1;
        for (int i = 0; i < n; i++) {
            float m = (float)((mraw[i * 3 + c] - mean) / sd) * M_SIGMA;
            A[i * 3 + c] = clampf(A[i * 3 + c] * expf(m), ALB_LO, ALB_HI);
        }
    }

    /* the split, measured from the finished albedo exactly as a user would */
    for (int i = 0; i < n * 3; i++) L[i] = logf(fmaxf(A[i], 1e-5f));
    float* tmp = calloc((size_t)n * 2, sizeof(float));
    blur_rgb(L, sc->LS, tmp, SIGMA_SPLIT);
    for (int i = 0; i < n * 3; i++) {
        sc->lm[i] = L[i] - sc->LS[i];              /* zero mean by construction */
        sc->A0[i] = expf(sc->LS[i] + sc->lm[i]);
    }

    free(tmp); free(A); free(L); free(elev); free(moist); free(strata); free(sn); free(mraw);
}

/* --------------------------------------------------------- layer faders */
static void layer_means(const Scene* sc, float mean[3]) {
    int n = N * N;
    for (int c = 0; c < 3; c++) mean[c] = 0;
    for (int i = 0; i < n; i++)
        for (int c = 0; c < 3; c++) mean[c] += sc->LS[i * 3 + c];
    for (int c = 0; c < 3; c++) mean[c] /= n;
}

/* the layered model: log-additive, each knob owns one band */
static void apply_layers(const Scene* sc, const float ks[3], const float km[3], float* out) {
    int n = N * N;
    float smean[3];
    layer_means(sc, smean);
    for (int i = 0; i < n; i++)
        for (int c = 0; c < 3; c++) {
            float s = smean[c] + ks[c] * (sc->LS[i * 3 + c] - smean[c]);
            s = clampf(s, LS_LO, LS_HI);
            float m = clampf(km[c] * sc->lm[i * 3 + c], -M_CLIP, M_CLIP);
            out[i * 3 + c] = clampf(expf(s + m), ALB_LO, ALB_HI);
        }
}

/* the merged single-layer scheme (status quo): one gain on the whole albedo.
   Every band scales together - there is no way to move detail without colour. */
static void apply_merged(const Scene* sc, const float kg[3], float* out) {
    int n = N * N;
    for (int i = 0; i < n; i++)
        for (int c = 0; c < 3; c++)
            out[i * 3 + c] = clampf(sc->A0[i * 3 + c] * kg[c], ALB_LO, ALB_HI);
}

/* linear (non-relative) micro layering, kept only as a counter-example */
static void apply_additive(const Scene* sc, const float ks[3], const float km[3], float* out) {
    int n = N * N;
    float smean[3];
    layer_means(sc, smean);
    for (int i = 0; i < n; i++)
        for (int c = 0; c < 3; c++) {
            float s0 = sc->LS[i * 3 + c];
            float s = clampf(smean[c] + ks[c] * (s0 - smean[c]), LS_LO, LS_HI);
            float base = expf(s);
            /* nominal absolute detail: the (0,1)-scaled difference is not
               right, so use exp(S') * (exp(m) - 1) evaluated at the *nominal*
               structure - this is what a merged pipeline bakes in and cannot
               rescale when the structure layer moves. */
            float detail = expf(clampf(s0, LS_LO, LS_HI)) * (expf(sc->lm[i * 3 + c]) - 1.0f);
            out[i * 3 + c] = clampf(base + km[c] * detail, ALB_LO, ALB_HI);
        }
}

/* --------------------------------------------------------------- weather */
typedef struct {
    const char* name;
    float sun[3];      /* all zero: no sun */
    float sky[3];
    float elev_deg, az_deg;
} Weather;

static const Weather WEATHER[3] = {
    { "clear",    { 3.05f, 2.80f, 2.40f }, { 0.33f, 0.47f, 0.76f }, 48.0f, 125.0f },
    { "overcast", { 0, 0, 0 },             { 0.92f, 0.97f, 1.03f },  0.0f,   0.0f },
    { "dusk",     { 1.25f, 0.58f, 0.26f }, { 0.11f, 0.17f, 0.36f },  7.0f, 240.0f },
};
/* neutral calibration light, D55-ish */
static const Weather NEUTRAL = { "neutral", { 2.30f, 2.20f, 2.05f }, { 0.52f, 0.55f, 0.62f }, 40.0f, 120.0f };
static const float LUMA[3] = { 0.2126f, 0.7152f, 0.0722f };

static void normal_at(const float* H, int x, int y, float n[3]) {
    int xm = x > 0 ? x - 1 : x, xp = x < N - 1 ? x + 1 : x;
    int ym = y > 0 ? y - 1 : y, yp = y < N - 1 ? y + 1 : y;
    float gx = (H[y * N + xp] - H[y * N + xm]) / (float)(xp - xm);
    float gy = (H[yp * N + x] - H[ym * N + x]) / (float)(yp - ym);
    float inv = 1.0f / sqrtf(gx * gx + gy * gy + 1.0f);
    n[0] = -gx * inv; n[1] = -gy * inv; n[2] = inv;
}

static float shadow_at(const float* H, int x, int y, float sx, float sy, float tan_e) {
    float hp = H[y * N + x];
    float inv = 1.0f / sqrtf(sx * sx + sy * sy);
    float dx = sx * inv, dy = sy * inv;
    for (float t = 1.25f; t < 130.0f; t += 1.25f) {
        float fx = x + dx * t, fy = y + dy * t;
        if (fx < 0 || fy < 0 || fx >= N - 1 || fy >= N - 1) break;
        int x0 = (int)fx, y0 = (int)fy;
        float tx = fx - x0, ty = fy - y0;
        float h = (H[y0 * N + x0] * (1 - tx) + H[y0 * N + x0 + 1] * tx) * (1 - ty)
                + (H[(y0 + 1) * N + x0] * (1 - tx) + H[(y0 + 1) * N + x0 + 1] * tx) * ty;
        if (h > hp + t * tan_e) return 0.0f;
    }
    return 1.0f;
}

static void irradiance(const float* H, const Weather* w, float* E) {
    int has_sun = w->sun[0] + w->sun[1] + w->sun[2] > 0;
    float sx = 0, sy = 0, sz = 0, tan_e = 0;
    if (has_sun) {
        float e = w->elev_deg * 3.14159265f / 180.0f, a = w->az_deg * 3.14159265f / 180.0f;
        sx = cosf(e) * cosf(a); sy = cosf(e) * sinf(a); sz = sinf(e);
        tan_e = tanf(e);
    }
    for (int y = 0; y < N; y++)
        for (int x = 0; x < N; x++) {
            int i = y * N + x;
            float n[3]; normal_at(H, x, y, n);
            float ndl = 0.0f, sh = 1.0f;
            if (has_sun) {
                ndl = fmaxf(0.0f, n[0] * sx + n[1] * sy + n[2] * sz);
                sh = shadow_at(H, x, y, sx, sy, tan_e);
            }
            float amb = has_sun ? (0.40f + 0.60f * n[2]) : (0.55f + 0.45f * n[2]);
            for (int c = 0; c < 3; c++)
                E[i * 3 + c] = (has_sun ? w->sun[c] * ndl * sh : 0.0f) + w->sky[c] * amb;
        }
}

static void wb_calibrate(const float* H, float wb[3]) {
    float* E = calloc((size_t)N * N * 3, sizeof(float));
    irradiance(H, &NEUTRAL, E);
    double mean[3] = { 0, 0, 0 };
    for (int i = 0; i < N * N; i++) for (int c = 0; c < 3; c++) mean[c] += E[i * 3 + c];
    for (int c = 0; c < 3; c++) mean[c] /= N * N;
    double y = mean[0] * LUMA[0] + mean[1] * LUMA[1] + mean[2] * LUMA[2];
    for (int c = 0; c < 3; c++) wb[c] = (float)(y / mean[c]);
    float g = wb[0] * LUMA[0] + wb[1] * LUMA[1] + wb[2] * LUMA[2];
    for (int c = 0; c < 3; c++) wb[c] /= g;
    free(E);
}

static float filmic(float x) {
    if (x <= 0) return 0;
    float v = (x * (2.51f * x + 0.03f)) / (x * (2.43f * x + 0.59f) + 0.14f);
    return v > 1 ? 1 : v;
}
static float srgb_encode(float x) {
    x = clampf(x, 0.0f, 1.0f);
    return x <= 0.0031308f ? 12.92f * x : 1.055f * powf(x, 1.0f / 2.4f) - 0.055f;
}

/* when > 0, forces the exposure (used to isolate tone-curve from auto-exposure) */
static float g_expo_override = 0.0f;

static void render_E(const float* A, const float* E, const float wb[3], float* D,
                     float* expo_out, float* clip_out) {
    float* Y = calloc((size_t)N * N, sizeof(float));
    double sumY = 0;
    for (int i = 0; i < N * N; i++) {
        float y = 0;
        for (int c = 0; c < 3; c++) y += A[i * 3 + c] * E[i * 3 + c] * LUMA[c];
        Y[i] = y; sumY += y;
    }
    float meanY = (float)(sumY / (N * N));
    float expo = meanY > 1e-9f ? 0.18f / meanY : 1.0f;
    if (g_expo_override > 0.0f) expo = g_expo_override;
    if (expo_out) *expo_out = expo;
    long clip = 0;
    for (int i = 0; i < N * N; i++) {
        float mag = Y[i] > 1e-9f ? filmic(Y[i] * expo) / (Y[i] * expo) : 0.0f;
        for (int c = 0; c < 3; c++) {
            float Lw = A[i * 3 + c] * E[i * 3 + c] * expo * wb[c];
            float v = srgb_encode(Lw * mag);
            if (v >= 0.999f) clip++;
            D[i * 3 + c] = v;
        }
    }
    if (clip_out) *clip_out = (float)clip / (N * N * 3);
    free(Y);
}

/* ---------------------------------------------------------------- stats */
static void channel_dev(const float* D, float dev[3]) {
    double mean[3] = { 0, 0, 0 };
    for (int i = 0; i < N * N; i++) for (int c = 0; c < 3; c++) mean[c] += D[i * 3 + c];
    for (int c = 0; c < 3; c++) mean[c] /= N * N;
    double y = mean[0] * LUMA[0] + mean[1] * LUMA[1] + mean[2] * LUMA[2];
    for (int c = 0; c < 3; c++) dev[c] = (float)(100.0 * (mean[c] / y - 1.0));
}

/* linear (post-WB, pre-tonemap) channel bias: the physical colour cast of the
   weather, uncontaminated by the highlight shoulder of the tone curve */
static void linear_dev(const float* A, const float* E, const float wb[3], float dev[3]) {
    int n = N * N;
    double mean[3] = { 0, 0, 0 };
    for (int i = 0; i < n; i++)
        for (int c = 0; c < 3; c++) mean[c] += A[i * 3 + c] * E[i * 3 + c] * wb[c];
    for (int c = 0; c < 3; c++) mean[c] /= n;
    double y = mean[0] * LUMA[0] + mean[1] * LUMA[1] + mean[2] * LUMA[2];
    for (int c = 0; c < 3; c++) dev[c] = (float)(100.0 * (mean[c] / y - 1.0));
}

typedef struct { float s_log[3]; float m_rel[3]; } BandStats;

/* log-domain band measurement of an albedo image:
   structure = std of gauss(log A) about its mean (log units)
   micro     = std of (log A - gauss(log A))            (log units)        */
static void albedo_bands(const float* A, float sigma, BandStats* bs) {
    int n = N * N;
    float* L = calloc((size_t)n * 3, sizeof(float));
    float* S = calloc((size_t)n * 3, sizeof(float));
    float* tmp = calloc((size_t)n * 2, sizeof(float));
    for (int i = 0; i < n * 3; i++) L[i] = logf(fmaxf(A[i], 1e-6f));
    blur_rgb(L, S, tmp, sigma);
    for (int c = 0; c < 3; c++) {
        double smean = 0, ss = 0, ms = 0;
        for (int i = 0; i < n; i++) smean += S[i * 3 + c];
        smean /= n;
        for (int i = 0; i < n; i++) {
            float s = S[i * 3 + c];
            float m = L[i * 3 + c] - s;
            ss += (s - smean) * (s - smean);
            ms += m * m;
        }
        bs->s_log[c] = (float)sqrt(ss / n);
        bs->m_rel[c] = (float)sqrt(ms / n);
    }
    free(L); free(S); free(tmp);
}

/* display-domain band measurement, for the observed (post-tonemap) picture */
static void display_bands(const float* D, float sigma, BandStats* bs) {
    int n = N * N;
    float* S = calloc((size_t)n * 3, sizeof(float));
    float* tmp = calloc((size_t)n * 2, sizeof(float));
    blur_rgb(D, S, tmp, sigma);
    for (int c = 0; c < 3; c++) {
        double smean = 0, ss = 0, ms = 0;
        for (int i = 0; i < n; i++) smean += S[i * 3 + c];
        smean /= n;
        for (int i = 0; i < n; i++) {
            float s = S[i * 3 + c];
            float m = (D[i * 3 + c] - s) / fmaxf(s, 1e-4f);
            ss += (s - smean) * (s - smean);
            ms += m * m;
        }
        bs->s_log[c] = (float)sqrt(ss / n);
        bs->m_rel[c] = (float)sqrt(ms / n);
    }
    free(S); free(tmp);
}

/* exact layer amplitudes, read from the two layer arrays themselves:
   structure = std of the (mean-locked, gain-scaled) structure layer, log units
   micro     = std of the (gain-scaled, clipped) micro layer, log units        */
static void layer_amplitudes(const Scene* sc, const float ks[3], const float km[3],
                             float s_amp[3], float m_amp[3]) {
    int n = N * N;
    float smean[3];
    layer_means(sc, smean);
    for (int c = 0; c < 3; c++) {
        double s2 = 0, m2 = 0;
        for (int i = 0; i < n; i++) {
            float sv = clampf(smean[c] + ks[c] * (sc->LS[i * 3 + c] - smean[c]), LS_LO, LS_HI);
            float mv = clampf(km[c] * sc->lm[i * 3 + c], -M_CLIP, M_CLIP);
            s2 += (sv - smean[c]) * (sv - smean[c]);
            m2 += mv * mv;
        }
        s_amp[c] = (float)sqrt(s2 / n);
        m_amp[c] = (float)sqrt(m2 / n);
    }
}

/* linear-reflectance band measurement.  This is the scheme-independent view:
   a single merged gain multiplies *every* band by the same factor, which is
   exactly what the split removes. */
static void linear_bands(const float* A, float sigma, float s_amp[3], float m_amp[3]) {
    int n = N * N;
    float* S = calloc((size_t)n * 3, sizeof(float));
    float* tmp = calloc((size_t)n * 2, sizeof(float));
    blur_rgb(A, S, tmp, sigma);
    for (int c = 0; c < 3; c++) {
        double smean = 0, ss = 0, ms = 0;
        for (int i = 0; i < n; i++) smean += S[i * 3 + c];
        smean /= n;
        for (int i = 0; i < n; i++) {
            float sv = S[i * 3 + c] - (float)smean;
            float mv = A[i * 3 + c] - S[i * 3 + c];
            ss += sv * sv;
            ms += mv * mv;
        }
        s_amp[c] = (float)sqrt(ss / n);
        m_amp[c] = (float)sqrt(ms / n);
    }
    free(S); free(tmp);
}

/* pair-mean cross-channel correlation of a high-pass band, in log albedo */
static float cross_corr_band(const float* A, float sigma) {
    int n = N * N;
    float* L = calloc((size_t)n * 3, sizeof(float));
    float* S = calloc((size_t)n * 3, sizeof(float));
    float* tmp = calloc((size_t)n * 2, sizeof(float));
    for (int i = 0; i < n * 3; i++) L[i] = logf(fmaxf(A[i], 1e-6f));
    blur_rgb(L, S, tmp, sigma);
    float* hp[3];
    for (int c = 0; c < 3; c++) hp[c] = calloc((size_t)n, sizeof(float));
    double mean[3] = { 0, 0, 0 }, sd[3] = { 0, 0, 0 };
    for (int c = 0; c < 3; c++) {
        for (int i = 0; i < n; i++) { hp[c][i] = L[i * 3 + c] - S[i * 3 + c]; mean[c] += hp[c][i]; }
        mean[c] /= n;
        for (int i = 0; i < n; i++) sd[c] += (hp[c][i] - mean[c]) * (hp[c][i] - mean[c]);
        sd[c] = sqrt(sd[c] / n);
    }
    double corr = 0;
    for (int a = 0; a < 3; a++)
        for (int b = a + 1; b < 3; b++) {
            double s = 0;
            for (int i = 0; i < n; i++)
                s += (hp[a][i] - mean[a]) * (hp[b][i] - mean[b]);
            corr += s / (n * sd[a] * sd[b]);
        }
    for (int c = 0; c < 3; c++) free(hp[c]);
    free(L); free(S); free(tmp);
    return (float)(corr / 3.0);
}

/* high-pass rms, in log units, at a given split scale */
static void highpass_sigma(const float* A, float sigma, float out[3]) {
    int n = N * N;
    float* L = calloc((size_t)n * 3, sizeof(float));
    float* S = calloc((size_t)n * 3, sizeof(float));
    float* tmp = calloc((size_t)n * 2, sizeof(float));
    for (int i = 0; i < n * 3; i++) L[i] = logf(fmaxf(A[i], 1e-6f));
    blur_rgb(L, S, tmp, sigma);
    for (int c = 0; c < 3; c++) {
        double s = 0;
        for (int i = 0; i < n; i++) {
            float h = L[i * 3 + c] - S[i * 3 + c];
            s += h * h;
        }
        out[c] = (float)sqrt(s / n);
    }
    free(L); free(S); free(tmp);
}

/* ---------------------------------------------------------------- output */
static void write_ppm(const char* path, const float* rgb, int w, int h) {
    FILE* f = fopen(path, "wb");
    if (!f) { fprintf(stderr, "cannot write %s\n", path); exit(1); }
    fprintf(f, "P6\n%d %d\n255\n", w, h);
    for (int i = 0; i < w * h * 3; i++) {
        int q = (int)(clampf(rgb[i], 0.0f, 1.0f) * 255.0f + 0.5f);
        fputc(q, f);
    }
    fclose(f);
}

static int cmp_float(const void* a, const void* b) {
    float x = *(const float*)a, y = *(const float*)b;
    return x < y ? -1 : (x > y ? 1 : 0);
}
static float percentile(float* sorted, int n, float p) {
    int i = (int)(p * (float)(n - 1));
    return sorted[i < 0 ? 0 : (i >= n ? n - 1 : i)];
}

/* --------------------------------------------------------------- tests */
static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static void bands_for(const Scene* sc, const float ks[3], const float km[3], int merged,
                      BandStats* bs) {
    int n = N * N;
    float* A = calloc((size_t)n * 3, sizeof(float));
    if (merged) apply_merged(sc, ks, A);
    else apply_layers(sc, ks, km, A);
    albedo_bands(A, SIGMA_SPLIT, bs);
    free(A);
}

static void test_scene(Scene* sc) {
    int n = N * N;
    float one[3] = { 1, 1, 1 };
    BandStats bs;
    bands_for(sc, one, one, 0, &bs);
    int bad = 0;
    double msum = 0;
    for (int i = 0; i < n; i++)
        for (int c = 0; c < 3; c++) {
            if (sc->LS[i * 3 + c] < LS_LO - 1e-5f || sc->LS[i * 3 + c] > LS_HI + 1e-5f) bad++;
            msum += sc->lm[i * 3 + c];
        }
    check("structure_within_clip", bad == 0);
    check("micro_zero_mean", fabs(msum / (n * 3)) < 1e-4);
    /* the measured micro band also contains the 28 px texture's knee, so it
       is slightly above the 4.5 % grain design sigma */
    check("micro_sigma_near_design", bs.m_rel[0] > 0.8f * M_SIGMA && bs.m_rel[0] < 1.5f * M_SIGMA);
    check("structure_sigma_exceeds_micro", bs.s_log[0] > 2.5f * bs.m_rel[0]);
}

static void test_knobs(Scene* sc) {
    int n = N * N;
    float one[3] = { 1, 1, 1 };
    float s0[3], m0[3];

    /* 1. the model claim: each knob owns exactly one layer array */
    layer_amplitudes(sc, one, one, s0, m0);
    {
        float ks[3] = { 1.20f, 1.20f, 1.20f }, s1[3], m1[3];
        layer_amplitudes(sc, ks, one, s1, m1);
        check("structure_knob_moves_structure",
              fabsf(s1[0] / s0[0] - 1.20f) < 1e-4f && fabsf(s1[2] / s0[2] - 1.20f) < 1e-4f);
        check("structure_knob_leaves_micro",
              fabsf(m1[0] / m0[0] - 1.0f) < 1e-6f && fabsf(m1[2] / m0[2] - 1.0f) < 1e-6f);
    }
    {
        float km[3] = { 1.20f, 1.20f, 1.20f }, s1[3], m1[3];
        layer_amplitudes(sc, one, km, s1, m1);
        check("micro_knob_moves_micro",
              fabsf(m1[0] / m0[0] - 1.20f) < 0.03f && fabsf(m1[2] / m0[2] - 1.20f) < 0.03f);
        check("micro_knob_leaves_structure",
              fabsf(s1[0] / s0[0] - 1.0f) < 1e-6f && fabsf(s1[2] / s0[2] - 1.0f) < 1e-6f);
    }
    {
        float ks_b[3] = { 1.0f, 1.0f, 0.85f }, s1[3], m1[3];
        layer_amplitudes(sc, ks_b, one, s1, m1);
        check("structure_knob_per_channel",
              fabsf(s1[0] / s0[0] - 1.0f) < 1e-6f && fabsf(s1[2] / s0[2] - 0.85f) < 1e-4f);
        float km_g[3] = { 1.0f, 1.20f, 1.0f }, s2[3], m2[3];
        layer_amplitudes(sc, one, km_g, s2, m2);
        check("micro_knob_per_channel",
              fabsf(m2[0] / m0[0] - 1.0f) < 1e-6f && fabsf(m2[1] / m0[1] - 1.20f) < 0.03f);
    }

    /* 2. the merged scheme, measured the same way in linear reflectance:
       one gain moves both bands by exactly that gain */
    {
        float* A = calloc((size_t)n * 3, sizeof(float));
        float sA[3], mA[3], sB[3], mB[3];
        apply_merged(sc, one, A);
        linear_bands(A, SIGMA_SPLIT, sA, mA);
        float kg[3] = { 1.20f, 1.0f, 1.0f };
        apply_merged(sc, kg, A);
        linear_bands(A, SIGMA_SPLIT, sB, mB);
        check("merged_scheme_couples_bands",
              fabsf(sB[0] / sA[0] - 1.20f) < 0.01f && fabsf(mB[0] / mA[0] - 1.20f) < 0.01f);
        free(A);
    }

    /* 3. composed albedo, same knob: the *relative* micro contrast is the
       invariant (the multiplicative layer scales absolute detail with local
       brightness by design, so linear micro grows with the structure band) */
    {
        float* A = calloc((size_t)n * 3, sizeof(float));
        BandStats ba, bb;
        apply_layers(sc, one, one, A);
        albedo_bands(A, SIGMA_SPLIT, &ba);
        float ks[3] = { 1.20f, 1.0f, 1.0f };
        apply_layers(sc, ks, one, A);
        albedo_bands(A, SIGMA_SPLIT, &bb);
        check("structure_knob_keeps_relative_micro_in_compose",
              fabsf(bb.m_rel[0] / ba.m_rel[0] - 1.0f) < 0.05f
              && fabsf(bb.m_rel[2] / ba.m_rel[2] - 1.0f) < 0.01f);
        check("structure_knob_moves_composed_structure",
              bb.s_log[0] / ba.s_log[0] > 1.10f && bb.s_log[2] / ba.s_log[2] > 0.98f
              && bb.s_log[2] / ba.s_log[2] < 1.02f);
        free(A);
    }

    /* 4. linear additive layering cross-talks: the detail is stored as an
       absolute amount, so raising the structure by +20 % dilutes the relative
       micro contrast by 1/1.2 in linear reflectance.  In log space the same
       error shows up as an unchanged relative micro (log division hides it),
       which is exactly why the log-domain split is the right home for it. */
    {
        float ks[3] = { 1.20f, 1.20f, 1.20f };
        float* A = calloc((size_t)n * 3, sizeof(float));
        float sA[3], mA[3], sB[3], mB[3];
        apply_additive(sc, one, one, A);
        linear_bands(A, SIGMA_SPLIT, sA, mA);
        apply_additive(sc, ks, one, A);
        linear_bands(A, SIGMA_SPLIT, sB, mB);
        /* absolute detail is unchanged, base grows: relative micro falls by 1/1.2 */
        float relA = mA[0] / sA[0], relB = mB[0] / sB[0];
        check("additive_layering_cross_talks", relB / relA < 0.87f && mB[0] / mA[0] > 0.97f);
        free(A);
    }
}

static void test_weather(Scene* sc) {
    int n = N * N;
    float wb[3];
    wb_calibrate(sc->H, wb);
    check("wb_sane", wb[0] > 0.8f && wb[0] < 1.3f && wb[2] > 0.8f && wb[2] < 1.3f);

    float* A = calloc((size_t)n * 3, sizeof(float));
    float* E = calloc((size_t)n * 3, sizeof(float));
    float* D = calloc((size_t)n * 3, sizeof(float));
    for (int i = 0; i < n * 3; i++) A[i] = 0.5f;
    irradiance(sc->H, &NEUTRAL, E);
    render_E(A, E, wb, D, NULL, NULL);
    float dev[3];
    channel_dev(D, dev);
    check("neutral_plane_gray", fabsf(dev[0]) < 0.5f && fabsf(dev[1]) < 0.5f && fabsf(dev[2]) < 0.5f);

    float one[3] = { 1, 1, 1 };
    apply_layers(sc, one, one, A);
    int finite = 1;
    for (int wi = 0; wi < 3; wi++) {
        irradiance(sc->H, &WEATHER[wi], E);
        render_E(A, E, wb, D, NULL, NULL);
        channel_dev(D, dev);
        for (int c = 0; c < 3; c++) if (!isfinite(dev[c])) finite = 0;
    }
    check("weather_outputs_finite", finite);
    irradiance(sc->H, &WEATHER[0], E);
    render_E(A, E, wb, D, NULL, NULL);
    float d1[3]; channel_dev(D, d1);
    render_E(A, E, wb, D, NULL, NULL);
    float d2[3]; channel_dev(D, d2);
    check("render_deterministic", d1[0] == d2[0] && d1[2] == d2[2]);
    free(A); free(E); free(D);
}

/* --------------------------------------------------------------- build */
static void build_outputs(Scene* sc) {
    int n = N * N;
    char path[256];
    FILE* f;
    float wb[3];
    wb_calibrate(sc->H, wb);
    const char* chname[3] = { "R", "G", "B" };
    float one[3] = { 1, 1, 1 }, zero[3] = { 0, 0, 0 };

    float* E[3];
    for (int wi = 0; wi < 3; wi++) {
        E[wi] = calloc((size_t)n * 3, sizeof(float));
        irradiance(sc->H, &WEATHER[wi], E[wi]);
    }

    /* --- 1. scale scan: where shared surface texture ends, grain begins --- */
    f = fopen("out/scale_scan.csv", "w");
    fprintf(f, "sigma_px,hp_cross_corr,hp_sigma_r,hp_sigma_g,hp_sigma_b\n");
    for (float sig = 0.5f; sig <= 12.05f; sig += 0.5f) {
        float hs[3];
        highpass_sigma(sc->A0, sig, hs);
        fprintf(f, "%.1f,%.4f,%.5f,%.5f,%.5f\n", sig, cross_corr_band(sc->A0, sig),
                hs[0], hs[1], hs[2]);
    }
    fclose(f);

    /* --- 1b. split-scale decision: analytic Gaussian transfer at the scene's
       own feature periods.  The criteria are stated in the table header:
         structure must keep >= 80 % of the finest *authored* structure feature
         (28 px period, the fine surface texture);
         micro must keep >= 99 % of the finest *grain* (1.5-3 px).
       Together those bound sigma to roughly 1.0 .. 2.75 px; 2.75 px is the
       largest value that still satisfies the structure criterion (82.7 %), so
       it is chosen: maximal separation, minimal leakage of authored texture. --- */
    {
        const float periods[5] = { 96.0f, 28.0f, 14.0f, 3.0f, 1.5f };
        const char* names[5] = { "material_96px", "fine_texture_28px",
                                 "texture_octave_14px", "grain_3px", "grain_fine_1.5px" };
        f = fopen("out/scale_decision.csv", "w");
        fprintf(f, "period_px,feature,sigma_px,structure_retained_pct,micro_retained_pct,"
                   "requirement,pass\n");
        for (int si = 0; si < 16; si++) {
            float sig = 1.0f + 0.25f * si;
            for (int k = 0; k < 5; k++) {
                double fr = 1.0 / periods[k];
                double hs = exp(-2.0 * 3.14159265358979 * 3.14159265358979 * sig * sig * fr * fr);
                double hm = 1.0 - hs;
                const char* req = k == 1 ? "structure>=80%" : (k >= 3 ? "micro>=99%" : "info");
                int pass = k == 1 ? hs >= 0.80 : (k >= 3 ? hm >= 0.99 : 1);
                fprintf(f, "%.1f,%s,%.2f,%.2f,%.2f,%s,%d\n", periods[k], names[k], sig,
                        100.0 * hs, 100.0 * hm, req, pass);
            }
        }
        fclose(f);
    }

    /* --- 2. layer maps --- */
    /* structure layer as a viewable reflectance image: exp(LS) */
    {
        float* Sv = malloc((size_t)n * 3 * sizeof(float));
        for (int i = 0; i < n * 3; i++) Sv[i] = expf(sc->LS[i]);
        write_ppm("out/structure.ppm", Sv, N, N);
        free(Sv);
    }
    float* Av = calloc((size_t)n * 3, sizeof(float));
    for (int i = 0; i < n * 3; i++) Av[i] = 0.5f + sc->lm[i] / (2.0f * M_CLIP);
    write_ppm("out/micro.ppm", Av, N, N);
    apply_layers(sc, one, one, Av);
    write_ppm("out/albedo.ppm", Av, N, N);

    /* histograms, exported at full float precision (the 8-bit PNGs quantise
       too coarsely to show the tails honestly) */
    {
        const int NB = 256;
        float* v = malloc((size_t)n * sizeof(float));
        f = fopen("out/hist_structure.csv", "w");
        fprintf(f, "channel,bin_lo,bin_hi,count\n");
        for (int c = 0; c < 3; c++) {
            float lo = 1e9f, hi = -1e9f;
            for (int i = 0; i < n; i++) {
                v[i] = expf(clampf(sc->LS[i * 3 + c], LS_LO, LS_HI));
                lo = fminf(lo, v[i]); hi = fmaxf(hi, v[i]);
            }
            lo = fminf(lo, S_REFL_LO) - 0.005f;
            hi = fmaxf(hi, S_REFL_HI) + 0.005f;
            float w = (hi - lo) / NB;
            int* h = calloc(NB, sizeof(int));
            for (int i = 0; i < n; i++) {
                int b = (int)((v[i] - lo) / w);
                h[b < 0 ? 0 : (b >= NB ? NB - 1 : b)]++;
            }
            for (int b = 0; b < NB; b++)
                fprintf(f, "%s,%.4f,%.4f,%d\n", chname[c], lo + b * w, lo + (b + 1) * w, h[b]);
            free(h);
        }
        fclose(f);

        f = fopen("out/hist_micro.csv", "w");
        fprintf(f, "channel,bin_lo,bin_hi,count\n");
        for (int c = 0; c < 3; c++) {
            float lo = 1e9f, hi = -1e9f;
            for (int i = 0; i < n; i++) {
                v[i] = 100.0f * (expf(clampf(sc->lm[i * 3 + c], -M_CLIP, M_CLIP)) - 1.0f);
                lo = fminf(lo, v[i]); hi = fmaxf(hi, v[i]);
            }
            float w = (hi - lo) / NB;
            int* h = calloc(NB, sizeof(int));
            for (int i = 0; i < n; i++) {
                int b = (int)((v[i] - lo) / w);
                h[b < 0 ? 0 : (b >= NB ? NB - 1 : b)]++;
            }
            for (int b = 0; b < NB; b++)
                fprintf(f, "%s,%.3f,%.3f,%d\n", chname[c], lo + b * w, lo + (b + 1) * w, h[b]);
            free(h);
        }
        fclose(f);
        free(v);
    }

    /* truncation masks at nominal gain: white where a layer is pinned at a
       bound.  Structure is 0 % at nominal, micro is a few tenths of a percent,
       so these are sparse - which is the point: truncation is local, not a
       global tone shift. */
    {
        float* m = calloc((size_t)n * 3, sizeof(float));
        for (int i = 0; i < n; i++)
            for (int c = 0; c < 3; c++) {
                int sclip = (sc->LS[i * 3 + c] <= LS_LO + 1e-6f) || (sc->LS[i * 3 + c] >= LS_HI - 1e-6f);
                m[i * 3 + c] = sclip ? 1.0f : 0.0f;
            }
        write_ppm("out/structure_clip.ppm", m, N, N);
        for (int i = 0; i < n; i++)
            for (int c = 0; c < 3; c++) {
                int mc = (sc->lm[i * 3 + c] <= -M_CLIP + 1e-6f) || (sc->lm[i * 3 + c] >= M_CLIP - 1e-6f);
                m[i * 3 + c] = mc ? 1.0f : 0.0f;
            }
        write_ppm("out/micro_clip.ppm", m, N, N);
        free(m);
    }

    /* --- 3. layer amplitude and truncation statistics --- */
    f = fopen("out/layer_stats.csv", "w");
    fprintf(f, "layer,channel,unit,min,p001,p01,p50,p99,p999,max,clip_frac,sigma\n");
    {
        float* v = malloc((size_t)n * sizeof(float));
        for (int c = 0; c < 3; c++) {
            /* structure layer as composed (gain 1, clipped at its bounds) */
            for (int i = 0; i < n; i++) v[i] = expf(clampf(sc->LS[i * 3 + c], LS_LO, LS_HI));
            long cl = 0;
            for (int i = 0; i < n; i++) if (sc->LS[i * 3 + c] <= LS_LO + 1e-6f || sc->LS[i * 3 + c] >= LS_HI - 1e-6f) cl++;
            qsort(v, n, sizeof(float), cmp_float);
            double mean = 0, sq = 0;
            for (int i = 0; i < n; i++) { mean += v[i]; sq += (double)v[i] * v[i]; }
            mean /= n;
            fprintf(f, "structure,%s,reflectance,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n",
                    chname[c], v[0], percentile(v, n, 0.001f), percentile(v, n, 0.01f),
                    percentile(v, n, 0.5f), percentile(v, n, 0.99f), percentile(v, n, 0.999f),
                    v[n - 1], (double)cl / n, (float)sqrt(sq / n - mean * mean));
            /* micro layer as composed (gain 1, clipped at +/-M_CLIP) */
            for (int i = 0; i < n; i++)
                v[i] = 100.0f * (expf(clampf(sc->lm[i * 3 + c], -M_CLIP, M_CLIP)) - 1.0f);
            cl = 0;
            for (int i = 0; i < n; i++) if (sc->lm[i * 3 + c] <= -M_CLIP + 1e-6f
                                          || sc->lm[i * 3 + c] >= M_CLIP - 1e-6f) cl++;
            qsort(v, n, sizeof(float), cmp_float);
            mean = sq = 0;
            for (int i = 0; i < n; i++) { mean += v[i]; sq += (double)v[i] * v[i]; }
            mean /= n;
            fprintf(f, "micro,%s,percent,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.5f,%.4f\n",
                    chname[c], v[0], percentile(v, n, 0.001f), percentile(v, n, 0.01f),
                    percentile(v, n, 0.5f), percentile(v, n, 0.99f), percentile(v, n, 0.999f),
                    v[n - 1], (double)cl / n, (float)sqrt(sq / n - mean * mean));
        }
        free(v);
    }
    fclose(f);

    /* --- 3b. truncation onset: at what knob position does each layer's clip
       start to bite?  clip_frac = fraction of pixels pinned at the bound --- */
    f = fopen("out/truncation.csv", "w");
    fprintf(f, "layer,channel,gain,clip_frac\n");
    for (int layer = 0; layer < 2; layer++)
        for (int c = 0; c < 3; c++)
            for (int gi = 0; gi < 21; gi++) {
                float g = 0.80f + 0.05f * gi;
                long cl = 0;
                float smean[3]; layer_means(sc, smean);
                for (int i = 0; i < N * N; i++) {
                    float sv = smean[c] + g * (sc->LS[i * 3 + c] - smean[c]);
                    if (layer == 0 && (sv <= LS_LO + 1e-6f || sv >= LS_HI - 1e-6f)) cl++;
                    float mv = g * sc->lm[i * 3 + c];
                    if (layer == 1 && (mv <= -M_CLIP + 1e-6f || mv >= M_CLIP - 1e-6f)) cl++;
                }
                fprintf(f, "%s,%s,%.2f,%.5f\n", layer == 0 ? "structure" : "micro",
                        chname[c], g, (double)cl / (N * N));
            }
    fclose(f);

    /* --- 4. per-weather baseline, renders, channel deviation --- */
    float* D = calloc((size_t)n * 3, sizeof(float));
    float base_dev[3][3], base_lin[3][3];
    BandStats base_ab[3], base_db[3];
    f = fopen("out/weather.csv", "w");
    fprintf(f, "weather,sun_elev,sun_az,channel,mean_display,dev_pct,dev_linear_pct,expo,clip_frac,"
               "albedo_structure_log,albedo_micro_log,display_structure,display_micro_rel\n");
    for (int wi = 0; wi < 3; wi++) {
        float expo = 0, clipf = 0;
        apply_layers(sc, one, one, Av);
        render_E(Av, E[wi], wb, D, &expo, &clipf);
        snprintf(path, sizeof path, "out/weather_%s.ppm", WEATHER[wi].name);
        write_ppm(path, D, N, N);
        channel_dev(D, base_dev[wi]);
        linear_dev(Av, E[wi], wb, base_lin[wi]);
        albedo_bands(Av, SIGMA_SPLIT, &base_ab[wi]);
        display_bands(D, SIGMA_SPLIT, &base_db[wi]);
        double mean[3] = { 0, 0, 0 };
        for (int i = 0; i < n; i++) for (int c = 0; c < 3; c++) mean[c] += D[i * 3 + c];
        for (int c = 0; c < 3; c++) mean[c] /= n;
        for (int c = 0; c < 3; c++)
            fprintf(f, "%s,%.1f,%.1f,%s,%.5f,%.3f,%.3f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n",
                    WEATHER[wi].name, WEATHER[wi].elev_deg, WEATHER[wi].az_deg, chname[c],
                    mean[c], base_dev[wi][c], base_lin[wi][c], expo, clipf,
                    base_ab[wi].s_log[c], base_ab[wi].m_rel[c],
                    base_db[wi].s_log[c], base_db[wi].m_rel[c]);
    }
    fclose(f);

    /* --- 5. knob independence.
       Two measurements are written, because they answer different questions:
         knob_independence.csv - the stored layer amplitudes (the model);
         resplit_check.csv     - the same knobs measured by re-splitting the
                                 recomposed albedo.  The second carries the
                                 split's own leakage (G*G != G), which is a
                                 measurement floor, not cross-talk. --- */
    {
        BandStats base;
        bands_for(sc, one, one, 0, &base);
        f = fopen("out/knob_independence.csv", "w");
        fprintf(f, "layer,gain,structure_ratio,micro_ratio\n");
        for (int layer = 0; layer < 2; layer++)
            for (int gi = 0; gi < 13; gi++) {
                float g = 0.70f + 0.05f * gi;
                float ks[3] = { 1, 1, 1 }, km[3] = { 1, 1, 1 };
                if (layer == 0) { ks[0] = ks[1] = ks[2] = g; } else { km[0] = km[1] = km[2] = g; }
                float s1[3], m1[3], s0[3], m0[3];
                layer_amplitudes(sc, one, one, s0, m0);
                layer_amplitudes(sc, ks, km, s1, m1);
                float sr = 0, mr = 0;
                for (int c = 0; c < 3; c++) { sr += s1[c] / s0[c]; mr += m1[c] / m0[c]; }
                fprintf(f, "%s,%.2f,%.4f,%.4f\n", layer == 0 ? "structure" : "micro", g,
                        sr / 3.0f, mr / 3.0f);
            }
        fclose(f);

        f = fopen("out/resplit_check.csv", "w");
        fprintf(f, "layer,gain,structure_ratio,micro_ratio\n");
        for (int layer = 0; layer < 2; layer++)
            for (int gi = 0; gi < 13; gi++) {
                float g = 0.70f + 0.05f * gi;
                float ks[3] = { 1, 1, 1 }, km[3] = { 1, 1, 1 };
                if (layer == 0) { ks[0] = ks[1] = ks[2] = g; } else { km[0] = km[1] = km[2] = g; }
                BandStats bs;
                bands_for(sc, ks, km, 0, &bs);
                float sr = 0, mr = 0;
                for (int c = 0; c < 3; c++) { sr += bs.s_log[c] / base.s_log[c]; mr += bs.m_rel[c] / base.m_rel[c]; }
                fprintf(f, "%s,%.2f,%.4f,%.4f\n", layer == 0 ? "structure" : "micro", g,
                        sr / 3.0f, mr / 3.0f);
            }
        fclose(f);
    }

    /* --- 5b. cross-term summary, same-channel, both layers on --- */
    {
        float* A = calloc((size_t)n * 3, sizeof(float));
        float* Dd = calloc((size_t)n * 3, sizeof(float));
        f = fopen("out/cross_terms.csv", "w");
        fprintf(f, "weather,space,max_cross_pct,worst_channel,worst_knob,definition\n");
        for (int wi = 0; wi < 3; wi++) {
            /* albedo space: read the stored layers exactly */
            BandStats ab0;
            apply_layers(sc, one, one, A);
            albedo_bands(A, SIGMA_SPLIT, &ab0);
            float worst_ab = 0; int wc = 0;
            for (int ci = 0; ci < 3; ci++) {
                float ks[3] = { 1, 1, 1 }; ks[ci] = 1.20f;
                BandStats ab1;
                apply_layers(sc, ks, one, A);
                albedo_bands(A, SIGMA_SPLIT, &ab1);
                for (int c = 0; c < 3; c++)
                    if (fabsf(ab1.m_rel[c] / ab0.m_rel[c] - 1.0f) > worst_ab) {
                        worst_ab = fabsf(ab1.m_rel[c] / ab0.m_rel[c] - 1.0f); wc = c;
                    }
            }
            fprintf(f, "%s,albedo,%.3f,%s,structure+20%%,micro band change when the structure knob moves\n",
                    WEATHER[wi].name, worst_ab * 100.0f, chname[wc]);

            /* display space: after weather + tone curve, same measurement */
            BandStats db0;
            apply_layers(sc, one, one, A);
            render_E(A, E[wi], wb, Dd, NULL, NULL);
            display_bands(Dd, SIGMA_SPLIT, &db0);
            float worst_db = 0; int wc2 = 0;
            for (int ci = 0; ci < 3; ci++) {
                float ks[3] = { 1, 1, 1 }; ks[ci] = 1.20f;
                BandStats db1;
                apply_layers(sc, ks, one, A);
                render_E(A, E[wi], wb, Dd, NULL, NULL);
                display_bands(Dd, SIGMA_SPLIT, &db1);
                for (int c = 0; c < 3; c++)
                    if (fabsf(db1.m_rel[c] / db0.m_rel[c] - 1.0f) > worst_db) {
                        worst_db = fabsf(db1.m_rel[c] / db0.m_rel[c] - 1.0f); wc2 = c;
                    }
            }
            fprintf(f, "%s,display,%.3f,%s,structure+20%%,same, after the weather tone curve\n",
                    WEATHER[wi].name, worst_db * 100.0f, chname[wc2]);

            /* same measurement with the exposure frozen at the nominal scene's
               value: isolates the tone-curve shoulder from auto-exposure */
            float expo0 = 0;
            apply_layers(sc, one, one, A);
            render_E(A, E[wi], wb, Dd, &expo0, NULL);
            g_expo_override = expo0;
            display_bands(Dd, SIGMA_SPLIT, &db0);
            float worst_fx = 0; int wc3 = 0;
            for (int ci = 0; ci < 3; ci++) {
                float ks[3] = { 1, 1, 1 }; ks[ci] = 1.20f;
                BandStats db1;
                apply_layers(sc, ks, one, A);
                render_E(A, E[wi], wb, Dd, NULL, NULL);
                display_bands(Dd, SIGMA_SPLIT, &db1);
                for (int c = 0; c < 3; c++)
                    if (fabsf(db1.m_rel[c] / db0.m_rel[c] - 1.0f) > worst_fx) {
                        worst_fx = fabsf(db1.m_rel[c] / db0.m_rel[c] - 1.0f); wc3 = c;
                    }
            }
            g_expo_override = 0.0f;
            fprintf(f, "%s,display_fixed_exposure,%.3f,%s,structure+20%%,same, exposure frozen\n",
                    WEATHER[wi].name, worst_fx * 100.0f, chname[wc3]);
            apply_layers(sc, one, one, A);
            render_E(A, E[wi], wb, Dd, NULL, NULL);
            display_bands(Dd, SIGMA_SPLIT, &db0);
            /* one-sided check: the micro band is allowed to *shrink* when the
               structure knob rises (message above); what would be a bug is it
               growing.  Report the largest positive change separately. */
            float grow = 0; int wcg = 0;
            for (int ci = 0; ci < 3; ci++) {
                float ks[3] = { 1, 1, 1 }; ks[ci] = 1.20f;
                BandStats db1;
                apply_layers(sc, ks, one, A);
                render_E(A, E[wi], wb, Dd, NULL, NULL);
                display_bands(Dd, SIGMA_SPLIT, &db1);
                for (int c = 0; c < 3; c++) {
                    float ch = (db1.m_rel[c] / db0.m_rel[c] - 1.0f) * 100.0f;
                    if (ch > grow) { grow = ch; wcg = c; }
                }
            }
            fprintf(f, "%s,display_micro_growth,%.3f,%s,structure+20%%,same, positive change only\n",
                    WEATHER[wi].name, grow, chname[wcg]);
            /* wide-range robustness: structure +50 % in R and G (B hits its
               upper clip at +50 %, so B is excluded from this test) */
            float worst50 = 0; int wc50 = 0;
            for (int ci = 0; ci < 2; ci++) {
                float ks[3] = { 1, 1, 1 }; ks[ci] = 1.50f;
                BandStats ab1;
                apply_layers(sc, ks, one, A);
                albedo_bands(A, SIGMA_SPLIT, &ab1);
                for (int c = 0; c < 3; c++) {
                    float ch = fabsf(ab1.m_rel[c] / ab0.m_rel[c] - 1.0f) * 100.0f;
                    if (ch > worst50) { worst50 = ch; wc50 = c; }
                }
            }
            fprintf(f, "%s,albedo,%.3f,%s,structure+50%%,micro band change, R/G only\n",
                    WEATHER[wi].name, worst50, chname[wc50]);
        }
        fclose(f);
        free(A); free(Dd);
    }

    /* --- 6. merged single-layer scheme: band coupling, other-channel drift --- */
    float base_lin_s[3], base_lin_m[3];
    {
        float* A = calloc((size_t)n * 3, sizeof(float));
        apply_merged(sc, one, A);
        linear_bands(A, SIGMA_SPLIT, base_lin_s, base_lin_m);
        free(A);
    }
    f = fopen("out/merged_scheme.csv", "w");
    fprintf(f, "weather,driven_channel,gain,channel,dev_change_pct,structure_ratio_linear,micro_ratio_linear\n");
    for (int wi = 0; wi < 3; wi++)
        for (int ci = 0; ci < 3; ci++)
            for (int gi = 0; gi < 13; gi++) {
                float g = 0.70f + 0.05f * gi;
                float kg[3] = { 1, 1, 1 };
                kg[ci] = g;
                apply_merged(sc, kg, Av);
                /* linear-reflectance bands: a merged gain scales both by g */
                float s1[3], m1[3];
                linear_bands(Av, SIGMA_SPLIT, s1, m1);
                render_E(Av, E[wi], wb, D, NULL, NULL);
                float dev[3];
                channel_dev(D, dev);
                for (int c = 0; c < 3; c++)
                    fprintf(f, "%s,%s,%.2f,%s,%.3f,%.4f,%.4f\n",
                            WEATHER[wi].name, chname[ci], g, chname[c],
                            dev[c] - base_dev[wi][c],
                            s1[c] / base_lin_s[c],
                            m1[c] / base_lin_m[c]);
            }
    fclose(f);

    /* --- 7. layered scheme sweeps.  Each knob is pushed one channel at a
       time.  Reported per row: the driven channel's display deviation, the
       same band amplitude in display space, and the *other* band's amplitude
       ratio - the crosstalk number the report quotes. --- */
    {
        float* A = calloc((size_t)n * 3, sizeof(float));
        float* Dd = calloc((size_t)n * 3, sizeof(float));
        for (int knob = 0; knob < 2; knob++) {
            f = fopen(knob == 0 ? "out/structure_sweep.csv" : "out/micro_sweep.csv", "w");
            fprintf(f, "weather,driven_channel,gain,channel,dev_pct,dev_change_pct,"
                       "display_structure,display_micro_rel,cross_band_change_pct\n");
            for (int wi = 0; wi < 3; wi++)
                for (int ci = 0; ci < 3; ci++)
                    for (int gi = 0; gi < 13; gi++) {
                        float g = 0.70f + 0.05f * gi;
                        float ks[3] = { 1, 1, 1 }, km[3] = { 1, 1, 1 };
                        if (knob == 0) ks[ci] = g; else km[ci] = g;
                        apply_layers(sc, ks, km, A);
                        render_E(A, E[wi], wb, Dd, NULL, NULL);
                        BandStats bs;
                        float dev[3];
                        channel_dev(Dd, dev);
                        display_bands(Dd, SIGMA_SPLIT, &bs);
                        for (int c = 0; c < 3; c++) {
                            /* the band the knob is *not* supposed to move */
                            float cross = knob == 0
                                ? (bs.m_rel[c] / base_db[wi].m_rel[c] - 1.0f) * 100.0f
                                : (bs.s_log[c] / base_db[wi].s_log[c] - 1.0f) * 100.0f;
                            fprintf(f, "%s,%s,%.2f,%s,%.3f,%.3f,%.5f,%.5f,%.3f\n",
                                    WEATHER[wi].name, chname[ci], g, chname[c], dev[c],
                                    dev[c] - base_dev[wi][c],
                                    bs.s_log[c], bs.m_rel[c], cross);
                        }
                    }
            fclose(f);
        }
        free(A); free(Dd);
    }

    /* --- 8. interference sample: both layers on, one layer pushed +20 %,
       the other band measured on the *same channel*.  Any change there is the
       cross term.  Measured in albedo space (the model) and in display space
       (after weather + tone curve), for each weather. --- */
    f = fopen("out/interference.csv", "w");
    fprintf(f, "weather,model,knob_layer,knob_channel,response_layer,response_channel,"
               "base,perturbed,cross_term_pct\n");
    {
        float* A = calloc((size_t)n * 3, sizeof(float));
        float* D = calloc((size_t)n * 3, sizeof(float));
        BandStats base_ab, base_db;
        apply_layers(sc, one, one, A);
        albedo_bands(A, SIGMA_SPLIT, &base_ab);
        for (int wi = 0; wi < 3; wi++) {
            render_E(A, E[wi], wb, D, NULL, NULL);
            display_bands(D, SIGMA_SPLIT, &base_db);

            for (int layer = 0; layer < 2; layer++)
                for (int ci = 0; ci < 3; ci++) {
                    float ks[3] = { 1, 1, 1 }, km[3] = { 1, 1, 1 };
                    if (layer == 0) ks[ci] = 1.20f; else km[ci] = 1.20f;
                    apply_layers(sc, ks, km, A);
                    BandStats ab, db;
                    albedo_bands(A, SIGMA_SPLIT, &ab);
                    render_E(A, E[wi], wb, D, NULL, NULL);
                    display_bands(D, SIGMA_SPLIT, &db);
                    /* structure knob -> micro band, and vice versa, same channel */
                    fprintf(f, "%s,albedo,%s,%s,micro,%s,%.6f,%.6f,%.4f\n",
                            WEATHER[wi].name, layer == 0 ? "structure" : "micro", chname[ci],
                            chname[ci], base_ab.m_rel[ci], ab.m_rel[ci],
                            (ab.m_rel[ci] / base_ab.m_rel[ci] - 1.0f) * 100.0f);
                    fprintf(f, "%s,albedo,%s,%s,structure,%s,%.6f,%.6f,%.4f\n",
                            WEATHER[wi].name, layer == 0 ? "structure" : "micro", chname[ci],
                            chname[ci], base_ab.s_log[ci], ab.s_log[ci],
                            (ab.s_log[ci] / base_ab.s_log[ci] - 1.0f) * 100.0f);
                    fprintf(f, "%s,display,%s,%s,micro,%s,%.6f,%.6f,%.4f\n",
                            WEATHER[wi].name, layer == 0 ? "structure" : "micro", chname[ci],
                            chname[ci], base_db.m_rel[ci], db.m_rel[ci],
                            (db.m_rel[ci] / base_db.m_rel[ci] - 1.0f) * 100.0f);
                    fprintf(f, "%s,display,%s,%s,structure,%s,%.6f,%.6f,%.4f\n",
                            WEATHER[wi].name, layer == 0 ? "structure" : "micro", chname[ci],
                            chname[ci], base_db.s_log[ci], db.s_log[ci],
                            (db.s_log[ci] / base_db.s_log[ci] - 1.0f) * 100.0f);
                }
        }

        /* Interference sample: does the micro layer change how the structure
           knob responds?  C = D(both, ks_R=1.2) - D(both, ks_R=1)
                             - [D(struct only, ks_R=1.2) - D(struct only, ks_R=1)]
           For a separable model C is flat.  Amplified x10 for the eye; the
           RMS ratio to the structure response is written to the table. */
        float* A0b = calloc((size_t)n * 3, sizeof(float));
        float* A1b = calloc((size_t)n * 3, sizeof(float));
        float* A0s = calloc((size_t)n * 3, sizeof(float));
        float* A1s = calloc((size_t)n * 3, sizeof(float));
        float* D0b = calloc((size_t)n * 3, sizeof(float));
        float* D1b = calloc((size_t)n * 3, sizeof(float));
        float* D0s = calloc((size_t)n * 3, sizeof(float));
        float* D1s = calloc((size_t)n * 3, sizeof(float));
        float* cross = calloc((size_t)n * 3, sizeof(float));
        float ks_p[3] = { 1.20f, 1, 1 };
        apply_layers(sc, one, one, A0b);
        apply_layers(sc, ks_p, one, A1b);
        apply_layers(sc, one, zero, A0s);
        apply_layers(sc, ks_p, zero, A1s);
        render_E(A0b, E[0], wb, D0b, NULL, NULL);
        render_E(A1b, E[0], wb, D1b, NULL, NULL);
        render_E(A0s, E[0], wb, D0s, NULL, NULL);
        render_E(A1s, E[0], wb, D1s, NULL, NULL);
        for (int c = 0; c < 3; c++) {
            double ec = 0, er = 0;
            for (int i = 0; i < n; i++) {
                double ct = (D1b[i * 3 + c] - D0b[i * 3 + c])
                          - (D1s[i * 3 + c] - D0s[i * 3 + c]);
                double resp = D1s[i * 3 + c] - D0s[i * 3 + c];
                ec += ct * ct; er += resp * resp;
                cross[i * 3 + c] = clampf(0.5f + 30.0f * (float)ct, 0.0f, 1.0f);
            }
            fprintf(f, "clear,display,structure_R,%s,micro_presence,%s,%.6f,%.6f,%.4f\n",
                    chname[0], chname[c], 0.0, sqrt(ec / n), sqrt(ec / er) * 100.0);
        }
        write_ppm("out/interference_sample.ppm", cross, N, N);
        free(A0b); free(A1b); free(A0s); free(A1s);
        free(D0b); free(D1b); free(D0s); free(D1s); free(cross);
        free(A); free(D);
    }
    fclose(f);

    /* --- 9. stability summary per weather.
       The question "which weather is the split most stable under" is a
       *display-space* question: the tone curve's shoulder decides how much of
       the micro band survives, and it compresses differently in each weather.
       Metrics:
         bias_{r,g,b}_linear_pct : the weather's own colour cast, physical
         bias_{r,g,b}_display_pct: what is left after the fixed tone curve
         layered_crosstalk_pct   : display micro-band drift when the structure
                                   knob is pushed +20 % (should be small)
         merged_crosstalk_pct    : display micro-band drift when the single
                                   merged gain is pushed +20 % (the status quo)
         display_micro_spread_pct: max/min of the per-channel display micro
                                   amplitude - how evenly the layer survives
       --- */
    f = fopen("out/stability.csv", "w");
    fprintf(f, "weather,bias_r_linear_pct,bias_g_linear_pct,bias_b_linear_pct,max_abs_lin_pct,"
               "bias_r_display_pct,bias_g_display_pct,bias_b_display_pct,max_abs_display_pct,"
               "layered_bias_drift_pct,merged_bias_drift_pct,"
               "display_micro_spread_pct,display_micro_r,display_micro_g,display_micro_b,"
               "albedo_micro_chroma_spread_pct,clip_frac,stability_index_pct\n");
    for (int wi = 0; wi < 3; wi++) {
        float maxlin = 0, maxdisp = 0;
        for (int c = 0; c < 3; c++) {
            maxlin = fmaxf(maxlin, fabsf(base_lin[wi][c]));
            maxdisp = fmaxf(maxdisp, fabsf(base_dev[wi][c]));
        }
        float* A = calloc((size_t)n * 3, sizeof(float));
        float* Dd = calloc((size_t)n * 3, sizeof(float));

        /* the user's complaint, measured directly: push ONE channel by
           +20 %, look at the *other* channels' display bias.  "Pressing one
           channel down makes another pop up" is exactly this number.  Both
           schemes are pushed per-channel, so the comparison is like-for-like.
           (A separate band-amplitude crosstalk number is in cross_terms.csv.) */
        float dev0[3], dev1[3];
        float lcross = 0, mcross = 0;
        apply_layers(sc, one, one, A);
        render_E(A, E[wi], wb, Dd, NULL, NULL);
        channel_dev(Dd, dev0);
        for (int ci = 0; ci < 3; ci++) {
            float ks1[3] = { 1, 1, 1 };
            ks1[ci] = 1.20f;
            apply_layers(sc, ks1, one, A);
            render_E(A, E[wi], wb, Dd, NULL, NULL);
            channel_dev(Dd, dev1);
            for (int c = 0; c < 3; c++)
                if (c != ci)
                    lcross = fmaxf(lcross, fabsf(dev1[c] - dev0[c]));
        }
        apply_merged(sc, one, A);
        render_E(A, E[wi], wb, Dd, NULL, NULL);
        channel_dev(Dd, dev0);
        for (int ci = 0; ci < 3; ci++) {
            float kg1[3] = { 1, 1, 1 };
            kg1[ci] = 1.20f;
            apply_merged(sc, kg1, A);
            render_E(A, E[wi], wb, Dd, NULL, NULL);
            channel_dev(Dd, dev1);
            for (int c = 0; c < 3; c++)
                if (c != ci)
                    mcross = fmaxf(mcross, fabsf(dev1[c] - dev0[c]));
        }
        free(A); free(Dd);

        float dl = fmaxf(fminf(base_db[wi].m_rel[0], fminf(base_db[wi].m_rel[1], base_db[wi].m_rel[2])), 1e-4f);
        float dh = fmaxf(base_db[wi].m_rel[0], fmaxf(base_db[wi].m_rel[1], base_db[wi].m_rel[2]));
        float al = fminf(base_ab[wi].m_rel[0], fminf(base_ab[wi].m_rel[1], base_ab[wi].m_rel[2]));
        float ah = fmaxf(base_ab[wi].m_rel[0], fmaxf(base_ab[wi].m_rel[1], base_ab[wi].m_rel[2]));
        float clipf = 0;
        apply_layers(sc, one, one, Av);
        render_E(Av, E[wi], wb, D, NULL, &clipf);
        fprintf(f, "%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.2f,%.4f,%.4f,%.4f,%.2f,%.5f,%.3f\n",
                WEATHER[wi].name,
                base_lin[wi][0], base_lin[wi][1], base_lin[wi][2], maxlin,
                base_dev[wi][0], base_dev[wi][1], base_dev[wi][2], maxdisp,
                lcross, mcross,
                (dh / dl - 1.0f) * 100.0f,
                base_db[wi].m_rel[0], base_db[wi].m_rel[1], base_db[wi].m_rel[2],
                (ah / al - 1.0f) * 100.0f, clipf, maxlin + lcross);
    }
    fclose(f);

    /* --- 10. the "press one, another pops" measurement.
       With a merged gain, zeroing channel c's bias is a *brightness* change:
       it moves every other channel's bias too.  With the layered scheme the
       equivalent correction is a contrast change on one band, which leaves
       all channel means (and therefore every other channel's bias) alone. --- */
    f = fopen("out/correction.csv", "w");
    fprintf(f, "weather,scheme,knob_channel,knob_gain,own_bias_after_pct,other_channel,"
               "other_bias_after_pct,other_bias_drift_pct,own_structure_change_pct,"
               "own_micro_change_pct\n");
    for (int wi = 0; wi < 3; wi++) {
        float base_lin[3];
        apply_layers(sc, one, one, Av);
        linear_dev(Av, E[wi], wb, base_lin);
        for (int ci = 0; ci < 3; ci++) {
            /* merged: bisect the gain that zeroes this channel's linear bias */
            float lo = 0.30f, hi = 3.0f, kg[3] = { 1, 1, 1 };
            float dev[3];
            for (int it = 0; it < 44; it++) {
                float mid = 0.5f * (lo + hi);
                kg[ci] = mid;
                apply_merged(sc, kg, Av);
                linear_dev(Av, E[wi], wb, dev);
                if (dev[ci] > 0) hi = mid; else lo = mid;
            }
            float gfix = 0.5f * (lo + hi);
            kg[ci] = gfix;
            apply_merged(sc, kg, Av);
            linear_dev(Av, E[wi], wb, dev);
            for (int c = 0; c < 3; c++)
                if (c != ci)
                    fprintf(f, "%s,merged,%s,%.4f,%.4f,%s,%.4f,%.4f,%.2f,%.2f\n",
                            WEATHER[wi].name, chname[ci], gfix, dev[ci], chname[c], dev[c],
                            dev[c] - base_lin[c],
                            (gfix - 1.0f) * 100.0f, (gfix - 1.0f) * 100.0f);

            /* layered: +20 % contrast on the same channel.  Its own band moves;
               the channel means - and so every other channel's bias - stay. */
            float ksl[3] = { 1, 1, 1 };
            BandStats b0, b1;
            apply_layers(sc, one, one, Av);
            albedo_bands(Av, SIGMA_SPLIT, &b0);
            ksl[ci] = 1.20f;
            apply_layers(sc, ksl, one, Av);
            linear_dev(Av, E[wi], wb, dev);
            albedo_bands(Av, SIGMA_SPLIT, &b1);
            for (int c = 0; c < 3; c++)
                if (c != ci)
                    fprintf(f, "%s,layered,%s,1.2000,%.4f,%s,%.4f,%.4f,%.2f,%.2f\n",
                            WEATHER[wi].name, chname[ci], dev[ci], chname[c], dev[c],
                            dev[c] - base_lin[c],
                            (b1.s_log[ci] / b0.s_log[ci] - 1.0f) * 100.0f,
                            (b1.m_rel[ci] / b0.m_rel[ci] - 1.0f) * 100.0f);
        }
    }
    fclose(f);

    for (int wi = 0; wi < 3; wi++) free(E[wi]);
    free(Av); free(D);
}

/* ---------------------------------------------------------------- main */
int main(int argc, char** argv) {
    const char* mode = argc > 1 ? argv[1] : "test";
    Scene sc;
    scene_build(&sc, 20261124u);
    if (strcmp(mode, "test") == 0) {
        test_scene(&sc);
        test_knobs(&sc);
        test_weather(&sc);
        printf("all_ok=%d\n", failures == 0);
        return failures == 0 ? 0 : 1;
    }
    if (strcmp(mode, "build") == 0) {
        build_outputs(&sc);
        printf("wrote out/ CSV reports and PPM images\n");
        return 0;
    }
    fprintf(stderr, "usage: %s [test|build]\n", argv[0]);
    return 2;
}
