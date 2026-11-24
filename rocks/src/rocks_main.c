/* rocks: seeded gravel-texture deviation lab.
 *
 * 1. Build a fixed pebble layout and tile it into a 12x8 chart.
 * 2. Render the chart under N texture seeds with two independent texture
 *    modules (value fBm vs gradient fBm).
 * 3. At every pixel take min/max luminance across seeds: the span is the
 *    "same coordinate, different seed" deviation. Fold it into a local
 *    deviation heatmap and flag sites whose span is absurd (> threshold).
 * 4. Inject a small Gaussian height bump, drape the chart over it, measure
 *    texture displacement + stretch, and crosshair-mark sites where the
 *    deformation exceeds a given scale.
 * 5. Re-run the identical terrain slice through both texture modules and
 *    verify bit-for-bit reproducibility of the measurements.
 *
 * Usage: rocks [--out DIR] [--seeds N] [--span-limit X] [--deform-scale X]
 *              [--drape-k X] [--bump CX,CY,SIGMA,AMP] [--modules A|B|AB]
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "analysis.h"
#include "gfx.h"
#include "stones.h"

#define PANEL_W 500
#define PANEL_CHART_W 480
#define PANEL_CHART_H 320
#define GAP 24
#define PAD 32
#define BOARD_W (PAD * 2 + PANEL_W * 3 + GAP * 2)
#define BOARD_H 1412

typedef struct {
    int module;
    Deviation dev;
    Drape drape;
    HeightSlice slice;
    uint32_t hash_flat;
    uint32_t hash_dev;
    uint32_t hash_drape;
    uint32_t hash_slice;
} Result;

/* mkdir -p, so --out can be any path */
static int failures = 0;
static void check(const char *name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static int mkdir_p(const char *path) {
    char buf[1024];
    size_t n = strlen(path);
    if (n == 0 || n >= sizeof(buf)) return -1;
    memcpy(buf, path, n + 1);
    for (size_t i = 1; i < n; i++) {
        if (buf[i] == '/') {
            buf[i] = '\0';
            mkdir(buf, 0755);
            buf[i] = '/';
        }
    }
    return mkdir(buf, 0755);
}

static void usage(void) {
    fprintf(stderr,
            "usage: rocks [--out DIR] [--seeds N] [--span-limit X] [--deform-scale X]\n"
            "             [--drape-k X] [--bump CX,CY,SIGMA,AMP] [--modules A|B|AB] [--check]\n");
}

static void upscale_chart(Canvas *dst, int x0, int y0, const float *lum, int w, int h, int scale) {
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            float v = lum[y * w + x];
            int g = (int)(v * 255.0f + 0.5f);
            if (g < 0) g = 0;
            if (g > 255) g = 255;
            fill_rect(dst, x0 + x * scale, y0 + y * scale, scale, scale, RGB(g, g, g));
        }
    }
}

/* Transfer from data value to ramp position.
 *
 * Two things are going on: the interesting values cluster well below the
 * peak, so a linear ramp spends most of its ink on the top few pixels; and
 * the far field is numerically ~1e-4 (about 0.1% of peak), which any
 * low-end lift turns into visible haze. So: a 1.5% dead zone (values below
 * it are exactly zero - that is the measurement noise floor, annotated on
 * the panel) followed by a gamma of 1.5 on the remainder. The legend strip
 * is drawn with this same function, so the ramp on the panel is the ramp in
 * the key. */
static float heat_t(float value, float domain) {
    const float dead = 0.015f;
    float t = domain > 0.0f ? value / domain : 0.0f;
    if (!(t > dead)) return 0.0f;
    if (t > 1.0f) t = 1.0f;
    return powf((t - dead) / (1.0f - dead), 0.66f);
}

static void upscale_heat(Canvas *dst, int x0, int y0, const float *val, int w, int h,
                         float domain, int scale) {
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            fill_rect(dst, x0 + x * scale, y0 + y * scale, scale, scale,
                      blue_ramp(heat_t(val[y * w + x], domain)));
        }
    }
}

static void panel_frame(Canvas *c, int x, int y, int w, int h, const char *title,
                        const char *sub) {
    fill_rect(c, x, y, w, h, COL_SURFACE);
    stroke_rect(c, x, y, w, h, COL_BORDER);
    draw_text(c, x + 16, y + 14, title, 2, COL_INK);
    if (sub) draw_text(c, x + 16, y + 38, sub, 1, COL_MUTED);
}

static void legend_bar(Canvas *c, int x, int y, int w, float lo, float hi, const char *lo_lab,
                       const char *hi_lab) {
    /* the strip is sampled in ramp position, so it shows the real transfer */
    for (int i = 0; i < w; i++) {
        float f = (float)i / (float)(w - 1);
        fill_rect(c, x + i, y, 1, 10, blue_ramp(f));
    }
    stroke_rect(c, x, y, w, 10, COL_BORDER);
    draw_text(c, x, y + 16, lo_lab, 1, COL_MUTED);
    int hw = text_width(hi_lab, 1);
    draw_text(c, x + w - hw, y + 16, hi_lab, 1, COL_MUTED);
    char mid[40];
    snprintf(mid, sizeof(mid), "%.2f", 0.5f * (lo + hi));
    int mw = text_width(mid, 1);
    draw_text(c, x + w / 2 - mw / 2, y + 16, mid, 1, COL_MUTED);
}

/* crosshair with a white halo so it survives dark and light texture */
static void marked_crosshair(Canvas *c, int cx, int cy, int arm) {
    draw_line(c, cx - arm - 1, cy, cx + arm + 1, cy, COL_WHITE);
    draw_line(c, cx, cy - arm - 1, cx, cy + arm + 1, COL_WHITE);
    draw_circle(c, cx, cy, 4, COL_WHITE);
    draw_line(c, cx - arm, cy, cx + arm, cy, COL_CRITICAL);
    draw_line(c, cx, cy - arm, cx, cy + arm, COL_CRITICAL);
    draw_circle(c, cx, cy, 3, COL_CRITICAL);
}

static void crop_chart(Canvas *dst, int x0, int y0, const float *lum, int w, int cx0, int cy0,
                       int cw, int ch, int scale) {
    for (int y = 0; y < ch; y++) {
        for (int x = 0; x < cw; x++) {
            float v = lum[(cy0 + y) * w + (cx0 + x)];
            int g = (int)(v * 255.0f + 0.5f);
            if (g < 0) g = 0;
            if (g > 255) g = 255;
            fill_rect(dst, x0 + x * scale, y0 + y * scale, scale, scale, RGB(g, g, g));
        }
    }
}

static void render_individual(const Config *cfg, const StoneSet *set, const Result *r,
                              const char *outdir) {
    int w = cfg->img_w, h = cfg->img_h, scale = 2;
    char path[1024];
    float *lum = (float *)malloc(sizeof(float) * (size_t)w * h);
    if (!lum) return;

    Canvas c;
    canvas_init(&c, w * scale, h * scale);

    render_lum_flat(cfg, set, r->module, cfg->display_seed, lum);
    upscale_chart(&c, 0, 0, lum, w, h, scale);
    snprintf(path, sizeof(path), "%s/texture_%c.png", outdir, r->module == 0 ? 'A' : 'B');
    canvas_write_png(&c, path);

    upscale_heat(&c, 0, 0, r->dev.span_smooth, w, h, r->dev.span_domain, scale);
    for (int i = 0; i < cfg->cols * cfg->rows; i++) {
        if (!r->dev.site_absurd[i]) continue;
        int gx = i % ROCKS_MAX_COLS, gy = i / ROCKS_MAX_COLS;
        stroke_rect(&c, gx * cfg->cell_px * scale + 2, gy * cfg->cell_px * scale + 2,
                    cfg->cell_px * scale - 4, cfg->cell_px * scale - 4, COL_CRITICAL);
    }
    snprintf(path, sizeof(path), "%s/span_heatmap_%c.png", outdir, r->module == 0 ? 'A' : 'B');
    canvas_write_png(&c, path);

    render_lum_draped_grid(cfg, &r->slice, set, r->module, cfg->display_seed, lum);
    upscale_chart(&c, 0, 0, lum, w, h, scale);
    for (int i = 0; i < r->drape.cross_count; i++) {
        int s = r->drape.cross_site[i];
        int gx = s % ROCKS_MAX_COLS, gy = s / ROCKS_MAX_COLS;
        int cx = (gx * cfg->cell_px + cfg->cell_px / 2) * scale;
        int cy = (gy * cfg->cell_px + cfg->cell_px / 2) * scale;
        marked_crosshair(&c, cx, cy, 14);
    }
    snprintf(path, sizeof(path), "%s/draped_%c.png", outdir, r->module == 0 ? 'A' : 'B');
    canvas_write_png(&c, path);

    upscale_heat(&c, 0, 0, r->drape.deform, w, h,
                 r->drape.deform_max > 0.0f ? r->drape.deform_max : 1.0f, scale);
    for (int i = 0; i < r->drape.cross_count; i++) {
        int s = r->drape.cross_site[i];
        int gx = s % ROCKS_MAX_COLS, gy = s / ROCKS_MAX_COLS;
        int cx = (gx * cfg->cell_px + cfg->cell_px / 2) * scale;
        int cy = (gy * cfg->cell_px + cfg->cell_px / 2) * scale;
        marked_crosshair(&c, cx, cy, 14);
    }
    snprintf(path, sizeof(path), "%s/deform_heatmap_%c.png", outdir,
             r->module == 0 ? 'A' : 'B');
    canvas_write_png(&c, path);

    free(lum);
    canvas_free(&c);
}

static void draw_metrics(Canvas *board, int x, int y, int w, const Config *cfg,
                         const Result *ra, const Result *rb, int bump_site,
                         const char *repro) {
    char line[256];
    int lh = 18;
    int yy = y;
    int colw = w / 2;

    draw_text(board, x, yy, "MEASUREMENT SUMMARY", 2, COL_INK);
    yy += 26;

    /* two metric columns: deviation on the left, deformation on the right */
    for (int m = 0; m < 2; m++) {
        const Result *r = m == 0 ? ra : rb;
        const char *tag = m == 0 ? "MODULE A" : "MODULE B";
        int cx = x + m * colw;
        int cy = yy;
        snprintf(line, sizeof(line), "%s  SPAN", tag);
        draw_text(board, cx, cy, line, 1, COL_INK2);
        cy += lh;
        snprintf(line, sizeof(line), "  MEAN %.3f   STD %.3f", r->dev.site_mean, r->dev.site_std);
        draw_text(board, cx, cy, line, 1, COL_MUTED);
        cy += lh;
        snprintf(line, sizeof(line), "  LIMIT %.3f (AUTO MEAN+2SD)   OVER %d/%d", 
                 r->dev.absurd_threshold, r->dev.absurd_count, cfg->cols * cfg->rows);
        draw_text(board, cx, cy, line, 1, COL_MUTED);
        cy += lh;
        snprintf(line, sizeof(line), "  PEAK SPAN %.3f   STRETCH MAX %.3f   MARKED %d",
                 r->dev.span_max, r->drape.deform_max, r->drape.cross_count);
        draw_text(board, cx, cy, line, 1, COL_MUTED);
        cy += lh;
        snprintf(line, sizeof(line), "  MAX DISPLACEMENT %.3f WORLD UNITS", r->drape.disp_max);
        draw_text(board, cx, cy, line, 1, COL_MUTED);
    }
    yy += lh * 5 + 12;

    draw_text(board, x, yy, "REPRODUCIBILITY AND CONTROL", 2, COL_INK);
    yy += 26;
    snprintf(line, sizeof(line),
             "  SEEDS %d   DISPLAY SEED %u   LAYOUT SEED %u   TERRAIN BUMP SITE %d",
             cfg->n_seeds, cfg->display_seed, cfg->layout_seed, bump_site);
    draw_text(board, x, yy, line, 1, COL_INK2);
    yy += lh;
    snprintf(line, sizeof(line), "  HASH A  FLAT %08X  DEV %08X  DRAPE %08X", ra->hash_flat,
             ra->hash_dev, ra->hash_drape);
    draw_text(board, x, yy, line, 1, COL_MUTED);
    yy += lh;
    snprintf(line, sizeof(line), "  HASH B  FLAT %08X  DEV %08X  DRAPE %08X", rb->hash_flat,
             rb->hash_dev, rb->hash_drape);
    draw_text(board, x, yy, line, 1, COL_MUTED);
    yy += lh + 4;

    char buf[1024];
    snprintf(buf, sizeof(buf), "%s", repro);
    char *save = NULL;
    for (char *tok = strtok_r(buf, "\n", &save); tok; tok = strtok_r(NULL, "\n", &save)) {
        draw_text(board, x, yy, tok, 1, COL_INK2);
        yy += 16;
    }
    (void)w;
}

static void draw_crosshairs(Canvas *board, const Config *cfg, const Drape *d, int ox, int oy) {
    for (int i = 0; i < d->cross_count; i++) {
        int s = d->cross_site[i];
        int gx = s % ROCKS_MAX_COLS, gy = s / ROCKS_MAX_COLS;
        int cx = ox + (gx * cfg->cell_px + cfg->cell_px / 2) * 2;
        int cy = oy + (gy * cfg->cell_px + cfg->cell_px / 2) * 2;
        marked_crosshair(board, cx, cy, 13);
    }
}

static void render_board(const Config *cfg, const StoneSet *set, const Result *ra,
                         const Result *rb, const float *lum_a, const float *lum_b,
                         int bump_site, const char *repro, const char *path) {
    Canvas board;
    if (canvas_init(&board, BOARD_W, BOARD_H) != 0) return;
    canvas_clear(&board, COL_PAGE);

    draw_text(&board, PAD, 28, "SEEDED GRAVEL TEXTURE - LOCAL DEVIATION LAB", 3, COL_INK);
    draw_text(&board, PAD, 64, "FIXED STONE LAYOUT / DIFFERENT TEXTURE SEEDS / SMALL TERRAIN BUMP",
              1, COL_INK2);

    int row1 = 100;
    int panel_h = 430;
    int row2 = row1 + panel_h + 24;

    /* panel 1: texture A */
    int x = PAD;
    panel_frame(&board, x, row1, PANEL_W, panel_h, "1  TEXTURE MODULE A", "VALUE-NOISE FBM");
    upscale_chart(&board, x + 10, row1 + 58, lum_a, cfg->img_w, cfg->img_h, 2);
    char sub[128];
    snprintf(sub, sizeof(sub), "SEED %u - FLAT / FULL CHART", cfg->display_seed);
    draw_text(&board, x + 16, row1 + panel_h - 24, sub, 1, COL_MUTED);

    /* panel 2: texture B */
    x += PANEL_W + GAP;
    panel_frame(&board, x, row1, PANEL_W, panel_h, "2  TEXTURE MODULE B", "GRADIENT-NOISE FBM");
    upscale_chart(&board, x + 10, row1 + 58, lum_b, cfg->img_w, cfg->img_h, 2);
    snprintf(sub, sizeof(sub), "SEED %u - FLAT / FULL CHART", cfg->display_seed);
    draw_text(&board, x + 16, row1 + panel_h - 24, sub, 1, COL_MUTED);

    /* panel 3: same site crop under 4 runs */
    x += PANEL_W + GAP;
    panel_frame(&board, x, row1, PANEL_W, panel_h, "3  SAME SITE, 4 RUNS",
                "FIXED LAYOUT - GRAIN ONLY");
    {
        const int crop_cells = 3;
        int cw = crop_cells * cfg->cell_px;
        /* centre the crop on the site that diverges most across seeds, so the
         * thumbnails actually show the phenomenon rather than a calm patch */
        int hot = 0;
        for (int i = 1; i < cfg->cols * cfg->rows; i++)
            if (ra->dev.site_span[i] > ra->dev.site_span[hot]) hot = i;
        int hgx = hot % cfg->cols, hgy = hot / cfg->cols;
        int cx0 = (hgx - 1) * cfg->cell_px;
        int cy0 = (hgy - 1) * cfg->cell_px;
        if (cx0 < 0) cx0 = 0;
        if (cy0 < 0) cy0 = 0;
        if (cx0 > cfg->img_w - cw) cx0 = cfg->img_w - cw;
        if (cy0 > cfg->img_h - cw) cy0 = cfg->img_h - cw;
        int cscale = 2;
        int thumb = cw * cscale;
        int tx = x + 16, ty = row1 + 60;
        unsigned seeds[4] = {cfg->seeds[0], cfg->seeds[1], cfg->seeds[0], cfg->seeds[1]};
        const char *labs[4] = {"A SEED 1", "A SEED 2", "B SEED 1", "B SEED 2"};
        float *buf = (float *)malloc(sizeof(float) * (size_t)cfg->img_w * cfg->img_h);
        for (int k = 0; k < 4 && buf; k++) {
            int mod = k < 2 ? 0 : 1;
            render_lum_flat(cfg, set, mod, seeds[k], buf);
            int px = tx + (k % 2) * (thumb + 10);
            int py = ty + (k / 2) * (thumb + 26);
            crop_chart(&board, px, py, buf, cfg->img_w, cx0, cy0, cw, cw, cscale);
            stroke_rect(&board, px - 1, py - 1, thumb + 2, thumb + 2, COL_BORDER);
            draw_text(&board, px, py + thumb + 5, labs[k], 1, COL_INK2);
        }
        free(buf);
        char cl[160];
        snprintf(cl, sizeof(cl), "WORST SITE %d,%d (SPAN %.2f) - SAME CROP WINDOW", hgx, hgy,
                 ra->dev.site_span[hot]);
        draw_text(&board, x + 16, row1 + panel_h - 24, cl, 1, COL_MUTED);
    }

    /* row 2 panel 4: span heatmap A */
    x = PAD;
    panel_frame(&board, x, row2, PANEL_W, panel_h, "4  LOCAL DEVIATION",
                "MAX-MIN LUMINANCE ACROSS SEEDS");
    upscale_heat(&board, x + 10, row2 + 58, ra->dev.span_smooth, cfg->img_w, cfg->img_h,
                 ra->dev.span_domain, 2);
    for (int i = 0; i < cfg->cols * cfg->rows; i++) {
        if (!ra->dev.site_absurd[i]) continue;
        int gx = i % ROCKS_MAX_COLS, gy = i / ROCKS_MAX_COLS;
        int px = x + 10 + gx * cfg->cell_px * 2;
        int py = row2 + 58 + gy * cfg->cell_px * 2;
        stroke_rect(&board, px + 2, py + 2, cfg->cell_px * 2 - 4, cfg->cell_px * 2 - 4,
                    COL_CRITICAL);
        fill_rect(&board, px + cfg->cell_px - 2, py + cfg->cell_px - 2, 5, 5, COL_CRITICAL);
    }
    char hi_lab[32];
    snprintf(hi_lab, sizeof(hi_lab), "%.2f", ra->dev.span_domain);
    legend_bar(&board, x + 16, row2 + 386, 190, 0.0f, ra->dev.span_domain, "0.00", hi_lab);
    {
        char note[160];
        snprintf(note, sizeof(note), "RED BOX  SPAN > %.2f   %d SITES   (NONLINEAR KEY)",
                 ra->dev.absurd_threshold, ra->dev.absurd_count);
        draw_text(&board, x + 16, row2 + 418, note, 1, COL_INK2);
    }

    /* panel 5: draped texture with crosshairs */
    x += PANEL_W + GAP;
    panel_frame(&board, x, row2, PANEL_W, panel_h, "5  TEXTURE OVER TERRAIN BUMP",
                "CROSSHAIR - LOCAL STRETCH > SCALE / GRID = 1 WORLD UNIT");
    {
        float *lum = (float *)malloc(sizeof(float) * (size_t)cfg->img_w * cfg->img_h);
        if (lum) {
            render_lum_draped_grid(cfg, &ra->slice, set, 0, cfg->display_seed, lum);
            upscale_chart(&board, x + 10, row2 + 58, lum, cfg->img_w, cfg->img_h, 2);
            free(lum);
        }
    }
    draw_crosshairs(&board, cfg, &ra->drape, x + 10, row2 + 58);
    {
        char note[160];
        snprintf(note, sizeof(note), "BUMP CENTRE %.1f,%.1f  SIGMA %.2f  AMP %.2f  K %.2f",
                 cfg->bump.cx, cfg->bump.cy, cfg->bump.sigma, cfg->bump.amp, cfg->drape_k);
        draw_text(&board, x + 16, row2 + 386, note, 1, COL_MUTED);
        draw_text(&board, x + 16, row2 + 418,
                  "GRID SPACING WIDENS ON SLOPES - TEXTURE IS STRETCHED", 1, COL_MUTED);
    }

    /* panel 6: deformation heatmap */
    x += PANEL_W + GAP;
    panel_frame(&board, x, row2, PANEL_W, panel_h, "6  DEFORMATION",
                "LOCAL STRETCH OF THE DRAPE FIELD");
    {
        float dom = ra->drape.deform_max > 0.0f ? ra->drape.deform_max : 1.0f;
        upscale_heat(&board, x + 10, row2 + 58, ra->drape.deform, cfg->img_w, cfg->img_h, dom, 2);
        draw_crosshairs(&board, cfg, &ra->drape, x + 10, row2 + 58);
        char lab_hi[32];
        snprintf(lab_hi, sizeof(lab_hi), "%.2f", dom);
        legend_bar(&board, x + 16, row2 + 386, 190, 0.0f, dom, "0.00", lab_hi);
        char note[160];
        snprintf(note, sizeof(note), "SCALE %.2f   %d SITES MARKED   (NONLINEAR KEY)",
                 cfg->deform_scale, ra->drape.cross_count);
        draw_text(&board, x + 16, row2 + 418, note, 1, COL_INK2);
    }

    /* metric band */
    int band_y = row2 + panel_h + 16;
    int band_h = 372;
    fill_rect(&board, PAD, band_y, BOARD_W - 2 * PAD, band_h, COL_SURFACE);
    stroke_rect(&board, PAD, band_y, BOARD_W - 2 * PAD, band_h, COL_BORDER);
    draw_metrics(&board, PAD + 20, band_y + 18, BOARD_W - 2 * PAD - 40, cfg, ra, rb, bump_site,
                 repro);

    draw_text(&board, PAD, BOARD_H - 20,
              "MINI-PROC-LAB / ROCKS - ALL MEASUREMENTS DETERMINISTIC FOR A FIXED CONFIG", 1,
              COL_MUTED);

    canvas_write_png(&board, path);
    canvas_free(&board);
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

static int run_once(const Config *cfg, const StoneSet *set, int module, Result *r) {
    memset(r, 0, sizeof(*r));
    r->module = module;
    slice_build(cfg, &r->slice);
    r->hash_slice = hash_floats(r->slice.h, r->slice.n * r->slice.n);
    deviation_analyse(cfg, set, module, &r->dev);
    drape_analyse(cfg, &r->slice, &r->drape);

    int npix = cfg->img_w * cfg->img_h;
    float *lum = (float *)malloc(sizeof(float) * (size_t)npix);
    if (!lum) return -1;
    render_lum_flat(cfg, set, module, cfg->display_seed, lum);
    r->hash_flat = hash_floats(lum, npix);
    render_lum_draped(cfg, &r->slice, set, module, cfg->display_seed, lum);
    r->hash_drape = hash_floats(lum, npix);
    free(lum);

    /* dev hash covers span field + site stats + absurd mask */
    r->hash_dev = hash_floats(r->dev.span_px, npix);
    r->hash_dev ^= hash_floats(r->dev.site_span, cfg->cols * cfg->rows) * 2654435761u;
    uint32_t mask_hash = 2166136261u;
    for (int i = 0; i < cfg->cols * cfg->rows; i++) {
        mask_hash ^= r->dev.site_absurd[i];
        mask_hash *= 16777619u;
    }
    r->hash_dev ^= mask_hash * 40503u;
    return 0;
}

int main(int argc, char **argv) {
    Config cfg;
    config_default(&cfg);
    const char *outdir = "artifacts";
    const char *modules = "AB";
    int check_mode = 0;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--out") && i + 1 < argc) {
            outdir = argv[++i];
        } else if (!strcmp(argv[i], "--seeds") && i + 1 < argc) {
            config_set_seed_count(&cfg, atoi(argv[++i]));
        } else if (!strcmp(argv[i], "--span-limit") && i + 1 < argc) {
            cfg.span_limit = (float)atof(argv[++i]);
        } else if (!strcmp(argv[i], "--deform-scale") && i + 1 < argc) {
            cfg.deform_scale = (float)atof(argv[++i]);
        } else if (!strcmp(argv[i], "--drape-k") && i + 1 < argc) {
            cfg.drape_k = (float)atof(argv[++i]);
        } else if (!strcmp(argv[i], "--bump") && i + 1 < argc) {
            float cx, cy, sg, am;
            if (sscanf(argv[++i], "%f,%f,%f,%f", &cx, &cy, &sg, &am) == 4) {
                cfg.bump.cx = cx; cfg.bump.cy = cy; cfg.bump.sigma = sg; cfg.bump.amp = am;
                cfg.bump_auto = 0; /* an explicit placement wins over auto */
            } else {
                fprintf(stderr, "bad --bump value, expected CX,CY,SIGMA,AMP\n");
                return 2;
            }
        } else if (!strcmp(argv[i], "--modules") && i + 1 < argc) {
            modules = argv[++i];
            if (!strchr(modules, 'A') && !strchr(modules, 'B')) {
                fprintf(stderr, "--modules expects A, B or AB\n");
                return 2;
            }
        } else if (!strcmp(argv[i], "--check")) {
            check_mode = 1;
        } else if (!strcmp(argv[i], "--help")) {
            usage();
            return 0;
        } else {
            usage();
            return 2;
        }
    }

    StoneSet set;
    stones_layout(&set, cfg.layout_seed);

    int bump_site = -1;
    if (cfg.bump_auto) {
        Deviation probe;
        deviation_analyse(&cfg, &set, 0, &probe);
        bump_site = config_place_bump_best_site(&cfg, &probe);
        printf("bump_auto site=%d centre=%.2f,%.2f\n", bump_site, cfg.bump.cx, cfg.bump.cy);
        deviation_free(&probe);
    }

    Result ra, rb;
    if (run_once(&cfg, &set, 0, &ra) != 0) return 1;
    if (run_once(&cfg, &set, 1, &rb) != 0) return 1;

    /* --- independent rerun for the reproducibility check --- */
    StoneSet set2;
    stones_layout(&set2, cfg.layout_seed);
    Result ra2, rb2;
    run_once(&cfg, &set2, 0, &ra2);
    run_once(&cfg, &set2, 1, &rb2);

    int layout_same = memcmp(&set, &set2, sizeof(StoneSet)) == 0;
    /* the two modules must have consumed the *same* discrete slice */
    int slice_same = memcmp(&ra.slice.h, &rb.slice.h, sizeof(ra.slice.h)) == 0 &&
                     ra.hash_slice == rb.hash_slice && ra2.hash_slice == ra.hash_slice;
    uint32_t slice_hash = ra.hash_slice;

    int run_a_same = ra.hash_flat == ra2.hash_flat && ra.hash_dev == ra2.hash_dev &&
                     ra.hash_drape == ra2.hash_drape;
    int run_b_same = rb.hash_flat == rb2.hash_flat && rb.hash_dev == rb2.hash_dev &&
                     rb.hash_drape == rb2.hash_drape;

    /* cross-module: geometry-only measurements must agree exactly */
    int cross_same = ra.drape.cross_count == rb.drape.cross_count;
    if (cross_same) {
        for (int i = 0; i < ra.drape.cross_count; i++) {
            if (ra.drape.cross_site[i] != rb.drape.cross_site[i]) { cross_same = 0; break; }
        }
    }
    int deform_same = ra.drape.deform_max > 0.0f
                          ? fabsf(ra.drape.deform_max - rb.drape.deform_max) < 1e-6f
                          : rb.drape.deform_max == 0.0f;

    float *flat_a = (float *)malloc(sizeof(float) * (size_t)cfg.img_w * cfg.img_h);
    float *flat_b = (float *)malloc(sizeof(float) * (size_t)cfg.img_w * cfg.img_h);
    float corr = 0.0f, seed_change = 0.0f;
    if (flat_a && flat_b) {
        render_lum_flat(&cfg, &set, 0, cfg.display_seed, flat_a);
        render_lum_flat(&cfg, &set, 1, cfg.display_seed, flat_b);
        corr = pearson_subsampled(flat_a, flat_b, cfg.img_w, cfg.img_h, 3);
        float *flat_a2 = (float *)malloc(sizeof(float) * (size_t)cfg.img_w * cfg.img_h);
        if (flat_a2) {
            render_lum_flat(&cfg, &set, 0, cfg.seeds[1], flat_a2);
            double acc = 0.0;
            for (int i = 0; i < cfg.img_w * cfg.img_h; i++)
                acc += fabs((double)flat_a[i] - flat_a2[i]);
            seed_change = (float)(acc / (cfg.img_w * cfg.img_h));
            free(flat_a2);
        }
    }

    /* ---------- machine-readable report ---------- */
    printf("config cols=%d rows=%d cell_px=%d seeds=%d layout_seed=%u\n", cfg.cols, cfg.rows,
           cfg.cell_px, cfg.n_seeds, cfg.layout_seed);
    printf("config deform_scale=%.4f drape_k=%.4f bump=%.3f,%.3f,%.3f,%.3f\n",
           cfg.deform_scale, cfg.drape_k, cfg.bump.cx, cfg.bump.cy, cfg.bump.sigma, cfg.bump.amp);
    printf("layout stones=%d identical_rerun=%d\n", set.n, layout_same);
    printf("span A mean=%.4f std=%.4f thresh=%.4f absurd=%d/%d max=%.4f\n", ra.dev.site_mean,
           ra.dev.site_std, ra.dev.absurd_threshold, ra.dev.absurd_count,
           cfg.cols * cfg.rows, ra.dev.span_max);
    printf("span B mean=%.4f std=%.4f thresh=%.4f absurd=%d/%d max=%.4f\n", rb.dev.site_mean,
           rb.dev.site_std, rb.dev.absurd_threshold, rb.dev.absurd_count,
           cfg.cols * cfg.rows, rb.dev.span_max);
    printf("deform scale=%.4f max_stretch=%.4f max_disp=%.4f crosshairs=%d\n", cfg.deform_scale,
           ra.drape.deform_max, ra.drape.disp_max, ra.drape.cross_count);
    printf("repro runA_same=%d runB_same=%d shared_slice_same=%d shared_slice_hash=%08x\n",
           run_a_same, run_b_same, slice_same, slice_hash);
    printf("repro crosshairs_A_B_identical=%d deform_max_A_B_equal=%d cross_module_lum_corr=%.4f\n",
           cross_same, deform_same, corr);
    printf("texture seed_change_mean_abs_lum=%.4f\n", seed_change);
    printf("hash A flat=%08x dev=%08x drape=%08x\n", ra.hash_flat, ra.hash_dev, ra.hash_drape);
    printf("hash B flat=%08x dev=%08x drape=%08x\n", rb.hash_flat, rb.hash_dev, rb.hash_drape);

    /* ---------- self checks (same format as the sibling modules) ---------- */
    int nsites = cfg.cols * cfg.rows;
    int absurd_mask_ok = 1, absurd_n = 0;
    for (int i = 0; i < nsites; i++) {
        if (ra.dev.site_absurd[i]) absurd_n++;
        if (ra.dev.site_absurd[i] != (ra.dev.site_span[i] > ra.dev.absurd_threshold))
            absurd_mask_ok = 0;
    }
    int span_range_ok = ra.dev.span_max <= 1.0f + 1e-6f && ra.dev.span_max >= 0.0f;
    int cross_ok = 1;
    for (int i = 0; i < ra.drape.cross_count; i++) {
        int site = ra.drape.cross_site[i];
        if (site < 0 || site >= nsites) cross_ok = 0;
        else if (!(ra.drape.site_deform[site] > cfg.deform_scale)) cross_ok = 0;
    }
    int cross_cap_ok = ra.drape.cross_count <= ROCKS_MAX_MARKS;
    int mark_complete = 1;
    for (int i = 0; i < nsites; i++) {
        int over = ra.drape.site_deform[i] > cfg.deform_scale;
        int found = 0;
        for (int k = 0; k < ra.drape.cross_count; k++)
            if (ra.drape.cross_site[k] == i) { found = 1; break; }
        /* every over-threshold site is marked unless the display cap clipped */
        if (over != found && ra.drape.cross_count < ROCKS_MAX_MARKS) mark_complete = 0;
    }

    check("layout_deterministic", layout_same);
    check("shared_slice_identical", slice_same);
    check("moduleA_run_deterministic", run_a_same);
    check("moduleB_run_deterministic", run_b_same);
    check("crosshair_sets_agree", cross_same);
    check("deform_max_modules_agree", deform_same);
    check("absurd_mask_consistent", absurd_mask_ok && absurd_n == ra.dev.absurd_count);
    check("span_in_unit_range", span_range_ok);
    check("span_threshold_positive", ra.dev.absurd_threshold > 0.0f);
    check("crosshairs_over_scale", cross_ok && cross_cap_ok);
    check("crosshair_mark_complete", mark_complete);
    check("seed_changes_texture", seed_change > 0.0f);
    check("modules_not_identical", corr < 0.9999f && corr > 0.5f);

    /* ---------- artifacts ---------- */
    char rocks_dir[512];
    char path[1100];
    snprintf(rocks_dir, sizeof(rocks_dir), "%s/rocks", outdir);
    if (!check_mode) mkdir_p(rocks_dir);
    snprintf(path, sizeof(path), "%s/report.png", rocks_dir);

    char repro[1024];
    snprintf(repro, sizeof(repro),
             "LAYOUT MEMCPY EQUAL: %s\n"
             "SHARED HEIGHT SLICE 33X33 EQUAL: %s (HASH %08X)\n"
             "RUN A BIT-IDENTICAL: %s\n"
             "RUN B BIT-IDENTICAL: %s\n"
             "CROSSHAIR SET A == B: %s\n"
             "DEFORM MAX A == B: %s (%.4f)\n"
             "A/B LUMINANCE CORR: %.3f\n"
             "SEED-TO-SEED MEAN |DL|: %.3f",
             layout_same ? "YES" : "NO", slice_same ? "YES" : "NO", slice_hash,
             run_a_same ? "YES" : "NO", run_b_same ? "YES" : "NO", cross_same ? "YES" : "NO",
             deform_same ? "YES" : "NO", ra.drape.deform_max, corr, seed_change);

    float *lum_a = (float *)malloc(sizeof(float) * (size_t)cfg.img_w * cfg.img_h);
    float *lum_b = (float *)malloc(sizeof(float) * (size_t)cfg.img_w * cfg.img_h);
    if (!check_mode && lum_a && lum_b) {
        render_lum_flat(&cfg, &set, 0, cfg.display_seed, lum_a);
        render_lum_flat(&cfg, &set, 1, cfg.display_seed, lum_b);
        render_board(&cfg, &set, &ra, &rb, lum_a, lum_b, bump_site, repro, path);
        printf("artifact %s\n", path);
    }

    if (!check_mode && strchr(modules, 'A')) render_individual(&cfg, &set, &ra, rocks_dir);
    if (!check_mode && strchr(modules, 'B')) render_individual(&cfg, &set, &rb, rocks_dir);

    /* machine-readable per-site records, anchored on the same data as the art */
    snprintf(path, sizeof(path), "%s/sites.csv", rocks_dir);
    FILE *csv = check_mode ? NULL : fopen(path, "w");
    if (csv) {
        fprintf(csv, "site,gx,gy,span_A,span_B,absurd_A,crosshair_A,stretch_A\n");
        const Drape *dr = &ra.drape;
        for (int i = 0; i < cfg.cols * cfg.rows; i++) {
            int marked = 0;
            for (int k = 0; k < dr->cross_count; k++)
                if (dr->cross_site[k] == i) { marked = 1; break; }
            fprintf(csv, "%d,%d,%d,%.4f,%.4f,%d,%d,%.4f\n", i, i % cfg.cols, i / cfg.cols,
                    ra.dev.site_span[i], rb.dev.site_span[i], ra.dev.site_absurd[i], marked,
                    dr->site_deform[i]);
        }
        fclose(csv);
        printf("artifact %s\n", path);
    }

    printf("all_ok=%d\n", failures == 0 ? 1 : 0);
    free(lum_a);
    free(lum_b);
    free(flat_a);
    free(flat_b);
    deviation_free(&ra.dev); drape_free(&ra.drape);
    deviation_free(&rb.dev); drape_free(&rb.drape);
    deviation_free(&ra2.dev); drape_free(&ra2.drape);
    deviation_free(&rb2.dev); drape_free(&rb2.drape);
    (void)modules;
    return failures == 0 ? 0 : 1;
}
