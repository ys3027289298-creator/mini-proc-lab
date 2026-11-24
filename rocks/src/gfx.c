#include "gfx.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <png.h>

int canvas_init(Canvas *c, int w, int h) {
    c->w = w;
    c->h = h;
    c->px = (uint8_t *)calloc((size_t)w * (size_t)h * 3u, 1);
    return c->px ? 0 : -1;
}

void canvas_free(Canvas *c) {
    free(c->px);
    c->px = NULL;
    c->w = c->h = 0;
}

void canvas_clear(Canvas *c, Color color) {
    fill_rect(c, 0, 0, c->w, c->h, color);
}

void canvas_copy(Canvas *dst, const Canvas *src) {
    memcpy(dst->px, src->px, (size_t)src->w * (size_t)src->h * 3u);
}

void put_px(Canvas *c, int x, int y, Color color) {
    if (x < 0 || y < 0 || x >= c->w || y >= c->h) return;
    uint8_t *p = c->px + ((size_t)y * (size_t)c->w + (size_t)x) * 3u;
    p[0] = (uint8_t)((color >> 16) & 0xff);
    p[1] = (uint8_t)((color >> 8) & 0xff);
    p[2] = (uint8_t)(color & 0xff);
}

void fill_rect(Canvas *c, int x, int y, int w, int h, Color color) {
    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++) put_px(c, x + i, y + j, color);
}

void stroke_rect(Canvas *c, int x, int y, int w, int h, Color color) {
    draw_line(c, x, y, x + w - 1, y, color);
    draw_line(c, x, y + h - 1, x + w - 1, y + h - 1, color);
    draw_line(c, x, y, x, y + h - 1, color);
    draw_line(c, x + w - 1, y, x + w - 1, y + h - 1, color);
}

void draw_line(Canvas *c, int x0, int y0, int x1, int y1, Color color) {
    int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
    int dy = -abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
    int err = dx + dy;
    for (;;) {
        put_px(c, x0, y0, color);
        if (x0 == x1 && y0 == y1) break;
        int e2 = 2 * err;
        if (e2 >= dy) { err += dy; x0 += sx; }
        if (e2 <= dx) { err += dx; y0 += sy; }
    }
}

void draw_circle(Canvas *c, int cx, int cy, int radius, Color color) {
    int x = radius, y = 0, err = 1 - radius;
    while (x >= y) {
        put_px(c, cx + x, cy + y, color); put_px(c, cx + y, cy + x, color);
        put_px(c, cx - y, cy + x, color); put_px(c, cx - x, cy + y, color);
        put_px(c, cx - x, cy - y, color); put_px(c, cx - y, cy - x, color);
        put_px(c, cx + y, cy - x, color); put_px(c, cx + x, cy - y, color);
        y++;
        if (err < 0) {
            err += 2 * y + 1;
        } else {
            x--;
            err += 2 * (y - x) + 1;
        }
    }
}

void draw_crosshair(Canvas *c, int cx, int cy, int arm, Color color) {
    draw_line(c, cx - arm, cy, cx + arm, cy, color);
    draw_line(c, cx, cy - arm, cx, cy + arm, color);
    draw_circle(c, cx, cy, 3, color);
}

Color blue_ramp(float t) {
    static const Color stops[] = {
        0xfcfcfb, 0xcde2fb, 0xb7d3f6, 0x9ec5f4, 0x86b6ef, 0x6da7ec,
        0x5598e7, 0x3987e5, 0x2a78d6, 0x256abf, 0x1c5cab, 0x184f95,
        0x104281, 0x0d366b,
    };
    int n = (int)(sizeof(stops) / sizeof(stops[0])) - 1;
    if (!(t > 0.0f)) t = 0.0f; /* also catches NaN */
    if (t > 1.0f) t = 1.0f;
    float f = t * (float)n;
    int i = (int)f;
    if (i >= n) i = n - 1;
    float u = f - (float)i;
    /* channel bytes must be signed before differencing - the packed Color is
     * unsigned, so truncating a channel on the way down wraps without this. */
    float ar = (float)((stops[i] >> 16) & 0xffu), br = (float)((stops[i + 1] >> 16) & 0xffu);
    float ag = (float)((stops[i] >> 8) & 0xffu), bg = (float)((stops[i + 1] >> 8) & 0xffu);
    float ab = (float)(stops[i] & 0xffu), bb = (float)(stops[i + 1] & 0xffu);
    int r = (int)(ar + (br - ar) * u + 0.5f);
    int g = (int)(ag + (bg - ag) * u + 0.5f);
    int bl = (int)(ab + (bb - ab) * u + 0.5f);
    return RGB(r, g, bl);
}

/* ------------------------------------------------------------------ */
/* 5x7 font                                                            */
/* ------------------------------------------------------------------ */

typedef struct {
    char ch;
    unsigned char rows[7]; /* low 5 bits per row */
} Glyph;

static const Glyph glyphs[] = {
    {'A', {0x0e, 0x11, 0x11, 0x1f, 0x11, 0x11, 0x11}},
    {'B', {0x1e, 0x11, 0x11, 0x1e, 0x11, 0x11, 0x1e}},
    {'C', {0x0e, 0x11, 0x10, 0x10, 0x10, 0x11, 0x0e}},
    {'D', {0x1e, 0x11, 0x11, 0x11, 0x11, 0x11, 0x1e}},
    {'E', {0x1f, 0x10, 0x10, 0x1e, 0x10, 0x10, 0x1f}},
    {'F', {0x1f, 0x10, 0x10, 0x1e, 0x10, 0x10, 0x10}},
    {'G', {0x0e, 0x11, 0x10, 0x17, 0x11, 0x11, 0x0f}},
    {'H', {0x11, 0x11, 0x11, 0x1f, 0x11, 0x11, 0x11}},
    {'I', {0x1f, 0x04, 0x04, 0x04, 0x04, 0x04, 0x1f}},
    {'J', {0x07, 0x02, 0x02, 0x02, 0x12, 0x12, 0x0c}},
    {'K', {0x11, 0x12, 0x14, 0x18, 0x14, 0x12, 0x11}},
    {'L', {0x10, 0x10, 0x10, 0x10, 0x10, 0x10, 0x1f}},
    {'M', {0x11, 0x1b, 0x15, 0x15, 0x11, 0x11, 0x11}},
    {'N', {0x11, 0x19, 0x15, 0x13, 0x11, 0x11, 0x11}},
    {'O', {0x0e, 0x11, 0x11, 0x11, 0x11, 0x11, 0x0e}},
    {'P', {0x1e, 0x11, 0x11, 0x1e, 0x10, 0x10, 0x10}},
    {'Q', {0x0e, 0x11, 0x11, 0x11, 0x15, 0x12, 0x0d}},
    {'R', {0x1e, 0x11, 0x11, 0x1e, 0x14, 0x12, 0x11}},
    {'S', {0x0f, 0x10, 0x10, 0x0e, 0x01, 0x01, 0x1e}},
    {'T', {0x1f, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04}},
    {'U', {0x11, 0x11, 0x11, 0x11, 0x11, 0x11, 0x0e}},
    {'V', {0x11, 0x11, 0x11, 0x11, 0x11, 0x0a, 0x04}},
    {'W', {0x11, 0x11, 0x11, 0x15, 0x15, 0x15, 0x0a}},
    {'X', {0x11, 0x11, 0x0a, 0x04, 0x0a, 0x11, 0x11}},
    {'Y', {0x11, 0x11, 0x0a, 0x04, 0x04, 0x04, 0x04}},
    {'Z', {0x1f, 0x01, 0x02, 0x04, 0x08, 0x10, 0x1f}},
    {'0', {0x0e, 0x11, 0x13, 0x15, 0x19, 0x11, 0x0e}},
    {'1', {0x04, 0x0c, 0x04, 0x04, 0x04, 0x04, 0x0e}},
    {'2', {0x0e, 0x11, 0x01, 0x02, 0x04, 0x08, 0x1f}},
    {'3', {0x1f, 0x02, 0x04, 0x02, 0x01, 0x11, 0x0e}},
    {'4', {0x02, 0x06, 0x0a, 0x12, 0x1f, 0x02, 0x02}},
    {'5', {0x1f, 0x10, 0x1e, 0x01, 0x01, 0x11, 0x0e}},
    {'6', {0x06, 0x08, 0x10, 0x1e, 0x11, 0x11, 0x0e}},
    {'7', {0x1f, 0x01, 0x02, 0x04, 0x08, 0x08, 0x08}},
    {'8', {0x0e, 0x11, 0x11, 0x0e, 0x11, 0x11, 0x0e}},
    {'9', {0x0e, 0x11, 0x11, 0x0f, 0x01, 0x02, 0x0c}},
    {' ', {0, 0, 0, 0, 0, 0, 0}},
    {'.', {0, 0, 0, 0, 0, 0x0c, 0x0c}},
    {',', {0, 0, 0, 0, 0x0c, 0x04, 0x08}},
    {':', {0, 0x0c, 0x0c, 0, 0x0c, 0x0c, 0}},
    {'-', {0, 0, 0, 0x1f, 0, 0, 0}},
    {'+', {0, 0x04, 0x04, 0x1f, 0x04, 0x04, 0}},
    {'=', {0, 0, 0x1f, 0, 0x1f, 0, 0}},
    {'/', {0x01, 0x02, 0x02, 0x04, 0x08, 0x08, 0x10}},
    {'%', {0x18, 0x19, 0x02, 0x04, 0x08, 0x13, 0x03}},
    {'(', {0x02, 0x04, 0x08, 0x08, 0x08, 0x04, 0x02}},
    {')', {0x08, 0x04, 0x02, 0x02, 0x02, 0x04, 0x08}},
    {']', {0x0e, 0x02, 0x02, 0x02, 0x02, 0x02, 0x0e}},
    {'[', {0x0e, 0x08, 0x08, 0x08, 0x08, 0x08, 0x0e}},
    {'<', {0x02, 0x04, 0x08, 0x10, 0x08, 0x04, 0x02}},
    {'>', {0x08, 0x04, 0x02, 0x01, 0x02, 0x04, 0x08}},
    {'_', {0, 0, 0, 0, 0, 0, 0x1f}},
    {'!', {0x04, 0x04, 0x04, 0x04, 0x04, 0, 0x04}},
    {'*', {0, 0x15, 0x0e, 0x1f, 0x0e, 0x15, 0}},
    {'#', {0x0a, 0x0a, 0x1f, 0x0a, 0x1f, 0x0a, 0x0a}},
    {'?', {0x0e, 0x11, 0x01, 0x02, 0x04, 0, 0x04}},
};

static const Glyph *find_glyph(char ch) {
    if (ch >= 'a' && ch <= 'z') ch = (char)(ch - 'a' + 'A');
    for (size_t i = 0; i < sizeof(glyphs) / sizeof(glyphs[0]); i++)
        if (glyphs[i].ch == ch) return &glyphs[i];
    return NULL;
}

int text_width(const char *s, int scale) {
    int n = (int)strlen(s);
    return n > 0 ? (n * 6 - 1) * scale : 0;
}

int text_height(const char *s, int scale) {
    (void)s;
    return 7 * scale;
}

void draw_text(Canvas *c, int x, int y, const char *s, int scale, Color color) {
    int pen = x;
    for (const char *p = s; *p; p++) {
        const Glyph *g = find_glyph(*p);
        if (g) {
            for (int row = 0; row < 7; row++) {
                unsigned bits = g->rows[row];
                for (int col = 0; col < 5; col++) {
                    if (bits & (1u << (4 - col)))
                        fill_rect(c, pen + col * scale, y + row * scale, scale, scale, color);
                }
            }
        }
        pen += 6 * scale;
    }
}

/* ------------------------------------------------------------------ */
/* PNG                                                                 */
/* ------------------------------------------------------------------ */

static void on_png_error(png_structp png, png_const_charp msg) {
    (void)png;
    fprintf(stderr, "png error: %s\n", msg);
}

int canvas_write_png(const Canvas *c, const char *path) {
    FILE *fp = fopen(path, "wb");
    if (!fp) {
        fprintf(stderr, "cannot open %s\n", path);
        return -1;
    }
    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, on_png_error, NULL);
    png_infop info = png ? png_create_info_struct(png) : NULL;
    if (!png || !info) {
        fclose(fp);
        return -1;
    }

    png_init_io(png, fp);
    png_set_IHDR(png, info, (png_uint_32)c->w, (png_uint_32)c->h, 8,
                 PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    png_set_sRGB(png, info, PNG_sRGB_INTENT_PERCEPTUAL);
    png_write_info(png, info);

    png_bytep *rows = (png_bytep *)malloc(sizeof(png_bytep) * (size_t)c->h);
    if (!rows) {
        png_destroy_write_struct(&png, &info);
        fclose(fp);
        return -1;
    }
    for (int y = 0; y < c->h; y++) rows[y] = (png_bytep)(c->px + (size_t)y * (size_t)c->w * 3u);

    png_write_image(png, rows);
    png_write_end(png, NULL);
    free(rows);
    png_destroy_write_struct(&png, &info);
    fclose(fp);
    return 0;
}
