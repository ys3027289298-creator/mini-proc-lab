#ifndef ROCKS_GFX_H
#define ROCKS_GFX_H

#include <stdint.h>

/* Packed 0x00RRGGBB. */
typedef uint32_t Color;

#define RGB(r, g, b) \
    ((((uint32_t)(r) & 0xffu) << 16) | (((uint32_t)(g) & 0xffu) << 8) | ((uint32_t)(b) & 0xffu))

#define COL_PAGE      RGB(249, 249, 247)
#define COL_SURFACE   RGB(252, 252, 251)
#define COL_INK       RGB(11, 11, 11)
#define COL_INK2      RGB(82, 81, 78)
#define COL_MUTED     RGB(137, 135, 129)
#define COL_HAIR      RGB(225, 224, 217)
#define COL_BORDER    RGB(211, 209, 200)
#define COL_CRITICAL  RGB(208, 59, 59)
#define COL_WHITE     RGB(255, 255, 255)
#define COL_BLACK     RGB(0, 0, 0)

typedef struct {
    int w, h;
    uint8_t *px; /* RGB, row-major, 3 bytes per pixel */
} Canvas;

int canvas_init(Canvas *c, int w, int h);
void canvas_free(Canvas *c);
void canvas_clear(Canvas *c, Color color);
void canvas_copy(Canvas *dst, const Canvas *src);

void put_px(Canvas *c, int x, int y, Color color);
void fill_rect(Canvas *c, int x, int y, int w, int h, Color color);
void stroke_rect(Canvas *c, int x, int y, int w, int h, Color color);
void draw_line(Canvas *c, int x0, int y0, int x1, int y1, Color color);
void draw_circle(Canvas *c, int cx, int cy, int radius, Color color);
void draw_crosshair(Canvas *c, int cx, int cy, int arm, Color color);

/* 5x7 bitmap font, uppercase + digits + common punctuation. */
int text_width(const char *s, int scale);
int text_height(const char *s, int scale);
void draw_text(Canvas *c, int x, int y, const char *s, int scale, Color color);

/* Sequential blue ramp from the dataviz reference palette (surface -> 700). */
Color blue_ramp(float t);

int canvas_write_png(const Canvas *c, const char *path);

#endif
