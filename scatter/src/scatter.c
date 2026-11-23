
/* scatter: deterministic jittered scatter with minimum distance rejection. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAXP 64

typedef struct { float x, y; } P;

typedef struct { P p[MAXP]; int n; } Set;

static unsigned h1(unsigned i) {
    unsigned h = i * 2654435761u;
    h ^= h >> 15;
    return h;
}

static float u01(unsigned i) { return (float)(h1(i) & 0xffff) / 65535.0f; }

static void set_init(Set* s) { memset(s, 0, sizeof(*s)); }

static int set_add(Set* s, float x, float y, float mind) {
    if (s->n >= MAXP) return -1;
    for (int i = 0; i < s->n; i++) {
        float dx = s->p[i].x - x, dy = s->p[i].y - y;
        if (dx * dx + dy * dy < mind * mind) return -2;
    }
    s->p[s->n].x = x; s->p[s->n].y = y;
    return s->n++;
}

static int jitter_grid(Set* s, int nx, int ny, float cell, float mind, unsigned seed) {
    int tries = 0;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            unsigned k = (unsigned)(i * 131 + j * 17) + seed * 7919u;
            float x = (i + u01(k)) * cell;
            float y = (j + u01(k + 1)) * cell;
            int r = set_add(s, x, y, mind);
            if (r == -2) { tries++; continue; }
            if (r == -1) return -1;
        }
    }
    return tries;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Set s;
    set_init(&s);
    int rejected = jitter_grid(&s, 4, 4, 1.0f, 0.0f, 1u);
    check("grid_count", s.n == 16 && rejected == 0);

    Set t;
    set_init(&t);
    rejected = jitter_grid(&t, 4, 4, 0.5f, 0.6f, 2u);
    check("min_distance_enforced", rejected > 0);
    int violations = 0;
    for (int i = 0; i < t.n; i++)
        for (int j = i + 1; j < t.n; j++) {
            float dx = t.p[i].x - t.p[j].x, dy = t.p[i].y - t.p[j].y;
            if (dx * dx + dy * dy < 0.36f - 1e-5f) violations++;
        }
    check("no_violation", violations == 0);

    check("deterministic", ({ Set a, b; set_init(&a); set_init(&b);
        jitter_grid(&a, 3, 3, 1.0f, 0.0f, 9u); jitter_grid(&b, 3, 3, 1.0f, 0.0f, 9u);
        a.n == b.n && fabsf(a.p[0].x - b.p[0].x) < 1e-9f; }));

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
