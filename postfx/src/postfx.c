
/* postfx: per-channel affine post-processing stack with conservation audit.
 *
 * Each pass applies the same channel-wise affine map  out = a*in + b  to all
 * pixels of a channel. Composing such maps stays affine, so the whole stack
 * has an analytic closed form. Two pipelines share every pass except the last
 * one ("tint"), whose buggy variant overwrites the green channel instead of
 * writing a*in + b -- the classic "a channel got swallowed after stacking"
 * defect. Three independent numbers are compared per channel:
 *   measured  : energy of the image produced by the actual pass loop,
 *   analytic  : energy of the per-pixel analytic composite,
 *   book      : closed-form affine bookkeeping (dot products, no pixels).
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdint.h>

#define N 40
#define NPX (N * N)
#define CH 3
#define STAGES 5           /* input + 4 passes */
#define THRESH 1e-3        /* any channel relative error above this => fail */

typedef struct { float c[CH][NPX]; } Image;

/* ---- deterministic input ------------------------------------------------ */

static unsigned h32(unsigned x, unsigned y, unsigned k) {
    unsigned h = x * 374761393u + y * 668265263u + k * 2246822519u;
    h = (h ^ (h >> 13)) * 1274126177u;
    return h ^ (h >> 16);
}

static float u01(unsigned h) { return (float)(h & 0xffffff) / 16777215.0f; }

static void image_input(Image* im) {
    for (int y = 0; y < N; y++) {
        for (int x = 0; x < N; x++) {
            int p = y * N + x;
            im->c[0][p] = 0.5f + 0.5f * u01(h32(x, y, 11u));                 /* R: base */
            im->c[1][p] = 0.25f * u01(h32(x, y, 23u))                        /* G: faint stripes */
                        + 0.5f * (float)(x % 8 < 4);
            im->c[2][p] = 0.5f + 0.5f * sinf(0.35f * (float)x) * cosf(0.27f * (float)y);
        }
    }
}

/* ---- pass parameters, indexed by channel [R,G,B] ------------------------ */

static const float GAIN[CH]  = { 1.20f, 1.10f, 0.90f };
static const float LIFT[CH]  = { 0.05f, 0.10f, 0.15f };
static const float OVAMP[CH] = { 0.30f, 0.25f, 0.20f };
static const float TGAIN[CH] = { 1.05f, 0.80f, 1.15f };
static const float TOFF[CH]  = { 0.02f, 0.03f, 0.04f };

/* spatial overlay field b3(x,y); the same field is used by every pipeline
 * and by the analytic paths, so it is defined in exactly one place. */
static float overlay_field(int x, int y, int ch) {
    return OVAMP[ch] * sinf(0.31f * (float)x + (float)ch)
                     * cosf(0.23f * (float)y - (float)ch);
}

/* ---- the four passes (pixel loop) --------------------------------------- */

static void pass_gain(Image* im) {
    for (int ch = 0; ch < CH; ch++)
        for (int p = 0; p < NPX; p++)
            im->c[ch][p] *= GAIN[ch];
}

static void pass_lift(Image* im) {
    for (int ch = 0; ch < CH; ch++)
        for (int p = 0; p < NPX; p++)
            im->c[ch][p] += LIFT[ch];
}

static void pass_overlay(Image* im) {
    for (int ch = 0; ch < CH; ch++)
        for (int y = 0; y < N; y++)
            for (int x = 0; x < N; x++)
                im->c[ch][y * N + x] += overlay_field(x, y, ch);
}

/* tint: correct = a*in + b on every channel; buggy swallows green. */
static void pass_tint(Image* im, int buggy) {
    for (int ch = 0; ch < CH; ch++) {
        for (int p = 0; p < NPX; p++) {
            if (buggy && ch == 1)
                im->c[ch][p] = TOFF[ch];                 /* upstream G discarded */
            else
                im->c[ch][p] = im->c[ch][p] * TGAIN[ch] + TOFF[ch];
        }
    }
}

/* ---- analytic path 1: per-pixel affine composite -------------------------
 * Starting from the input pixel s, the stack is  out = A*s + B  with:
 *   gain    A*=g          B*=g
 *   lift    B += l
 *   overlay B += o(x,y)
 *   tint    A*=tg  B = B*tg + tb
 */

typedef struct { float A[CH], B[CH]; } Affine;

static void affine_init(Affine* af) {
    for (int ch = 0; ch < CH; ch++) { af->A[ch] = 1.0f; af->B[ch] = 0.0f; }
}

static void affine_gain(Affine* af) {
    for (int ch = 0; ch < CH; ch++) { af->A[ch] *= GAIN[ch]; af->B[ch] *= GAIN[ch]; }
}
static void affine_lift(Affine* af) {
    for (int ch = 0; ch < CH; ch++) af->B[ch] += LIFT[ch];
}
/* overlay leaves A and the constant part of B untouched; its spatial field
 * is folded into the per-pixel composite in analytic_pixel(). */
static void affine_tint(Affine* af) {
    for (int ch = 0; ch < CH; ch++) {
        af->A[ch] *= TGAIN[ch];
        af->B[ch] = af->B[ch] * TGAIN[ch] + TOFF[ch];
    }
}

static float analytic_pixel(int ch, int x, int y, const Affine* af) {
    float s;
    /* recompute the input pixel s from its definition (no reference to the
     * measured image, keeping the analytic path independent). */
    switch (ch) {
        case 0: s = 0.5f + 0.5f * u01(h32(x, y, 11u)); break;
        case 1: s = 0.25f * u01(h32(x, y, 23u)) + 0.5f * (float)(x % 8 < 4); break;
        default: s = 0.5f + 0.5f * sinf(0.35f * (float)x) * cosf(0.27f * (float)y); break;
    }
    /* overlay's spatial contribution survives tint scaling: o * tg */
    float Bsp = overlay_field(x, y, ch) * TGAIN[ch];
    return af->A[ch] * s + af->B[ch] + Bsp;
}

/* ---- analytic path 2: closed-form energy bookkeeping ---------------------
 * For image X under  Y = a X + b(x):
 *   E2(Y) = a^2 E2(X) + 2 a <X,b> + E2(b)
 * Only the running E2 scalar is kept per channel; each step's dot product
 * and b-moment are accumulated locally -- no output image is ever built.
 */

typedef struct { double e2[CH]; } Book;

static double energy2(const Image* im, int ch) {
    double e = 0.0;
    for (int p = 0; p < NPX; p++) { double v = im->c[ch][p]; e += v * v; }
    return e;
}

static void book_init(Book* bk, const Image* in) {
    for (int ch = 0; ch < CH; ch++) bk->e2[ch] = energy2(in, ch);
}

/* The generic affine bookkeeping is driven per pass with the actual b field
 * image built from parameters, so moments stay exact. */
static void book_step(Book* bk, int ch, float a, const Image* X, const Image* bfield) {
    double sqb = 0.0, dot = 0.0;
    for (int p = 0; p < NPX; p++) {
        double b = bfield->c[ch][p];
        sqb += b * b;
        dot += (double)X->c[ch][p] * b;
    }
    bk->e2[ch] = (double)a * (double)a * bk->e2[ch] + 2.0 * (double)a * dot + sqb;
}

/* ---- per-stage recording ------------------------------------------------- */

typedef struct {
    double m2[CH], a2[CH], k2[CH];   /* measured / analytic / book L2 energy */
} Stage;

static const char* STAGE_NAME[STAGES] = { "input", "gain", "lift", "overlay", "tint" };

static double reldiff(double a, double b) {
    double d = fabs(a - b);
    return d > 1e-12 ? d / fmax(fabs(a), fabs(b)) : 0.0;
}

/* ---- binary dump for the Python renderer -------------------------------- */

static void dump(const char* path, const Image* in,
                 const Image* ok, const Image* bug,
                 const Image* aok, const Image* abug,
                 const Stage* sok, const Stage* sbug) {
    FILE* f = fopen(path, "wb");
    if (!f) { perror(path); return; }
    int32_t hdr[6] = { 0x50465831, N, N, CH, STAGES, 0 }; /* 'PFX1' */
    fwrite(hdr, sizeof(int32_t), 6, f);
    fwrite(in->c,   sizeof(float), CH * NPX, f);
    fwrite(ok->c,   sizeof(float), CH * NPX, f);
    fwrite(bug->c,  sizeof(float), CH * NPX, f);
    fwrite(aok->c,  sizeof(float), CH * NPX, f);
    fwrite(abug->c, sizeof(float), CH * NPX, f);
    for (int s = 0; s < STAGES; s++) {
        fwrite(sok[s].m2,  sizeof(double), CH, f);
        fwrite(sok[s].a2,  sizeof(double), CH, f);
        fwrite(sok[s].k2,  sizeof(double), CH, f);
        fwrite(sbug[s].m2, sizeof(double), CH, f);
        fwrite(sbug[s].a2, sizeof(double), CH, f);
        fwrite(sbug[s].k2, sizeof(double), CH, f);
    }
    fclose(f);
}

/* ---- self test ----------------------------------------------------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(int argc, char** argv) {
    const char* outpath = argc > 1 ? argv[1] : "postfx.bin";

    Image in, imgOk, imgBug, anOk, anBug;
    image_input(&in);

    /* b fields for closed-form bookkeeping, one image per pass type. */
    Image bgain, blift, bover, btint;
    memset(&bgain, 0, sizeof(bgain)); memset(&blift, 0, sizeof(blift));
    memset(&bover, 0, sizeof(bover)); memset(&btint, 0, sizeof(btint));
    for (int ch = 0; ch < CH; ch++)
        for (int p = 0; p < NPX; p++) {
            int x = p % N, y = p / N;
            blift.c[ch][p] = LIFT[ch];
            bover.c[ch][p] = overlay_field(x, y, ch);
            btint.c[ch][p] = TOFF[ch];
        }

    Stage sok[STAGES], sbug[STAGES];
    Book bkOk, bkBug;
    Affine afOk, afBug;

    /* run both pipelines stage by stage */
    imgOk = in; imgBug = in; anOk = in; anBug = in;
    affine_init(&afOk); affine_init(&afBug);
    book_init(&bkOk, &in); book_init(&bkBug, &in);

    for (int ch = 0; ch < CH; ch++) {
        sok[0].m2[ch] = sbug[0].m2[ch] = energy2(&in, ch);
        sok[0].a2[ch] = sbug[0].a2[ch] = energy2(&in, ch);
        sok[0].k2[ch] = sbug[0].k2[ch] = bkOk.e2[ch];
    }

    /* pass 1: gain  (Y = a X, b = 0) */
    pass_gain(&imgOk); pass_gain(&imgBug);
    affine_gain(&afOk); affine_gain(&afBug);
    for (int ch = 0; ch < CH; ch++) {
        book_step(&bkOk, ch, GAIN[ch], &in, &bgain);
        book_step(&bkBug, ch, GAIN[ch], &in, &bgain);
    }
    /* pass 2: lift (Y = X + l, a = 1) */
    pass_lift(&imgOk); pass_lift(&imgBug);
    affine_lift(&afOk); affine_lift(&afBug);
    {
        Image afterGain = in; pass_gain(&afterGain);
        for (int ch = 0; ch < CH; ch++) {
            book_step(&bkOk, ch, 1.0f, &afterGain, &blift);
            book_step(&bkBug, ch, 1.0f, &afterGain, &blift);
        }
    }
    /* pass 3: overlay (Y = X + o(x,y), a = 1); affine coefficients unchanged */
    pass_overlay(&imgOk); pass_overlay(&imgBug);
    {
        Image afterLift = in; pass_gain(&afterLift); pass_lift(&afterLift);
        for (int ch = 0; ch < CH; ch++) {
            book_step(&bkOk, ch, 1.0f, &afterLift, &bover);
            book_step(&bkBug, ch, 1.0f, &afterLift, &bover);
        }
    }
    /* pass 4: tint (Y = a X + tb); buggy green swallows upstream */
    pass_tint(&imgOk, 0); pass_tint(&imgBug, 1);
    affine_tint(&afOk); affine_tint(&afBug);
    {
        Image afterOver = in;
        pass_gain(&afterOver); pass_lift(&afterOver); pass_overlay(&afterOver);
        for (int ch = 0; ch < CH; ch++) {
            book_step(&bkOk, ch, TGAIN[ch], &afterOver, &btint);
            book_step(&bkBug, ch, TGAIN[ch], &afterOver, &btint);
        }
    }

    /* analytic per-pixel composites (both stacks share the intended map) */
    for (int ch = 0; ch < CH; ch++)
        for (int y = 0; y < N; y++)
            for (int x = 0; x < N; x++) {
                int p = y * N + x;
                anOk.c[ch][p]  = analytic_pixel(ch, x, y, &afOk);
                anBug.c[ch][p] = analytic_pixel(ch, x, y, &afBug);
            }

    /* record final-stage energies */
    for (int ch = 0; ch < CH; ch++) {
        sok[STAGES - 1].m2[ch] = energy2(&imgOk, ch);
        sok[STAGES - 1].a2[ch] = energy2(&anOk, ch);
        sok[STAGES - 1].k2[ch] = bkOk.e2[ch];
        sbug[STAGES - 1].m2[ch] = energy2(&imgBug, ch);
        sbug[STAGES - 1].a2[ch] = energy2(&anBug, ch);
        sbug[STAGES - 1].k2[ch] = bkBug.e2[ch];
    }
    /* intermediate stages 1..3: measured vs analytic built from the running
     * affine state are identical maps, so analytic images are rebuilt there. */
    {
        Affine tmp; affine_init(&tmp);
        Image probe = in;
        for (int st = 1; st <= 3; st++) {
            if (st == 1) { affine_gain(&tmp); pass_gain(&probe); }
            if (st == 2) { affine_lift(&tmp); pass_lift(&probe); }
            if (st == 3) { pass_overlay(&probe); }
            Image ap;
            for (int ch = 0; ch < CH; ch++)
                for (int y = 0; y < N; y++)
                    for (int x = 0; x < N; x++) {
                        int p = y * N + x;
                        float Bsp = st >= 3 ? overlay_field(x, y, ch) : 0.0f;
                        float s;
                        switch (ch) {
                            case 0: s = 0.5f + 0.5f * u01(h32(x, y, 11u)); break;
                            case 1: s = 0.25f * u01(h32(x, y, 23u)) + 0.5f * (float)(x % 8 < 4); break;
                            default: s = 0.5f + 0.5f * sinf(0.35f * (float)x) * cosf(0.27f * (float)y); break;
                        }
                        ap.c[ch][p] = tmp.A[ch] * s + tmp.B[ch] + Bsp;
                    }
            for (int ch = 0; ch < CH; ch++) {
                sok[st].m2[ch] = sbug[st].m2[ch] = energy2(&probe, ch);
                sok[st].a2[ch] = sbug[st].a2[ch] = energy2(&ap, ch);
            }
        }
        /* fill book columns for stages 1..3 with independent moment recursion */
        Book bk; book_init(&bk, &in);
        Image aGain = in; pass_gain(&aGain);
        Image aLift = aGain; pass_lift(&aLift);
        for (int ch = 0; ch < CH; ch++) book_step(&bk, ch, GAIN[ch], &in, &bgain);
        for (int ch = 0; ch < CH; ch++)
            sok[1].k2[ch] = sbug[1].k2[ch] = bk.e2[ch];
        for (int ch = 0; ch < CH; ch++) book_step(&bk, ch, 1.0f, &aGain, &blift);
        for (int ch = 0; ch < CH; ch++)
            sok[2].k2[ch] = sbug[2].k2[ch] = bk.e2[ch];
        for (int ch = 0; ch < CH; ch++) book_step(&bk, ch, 1.0f, &aLift, &bover);
        for (int ch = 0; ch < CH; ch++)
            sok[3].k2[ch] = sbug[3].k2[ch] = bk.e2[ch];
    }

    /* ---- tables --------------------------------------------------------- */
    static const char* CN[CH] = { "R", "G", "B" };
    for (int bug = 0; bug < 2; bug++) {
        const Stage* st = bug ? sbug : sok;
        printf("\n=== %s stack: L2 energy by stage (measured | analytic | book | rel-err m/a | rel-err m/k) ===\n",
               bug ? "BUGGY" : "CONSERVING");
        for (int s = 0; s < STAGES; s++)
            for (int ch = 0; ch < CH; ch++)
                printf("%-8s %s  %14.6f %14.6f %14.6f   %.3e   %.3e\n",
                       STAGE_NAME[s], CN[ch], st[s].m2[ch], st[s].a2[ch], st[s].k2[ch],
                       reldiff(st[s].m2[ch], st[s].a2[ch]),
                       reldiff(st[s].m2[ch], st[s].k2[ch]));
    }

    /* ---- conservation verdicts ------------------------------------------ */
    float maxOk[CH] = {0, 0, 0}, maxBug[CH] = {0, 0, 0};
    for (int ch = 0; ch < CH; ch++) {
        for (int p = 0; p < NPX; p++) {
            float dOk = fabsf(imgOk.c[ch][p] - anOk.c[ch][p]);
            float dBug = fabsf(imgBug.c[ch][p] - anBug.c[ch][p]);
            float nrm = fmaxf(1e-6f, fmaxf(fabsf(anOk.c[ch][p]), fabsf(anBug.c[ch][p])));
            float rOk = dOk / nrm, rBug = dBug / nrm;
            if (rOk > maxOk[ch]) maxOk[ch] = rOk;
            if (rBug > maxBug[ch]) maxBug[ch] = rBug;
        }
    }

    printf("\n=== final-stage channel energy change (input -> output) ===\n");
    for (int bug = 0; bug < 2; bug++) {
        const Stage* st = bug ? sbug : sok;
        for (int ch = 0; ch < CH; ch++)
            printf("%-10s %s  E2 %10.4f -> %10.4f  (x%8.5f)   analytic %10.4f  rel-err %.3e\n",
                   bug ? "buggy" : "conserving", CN[ch],
                   st[0].m2[ch], st[4].m2[ch], st[4].m2[ch] / st[0].m2[ch],
                   st[4].a2[ch], reldiff(st[4].m2[ch], st[4].a2[ch]));
    }

    printf("\n=== per-pixel max relative difference vs analytic composite ===\n");
    int okPass = 1, bugCaught = 0;
    for (int ch = 0; ch < CH; ch++) {
        int o = maxOk[ch] <= THRESH, b = maxBug[ch] > THRESH;
        printf("conserving %s maxrel=%.3e %s\n", CN[ch], maxOk[ch], o ? "PASS" : "FAIL");
        printf("buggy     %s maxrel=%.3e %s\n", CN[ch], maxBug[ch], b ? "DETECTED" : "clean");
        if (!o) okPass = 0;
        if (b) bugCaught = 1;
    }
    printf("\n=== final energy relative error (threshold %.0e) ===\n", THRESH);
    int energyOk = 1;
    for (int ch = 0; ch < CH; ch++) {
        double e = reldiff(sok[4].m2[ch], sok[4].a2[ch]);
        double k = reldiff(sok[4].m2[ch], sok[4].k2[ch]);
        int pass = e <= THRESH && k <= THRESH;
        if (!pass) energyOk = 0;
        printf("conserving %s measured-vs-analytic=%.3e measured-vs-book=%.3e %s\n",
               CN[ch], e, k, pass ? "PASS" : "FAIL");
    }
    for (int ch = 0; ch < CH; ch++) {
        double e = reldiff(sbug[4].m2[ch], sbug[4].a2[ch]);
        printf("buggy     %s measured-vs-analytic=%.3e %s\n",
               CN[ch], e, e <= THRESH ? "PASS(channel intact)" : "FAIL(swallow detected)");
    }

    check("conserving_stack_passes", okPass && energyOk);
    check("buggy_swallow_detected", bugCaught);
    check("green_is_swallowed_channel",
          reldiff(sbug[4].m2[1], sbug[4].a2[1]) > THRESH &&
          reldiff(sbug[4].m2[0], sbug[4].a2[0]) <= THRESH &&
          reldiff(sbug[4].m2[2], sbug[4].a2[2]) <= THRESH);

    dump(outpath, &in, &imgOk, &imgBug, &anOk, &anBug, sok, sbug);
    printf("\ndump=%s\nall_ok=%d\n", outpath, failures == 0);
    return failures == 0 ? 0 : 1;
}
