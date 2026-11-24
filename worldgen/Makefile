CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -ffp-contract=off -fno-fast-math
CPPFLAGS ?= -Iinclude -I../common
LDLIBS ?= -lm -pthread

COMMON := ../common/seed.c ../common/rng.c ../common/grid.c ../common/parallel.c
LAYERS := src/noise.c src/erosion.c src/biome.c src/hydro.c src/scatter.c \
          src/common.c
HEADERS := include/worldgen.h ../common/seed.h ../common/rng.h \
           ../common/grid.h ../common/parallel.h

all: worldgen

worldgen: src/main.c $(LAYERS) $(COMMON) $(HEADERS)
	$(CC) $(CPPFLAGS) $(CFLAGS) src/main.c $(LAYERS) $(COMMON) -o $@ $(LDLIBS)

# Negative control for scripts/nondeterminism_demo.sh.  The only difference is
# -DSP_SELFTEST_NONDET, which makes hydro scale its output by a mean computed
# with a thread-order-dependent reduction.  Do not ship this binary.
worldgen-nondet: src/main.c $(LAYERS) src/selftest_nondet.c $(COMMON) $(HEADERS)
	$(CC) $(CPPFLAGS) $(CFLAGS) -DSP_SELFTEST_NONDET src/main.c $(LAYERS) \
	      src/selftest_nondet.c $(COMMON) -o $@ $(LDLIBS)

# A world large enough that a 64x64 window plus its halo does not touch the
# world edge, so the seam test really exercises the halo instead of being
# masked by edge clamping.
worldgen-512: src/main.c $(LAYERS) $(COMMON) $(HEADERS)
	$(CC) $(CPPFLAGS) $(CFLAGS) -DWG_WORLD=512 src/main.c $(LAYERS) $(COMMON) \
	      -o $@ $(LDLIBS)

# Negative control for seam detection: correct world size, deliberately short
# halos.  Must FAIL the seam check.  Test-only.
worldgen-shorthalo-512: src/main.c $(LAYERS) $(COMMON) $(HEADERS)
	$(CC) $(CPPFLAGS) $(CFLAGS) -DWG_WORLD=512 -DSP_SELFTEST_SHORT_HALO \
	      src/main.c $(LAYERS) $(COMMON) -o $@ $(LDLIBS)

test: worldgen
	./tests/run_unit.sh

clean:
	rm -f worldgen worldgen-nondet worldgen-512 worldgen-shorthalo-512 tests/unit

.PHONY: all test clean
