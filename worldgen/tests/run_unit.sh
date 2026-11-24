#!/bin/sh
# Build and run the layer unit tests, then the documentation seed check.
set -eu

HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
WG=$(CDPATH= cd -- "$HERE/.." && pwd)
ROOT=$(CDPATH= cd -- "$WG/.." && pwd)

CC=${CC:-cc}
CFLAGS=${CFLAGS:--O2 -std=c11 -Wall -Wextra -ffp-contract=off -fno-fast-math}

BIN="$WG/tests/unit"
$CC -I"$WG/include" -I"$ROOT/common" $CFLAGS \
    "$WG/tests/unit.c" \
    "$WG/src/noise.c" "$WG/src/erosion.c" "$WG/src/biome.c" \
    "$WG/src/hydro.c" "$WG/src/scatter.c" "$WG/src/common.c" \
    "$ROOT/common/seed.c" "$ROOT/common/rng.c" "$ROOT/common/grid.c" \
    "$ROOT/common/parallel.c" \
    -o "$BIN" -lm -pthread

"$BIN"
"$ROOT/scripts/verify_seeds.sh"
