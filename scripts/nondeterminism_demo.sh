#!/bin/sh
# Negative control: prove the determinism checks actually detect a
# thread-order-dependent reduction, and that they localize it to hydro.
#
# worldgen-nondet is worldgen plus one deliberate bug: the hydro layer scales
# its output by a mean computed as "one float partial per thread, added in
# completion order" (src/selftest_nondet.c).  Every other layer is untouched.
#
# Expected result: noise/erosion/biome/scatter identical across 1/4/8 threads,
# hydro differs.  This demonstrates that a failing check points at the layer
# that introduced the reduction.
set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
SEED=${1:-20261124}
OUT="$ROOT/reports/nondeterminism-demo"
BIN="$ROOT/worldgen/worldgen-nondet"
GOOD="$ROOT/worldgen/worldgen"

if [ ! -x "$BIN" ]; then
    echo "demo: build first: make -C worldgen worldgen-nondet" >&2
    exit 2
fi
if [ ! -x "$GOOD" ]; then
    echo "demo: build first: make -C worldgen" >&2
    exit 2
fi

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT INT TERM
mkdir -p "$OUT"

echo "== negative control: one deliberate non-deterministic reduction in hydro =="
echo "binary: worldgen-nondet  seed: $SEED"
echo ""
printf '%-10s %-8s %-16s %s\n' layer threads hash identical_to_1
: > "$OUT/threads.csv"
echo "layer,threads,sha256,identical_to_1" >> "$OUT/threads.csv"

suspect=""
for layer in noise erosion biome hydro scatter; do
    case "$layer" in
        noise) f=noise.f32 ;;
        erosion) f=erosion.f32 ;;
        biome) f=biome.u8 ;;
        hydro) f=hydro.f32 ;;
        scatter) f=scatter.f32 ;;
    esac
    ref=""
    for t in 1 4 8; do
        d="$WORK/${layer}_$t"
        "$BIN" --seed "$SEED" --layer "$layer" --threads "$t" --out "$d" >/dev/null
        h=$(sha256sum "$d/$f" | awk '{print $1}')
        if [ -z "$ref" ]; then
            ref=$h
            printf '%-10s %-8s %-16s %s\n' "$layer" "$t" "$(printf %.16s "$h")" "yes (reference)"
            echo "$layer,$t,$h,yes" >> "$OUT/threads.csv"
        else
            if [ "$h" = "$ref" ]; then
                printf '%-10s %-8s %-16s %s\n' "$layer" "$t" "$(printf %.16s "$h")" "yes"
                echo "$layer,$t,$h,yes" >> "$OUT/threads.csv"
            else
                printf '%-10s %-8s %-16s %s\n' "$layer" "$t" "$(printf %.16s "$h")" "NO"
                echo "$layer,$t,$h,no" >> "$OUT/threads.csv"
                suspect="$layer"
            fi
        fi
    done
done

echo ""
if [ "$suspect" = "hydro" ]; then
    echo "RESULT: hydro is the only layer that varies with thread count."
    echo "        Localization correct: hydro is where the unsafe reduction was linked."
elif [ -z "$suspect" ]; then
    echo "RESULT: FAILED to reproduce -- the planted bug is no longer detectable."
    exit 1
else
    echo "RESULT: unexpected suspect: $suspect"
    exit 1
fi

# Also confirm the shipped binary is clean on the same input.
echo ""
echo "== control: shipped worldgen must be identical across thread counts =="
: > "$OUT/shipped.csv"
echo "layer,threads,sha256,identical_to_1" >> "$OUT/shipped.csv"
clean=1
for layer in noise erosion biome hydro scatter; do
    case "$layer" in
        noise) f=noise.f32 ;;
        erosion) f=erosion.f32 ;;
        biome) f=biome.u8 ;;
        hydro) f=hydro.f32 ;;
        scatter) f=scatter.f32 ;;
    esac
    ref=""
    for t in 1 4 8; do
        d="$WORK/good_${layer}_$t"
        "$GOOD" --seed "$SEED" --layer "$layer" --threads "$t" --out "$d" >/dev/null
        h=$(sha256sum "$d/$f" | awk '{print $1}')
        [ -z "$ref" ] && ref=$h
        if [ "$h" = "$ref" ]; then
            echo "$layer,$t,$h,yes" >> "$OUT/shipped.csv"
        else
            echo "$layer,$t,$h,no" >> "$OUT/shipped.csv"
            clean=0
        fi
    done
done
if [ "$clean" = 1 ]; then
    echo "shipped worldgen: all layers identical at 1/4/8 threads"
else
    echo "shipped worldgen: NON-DETERMINISTIC -- investigate"
    exit 1
fi
