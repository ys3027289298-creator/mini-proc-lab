#!/bin/sh
# Determinism, seam and seed-derivation checks for worldgen.
#
#   1. per layer: same seed run twice, output must be bit-identical
#   2. thread counts 1/4/8: every layer must be bit-identical
#   3. seams: 4x4 tiles of 64x64 must stitch into the exact full world
#   4. master seed: all five outputs change, and so do all five sub-seeds
#   5. sub-seed derivation: C output must match tools/seed_reference.py
#
# Usage: scripts/verify.sh [--seed N] [--out DIR] [--bin PATH]
# Exit status 0 iff all checks pass.  Writes CSV/Markdown reports to DIR.
set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
SEED=20261124
OUT="$ROOT/reports/verify"
BIN="$ROOT/worldgen/worldgen"

while [ $# -gt 0 ]; do
    case "$1" in
        --seed) SEED=$2; shift 2 ;;
        --out) OUT=$2; shift 2 ;;
        --bin) BIN=$2; shift 2 ;;
        -h|--help)
            sed -n '2,13p' "$0" | sed 's/^# \{0,1\}//'
            exit 0 ;;
        *) echo "verify.sh: unknown argument: $1" >&2; exit 2 ;;
    esac
done

if [ ! -x "$BIN" ]; then
    echo "verify.sh: build first: make -C worldgen" >&2
    exit 2
fi

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT INT TERM
mkdir -p "$OUT"

LAYERS="noise erosion biome hydro scatter"

file_for() {
    case "$1" in
        noise) echo noise.f32 ;;
        erosion) echo erosion.f32 ;;
        biome) echo biome.u8 ;;
        hydro) echo hydro.f32 ;;
        scatter) echo scatter.f32 ;;
    esac
}
say() { printf '%s\n' "$*"; }

fail=0
note_fail() { fail=1; }

#############################
# 1. same seed, two runs    #
#############################
say "== 1. per-layer repeat determinism (seed $SEED, threads 4, 2 runs) =="
printf 'layer,run_a_sha256,run_b_sha256,bytes,identical\n' > "$OUT/repeat.csv"
for layer in $LAYERS; do
    f=$(file_for "$layer")
    a="$WORK/${layer}_a"
    b="$WORK/${layer}_b"
    "$BIN" --seed "$SEED" --layer "$layer" --threads 4 --out "$a" >/dev/null
    "$BIN" --seed "$SEED" --layer "$layer" --threads 4 --out "$b" >/dev/null
    ha=$(sha256sum "$a/$f" | awk '{print $1}')
    hb=$(sha256sum "$b/$f" | awk '{print $1}')
    n=$(wc -c < "$a/$f")
    if cmp -s "$a/$f" "$b/$f"; then
        say "  $layer: IDENTICAL ($n bytes) $ha"
        echo "$layer,$ha,$hb,$n,yes" >> "$OUT/repeat.csv"
    else
        off=$(cmp "$a/$f" "$b/$f" 2>/dev/null | sed -n 's/.*byte \([0-9]*\).*/\1/p')
        say "  $layer: DIFFERS first byte ${off:-?} ($n bytes)"
        echo "$layer,$ha,$hb,$n,no" >> "$OUT/repeat.csv"
        note_fail
    fi
done

################################
# 2. thread counts 1 / 4 / 8   #
################################
say ""
say "== 2. thread-count invariance (1, 4, 8), each layer and combined 'all' =="
printf 'layer,threads,sha256,identical_to_threads_1\n' > "$OUT/threads.csv"
for layer in $LAYERS; do
    f=$(file_for "$layer")
    ref=""
    for t in 1 4 8; do
        d="$WORK/th_${layer}_$t"
        "$BIN" --seed "$SEED" --layer "$layer" --threads "$t" --out "$d" >/dev/null
        h=$(sha256sum "$d/$f" | awk '{print $1}')
        if [ -z "$ref" ]; then
            ref=$h
            say "  $layer threads=1  $h  (reference)"
            echo "$layer,$t,$h,yes" >> "$OUT/threads.csv"
        elif [ "$h" = "$ref" ]; then
            say "  $layer threads=$t  $h  identical"
            echo "$layer,$t,$h,yes" >> "$OUT/threads.csv"
        else
            say "  $layer threads=$t  $h  DIFFERS  <-- suspect layer"
            echo "$layer,$t,$h,no" >> "$OUT/threads.csv"
            note_fail
        fi
    done
done

for t in 1 4 8; do
    "$BIN" --seed "$SEED" --layer all --threads "$t" --out "$WORK/all_$t" >/dev/null
done
for layer in $LAYERS; do
    f=$(file_for "$layer")
    for t in 4 8; do
        if ! cmp -s "$WORK/all_1/$f" "$WORK/all_$t/$f"; then
            say "  all threads=$t: $f differs from threads=1  <-- suspect layer"
            note_fail
        fi
    done
done
say "  combined 'all': checked at 1/4/8"

#############################
# 3. seams: 4x4 tile stitch #
#############################
say ""
SEAM_BIN="$BIN"
if [ -x "$ROOT/worldgen/worldgen-512" ]; then
    SEAM_BIN=$ROOT/worldgen/worldgen-512
    say "== 3. seam check: 8x8 tiles of 64x64 stitched vs full 512x512 world =="
    say "   (a 64x64 window plus the 116-cell hydro halo stays inside the world,"
    say "    so this tests the halo itself rather than world-edge clamping)"
else
    say "== 3. seam check: 4x4 tiles of 64x64 stitched vs full 256x256 world =="
    say "   (build worldgen-512 for the stronger version of this test)"
fi
WORLD=$("$SEAM_BIN" --dump-config | awk '$1=="world" {print $2}')
NT=$((WORLD / 64))

printf 'layer,tile_bytes,mismatch_bytes,identical\n' > "$OUT/seams.csv"
"$SEAM_BIN" --seed "$SEED" --layer all --threads 1 --out "$WORK/full" >/dev/null
TILES="$WORK/tiles"
mkdir -p "$TILES"
ty=0
while [ "$ty" -lt "$NT" ]; do
    tx=0
    while [ "$tx" -lt "$NT" ]; do
        "$SEAM_BIN" --seed "$SEED" --layer all --threads 4 \
               --window "$((tx * 64)),$((ty * 64)),64,64" \
               --out "$TILES/${tx}_${ty}" >/dev/null
        tx=$((tx + 1))
    done
    ty=$((ty + 1))
done

python3 -I - "$WORK/full" "$TILES" "$OUT/seams.csv" "$WORLD" "$NT" \
    <<'PY' || note_fail
import struct, sys

full_dir, tiles_dir, csv_path = sys.argv[1], sys.argv[2], sys.argv[3]
WORLD, NT = int(sys.argv[4]), int(sys.argv[5])
TILE = WORLD // NT
failed = False
rows = []


def grid(path, size):
    return open(path, "rb").read()


for name, size in (("noise.f32", 4), ("erosion.f32", 4),
                   ("biome.u8", 1), ("hydro.f32", 4)):
    full = grid(f"{full_dir}/{name}", size)
    stitched = bytearray(WORLD * WORLD * size)
    for ty in range(NT):
        for tx in range(NT):
            tile = grid(f"{tiles_dir}/{tx}_{ty}/{name}", size)
            for y in range(TILE):
                src = (y * TILE) * size
                dst = ((ty * TILE + y) * WORLD + tx * TILE) * size
                stitched[dst:dst + TILE * size] = tile[src:src + TILE * size]
    bad = sum(1 for a, b in zip(full, stitched) if a != b)
    rows.append((name, len(stitched), bad, "yes" if bad == 0 else "no"))
    if bad:
        failed = True


def points(path):
    data = open(path, "rb").read()
    n = struct.unpack("<I", data[:4])[0]
    return [struct.unpack("<ff", data[4 + i * 8:12 + i * 8]) for i in range(n)]


full_pts = sorted(points(f"{full_dir}/scatter.f32"))
tile_pts = []
for ty in range(NT):
    for tx in range(NT):
        tile_pts.extend(points(f"{tiles_dir}/{tx}_{ty}/scatter.f32"))
tile_pts.sort()
bad = 0 if full_pts == tile_pts else abs(len(full_pts) - len(tile_pts)) + 1
rows.append(("scatter.f32", 4 + 8 * len(tile_pts), bad,
             "yes" if full_pts == tile_pts else "no"))
if full_pts != tile_pts:
    failed = True

with open(csv_path, "w") as f:
    f.write("layer,tile_bytes,mismatch_bytes,identical\n")
    for name, nbytes, bad, ok in rows:
        detail = f"{len(full_pts)} points" if name == "scatter.f32" else f"{nbytes} bytes"
        print(f"  {name}: {detail}, {bad} mismatching, "
              f"{'IDENTICAL' if ok == 'yes' else 'DIFFERS'}")
        f.write(f"{name},{nbytes},{bad},{ok}\n")
sys.exit(1 if failed else 0)
PY

################################
# 4. master seed changes all   #
################################
say ""
say "== 4. different master seed moves all five layers =="
printf 'layer,seed_a,seed_b,sha256_a,sha256_b,changed\n' > "$OUT/seeda.csv"
OTHER=$((SEED ^ 1))
"$BIN" --seed "$SEED"  --layer all --threads 4 --out "$WORK/seed_a" >/dev/null
"$BIN" --seed "$OTHER" --layer all --threads 4 --out "$WORK/seed_b" >/dev/null
for layer in $LAYERS; do
    f=$(file_for "$layer")
    ha=$(sha256sum "$WORK/seed_a/$f" | awk '{print $1}')
    hb=$(sha256sum "$WORK/seed_b/$f" | awk '{print $1}')
    if [ "$ha" != "$hb" ]; then
        say "  $layer: changed"
        echo "$layer,$SEED,$OTHER,$ha,$hb,yes" >> "$OUT/seeda.csv"
    else
        say "  $layer: UNCHANGED  <-- sub-seed not wired"
        echo "$layer,$SEED,$OTHER,$ha,$hb,no" >> "$OUT/seeda.csv"
        note_fail
    fi
done

#############################
# 5. seed reference match   #
#############################
say ""
say "== 5. sub-seed derivation: C vs Python reference =="
"$BIN" --dump-seeds "$SEED" "$OTHER" 0 1 123456789 0xDEADBEEFCAFEF00D \
    > "$WORK/c_seeds.txt"
python3 -I "$ROOT/tools/seed_reference.py" \
    "$SEED" "$OTHER" 0 1 123456789 0xDEADBEEFCAFEF00D > "$WORK/py_seeds.txt"
if cmp -s "$WORK/c_seeds.txt" "$WORK/py_seeds.txt"; then
    say "  IDENTICAL"
    cp "$WORK/c_seeds.txt" "$OUT/seeds.txt"
else
    say "  DIFFERS"
    diff -u "$WORK/py_seeds.txt" "$WORK/c_seeds.txt" || true
    note_fail
fi

# ---------------------------------------------------------------- summary
say ""
if [ "$fail" = 0 ]; then
    say "RESULT: PASS"
else
    say "RESULT: FAIL"
fi

# One report per layer: two runs of the same seed, plus the 1/4/8 thread
# comparison.  This is the "one comparison result per layer" deliverable.
for layer in $LAYERS; do
    f=$(file_for "$layer")
    rep="$OUT/layer-$layer.md"
    ha=$(awk -F, -v l="$layer" '$1==l {print $2}' "$OUT/repeat.csv")
    hb=$(awk -F, -v l="$layer" '$1==l {print $3}' "$OUT/repeat.csv")
    rn=$(awk -F, -v l="$layer" '$1==l {print $4}' "$OUT/repeat.csv")
    rok=$(awk -F, -v l="$layer" '$1==l {print $5}' "$OUT/repeat.csv")
    {
        echo "# Layer: $layer"
        echo ""
        echo "Master seed: \`$SEED\`  |  binary: \`worldgen\`  |  output: \`$f\`"
        echo ""
        echo "## Two runs, same seed"
        echo ""
        echo "| run | bytes | sha256 |"
        echo "|---|---:|---|"
        echo "| A | $rn | \`$ha\` |"
        echo "| B | $rn | \`$hb\` |"
        echo ""
        if [ "$rok" = "yes" ]; then
            echo "**Verdict: bit-identical** (\`cmp\` clean)."
        else
            echo "**Verdict: DIFFERS.**"
        fi
        echo ""
        echo "## Thread counts"
        echo ""
        echo "| threads | sha256 | identical to threads=1 |"
        echo "|---|---|---|"
        for t in 1 4 8; do
            h=$(awk -F, -v l="$layer" -v t="$t" '$1==l && $2==t {print $3}' "$OUT/threads.csv")
            ok=$(awk -F, -v l="$layer" -v t="$t" '$1==l && $2==t {print $4}' "$OUT/threads.csv")
            echo "| $t | \`$h\` | $ok |"
        done
        echo ""
        if awk -F, -v l="$layer" '$1==l && $4=="no"' "$OUT/threads.csv" | grep -q .; then
            echo "**Verdict: NON-DETERMINISTIC.** This layer changes with thread"
            echo "count. Look for a cross-thread floating-point or integer"
            echo "reduction (partial sums combined in completion order, atomic"
            echo "adds, or work-stealing that reorders accumulation)."
        else
            echo "**Verdict: bit-identical at 1, 4 and 8 threads.**"
        fi
    } > "$rep"
done
say "per-layer reports: $OUT/layer-*.md"

{
    echo "# worldgen verification report"
    echo ""
    echo "- master seed: \`$SEED\`"
    echo "- binary: \`$BIN\`"
    echo "- result: **$([ "$fail" = 0 ] && echo PASS || echo FAIL)**"
    echo ""
    echo "## Per-layer repeat (same seed, two runs, threads 4)"
    echo ""
    echo "| layer | identical | sha256 |"
    echo "|---|---|---|"
    awk -F, 'NR>1 {printf "| %s | %s | `%s` |\n", $1, ($5=="yes"?"yes":"**NO**"), $2}' "$OUT/repeat.csv"
    echo ""
    echo "## Thread counts 1 / 4 / 8"
    echo ""
    echo "| layer | threads=1 | threads=4 | threads=8 |"
    echo "|---|---|---|---|"
    for layer in $LAYERS; do
        a=$(awk -F, -v l="$layer" '$1==l && $2==1 {print $4}' "$OUT/threads.csv")
        b=$(awk -F, -v l="$layer" '$1==l && $2==4 {print $4}' "$OUT/threads.csv")
        c=$(awk -F, -v l="$layer" '$1==l && $2==8 {print $4}' "$OUT/threads.csv")
        echo "| $layer | $a | $b | $c |"
    done
    echo ""
    echo "## Seams (64x64 tiles stitched vs full world)"
    echo ""
    echo "| layer | identical | mismatching bytes |"
    echo "|---|---|---|"
    awk -F, 'NR>1 {printf "| %s | %s | %s |\n", $1, ($4=="yes"?"yes":"**NO**"), $3}' "$OUT/seams.csv"
    echo ""
    echo "## Master seed changes every layer"
    echo ""
    echo "| layer | changed |"
    echo "|---|---|"
    awk -F, 'NR>1 {printf "| %s | %s |\n", $1, ($6=="yes"?"yes":"**NO**")}' "$OUT/seeda.csv"
} > "$OUT/summary.md"

say "reports written to $OUT"
exit "$fail"
