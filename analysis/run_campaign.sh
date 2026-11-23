#!/usr/bin/env bash
# Run the three iteration tiers (plus repeats, a size-scaling spot check and a
# second seed) for the erosion convergence study.
#
#   analysis/run_campaign.sh <outdir>
#
# Every run: n=256, drops=65536 (one droplet per cell), the tuned stable
# parameters from erosion/README.md, fixed droplet sample (--resample off) so
# the residual measures the pass operator rather than resampling noise.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$ROOT/erosion/erosion"
OUT="${1:?usage: run_campaign.sh <outdir>}"
mkdir -p "$OUT/logs" "$OUT/fields"

N=256
DROPS=65536
COMMON=(--n "$N" --drops "$DROPS" --seed 7
        --max-steps 48 --capacity 0.05 --evaporate 0.01
        --radius 2 --inertia 0.05 --erode 0.3 --deposit 0.3 --gravity 4)

run() {  # run <tag> <iters> [extra...]
    local tag="$1" iters="$2"; shift 2
    echo "run $tag iters=$iters" >&2
    "$BIN" "${COMMON[@]}" --iters "$iters" --out "$OUT/fields/$tag" "$@" \
        > "$OUT/logs/$tag.json" 2> "$OUT/logs/$tag.err"
}

# --- three iteration tiers: 10 / 30 / 90 (main measurement) ---
run t10 10
run t30 30
run t90 90

# --- repeatability: same tiers, fresh process, same seed => identical fields ---
run t30_rep 30
run t90_rep 90

# --- second seed at the top tier, to check the shape of the curve is not seed luck ---
"$BIN" --n "$N" --drops "$DROPS" --seed 11 \
        --max-steps 48 --capacity 0.05 --evaporate 0.01 \
        --radius 2 --inertia 0.05 --erode 0.3 --deposit 0.3 --gravity 4 \
        --iters 90 --out "$OUT/fields/t90_s11" \
        > "$OUT/logs/t90_s11.json" 2> "$OUT/logs/t90_s11.err"

# --- single-iteration cost vs field size (drops fixed at 1/cell) ---
for n in 128 256 512; do
    "$BIN" --n "$n" --drops $((n*n)) --seed 7 \
            --max-steps 48 --capacity 0.05 --evaporate 0.01 \
            --radius 2 --inertia 0.05 --erode 0.3 --deposit 0.3 --gravity 4 \
            --iters 5 --out "$OUT/fields/size$n" \
            > "$OUT/logs/size$n.json" 2> "$OUT/logs/size$n.err"
done

# --- cost vs drops-per-cell at n=256 ---
for d in 1 2 4; do
    "$BIN" --n "$N" --drops $((DROPS*d)) --seed 7 \
            --max-steps 48 --capacity 0.05 --evaporate 0.01 \
            --radius 2 --inertia 0.05 --erode 0.3 --deposit 0.3 --gravity 4 \
            --iters 5 --out "$OUT/fields/drops${d}x" \
            > "$OUT/logs/drops${d}x.json" 2> "$OUT/logs/drops${d}x.err"
done

echo "done -> $OUT" >&2
