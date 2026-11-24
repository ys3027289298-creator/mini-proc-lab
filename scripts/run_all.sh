#!/bin/sh
# Runs every check and writes a consolidated report.
#
#   scripts/run_all.sh
#
# Produces:
#   reports/verify/                 per-layer repeat, threads 1/4/8, seams, seeds
#   reports/nondeterminism-demo/    negative control: planted bad reduction
set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"

echo "=== build ==="
make -C worldgen worldgen worldgen-nondet worldgen-512 worldgen-shorthalo-512
echo ""

echo "=== determinism / seams / seed derivation ==="
scripts/verify.sh --out "$ROOT/reports/verify"
echo ""

echo "=== negative control (planted non-deterministic reduction) ==="
scripts/nondeterminism_demo.sh
echo ""

echo "=== negative control (halos shortened to prove seam check is real) ==="
scripts/halo_demo.sh
echo ""

echo "reports:"
ls -1 "$ROOT/reports/verify" "$ROOT/reports/nondeterminism-demo" \
      "$ROOT/reports/halo-demo"
