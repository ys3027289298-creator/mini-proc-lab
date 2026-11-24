#!/bin/sh
# Checks docs/seed-derivation.md against the two implementations:
# the worked example table for master = 123456789 must appear verbatim in the
# C and Python outputs.
set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
BIN="$ROOT/worldgen/worldgen"

if [ ! -x "$BIN" ]; then
    echo "verify_seeds.sh: build first: make -C worldgen" >&2
    exit 2
fi

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT INT TERM

"$BIN" --dump-seeds 123456789 > "$WORK/c.txt"
python3 -I "$ROOT/tools/seed_reference.py" 123456789 > "$WORK/py.txt"

if ! cmp -s "$WORK/c.txt" "$WORK/py.txt"; then
    echo "C and Python disagree:" >&2
    diff -u "$WORK/py.txt" "$WORK/c.txt" >&2 || true
    exit 1
fi

# Expected values from docs/seed-derivation.md.
cat > "$WORK/expected.txt" <<'EOF'
master 123456789 0x00000000075bcd15
  noise     9577184121432529850  0x84e8fe37627f5bba
  scatter  13655677282127860315  0xbd82b623209a525b
  erosion   7746194062884428957  0x6b80027587775c9d
  biome      114714165651301812  0x01978bed0659e5b4
  hydro    13156795571732877193  0xb69653dd63920f89
EOF

if cmp -s "$WORK/c.txt" "$WORK/expected.txt"; then
    echo "seed derivation matches docs/seed-derivation.md"
else
    echo "seed derivation no longer matches the documented example:" >&2
    diff -u "$WORK/expected.txt" "$WORK/c.txt" >&2 || true
    exit 1
fi
