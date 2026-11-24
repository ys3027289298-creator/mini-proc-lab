#!/bin/sh
# Negative control for the seam test: prove it detects a too-small halo.
#
# worldgen-shorthalo-512 is the 512x512 build with intentionally tiny halos
# (3..6 cells instead of 66..116).  Independently generated windows must then
# disagree with the full world near their edges.  The real build must agree.
#
# Expected: shipped 0 mismatching bytes, shorthalo > 0 mismatching bytes.
set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
SEED=${1:-20261124}
OUT="$ROOT/reports/halo-demo"
GOOD="$ROOT/worldgen/worldgen-512"
BAD="$ROOT/worldgen/worldgen-shorthalo-512"

if [ ! -x "$GOOD" ] || [ ! -x "$BAD" ]; then
    echo "demo: build first: make -C worldgen worldgen-512 worldgen-shorthalo-512" >&2
    exit 2
fi

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT INT TERM
mkdir -p "$OUT"

check_binary() {
    # $1 = binary, $2 = label; writes "$2.csv"; prints verdict; returns 0 if
    # every window matched the full world.
    err=0
    printf 'layer,mismatch_bytes,identical\n' > "$OUT/$2.csv"
    "$1" --seed "$SEED" --layer all --threads 1 --out "$WORK/$2-full" >/dev/null
    for spec in "0,0,64,64" "200,200,64,64" "260,300,64,64" "96,420,64,64"; do
        tx=$(printf '%s' "$spec" | tr ',' '_')
        "$1" --seed "$SEED" --layer all --threads 4 --window "$spec" \
             --out "$WORK/$2-$tx" >/dev/null
    done
    python3 -I - "$WORK" "$2" "$OUT/$2.csv" "$1" <<'PY' || err=1
import struct, subprocess, sys

work, label, csv_path, binary = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
world = int(subprocess.check_output([binary, "--dump-config"]).decode()
            .split("world ")[1].split()[0])
windows = [(0, 0, 64, 64), (200, 200, 64, 64), (260, 300, 64, 64),
           (96, 420, 64, 64)]
failed = False
rows = []
for name, size in (("noise.f32", 4), ("erosion.f32", 4),
                   ("biome.u8", 1), ("hydro.f32", 4)):
    full = open(f"{work}/{label}-full/{name}", "rb").read()
    bad = 0
    for (x0, y0, w, h) in windows:
        tag = f"{x0}_{y0}_{w}_{h}"
        tile = open(f"{work}/{label}-{tag}/{name}", "rb").read()
        for y in range(h):
            off = ((y0 + y) * world + x0) * size
            if full[off:off + w * size] != tile[y * w * size:(y + 1) * w * size]:
                bad += w * size
    rows.append((name, bad, bad == 0))
    if bad:
        failed = True


def points(path):
    data = open(path, "rb").read()
    n = struct.unpack("<I", data[:4])[0]
    return [struct.unpack("<ff", data[4 + i * 8:12 + i * 8]) for i in range(n)]


bad = 0
for (x0, y0, w, h) in windows:
    tag = f"{x0}_{y0}_{w}_{h}"
    full = points(f"{work}/{label}-full/scatter.f32")
    inside = [p for p in full if x0 <= p[0] < x0 + w and y0 <= p[1] < y0 + h]
    if inside != points(f"{work}/{label}-{tag}/scatter.f32"):
        bad += 1
rows.append(("scatter.f32", bad, bad == 0))
if bad:
    failed = True

with open(csv_path, "w") as f:
    f.write("layer,mismatch_bytes,identical\n")
    for name, bad, ok in rows:
        print(f"  {name:12s} {bad:>8d} mismatching  {'ok' if ok else 'SEAM'}")
        f.write(f"{name},{bad},{'yes' if ok else 'no'}\n")
sys.exit(1 if failed else 0)
PY
    return "$err"
}

echo "== seam negative control: 512x512 world, 64x64 windows =="
echo ""
echo "-- shipped worldgen-512 (real halos) --"
good_ok=1
check_binary "$GOOD" shipped || good_ok=0

echo ""
echo "-- worldgen-shorthalo-512 (halos 1..5 cells, deliberately wrong) --"
bad_ok=1
check_binary "$BAD" shorthalo || bad_ok=0

echo ""
if [ "$good_ok" = 1 ] && [ "$bad_ok" = 0 ]; then
    echo "RESULT: PASS -- real build has no seam; short-halo build does,"
    echo "        so the seam check is not vacuous."
    exit 0
fi
if [ "$good_ok" = 0 ]; then
    echo "RESULT: FAIL -- the real build shows a seam."
    exit 1
fi
echo "RESULT: FAIL -- the short-halo build produced no seam;"
echo "        the seam check is not actually testing anything."
exit 1
