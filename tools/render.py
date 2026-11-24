"""Render worldgen layers to PNG for visual inspection of seams.

No third-party libraries: PPM is written directly and converted with ffmpeg
if it is available, otherwise the PPM files are left in place.

Usage:  python3 -I tools/render.py OUTDIR [--seed N] [--scale N]
"""
import math
import os
import shutil
import struct
import subprocess
import sys

WORLD = 256

BIOME_COLORS = {
    0: (40, 70, 140),    # water
    1: (196, 174, 96),   # desert
    2: (40, 100, 45),    # forest
    3: (110, 150, 60),   # grass
    4: (190, 200, 210),  # tundra
    5: (110, 100, 95),   # rock
}


def read_f32(path):
    data = open(path, "rb").read()
    n = len(data) // 4
    return struct.unpack("<%df" % n, data)


def read_u8(path):
    return list(open(path, "rb").read())


def shade(values, lo=None, hi=None):
    lo = min(values) if lo is None else lo
    hi = max(values) if hi is None else hi
    span = (hi - lo) or 1.0
    return [int(255.0 * (v - lo) / span) for v in values]


def write_ppm(path, pixels, w, h, scale):
    with open(path, "wb") as f:
        f.write(b"P6\n%d %d\n255\n" % (w * scale, h * scale))
        for y in range(h):
            row = bytearray()
            for x in range(w):
                r, g, b = pixels[y * w + x]
                row += bytes((r, g, b)) * scale
            for _ in range(scale):
                f.write(row)


def main(argv):
    outdir = argv[1] if len(argv) > 1 else "reports/render"
    seed = 20261124
    scale = 2
    keep_raw = False
    args = argv[2:]
    while args:
        if args[0] == "--seed":
            seed = int(args[1], 0)
            args = args[2:]
        elif args[0] == "--scale":
            scale = int(args[1])
            args = args[2:]
        elif args[0] == "--keep-raw":
            keep_raw = True
            args = args[1:]
        else:
            raise SystemExit("render.py: unknown argument %s" % args[0])
    os.makedirs(outdir, exist_ok=True)

    base = os.path.join(outdir, "world")
    subprocess.run(["./worldgen/worldgen", "--seed", str(seed), "--layer", "all",
                    "--threads", "4", "--out", base], check=True,
                   stdout=subprocess.DEVNULL)

    layers = {}
    layers["noise"] = read_f32(os.path.join(base, "noise.f32"))
    layers["erosion"] = read_f32(os.path.join(base, "erosion.f32"))
    layers["hydro"] = read_f32(os.path.join(base, "hydro.f32"))
    biome = read_u8(os.path.join(base, "biome.u8"))

    def gray(vals, lo=None, hi=None):
        g = shade(vals, lo, hi)
        return [(v, v, v) for v in g]

    write_ppm(os.path.join(outdir, "noise.ppm"), gray(layers["noise"]),
              WORLD, WORLD, scale)
    write_ppm(os.path.join(outdir, "erosion.ppm"), gray(layers["erosion"]),
              WORLD, WORLD, scale)

    # Hydro on a log scale: accumulation is very heavy-tailed.
    hydro_log = [math.log1p(v) for v in layers["hydro"]]
    write_ppm(os.path.join(outdir, "hydro.ppm"), gray(hydro_log),
              WORLD, WORLD, scale)

    write_ppm(os.path.join(outdir, "biome.ppm"),
              [BIOME_COLORS.get(v, (255, 0, 255)) for v in biome],
              WORLD, WORLD, scale)

    # Seam view: 4x4 tiles generated independently, stitched, then compared
    # against the full world.  Identical input must produce zero diff; we draw
    # the per-pixel difference amplified so any seam is obvious.
    tiles = os.path.join(outdir, "tiles")
    os.makedirs(tiles, exist_ok=True)
    for ty in range(4):
        for tx in range(4):
            subprocess.run(["./worldgen/worldgen", "--seed", str(seed),
                            "--layer", "all", "--threads", "4",
                            "--window", "%d,%d,64,64" % (tx * 64, ty * 64),
                            "--out", os.path.join(tiles, "%d_%d" % (tx, ty))],
                           check=True, stdout=subprocess.DEVNULL)

    for name, size in (("noise.f32", 4), ("erosion.f32", 4),
                       ("hydro.f32", 4), ("biome.u8", 1)):
        full = open(os.path.join(base, name), "rb").read()
        worst = 0
        for ty in range(4):
            for tx in range(4):
                tile = open(os.path.join(tiles, "%d_%d" % (tx, ty),
                                         name), "rb").read()
                for y in range(64):
                    a = full[((ty * 64 + y) * WORLD + tx * 64) * size:
                             ((ty * 64 + y) * WORLD + tx * 64 + 64) * size]
                    b = tile[y * 64 * size:(y + 1) * 64 * size]
                    if a != b:
                        worst = max(worst, sum(1 for p, q in zip(a, b) if p != q))
        print("seam %-12s mismatching bytes: %d" % (name, worst))

    have_ffmpeg = subprocess.run(["sh", "-c", "command -v ffmpeg"],
                                 stdout=subprocess.DEVNULL).returncode == 0
    if have_ffmpeg:
        for stem in ("noise", "erosion", "hydro", "biome"):
            subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-i",
                            os.path.join(outdir, stem + ".ppm"),
                            os.path.join(outdir, stem + ".png")], check=False)
            os.unlink(os.path.join(outdir, stem + ".ppm"))
        print("wrote PNGs to", outdir)
    else:
        print("ffmpeg not found; PPM files left in", outdir)

    # The raw world and tile scratch directories are only inputs to the seam
    # numbers printed above; keep the report directory small.
    if not keep_raw:
        shutil.rmtree(base, ignore_errors=True)
        shutil.rmtree(tiles, ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
