#!/usr/bin/env python3
"""Assemble a 3 speeds x 3 directions montage of PGM textures into one PNG.

Run under: python3 -I montage.py MODE
where MODE is the output subdirectory name (fixed|coupled).
"""
import os
import sys
import site

_usp = site.getusersitepackages()
if _usp not in sys.path:
    sys.path.append(_usp)

from PIL import Image, ImageDraw

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

SPEEDS = [4, 7, 10]
DIRS = [0, 45, 90]
CELL = 256
MARGIN = 64
LABEL = 22


def main():
    mode = sys.argv[1]
    d = os.path.join(ROOT, "out", mode)
    canvas = Image.new("L", (MARGIN + 3 * CELL + 2 * LABEL,
                             MARGIN + 3 * CELL + 2 * LABEL), 255)
    draw = ImageDraw.Draw(canvas)
    draw.text((8, 4), "mode=%s (rows: wind speed m/s, cols: direction deg)" % mode, fill=0)

    for j, ang in enumerate(DIRS):
        x = MARGIN + j * (CELL + LABEL)
        draw.text((x + CELL // 2 - 14, 26), "dir=%d" % ang, fill=0)
    for i, spd in enumerate(SPEEDS):
        y = MARGIN + i * (CELL + LABEL)
        draw.text((6, y + CELL // 2), "v=%d" % spd, fill=0)

    for i, spd in enumerate(SPEEDS):
        for j, ang in enumerate(DIRS):
            img = Image.open(os.path.join(d, "s%02d_d%03d.pgm" % (spd, ang)))
            x = MARGIN + j * (CELL + LABEL)
            y = MARGIN + i * (CELL + LABEL)
            canvas.paste(img, (x, y))
            draw.rectangle([x - 1, y - 1, x + CELL, y + CELL], outline=128)

    out = os.path.join(d, "montage.png")
    canvas.save(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
