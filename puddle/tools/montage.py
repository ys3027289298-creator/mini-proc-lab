#!/usr/bin/env python3
"""Montage the puddle LOD PPM dumps into one comparison PNG.

Usage: python3 -I montage.py <ppm_dir> <out_png>
Rows: near / mid / far (same world region, same panel size).
Cols: high-res reference, existence (majority) path, coverage-blend path.
"""
import os
import sys

from PIL import Image, ImageDraw

PANEL = 256
PAD = 12
HEADER = 28
LABEL = 18
TIERS = ["near", "mid", "far"]
COLS = ["reference (512)", "existence / majority", "coverage blend"]


def load(ppm_dir, name):
    im = Image.open(os.path.join(ppm_dir, name))
    return im.resize((PANEL, PANEL), Image.NEAREST).convert("RGB")


def main():
    ppm_dir, out_png = sys.argv[1], sys.argv[2]
    w = PAD + 3 * (PANEL + PAD)
    h = HEADER + PAD + 3 * (PANEL + LABEL + PAD)
    canvas = Image.new("RGB", (w, h), (24, 24, 24))
    draw = ImageDraw.Draw(canvas)

    for c, title in enumerate(COLS):
        x = PAD + c * (PANEL + PAD)
        draw.text((x + 4, 8), title, fill=(230, 230, 230))

    for r, tier in enumerate(TIERS):
        y0 = HEADER + PAD + r * (PANEL + LABEL + PAD)
        panels = [
            load(ppm_dir, "reference.ppm"),
            load(ppm_dir, f"{tier}_hard.ppm"),
            load(ppm_dir, f"{tier}_blend.ppm"),
        ]
        for c, im in enumerate(panels):
            x = PAD + c * (PANEL + PAD)
            canvas.paste(im, (x, y0))
            draw.text((x + 4, y0 + PANEL + 2),
                      tier if c == 0 else "", fill=(230, 230, 230))

    canvas.save(out_png)
    print("wrote", out_png)


if __name__ == "__main__":
    main()
