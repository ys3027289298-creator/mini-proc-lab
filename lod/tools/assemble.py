import os
import sys

from PIL import Image, ImageDraw

def load(d, n):
    return Image.open(os.path.join(d, n)).convert("RGB")

def frame_diff(a, b, thresh=24, gain=4):
    out = Image.new("RGB", a.size)
    pa, pb, po = a.load(), b.load(), out.load()
    w, h = a.size
    cnt = 0
    for j in range(h):
        for i in range(w):
            d = sum(abs(pa[i, j][k] - pb[i, j][k]) for k in range(3))
            if d > thresh:
                cnt += 1
            v = min(255, d * gain)
            po[i, j] = (v, v // 3, v // 3)
    return out, cnt

def side_by_side(left, right, label_l, label_r, path, gap=12, header=26):
    w, h = left.size
    canvas = Image.new("RGB", (w * 2 + gap, h + header), (18, 18, 18))
    canvas.paste(left, (0, header))
    canvas.paste(right, (w + gap, header))
    dr = ImageDraw.Draw(canvas)
    dr.text((8, 7), label_l, fill=(255, 200, 120))
    dr.text((w + gap + 8, 7), label_r, fill=(140, 255, 160))
    canvas.save(path)

def main():
    d = sys.argv[1]
    for t in ("tier_near", "tier_mid", "tier_far"):
        load(d, t + ".ppm").save(os.path.join(d, t + ".png"))

    nb, nt = load(d, "near_base.ppm"), load(d, "near_tuned.ppm")
    nd, nc = frame_diff(nb, nt)
    w, h = nb.size
    canvas = Image.new("RGB", (w * 3 + 24, h + 26), (18, 18, 18))
    canvas.paste(nb, (0, 26))
    canvas.paste(nt, (w + 12, 26))
    canvas.paste(nd, (2 * w + 24, 26))
    dr = ImageDraw.Draw(canvas)
    dr.text((8, 7), "baseline tau=4.0 -> L1 (near too coarse)", fill=(255, 200, 120))
    dr.text((w + 20, 7), "tuned tau=1.2 -> L0", fill=(140, 255, 160))
    dr.text((2 * w + 32, 7), "diff x4 (%d px)" % nc, fill=(200, 200, 200))
    canvas.save(os.path.join(d, "compare_near.png"))
    print("near_diff_px=%d" % nc)

    db, cb = frame_diff(load(d, "far_b0.ppm"), load(d, "far_b1.ppm"))
    dt, ct = frame_diff(load(d, "far_t0.ppm"), load(d, "far_t1.ppm"))
    side_by_side(
        db, dt,
        "baseline frame-diff x4 (%d px)" % cb, "tuned frame-diff x4 (%d px)" % ct,
        os.path.join(d, "compare_far.png"),
    )
    print("flicker_px_base=%d flicker_px_tuned=%d" % (cb, ct))

    for f in os.listdir(d):
        if f.endswith(".ppm"):
            os.remove(os.path.join(d, f))

if __name__ == "__main__":
    main()
