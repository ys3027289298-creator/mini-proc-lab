#!/usr/bin/env python3
# 河谷五段地表色过渡：宽度公式 + 敏感度公式 -> 定宽 -> 仿真校验 -> 修正 -> 出剖面 SVG
# 纯 stdlib。
#
# 式1(地形驱动):  w_terrain = W0 + ALPHA*dH + BETA*dS
# 式2(颜色敏感度): 串色 B = dC*w/(6L) <= RHO_B*delta  ->  w <= w_hi = 6*RHO_B*L*delta/dC
#                 跳变 J = 1.5*dC/w <= RHO_J*delta    ->  w >= w_lo = 1.5*dC/(RHO_J*delta)
# 两式互制: w = clamp(w_terrain, w_lo, w_hi); 敏感度量 S = dB/dw = dC/(6L), 段越短、邻段色差越大越敏感。
import math, os

SEGS = [
    dict(name="河道",   x0=0,   x1=20,  H=1.0,  s=0.100, C=( 88,120,80), d_est=12.0, amp=3.5),
    dict(name="河漫滩", x0=20,  x1=80,  H=4.0,  s=0.067, C=(110,150,70), d_est=18.0, amp=2.6),
    dict(name="阶地",   x0=80,  x1=110, H=12.0, s=0.400, C=(150,160,90), d_est=25.0, amp=7.2),
    dict(name="坡地",   x0=110, x1=200, H=39.0, s=0.467, C=(120,140,95), d_est=30.0, amp=8.7),
    dict(name="山脊",   x0=200, x1=260, H=67.5, s=0.250, C=(110,115,85), d_est=15.0, amp=4.3),
]
W0, ALPHA, BETA = 2.0, 0.5, 15.0
RHO_B, RHO_J = 0.4, 0.5
DX = 0.5

def smoothstep(t):
    t = min(1.0, max(0.0, t))
    return t * t * (3 - 2 * t)

def cdist(a, b):
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))

def seg_len(k): return SEGS[k]["x1"] - SEGS[k]["x0"]

def compute_bounds(dkey):
    """按指定的段内色差字段(d_est 或 d_meas)推各边界宽度。"""
    out = []
    for i in range(len(SEGS) - 1):
        a, b = SEGS[i], SEGS[i + 1]
        dH = abs(b["H"] - a["H"]); dS = abs(b["s"] - a["s"]); dC = cdist(a["C"], b["C"])
        w_t = W0 + ALPHA * dH + BETA * dS
        w_lo = 1.5 * dC / (RHO_J * min(a[dkey], b[dkey]))
        w_hi = min(6 * RHO_B * seg_len(i) * a[dkey] / dC,
                   6 * RHO_B * seg_len(i + 1) * b[dkey] / dC)
        out.append(dict(i=i, x=a["x1"], dH=dH, dS=dS, dC=dC,
                        w_terrain=w_t, w_lo=w_lo, w_hi=w_hi,
                        w=min(max(w_t, w_lo), w_hi)))
    return out

def jitter(k, x):
    ph = k * 1.7 + 0.9
    return SEGS[k]["amp"] * (0.6 * math.sin(0.9 * x + ph) + 0.4 * math.sin(2.3 * x + 1.7 * ph))

def field(k, x):
    j = jitter(k, x)
    return tuple(SEGS[k]["C"][ch] + j for ch in range(3))

def weights(bounds, x):
    ws = []
    for k, sg in enumerate(SEGS):
        w = 1.0
        if k > 0:
            b = bounds[k - 1]
            w *= smoothstep((x - (b["x"] - b["w"] / 2)) / b["w"])
        if k < len(SEGS) - 1:
            b = bounds[k]
            w *= 1.0 - smoothstep((x - (b["x"] - b["w"] / 2)) / b["w"])
        ws.append(w)
    s = sum(ws)
    return [v / s for v in ws]

def color_after(bounds, x, base_only=False):
    ws = weights(bounds, x)
    c = [0.0, 0.0, 0.0]
    for k in range(len(SEGS)):
        f = SEGS[k]["C"] if base_only else field(k, x)
        for ch in range(3):
            c[ch] += ws[k] * f[ch]
    return tuple(c)

def color_before(x):
    for k, sg in enumerate(SEGS):
        if sg["x0"] <= x < sg["x1"] or (x == SEGS[-1]["x1"] and k == len(SEGS) - 1):
            return field(k, x)
    return field(len(SEGS) - 1, x)

# ---- 实测段内色差(对扰动采样, 标量扰动 x√3 -> RGB 距离) ----
xs_all = [i * DX for i in range(int(SEGS[-1]["x1"] / DX) + 1)]
for k, sg in enumerate(SEGS):
    js = [jitter(k, x) for x in xs_all if sg["x0"] + 1 <= x <= sg["x1"] - 1]
    sg["d_meas"] = math.sqrt(3) * (max(js) - min(js))

def run_pass(bounds, tag):
    """对每条边界实测跳变(基色混合趋势, 每米), 与弱侧段内色差实测值求比。"""
    print(f"\n== {tag} ==")
    print(f"{'边界':<13}{'dH':>5}{'dS':>7}{'dC':>6}{'w地形':>7}{'w下界':>7}{'w上界':>7}{'w定':>7}"
          f"{'J/m':>7}{'δ弱侧':>7}{'比值':>6}  判定")
    all_ok = True
    for b in bounds:
        xs = [b["x"] - b["w"] / 2 + i * DX for i in range(int(b["w"] / DX) + 1)]
        J = max(cdist(color_after(bounds, a, True), color_after(bounds, c, True)) / DX
                for a, c in zip(xs, xs[1:]))
        d_weak = min(SEGS[b["i"]]["d_meas"], SEGS[b["i"] + 1]["d_meas"])
        r = J / d_weak
        ok = r <= RHO_J
        all_ok &= ok
        nm = SEGS[b["i"]]["name"] + "|" + SEGS[b["i"] + 1]["name"]
        print(f"{nm:<13}{b['dH']:>5.1f}{b['dS']:>7.3f}{b['dC']:>6.1f}{b['w_terrain']:>7.1f}"
              f"{b['w_lo']:>7.1f}{b['w_hi']:>7.1f}{b['w']:>7.1f}{J:>7.2f}{d_weak:>7.1f}{r:>6.2f}"
              f"  {'OK' if ok else '超差'}")
    return all_ok

print("段内色差 估值 vs 实测:")
for sg in SEGS:
    flag = "  <-- 估值偏差大" if abs(sg["d_est"] - sg["d_meas"]) / sg["d_meas"] > 0.3 else ""
    print(f"  {sg['name']:<5} est={sg['d_est']:>5.1f}  meas={sg['d_meas']:>5.1f}{flag}")

b1 = compute_bounds("d_est")
ok1 = run_pass(b1, "第一轮: 按估值定宽")

if not ok1:
    # 诊断: 超差边界共邻的、估值偏差大的段 -> 用实测值修正其敏感度, 重定宽
    bad = [k for k, sg in enumerate(SEGS)
           if abs(sg["d_est"] - sg["d_meas"]) / sg["d_meas"] > 0.3]
    print("\n诊断: 段内色差估值失准的段 ->", "、".join(SEGS[k]["name"] for k in bad))
    for k in bad:
        SEGS[k]["d_est"] = SEGS[k]["d_meas"]
    b2 = compute_bounds("d_est")
    print("修正后宽度:", {SEGS[b['i']]['name'] + '|' + SEGS[b['i']+1]['name']: round(b['w'], 1) for b in b2})
    run_pass(b2, "第二轮: 修正敏感度后定宽")
    bounds = b2
else:
    bounds = b1

# ---------- 五段剖面 SVG(改前/改后各一条, 用最终宽度) ----------
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
os.makedirs(OUT, exist_ok=True)

def lum(c): return 0.299 * c[0] + 0.587 * c[1] + 0.114 * c[2]
def rgb(c): return "#%02x%02x%02x" % tuple(min(255, max(0, int(round(v)))) for v in c)

PW, PH = 860, 460
ML, MR, MT, MB = 60, 20, 44, 96
for k, sg in enumerate(SEGS):
    wl = bounds[k - 1]["w"] if k > 0 else 0.0
    wr = bounds[k]["w"] if k < len(bounds) else 0.0
    margin = max(14.0, 1.3 * max(wl, wr))
    xa = max(0.0, sg["x0"] - margin); xb = min(SEGS[-1]["x1"], sg["x1"] + margin)
    xs = [xa + i * DX for i in range(int((xb - xa) / DX) + 1)]
    yb = [lum(color_before(x)) for x in xs]
    ya = [lum(color_after(bounds, x)) for x in xs]
    y0, y1 = min(yb + ya) - 6, max(yb + ya) + 6
    X = lambda x: ML + (x - xa) / (xb - xa) * (PW - ML - MR)
    Y = lambda y: MT + (y1 - y) / (y1 - y0) * (PH - MT - MB)
    poly = lambda ys: " ".join(f"{X(x):.1f},{Y(y):.1f}" for x, y in zip(xs, ys))

    p = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{PW}" height="{PH}"'
         f' font-family="monospace" style="background:#fafaf7">']
    p.append(f'<text x="{ML}" y="26" font-size="16">段{k+1} {sg["name"]}  x∈[{sg["x0"]},{sg["x1"]}]m'
             f'  过渡宽 左{wl:.1f}m / 右{wr:.1f}m</text>')
    for b in bounds:
        if xa <= b["x"] <= xb:
            p.append(f'<rect x="{X(b["x"]-b["w"]/2):.1f}" y="{MT}"'
                     f' width="{X(b["x"]+b["w"]/2)-X(b["x"]-b["w"]/2):.1f}" height="{PH-MT-MB}"'
                     f' fill="#000" opacity="0.06"/>')
            p.append(f'<line x1="{X(b["x"]):.1f}" y1="{MT}" x2="{X(b["x"]):.1f}" y2="{PH-MB}"'
                     f' stroke="#888" stroke-dasharray="3 3"/>')
    p.append(f'<line x1="{ML}" y1="{PH-MB}" x2="{PW-MR}" y2="{PH-MB}" stroke="#333"/>')
    p.append(f'<line x1="{ML}" y1="{MT}" x2="{ML}" y2="{PH-MB}" stroke="#333"/>')
    for gx in range(int(xa) - int(xa) % 20, int(xb) + 1, 20):
        p.append(f'<line x1="{X(gx):.1f}" y1="{PH-MB}" x2="{X(gx):.1f}" y2="{PH-MB+5}" stroke="#333"/>')
        p.append(f'<text x="{X(gx):.1f}" y="{PH-MB+18}" font-size="11" text-anchor="middle">{gx}m</text>')
    p.append(f'<text x="16" y="{(MT+PH-MB)/2:.0f}" font-size="11"'
             f' transform="rotate(-90 16 {(MT+PH-MB)/2:.0f})">亮度</text>')
    p.append(f'<polyline points="{poly(yb)}" fill="none" stroke="#b04a3a"'
             f' stroke-width="1.6" stroke-dasharray="5 3"/>')
    p.append(f'<polyline points="{poly(ya)}" fill="none" stroke="#2e6b3e" stroke-width="2.2"/>')
    for row, fn, lab in ((0, color_before, "改前"), (1, lambda x: color_after(bounds, x), "改后")):
        y = PH - MB + 30 + row * 22
        p.append(f'<text x="{ML-34}" y="{y+13}" font-size="11">{lab}</text>')
        for x in xs:
            p.append(f'<rect x="{X(x):.1f}" y="{y}" width="{X(x+DX)-X(x)+0.6:.1f}" height="16"'
                     f' fill="{rgb(fn(x))}"/>')
    p.append(f'<text x="{ML}" y="{PH-6}" font-size="11" fill="#b04a3a">— — 改前(硬切)</text>')
    p.append(f'<text x="{ML+150}" y="{PH-6}" font-size="11" fill="#2e6b3e">—— 改后(式1/式2定宽)</text>')
    p.append('</svg>')
    fp = os.path.join(OUT, f"seg{k+1}_{sg['name']}.svg")
    with open(fp, "w") as f:
        f.write("\n".join(p))
    print("写出", fp)
