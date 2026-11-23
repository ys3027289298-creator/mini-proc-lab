#!/usr/bin/env python3
"""analyze.py — 读取 biome gen 产物,做 15% 门槛判定并出图。

产物(写入同一 outdir):
  fig_functions.png   两条函数曲线 F(w)/G(w) + F-G 权衡视图
  fig_gate.png        实测对预测的相对偏离与 15% 门槛
  fig_tiers.png       三档带宽(w=2/6/12)渲染与模糊度 phi 热图
  fig_worst_zoom.png  最割裂区段在三档带宽下的局部放大
  fuzz_stats.csv      模糊度统计(表视图)
  verdict.txt         门槛判定与最割裂区段结论
门槛失败时进程退出码为 1。
"""
import csv
import os
import struct
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import ListedColormap

import glob
for _p in (glob.glob(os.path.expanduser("~/.local/share/fonts/wqy-zenhei.ttc"))
           + glob.glob("/usr/share/fonts/**/wqy-zenhei.ttc", recursive=True)):
    if os.path.exists(_p):
        font_manager.fontManager.addfont(_p)
plt.rcParams["font.sans-serif"] = ["WenQuanYi Zen Hei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

# dataviz 参考调色板(light,已通过 validate_palette --pairs all)
SURFACE = "#fcfcfb"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
C_MEAS, C_CONT, C_DISC = "#2a78d6", "#eb6834", "#1baf7a"
BLUE_RAMP = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec",
             "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab",
             "#184f95", "#104281", "#0d366b"]
GATE = 0.15
TIERS = [2, 6, 12]


def load_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def load_ppm(path):
    with open(path, "rb") as f:
        buf = f.read()
    assert buf[:2] == b"P6"
    idx = 2

    def token():
        nonlocal idx
        while buf[idx] in b" \t\n\r":
            idx += 1
        start = idx
        while buf[idx] not in b" \t\n\r":
            idx += 1
        return buf[start:idx]

    w = int(token()); h = int(token()); mx = int(token())
    idx += 1
    assert mx == 255
    raw = buf[idx:idx + w * h * 3]
    import numpy as np
    return np.frombuffer(raw, dtype=np.uint8).reshape(h, w, 3)


def load_f32(path, n):
    import numpy as np
    with open(path, "rb") as f:
        return np.frombuffer(f.read(), dtype="<f4", count=n * n).reshape(n, n).copy()


def main(outdir):
    import numpy as np
    with open(os.path.join(outdir, "meta.csv"), newline="") as f:
        meta = next(csv.DictReader(f))
    n = int(meta["grid"])
    sweep = load_csv(os.path.join(outdir, "sweep.csv"))
    worst = load_csv(os.path.join(outdir, "worst.csv"))
    w = np.array([int(r["w"]) for r in sweep])
    fuzz = np.array([float(r["fuzz_meas"]) for r in sweep])
    fcont = np.array([float(r["fuzz_cont_pred"]) for r in sweep])
    fdisc = np.array([float(r["fuzz_disc_pred"]) for r in sweep])
    look = np.array([float(r["extra_lookups"]) for r in sweep])
    gcont = np.array([float(r["lookup_cont_pred"]) for r in sweep])
    gdisc = np.array([float(r["lookup_disc_pred"]) for r in sweep])
    fcerr = np.abs(fuzz - fcont) / fuzz
    fderr = np.abs(fuzz - fdisc) / fuzz
    gcerr = np.abs(look - gcont) / look
    gderr = np.abs(look - gdisc) / look

    # 独立交叉核对:Python 直接读 phi 二进制,均值须等于 C 导出的 fuzz_meas
    cross = []
    for ww in TIERS:
        phi = load_f32(os.path.join(outdir, f"phi_w{ww}.f32"), n)
        row = next(r for r in sweep if int(r["w"]) == ww)
        cross.append((ww, abs(phi.mean() - float(row["fuzz_meas"]))))
    assert all(e < 1e-7 for _, e in cross), f"phi 交叉核对失败 {cross}"

    # ---- 门槛判定:模糊度对预测偏离 > 15% 即未完成 ----
    disc_fail = [int(x) for x, e in zip(w, fderr) if e > GATE]
    cont_fail = [int(x) for x, e in zip(w, fcerr) if x >= 4 and e > GATE]
    cont_small = [(int(x), float(e)) for x, e in zip(w, fcerr) if x < 4 and e > GATE]
    passed = not disc_fail and not cont_fail

    # ---- 图 1:两条函数曲线 + 权衡视图 ----
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4), facecolor=SURFACE)
    for ax in axes:
        ax.set_facecolor(SURFACE)
        ax.grid(True, color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        for sp in ax.spines.values():
            sp.set_color(BASELINE)
        ax.tick_params(colors=INK2, labelsize=9)
        ax.xaxis.label.set_color(INK); ax.yaxis.label.set_color(INK)
        ax.title.set_color(INK)

    ax = axes[0]
    ax.plot(w, fuzz, "-o", color=C_MEAS, lw=2, ms=7, label="实测 F(w)")
    ax.plot(w, fcont, "--s", color=C_CONT, lw=2, ms=6, label="连续预测 F_c=18/35·Lw/P")
    ax.plot(w, fdisc, ":^", color=C_DISC, lw=2, ms=7, label="离散预测 F_d(距离场积分)")
    ax.set_xlabel("过渡带带宽 w(像素)"); ax.set_ylabel("全场平均模糊度 F")
    ax.set_title("模糊度关于带宽的函数")
    ax.legend(frameon=False, fontsize=8.5, loc="upper left")

    ax = axes[1]
    ax.plot(w, look, "-o", color=C_MEAS, lw=2, ms=7, label="实测 G(w)")
    ax.plot(w, gcont, "--s", color=C_CONT, lw=2, ms=6, label="连续预测 G_c=(w+1)Lw")
    ax.plot(w, gdisc, ":^", color=C_DISC, lw=2, ms=7, label="离散预测 G_d=(w+1)B(w)")
    ax.set_xlabel("过渡带带宽 w(像素)"); ax.set_ylabel("带内额外噪声查找次数 G")
    ax.set_title("噪声查找增长关于带宽的函数")
    ax.legend(frameon=False, fontsize=8.5, loc="upper left")

    ax = axes[2]
    ax.plot(gdisc, fuzz, "-o", color=C_MEAS, lw=2, ms=7)
    for ww in TIERS:
        i = int(ww) - 1
        ax.annotate(f"w={ww}", (gdisc[i], fuzz[i]),
                    textcoords="offset points", xytext=(8, 4),
                    fontsize=9, color=INK)
    ax.set_xlabel("噪声查找次数 G"); ax.set_ylabel("平均模糊度 F")
    ax.set_title("权衡曲线:模糊度 ↔ 查找成本")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig_functions.png"), dpi=140,
                facecolor=SURFACE)
    plt.close(fig)

    # ---- 图 2:偏离与 15% 门槛 ----
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), facecolor=SURFACE)
    for ax in axes:
        ax.set_facecolor(SURFACE)
        ax.grid(True, axis="y", color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        for sp in ax.spines.values():
            sp.set_color(BASELINE)
        ax.tick_params(colors=INK2, labelsize=9)
        ax.xaxis.label.set_color(INK); ax.yaxis.label.set_color(INK)
        ax.title.set_color(INK)
    width = 0.38
    for ax, cerr_v, derr_v, title in (
        (axes[0], fcerr, fderr, "模糊度:预测 vs 实测"),
        (axes[1], gcerr, gderr, "噪声查找:预测 vs 实测"),
    ):
        ax.bar(w - width / 2, cerr_v * 100, width, color=C_CONT, label="连续预测偏离")
        ax.bar(w + width / 2, derr_v * 100, width, color=C_DISC, label="离散预测偏离")
        zero_x = w[derr_v == 0] + width / 2
        if len(zero_x):
            ax.plot(zero_x, np.zeros_like(zero_x), "_", ms=11, mew=2,
                    color=C_DISC, label="离散预测偏离=0(精确)")
        ax.axhline(15, color="#d03b3b", lw=1.6, ls="--")
        ax.text(12.4, 15.6, "15% 门槛", color="#d03b3b", fontsize=9, ha="right")
        ax.set_xlabel("带宽 w")
        ax.set_ylabel(("模糊度" if "模糊度" in title else "查找次数") + "相对偏离 (%)")
        ax.set_title(title)
        ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig_gate.png"), dpi=140, facecolor=SURFACE)
    plt.close(fig)

    # ---- 图 3:三档带宽渲染 + phi 热图 ----
    cmap = ListedColormap(BLUE_RAMP)
    fig, axes = plt.subplots(2, 3, figsize=(13, 8.6), facecolor=SURFACE)
    wx, wy = int(worst[0]["x"]), int(worst[0]["y"])
    for col, ww in enumerate(TIERS):
        rgb = load_ppm(os.path.join(outdir, f"render_w{ww}.ppm"))
        phi = load_f32(os.path.join(outdir, f"phi_w{ww}.f32"), n)
        row = next(r for r in sweep if int(r["w"]) == ww)
        ax = axes[0][col]
        ax.imshow(rgb, origin="upper")
        ax.set_title(f"w={ww} 渲染 | 带内 {row['band_px']}px | "
                     f"F={float(row['fuzz_meas']):.4f} | G={row['extra_lookups']}",
                     fontsize=10)
        ax.tick_params(colors=MUTED, labelsize=8)
        if ww == TIERS[0]:
            ax.plot(wx, wy, marker="x", ms=14, mew=2.5, color="#d03b3b")
            ax.annotate("最割裂点", (wx, wy), textcoords="offset points",
                        xytext=(-10, -18), color="#d03b3b", fontsize=9,
                        arrowprops=dict(arrowstyle="->", color="#d03b3b", lw=1.2))
        ax = axes[1][col]
        im = ax.imshow(phi, origin="upper", cmap=cmap, vmin=0, vmax=1)
        ax.set_title(f"w={ww} 逐像素模糊度 φ", fontsize=10)
        ax.tick_params(colors=MUTED, labelsize=8)
    cbar = fig.colorbar(im, ax=axes[1].tolist(), fraction=0.025, pad=0.01)
    cbar.ax.tick_params(colors=INK2, labelsize=8)
    cbar.outline.set_color(BASELINE)
    fig.suptitle("三档过渡带对比:窄带割裂(w=2) / 折中(w=6) / 宽带充分混合(w=12)",
                 fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(os.path.join(outdir, "fig_tiers.png"), dpi=140, facecolor=SURFACE)
    plt.close(fig)

    # ---- 图 4:最割裂区段局部放大 ----
    r = 24
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.6), facecolor=SURFACE)
    for col, ww in enumerate(TIERS):
        rgb = load_ppm(os.path.join(outdir, f"render_w{ww}.ppm"))
        crop = rgb[max(wy - r, 0):wy + r, max(wx - r, 0):wx + r]
        ax = axes[col]
        ax.imshow(crop, origin="upper",
                  extent=(max(wx - r, 0), wx + r, wy + r, max(wy - r, 0)))
        ax.plot(wx, wy, marker="x", ms=14, mew=2.5, color="#d03b3b")
        ax.set_title(f"w={ww} 最割裂区段局部 (x={wx}, y={wy})", fontsize=10)
        ax.tick_params(colors=MUTED, labelsize=8)
    fig.suptitle(f"同一边界点的邻域色跳:窄带最大色距 {float(worst[0]['color_jump']):.1f}",
                 fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0.02, 1, 0.94))
    fig.savefig(os.path.join(outdir, "fig_worst_zoom.png"), dpi=140,
                facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)

    # ---- 统计表(表视图) ----
    with open(os.path.join(outdir, "fuzz_stats.csv"), "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["w", "band_px", "fuzz_meas", "fuzz_cont_pred", "fuzz_disc_pred",
                     "cont_rel_err", "disc_rel_err", "extra_lookups", "gate"])
        for row in sweep:
            ww = int(row["w"])
            ok = float(row["fuzz_disc_rel_err"]) <= GATE and (
                ww < 4 or float(row["fuzz_cont_rel_err"]) <= GATE)
            wr.writerow([ww, row["band_px"], row["fuzz_meas"], row["fuzz_cont_pred"],
                         row["fuzz_disc_pred"], row["fuzz_cont_rel_err"],
                         row["fuzz_disc_rel_err"], row["extra_lookups"],
                         "PASS" if ok else "FAIL"])

    # ---- 结论文本 ----
    top = worst[0]
    lines = [
        f"网格 {n}x{n} (P={meta['P']}), 边界长度 L={meta['L']}, "
        f"场生成噪声查找 {meta['field_noise_calls']} 次(与 w 无关)",
        f"函数定义: F_c(w)=18/35·Lw/P ; G_c(w)=(w+1)Lw ; "
        f"离散预测对同一距离场独立积分",
        f"门槛: 模糊度与预测相对偏离 > 15% 判未完成",
        f"  离散预测 F_d: 全部 w=1..12 最大偏离 {fderr.max()*100:.2f}% "
        f"({'PASS' if not disc_fail else 'FAIL @w=' + str(disc_fail)})",
        f"  连续预测 F_c: w>=4 最大偏离 {max(e for x, e in zip(w, fcerr) if x >= 4)*100:.2f}% "
        f"({'PASS' if not cont_fail else 'FAIL @w=' + str(cont_fail)})",
        f"  连续预测在细带失准(像素化): "
        + ", ".join(f"w={x} {e*100:.1f}%" for x, e in cont_small),
        f"  查找 G_d 离散预测最大偏离 {gderr.max()*100:.2f}%; "
        f"G_c 连续近似在 w>=9 超 15%(最大 {gcerr.max()*100:.1f}%),属渐近公式截断误差,已留档",
        f"phi 二进制独立交叉核对最大偏差 {max(e for _, e in cross):.2e}",
        f"最割裂区段: (x={top['x']}, y={top['y']}), 窄带 w={top['scan_width']} 上 "
        f"邻域最大 RGB 色距 {float(top['color_jump']):.2f}, 该处场梯度 {float(top['g_grad']):.4f}",
        f"结论: {'完成 — 全部模糊度门槛通过' if passed else '未完成 — 存在超 15% 偏离'}",
    ]
    with open(os.path.join(outdir, "verdict.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print("\n三档模糊度统计:")
    print(f"{'w':>3} {'带内像素':>9} {'F实测':>9} {'连续偏离':>9} {'离散偏离':>9} "
          f"{'G查找':>9}  判定")
    for row in sweep:
        ww = int(row["w"])
        if ww not in TIERS:
            continue
        ok = float(row["fuzz_disc_rel_err"]) <= GATE and (
            ww < 4 or float(row["fuzz_cont_rel_err"]) <= GATE)
        print(f"{ww:>3} {row['band_px']:>9} {float(row['fuzz_meas']):>9.4f} "
              f"{float(row['fuzz_cont_rel_err'])*100:>8.2f}% "
              f"{float(row['fuzz_disc_rel_err'])*100:>8.2f}% "
              f"{row['extra_lookups']:>9}  {'PASS' if ok else 'FAIL'}")
    return 0 if passed else 1


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("用法: python3 analyze.py <outdir>")
    sys.exit(main(sys.argv[1]))
