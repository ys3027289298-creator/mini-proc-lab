#!/usr/bin/env python3
"""Aggregate vibfield 10-run output into OOB-cell repeat-rate statistics.

Reads the machine-readable records emitted by `vibfield runs <mode> 10`:
  RUN,run,seed,n_oob,n_pseudo
  OOB,run,r,c,side,dist,rho,freq,src_state,s1,s2,s3,s4,risk

Usage: python3 -I repeat_report.py <filtered.csv> <raw.csv> <out.md>
"""
import csv
import sys
from collections import Counter, defaultdict

RUNS = 10


def load(path):
    runs = {}
    oobs = defaultdict(list)
    with open(path, newline="") as fh:
        for row in csv.reader(fh):
            if not row:
                continue
            if row[0] == "RUN":
                runs[int(row[1])] = {
                    "seed": int(row[2]),
                    "n_oob": int(row[3]),
                    "n_pseudo": int(row[4]),
                }
            elif row[0] == "OOB":
                rec = {
                    "run": int(row[1]), "r": int(row[2]), "c": int(row[3]),
                    "side": row[4], "dist": int(row[5]),
                    "rho": float(row[6]), "freq": float(row[7]),
                    "src": row[8],
                    "s1": int(row[9]), "s2": int(row[10]),
                    "s3": int(row[11]), "s4": int(row[12]),
                    "risk": int(row[13]),
                }
                oobs[(rec["r"], rec["c"])].append(rec)
    return runs, oobs


def rate_bucket(k):
    if k == RUNS:
        return "10/10 (结构性，每轮必现)"
    if k >= 7:
        return "7-9/10 (高频，近结构)"
    if k >= 4:
        return "4-6/10 (中频)"
    if k >= 2:
        return "2-3/10 (低频，纹理驱动)"
    return "1/10 (偶发，纹理驱动)"


BUCKET_ORDER = [
    (lambda k: k == RUNS, "10/10 (结构性，每轮必现)"),
    (lambda k: 7 <= k < RUNS, "7-9/10 (高频，近结构)"),
    (lambda k: 4 <= k < 7, "4-6/10 (中频)"),
    (lambda k: 2 <= k < 4, "2-3/10 (低频，纹理驱动)"),
    (lambda k: k == 1, "1/10 (偶发，纹理驱动)"),
]


def summarize(name, runs, oobs):
    cells = [(rc, len(recs)) for rc, recs in oobs.items()]
    total_distinct = len(cells)
    total_visits = sum(k for _, k in cells)
    hist = Counter(k for _, k in cells)
    # repeat rate = average per-cell appearance probability over 10 runs
    repeat_rate = total_visits / (total_distinct * RUNS)
    risk_cells = {rc for rc, recs in oobs.items() if recs[0]["risk"]}
    risk_hist = Counter(len(oobs[rc]) for rc in risk_cells)
    side_hist = Counter(recs[0]["side"] for recs in oobs.values())
    freq_vals = [rec["freq"] for recs in oobs.values() for rec in recs]
    rho_vals = [rec["rho"] for recs in oobs.values() for rec in recs]
    return {
        "name": name, "runs": runs, "oobs": oobs,
        "total_distinct": total_distinct, "total_visits": total_visits,
        "hist": hist, "repeat_rate": repeat_rate,
        "risk_cells": risk_cells, "risk_hist": risk_hist,
        "side_hist": side_hist,
        "freq_min": min(freq_vals), "freq_max": max(freq_vals),
        "rho_min": min(rho_vals), "rho_max": max(rho_vals),
    }


def md_table_oob(s, limit=18):
    """Top OOB cells by repeat rate, then highest frequency."""
    rows = sorted(
        s["oobs"].items(),
        key=lambda kv: (-len(kv[1]), -kv[1][0]["freq"]),
    )
    lines = [
        "| 格子 (r,c) | 边界 | 外推距离 | 出现轮次 | 重复率 | 密度 ρ | 频率 Hz | 源状态 | 落石 |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for (r, c), recs in rows[:limit]:
        x = recs[0]
        lines.append(
            f"| ({r},{c}) | {x['side']} | {x['dist']} | {len(recs)}/{RUNS} | "
            f"{len(recs) / RUNS:.0%} | {x['rho']:.3f} | {x['freq']:.2f} | "
            f"{x['src']} | {'是' if x['risk'] else '否'} |"
        )
    return "\n".join(lines)


def main():
    f_runs, f_oobs = load(sys.argv[1])
    r_runs, r_oobs = load(sys.argv[2])
    fs = summarize("FILTERED（抑制伪震动）", f_runs, f_oobs)
    rs = summarize("RAW（不抑制，对照）", r_runs, r_oobs)
    out = sys.argv[3]

    L = []
    L.append("# 越界格子十轮重复率报告\n")
    L.append("固定地形（缓坡 × 斜向陡崖交叉），每轮换低海拔纹理种子（100–109），"
             "缓冲区宽 4 格，频率告警阈值 4.50 Hz。\n")

    L.append("## 1. 十轮总览\n")
    L.append("| 模式 | 每轮 OOB 数 | 抑制伪震动格/轮 | 不同 OOB 格总数 | 格-轮次总数 | 平均重复率 |")
    L.append("|---|---|---|---|---|---|")
    for s in (fs, rs):
        counts = [v["n_oob"] for v in s["runs"].values()]
        pseudo = [v["n_pseudo"] for v in s["runs"].values()]
        L.append(
            f"| {s['name']} | {min(counts)}–{max(counts)} | "
            f"{min(pseudo)}–{max(pseudo)} | {s['total_distinct']} | "
            f"{s['total_visits']} | {s['repeat_rate']:.1%} |"
        )
    L.append("")
    L.append("平均重复率 = 所有出现过的越界格在 10 轮中被观测到的平均概率；"
             "100% 表示每个越界格十轮全现。\n")

    L.append("## 2. 重复率分布（FILTERED）\n")
    L.append("| 重复档位 | 格子数 | 占比 |")
    L.append("|---|---|---|")
    for pred, label in BUCKET_ORDER:
        n = sum(cnt for k, cnt in fs["hist"].items() if pred(k))
        if n:
            L.append(f"| {label} | {n} | {n / fs['total_distinct']:.1%} |")
    L.append("")

    L.append("## 3. 落石格 vs 普通越界格的稳定性\n")
    L.append("| 类别 | 格子数 | 10/10 稳定 | 非全现 |")
    L.append("|---|---|---|---|")
    fr = fs["risk_cells"]
    fr_stable = sum(1 for rc in fr if len(fs["oobs"][rc]) == RUNS)
    non = set(fs["oobs"]) - fr
    non_stable = sum(1 for rc in non if len(fs["oobs"][rc]) == RUNS)
    L.append(f"| 落石格（rockfall） | {len(fr)} | {fr_stable} | {len(fr) - fr_stable} |")
    L.append(f"| 普通越界格 | {len(non)} | {non_stable} | {len(non) - non_stable} |")
    L.append("")
    L.append("落石检查的前三道门基于无纹理基准地形，因此落石格与纹理种子无关、十轮全现；"
             "普通越界格由带纹理的密度场外推产生，随种子漂移。\n")

    L.append("## 4. 重复率最高的越界格（FILTERED，前 18）\n")
    L.append(md_table_oob(fs))
    L.append("")

    L.append("## 5. 边界与频率范围\n")
    L.append("| 模式 | N | S | W | E | 频率 min–max Hz | 密度 min–max |")
    L.append("|---|---|---|---|---|---|---|")
    for s in (fs, rs):
        h = s["side_hist"]
        L.append(
            f"| {s['name']} | {h.get('N', 0)} | {h.get('S', 0)} | "
            f"{h.get('W', 0)} | {h.get('E', 0)} | "
            f"{s['freq_min']:.2f}–{s['freq_max']:.2f} | "
            f"{s['rho_min']:.3f}–{s['rho_max']:.3f} |"
        )
    L.append("")

    L.append("## 6. FILTERED 与 RAW 对照\n")
    common = set(fs["oobs"]) & set(rs["oobs"])
    only_raw = set(rs["oobs"]) - set(fs["oobs"])
    only_f = set(fs["oobs"]) - set(rs["oobs"])
    L.append(f"- 两种模式共同越界格：{len(common)}")
    L.append(f"- 仅 RAW 越界（伪震动未抑制、外推出的假警报）：{len(only_raw)}")
    L.append(f"- 仅 FILTERED 越界：{len(only_f)}")
    L.append(f"- RAW 平均重复率 {rs['repeat_rate']:.1%} vs "
             f"FILTERED {fs['repeat_rate']:.1%}")
    L.append("")

    with open(out, "w") as fh:
        fh.write("\n".join(L) + "\n")

    # full per-cell detail table alongside the markdown report
    detail = out.rsplit(".", 1)[0] + "_cells.csv"
    with open(detail, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["r", "c", "side", "dist", "appearances", "repeat_rate",
                    "rho", "freq_hz", "src_state", "s1_cliff", "s2_outward",
                    "s3_channel", "s4_alarm", "rockfall"])
        for (r, c), recs in sorted(fs["oobs"].items(),
                                   key=lambda kv: (-len(kv[1]), kv[0])):
            x = recs[0]
            w.writerow([r, c, x["side"], x["dist"], len(recs),
                        f"{len(recs) / RUNS:.2f}", f"{x['rho']:.4f}",
                        f"{x['freq']:.4f}", x["src"], x["s1"], x["s2"],
                        x["s3"], x["s4"], x["risk"]])

    print(f"wrote {out}")
    print(f"wrote {detail}")
    print(f"FILTERED distinct={fs['total_distinct']} repeat={fs['repeat_rate']:.3f} "
          f"rockfall={len(fr)} stable={fr_stable}")
    print(f"RAW distinct={rs['total_distinct']} repeat={rs['repeat_rate']:.3f}")


if __name__ == "__main__":
    main()
