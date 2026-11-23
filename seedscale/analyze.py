#!/usr/bin/env python3
"""
analyze.py — 读取 fieldstats 的 CSV，产出三档特征统计与独立性检验报告。

用法: python3 -I analyze.py out [report.txt]

检验分三层:
  H1 跨种子独立性: 相邻种子场互相关 xcorr 必须 ~0（t 检验，Bonferroni 校正）
  H2 尺度独立性(跨档): 相关长度、方差等特征在三档种子间无显著差异
                    （Kruskal-Wallis + 跨档均值漂移的效应量）
  H3 尺度独立性(档内): 按种子低 3 位分 8 组，相关长度/高频比组间无显著差异。
                    H2 只能抓"档位差异"，H3 抓与种子位系统性耦合的尺度变化
                    （sx 模式的格距乘数只依赖种子低 3 位，三档分布会被平均掉）
另外 good_nn 对照演示采样端混叠：H1/H2/H3 全过，但高频能量比 hf 异常升高。
"""
import csv
import math
import os
import sys

LAGP = 33
PPC = 8
TIERS = ["S", "M", "L"]
TIER_LABEL = {"S": "S_small(0..127)", "M": "M_mid(2^16+)", "L": "L_full32"}


def load(path):
    rows = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            for k, v in list(r.items()):
                if k in ("tier", "seed", "octaves"):
                    r[k] = int(v) if k != "tier" else v
                else:
                    r[k] = float(v)
            rows.append(r)
    return rows


def mean(xs):
    return sum(xs) / len(xs)


def sd(xs):
    m = mean(xs)
    return math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1))


def tstat_onesample(xs, mu0=0.0):
    n = len(xs)
    m, s = mean(xs), sd(xs)
    return (m - mu0) / (s / math.sqrt(n))


# 双侧 t 尾概率（n 不大时用正态近似足够；补一个 t 分布的展开）
def t_sf_two_sided(t, df):
    # 正态近似 + 一阶 t 修正
    z = abs(t)
    p_norm = math.erfc(z / math.sqrt(2.0))
    # t 比正态重尾：(z^3+z)/(4df) 一阶修正
    p_t = p_norm + (z ** 3 + z) / (4.0 * df) * math.exp(-0.5 * z * z) / math.sqrt(2 * math.pi)
    return min(1.0, max(0.0, p_t))


def rankdata(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    ranks = [0.0] * len(xs)
    i = 0
    while i < len(xs):
        j = i
        while j + 1 < len(xs) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def kruskal_wallis(groups):
    """三档非参数比较；返回 H 统计量与上侧 chi^2(2) 近似 p。
    全样本同值（退化，例如坏哈希让所有种子产出同一场）时返回 H=0, p=0。"""
    allv = [x for g in groups for x in g]
    n = len(allv)
    if len(set(allv)) <= 1:
        return 0.0, 0.0
    ranks = rankdata(allv)
    # 结修正
    ties = {}
    for x in allv:
        ties[x] = ties.get(x, 0) + 1
    tie_corr = 1.0 - sum(t ** 3 - t for t in ties.values()) / (n ** 3 - n)
    if tie_corr <= 0.0:
        return 0.0, 0.0
    H = 0.0
    idx = 0
    for g in groups:
        ng = len(g)
        rsum = sum(ranks[idx:idx + ng])
        H += rsum ** 2 / ng
        idx += ng
    H = (12.0 / (n * (n + 1)) * H - 3 * (n + 1)) / tie_corr
    # chi^2(df=2) 上侧尾：p = exp(-H/2)
    p = math.exp(-0.5 * H)
    return H, p


def corr_len(ac):
    """1/e 相关长度（像素）。"""
    for d in range(1, len(ac)):
        if ac[d] < math.exp(-1):
            return d - 1 + (ac[d - 1] - math.exp(-1)) / (ac[d - 1] - ac[d])
    return float("nan")


def per_seed_features(rows):
    out = []
    for r in rows:
        acx = [r[f"acx{d}"] for d in range(LAGP)]
        acy = [r[f"acy{d}"] for d in range(LAGP)]
        acm = [0.5 * (acx[d] + acy[d]) for d in range(LAGP)]
        out.append({
            "tier": r["tier"],
            "seed": r["seed"],
            "k3": r["seed"] & 7,
            "mean": r["mean"],
            "var": r["var"],
            "hf": r["hf"],
            "range": r["max"] - r["min"],
            "crossings": r["crossings"],
            "runlen": r["runlen"],
            "iso": r["iso"],
            "clen_px": corr_len(acm),
            "clen_cells": corr_len(acm) / PPC,
            "acx": acx, "acy": acy, "acm": acm,
            "xcorr": r["xcorr_next"],
        })
    return out


def h1_test(feats, alpha=0.05):
    """H1: 档内相邻种子场互相关 ~0。"""
    res = {}
    m_tests = 3
    for t in TIERS:
        xs = [f["xcorr"] for f in feats if f["tier"] == t and not math.isnan(f["xcorr"])]
        if sd(xs) == 0.0:
            tt, p = float("inf"), 0.0
        else:
            tt = tstat_onesample(xs)
            p = t_sf_two_sided(tt, len(xs) - 1)
        res[t] = {"mean": mean(xs), "sd": sd(xs), "t": tt,
                  "p": p, "pass": p > alpha / m_tests, "n": len(xs)}
    return res


def h2_test(feats, metrics, alpha=0.05):
    """H2: 每个特征三档分布一致（Kruskal-Wallis）+ 效应量（极差/合并sd）。"""
    res = {}
    m_tests = len(metrics)
    for name in metrics:
        groups = [[f[name] for f in feats if f["tier"] == t] for t in TIERS]
        # 退化检测：档内零方差说明种子完全不起作用，直接判失败
        zero_var = any(sd(g) == 0.0 for g in groups)
        H, p = kruskal_wallis(groups)
        sds = [sd(g) for g in groups]
        if sum(sds) == 0.0:
            drift = float("inf")
        else:
            pooled_sd = math.sqrt(sum((len(g) - 1) * s_ ** 2
                                      for g, s_ in zip(groups, sds))
                                  / sum(len(g) - 1 for g in groups))
            drift = (max(mean(g) for g in groups) - min(mean(g) for g in groups)) / pooled_sd
        res[name] = {"H": H, "p": p, "drift": drift,
                     "means": {t: mean(g) for t, g in zip(TIERS, groups)},
                     "sds": {t: sd(g) for t, g in zip(TIERS, groups)},
                     "pass": (not zero_var) and p > alpha / m_tests and drift < 0.5}
    return res


def h3_test(feats, metrics, alpha=0.05):
    """H3: 档内按种子低 3 位分 8 组，尺度类指标组间无显著差异。
    三档分别做 Kruskal-Wallis；任一挂 => 种子位与特征尺度耦合（违反 C2）。"""
    res = {}
    for name in metrics:
        per_tier = {}
        any_fail = False
        for t in TIERS:
            groups = [[f[name] for f in feats if f["tier"] == t and f["k3"] == g]
                      for g in range(8)]
            if any(len(g) == 0 for g in groups):
                per_tier[t] = None
                continue
            H, p = kruskal_wallis(groups)
            ok = p > alpha / (len(metrics) * 3)
            per_tier[t] = {"H": H, "p": p, "pass": ok,
                           "group_means": [mean(g) for g in groups]}
            any_fail = any_fail or (not ok)
        res[name] = {"per_tier": per_tier, "pass": not any_fail}
    return res


CONFIGS = [
    ("good_vn.csv",       "good + 值噪声(合格基线)"),
    ("good_vn_fbm.csv",   "good + fBm 4 八度"),
    ("good_nn.csv",       "good + 最近邻(采样坏法, 混叠)"),
    ("good_sx.csv",       "good + 种子耦合频率(管线坏法, 违反 C2)"),
    ("bad_linear_vn.csv", "bad-linear + 值噪声(哈希坏法)"),
]

METRICS = ["var", "clen_cells", "iso", "runlen", "crossings", "hf"]
H3_METRICS = ["clen_cells", "hf"]


def main():
    outdir = sys.argv[1]
    report_path = sys.argv[2] if len(sys.argv) > 2 else None
    lines = []

    def emit(s=""):
        lines.append(s)
        print(s)

    all_res = {}
    for fname, label in CONFIGS:
        path = os.path.join(outdir, fname)
        feats = per_seed_features(load(path))
        emit("=" * 78)
        emit(f"配置: {label}  ({fname})")
        emit("=" * 78)

        emit("\n[三档特征统计] 每档 128 个种子, mean ± sd")
        hdr = f"{'特征':<12}" + "".join(f"{TIER_LABEL[t]:>22}" for t in TIERS)
        emit(hdr)
        for name in METRICS:
            cells = []
            for t in TIERS:
                xs = [f[name] for f in feats if f["tier"] == t]
                cells.append(f"{mean(xs):>10.5f} ± {sd(xs):.5f}")
            emit(f"{name:<12}" + "".join(f"{c:>22}" for c in cells))

        emit("\n[H1 跨种子独立性] 相邻种子场互相关应 ~0")
        h1 = h1_test(feats)
        for t in TIERS:
            r = h1[t]
            emit(f"  {t}: xcorr={r['mean']:+.5f} ± {r['sd']:.5f}  "
                 f"t={r['t']:+.2f} p={r['p']:.4f}  "
                 f"{'PASS' if r['pass'] else 'FAIL'}")

        emit("\n[H2 尺度独立性(跨档)] 三档特征分布 Kruskal-Wallis + 跨档漂移(单位:合并sd)")
        h2 = h2_test(feats, METRICS)
        for name, r in h2.items():
            emit(f"  {name:<12} H={r['H']:7.3f} p={r['p']:.4f}  "
                 f"drift={r['drift']:.3f}  {'PASS' if r['pass'] else 'FAIL'}")

        emit("\n[H3 尺度独立性(档内)] 按种子低3位分8组的 Kruskal-Wallis")
        h3 = h3_test(feats, H3_METRICS)
        for name, r in h3.items():
            for t in TIERS:
                pt = r["per_tier"][t]
                if pt is None:
                    continue
                emit(f"  {name:<12} {t}: H={pt['H']:8.3f} p={pt['p']:.4g}  "
                     f"{'PASS' if pt['pass'] else 'FAIL'}")

        # 场级"指纹"：平均自相关曲线的跨档最大差异（同一配置内）
        ac_means = {}
        for t in TIERS:
            sel = [f for f in feats if f["tier"] == t]
            ac_means[t] = [mean([f["acm"][d] for f in sel]) for d in range(LAGP)]
        maxdev = 0.0
        for d in range(1, LAGP):
            for a in TIERS:
                for b in TIERS:
                    maxdev = max(maxdev, abs(ac_means[a][d] - ac_means[b][d]))
        emit(f"\n  平均自相关曲线跨档最大偏差(滞后1..32px): {maxdev:.5f}")
        hf_mean = mean([f["hf"] for f in feats])
        all_res[fname] = {"h1": h1, "h2": h2, "h3": h3,
                          "ac_maxdev": maxdev, "hf": hf_mean}

    # 混叠阈值按八度数分组：fBm 合法地含有更高频能量，不能跟单八度基线比。
    # 同哈希、同八度、合格采样(vn)的 hf 才是参照系。
    base_files = {1: "good_vn.csv", 4: "good_vn_fbm.csv"}
    hf_thresh = {}
    for octv, bf in base_files.items():
        sel = per_seed_features(load(os.path.join(outdir, bf)))
        m, s = mean([f["hf"] for f in sel]), sd([f["hf"] for f in sel])
        hf_thresh[octv] = (m, s, m + 5.0 * s)

    emit("\n" + "#" * 78)
    emit("# 归因汇总")
    emit("#" * 78)
    emit("# 判据:")
    emit("#   哈希级 hash_test 先筛；场级 H1 挂 => 哈希引入跨格点/跨种子相关(哈希问题)；")
    emit("#   H1 过、H3 挂且尺度类指标随种子位变化 => 种子耦合进尺度(管线/哈希构造问题)；")
    emit("#   H1/H2/H3 全过、仅 hf 超同八度基线阈值 => 采样端混叠(采样问题)。")
    for octv, (m, s, th) in hf_thresh.items():
        emit(f"# hf 基线(oct={octv})={m:.5f} (sd={s:.5f}), 混叠阈值={th:.5f}")
    for fname, label in CONFIGS:
        r = all_res[fname]
        octv = 4 if "fbm" in fname else 1
        _, _, th = hf_thresh[octv]
        h1_pass = all(v["pass"] for v in r["h1"].values())
        h2_pass = all(v["pass"] for v in r["h2"].values())
        h3_pass = all(v["pass"] for v in r["h3"].values())
        hf_alias = r["hf"] > th
        overall = h1_pass and h2_pass and h3_pass and not hf_alias
        emit(f"  {label:<34} H1={'PASS' if h1_pass else 'FAIL':<4} "
             f"H2={'PASS' if h2_pass else 'FAIL':<4} "
             f"H3={'PASS' if h3_pass else 'FAIL':<4} "
             f"hf={r['hf']:.4f}{'*' if hf_alias else ' '} "
             f"[{'PASS' if overall else 'FAIL'}]")

    if report_path:
        with open(report_path, "w") as fh:
            fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
