#!/usr/bin/env python3
"""展开图结构校验：XML 合法性、元素边界、pattern 引用、超限面标注。
用法：python3 check_svg.py unfolding.svg
"""
import sys
import xml.dom.minidom as minidom


def attr(e, k):
    v = e.getAttribute(k)
    return v if v else ""


def main(path):
    d = minidom.parse(path)
    root = d.getElementsByTagName("svg")[0]
    cw, ch = (float(v) for v in root.getAttribute("viewBox").split()[2:])
    rects = d.getElementsByTagName("rect")
    texts = d.getElementsByTagName("text")
    patterns = d.getElementsByTagName("pattern")

    def num(r, k):
        v = r.getAttribute(k)
        return float(v) if v else 0.0

    # 1. 所有几何元素在 viewBox 内
    out = []
    for r in rects:
        if not r.getAttribute("width"):
            continue
        x, y, w, h = (num(r, k) for k in ("x", "y", "width", "height"))
        if x < -0.5 or y < -0.5 or x + w > cw + 0.5 or y + h > ch + 0.5:
            out.append((round(x, 1), round(y, 1), round(w, 1), round(h, 1)))
    assert not out, f"rects out of viewBox: {out}"

    # 2. pattern 引用必须有定义且无重复 id
    ids = {attr(e, "id") for e in d.getElementsByTagName("*") if attr(e, "id")}
    refs = {attr(r, "fill") for r in rects if attr(r, "fill").startswith("url(#")}
    missing = [u for u in refs if u[5:-1] not in ids]
    assert not missing, f"missing pattern defs: {missing}"
    pids = [attr(p, "id") for p in patterns]
    assert len(pids) == len(set(pids)), f"duplicate pattern ids: {pids}"

    # 3. 三档与超限面 E 的文字标注必须存在
    alltext = " ".join(t.firstChild.data for t in texts if t.firstChild)
    for needle in ("OVER", "face E", "东墙 E", "H/L=3.00", "北墙 N", "西墙 W"):
        assert needle in alltext, f"missing label: {needle!r}"

    # 4. 文字不越界
    for t in texts:
        x, y = float(t.getAttribute("x")), float(t.getAttribute("y"))
        assert 0 <= x <= cw and 0 <= y <= ch, f"text out of bounds: {x},{y}"

    print("svg_check=1 rects=%d texts=%d patterns=%d" % (len(rects), len(texts), len(patterns)))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "unfolding.svg")
