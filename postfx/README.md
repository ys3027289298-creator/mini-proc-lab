# postfx

多后处理叠加的逐通道守恒校验。四个仿射 pass（gain / lift / overlay / tint）叠在 40×40
三通道图像上，对比两条管线：

- **守恒栈**：每个 pass 都写 `a*in + b`；
- **缺陷栈**：tint 的 G 通道被常量覆盖（吞通道），其余完全相同。

每条管线的输出都与两条独立解析路径交叉校验：

1. 逐像素仿射复合 $\hat X=A s+B$（推导见 [DERIVATION.md](DERIVATION.md)）；
2. 闭合式能量簿记 $E_2(aX+b)=a^2E_2(X)+2a\langle X,b\rangle+E_2(b)$，只走矩递推、不碰像素。

`make test` 编译运行 C 自检并调用 `tools/render.py`（`python3 -I` 隔离模式）渲染，
产物输出到 `../artifacts/postfx/`：

| 文件 | 内容 |
|---|---|
| `report.txt` | 分阶段三通道实测/解析/簿记能量全表与判定 |
| `energy_table.md` | Markdown 三通道能量变化表 |
| `diffmap_conserving.png` | 守恒栈逐像素 \|实测−解析\| 差异图（三通道） |
| `diffmap_buggy.png` | 缺陷栈有符号逐像素差异图（三通道，独立色标） |
| `stats.json` | 逐通道最大绝对/相对差异与 PASS/DETECTED 判定 |
| `postfx.bin` | PFX1 二进制 bundle（输入、两套实测、两套解析、阶段能量） |

判据：任一通道相对差超过 $10^{-3}$ 即未完成。当前守恒栈最差 $8.6\times10^{-8}$；
缺陷栈 G 通道实测 $9.97\times10^{-1}$，吞通道被定位，R/B 保持干净。
