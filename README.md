# mini-proc-lab — 五层世界生成

同一主种子下，噪声、散布、侵蚀、群系、水系五层从同一个种子派生各自的子种子，
所有随机量都按**全局世界坐标**求值，因此：

- 同一个种子跑两次，五层输出逐位相同；
- 线程数 1 / 4 / 8，五层输出逐位相同；
- 任意窗口切块生成再拼回去，与整幅世界逐位相同（接缝不断）。

## 目录

| 路径 | 内容 |
|---|---|
| `common/seed.{c,h}` | 五层子种子派生（规则见 `docs/seed-derivation.md`） |
| `common/rng.{c,h}` | 坐标哈希、值噪声、fBm |
| `common/grid.{c,h}` | 带世界坐标的浮点栅格、双线性采样 |
| `common/parallel.{c,h}` | 按行静态切分的线程池；无跨线程归约 |
| `worldgen/src/noise.c` | 基础高度场 |
| `worldgen/src/erosion.c` | 热侵蚀（单目的地、双缓冲） |
| `worldgen/src/biome.c` | 群系分类（高度 + 坡度 + 两个气候场） |
| `worldgen/src/hydro.c` | D8 流向 + 汇流累加（Jacobi 迭代） |
| `worldgen/src/scatter.c` | 抖动网格 + 最小间距（邻域最优者胜） |
| `noise/ terrain/ scatter/` | 早期的独立原型，自带自测，保留不动 |
| `docs/` | 种子派生规则与输出格式 |
| `scripts/` | 复现/比对脚本 |
| `reports/` | 生成的结果（可删除重建） |

## 构建与测试

```sh
make -C worldgen              # 生成 worldgen/worldgen
make -C worldgen test         # 各层单元测试 + 文档种子样例核对
scripts/run_all.sh            # 全部比对，产出 reports/
```

`scripts/run_all.sh` 等价于：

```sh
scripts/verify.sh                 # 复现 / 线程数 / 接缝 / 换种子 / 种子派生
scripts/nondeterminism_demo.sh    # 负向对照：故意放一个非确定性归约
scripts/halo_demo.sh              # 负向对照：故意缩短 halo
scripts/verify_seeds.sh           # 文档样例 vs C vs Python
python3 -I tools/render.py reports/render   # 四层 PNG 预览
```

结果汇总见 `reports/RESULTS.md`。

## 用法

```sh
worldgen/worldgen --seed 20261124 --layer all --threads 4 --out out/world
worldgen/worldgen --seed 20261124 --layer hydro --threads 1 --out out/h1
worldgen/worldgen --seed 20261124 --layer all --threads 8 \
                  --window 96,32,64,64 --out out/tile
worldgen/worldgen --dump-seeds 20261124
```

`--layer` 取 `noise|erosion|biome|hydro|scatter|all`。输出格式见
`docs/output-format.md`。

## 确定性约定

1. 五层的随机量只来自各自的子种子，派生规则是一段可以照着复算的文本。
2. 噪声与阈值扰动按世界坐标求值，绝不使用“块内计数器”。
3. 侵蚀与汇流都是双缓冲 Jacobi：本轮只读上一轮、只写本轮。
4. `sp_parallel_rows` 按 `row % threads` 静态分配，回调只写自己那一行。
5. 散布用“邻域内排名最优先者胜”，结果与遍历顺序无关。
6. 编译带 `-ffp-contract=off -fno-fast-math`，不让编译器引入 FMA 造成末位差异。

`scripts/nondeterminism_demo.sh` 会把一个“按线程完成顺序合并部分和”的
归约接进 hydro，验证比对流程确实能把它抓出来并定位到该层。
`scripts/halo_demo.sh` 则把 halo 故意缩短，验证接缝检查不是空转。
