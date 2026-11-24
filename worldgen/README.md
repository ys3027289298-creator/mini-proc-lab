# 世界生成层（worldgen）

五层共用一份主种子，各自持有派生出的子种子。窗口、线程数、重复运行都不改变
输出。详细规则见 `../docs/seed-derivation.md` 和 `../docs/output-format.md`。

```sh
make            # 生成 ./worldgen
make test       # 单元测试 + 文档种子样例
```

```sh
./worldgen --seed 20261124 --layer all --threads 4 --out out/world
./worldgen --dump-seeds 20261124
```

## 每层的范围

| 层 | 依赖 | 依赖的子种子 | 说明 |
|---|---|---|---|
| noise | — | noise | 6+3 个八度的 fBm，基础高度 |
| erosion | noise | erosion | 64 轮热侵蚀，单目的地搬运 |
| biome | noise, erosion | biome | 高度/坡度 + 温度/湿度两个 fBm |
| hydro | noise, erosion | hydro | D8 流向 + 48 轮汇流 Jacobi |
| scatter | noise, erosion | scatter | 2.0 格抖动网格，最小间距 1.1 |

`--layer <name>` 只写该层文件；`--layer all` 依次写全部五个。每层都从主种子
重算所需的上游层，不读取其他层的输出，因而单层可比对、可并行。

## 为什么接缝不会断

每个层在请求窗口外多算一圈 halo，大小覆盖该层信息传播的最远距离
（`wg_halo_for_layer`）。halo 之外的一切都按世界坐标求值，窗口内部结果与整幅
世界在对应区域逐位相同。`scripts/verify.sh` 用 4×4 个 64×64 的独立窗口拼接
后与整幅世界逐字节比较。

## 负数对照

`make worldgen-nondet` 会额外定义 `SP_SELFTEST_NONDET`，让 hydro 用一个
“按线程完成顺序合并部分和”的均值去缩放输出。它**只**用于
`scripts/nondeterminism_demo.sh`，证明比对流程能定位到具体层；正常构建不包含
这段代码。
