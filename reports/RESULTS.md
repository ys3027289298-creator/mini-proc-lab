# 验证结果

所有数字由 `scripts/run_all.sh` 生成，可随时重跑覆盖。主种子
`20261124`，世界 256×256。

## 一、每层同种子跑两次

`reports/verify/repeat.csv`，逐层详情见 `reports/verify/layer-<name>.md`。

| 层 | 输出 | 两次运行 |
|---|---|---|
| noise | `noise.f32` (262144 B) | 逐位相同 |
| erosion | `erosion.f32` (262144 B) | 逐位相同 |
| biome | `biome.u8` (65536 B) | 逐位相同 |
| hydro | `hydro.f32` (262144 B) | 逐位相同 |
| scatter | `scatter.f32` (40652 B) | 逐位相同 |

每层的 sha256 见对应 `layer-<name>.md`。

## 二、线程数 1 / 4 / 8

`reports/verify/threads.csv`。

| 层 | 1 | 4 | 8 |
|---|---|---|---|
| noise | 基准 | 相同 | 相同 |
| erosion | 基准 | 相同 | 相同 |
| biome | 基准 | 相同 | 相同 |
| hydro | 基准 | 相同 | 相同 |
| scatter | 基准 | 相同 | 相同 |

另外用线程数 2、3、5、7、16、32 各跑一遍，五层仍只有一个哈希值。合并模式
`--layer all` 在 1/4/8 下五个文件也全部一致。

### 非确定性归约的定位验证

`scripts/nondeterminism_demo.sh` 用 `worldgen-nondet` 构建（唯一区别是 hydro
里接入了一个“每线程一个部分和、按线程完成顺序相加”的归约）：

| 层 | 1 线程 | 4 线程 | 8 线程 |
|---|---|---|---|
| noise | 基准 | 相同 | 相同 |
| erosion | 基准 | 相同 | 相同 |
| biome | 基准 | 相同 | 相同 |
| **hydro** | 基准 | **不同** | **不同** |
| scatter | 基准 | 相同 | 相同 |

结论：比对流程能锁定到引入非确定性归约的那一层（hydro），而不是只报“整体
不一致”。同一台机器上用正常构建跑同样输入，五层在 1/4/8 下完全一致。

## 三、接缝

`reports/verify/seams.csv`：在 512×512 世界上切成 8×8 个 64×64 的窗口分别
生成，再拼回整幅世界逐字节比较（用 512 世界的原因见第六节）。

| 层 | 比较字节 | 不一致字节 |
|---|---:|---:|
| noise | 1048576 | 0 |
| erosion | 1048576 | 0 |
| biome | 262144 | 0 |
| hydro | 1048576 | 0 |
| scatter | 18861 个点 | 0 |

在 256 世界上另对 4 个种子 × 10 个窗口（含贴边、细长、跨角）做过同类比较，
40 组全部通过。

## 四、换主种子

主种子改 1 位后，五层输出全部改变，说明每层的子种子都真正接入了生成过程，
没有哪层仍在用常量。

## 五、子种子派生

`worldgen --dump-seeds` 与 `tools/seed_reference.py`（独立 Python 实现）对
`20261124`、`20261123`、`0`、`1`、`123456789`、`0xDEADBEEFCAFEF00D` 的输出
逐字符相同；`123456789` 的样例值写在 `docs/seed-derivation.md` 表格中，由
`scripts/verify_seeds.sh` 每次回归时核对。

## 六、接缝检查不是空转

256 世界下 hydro 的 halo 是 116 格，任何 64×64 窗口都会被撑到整幅世界，接缝
比较就退化成“整幅重算”，说明不了 halo 够不够。因此接缝测试改用
`worldgen-512`（512×512 世界），窗口加 halo 后仍留有边界余量。

再做一个负向对照 `worldgen-shorthalo-512`：同样的 512 世界，但把 halo 故意
改成 1～5 格。用 `scripts/halo_demo.sh` 比较：

| 层 | 正常 halo | 短 halo |
|---|---:|---:|
| noise | 0 | 0（该层无 halo，符合预期） |
| erosion | 0 | 37120 字节不一致 |
| biome | 0 | 320 字节不一致 |
| hydro | 0 | 41984 字节不一致 |
| scatter | 0 | 1 个窗口点集不一致 |

正常构建零差异、短 halo 全部出现断裂，说明接缝检查确实能抓到问题。

## 七、工具链检查

- 单元测试 36 项全部通过（`make -C worldgen test`）。
- ThreadSanitizer：8 线程跑 `--layer all`，无数据竞争报告。
- AddressSanitizer + UBSan：整幅世界与窗口模式各跑一次，无越界、无未定义行为。

## 产物

- `reports/verify/` — CSV 与 Markdown 报告
- `reports/nondeterminism-demo/` — 负向对照结果
- `reports/render/` — 四层的 PNG 预览（`tools/render.py` 生成）
- `reports/halo-demo/` — 接缝负向对照结果

预览图仅用于目测地貌是否连续；判断一致性以逐字节比较为准。
