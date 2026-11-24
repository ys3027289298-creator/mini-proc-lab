# 子种子派生规则（v1）

一份世界只有一个 **主种子**（64 位无符号整数）。五个层各自使用一个从主种子
派生的 **子种子**，层与层之间不共享、不叠加。任何按本文实现的程序，对同一个
主种子都会算出完全相同的五个子种子。

## 规则

记主种子为 `m`（uint64），层名为 `L`，取值依次为：

```
noise  scatter  erosion  biome  hydro
```

对每个层分别执行：

1. 拼出字节串 `S = "mini-proc-lab/v1/" || L`（UTF-8/ASCII，**不含**结尾的
   NUL；共用前缀共 17 字节）。
2. 对 `S` 计算 64 位 FNV-1a 哈希，得到 `h`：
   - 初值 `h = 0xCBF29CE484222325`
   - 对每个字节 `b`：`h = ((h XOR b) * 0x100000001B3) mod 2^64`
3. 计算 `x = (m XOR h) mod 2^64`。
4. 对 `x` 做 SplitMix64 终结混合，结果即为子种子：
   - `x = (x + 0x9E3779B97F4A7C15) mod 2^64`
   - `x = ((x XOR (x >> 30)) * 0xBF58476D1CE4E5B9) mod 2^64`
   - `x = ((x XOR (x >> 27)) * 0x94D049BB133111EB) mod 2^64`
   - `子种子 = x XOR (x >> 31)`
5. 该层后续所有噪声、抖动、阈值扰动都只用这个子种子作为根；不得再混入主种子
   或其他层的子种子。

所有运算均为 64 位无符号模运算（`>>` 为逻辑右移），字节序只在拼接字节串时
按书写顺序处理，与机器无关。

## 可复算样例（主种子 m = 123456789 = 0x00000000075BCD15）

| 层 | `S` 的 FNV-1a 值 `h` | `x = m XOR h` | 子种子 `splitmix64(x)` |
|---|---:|---:|---:|
| noise | `0x019272F6A9B88BA0` | `0x019272F6AE6346B5` | `0x84E8FE37627F5BBA` |
| scatter | `0x44A3EBDA695D8C94` | `0x44A3EBDA6E064181` | `0xBD82B623209A525B` |
| erosion | `0xCAEF57DF0B1CA25F` | `0xCAEF57DF0C476F4A` | `0x6B80027587775C9D` |
| biome | `0xCB7DAC490AFEEDA6` | `0xCB7DAC490DA520B3` | `0x01978BED0659E5B4` |
| hydro | `0x70BF097403C63350` | `0x70BF0974049DFE45` | `0xB69653DD63920F89` |

十进制形式（顺序同上）：

```
noise    9577184121432529850
scatter 13655677282127860315
erosion  7746194062884428957
biome      114714165651301812
hydro   13156795571732877193
```

## 自行核对

```sh
make -C worldgen worldgen
./worldgen/worldgen --dump-seeds 123456789      # C 实现
python3 -I tools/seed_reference.py 123456789    # 独立 Python 实现，应逐字符相同
```

`scripts/verify_seeds.sh` 会比较两份实现和文档里的样例值。

## 兼容性

派生规则一旦发布就不能改动。需要改变规则时，必须同时修改
`SP_SEED_PREFIX`（例如改成 `mini-proc-lab/v2/`），这样旧世界不会在无声无息
中变成另一个样子。
