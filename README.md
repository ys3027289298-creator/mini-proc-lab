# mini-proc-lab

程序化地形与散布基础模块，C 实现，每目录独立构建可自测。

| 目录 | 内容 |
|---|---|
| `noise/` | 值噪声、fBm 叠加、确定性哈希 |
| `terrain/` | 高度场双线性采样、坡度 |
| `scatter/` | 抖动网格散布、最小间距拒绝、坡度密度衰减 |
| `slope/` | 坡度–密度衰减曲线拟合、置信区间、分层统计、上线校验 |

## 坡度与散布密度

`slope/tools/slope_fit_report` 把散布实测数据拟合成衰减曲线，给出参数置信区间与按坡度分层的统计，并对一组部署参数做上线判定。详见 `slope/README.md`。

```sh
make -C noise && make -C terrain && make -C scatter && make -C slope
make -C slope report
```

构建产物（`*.o`、`*.a`、测试与工具可执行文件）未纳入版本管理。
