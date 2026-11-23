# terrain

`make test` 跑自检。

## 模块

- `field_sample`：双线性高度场采样（越界钳位到边界值）。
- `field_grad` / `field_slope`：边界感知梯度与坡度。内部节点二阶中心差分
  （截断误差 ≤ M₃h²/6），边界节点三点二阶单侧差分（≤ M₃h²/3），
  全网格保持 O(h²) 同阶。
- 法线：`n = (−fx, −fy, 1)/√(1+fx²+fy²)`，夹角上界
  `α = asin(ε/(1−sε))`。

边界坡度/法线的偏差上界推导、适用网格尺度与四边界验证见
[BOUNDARY_NORMALS.md](BOUNDARY_NORMALS.md)。

## 独立复核

```
python3 -I verify/verify_boundary_normals.py
```

纯 Python（float64）独立重推模板常数与四边界夹角，不与 C 实现共享代码；
输出存档于 `verify/results.txt`。
