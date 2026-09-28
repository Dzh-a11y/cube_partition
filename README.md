# 49 个小正方体：最大化尺寸种数

本模型把大正方体固定为 `[0,N]^3`，只考虑与坐标轴平行、位置和边长均为整数的小正方体。归一化到单位大正方体后，坐标和边长都是 `1/N` 的整数倍。

## 数学模型

候选集合 `C` 中的元素为 `c=(i,j,k,s)`：`i,j,k>=0`，`s>=1`，全部为整数，且 `i+s,j+s,k+s<=N`。变量 `x_c` 表示选不选候选，变量 `y_s` 表示用不用边长 `s`；全部是二元变量。

目标：`maximize sum_s y_s`。

约束：

1. 对每个单位格 `v`，`sum_{c: v in c} x_c = 1`，保证无重叠、无空隙并恰好填满大正方体。
2. `sum_c x_c = K`，默认 `K=49`。
3. 对每种尺寸 `s`，`y_s <= sum_{c: size(c)=s} x_c <= M_s*y_s`，其中 `M_s=min(K,floor(N^3/s^3))`。

候选中排除 `s^3+(K-1)>N^3` 的尺寸：其余每个小正方体至少占一个单位格，因此这些候选不可能被使用。无需另加总体积约束。

SciPy 的 `milp` 调用 HiGHS 求解器。由于 API 做最小化，代码传入的是 `-sum_s y_s`。设定相对 gap 为零；达到时间上限时会保留并验证可行解，但不会标成已证明最优。

## 运行

本项目现有的根目录 `.venv` 已安装 SciPy。从项目根目录可直接运行：

```sh
.venv/bin/python exercises/cube_partition_mip/solve.py --n 6 --count 49 --time-limit 60 --output outputs/cube_partition_mip/result_n6.json
.venv/bin/python -m unittest discover -s exercises/cube_partition_mip -v
```

已验证的 N=6、K=49 实例最优值为 **3 种尺寸**：边长 1 有 36 个，边长 2 有 9 个，边长 3 有 4 个。总数 `36+9+4=49`，总体积 `36+9*8+4*27=216=6^3`；程序另外逐格验证几何覆盖。完整坐标保存在项目根目录下的 `outputs/cube_partition_mip/result_n6.json`。

项目路径包含冒号 `:`，`virtualenv` 会拒绝在该路径下新建环境。直接使用上面的项目根目录解释器即可。IDE 中选择已有解释器 `/Users/douzihao/I2DL2026:2027/.venv/bin/python`，无需在本脚本目录下创建新的 `.venv`。

## 枚举全部可行分割（N=6，总数固定 49）

`enumerate_solutions.py` 保留逐格覆盖与总数约束，把目标设为常数 0。每找到一组 S，就加入 `sum_{c in S} x_c <= 48` 再求解。总数始终等于 49；求解器可重新选择全部候选。原 `solve.py` 的最大化尺寸种数功能保持独立。

从项目根目录运行一轮，时间单位为秒：

```sh
.venv/bin/python exercises/cube_partition_mip/enumerate_solutions.py --n 6 --count 49 --time-limit 60 --output outputs/cube_partition_mip/feasible_n6.jsonl
```

接着用同一个文件续跑（可反复运行）：

```sh
.venv/bin/python exercises/cube_partition_mip/enumerate_solutions.py --n 6 --count 49 --time-limit 60 --output outputs/cube_partition_mip/feasible_n6.jsonl --resume
```

每个解立即进行整数几何验证并追加保存，文件第一行为参数，其余每行是一个解。续跑会校验旧解并重建全部排除约束，不会重复保存。已有文件必须显式指定 `--resume`，否则拒绝覆盖。可选 `--max-solutions 100` 限制累计解数（包含之前保存的解）；默认不限制解数。

同目录的 `feasible_n6.summary.json` 记录状态。只有剩余模型被证明不可行时，才会出现 `status=exhausted, complete=true`；时限、解数上限或 Ctrl+C 都表示未确认找全。时间预算按每次启动独立计算，不包含加载并验证旧文件的时间。每次启动均会保存已获得的解；硬终止时摘要可能滞后，以 JSONL 中完整写入的记录为准。

枚举对象是不同的空间摆放，旋转和镜像分别计数；同一组候选的排列顺序不重复计数。尺寸数量组合单独汇总在摘要的 `size_histograms` 中。SciPy 每轮重新求解，累计排除约束较多时会变慢。

`--n` 可改为更细网格，`--count` 可改为其他个数。候选数随 N 很快增长，时间上限只约束求解阶段，不包含候选生成和矩阵构建。默认 N=6 用于提供一个小规模可运行实例，并不是对原题的额外限制。

JSON 中 `cubes` 给出全部选中候选的位置及边长，`size_counts` 给出各尺寸数量。`verified=true` 表示已用整数逐格计数独立检查数量、边界、覆盖和尺寸指示变量。`optimal_on_grid=true` 仅表示当前 N 下已证明最优。`size_count_upper_bound` 是求解器给出的尺寸种数上界；达到时间上限且没有可行解时，不会输出虚构的分割方案。

固定网格的最优值不能直接作为允许任意实数尺寸和位置时的全局最优值。不同 N 的网格未必相互包含；N 变为原来的整数倍时，原有分割才一定仍可表达。
