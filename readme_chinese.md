# 短程有序（SRO）参数计算用合金筛选系统

## 概述

该计算框架用于对多组元合金成分进行系统化筛选，以计算并分析短程有序（SRO）参数、有效对相互作用（EPI）、失配体积以及归一化 SRO 参数。系统将分子动力学模拟（通过 LAMMPS）与统计热力学相结合，用于预测复杂合金在一系列温度下的有序化行为。

## 项目结构

```
Project/
├── main_screening.py                   # 主执行脚本（交互式 CLI）
├── calculators/                        # 核心计算模块
│   ├── __init__.py
│   ├── EPI_calculator.py               # 有效对相互作用（EPI）计算
│   ├── calculate_misfit_volume.py      # 失配体积计算（来自 database.json）
│   ├── calculate_sro_parameters.py     # Warren–Cowley SRO 参数计算
│   └── normalized_sro_parameters.py    # 归一化 SRO 参数 β_i/β_0
├── utils/                              # 辅助与工具模块
│   ├── __init__.py
│   ├── composition.py                  # 成分解析与命名
│   ├── create_calculation_structure.py # 计算目录结构创建 + CLI 辅助
│   ├── generate_calculation_files.py   # LAMMPS/求解脚本模板（按每个对/每个壳层生成）
│   ├── create_average_potential.py     # EAM 势函数加权平均工具
│   └── normalized_sro_screening.py     # 归一化 SRO 筛选 + 报告生成
└── data/
    ├── input/
    │   ├── database.json          # 元素性质数据库
    │   ├── potentials/
    │   │   ├── eam/              # EAM 势函数文件
    │   │   └── average_potentials/ # 生成的平均势函数
    │   └── *.xls                 # 成分输入文件
    └── output/
        └── results/
            └── {alloy_system}/   # 按合金体系组织结果
                └── {alloy_name}/ # 单个合金结果
                    ├── {alloy_name}.json              # 配置文件
                    ├── {alloy_name}_epi_results.txt   # EPI 结果
                    ├── {alloy_name}_misfit_volumes_results.txt # 失配体积
                    ├── {alloy_name}_sro_results.json  # SRO 参数
                    ├── {alloy_name}_fugacity_results.txt # 逸度值
                    ├── {alloy_name}_normalized_sro_parameters.txt # 归一化 SRO
                    ├── {pair_name}/                   # 成对目录
                    │   └── {nn}/                     # 近邻壳层
                    │       ├── calculate_{nn}_S-S_interaction_by_pressure_relaxation.py
                    │       ├── run_calculation_{nn}.sh
                    │       └── log*.lammps
                    └── all_alloys_normalized_sro_report.txt # 汇总报告
```


## 快速开始（最简）

1. 确保输入文件存在：
   - `data/input/database.json`
   - `data/input/potentials/eam/` 中的 EAM 势函数
2. 运行交互式 CLI：

```bash
python main_screening.py
```

## 方法学

### 1. 成分处理

系统接受合金成分的方式包括：
- **交互式输入**：手动输入元素含量
- **Excel 文件输入**：从包含成分数据的 `.xls` 文件批量处理

成分会被解析与归一化，生成：
- 合金名称（例如 `Fe57Ni19Cr24`）
- 组元列表与浓度
- 组元数与唯一元素对数量

### 2. 势函数生成

对于每个成分，系统：
1. **定位 EAM 势函数**：搜索与合金组元匹配的嵌入原子法（EAM）势函数文件
2. **生成平均势函数**：使用 `generate_averaged_potential` 生成按浓度加权的平均势函数
3. **保存平均势函数**：将平均势函数存到 `data/input/potentials/average_potentials/`

### 3. 计算结构搭建

系统创建分层目录结构：
```
{alloy_name}/
├── {pair_name}/          # 每个唯一元素对（例如 CrFe、CrNi、FeNi）
    └── {nn}/             # 每个近邻壳层（1st、2nd、3rd、4th）
        ├── calculate_{nn}_S-S_interaction_by_pressure_relaxation.py
        └── run_calculation_{nn}.sh
```

### 4. LAMMPS 计算

对于每个元素对与近邻壳层，系统：
1. **生成 Python 脚本**：生成用于溶质-溶质相互作用的 LAMMPS 计算脚本
2. **生成 shell 脚本**：生成用于执行所有对组合（11、12、22）的 bash 脚本
3. **执行计算**：通过 subprocess 调用运行 LAMMPS 模拟
4. **提取结果**：解析 log 文件提取相互作用能量


### 5. 有效对相互作用（EPI）

`EPI_calculator` 模块：
- 解析每个元素对/壳层组合的 LAMMPS log 文件
- 计算有效对相互作用能：`V_eff = E_12 - (E_11 + E_22)/2`
- 汇总所有元素对与壳层结果
- 保存到 `{alloy_name}_epi_results.txt`

### 6. 失配体积计算

`misfit_volume_calculator` 模块：
- 从 `database.json` 读取原子体积
- 计算平均原子体积：`V_avg = Σ(c_i × V_i)`
- 计算失配体积：`ΔV_i = V_i - V_avg`
- 保存到 `{alloy_name}_misfit_volumes_results.txt`

**备注**：使用按成分加权的平均体积与（类 Vegard 定律的）线性混合近似一致，这在溶质强化模型中很常见，
例如 Varvenne 等（2016），DOI：`10.1016/j.actamat.2016.07.040`。

### 7. SRO 参数计算

`calculate_sro_parameters` 模块：
- **逸度计算**：基于 EPI 计算随温度变化的逸度
- **SRO 参数计算**：为每个元素对与壳层计算 Warren-Cowley SRO 参数（α_nm）：
  - 二元体系：由逸度直接计算
  - 三元体系：通过矩阵求解耦合方程
  - 四元/五元体系：扩展矩阵形式
- **温度依赖**：在指定温度范围内计算 SRO 参数
- 保存到 `{alloy_name}_sro_results.json` 和 `{alloy_name}_fugacity_results.txt`

### 8. 归一化 SRO 参数

`normalized_sro_parameters` 模块：
- **Beta 计算**：结合失配体积计算归一化 SRO 参数（β_i）：
  ```
  β_i = Σ_n,m c_n × c_m × (ΔV_n - ΔV_m)² × α_nm^(i)
  ```
- **Beta_0 归一化**：计算参考值：
  ```
  β_0 = Σ_i c_i × (ΔV_i)²
  ```
- **相依赖平均**：
  - **FCC**：`β_1/β_0` 与 `(β_2 + 4β_3 + 2β_4)/β_0`
  - **BCC**：`(1.5β_1 + β_2)/β_0` 与 `(β_3 + 1.8β_4)/β_0`
- 保存到 `{alloy_name}_normalized_sro_parameters.txt`

### 9. 报告生成

系统生成综合报告：
- **单个合金报告**：每个成分的详细结果
- **汇总报告**：合金体系内所有合金的聚合结果（`all_alloys_normalized_sro_report.txt`）
- 报告包含归一化 SRO 参数、条件满足标记与温度依赖行为

## 使用方法

### 基本运行

```bash
python main_screening.py
```

### 交互式流程

1. **模式选择**：
   - `1`：计算新成分（完整计算流水线）
   - `2`：从已有文件加载（后处理/报告）

2. **合金体系**：输入合金体系名称（例如 `FeNiCr`）

3. **晶体相**：选择 `fcc` 或 `bcc`

4. **近邻壳层**：指定要计算的壳层数（1-4）. 不过，目前仅实现了第一最近邻（1st NN）的计算。

5. **晶胞数**：用于 LAMMPS 模拟的晶胞数量（默认：5）. 晶胞数量越多，计算成本越高。

6. **晶格常数**：（仅模式 1）为所有成分输入统一晶格常数

7. **成分输入**：（仅模式 1）提供包含成分的 Excel 文件路径

### 输出文件

#### 配置文件（`{alloy_name}.json`）
包含所有计算参数：
- 合金名称与成分
- 组元浓度
- 近邻壳层
- 晶格常数与晶体相
- 势函数类型

#### EPI 结果（`{alloy_name}_epi_results.txt`）
每个壳层与元素对组合的有效对相互作用能。

#### 失配体积（`{alloy_name}_misfit_volumes_results.txt`）
合金中各元素的失配体积。

#### SRO 参数（`{alloy_name}_sro_results.json`）
所有元素对与壳层的、随温度变化的 Warren-Cowley SRO 参数（α_nm）。

#### 逸度结果（`{alloy_name}_fugacity_results.txt`）
用于 SRO 计算的、随温度变化的逸度值。

#### 归一化 SRO 参数（`{alloy_name}_normalized_sro_parameters.txt`）
归一化 SRO 参数（β_i/β_0）与相依赖平均值。

#### 汇总报告（`all_alloys_normalized_sro_report.txt`）
体系内所有合金的聚合结果，并带有条件满足指示。

## 依赖

### Python 包
- `numpy`：数值计算
- `pandas`：Excel 文件读取
- `pathlib`：路径处理
- `json`：配置与数据存储
- `subprocess`：LAMMPS 执行
- `matplotlib`：可视化（可选）
- `scipy`：优化求解（用于 SRO 求解；也用于生成脚本）

### Python 版本
- 推荐：**Python 3.8+**（SciPy/matplotlib 兼容性）。

### 外部软件
- **LAMMPS**：分子动力学模拟器（计算必需）
- **mpi4py**：并行支持（可选，生成脚本中使用）

### 数据需求
- **EAM 势函数**：LAMMPS 格式的嵌入原子法势函数文件
- **数据库文件**：`database.json`，包含元素性质（原子体积、晶格常数等）

## 关键函数

### 主函数（`main_screening.py`）

- `process_composition()`：对单个成分编排完整计算流水线
- `main()`：用户交互与批处理协调

### 计算模块（`calculators/`）

- `calculate_epis()`：提取并计算有效对相互作用
- `misfit_volume_calculator()`：从原子体积计算失配体积
- `calculate_sro_parameters()`：计算 Warren-Cowley SRO 参数
- `calculate_normalized_sro_parameters()`：计算归一化 SRO 参数

### 工具模块（`utils/`）

- `get_alloy_name_and_concentrations()`：解析成分并生成命名
- `create_folder_structure()`：创建计算目录层级
- `find_potential()`：定位合适的 EAM 势函数
- `generate_averaged_potential()`：生成按浓度加权的平均势函数
- `create_shell_script()` / `create_python_script()`：生成 LAMMPS 计算脚本
- `save_normalized_sro_parameters()` / `load_normalized_sro_parameters()`：归一化 SRO 数据 I/O
- `generate_alloy_report()`：生成汇总报告

## 科学背景

### 短程有序（SRO）

SRO 参数（α_nm）量化偏离随机混合的程度：
- α_nm = 0：随机分布
- α_nm > 0：偏向团聚（clustering）
- α_nm < 0：偏向有序（ordering）

### 归一化 SRO 参数

归一化参数（β_i/β_0）同时包含化学与尺寸效应：
- **β_i**：壳层相关的归一化 SRO（包含失配体积）
- **β_0**：参考归一化因子
- 相依赖平均反映晶体结构效应

**文献依据**：
- Rao, Y. 和 Curtin, W. A.（2022），*Analytical models of short-range order in FCC and BCC alloys*，Acta Materialia 226, 117621。DOI：`10.1016/j.actamat.2022.117621`
- Nag, S. 和 Curtin, W. A.（2024），*Solute-strengthening in metal alloys with short-range order*，Acta Materialia 263, 119472。DOI：`10.1016/j.actamat.2023.119472`

### 有效对相互作用

EPI（V_nm）表示混合与分离对之间的能量差：
- 与有序/团聚倾向直接相关
- 通过逸度关系体现温度依赖
- 用于通过统计热力学预测 SRO 行为

## 备注

- 系统支持二元、三元、四元与五元合金
- 结果会缓存；若已有计算文件将被复用
- 系统会自动检查文件是否存在，以避免重复计算

## 作者与贡献（代码）

- **osintsevkirill**
  - 实现核心流水线模块，包括失配体积计算（`calculators/calculate_misfit_volume.py`）、归一化 SRO 参数（`calculators/normalized_sro_parameters.py`）以及 SRO 计算（`calculators/calculate_sro_parameters.py`）。
- **Xin Liu（EPFL）** — `xin.liu@epfl.ch`
  - 在 `utils/generate_calculation_files.py` 使用的生成运行脚本模板（shell 脚本头部注释）中被署名/致谢。

## 第三方代码与许可证

- `utils/create_average_potential.py` 包含改编自 `matscipy`（`average_atom.py`）的代码：
  - Lars Pastewka（U. Freiburg），2021
  - Wolfram G. Nöhring（U. Freiburg），2020
  - 项目：`matscipy` — `https://github.com/libAtoms/matscipy`
  - 许可证：GNU 通用公共许可证（GPL），v2 或更高版本（见 `http://www.gnu.org/licenses/`）

## 参考

该实现遵循以下领域的既有方法：
- Warren-Cowley SRO 参数计算
- EAM 势函数平均技术
- 合金有序化的统计热力学
- 多组元体系的失配体积计算

实现中引用的关键论文：
- Rao, Y. 和 Curtin, W. A.（2022），*Analytical models of short-range order in FCC and BCC alloys*，Acta Materialia 226, 117621。DOI：`10.1016/j.actamat.2022.117621`
- Nag, S. 和 Curtin, W. A.（2024），*Solute-strengthening in metal alloys with short-range order*，Acta Materialia 263, 119472。DOI：`10.1016/j.actamat.2023.119472`
- Varvenne, C., Luque, A., Curtin, W. A.（2016），*Theory of strengthening in fcc high entropy alloys*，Acta Materialia 118, 164–176。DOI：`10.1016/j.actamat.2016.07.040`

---

**版本**：1.0  
**最后更新**：2025 年 1 月

