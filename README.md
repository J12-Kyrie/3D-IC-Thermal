# 3D-IC Analytical Thermal Simulation Framework (3D-IC热分析解析框架)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![Paper](https://img.shields.io/badge/paper-TVLSI'25-b31b1b.svg)](https://ieeexplore.ieee.org/xpl/RecentIssue.jsp?punumber=9444)
[![Stars](https://img.shields.io/github/stars/J12-Kyrie/3D-IC-Thermal?style=social)](https://github.com/J12-Kyrie/3D-IC-Thermal/stargazers)

一种针对三维集成电路（3D-IC）的高效、高精度稳态热分析解析框架，能够精确处理任意形状的功率密度分布。

该仓库是论文 **"An Analytical 3D-IC Thermal Simulation Framework Using Adaptive Rectangular Approximation and Conformal Mesh Method"** 的官方代码实现。

---

## ✨ 核心特性

**处理任意功率密度图**: 创新性地提出 **自适应矩形近似算法**，能够将任意曲线形状的功率模块高效、精确地离散化为矩形热源集合。
**高效解析求解**: 基于 **域分解法 (Domain Decomposition Method)** 和 **傅里叶级数展开** ，对多层3D-IC结构进行全解析求解，避免了传统数值方法的离散化误差。
**共形网格与步进积分**: 结合 **共形网格生成算法 (Conformal Mesh Generation)** 与解析步进积分，极大地提高了傅里叶系数的计算效率和精度。
**卓越的性能**: 与商业有限元软件（如COMSOL）相比，本框架在保持高精度（最大绝对误差低于0.5 K）的同时，实现了 **高达60倍** 的计算速度提升。
**支持各向异性材料**: 模型在求解过程中考虑了材料热导率的各向异性，更贴近真实物理情况。

## 🚀 安装指南

本项目使用Python实现，并依赖于NumPy等科学计算库。

1.  **克隆仓库**
    ```bash
    git clone [https://github.com/J12-Kyrie/3D-IC-Thermal.git](https://github.com/J12-Kyrie/3D-IC-Thermal.git)
    cd 3D-IC-Thermal
    ```

2.  **创建虚拟环境 (推荐)**
    ```bash
    python -m venv venv
    source venv/bin/activate  # On Windows, use `venv\Scripts\activate`
    ```

3.  **安装依赖**
    ```bash
    pip install -r requirements.txt
    ```

## 💡 快速开始

以下是一个如何使用本框架进行3D-IC热仿真的基本示例。

```python
import numpy as np
from thermal_solver import AnalyticalThermalSolver, PowerMap, ChipStructure

# 1. 定义芯片结构
# 定义各层材料、厚度和热导率
chip_layers = [
    {'name': 'HeatSink', 'thickness': 5e-3, 'conductivity': [400, 400, 400]},
    {'name': 'TIM1', 'thickness': 0.5e-3, 'conductivity': [2, 2, 2]},
    {'name': 'Chip', 'thickness': 0.5e-3, 'conductivity': [130, 130, 130]},
    {'name': 'TIM2', 'thickness': 0.5e-3, 'conductivity': [2, 2, 2]},
]
structure = ChipStructure(layers=chip_layers, width=0.024, length=0.024)

# 2. 加载并近似功率密度图
# power_map_data 可以是从EDA工具导出的不规则形状功率数据
power_map = PowerMap.from_arbitrary_shape(power_map_data)
approximated_rects = power_map.adaptive_rectangular_approximation(error_threshold=6e-10)

# 3. 初始化并运行求解器
solver = AnalyticalThermalSolver(
    structure=structure,
    power_rects=approximated_rects,
    eigenvalues=(30, 30) # 设置傅里叶级数截断项数
)
temperature_distribution = solver.solve()

# 4. 可视化或分析结果
solver.plot_temperature_map(layer_name='Chip')
max_temp = temperature_distribution.max()
print(f"芯片最高温度为: {max_temp:.2f} K")
