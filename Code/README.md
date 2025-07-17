# 3D-IC Analytical Thermal Simulation Framework

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![Paper](https://img.shields.io/badge/paper-TVLSI'25-b31b1b.svg)](https://ieeexplore.ieee.org/xpl/RecentIssue.jsp?punumber=9444)
[![Stars](https://img.shields.io/github/stars/leonchenucr/AnalyticalThermal3DIC?style=social)](https://github.com/leonchenucr/AnalyticalThermal3DIC/stargazers) 

An efficient and high-accuracy analytical framework for steady-state thermal analysis of three-dimensional integrated circuits (3D-ICs), capable of precisely handling power density distributions of arbitrary shapes.

This repository is the official code implementation of the paper "An Analytical 3D-IC Thermal Simulation Framework Using Adaptive Rectangular Approximation and Conformal Mesh Method".

---

## ✨ Key Features  

**Handles Arbitrary Power Maps**: Proposes an innovative **Adaptive Rectangular Approximation Algorithm** to efficiently and accurately discretize arbitrarily shaped power modules into rectangular heat source sets.  

**Efficient Analytical Solution**: Employs **Domain Decomposition Method (DDM)** and **Fourier Series Expansion** for fully analytical thermal modeling of multi-tier 3D-ICs, eliminating discretization errors inherent in traditional numerical methods.  

**Conformal Mesh & Stepwise Integration**: Combines **Conformal Mesh Generation** with analytical stepwise integration to significantly enhance the computational efficiency and accuracy of Fourier coefficient calculations.  

**Outstanding Performance**: Compared to commercial FEM tools (e.g., COMSOL), this framework achieves **up to 60× speedup** while maintaining high accuracy (maximum absolute error < 0.5 K).  

**Supports Anisotropic Materials**: The model accounts for anisotropic thermal conductivity, better reflecting real-world physical scenarios.


## 🚀 Installation Guide  

This project is implemented in Python and relies on scientific computing libraries such as NumPy.  

1. **Clone the Repository**  
    ```bash
    git clone https://github.com/leonchenucr/AnalyticalThermal3DIC.git  
    cd AnalyticalThermal3DIC  
    ```  

2. **Create a Virtual Environment (Recommended)**  
    ```bash
    python -m venv venv  
    source venv/bin/activate  # On Windows, use `venv\Scripts\activate`  
    ```  

3. **Install Dependencies**  
    ```bash
    pip install -r requirements.txt  
    ```  

4. **Verify Installation**  
    Run the demo script to ensure everything works correctly:  
    ```bash
    python demo.py  
    ```  

## 💡 Quick Start

Here's a basic example of how to use this framework for 3D-IC thermal simulation:

### Basic Workflow
1. **Curve Recognition** (Process irregular power maps)
   ```bash
   python Curve_recognition.py
   ```

2. **Rectangle Filling** (Discretize power modules)
   ```bash
   python Rectangle_filling.py
   ```

3. **Power Mapping** (Assign power values)
   ```bash
   python Rectangular_corresponding_power.py
   ```

4. **Mesh Generation** (Create conformal mesh)
   ```bash
   python mesh_generation.py
   ```

5. **Thermal Simulation** (Run 11-layer 3D-IC analysis)
   ```bash
   python 3d_11layer.py
   ```
