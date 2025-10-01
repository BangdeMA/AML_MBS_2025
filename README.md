**Project Overview**

This project integrates clinical data and computational models to analyze AML-related dynamics.
The workflow includes:

* Parameter optimization using Julia (SDE_Optimization.jl);

* Numerical simulation using C++ (main.cpp);

* Visualization using MATLAB (plotFigX.m);

**File Description**

Clinical Data (TCGA_OS.csv,TCGA_clinical_OS.csv): Overall survival (OS) data, obtained from The Cancer Genome Atlas (TCGA).

**Julia Script**

* SDE_Optimization.jl

* Performs parameter optimization based on a stochastic differential equations (SDEs) model.

* Input: TCGA data / initial parameter guesses.

* Output: optimized parameter set.

**C++ Program**

* main.cpp with optimized parameters (from Julia).

* Main numerical simulation code.

* Output: result files in .txt format.

**MATLAB Scripts**

* plotFig2.m, plotFig5.m, plotFig8.m, plotFig9.m, plotFig10.m, plotFig34.m, plotFig67.m

* Plotting scripts corresponding to Figures 2, 5, 8, 9, 10, 3&4, and 6&7 in the paper.

* They read the .txt output from main.cpp and generate figures.

**Workflow**
flowchart TD

    A[Clinical Data\n(TCGA_OS.csv)] --> B[Julia\nSDE_Optimization.jl\n(Parameter Optimization)]
    B --> C[C++\nmain.cpp\n(Numerical Simulation)]
    C --> D[Text Results\n(result_*.txt)]
    D --> E[MATLAB\nplotFigX.m\n(Figure Generation)]
    E --> F[Final Figures]

**How to Run**

**Prepare Data**
* Place TCGA_OS.csv in the project root folder.

**Run Julia Optimization**
* julia SDE_Optimization.jl
→ Generates optimized parameters.

**Run C++ Simulation**
→ Outputs simulation results (result_*.txt).

**Run MATLAB Plotting**
* Open MATLAB and execute, for example:plotFig2
→ Generates corresponding figures.

**Requirements**

* Julia ≥ 1.9 (with packages: DifferentialEquations.jl, Optim.jl, et al.)

* C++11 or later compiler

* MATLAB ≥ R2024b

**Acknowledgment**

* Clinical data were obtained from The Cancer Genome Atlas (TCGA).
