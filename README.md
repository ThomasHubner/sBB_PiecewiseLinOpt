# sBB_PiecewiseLinOpt
Spatial branch-and-bound algorithm for piecewise-linear optimization problems.

This repository is an accompaniment to the paper titled *Spatial branch-and-bound for nonconvex separable piecewise-linear optimization* by Thomas Hübner, Akshay Gupte and Steffen Rebennack.
It contains the code used to conduct computational experiments. 
Details and results can be found in Section 5 of the above-mentioned paper.

The code can be found in the folder *Algorithms*. It contains:
  - main.py
  - instance_generation.py
  - sBB_main.py
  - sBB_functions.py
  - gurobi_solver.py
  - main_approximation.py
  - folder "Julia-MIP" which contains
    - main.jl
    - MIP_solver.jl
    - folder "results"
      - writeLatexTable.py

The code works as follows. 
The file *main.py* uses functions contained in *instance_generation.py* to build random instances of the knapsack or network flow problems.
These instances are then handed over to the sBB algorithm in *sBB_main.py*, which uses functions from *sBB_functions.py* to solve the problem.
Alternatively, those random instances are solved by Gurobi's built-in piecewise-linear function handler or Gurobi's built-in MINLP solver. Both are called by the functions within *gurobi_solver.py*.
The problem information and the solve time of the sBB and Gurobis built-in methods are then written in a CSV file in the folder *Julia-MIP*.
The code in *main.jl* reads the CSV file and hands the problem information to a function in *MIP_solver.jl*, which solves the problem by using predefined logarithmic MIP models.
Finally, the solve times of all methods used are written in a CSV file in the folder *results*.
