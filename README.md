# sBB_PiecewiseLinOpt
Spatial branch-and-bound (sBB) algorithm for piecewise-linear optimization problems.

This repository is an accompaniment to the paper titled *Spatial branch-and-bound for nonconvex piecewise-linear optimization* by Thomas Hübner.
It contains the code used to conduct computational experiments. 
Details and results can be found in Section 4 and the Appendix A3 of the above mentioned paper.

The code can be found in the folder algorithms. It contains:
  - main.py
  - instance_generation.py
  - sBB_main.py
  - sBB_functions.py
  - folder "Julia-MIP" which contains
    - main.jl
    - MIP_solver.jl
    - folder "results"

The code works as follows. 
The file *main.py* uses functions contained in *instance_generation.py* to build instances of the knapsack or network flow problems.
These instances are then handed over to the sBB algorithm in *sBB_main.py* which uses functions from *sBB_functions.py* to solve the problem.
The problem information together with the solve time of the sBB are then written in an csv-file located in the folder *Julia-MIP*.
The code in *main.jl* reads the csv-file and hands the problem information to a function in *MIP_solver.jl* which solves the problem.
Finally, the solve times of all used methods are written in a csv-file in the folder *results*.


