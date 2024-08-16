"""
Main file. Generates problem instances and solves them by the sBB algorithm as well as Gurobi's built-in PLF and global solvers.
Logarithmic models are impmemented in Julia (there exist a package doing this).
"""

#import statements
import numpy as np
import time 
import bisect
import random
import csv
from gurobipy import *
import io #Necessary for suppressing print
from contextlib import redirect_stdout #Necessary for suppressing print
import sBB_functions as sBB
from sBB_main import *
from instance_generation import *
import gurobi_solver

#-------------------------------
# Set parameters for computations
#-------------------------------
# Problem type
problem = "network flow" #knapsack, concave-knapsack, global-knapsack, network flow, discontinuous network flow 

# Dimension of problems
n_list =  [100] # [5,10,15,20] for global-knapsack, [100] else 

# Number of segments which should be considered
K_list = [10, 100, 500, 1000, 5000, 10000] # [10000] for global-knapsack, [10, 100, 500, 1000, 5000, 10000] else

# Number of random instances
num_instances = 50 # 50 for all problems

# Set the seed for reproducibility
random.seed(1)  

# Timelimit in s
timelimit = 1800 # 1800

# Relative optimality gap
epsilon = 0.00001 

#-------------------
# "Heat up" solver
#-------------------
#Suppressing Gurobi output
trap = io.StringIO() # set a trap and redirect stdout
with redirect_stdout(trap):# set a trap and redirect stdout
    m = Model()
    x = m.addVar(name = 'x', lb = 0)
    y = m.addVar(name = 'y', lb=0, ub=3)
    m.setObjective(12*x+20*y, GRB.MINIMIZE)
    m.addConstr(6*x+8*y >= 100)
    m.addConstr(7*+12*y >= 120)
    m.optimize()

#-------------------
# Run computations
#------------------- 

# Iterate over dimension of the problem and number of segments

for n in n_list:
    
    for K in K_list:
        
        #open csv
        file = f"{problem}_K={K}.csv" if problem != "global-knapsack" else f"{problem}_n={n}.csv"
        path = "Julia-MIP/"+file
        f = open(path, 'w', newline='')
        writer = csv.writer(f) 
        
        #Generate random instance
        for instance in range(1,num_instances+1):
            
            print("")
            print("--------------------------------------------")
            print(f"n={n}, K={K}, instance={instance}")
            print("--------------------------------------------")   
            
            #Knapsack problem
            if problem == "knapsack":
                func_numbers = [random.randint(1,20) for i in range(0,n)] 
                PLFs_breakpoints_x, PLFs_breakpoints_y = getPLFs(n, K, func_numbers, False)
                rhs = getRHS(n, K, problem, PLFs_breakpoints_x)
            
            #Concave Knapsack problem
            elif problem == "concave-knapsack":
                func_numbers = [random.randint(1,20) for i in range(0,n)]
                PLFs_breakpoints_x, PLFs_breakpoints_y = getPLFs(n, K, func_numbers, True)
                rhs = getRHS(n, K, problem, PLFs_breakpoints_x)
            
            # Knapsack problem solved by Gurobi's sBB solver    
            elif problem == "global-knapsack":
                func_numbers = [random.choice([2,9,11,12,13,14,15,20]) for i in range(0,n)] # restricted set of functions due to Gurobi's limited function support
                PLFs_breakpoints_x, PLFs_breakpoints_y = getPLFs(n, K, func_numbers, False)
                rhs = getRHS(n, K, problem, PLFs_breakpoints_x)
            
            #network flow problem
            elif problem == "network flow":
                PLFs_breakpoints_x, PLFs_breakpoints_y = getNetworkPLFs(n, K)
                rhs = getRHS(n, K, problem, PLFs_breakpoints_x)
                
            # Discontinuous network flow problem
            elif problem == "discontinuous network flow":
                PLFs_breakpoints_x, PLFs_breakpoints_y = getDiscontinuousNetworkPLFs(n, K)
                rhs = getRHS(n, K, problem, PLFs_breakpoints_x)
                
            else:
                print("Problem not recogonised.") 
                
            
            #-------------------
            # Solve by sBB method
            #--------------------
                
            #Solve problem by sBB
            (sBB_time, grb_time, lp_time, env_time, plf_time, sBB_nodes, sBB_ub, sBB_lb, sBB_bestpoint, sBB_root_node_lb) = spatialBB(PLFs_breakpoints_x,PLFs_breakpoints_y, n, K, problem, rhs, timelimit, epsilon)
        
            #If generated instance is infeasible, generate a new instance (can only happen in corner cases)
            if sBB_bestpoint == 0:
                print("infeasible instance")
                instance = instance - 1
                continue
            
            # Print runtimes, lower and upper bounds of this instance
            print("sBB")
            print(f"lb: {sBB_lb}, ub: {sBB_ub}")
            print("time: ", sBB_time)
            
            #--------------------------------------------------
            # Solve by other approaches depending on the problem
            #--------------------------------------------------
            
            if problem == "global-knapsack":
                
                #Solve problem by Gurobi sB&B 
                node_limit = "none"
                (GRB_sBB_time, GRB_sBB_ub, GRB_sBB_lb, GRB_sBB_incumbent) = gurobi_solver.gurobi_sBB(func_numbers, n, K, problem, rhs, timelimit, epsilon, node_limit)
                
                #Solve problem by Gurobi sB&B but stop after root node to retrieve root node relaxation value
                node_limit = 0
                (GRB_sBB_time_root, GRB_sBB_ub_root, GRB_sBB_lb_root, GRB_sBB_incumbent_root) = gurobi_solver.gurobi_sBB(func_numbers, n, K, problem, rhs, timelimit, epsilon, node_limit)
                
                
                # Print runtimes, lower and upper bounds of this instance
                print("Gurobi sBB")
                print(f"lb: {GRB_sBB_lb}, ub: {GRB_sBB_ub}")
                print("time: ", GRB_sBB_time)
                
                
                # Print root node comparison
                print(f"root lb sBB: {sBB_root_node_lb}, root lb GRB sBB: {GRB_sBB_lb_root}")
                
                #Write problem type and solution times in the first row
                writer.writerow([problem,n,K, sBB_time, GRB_sBB_time, sBB_root_node_lb, GRB_sBB_lb_root])
                
                #Write selected functions in second row
                writer.writerow(func_numbers)
                
                
            elif problem == "discontinuous network flow":
                
                #Solve problem by Gurobi PLF
                (GRB_PLF_time, GRB_PLF_ub, GRB_PLF_lb, GRB_PLF_incumbent) = gurobi_solver.gurobi_PLF(PLFs_breakpoints_x, PLFs_breakpoints_y, n, K, problem, rhs, timelimit, epsilon, "none")    
        
                # Print runtimes, lower and upper bounds of this instance
                print("Gurobi PLF")
                print(f"lb: {GRB_PLF_lb}, ub: {GRB_PLF_ub}")
                print("time: ", GRB_PLF_time)
                
                #Write problem type and solution times in the first row
                writer.writerow([problem,n,K,sBB_time, GRB_PLF_time])
                
            else:
                
                #Solve problem by Gurobi PLF
                (GRB_PLF_time, GRB_PLF_ub, GRB_PLF_lb, GRB_PLF_incumbent) = gurobi_solver.gurobi_PLF(PLFs_breakpoints_x, PLFs_breakpoints_y, n, K, problem, rhs, timelimit, epsilon, "none")    
        
                # Print runtimes, lower and upper bounds of this instance
                print("Gurobi PLF")
                print(f"lb: {GRB_PLF_lb}, ub: {GRB_PLF_ub}")
                print("time: ", GRB_PLF_time)
                
                #Write problem type and solution times in the first row
                writer.writerow([problem,n,K,sBB_time,grb_time,lp_time,env_time,plf_time,GRB_PLF_time])
            
            
            #-------------------------------------------------------------------------------------------------
            # Write problem data to csv (might be needed to hand over to Julia for logarithmic PLF models)
            #-----------------------------------------------------------------------------------------------
            
            #Write right-hand-side values in the second row
            if isinstance(rhs, (np.floating, float)):
                rhs = [rhs] #convert to list if rhs is no list
            writer.writerow(rhs)
            
            #Write PLF breakpoints x values
            for i in range(0,n):
                writer.writerow(PLFs_breakpoints_x[i])
            
            #Write PLF breakpoints y values
            for i in range(0,n):
                writer.writerow(PLFs_breakpoints_y[i])

        f.close()
    
print("Computations are done.")

##################################################################