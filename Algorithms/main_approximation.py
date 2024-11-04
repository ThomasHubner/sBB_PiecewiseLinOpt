"""
Computes the benfit of approximating a nonliner function by more segments.
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
problem = "knapsack" # only knapsack uses an approximation of nonlinear functions

# Dimension of problems
n = 100 # use n=100 variables 

# Number of segments which should be considered
K_list = [10, 20, 50, 100, 500, 1000, 5000, 10000] 

# Number of random instances
num_instances = 50 #50

# Set the seed for reproducibility
random.seed(1)  

# Timelimit in s
timelimit = 1800 #1800

# Relative optimality gap
epsilon = 0.00001 


#-------------------
# Run computations
#------------------- 

#open csv
file = f"approximation_quality.csv" 
f = open(file, 'w', newline='')
writer = csv.writer(f)  

#Generate random instance
instance = 0
while instance < num_instances:
    instance += 1
    
    #list to store objective function values
    optimal_values = []
    
    #Check whether instance can be solved for K=10,000
    K=10000

    #Knapsack problem
    if problem == "knapsack":
        func_numbers = [random.randint(1,20) for i in range(0,n)] 
        PLFs_breakpoints_x, PLFs_breakpoints_y = getPLFs(n, K, func_numbers, False)
        rhs = getRHS(n, K, problem, PLFs_breakpoints_x)
    
    #Only knapsack problem is considered
    else:
        print("Problem not recogonised. Only knpasack problem uses approximated nonlinear functions.") 
          
    #Solve problem by sBB
    (sBB_time, grb_time, lp_time, env_time, plf_time, sBB_nodes, sBB_ub, sBB_lb, sBB_bestpoint, sBB_root_node_lb) = spatialBB(PLFs_breakpoints_x,PLFs_breakpoints_y, n, K, problem, rhs, timelimit, epsilon)

    #If generated instance is infeasible or cannot be solved for K=10000, generate a new instance 
    if sBB_bestpoint == 0:
        print("infeasible instance")
        instance = instance - 1
        continue
    elif sBB_time >= timelimit:
        print("unsolvable instance")
        instance = instance - 1
        continue

    #--------------------------------
    # Iterate over number of segments
    #--------------------------------
    for K in K_list:

        print("")
        print("--------------------------------------------")
        print(f"n={n}, K={K}, instance={instance}")
        print("--------------------------------------------")   
    
        #Get refined approximation
        PLFs_breakpoints_x, PLFs_breakpoints_y = getPLFs(n, K, func_numbers, False)
        
        #Solve problem by sBB
        (sBB_time, grb_time, lp_time, env_time, plf_time, sBB_nodes, sBB_ub, sBB_lb, sBB_bestpoint, sBB_root_node_lb) = spatialBB(PLFs_breakpoints_x,PLFs_breakpoints_y, n, K, problem, rhs, timelimit, epsilon)
    
        # Print runtimes, lower and upper bounds of this instance
        print("sBB")
        print(f"lb: {sBB_lb}, ub: {sBB_ub}")
        print("time: ", sBB_time)
        
        # Optimal value of original not-approximated problem
        # Insert solution found by sBB in the real nonlinear objective function
        # As objective function is separable, this is computed by a for loop
        optimal_value = 0
        for i in range(n):
            num = func_numbers[i] # retrieve function number
            function = globals()["func"+str(num)] # retrieve nonlinear function
            value = function(sBB_bestpoint[i]) # value of the nonlinear function at the point
            optimal_value += value
    
        # Add value to list
        optimal_values.append(optimal_value)
    
    # Get relative improvement of the approximation refinement
    rel_improvments = []
    for i in range(len(optimal_values) - 1):
        diff = optimal_values[i + 1] - optimal_values[i] # absolute improvment
        rel_improvment = diff / optimal_values[i] # relative improvement
        rel_improvments.append(rel_improvment)
    
    # Write relative improvement of this instance to row
    writer.writerow(rel_improvments)
 
f.close()
    
print("Computations are done.")

##################################################################