"""
Main file. Contains the computations performed with the sBB.
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

#-------------------------------
#Set parameters for computations
#-------------------------------

n = 100 #Dimension of problems
K_list = [10, 100, 500, 1000, 5000, 10000] #Number of segments which should be considered
num_instances = 100 #Number of instances solved for each K
problem = "network flow" #knapsack, concave-knapsack, network flow 
timelimit = 1800 #timelimit
epsilon = 0.00001 #relative optimality gap

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

#Start computations
for K in K_list:
    
    #open csv
    file = f"problem_info_K={K}.csv"
    path = "Julia-MIP/"+file
    f = open(path, 'w', newline='')
    writer = csv.writer(f) 
    
    #Generate random instance
    for instance in range(1,num_instances+1):
        
        feasible = False
        
        #Check thtat feasible instance was generated
        while not feasible:
            
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
            
            #network flow problem
            elif problem == "network flow":
                PLFs_breakpoints_x, PLFs_breakpoints_y = getNetworkPLFs(n, K)
                rhs = getRHS(n, K, problem, PLFs_breakpoints_x)
                
            else:
                print("Problem not recogonised. Run again with either 'knapsack' or 'network flow'") 
                
            #Solve problem
            (sBB_time, grb_time, lp_time, env_time, plf_time, sBB_nodes, sBB_ub, sBB_lb,
             sBB_bestpoint) = spatialBB(PLFs_breakpoints_x,PLFs_breakpoints_y,
                                        n, K, problem, rhs, timelimit, epsilon)
        
            #Check if feasible
            if sBB_bestpoint == 0:
                print("infeasible instance")
            else:
                feasible = True
                
        
        #Write problem type in the first row
        writer.writerow([problem,n,K,sBB_time,grb_time,lp_time,env_time,plf_time])
        
        #Write rhs values in the second row
        if isinstance(rhs, (np.floating, float)):
            rhs = [rhs]
        writer.writerow(rhs)
        
        #Write PLF breakpoints x values
        for i in range(0,n):
            writer.writerow(PLFs_breakpoints_x[i])
        
        #Write PLF breakpoints y values
        for i in range(0,n):
            writer.writerow(PLFs_breakpoints_y[i])
            
        print(f"lb: {sBB_lb}, ub: {sBB_ub}")
        print("time: ", sBB_time)

    f.close()
print("Computations are done.")

##################################################################