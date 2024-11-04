"""
Contains functions which call the built-in Gurobi solvers to solve the test problems.
    - gurobi_PLF (Gurobi's built-in piecewise-linear function solver.)
    - gurobi_sBB (Gurobi's built-in spatial B&B solver for nonlinear nonconvex functions (other than PLFs).
"""

import numpy as np
from gurobipy import *
import time 


def gurobi_PLF(PLFs_breakpoints_x, PLFs_breakpoints_y, n, K, problem, rhs, timelimit, epsilon, node_limit):
    '''
    Gurobi's built-in piecewise-linear function solver.

    Parameters
    ----------
    PLFs_breakpoints_x : list
        A list of the x-values of the breakpoints for each PLF.
    PLFs_breakpoints_y : list
        A list of the y-values of the breakpoints for each PLF.
    n : integer
        Dimension of the problem.
    K : integer
        Number of segments of each PLF.
    problem : string
        Problem which is to be solved: knapsack or network flow
    rhs : list of floats
        the right hand side values of the problem
    timelimit : float
        Time limit for solver.
    epsilon : float
        sufficient relative optimality gap for termination
    node_limit : integer
        max number of nodes to be explored (we can set this to 0 to retrieve the root node relaxation)

    Returns
    -------
    
    sol_time: float
        total solution time
    upper_bound: float
        upper bound
    lower_bound: float
        lower bound
    incumbent: array
        solution point
    '''
    
    #Extracting lower and upper bounds for variables
    lower_bounds = []
    upper_bounds = []
    for breakpoints in PLFs_breakpoints_x:
        lower_bounds.append(breakpoints[0])
        upper_bounds.append(breakpoints[-1])
    
    #Starttime
    starttime = time.time()
    
    #Generating model
    m = Model("P")
    
    # Disable gurobi log output to console:: 0(no), 1(yes) 
    m.Params.LogToConsole = 0
    
    # Set MIP gap
    m.Params.MIPGap = epsilon
    
    # Set timelimit
    m.Params.TimeLimit = timelimit
    
    # Node limit
    if node_limit != "none":
        m.Params.NodeLimit = node_limit

    #Variables
    x = m.addMVar(shape=n, lb=lower_bounds, ub=upper_bounds,
                  vtype = GRB.CONTINUOUS, name = "x")

    #Constraints
    if problem == "knapsack" or problem == "concave-knapsack" or problem == "global-knapsack":
        m.addConstr(np.ones(n)@x == rhs)
    elif problem == "network flow" or problem == "discontinuous network flow":
        nr = int(np.sqrt(n))
        m.addConstrs(quicksum([x[j] for j in range(i*nr,(i+1)*nr)])
                     -quicksum([x[nr*j+i] for j in range(0,nr)])== rhs[i]
                     for i in range(0,nr))
    else:
        print("Problem not recogonised.")  
        return 
    m.update()
    
    #Objective function
    if problem == "discontinuous network flow": # Add fixed charge jump (y=0 at x=0)
        for i in range(n):
            m.setPWLObj(m.getVarByName(f'x[{i}]'), [0] + PLFs_breakpoints_x[i], [0] + PLFs_breakpoints_y[i])
    else:
        for i in range(n):
            m.setPWLObj(m.getVarByName(f'x[{i}]'), PLFs_breakpoints_x[i], PLFs_breakpoints_y[i])
    m.update()
    
    #Solver
    m.optimize() 
    
    #Getting time 
    sol_time = time.time() - starttime
    
    # Getting upper bound
    try:
        upper_bound = m.ObjVal
    except:
        upper_bound = 1e09
    
    # Getting lower bound
    lower_bound = m.ObjBound
    
    # Getting solution
    try:
        incumbent = m.X
    except:
        incumbent = [0 for i in range(n)]
        
    # If GRB reports problem to be infeasible return 1800 (problem is feasible by construction but GRB encounters numerical issues)
    if m.Status == 3:
        sol_time = 1800
        
    return sol_time, upper_bound, lower_bound, incumbent


#---------------------------------------------------------------------------------------------------------------


def gurobi_sBB(func_numbers, n, K, problem, rhs,timelimit, epsilon, node_limit):
    '''
    Gurobi's built-in spatial B&B solver for nonlinear nonconvex functions (other than PLFs).
    Can be used to "solve" knapsack problems (only 8 out of the 20 given functions are implemented).
    FOther functions cause numerical issues and can lead to infeasibility.

    Parameters
    ----------
    func_numbers : list
        List of integers which define which function should be slected
    n : integer
        Dimension of the problem.
    K : integer
        Number of segments of each PLF.
    problem : string
        Problem which is to be solved: knapsack or network flow
    rhs : list of floats
        the right hand side values of the problem
    timelimit : float
        Time limit for solver.
    epsilon : float
        sufficient relative optimality gap for termination
    node_limit : integer
        max number of nodes to be explored (we can set this to 0 to retrieve the root node relaxation)

    Returns
    -------
    
    sol_time: float
        total solution time
    upper_bound: float
        upper bound
    lower_bound: float
        lower bound
    incumbent: array
        solution point
    '''
    
    # List of function variable bounds
    function_bounds = [ [-5,5], [-5,5], [-5,5], [-10,10], [0.1,10], [0.1,10], [-10,10], [-5,5], [-4,4], [-10,10], [-1,7], [-4,4], [0,5], [-1,4], [-4,4], [-5,5], [-10,10], [-1.5,3], [-5,5], [0.1, 10]]

    #Extracting lower and upper bounds for variables
    lower_bounds = []
    upper_bounds = []
    for i in func_numbers:
        lower_bounds.append(function_bounds[i-1][0])
        upper_bounds.append(function_bounds[i-1][1])
    
    #Starttime
    starttime = time.time()
    
    #Generating model
    m = Model("P")
    
    # Disable gurobi log output to console:: 0(no), 1(yes) 
    m.Params.LogToConsole = 0 
    
    # Set MIP gap
    m.Params.MIPGap = epsilon
    
    # Set timelimit
    m.Params.TimeLimit = timelimit
    
    # Set solution method to sB&B 
    m.Params.FuncNonlinear = 1
    
    # Node limit
    if node_limit != "none":
        m.Params.NodeLimit = node_limit

    #Variables
    x = m.addMVar(shape=n, lb=lower_bounds, ub=upper_bounds,
                  vtype = GRB.CONTINUOUS, name = "x")
    
    #Epigraph variables - Gurobi supports general nonconvex functions only in constrains
    e = m.addMVar(shape=n, lb=-GRB.INFINITY,
              vtype = GRB.CONTINUOUS, name = "e")
    m.update()

    #Constraints
    if problem == "global-knapsack":
        m.addConstr(np.ones(n)@x == rhs)
    else:
        print("Only knapsack problem is supported.")  
        return 
    m.update()
    
    # Objective function
    m.setObjective(np.ones(n)@e, GRB.MINIMIZE)
    
    # Add epigraph constraints
    for i in range(n):
        x_i = m.getVarByName(f'x[{i}]')
        e_i = m.getVarByName(f'e[{i}]')
        if func_numbers[i] == 2: # -0.2*exp(-x) + x**2
            a1 = m.addVar(lb=-GRB.INFINITY, vtype = GRB.CONTINUOUS)
            b1 = m.addVar(lb=-GRB.INFINITY, vtype = GRB.CONTINUOUS)
            c1 = m.addVar(lb=-GRB.INFINITY, vtype = GRB.CONTINUOUS)
            m.addConstr(a1 == -x_i)
            m.addGenConstrExp(a1,b1) # b1 = exp(a1)
            m.addGenConstrPoly(x_i,c1,[1,0,0])
            m.addConstr(-0.2*b1 + c1 <= e_i)
        elif func_numbers[i] == 9: #-x**7 / 5040 + x**5 / 120 - x**3 / 3 + x
            a1 = m.addVar(lb=-GRB.INFINITY, vtype = GRB.CONTINUOUS)
            m.addGenConstrPoly(x_i,a1,[-1/5040, 0, 1/120, 0, -1/3, 0, 1, 0])
            m.addConstr(a1 <= e_i)
        elif func_numbers[i] == 11: # x**4 - 12*x**3 + 47*x**2 - 60*x
            a1 = m.addVar(lb=-GRB.INFINITY, vtype = GRB.CONTINUOUS)
            m.addGenConstrPoly(x_i,a1,[1, -12, 47, -60, 0])
            m.addConstr(a1 <= e_i)
        elif func_numbers[i] == 12: # x**6 -15*x**4 + 27*x**2 +250
            a1 = m.addVar(lb=-GRB.INFINITY, vtype = GRB.CONTINUOUS)
            m.addGenConstrPoly(x_i,a1,[1, 0, -15, 0, 27, 0, 250])
            m.addConstr(a1 <= e_i)
        elif func_numbers[i] == 13: # x**4 - 10*x**3 +35*x**2 -50*x + 24
            a1 = m.addVar(lb=-GRB.INFINITY, vtype = GRB.CONTINUOUS)
            m.addGenConstrPoly(x_i,a1,[1, -10, 35, -50, 24])
            m.addConstr(a1 <= e_i)        
        elif func_numbers[i] == 14: # 0.2 *x**5 -1.25*x**4 +2.33*x**3 -2.5*x**2 +6*x
            a1 = m.addVar(lb=-GRB.INFINITY, vtype = GRB.CONTINUOUS)
            m.addGenConstrPoly(x_i,a1,[0.2, -1.25, 2.33, -2.5, 6, 0])
            m.addConstr(a1 <= e_i)  
        elif func_numbers[i] == 15: # x**3 -7*x + 7
            a1 = m.addVar(lb=-GRB.INFINITY, vtype = GRB.CONTINUOUS)
            m.addGenConstrPoly(x_i,a1,[1, 0, -7, 7])
            m.addConstr(a1 <= e_i) 
        elif func_numbers[i] == 20: # 1/x + 2*log(x) -2
            a1 = m.addVar(lb=-GRB.INFINITY, vtype = GRB.CONTINUOUS)
            b1 = m.addVar(lb=-GRB.INFINITY, vtype = GRB.CONTINUOUS)
            m.addGenConstrPow(x_i,a1,-1) # a1 = 1/x
            m.addGenConstrLog(x_i,b1)
            m.addConstr(a1 + 2*b1 -2 <= e_i)
        else:
            print("Function not supported by Gurobi sBB solver.")
                 
        m.update()
    
    #Solver
    m.optimize() 
    
    #Getting time 
    sol_time = time.time() - starttime
    
    # Getting upper bound
    try:
        upper_bound = m.ObjVal
    except:
        upper_bound = 1e09
    
    # Getting lower bound
    lower_bound = m.ObjBound
    
    # Getting solution
    try:
        incumbent = m.X
    except:
        incumbent = [0 for i in range(n)]
    
    return sol_time, upper_bound, lower_bound, incumbent

