"""
Spatial branch-and-bound algorithm to solve the piecewise linear optimization problem.
"""

import numpy as np
from gurobipy import *
import io #Necessary for suppressing print
from contextlib import redirect_stdout #Necessary for suppressing print
import time 
import bisect as bisect
import sBB_functions as sBB

def spatialBB(PLFs_breakpoints_x, PLFs_breakpoints_y,
              n, K, problem, rhs,timelimit, epsilon):
    '''
    sBB algorithm for global solution of problems with PLF.

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

    Returns
    -------
    
    sol_time: float
        total solution time
    grb_time: float
        total time spend at the gurobi-python interface to build the LPs
    lp_time: float
        total time needed to solve the LPs
    env_time: float
        total time needed to construct the lower convex envelopes
    node_count: integer
        number of nodes explored
    upper_bound: float
        upper bound
    lower_bound: float
        lower bound
    incumbent: array
        solution point
    '''
    
    #Getting time
    sbb_starttime = time.time()
    
    #Extracting lower and upper bounds for variables
    lower_bounds = []
    upper_bounds = []
    for breakpoints in PLFs_breakpoints_x:
        lower_bounds.append(breakpoints[0])
        upper_bounds.append(breakpoints[-1])
        
    #Getting lower convex envelopes for each PLF
    env_breakpoints_x = []
    env_breakpoints_y = []
    env_time = 0
    for i in range(0,n):
        lower_x, lower_y, env_runtime = sBB.getEnvelope(PLFs_breakpoints_x[i], PLFs_breakpoints_y[i],
                                       [lower_bounds[i],upper_bounds[i]]) 
        env_breakpoints_x.append(lower_x)
        env_breakpoints_y.append(lower_y)
        env_time += env_runtime
        
    #Suppressing Gurobi output
    trap = io.StringIO() # set a trap and redirect stdout
    with redirect_stdout(trap):# set a trap and redirect stdout
        
        grb_starttime = time.time()
        
        #Generating model
        m = Model("P")
    
        #Variables
        x = m.addMVar(shape=n, lb=lower_bounds, ub=upper_bounds,
                      vtype = GRB.CONTINUOUS, name = "x")
        #Epigraph variables
        e = m.addMVar(shape=n, lb=-GRB.INFINITY,
                  vtype = GRB.CONTINUOUS, name = "e")
        m.update()
    
        #Constraints
        if problem == "knapsack" or problem == "concave-knapsack":
            m.addConstr(np.ones(n)@x == rhs)
        elif problem == "network flow":
            nr = int(np.sqrt(n))
            m.addConstrs(quicksum([x[j] for j in range(i*nr,(i+1)*nr)])
                         -quicksum([x[nr*j+i] for j in range(0,nr)])== rhs[i]
                         for i in range(0,nr))
        else:
            print("Problem not recogonised. Run again with either 'knapsack' or 'network flow'")  
            return 
        m.update()
        
        #Objective function
        m.setObjective(np.ones(n)@e, GRB.MINIMIZE)
        for i in range(0,n):
            lower_x = env_breakpoints_x[i]
            lower_y = env_breakpoints_y[i]
            for k in range(1,len(lower_x)):
                slope = (lower_y[k] - lower_y[k-1]) / (lower_x[k] - lower_x[k-1])
                m.addConstr(e[i] >= slope * (x[i] - lower_x[k-1]) + lower_y[k-1]) 
        m.update()
        
        grb_runtime = time.time() - grb_starttime
        
        #Solver- Computing bounds
        m.optimize() 
        
        #Check if feasible
        if m.status !=2:
            print("Problem is infeasible!")
            return (0,0,0,0,0,0,0,0,0)
    
    #Getting time 
    lp_time = m.Runtime
    grb_time = grb_runtime
    
    #Calculating real PLF values (upper bound)
    plf_starttime = time.time()
    plf_values = [sBB.getPLFvalue(PLFs_breakpoints_x[i], PLFs_breakpoints_y[i] , m.x[i]) 
               for i in range(0,n)]
    plf_time = time.time() - plf_starttime
    
    #Intializing global upper bound and incumbent
    global_ub = sum(plf_values)
    incumbent = m.x[0:n]
    
    #List of nodes, i.e. the model and the envelopes at these nodes
    nodes_lb = [m.objval]
    nodes_model = [m]
    nodes_envelope_x = [env_breakpoints_x]
    nodes_envelope_y = [env_breakpoints_y]
    nodes_plf_values = [plf_values] #saved to avoid re-evaluation
    
    #Node count
    node_count = 1
    
    #------------------------------------
    # Starting branch-and-bound iteration
    #------------------------------------ 
    
    while( len(nodes_lb) != 0 
          and (global_ub - nodes_lb[0]) / abs(global_ub) > epsilon 
          and time.time() - sbb_starttime < timelimit):
        
        
        #------------------------------------------------------
        #Best first search: first node of sorted list is selected
        #-------------------------------------------------------
        m = nodes_model[0]
        env_breakpoints_x = nodes_envelope_x[0]
        env_breakpoints_y = nodes_envelope_y[0]
        plf_values = nodes_plf_values[0]


        #---------------------------------------------------------------------
        #Branching and Bounding
        #--------------------------------------------------------------------
        child_nodes, comp_times = sBB.branch_and_bound(m,n, plf_values, env_breakpoints_x, env_breakpoints_y,
                                                    PLFs_breakpoints_x, PLFs_breakpoints_y)
        
        #Extract information
        left_model, left_env_breakpoints_x,left_env_breakpoints_y, left_plf_values = child_nodes[0]
        right_model,right_env_breakpoints_x,right_env_breakpoints_y,right_plf_values = child_nodes[1] 
        
        #Increase time
        env_time += comp_times[0]
        grb_time += comp_times[1]
        lp_time += comp_times[2]
        plf_time += comp_times[3]    
        
        #------------------------------------------------------------------------------------------------    
        # Branching&Bounding ends here: child nodes are determined and solved -----------------------------------
        #------------------------------------------------------------------------------------------------    

        #Upper bounds
        ub_left = sum(left_plf_values)
        ub_right = sum(right_plf_values)
        
        #Updating global upper bound
        if ub_left < global_ub and ub_left < ub_right:
            global_ub = ub_left
            incumbent = left_model.x[0:n]
        elif ub_right < global_ub:
            global_ub = ub_right
            incumbent = right_model.x[0:n]     
            
        #Removing parent node
        nodes_lb.pop(0)
        nodes_model.pop(0)
        nodes_envelope_x.pop(0)
        nodes_envelope_y.pop(0)
        nodes_plf_values.pop(0)
        
        #Inserting new nodes in list of nodes
        if left_model.objval < global_ub:
            left_pos = bisect.bisect(nodes_lb,left_model.objval)
            nodes_lb.insert(left_pos, left_model.objval)
            nodes_model.insert(left_pos, left_model)
            nodes_envelope_x.insert(left_pos, left_env_breakpoints_x)
            nodes_envelope_y.insert(left_pos, left_env_breakpoints_y)
            nodes_plf_values.insert(left_pos, left_plf_values)
        
        if right_model.objval < global_ub:
            right_pos = bisect.bisect(nodes_lb,right_model.objval)
            nodes_lb.insert(right_pos, right_model.objval)
            nodes_model.insert(right_pos, right_model)
            nodes_envelope_x.insert(right_pos, right_env_breakpoints_x)
            nodes_envelope_y.insert(right_pos, right_env_breakpoints_y)
            nodes_plf_values.insert(right_pos, right_plf_values)
        
        #Fathoming/Pruning
        del_pos = bisect.bisect_left(nodes_lb, global_ub)
        nodes_lb = nodes_lb[0:del_pos]
        nodes_model = nodes_model[0:del_pos]
        nodes_envelope_x = nodes_envelope_x[0:del_pos]
        nodes_envelope_y = nodes_envelope_y[0:del_pos]
        nodes_plf_values = nodes_plf_values[0:del_pos]
        
        node_count += 2
        
    return (time.time() - sbb_starttime, grb_time, lp_time, env_time, plf_time, node_count,
            global_ub, nodes_lb[0] if len(nodes_lb)!=0 else global_ub, incumbent)
            
    

