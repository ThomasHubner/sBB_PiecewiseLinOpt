"""
Defines auxiliary functions used in the sBB. It includes:
    - getPLFvalue (computes the value of a PLF at  given point)
    - getEnvelope (computes the lower convex envelope)
    - branch_and_bound (branches current parent node and bounds child nodes)
    - solve_node (solves the child node, is a subfunction of branch_and_bound)
"""

import numpy as np
from gurobipy import *
import io #Necessary for suppressing print
from contextlib import redirect_stdout #Necessary for suppressing print
import time 
import bisect as bisect


#--------------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------
# Function to return the value of a PLF at a certain point x
#---------------------------------------------------------

def getPLFvalue(breakpoints_x, breakpoints_y, point, problem):
    '''
    Returns the value of the PLF defined by "breakpoints_x" (list) and "breakpoints_y" (list) at the point "point"(float).
    If problem = discontinuous network flow, it considers the fixed-charge jump at point 0.
    '''
    
    
    if point <= breakpoints_x[0]:
        return 0 if problem == "discontinuous network flow" else breakpoints_y[0] # fixed charge jump in discontinuous network flow problem
    elif point >= breakpoints_x[-1]:   
        return breakpoints_y[-1]
    else:
        pos = bisect.bisect(breakpoints_x, point)
        value = ((breakpoints_y[pos] - breakpoints_y[pos-1]) 
                 / (breakpoints_x[pos]-breakpoints_x[pos-1]) 
                 ) * (point-breakpoints_x[pos-1]) + breakpoints_y[pos-1]
    
        return value
    
    
#--------------------------------------------------------------------------------------------------------------------------------
    
    
#--------------------------------------------------------
# Function which computes the lower convex envelope of a PLF
#-------------------------------------------------------

def getEnvelope(breakpoints_x, breakpoints_y, interval, problem):
    '''
    Computes the lower convex envelope of a PLF over a given interval.
    If problem = discontinuous network flow, it considers the fixed-charge jump at point 0.
    
    Parameters
    ----------
    breakpoints_x : list
        X-value of breakpoints of PLF.
    breakpoints_y : list
        Y-value of breakpoints of PLF.
    interval : array
        Interval over which the envelope should be calculated.
    problem : string
        Test problem to be solved.

    Returns
    -------
    lower_x : list
        X-value of breakpoints of lower Envelope.
    lower_y : list
        Y-value of breakpoints of lower Envelope.
    runtime: float
        runtime of the algorithm
    '''
    
    # Measure time
    starttime = time.time()
    
    # Get copies of lists
    breakpoints_x = breakpoints_x.copy()
    breakpoints_y = breakpoints_y.copy()
    
    #1) Get PLF over sub-interval
    
    #Interval is full range of the functions domain
    if interval[0]==breakpoints_x[0] and interval[1]==breakpoints_x[-1]:
        points_x = breakpoints_x
        points_y = breakpoints_y
    
    #Not the full interval is considered
    else:
    
        #Get y value of the interval edges
        pos_left = bisect.bisect_right(breakpoints_x,interval[0])
        pos_right = bisect.bisect_left(breakpoints_x,interval[1])
        interval_y = np.array([
            ((breakpoints_y[pos_left] - breakpoints_y[pos_left-1]) /
            (breakpoints_x[pos_left] - breakpoints_x[pos_left-1])) * 
            (interval[0] - breakpoints_x[pos_left-1] ) + breakpoints_y[pos_left-1],
            ((breakpoints_y[pos_right] - breakpoints_y[pos_right-1]) /
            (breakpoints_x[pos_right] - breakpoints_x[pos_right-1])) * 
            (interval[1] - breakpoints_x[pos_right-1] ) + breakpoints_y[pos_right-1]
            ])
        
        #Remove all points which are not in interval and add interval edges
        points_x = breakpoints_x[pos_left:pos_right]
        points_y = breakpoints_y[pos_left:pos_right]
        points_x.insert(0,interval[0])
        points_x.append(interval[1])
        points_y.insert(0,interval_y[0])
        points_y.append(interval_y[1])
        
    #1a) If fixed-charge jump at 0 -> replace breakpoint by its lower limit
    if problem == "discontinuous network flow":
        if interval[0] <= 0: # If 0 is part of the interval
            points_y[0] = 0 
    
    #2) Get lower convex envelope
    
    lower_x = []
    lower_y = []
    
    for i in range(0, len(points_x)):
        
        while (len(lower_x) >= 2 and 
               ((points_y[i] - lower_y[-1]) / (points_x[i] - lower_x[-1]) < 
                (lower_y[-1] - lower_y[-2]) / (lower_x[-1] - lower_x[-2]))):
            lower_x.pop()
            lower_y.pop()
            
        lower_x.append(points_x[i])
        lower_y.append(points_y[i])
    
    return lower_x, lower_y, time.time() - starttime


#---------------------------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------------------------------------------
# Function which implements the largest-error branching rule and calls a subfunction to solve the resulting child nodes by warm start
#------------------------------------------------------------------------------------------------------------------------------------

def branch_and_bound(m,n,plf_values,env_breakpoints_x,env_breakpoints_y,
                     PLFs_breakpoints_x,PLFs_breakpoints_y, problem):
    '''
    Implements the greatest error branching and calls a subfunction to solve the resulting child nodes by warm start.

    Parameters
    ----------
    m : GRB model
        model of parent node
    plf_values : list
        value of the plfs at the solution
    env_breakpoints_x : list
        convex envelope of the parent node. x-values
    env_breakpoints_y : list
        convex envelope of the parent node. y-values
    PLFs_breakpoints_x : list
        PLFs of the objective function. x values
    PLFs_breakpoints_y : list
        plfs of the objective function. y-values
    problem : string
        Test problem to be solved.

    Returns
    -------
    child_nodes : list
        contains all the solution information of the two child nodes
    comp_times : list
        contains the three computation times for grb, lp, and env

    '''
    #Return-objects: solved child nodes and computation times
    child_nodes = []
    comp_times = [0,0,0,0]
    
    #select variable with "worst" convex envelope
    approx_error = np.array(plf_values) - np.array(m.x[n:2*n])
    branch_index = np.argmax(approx_error)
    branch_point = m.x[branch_index]
    
    #Get posititon where to divide the old convex envelope
    pos = bisect.bisect(env_breakpoints_x[branch_index],branch_point)
    
    #Get primal dual basis of parent node for warm start
    var = m.VBasis #primal optimal basis of parent node
    con = m.CBasis #dual optimal basis of parent node
    
    #Solve right and left child nodes
    for side in ["left","right"]:

        #Solve the node by using function solve_node
        res,times = solve_node(side, m.copy(), var, con, branch_index, branch_point, n,
                              PLFs_breakpoints_x, PLFs_breakpoints_y,
                              pos, env_breakpoints_x.copy(), env_breakpoints_y.copy(), problem)
        #Increase time
        comp_times[0] += times[0]
        comp_times[1] += times[1]
        comp_times[2] += times[2]
        comp_times[3] += times[3]
        
        #Add child node information to list 
        child_nodes.append(res)
    
    return child_nodes, comp_times



#--------------------------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------
# Subfunction of branch-and-bound function to sole a child node by using warm start from parent node
#---------------------------------------------------------------------------------------------------

def solve_node(side, m, var, con, branch_index, branch_point,
             n, PLFs_breakpoints_x, PLFs_breakpoints_y,
             pos, env_breakpoints_x, env_breakpoints_y, problem):
    '''
    Subfunction of branch-and-bound function to sole a child node by using warm start from parent node.

    Parameters
    ----------
    side: string
        whether it is the left or the right child node
    m : Gurobi model
        Gurobi model of the parent node.
    var : list
        Primal optimal basis of parent node.
    con : list
        Dual optimal basis of parent node.
    n : integer
        Dimension of the problem.
    PLFs_breakpoints_x : list
        Contains the x-values of the PLFs breakpoints.
    PLFs_breakpoints_y : list
        Contains the y-values of the PLFs breakpoints.
    env_breakpoints_x : list
        Contains the x-values of the convex envelope's breakpoints of the parent node.
    env_breakpoints_y : list
        Contains the y-values of the convex envelope's breakpoints of the parent node.
    branch_index : integer
        The variables index on which the brnahcing takes place.
    branch_point : float
        Point at which the interval of the "branch_index" is divided.
    pos : integer
        Position on which the branch point divides the convex envelope of the 
        parent node.
    problem : string
        Test problem to be solved.
        
    Returns
    -------
    First 4 items: model, plf values, envelope
    Second three items: times of grb,lp and env operations.

    '''
    
    #Lower convex envelope of the parent node
    lower_x = env_breakpoints_x[branch_index]
    lower_y = env_breakpoints_y[branch_index]
    
    #Get narrowed interval over which the convex envlope can be tightened
    if side == "left":
        interval = [lower_x[pos-1], branch_point]
    else:
        interval = [branch_point, lower_x[pos]]
    
    #Calculate tighter convex envelope
    new_lower_x, new_lower_y, env_runtime = getEnvelope(
        PLFs_breakpoints_x[branch_index],
        PLFs_breakpoints_y[branch_index], interval, problem)
    
        
    # Counting time used by Gurobi interace
    grb_starttime = time.time()
    
    # Disable gurobi log output to console:: 0(no), 1(yes) 
    m.Params.LogToConsole = 0 
    
    #Add linear inequalitites to the model of the parent node
    x_i = m.getVars()[branch_index]       
    e_i = m.getVars()[n+branch_index]
    
    #Add branching constraint
    if side == "left":
        m.addConstr(x_i <= branch_point)
    else:
        m.addConstr(x_i >= branch_point)
    m.update()

    #Add constriants for tighter convex envelope
    for k in range(1,len(new_lower_x)):
        slope = (new_lower_y[k] - new_lower_y[k-1]) / (new_lower_x[k] - new_lower_x[k-1])
        m.addConstr(e_i >= slope * (x_i - new_lower_x[k-1]) + new_lower_y[k-1]) 
    m.update()
    
    #Warm Start
    x = m.getVars()
    c = m.getConstrs()
    for h in range(0, len(var)):
        x[h].VBasis = var[h]
    for h in range(0, len(con)):
        c[h].CBasis = con[h]
    for k in range(0,len(new_lower_x)-1): #basis for new constraints
        c[-(k+1)].CBasis = 0
    c[-(k+2)].CBasis = 0 
    m.update()
    
    #Time to construct model
    grb_runtime = time.time() - grb_starttime
    
    #Resolve model by dual simplex
    m.optimize()
    
    #Time to solve model
    lp_runtime = m.runtime
    
    #New convex envelope for all variables
    if side == "left":
        env_breakpoints_x[branch_index] = lower_x[0:pos-1] + new_lower_x
        env_breakpoints_y[branch_index] = lower_y[0:pos-1] + new_lower_y
    else:
        env_breakpoints_x[branch_index] = new_lower_x + lower_x[pos+1:]
        env_breakpoints_y[branch_index] = new_lower_y + lower_y[pos+1:]
        
    #Upper bound: PLF values
    plf_starttime = time.time()
    plf_values = [getPLFvalue(PLFs_breakpoints_x[i], PLFs_breakpoints_y[i], m.x[i], problem) 
               for i in range(0, n)]
    plf_runtime = time.time() - plf_starttime

    return (m,env_breakpoints_x, env_breakpoints_y, plf_values),(env_runtime, grb_runtime, lp_runtime, plf_runtime)

