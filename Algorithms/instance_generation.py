"""
Contains the functions to generate test problems and run computations. It includes:
    - getRHS (generates the right hand side values of the problems)
    - getconcavePLFs (generates concave PLFs for the network flow problem)
    - getPLFs (generates PLFs by approximating nonconvex functions)
"""

#import statements
import numpy as np
import time 
import bisect
import random

#-----------------------------------------------------------------
# Function to generate the right hand side values for the problems
#-----------------------------------------------------------------
def getRHS(n,K,problem,PLFs_breakpoints_x):
    
    #Extracting lower and upper bounds for variables
    lower_bounds = []
    upper_bounds = []
    for breakpoints in PLFs_breakpoints_x:
        lower_bounds.append(breakpoints[0])
        upper_bounds.append(breakpoints[-1])
    
    #Compute right hand side value for knapsack problem
    if problem == "knapsack" or problem == "concave-knapsack":
        l = sum(lower_bounds)
        r = sum(upper_bounds)
        length = r-l
        rhs = random.uniform(l+0.25*length,r-0.25*length)
    elif problem == "network flow":
        rhs = []
        for i in range(0,int(np.sqrt(n))-1):
            rv1 = random.randint(1, 3)
            if rv1 == 1: #transshipment node
                rhs.append(0)
            elif rv1 == 2: #supply node
                rhs.append(random.uniform(5,50))
            else: #demand node
                rhs.append(-random.uniform(5,50))
        #Balance network
        rhs.append(-sum(rhs))  
    else:
        print("Problem not recogonised. Run again with either 'knapsack' or 'network flow'")  
           
    return rhs

#--------------------------------------------------------------
# Funtion to generate concave PLFs for the network flow problem
#--------------------------------------------------------------

def getNetworkPLFs(n,K):
    
    PLFs_breakpoints_x = []
    PLFs_breakpoints_y = []
    
    for i in range(0,n):
        
        breakpoints_x = [0]
        breakpoints_y = [0]
        
        breakpoints_x.append(random.uniform(5,50))
        for k in range(1,K):
            new = random.uniform(breakpoints_x[0],breakpoints_x[-1])
            while new in breakpoints_x:
                new = random.uniform(breakpoints_x[0],breakpoints_x[-1])
            bisect.insort(breakpoints_x, random.uniform(breakpoints_x[0],breakpoints_x[-1]))
        
        #Generate and sort slopes
        slopes = []
        for k in range(0,K):
            slopes.append(random.uniform(1,2000)/1000)
        slopes.sort(reverse=True)
        for k in range(0,K):
            breakpoints_y.append( slopes[k] * (breakpoints_x[k+1] - breakpoints_x[k])
                                 + breakpoints_y[-1])        
    
        #Append breakpoints to PLF list
        PLFs_breakpoints_x.append(breakpoints_x)
        PLFs_breakpoints_y.append(breakpoints_y)
        
    return PLFs_breakpoints_x, PLFs_breakpoints_y

#-------------------------------------------------------------------
# Function to generate n nonconvex PLFs out of randomly selected functions
#-----------------------------------------------------------------------

def getPLFs(n, K, func_number, concave):
    ''' 
    Generates n PLFs with K segments according to the selected functions 
    listed in the array func_number.
    
    Input:
    n: number of functions (int)
    K: number of segments of each PLF (int)
    func_number: array of integers which define which function should be slected
    concave: If True sgements are sorted so that concave function arise
    
    Output:
    PLFs_breakpoints_x: List which contains n lists where ach contain the breakpoints 
    of a PLF projected on the x-axis
    PLFs_breakpoints_y: List which contains n lists where ach contain the breakpoints 
    of a PLF projected on the y-axis
    
    '''
    
    PLFs_breakpoints_x = []
    PLFs_breakpoints_y = []
    
    for num in func_number:
        
        #Select function
        interval = globals()["X"+str(num)]
        function = globals()["func"+str(num)]
        
        #Generate x breakpoints
        breakpoints_x = [random.uniform(interval[0],interval[1]) for k in range(K-1)]
        breakpoints_x.append(interval[0])
        breakpoints_x.append(interval[1])
        breakpoints_x.sort()
        #Generate y breakpoints
        breakpoints_y = []
        for point in breakpoints_x:
            breakpoints_y.append(function(point))
            
        #Make function concave if wanted
        if concave == True:
            slopes = []
            for k in range(K):
                slope = (breakpoints_y[k+1] - breakpoints_y[k]) / (breakpoints_x[k+1] -breakpoints_x[k])
                slopes.append(slope)
            slopes.sort(reverse=True)   
            for k in range(K):
                breakpoints_y[k+1] = slopes[k]*(breakpoints_x[k+1]-breakpoints_x[k]) + breakpoints_y[k]
        
        #Append breakpoints to PLF list
        PLFs_breakpoints_x.append(breakpoints_x)
        PLFs_breakpoints_y.append(breakpoints_y)
        

    return PLFs_breakpoints_x, PLFs_breakpoints_y

#----------------------------------------------------------
# Set of univariate continuous functions (def func) over their support (X)
#----------------------------------------------------------

#Function 1
X1 = np.array([-5,5])
def func1(x): return np.exp(-3*x-12) - x**2 + 20

#Function 2
X2 = np.array([-5,5])
def func2(x): return -0.2*np.exp(-x) + x**2

#Function 3
X3 = np.array([-5,5])
def func3(x): return x**3 *np.exp(-x**2)

#Function 4
X4 = np.array([-10,10])
def func4(x): return (x**5 - 20*x**2 +5) / (x**4 +1)

#Function 5
X5 = np.array([0.1,10])
def func5(x): return np.log(3*x) * np.log(2*x) - 1

#Function 6
X6 = np.array([0.1,10])
def func6(x): return 10*np.log(x) - 3*x + (x-5)**2

#Function 7
X7 = np.array([-10,10])
def func7(x): return (-x**5 - 10*x**2)/(x**6 + 5)

#Function 8
X8 = np.array([-5,5])
def func8(x): return  x * np.exp(-x**2)

#Function 9
X9 = np.array([-4,4])
def func9(x): return -x**7 / 5040 + x**5 / 120 - x**3 / 3 + x

#Function 10
X10 = np.array([-10,10])
def func10(x): return (x**2 - 5*x +6) / (x**2 +1) -1

#Function 11
X11 = np.array([-1,7])
def func11(x): return x**4 - 12*x**3 + 47*x**2 - 60*x

#Function 12
X12 = np.array([-4,4])
def func12(x): return x**6 -15*x**4 + 27*x**2 +250

#Function 13
X13 = np.array([0,5])
def func13(x): return x**4 - 10*x**3 +35*x**2 -50*x + 24

#Function 14
X14 = np.array([-1,4])
def func14(x): return 0.2 *x**5 -1.25*x**4 +2.33*x**3 -2.5*x**2 +6*x

#Function 15
X15 = np.array([-4,4])
def func15(x): return x**3 -7*x + 7

#Function 16
X16 = np.array([-5,5])
def func16(x): return (x**4 - 4*x +10) / (x**2 +1) -1

#Function 17
X17 = np.array([-10,10])
def func17(x): return -x**5 *np.exp(-x**2) 

#Function 18
X18 = np.array([-1.5,3])
def func18(x): return x**5 -3*x**4 +4*x**3 + 2*x**2 -10*x -4

#Function 19
X19 = np.array([-5,5])
def func19(x): return (x**3 - 5*x +6) / (x**2 +1) -1

#Function 20
X20 = np.array([0.1,10])
def func20(x): return 1/x + 2*np.log(x) -2
