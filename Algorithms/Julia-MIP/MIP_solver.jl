#Solves the piecewise linear optimization problem casted as MIP by using a commercial solver.

function solveMIP(PLFs_breakpoints_x, PLFs_breakpoints_y, rhs, n, problem, method, timelimit)
    
    #Extracting lower and upper bounds for variables
    lower_bounds = []
    upper_bounds = []
    for breakpoints in PLFs_breakpoints_x
        push!(lower_bounds, breakpoints[1])
        push!(upper_bounds, breakpoints[end])
    end 

    # Define the model
    m = Model(Gurobi.Optimizer) 
    @variable(m, lower_bounds[i] <= x[i = 1:n] <= upper_bounds[i])
    if problem == "knapsack" || problem == "concave-knapsack"
        @constraint(m, sum(x) == rhs)
    elseif problem == "network flow"
        nr = trunc(Int,sqrt(n))
        @constraint(m, [i=1:nr], sum(x[j] for j in (i-1)*nr+1:i*nr) - sum(x[nr*j+i] for j in 0:nr-1) == rhs[i])
    end
    
    # Define the objective function
    z = []
    for i in 1:n
        push!(z,piecewiselinear(m, x[i], PLFs_breakpoints_x[i], PLFs_breakpoints_y[i],method=method))
    end
    @objective(m, Min, sum(z))
    
    # Run the solver
    set_silent(m)
    set_optimizer_attribute(m, "MIPGap", 0.00001)
    set_optimizer_attribute(m, "TimeLimit", timelimit)
    optimize!(m)

    println("Method: ", method)
    try
        println("lb: ", objective_bound(m), " ub: ", objective_value(m))
        return objective_value(m), value.(x)
    catch
        println("lb: ", objective_bound(m), " No upper bound found so far.")
        return 0, zeros(n)
    end
end