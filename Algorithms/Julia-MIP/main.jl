#Reads the problem from csv file, solves it with MIP solver and writes solve times to csv

#Import
ENV["GUROBI_HOME"] = "C:\\ProgramData\\Microsoft\\Windows\\Start Menu\\Programs\\Gurobi 9.5.0 (win64)"
import Pkg
Pkg.add("Gurobi")
Pkg.build("Gurobi")
Pkg.add("JuMP")
Pkg.add("PiecewiseLinearOpt")
Pkg.add("CSV")
Pkg.add("DataFrames")
using JuMP, Gurobi, PiecewiseLinearOpt
using DataFrames, CSV, Statistics
include("MIP_solver.jl")

#------------
# Parameters 
#------------
solver_list = ["sBB", :Logarithmic, :DisaggLogarithmic, :ZigZag, :ZigZagInteger]
timelimit = 1800

#-----------------
# "Heat up" solver
#-----------------

model = Model(Gurobi.Optimizer)
@variable(model, x >= 0)
@variable(model, 0 <= y <= 3)
@objective(model, Min, 12x + 20y)
@constraint(model, c1, 6x + 8y >= 100)
@constraint(model, c2, 7x + 12y >= 120)
set_silent(model)
optimize!(model)

#--------------------
# Run computations
#--------------------

K_list = [10,100,500,1000,5000,10000]

for K in K_list

    #list of computation times
    times_list = [[] for method in solver_list]
    grb_time_list = []
    lp_time_list = []
    env_time_list = []
    plf_time_list = []

    #Read problem from csv
    csv_file = DataFrame(CSV.File("problem_info_K="*string(K)*".csv",header=0))
    row = 1
    while row <= nrow(csv_file)

        #Get problem information and sBB time
        problem = csv_file[row,1]
        n = trunc(Int,csv_file[row,2])
        K = trunc(Int,csv_file[row,3])

        #Get rhs
        if problem == "knapsack" || problem == "concave-knapsack"
            rhs = parse(Float64,csv_file[row+1,1])
        elseif problem == "network flow"
            rhs = [parse(Float64,csv_file[row+1,1])]
            for i in 2:trunc(Int,sqrt(n))
                push!(rhs, csv_file[row+1,i])
            end
        end

        #Get PLF x breakpoints
        PLFs_breakpoints_x=[]
        for i in 1:n
            breakpoints = [parse(Float64,csv_file[trunc(Int,row+i+1),1])]
            for j in 2:trunc(Int,K+1)
                push!(breakpoints,csv_file[trunc(Int,row+i+1),j])
            end
            push!(PLFs_breakpoints_x, breakpoints)
        end

        #Get PLF y breakpoints
        PLFs_breakpoints_y=[]
        for i in 1:n
            breakpoints = [parse(Float64,csv_file[trunc(Int,i+row+n+1),1])]
            for j in 2:trunc(Int,K+1)
                push!(breakpoints,csv_file[trunc(Int,i+row+n+1),j])
            end
            push!(PLFs_breakpoints_y, breakpoints)
        end

        #Solve problem
        push!(times_list[1],csv_file[row,4]) #sBB time
        for i in 2:length(solver_list)
            method = solver_list[i]
            sol_time = @elapsed obj_val, solution = solveMIP(PLFs_breakpoints_x, PLFs_breakpoints_y, rhs, n, problem, method, timelimit)
            push!(times_list[i],sol_time)
            println("Solve time: ", sol_time)
        end

        #Get sBB operation times
        push!(grb_time_list,csv_file[row,5]/csv_file[row,4])
        push!(lp_time_list,csv_file[row,6]/csv_file[row,4])
        push!(env_time_list,csv_file[row,7]/csv_file[row,4])
        push!(plf_time_list,csv_file[row,8]/csv_file[row,4])

        #Go to next instance
        row += 2+2*n
    end

    #Write to excel if instances for K are completed
    col_names = [string(method) for method in solver_list]
    df = DataFrame(times_list,col_names) 
    CSV.write("results/"*string(K)*".csv",df)
    df = DataFrame(grb=grb_time_list, lp=lp_time_list, env=env_time_list, plf=plf_time_list) 
    CSV.write("results/"*string(K)*"_sBB_operation.csv",df)

end