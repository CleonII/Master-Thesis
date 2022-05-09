using ModelingToolkit, DifferentialEquations, BenchmarkTools, DataFrames, CSV, LSODA, Sundials, ODEInterface, ODEInterfaceDiffEq


include(joinpath(pwd(), "Pipeline_ModelSolver", "BigFloatODEProblem.jl"))
include(joinpath(pwd(), "Additional_functions", "additional_tools.jl"))
include(joinpath(pwd(), "Additional_functions", "Solver_info.jl"))

allModelFunctionVector = includeAllModels(getModelFiles(joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels")), 
        joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels"))


function modelSolver(modelFunction, modelFile, timeEnd, solvers, hiAccSolvers, Tols, iterations, writefile)
    nonStiffHiAccSolver, stiffHiAccSolver = hiAccSolvers
    println(modelFile)

    sys, initialSpeciesValues, trueParameterValues = modelFunction()
    new_sys = ode_order_lowering(sys)
    u0 = initialSpeciesValues
    pars = trueParameterValues 
    tspan = (0.0, timeEnd)

    # Estimate ground truth
    println("Obtaining high accuracy solution")
    bfProb = BigFloatODEProblem(new_sys, u0, tspan, pars)
    local hiAccSol
    try 
        hiAccSol = solve(bfProb, stiffHiAccSolver, reltol=1e-15, abstol=1e-15) 
    catch 
        hiAccSol = solve(bfProb, nonStiffHiAccSolver, reltol=1e-15, abstol=1e-15) 
    end

    if hiAccSol.retcode != :Success
        println("High accuracy solver failed for $modelFile")
        open(joinpath(pwd(), "Pipeline_ModelSolver", "Log.txt"), "a") do io
            println(io, "Failed with high accuracy solution for $modelFile")
        end
        return 
    else
        println("Done with high accuracy solver")
    end

    ts = hiAccSol.t
    
    # Solve with different solvers
    prob = ODEProblem(new_sys,u0,tspan,pars,jac=true)

    alg_solvers, alg_hints = solvers
    for alg_solver in alg_solvers
        println("Alg_solver = ", alg_solver)
        if ~((modelFile == "model_Crauste_CellSystems2017.jl") && alg_solver == AutoTsit5(Rosenbrock23())) 
            for tol in Tols
            
                benchRunTime = Vector{Float64}(undef, iterations)
                benchMemory = Vector{Float64}(undef, iterations)
                benchAllocs = Vector{Float64}(undef, iterations)
                sqDiff = Float64
                success = true

                try
                    sol = solve(prob, alg_solver, reltol = tol, abstol = tol, saveat = ts)
                    if sol.t[end] == timeEnd && sol.retcode == :Success
                        if length(sol.t) != length(ts)
                            index = findall(sol.t .!= setdiff(sol.t, ts))
                            sqDiff = sum((sol[:,index] - hiAccSol[:,:]).^2)
                        else
                            sqDiff = sum((sol[:,:] - hiAccSol[:,:]).^2)
                        end
                    else
                        sqDiff = NaN
                        success = false
                    end
                catch 
                    sqDiff = NaN
                    success = false
                end
                    
                if success
                    for i in 1:iterations
                        b = @benchmark solve($prob, $alg_solver, reltol = $tol, abstol = $tol) samples=1 evals=1
                        bMin = minimum(b)
                        benchRunTime[i] = bMin.time # microsecond
                        benchMemory[i] = bMin.memory # bytes
                        benchAllocs[i] = bMin.allocs # number of allocations
                    end
                else
                    benchRunTime .= NaN
                    benchMemory .= NaN
                    benchAllocs .= NaN
                end
                data = DataFrame(model = modelFile, solver = alg_solver, reltol = tol, abstol = tol, 
                                 success = success, runTime = benchRunTime, memory = benchMemory, allocs = benchAllocs,
                                 sqDiff = sqDiff, iteration = 1:iterations)
                if isfile(writefile)
                    CSV.write(writefile, data, append = true)
                else
                    CSV.write(writefile, data)
                end
            end
        else
            for tol in Tols
                data = DataFrame(model = modelFile, solver = alg_solver, reltol = tol, abstol = tol, 
                                    success = false, runTime = NaN, memory = NaN, allocs = NaN,
                                    sqDiff = NaN, iteration = 1:iterations)
                if isfile(writefile)
                    CSV.write(writefile, data, append = true)
                else
                    CSV.write(writefile, data)
                end
            end
        end
    end
    for alg_hint in alg_hints
        println("Alg_hint = ", alg_hint[1])
        if ~(modelFile == "model_Crauste_CellSystems2017.jl")
            for tol in Tols
                benchRunTime = Vector{Float64}(undef, iterations)
                benchMemory = Vector{Float64}(undef, iterations)
                benchAllocs = Vector{Float64}(undef, iterations)
                sqDiff = Float64
                success = true

                try
                    sol = solve(prob, alg_hint = alg_hint, reltol = tol, abstol = tol, saveat = ts)
                    if sol.t[end] == timeEnd && sol.retcode == :Success
                        if length(sol.t) != length(ts)
                            index = findall(sol.t .!= setdiff(sol.t, ts))
                            sqDiff = sum((sol[:,index] - hiAccSol[:,:]).^2)
                        else
                            sqDiff = sum((sol[:,:] - hiAccSol[:,:]).^2)
                        end
                    else
                        sqDiff = NaN
                        success = false
                    end
                catch 
                    sqDiff = NaN
                    success = false
                end

                if success
                    for i in 1:iterations
                        b = @benchmark solve($prob, alg_hint = $alg_hint, reltol = $tol, abstol = $tol) samples=1 evals=1
                        bMin = minimum(b)
                        benchRunTime[i] = bMin.time # microsecond
                        benchMemory[i] = bMin.memory # bytes
                        benchAllocs[i] = bMin.allocs # number of allocations
                    end
                else
                    benchRunTime .= NaN
                    benchMemory .= NaN
                    benchAllocs .= NaN
                end

                data = DataFrame(model = modelFile, solver = alg_hint[1], reltol = tol, abstol = tol, 
                                 success = success, runTime = benchRunTime, memory = benchMemory, allocs = benchAllocs, 
                                 sqDiff = sqDiff, iteration = 1:iterations)
                if isfile(writefile)
                    CSV.write(writefile, data, append = true)
                else
                    CSV.write(writefile, data)
                end
            end
        else
            for tol in Tols
                data = DataFrame(model = modelFile, solver = alg_hint[1], reltol = tol, abstol = tol, 
                                    success = false, runTime = NaN, memory = NaN, allocs = NaN, 
                                    sqDiff = NaN, iteration = 1:iterations)
                if isfile(writefile)
                    CSV.write(writefile, data, append = true)
                else
                    CSV.write(writefile, data)
                end
            end
        end
    end
    
end


function modelSolverIterator(modelFunctionVector, modelFiles, timeEnds, solvers, hiAccSolvers, Tols, iterations, writefile)
    for (i, modelFile) in enumerate(modelFiles)
        timeEnd = timeEnds[timeEnds[:,1] .== modelFile, 2][1]
        modelFunction = modelFunctionVector[i]
        modelSolver(modelFunction, modelFile, timeEnd, solvers, hiAccSolvers, Tols, iterations, writefile)
    end
end


function main(; modelFiles=["all"], modelsExclude=[""])
    readPath = joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels")
    writePath = joinpath(pwd(), "Pipeline_ModelSolver", "IntermediaryResults")
    fixDirectories(writePath)
    writefile = joinpath(writePath, "benchmark_" * string(getNumberOfFiles(writePath) + 1) * ".csv")

    allModelFiles = getModelFiles(readPath)

    if modelFiles[1] == "all"
        modelFiles = getModelFiles(readPath)
    end
    for modelExclude in modelsExclude
        modelFiles = modelFiles[.!(occursin.(modelFiles, modelExclude))]
    end

    usedModelFunctionVector = allModelFunctionVector[[allModelFile in modelFiles for allModelFile in allModelFiles]]
    
    timeEnds = CSV.read(joinpath(writePath, "timeScales.csv"), DataFrame)
    solvers = getSolvers()
    hiAccSolvers = getHiAccSolver()
    tolList = getTolerances(onlyMaxTol = false) # Get tolerances = [1e-6, 1e-9, 1e-12]
    iterations = 15
    
    modelSolverIterator(usedModelFunctionVector, modelFiles, timeEnds, solvers, hiAccSolvers, tolList, iterations, writefile)
end


main(modelFiles=["model_Laske_PLOSComputBiol2019.jl"], modelsExclude=["model_Chen_MSB2009.jl"])
