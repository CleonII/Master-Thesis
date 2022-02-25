using ModelingToolkit, DifferentialEquations, BenchmarkTools, DataFrames, CSV, LSODA, Sundials, ODEInterface, ODEInterfaceDiffEq

include(pwd() * "/Pipeline_ModelSolver/BigFloatODEProblem.jl")

function getModelFiles(path)
    if isdir(path)
        modelFiles = readdir(path)
        return modelFiles
    else
        println("No such directory")
        return nothing
    end
end

function getSolvers()
    algs_DiffEq = [Vern6(), Vern7(), Vern8(), Tsit5(), DP5(), DP8(), Feagin14(), VCABM(),
                   Rosenbrock23(), TRBDF2(), Rodas4(), Rodas4P(), Rodas4P2(), Rodas5(), QNDF(), FBDF(), 
                   Trapezoid(), KenCarp4(), Kvaerno5(), RadauIIA3(), RadauIIA5(), 
                   AutoTsit5(Rosenbrock23()), AutoVern7(Rodas5()), AutoVern9(Rodas4P()), AutoVern9(Rodas5())]
    algs_LSODA = [lsoda()]
    algs_Sundials = [CVODE_BDF(linear_solver=:Dense), CVODE_BDF(linear_solver=:LapackDense), CVODE_BDF(linear_solver=:GMRES), 
                     CVODE_Adams(linear_solver=:Dense), CVODE_Adams(linear_solver=:LapackDense), 
                     ARKODE(Sundials.Explicit(), order=4), ARKODE(Sundials.Explicit(), order=8), 
                     ARKODE(Sundials.Implicit(), order = 3), ARKODE(Sundials.Implicit(), order = 5)]
    algs_ODEInterface = [dopri5(), dop853(), radau(), radau5(), rodas()]
    algs = [algs_DiffEq; algs_LSODA; algs_Sundials; algs_ODEInterface]
    alg_hints = [[:auto], [:nonstiff], [:stiff]]
    return algs, alg_hints
end

function getHiAccSolver()
    algs = [AutoVern9(Rodas5()), Rodas4P()]
    return algs
end

function getTolerances()
    tolList = [1e-6, 1e-9, 1e-12, 1e-15] # 1e-6 is standard
    return [tolList[1]]
end

function fixDirectories(path)
    if ~isdir(path)
        println("Create directory: " * path)
        mkdir(path)
    end
end

function getNumberOfFiles(path)
    if isdir(path)
        files = readdir(path)
        return length(files)
    else
        mkdir(path)
        return 0
    end
end

function modelSolver(modelFile, timeEnd, solvers, hiAccSolvers, Tols, iterations, readPath, writefile)
    nonStiffHiAccSolver, stiffHiAccSolver = hiAccSolvers
    println(modelFile)
    modelPath = readPath * "/" * modelFile
    include(modelPath)

    new_sys = ode_order_lowering(sys)
    u0 = initialSpeciesValues

    p = trueParameterValues

    c = trueConstantsValues

    tspan = (0.0, timeEnd)

    # Estimate ground truth
    println("Obtaining high accuracy solution")
    bfProb = BigFloatODEProblem(new_sys, u0, tspan, [p;c])
    local hiAccSol
    try 
        hiAccSol = solve(bfProb, stiffHiAccSolver, reltol = 1e-15, abstol = 1e-15) 
    catch 
        hiAccSol = solve(bfProb, nonStiffHiAccSolver, reltol = 1e-15, abstol = 1e-15) 
    end
    ts = hiAccSol.t
    println("Done")
    
    # Solve with different solvers
    prob = ODEProblem(new_sys,u0,tspan,[p;c],jac=true)

    alg_solvers, alg_hints = solvers
    for alg_solver in alg_solvers
        println("Alg_solver = ", alg_solver)
        if ~((modelFile == "model_Crauste_CellSystems2017.jl") && alg_solver == AutoTsit5(Rosenbrock23())) 
            println(alg_solver)
            for tol in Tols
            
                benchRunTime = Vector{Float64}(undef, iterations)
                benchMemory = Vector{Float64}(undef, iterations)
                benchAllocs = Vector{Float64}(undef, iterations)
                sqDiff = Float64
                success = true

                try
                    sol = solve(prob, alg_solver, reltol = tol, abstol = tol, saveat = ts)
                    if sol.t[end] == timeEnd && sol.retcode == :Success
                        sqDiff = sum((sol[:,:] - hiAccSol[:,:]).^2)
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
        end
    end
    for alg_hint in alg_hints
        if ~(modelFile == "model_Crauste_CellSystems2017.jl") && ~(modelFile == "model_Chen_MSB2009.jl")
            println(alg_hint)
            for tol in Tols
                benchRunTime = Vector{Float64}(undef, iterations)
                benchMemory = Vector{Float64}(undef, iterations)
                benchAllocs = Vector{Float64}(undef, iterations)
                sqDiff = Float64
                success = true

                try
                    sol = solve(prob, alg_hint = alg_hint, reltol = tol, abstol = tol, saveat = ts)
                    if sol.t[end] == timeEnd && sol.retcode == :Success
                        sqDiff = sum((sol[:,:] - hiAccSol[:,:]).^2)
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
        end
    end
    
end


function modelSolverIterator(modelFiles, timeEnds, solvers, hiAccSolvers, Tols, iterations, readPath, writefile)
    for modelFile in modelFiles
        timeEnd = timeEnds[timeEnds[:,1] .== modelFile, 2][1]
        modelSolver(modelFile, timeEnd, solvers, hiAccSolvers, Tols, iterations, readPath, writefile)
    end
end


function getBigFloatProb(modelFile)
    
    writePath = pwd() * "/Pipeline_ModelSolver/IntermediaryResults"
    readPath = pwd() * "/Pipeline_SBMLImporter/JuliaModels"
    fixDirectories(writePath)
    timeEnds = CSV.read(writePath * "/timeScales.csv", DataFrame)
    timeEnds = CSV.read(writePath * "/timeScales.csv", DataFrame)
    timeEnd = timeEnds[timeEnds[:,1] .== modelFile, 2][1]

    modelPath = readPath * "/" * modelFile
    include(modelPath)
    new_sys = ode_order_lowering(sys)
    u0 = initialSpeciesValues
    p = trueParameterValues
    c = trueConstantsValues
    tspan = (0.0, timeEnd)

    # Estimate ground truth
    bfProb = BigFloatODEProblem(new_sys, u0, tspan, [p;c])

    return bfProb
end


function main(;modelFiles=["all"], modelsExclude=[""])
    readPath = pwd() * "/Pipeline_SBMLImporter/JuliaModels"
    writePath = pwd() * "/Pipeline_ModelSolver/IntermediaryResults"
    fixDirectories(writePath)
    writefile = writePath * "/benchmark_" * string(getNumberOfFiles(writePath) + 1) * ".csv"

    if modelFiles[1] == "all"
        modelFiles = getModelFiles(readPath)
    end
    for modelExclude in modelsExclude
        modelFiles = modelFiles[.!(occursin.(modelFiles, modelExclude))]
    end
    
    timeEnds = CSV.read(writePath * "/timeScales.csv", DataFrame)
    solvers = getSolvers()
    hiAccSolvers = getHiAccSolver()
    tolList = getTolerances()
    iterations = 15
    
    modelSolverIterator(modelFiles, timeEnds, solvers, hiAccSolvers, tolList, iterations, readPath, writefile)
end


modelFiles = main(modelFiles=["model_Elowitz_Nature2000.jl"], modelsExclude=["model_Chen_MSB2009.jl", "model_Rahman_MBS2016.jl", "model_SalazarCavazos_MBoC2020.jl"])
