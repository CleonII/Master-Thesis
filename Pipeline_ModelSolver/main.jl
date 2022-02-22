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
    relTols = [1e-3, 1e-6, 1e-9, 1e-12, 1e-15] # 1e-3 is standard
    absTols = [1e-6, 1e-9, 1e-12, 1e-15] # 1e-6 is standard
    return [relTols[1]], [absTols[1]]
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

function modelSolver(modelFile, timeEnd, solvers, hiAccSolvers, relTols, absTols, iterations, readPath, writefile)
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
    minRelTol = minimum(relTols)
    minAbsTol = minimum(absTols)
    bfProb = BigFloatODEProblem(new_sys, u0, tspan, [p;c])
    local hiAccSol
    try 
        hiAccSol = solve(bfProb, nonStiffHiAccSolver, relTol = minRelTol, absTol = minAbsTol) 
    catch 
        hiAccSol = solve(bfProb, stiffHiAccSolver, relTol = minRelTol, absTol = minAbsTol) 
    end
    ts = hiAccSol.t
    
    # Solve with different solvers
    prob = ODEProblem(new_sys,u0,tspan,[p;c],jac=true)

    alg_solvers, alg_hints = solvers
    for alg_solver in alg_solvers
        if ~((modelFile == "model_Crauste_CellSystems2017.jl") && alg_solver == AutoTsit5(Rosenbrock23())) && ~(modelFile == "model_Chen_MSB2009.jl")
            println(alg_solver)
            for relTol in relTols
                for absTol in absTols
                    benchRunTime = Vector{Float64}(undef, iterations)
                    benchMemory = Vector{Float64}(undef, iterations)
                    benchAllocs = Vector{Float64}(undef, iterations)
                    sqDiff = Float64
                    success = true

                    try
                        sol = solve(prob, alg_solver, relTol = relTol, absTol = absTol, saveat = ts)
                        if sol.t[end] == timeEnd
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
                            b = @benchmark solve($prob, $alg_solver, relTol = $relTol, absTol = $absTol) samples=1 evals=1
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

                    data = DataFrame(model = modelFile, solver = alg_solver, relTol = relTol, absTol = absTol, 
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
    for alg_hint in alg_hints
        if ~(modelFile == "model_Crauste_CellSystems2017.jl") && ~(modelFile == "model_Chen_MSB2009.jl")
            println(alg_hint)
            for relTol in relTols
                for absTol in absTols
                    benchRunTime = Vector{Float64}(undef, iterations)
                    benchMemory = Vector{Float64}(undef, iterations)
                    benchAllocs = Vector{Float64}(undef, iterations)
                    sqDiff = Float64
                    success = true

                    try
                        sol = solve(prob, alg_hint = alg_hint, relTol = relTol, absTol = absTol, saveat = ts)
                        if sol.t[end] == timeEnd
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
                            b = @benchmark solve($prob, alg_hint = $alg_hint, relTol = $relTol, absTol = $absTol) samples=1 evals=1
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

                    data = DataFrame(model = modelFile, solver = alg_hint[1], relTol = relTol, absTol = absTol, 
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
    
end

function modelSolverIterator(modelFiles, timeEnds, solvers, hiAccSolvers, relTols, absTols, iterations, readPath, writefile)
    for modelFile in modelFiles
        timeEnd = timeEnds[timeEnds[:,1] .== modelFile, 2][1]
        modelSolver(modelFile, timeEnd, solvers, hiAccSolvers, relTols, absTols, iterations, readPath, writefile)
    end
end


function main()
    readPath = pwd() * "/Pipeline_SBMLImporter/JuliaModels"
    writePath = pwd() * "/Pipeline_ModelSolver/IntermediaryResults"
    fixDirectories(writePath)
    writefile = writePath * "/benchmark_" * string(getNumberOfFiles(writePath) + 1) * ".csv"

    modelFiles = getModelFiles(readPath)
    timeEnds = CSV.read(writePath * "/timeScales.csv", DataFrame)
    solvers = getSolvers()
    hiAccSolvers = getHiAccSolver()
    relTols, absTols = getTolerances()
    iterations = 15
    
    modelSolverIterator(modelFiles, timeEnds, solvers, hiAccSolvers, relTols, absTols, iterations, readPath, writefile)

end

main()