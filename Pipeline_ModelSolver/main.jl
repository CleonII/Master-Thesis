using ModelingToolkit, DifferentialEquations, BenchmarkTools, DataFrames, CSV, LSODA, Sundials, ODEInterface, ODEInterfaceDiffEq


include(joinpath(pwd(), "Pipeline_ModelSolver", "BigFloatODEProblem.jl"))
include(joinpath(pwd(), "Additional_functions", "additional_tools.jl"))
include(joinpath(pwd(), "Additional_functions", "Solver_info.jl"))
include(joinpath(pwd(), "Additional_functions", "benchmarkSolvers.jl"))

allModelFunctionVector = includeAllModels(getModelFiles(joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels")), 
        joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels"))


function modelSolver(modelFunction, modelFile, solvers, hiAccSolvers, Tols, iterations, writefile)
    
    nonStiffHiAccSolver, stiffHiAccSolver = hiAccSolvers
    println(modelFile)

    sys, initialSpeciesValues, trueParameterValues = modelFunction()
    new_sys = ode_order_lowering(sys)
    # Read experimental conditions and set up simulation conditions 
    modelName = replace.(modelFile, ".jl" => "")
    experimentalConditions, measurementData, parameterBounds = readDataFiles(modelName)
    stateMap = initialSpeciesValues
    paramMap = trueParameterValues     
    # Set stateMap and paramMap non condition parameter to the nominal reporeted values n parameter-files
    setParamToParamFileVal!(paramMap, stateMap, parameterBounds) 
    # Get data on experimental conditions 
    firstExpIds, shiftExpIds, simulateSS, parameterNames, stateNames = getSimulationInfo(measurementData, sys)
    # Set up for first experimental condtition 
    changeToCondUse! = (pVec, u0Vec, expID) -> changeToCond!(pVec, u0Vec, expID, parameterBounds, experimentalConditions, parameterNames, stateNames, paramMap, stateMap)
    prob = ODEProblem(new_sys, stateMap, (0.0, 5e3), paramMap, jac=true)
    prob = remake(prob, p = convert.(Float64, prob.p), u0 = convert.(Float64, prob.u0)) # Ensure everything is of proper data type (Float64)

    # Estimate ground truth
    println("Obtaining high accuracy solution")
    tspan = (0.0, 5.0e3)
    bfProb = BigFloatODEProblem(new_sys, stateMap, tspan, paramMap)
    local hiAccSolArr
    local successHighAcc
    try 
        hiAccSolArr, successHighAcc = solveOdeModelAllCond(bfProb, changeToCondUse!, simulateSS, measurementData, firstExpIds, shiftExpIds, 1e-15, nonStiffHiAccSolver, nTSave=100)
    catch 
        hiAccSolArr, successHighAcc = solveOdeModelAllCond(bfProb, changeToCondUse!, simulateSS, measurementData, firstExpIds, shiftExpIds, 1e-15, stiffHiAccSolver, nTSave=100)
    end

    if successHighAcc != true
        println("High accuracy solver failed for $modelFile")
        open(joinpath(pwd(), "Pipeline_ModelSolver", "Log.txt"), "a") do io
            println(io, "Failed with high accuracy solution for $modelFile")
        end
        return 
    else
        println("Done with high accuracy solver")
    end

    alg_solvers, alg_hints = solvers
    alg_solvers = [Rosenbrock23()]
    for alg_solver in alg_solvers
        println("Alg_solver = ", alg_solver)
        if ~((modelFile == "model_Crauste_CellSystems2017.jl") && alg_solver == AutoTsit5(Rosenbrock23())) # Crashes 
            for tol in Tols
            
                benchRunTime = Vector{Float64}(undef, iterations)
                benchMemory = Vector{Float64}(undef, iterations)
                benchAllocs = Vector{Float64}(undef, iterations)
                local sqDiff = Float64

                try
                    sqDiff = calcSqErr(prob, changeToCondUse!, hiAccSolArr, simulateSS, measurementData, firstExpIds, shiftExpIds, tol, alg_solver)
                catch 
                    sqDiff = Inf
                end
                success = isinf(sqDiff) ? false : true
                    
                if success
                    for i in 1:iterations
                        b = @benchmark solveOdeModelAllCond($prob, $changeToCondUse!, $simulateSS, $measurementData, $firstExpIds, $shiftExpIds, $tol, $alg_solver) samples=1 evals=1
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
        GC.gc()
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
                    sqDiff = calcSqErr(prob, changeToCondUse!, hiAccSolArr, simulateSS, measurementData, firstExpIds, shiftExpIds, tol, alg_hint)
                catch 
                    sqDiff = Inf
                end
                success = isinf(sqDiff) ? false : true

                if success
                    for i in 1:iterations
                        b = @benchmark solveOdeModelAllCond($prob, $changeToCondUse!, $simulateSS, $measurementData, $firstExpIds, $shiftExpIds, $tol, $alg_hint) samples=1 evals=1
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

    GC.gc()
end


function modelSolverIterator(modelFunctionVector, modelFiles, solvers, hiAccSolvers, Tols, iterations, writefile)
    for (i, modelFile) in enumerate(modelFiles)
        modelFunction = modelFunctionVector[i]
        return modelSolver(modelFunction, modelFile, solvers, hiAccSolvers, Tols, iterations, writefile)
    end
end


function main(; modelFiles=["all"], modelsExclude=[""], writefile::String="")
    readPath = joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels")
    writePath = joinpath(pwd(), "Pipeline_ModelSolver", "IntermediaryResults")

    if writefile == ""
        fixDirectories(writePath)
        writefile = joinpath(writePath, "benchmark_" * string(getNumberOfFiles(writePath) + 1) * ".csv")
    end

    allModelFiles = getModelFiles(readPath)

    if modelFiles[1] == "all"
        modelFiles = getModelFiles(readPath)
    end
    for modelExclude in modelsExclude
        modelFiles = modelFiles[.!(occursin.(modelFiles, modelExclude))]
    end

    usedModelFunctionVector = allModelFunctionVector[[allModelFile in modelFiles for allModelFile in allModelFiles]]
        
    solvers = getSolvers()
    hiAccSolvers = getHiAccSolver()
    tolList = getTolerances(onlyMaxTol = false) # Get tolerances = [1e-6, 1e-9, 1e-12]
    iterations = 3
    
    return modelSolverIterator(usedModelFunctionVector, modelFiles, solvers, hiAccSolvers, tolList, iterations, writefile)
end

dirSave = pwd() * "/Intermediate/ODESolvers/"
if !isdir(dirSave)
    mkpath(dirSave)
end
fileSave = dirSave * "Benchmark1.csv"

modelListUse = ["model_Beer_MolBioSystems2014.jl", "model_Weber_BMC2015.jl", "model_Schwen_PONE2014.jl", "model_Alkan_SciSignal2018.jl", 
    "model_Bachmann_MSB2011.jl", "model_Bertozzi_PNAS2020.jl", "model_Blasi_CellSystems2016.jl", "model_Boehm_JProteomeRes2014.jl", 
    "model_Borghans_BiophysChem1997.jl", "model_Brannmark_JBC2010.jl", "model_Bruno_JExpBot2016.jl", "model_Crauste_CellSystems2017.jl", 
    "model_Elowitz_Nature2000.jl", "model_Fiedler_BMC2016.jl", "model_Fujita_SciSignal2010.jl", "model_Giordano_Nature2020.jl", 
    "model_Isensee_JCB2018.jl", "model_Laske_PLOSComputBiol2019.jl", "model_Lucarelli_CellSystems2018.jl", "model_Okuonghae_ChaosSolitonsFractals2020.jl", 
    "model_Oliveira_NatCommun2021.jl", "model_Perelson_Science1996.jl", "model_Rahman_MBS2016.jl", "model_Raimundez_PCB2020.jl", 
    "model_SalazarCavazos_MBoC2020.jl", "model_Sneyd_PNAS2002.jl", "model_Zhao_QuantBiol2020.jl", "model_Zheng_PNAS2012.jl"]

modelListUse = ["model_Alkan_SciSignal2018.jl"]

for i in eachindex(modelListUse)
    println("Starting with a new model")
    main(modelFiles=[modelListUse[i]], modelsExclude=["model_Chen_MSB2009.jl"], writefile=fileSave)
    GC.gc()
end
