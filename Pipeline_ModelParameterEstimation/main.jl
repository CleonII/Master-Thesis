using ModelingToolkit, DifferentialEquations, DataFrames, CSV, Random
using Sundials
using JuMP, NLopt, Ipopt, LinearAlgebra, DiffEqSensitivity, ForwardDiff
using ModelingToolkit: varmap_to_vars
using ForwardDiff: GradientConfig, Chunk
using LatinHypercubeSampling
using DelimitedFiles
using Optim 


println("Done loading modules")

include(joinpath(pwd(), "Additional_functions", "additional_tools.jl"))
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "LatinHyperCubeParameters", "LatinHyperCubeSampledParameters.jl"))
include(joinpath(pwd(), "Additional_functions", "importModelInfo.jl"))
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "optAndSave.jl"))

include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "CommonParameterEstimationMethods", "adjointSensitivities.jl"))
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "CommonParameterEstimationMethods", "forwardAutomaticDifferentiation.jl"))
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "CommonParameterEstimationMethods", "forwardGradient.jl"))
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "CommonParameterEstimationMethods", "CreateIpoptProb.jl"))

allModelFunctionVector = includeAllModels(getModelFiles(joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels")), 
        joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels"))
modelsWithWorkingMethods = ["model_Alkan_SciSignal2018", "model_Bachmann_MSB2011", "model_Boehm_JProteomeRes2014"]
allMethodFunctionArray = includeAllMethods(joinpath.(pwd(), "Pipeline_ModelParameterEstimation", modelsWithWorkingMethods))

# Include with hessian (currently testing)
include(pwd() * "/Pipeline_ModelParameterEstimation/model_Boehm_JProteomeRes2014/forwardAutomaticDifferentiation_with_hessian.jl")

function benchMethod_proto(usedModelFunction, usedMethodFunctions, iStartPar, filesAndPaths, timeEnd, experimentalConditions, measurementData, observables, parameterBounds,
        optAlg, method)
    write = joinpath(filesAndPaths.writePath, filesAndPaths.writeFile)
    
    solver = QNDF()
    
    if method == "adjointSensitivities" || method == "all"
        println("Running: adjointSensitivities \n")
        senseAlgs = getSenseAlgs()
        solver_adjoint = KenCarp4()

        for senseAlg in senseAlgs
            model, p, doLogSearch = usedMethodFunctions[1](usedModelFunction, iStartPar, senseAlg, optAlg, solver, solver_adjoint, 
                    timeEnd, experimentalConditions, measurementData, observables, parameterBounds)
            optModelSaveResults(model, p, doLogSearch, solver, optAlg, senseAlg, solver_adjoint, iStartPar, "adjointSensitivities", write)
        end
    end
    if method == "forwardAutomaticDifferentiation" || method == "all"
        println("Running: forwardAutomaticDifferentiation \n")

        model, p, doLogSearch = usedMethodFunctions[2](usedModelFunction, iStartPar, optAlg, solver, 
                timeEnd, experimentalConditions, measurementData, observables, parameterBounds)
        optModelSaveResults(model, p, doLogSearch, solver, optAlg, "-", "-", iStartPar, "forwardAutomaticDifferentiation", write)

    end
    if method == "forwardGradient" || method == "all"
        println("Running: forwardGradient \n")
        stepRanges = [[0.0, -3.0], [-1.0, -4.0], [-2.0, -5.0], [-3.0, -3.0]]
        b2s = [0.999, 0.9]
        n_its = [100, 500, 1000, 5000, 10000]
        iterations = 4

        for n_it in n_its
            for b2 in b2s
                for stepRange in stepRanges

                    step, adam_opt, doLogSearch = usedMethodFunctions[3](usedModelFunction, iStartPar, n_it, b2, stepRange, solver, 
                            timeEnd, experimentalConditions, measurementData, observables, parameterBounds)
                    optAdamSaveResults(step, adam_opt, doLogSearch, solver, n_it, b2, stepRange, iStartPar, iterations, "forwardGradient", write)

                end
            end
        end
    end

    nothing

end


function iterateMethods(benchMethod, numStartGuesses)
    
    for iStartPar=1:numStartGuesses
        benchMethod(iStartPar)
    end

    nothing

end


function createOptimProb(fPre::Function, fGradPre::Function, fHessPre::Function, lowerBounds, upperBounds)

    nParam = length(lowerBounds)
    p0 = zeros(nParam)
    fUse = (x) -> fPre(x...)
    fGradUse = (grad, x) -> fGradPre(grad, x...)
    fHessUse = (hess, x) -> fHessPre(hess, x...)
    df = TwiceDifferentiable(fUse, fGradUse, fHessUse, p0)
    dfc = TwiceDifferentiableConstraints(lowerBounds .- 0.01, upperBounds .+ 0.01)
    
    return df, dfc
end


function main(; modelName = "model_Bachmann_MSB2011", optAlg = :Ipopt, method = "all")
    
    modelPath = joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels")
    modelNameShort = modelName[7:end]
    modelFile = modelName * ".jl" 
    methodPath = joinpath(pwd(), "Pipeline_ModelParameterEstimation", modelName)
    writePath = joinpath(pwd(), "Pipeline_ModelParameterEstimation", "IntermediaryResults")
    writeFile = joinpath("benchmark_" * modelName * "_" * string(getNumberOfFiles(writePath) + 1) * ".csv")

    filesAndPaths = FilesAndPaths(modelName, modelPath, modelFile, methodPath, writePath, writeFile)

    fixDirectories(writePath)

    # Reading data for model
    readDataPath = joinpath(pwd(), "Pipeline_ModelParameterEstimation", "Data", modelName)
    dataEnding = modelNameShort * ".tsv"
    experimentalConditions = CSV.read(joinpath(readDataPath, "experimentalCondition_" * dataEnding), DataFrame)
    measurementData = CSV.read(joinpath(readDataPath, "measurementData_" * dataEnding), DataFrame)
    timeEnd = maximum(measurementData[!, 5])
    observables = CSV.read(joinpath(readDataPath, "observables_" * dataEnding), DataFrame)
    parameterBounds = CSV.read(joinpath(readDataPath, "parameters_" * dataEnding), DataFrame)
    
    allModelFiles = getModelFiles(modelPath)
    usedModelFunction = allModelFunctionVector[[allModelFile in [modelFile] for allModelFile in allModelFiles]][1]

    usedMethodFunctions = allMethodFunctionArray[findfirst(modelsWithWorkingMethods .== modelName), :]

    benchMethod = (iStartPar) -> benchMethod_proto(usedModelFunction, usedMethodFunctions, iStartPar, filesAndPaths, timeEnd, experimentalConditions, 
            measurementData, observables, parameterBounds, optAlg, method)

    println("Starting benchmark")
    
    numStartGuesses = 10
    iterateMethods(benchMethod, numStartGuesses)

    nothing

end


function testOptimizerHessian(modelName, dirSave, solver, nEvals)


    # Set up directory where to save results 
    if !isdir(dirSave)
        mkpath(dirSave)
    end

    # Set up paths for data files 
    modelPath = joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels")
    modelNameShort = modelName[7:end]
    modelFile = modelName * ".jl" 
    methodPath = joinpath(pwd(), "Pipeline_ModelParameterEstimation", modelName)
    # Reading data for model
    readDataPath = joinpath(pwd(), "Pipeline_ModelParameterEstimation", "Data", modelName)
    dataEnding = modelNameShort * ".tsv"
    experimentalConditions = CSV.read(joinpath(readDataPath, "experimentalCondition_" * dataEnding), DataFrame)
    measurementData = CSV.read(joinpath(readDataPath, "measurementData_" * dataEnding), DataFrame)
    observables = CSV.read(joinpath(readDataPath, "observables_" * dataEnding), DataFrame)
    parameterBounds = CSV.read(joinpath(readDataPath, "parameters_" * dataEnding), DataFrame)
    allModelFiles = getModelFiles(modelPath)
    usedModelFunction = allModelFunctionVector[[allModelFile in [modelFile] for allModelFile in allModelFiles]][1]

    # Generate a cost function, inplace gradient function, in place hessian and lower and upper bounds 
    fPre, fGradPre, fHessPre, lowerBounds, upperBounds = forwardAutomaticDifferentiation_hessian_proto_model_Boehm_JProteomeRes2014(usedModelFunction, solver, experimentalConditions, measurementData, observables, parameterBounds)
    nParam = length(lowerBounds)

    # Set up functions into a usable format for hupercube and get start-guesses
    fCube = (x) -> fPre(x...)
    fileSaveCube = dirSave * "Cube.csv"
    createCube(1000, lowerBounds, upperBounds, fileSaveCube, fCube)
    # Read cube 
    cube = Matrix(CSV.read(fileSaveCube, DataFrame))

    # Ipopt and Optim problem 
    ipoptProb, iterArr = createIpoptProb(fPre, fGradPre, fHessPre, lowerBounds, upperBounds)
    Ipopt.AddIpoptIntOption(ipoptProb, "print_level", 0)
    Ipopt.AddIpoptIntOption(ipoptProb, "max_iter", 1000)
    dfOpt, dfcOpt = createOptimProb(fPre, fGradPre, fHessPre, lowerBounds, upperBounds)

    # Evaluate gradient + hessian (allow precompilation)
    dfOpt.f(cube[1, :])
    dfOpt.df(zeros(nParam), cube[1, :])
    dfOpt.h(zeros(nParam, nParam), cube[1, :])

    # Where to save results 
    # Optim 
    costOptim = zeros(Float64, nEvals)
    nIterOptim = zeros(Int, nEvals)
    retCodeOptim = Array{Any, 1}(undef, nEvals)
    runTimeOptim = Array{Any, 1}(undef, nEvals)
    # Ipopt 
    costIpopt = zeros(Float64, nEvals)
    nIterIpopt = zeros(Int, nEvals)
    retCodeIpopt = Array{Any, 1}(undef, nEvals)
    runTimeIpopt = Array{Any, 1}(undef, nEvals)
    # Ipopt no hessian 
    costIpoptBFGS = zeros(Float64, nEvals)
    nIterIpoptBFGS = zeros(Int, nEvals)
    retCodeIpoptBFGS = Array{Any, 1}(undef, nEvals)
    runTimeIpoptBFGS = Array{Any, 1}(undef, nEvals)

    fileSave = dirSave * "OptRes.csv"
    for i in 1:nEvals

        println("I = $i of $nEvals")

        p0 = cube[i, :] # Sample from hypercube 
        # Ipopt with hessian 
        Ipopt.AddIpoptStrOption(ipoptProb, "hessian_approximation", "exact")
        ipoptProb.x = deepcopy(p0)
        benchRunTime = @elapsed  sol_opt = Ipopt.IpoptSolve(ipoptProb)
        # Save relevant values 
        costIpopt[i] = ipoptProb.obj_val
        nIterIpopt[i] = iterArr[1]
        retCodeIpopt[i] = ipoptProb.status
        runTimeIpopt[i] = benchRunTime

        # Ipopt without hessian 
        Ipopt.AddIpoptStrOption(ipoptProb, "hessian_approximation", "limited-memory")
        ipoptProb.x = deepcopy(p0)
        benchRunTime = @elapsed  sol_opt = Ipopt.IpoptSolve(ipoptProb)
        # Save relevant values 
        costIpoptBFGS[i] = ipoptProb.obj_val
        nIterIpoptBFGS[i] = iterArr[1]
        retCodeIpoptBFGS[i] = ipoptProb.status
        runTimeIpoptBFGS[i] = benchRunTime

        # Optim with Hessian 
        res = Optim.optimize(dfOpt, dfcOpt, p0, IPNewton(), Optim.Options(iterations = 1000, show_trace = false))
        costOptim[i] = res.minimum
        nIterOptim[i] = res.iterations
        retCodeOptim[i] = res.f_converged
        runTimeOptim[i] = res.time_run

        # Save after each iteration (do not loose data)
        dataSave = ["IpoptHess" costIpopt[i] runTimeIpopt[i] retCodeIpopt[i] nIterIpopt[i] i;
                "IpoptBGFS" costIpoptBFGS[i] runTimeIpoptBFGS[i] retCodeIpoptBFGS[i] nIterIpoptBFGS[i] i;
                "OptimInt" costOptim[i] runTimeOptim[i] retCodeOptim[i] nIterOptim[i] i]
        dataSave = DataFrame(dataSave, ["Alg", "Cost", "Run_time", "Ret_code", "N_iter", "Start_guess"])
        shouldAppend = isfile(fileSave) ? true : false
        CSV.write(fileSave, dataSave, append=shouldAppend)
    end
end


function genIpoptHess(modelName, solver)

    # Set up paths for data files 
    modelPath = joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels")
    modelNameShort = modelName[7:end]
    modelFile = modelName * ".jl" 
    methodPath = joinpath(pwd(), "Pipeline_ModelParameterEstimation", modelName)
    # Reading data for model
    readDataPath = joinpath(pwd(), "Pipeline_ModelParameterEstimation", "Data", modelName)
    dataEnding = modelNameShort * ".tsv"
    experimentalConditions = CSV.read(joinpath(readDataPath, "experimentalCondition_" * dataEnding), DataFrame)
    measurementData = CSV.read(joinpath(readDataPath, "measurementData_" * dataEnding), DataFrame)
    observables = CSV.read(joinpath(readDataPath, "observables_" * dataEnding), DataFrame)
    parameterBounds = CSV.read(joinpath(readDataPath, "parameters_" * dataEnding), DataFrame)
    allModelFiles = getModelFiles(modelPath)
    usedModelFunction = allModelFunctionVector[[allModelFile in [modelFile] for allModelFile in allModelFiles]][1]

    # Generate a cost function, inplace gradient function, in place hessian and lower and upper bounds 
    fPre, fGradPre, fHessPre, lowerBounds, upperBounds = forwardAutomaticDifferentiation_hessian_proto_model_Boehm_JProteomeRes2014(usedModelFunction, solver, experimentalConditions, measurementData, observables, parameterBounds)
    nParam = length(lowerBounds)

    ipoptProb, iterArr = createIpoptProb(fPre, fGradPre, fHessPre, lowerBounds, upperBounds)
    return ipoptProb, iterArr
end


dirSave = pwd() * "/Intermediate/Bohem/"
modelName = "model_Boehm_JProteomeRes2014"
solver = QNDF()
testOptimizerHessian(modelName, dirSave, solver, 2)
