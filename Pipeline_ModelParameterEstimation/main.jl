using ModelingToolkit, DifferentialEquations, BenchmarkTools, DataFrames, CSV
using JuMP, NLopt, Ipopt, LinearAlgebra, DiffEqSensitivity, ForwardDiff
using ModelingToolkit: varmap_to_vars
using ForwardDiff: GradientConfig, Chunk
using LatinHypercubeSampling

println("Done loading modules")

include(joinpath(pwd(), "Additional_functions", "additional_tools.jl"))
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "LatinHyperCubeParameters", "LatinHyperCubeSampledParameters.jl"))
include(joinpath(pwd(), "Additional_functions", "importModelInfo.jl"))

include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "optAndSave.jl"))
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "model_Bachmann_MSB2011", "adjointSensitivities.jl"))
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "model_Bachmann_MSB2011", "forwardAutomaticDifferentiation.jl"))
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "model_Bachmann_MSB2011", "forwardGradient.jl"))
include(joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels", "model_Bachmann_MSB2011.jl"))

function orderParameters(modelParameters, parameterBounds;
    scaleDeterminer = "scale", offsetDeterminer = "offset", varianceDeterminer = "sd_")

    numScale = modelParameters.numScale
    numOffset = modelParameters.numOffset
    numVariance = modelParameters.numVariance
    numOptParameters = modelParameters.numOptParameters

    scaleNames = modelParameters.scaleNames
    offsetNames = modelParameters.offsetNames
    varianceNames = modelParameters.varianceNames
    optParameterNames = modelParameters.optParameterNames

    newScaleVector = Vector{Float64}(undef, numScale)
    newOffsetVector = Vector{Float64}(undef, numOffset)
    newVarianceVector = Vector{Float64}(undef, numVariance)
    newDynamicParameterVector = Vector{Float64}(undef, numOptParameters)

    for (i, parId) in enumerate(parameterBounds[!, 1])
        if occursin(scaleDeterminer, parId)
            sIndex = findfirst(lowercase(parId) .== lowercase.(scaleNames))
            newScaleVector[sIndex] = parameterBounds[i, :nominalValue]
        elseif occursin(offsetDeterminer, parId)
            oIndex = findfirst(lowercase(parId) .== lowercase.(offsetNames))
            newOffsetVector[oIndex] = parameterBounds[i, :nominalValue]
        elseif occursin(varianceDeterminer, parId)
            vIndex = findfirst(lowercase(parId) .== lowercase.(varianceNames))
            newVarianceVector[vIndex] = parameterBounds[i, :nominalValue]^2
        else
            pIndex = findfirst(lowercase(parId) .== lowercase.(optParameterNames))
            newDynamicParameterVector[pIndex] = parameterBounds[i, :nominalValue]
        end
    end

    return vcat(newScaleVector, newOffsetVector, newVarianceVector, newDynamicParameterVector)
end

function calcCorrectCost(timeEnd, experimentalConditions, measurementData, observables, parameterBounds)

    sys, initialSpeciesValues, trueParameterValues = getODEModel()
    new_sys = ode_order_lowering(sys)
    u0 = initialSpeciesValues
    pars = trueParameterValues 
    tspan = (0.0, timeEnd)
    prob = ODEProblem(new_sys,u0,tspan,pars,jac=true)


    initVariableNames = ["EpoRJAK2_CIS", "SHP1", "STAT5", "EpoRJAK2", "SOCS3", "CIS"]

    observableVariableNames = ["CISRNA", "CIS", "SHP1", "SHP1Act", "SOCS3RNA", "SOCS3", "STAT5", 
            "p12EpoRpJAK2", "p1EpoRpJAK2", "p2EpoRpJAK2", "EpoRpJAK2", "pSTAT5"]

    parameterInU0Names = ["init_EpoRJAK2_CIS", "init_SHP1", "init_SHP1_multiplier", "SHP1ProOE", "init_STAT5", 
            "init_EpoRJAK2", "init_SOCS3_multiplier", "SOCS3EqcOE", "SOCS3Eqc", "init_CIS_multiplier", "CISEqc", "CISEqcOE"]

    parameterInObservableNames = ["CISRNAEqc", "CISEqc", "SOCS3RNAEqc", "SOCS3Eqc", "init_EpoRJAK2", "init_STAT5", "init_SHP1"]

    # Initialize structs

    modelData = ModelData(new_sys, prob, observables, experimentalConditions, 
            initVariableNames, observableVariableNames, parameterInU0Names, parameterInObservableNames)

    experimentalData = ExperimentalData(observables, experimentalConditions, measurementData, modelData)

    modelParameters = ModelParameters(new_sys, prob, parameterBounds, experimentalConditions, measurementData, observables, experimentalData)

    parameterSpace, numAllStartParameters, lowerBounds, upperBounds = ParameterSpace(modelParameters, parameterBounds)

    modelOutput = ModelOutput(Float64, experimentalData, modelParameters)

    # Get correct parameters

    corrParameters = orderParameters(modelParameters, parameterBounds)

    #

    updateAllParameterVectors = () -> updateAllParameterVectors_proto(modelParameters, modelData)

    solveODESystem = (iCond) -> solveODESystem_AdjSens_proto(prob, modelParameters.dynamicParametersVector, modelParameters.u0Vector, modelData, modelOutput, iCond)

    g_unscaledObservables = (type, u, dynPar, i, iCond) -> g_unscaledObservables_AdjSens_proto(type, u, dynPar, i, iCond, modelData, experimentalData)

    g_scaledObservationFunctions = (type, h_bar, scaleVector, offsetVector, iCond) -> g_scaledObservationFunctions_AdjSens_proto(type, h_bar, scaleVector, offsetVector, iCond, modelParameters, experimentalData, modelData)

    g_cost = (type, h_hat, varianceVector, i, iCond) -> g_cost_AdjSens_proto(type, h_hat, varianceVector, i, iCond, modelParameters, experimentalData)

    g = (u, dynPar, scale, offset, variance, i, iCond; type = get_type([u, dynPar, scale, offset, variance])) -> g_AdjSens_proto(u, dynPar, scale, offset, variance, i, iCond, g_unscaledObservables, g_scaledObservationFunctions, g_cost; type = get_type([u, dynPar, scale, offset, variance]))

    G = (iCond) -> G_AdjSens_proto(solveODESystem, g, iCond, modelParameters, experimentalData, modelOutput)

    allConditionsCost = (p...) -> allConditionsCost_AdjSens_proto(parameterSpace, modelParameters, experimentalData, modelData, 
    updateAllParameterVectors, G, p...)

    view(corrParameters, parameterSpace.doLogSearch) .= log10.(view(corrParameters, parameterSpace.doLogSearch))
    cost = allConditionsCost(corrParameters...)

    return cost
    
end





function benchMethod(methodFile, iStartPar, iterations, 
        filesAndPaths, timeEnd, experimentalConditions, measurementData, observables, parameterBounds)
    write = joinpath(filesAndPaths.writePath, filesAndPaths.writeFile)
    
    # runtime, objective_value (at each step), end parameters, iterations 
    senseAlgs = getSenseAlgs()
    optAlgs = [:LD_MMA, :LD_LBFGS, :Ipopt]
    chunkSizes = [23, 60, 115]
    stepRanges = [[0.0, -3.0], [-1.0, -4.0], [-2.0, -5.0], [-3.0, -3.0]]
    b2s = [0.999, 0.9]
    n_its = [100, 500, 1000, 5000, 10000]
    # r?
    # c?

    if methodFile == "adjointSensitivities.jl"
        println("Running: adjointSensitivities.jl \n")

        for senseAlg in senseAlgs
            model, p, doLogSearch = adjointSensitivities(iStartPar, senseAlg, optAlgs[1], 
                    timeEnd, experimentalConditions, measurementData, observables, parameterBounds)
            optModelSaveResults(model, p, doLogSearch, senseAlg, optAlgs[1], "-", iStartPar, iterations, methodFile, write)
        end
        for optAlg in optAlgs
            model, p, doLogSearch = adjointSensitivities(iStartPar, senseAlgs[1], optAlg, 
                    timeEnd, experimentalConditions, measurementData, observables, parameterBounds)
            optModelSaveResults(model, p, doLogSearch, senseAlgs[1], optAlg, "-", iStartPar, iterations, methodFile, write)
        end
        
    elseif methodFile == "forwardAutomaticDifferentiation.jl"
        println("Running: forwardAutomaticDifferentiation.jl \n")

        for chunkSize in chunkSizes
            model, p, doLogSearch = forwardAutomaticDifferentiation(iStartPar, chunkSize, optAlgs[1], 
                    timeEnd, experimentalConditions, measurementData, observables, parameterBounds)
            optModelSaveResults(model, p, doLogSearch, chunkSize, optAlgs[1], "-", iStartPar, iterations, methodFile, write)
        end
        for optAlg in optAlgs
            model, p, doLogSearch = forwardAutomaticDifferentiation(iStartPar, chunkSizes[1], optAlg, 
                    timeEnd, experimentalConditions, measurementData, observables, parameterBounds)
            optModelSaveResults(model, p, doLogSearch, chunkSizes[1], optAlg, "-", iStartPar, iterations, methodFile, write)
        end
        
    elseif methodFile == "forwardGradient.jl"
        println("Running: forwardGradient.jl \n")

        for n_it in n_its
            step, adam_opt, doLogSearch = forwardGradient(iStartPar, n_it, b2s[1], stepRanges[2], 
                    timeEnd, experimentalConditions, measurementData, observables, parameterBounds)
            optAdamSaveResults(step, adam_opt, doLogSearch, n_it, b2s[1], stepRanges[2], iStartPar, iterations, methodFile, write)
        end

        for b2 in b2s
            step, adam_opt, doLogSearch = forwardGradient(iStartPar, n_its[3], b2, stepRanges[2], 
                    timeEnd, experimentalConditions, measurementData, observables, parameterBounds)
            optAdamSaveResults(step, adam_opt, doLogSearch, n_its[3], b2, stepRanges[2], iStartPar, iterations, methodFile, write)
        end

        for stepRange in stepRanges
            step, adam_opt, doLogSearch = forwardGradient(iStartPar, n_its[3], b2s[1], stepRange, 
                    timeEnd, experimentalConditions, measurementData, observables, parameterBounds)
            optAdamSaveResults(step, adam_opt, doLogSearch, n_its[3], b2s[1], stepRange, iStartPar, iterations, methodFile, write)
        end

    end

    nothing

end





function iterateMethods(methodFiles, iterations, 
        filesAndPaths, timeEnd, experimentalConditions, measurementData, observables, parameterBounds)
    
    for iStartPar=1:50
        for methodFile in methodFiles
            benchMethod(methodFile, iStartPar, iterations, 
                    filesAndPaths, timeEnd, experimentalConditions, measurementData, observables, parameterBounds)
        end
    end

    nothing

end



function main(; modelName = "model_Bachmann_MSB2011")
    
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
    #

    #correctCost = calcCorrectCost(filesAndPaths, timeEnd, experimentalConditions, measurementData, observables, parameterBounds)
    #println("Cost with correct parameters: ", correctCost) # 

    methodFiles = readdir(methodPath)

    println("Starting benchmark")
    
    iterations = 5
    iterateMethods(methodFiles, iterations, 
            filesAndPaths, timeEnd, experimentalConditions, measurementData, observables, parameterBounds)

    nothing

end

main()