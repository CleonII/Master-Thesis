using ModelingToolkit, DifferentialEquations, DataFrames, CSV, Random
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


allModelFunctionVector = includeAllModels(getModelFiles(joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels")), 
        joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels"))



function benchMethod_proto(usedModelFunction, iStartPar, filesAndPaths, timeEnd, experimentalConditions, measurementData, observables, parameterBounds,
        optAlg, method)
    write = joinpath(filesAndPaths.writePath, filesAndPaths.writeFile)
    
    
    if method == "adjointSensitivities" || method == "all"
        println("Running: adjointSensitivities \n")
        senseAlgs = getSenseAlgs()

        for senseAlg in senseAlgs
            model, p, doLogSearch = adjointSensitivities(usedModelFunction, iStartPar, senseAlg, optAlg, 
                    timeEnd, experimentalConditions, measurementData, observables, parameterBounds)
            optModelSaveResults(model, p, doLogSearch, senseAlg, optAlg, "-", iStartPar, "adjointSensitivities", write)
        end
    end
    if method == "forwardAutomaticDifferentiation" || method == "all"
        println("Running: forwardAutomaticDifferentiation \n")
        chunkSizes = [23, 60, 115]

        for chunkSize in chunkSizes
            model, p, doLogSearch = forwardAutomaticDifferentiation(usedModelFunction, iStartPar, chunkSize, optAlg, 
                    timeEnd, experimentalConditions, measurementData, observables, parameterBounds)
            optModelSaveResults(model, p, doLogSearch, chunkSize, optAlg, "-", iStartPar, "forwardAutomaticDifferentiation", write)
        end
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

                    step, adam_opt, doLogSearch = forwardGradient(usedModelFunction, iStartPar, n_it, b2, stepRange, 
                            timeEnd, experimentalConditions, measurementData, observables, parameterBounds)
                    optAdamSaveResults(step, adam_opt, doLogSearch, n_it, b2, stepRange, iStartPar, iterations, "forwardGradient", write)

                end
            end
        end
    end

    nothing

end





function iterateMethods(benchMethod)
    
    for iStartPar=1:50
        benchMethod(iStartPar)
    end

    nothing

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
    usedModelFunction = allModelFunctionVector[[allModelFile in [modelFile] for allModelFile in allModelFiles]]

    benchMethod = (iStartPar) -> benchMethod_proto(usedModelFunction, iStartPar, filesAndPaths, timeEnd, experimentalConditions, 
            measurementData, observables, parameterBounds, optAlg, method)

    println("Starting benchmark")
    
    iterateMethods(benchMethod)

    nothing

end

main(modelName = "model_Raimundez_PCB2020", method = "forwardAutomaticDifferentiation")

#optAlgs = [:LD_MMA, :LD_LBFGS, :Ipopt]
#methods = ["adjointSensitivities", "forwardAutomaticDifferentiation", "forwardGradient", "all"]