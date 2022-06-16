using ModelingToolkit, DifferentialEquations, DataFrames, CSV, Random
using JuMP, NLopt, Ipopt, LinearAlgebra, DiffEqSensitivity, ForwardDiff
using ModelingToolkit: varmap_to_vars
using ForwardDiff: GradientConfig, Chunk
using LatinHypercubeSampling
using DelimitedFiles
using Optim 

# Make sure to load all correct libraries and files etc
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


# Generate the cost function, gradient function and hessian function for a model (only currently supports Boehm). 
# The cost and gradient (etc) will then be calculated by the provided solver (below I use Rodas4P) as it often 
# is a good solver. Function also return number of parameters 
function genarateCostGradHess(modelName, solver)

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

    fUse = (x) -> fPre(x...)
    fGradUse = (grad, x) -> fGradPre(grad, x...)
    fHessUse = (hess, x) -> fHessPre(hess, x...)
 
    return fUse, fGradUse, fHessUse, nParam
end

solver = Rodas4P()
modelName = "model_Boehm_JProteomeRes2014"
fUse, fGradUse, fHessUse, nParam = genarateCostGradHess(modelName, solver)

# Reading a parameter vector from the stored hypercube
fileSaveCube = pwd() * "/Intermediate/Bohem/Cube.csv"
cube = Matrix(CSV.read(fileSaveCube, DataFrame))
p = cube[1, :] # A random parameter vector 

cost = fUse(p)
# Note gradients and hessian are modified in place 
gradVec = zeros(nParam)
hessMat = zeros(nParam, nParam)
fGradUse(gradVec, p)
fHessUse(hessMat, p)

println("Cost = ", cost)
println("Grad = ", gradVec)
println("Hess = ", hessMat)