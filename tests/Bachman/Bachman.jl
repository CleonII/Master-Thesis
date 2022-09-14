using ModelingToolkit, DifferentialEquations, BenchmarkTools, DataFrames, CSV, LSODA, Sundials, ODEInterface, ODEInterfaceDiffEq
using ForwardDiff
using ReverseDiff
using Plots
using StatsBase
using Random
using LinearAlgebra
using Calculus
using Ipopt
using Optim
using Distributions
using Printf

include(joinpath(pwd(), "Additional_functions", "additional_tools.jl"))
include(joinpath(pwd(), "Additional_functions", "benchmarkSolvers.jl"))

# Ipopt wrapper 
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "CommonParameterEstimationMethods", "CreateIpoptProb.jl"))

# PeTab importer 
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "PeTabImporter.jl"))

# Observation function, needs to be generalised 
include(pwd() * "/Pipeline_ModelParameterEstimation/Data/model_Boehm_JProteomeRes2014/Boehm_JProteomeRes2014Obs.jl")
include(pwd() * "/Pipeline_ModelParameterEstimation/Data/model_Bachmann_MSB2011/Bachmann_MSB2011Obs.jl")

allModelFunctionVector = includeAllModels(getModelFiles(joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels")), 
        joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels"))
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "LatinHyperCubeParameters", "LatinHyperCubeSampledParameters.jl"))


solver, tol = Rodas4P(), 1e-9
modelName = "model_Bachmann_MSB2011"
evalObs = Bachmann_MSB2011
evalU0 = Bachmann_MSB2011_t0!
evalSd = Bachmann_MSB2011_sd!

paramVals = CSV.read(pwd() * "/tests/Bachman/Params.csv", DataFrame)
paramMat = paramVals[!, Not([:Id, :SOCS3RNAEqc, :CISRNAEqc])]


evalF, evalGradF, evalHessianApproxF, paramVecEstTmp, lowerBounds, upperBounds, idParam = setUpCostFunc(modelName, evalObs, evalU0, evalSd, solver, tol)
# "Exact" hessian via autodiff (should only be used for smaller models)
evalH = (hessianMat, paramVec) -> begin hessianMat .= Symmetric(ForwardDiff.hessian(evalF, paramVec)) end
a = 1
#cost = evalF(paramVecEstTmp)

p0 = collect(paramMat[1, :])
cost = evalF(p0)




experimentalConditions, measurementData, parameterBounds, observableData = readDataFiles(modelName, readObs=true)

# Read file with observed values and put into struct 
obsData = processObsData(measurementData, observableData)

i_relevant = findall(x -> x == "observable_pSTAT5B_rel", obsData.observebleDd)
print(obsData.yObs[i_relevant])


# Error is in paramMap, I always copy from the reference paramMap so changes are not captured accurately 
# when going between conditions. 


#conditionIdSol = ["model1_data2", "model1_data1", "model1_data11", "model1_data12", "model1_data13", "model1_data14", "model1_data4", "model1_data5", "model1_data7", "model1_data8", "model1_data38", "model1_data39", "model1_data40", "model1_data41", "model1_data42", "model1_data3", "model1_data10", "model1_data15", "model1_data16", "model1_data17", "model1_data18", "model1_data19", "model1_data20", "model1_data21", "model1_data22", "model1_data23", "model1_data24", "model1_data25", "model1_data6", "model1_data9", "model1_data26", "model1_data27", "model1_data30", "model1_data31", "model1_data32", "model1_data33"]