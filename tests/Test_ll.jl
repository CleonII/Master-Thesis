#= 
    Check the accruacy of the PeTab importer by checking the log-likelihood value against known values for several 
    models.
=#


using ModelingToolkit 
using DifferentialEquations
using DataFrames
using CSV 
using ForwardDiff
using ReverseDiff
using StatsBase
using Random
using LinearAlgebra
using Distributions
using Printf
using Calculus


# Relevant PeTab structs for compuations 
include(joinpath(pwd(), "src", "PeTab_structs.jl"))

# Functions for solving ODE system 
include(joinpath(pwd(), "src", "Solve_ODE_model", "Solve_ode_model.jl"))

# PeTab importer to get cost, grad etc 
include(joinpath(pwd(), "src", "PeTab_importer", "Create_cost_grad_hessian.jl"))

# HyperCube sampling 
include(joinpath(pwd(), "src", "Optimizers", "Lathin_hypercube.jl"))

# For converting to SBML 
include(joinpath(pwd(), "src", "SBML", "SBML_to_ModellingToolkit.jl"))


passTest = true


# Bachman model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Bachmann_MSB2011/"
peTabModel = setUpPeTabModel("model_Bachmann_MSB2011", dirModel)
peTabOpt = setUpCostGradHess(peTabModel, Rodas4P(), 1e-12)
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed)
diff = cost + (418.40573341425295)
if diff > 1e-3
    println("Does not pass ll-test for Bachman model")
    passTest = false
end


# Beer model - Currently we do not support parameter estimation for models with condition specific parameters to estimate
dirModel = pwd() * "/Intermediate/PeTab_models/model_Beer_MolBioSystems2014/"
peTabModel = setUpPeTabModel("model_Beer_MolBioSystems2014", dirModel)
peTabOpt = setUpCostGradHess(peTabModel, Rodas4P(), 1e-12)
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed)
diff = cost + (58622.9145631413)
if diff > 1e-3
    println("Does not pass ll-test for Beer model")
    passTest = false
end


# Boehm model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Boehm_JProteomeRes2014/"
peTabModel = setUpPeTabModel("model_Boehm_JProteomeRes2014", dirModel)
peTabOpt = setUpCostGradHess(peTabModel, Rodas4P(), 1e-12)
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed)
diff = cost + (-138.22199693517703)
if diff > 1e-3
    println("Does not pass ll-test for Boehm model")
    passTest = false
end


# Brännmark model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Brannmark_JBC2010/"
peTabModel = setUpPeTabModel("model_Brannmark_JBC2010", dirModel)
peTabOpt = setUpCostGradHess(peTabModel, Rodas4P(), 1e-12)
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed)
diff = abs(cost + (-141.889113770537))
if diff > 1e-3
    println("Does not pass ll-test for Brännmark model")
    passTest = false
end


# Bruno model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Bruno_JExpBot2016/"
peTabModel = setUpPeTabModel("model_Bruno_JExpBot2016", dirModel)
peTabOpt = setUpCostGradHess(peTabModel, Rodas4P(), 1e-12)
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed)
diff = abs(cost - (-46.688176988431806))
if diff > 1e-3
    println("Does not pass ll-test for Bruno model")
    passTest = false
end
# Check that we can compute gradients 
gradBrunoForward = ForwardDiff.gradient(peTabOpt.evalF, peTabOpt.paramVecTransformed) 
gradBrunoCalc = Calculus.gradient(peTabOpt.evalF, peTabOpt.paramVecTransformed) 
gradImporter = zeros(length(peTabOpt.paramVecTransformed))
peTabOpt.evalGradF(gradImporter, peTabOpt.paramVecTransformed)
if sum((gradImporter - gradBrunoCalc).^2) > 1e-6
    println("Does not pass gradient test for Bruno model")
    passTest = false
end


# Crauste model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Crauste_CellSystems2017/"
peTabModel = setUpPeTabModel("model_Crauste_CellSystems2017", dirModel)
peTabOpt = setUpCostGradHess(peTabModel, Rodas4P(), 1e-12)
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed)
diff = abs(cost + (-190.96521897435176))
if diff > 1e-2
    println("Does not pass ll-test for Crauste model")
    passTest = false
end


# Elowitz model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Elowitz_Nature2000/"
peTabModel = setUpPeTabModel("model_Elowitz_Nature2000", dirModel)
peTabOpt = setUpCostGradHess(peTabModel, Rodas4P(), 1e-12)
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed)
diff = abs(cost + (63.20279991419332))
if diff > 1e-3
    println("Does not pass ll-test for Elowitz model")
    passTest = false
end


# Fiedler model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Fiedler_BMC2016/"
peTabModel = setUpPeTabModel("model_Fiedler_BMC2016", dirModel)
peTabOpt = setUpCostGradHess(peTabModel, Rodas4P(), 1e-9)
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed)
diff = abs(cost + (58.58390161681))
if diff > 1e-3
    println("Does not pass ll-test for Fiedler model")
    passTest = false
end


# Fujita model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Fujita_SciSignal2010/"
peTabModel = setUpPeTabModel("model_Fujita_SciSignal2010", dirModel)
peTabOpt = setUpCostGradHess(peTabModel, Rodas4P(), 1e-12)
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed)
diff = abs(cost + (53.08377736998929))
if diff > 1e-3
    println("Does not pass ll-test for Fujita model")
    passTest = false
end


# Isensee model - Break code check why 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Isensee_JCB2018/"
peTabModel = setUpPeTabModel("model_Isensee_JCB2018", dirModel)
peTabOpt = setUpCostGradHess(peTabModel, Rodas4P(), 1e-9)
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed)
diff = abs(cost + (-3949.375966548649))
if diff > 1e-2
    println("Does not pass ll-test for Isensee model")
    passTest = false
end



# Lucarelli - Breaks code by the same reason as above 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Lucarelli_CellSystems2018/"
peTabModel = setUpPeTabModel("model_Lucarelli_CellSystems2018", dirModel)
peTabOpt = setUpCostGradHess(peTabModel, Rodas4P(), 1e-12)
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed)
diff = abs(cost + (-1681.6059879426584))
if diff > 1e-3
    println("Does not pass ll-test for Lucarelli model")
    passTest = false
end


# Schwen model - We do not handle Priors yet
#=
dirModel = pwd() * "/Intermediate/PeTab_models/model_Schwen_PONE2014/"
peTabModel = setUpPeTabModel("model_Schwen_PONE2014", dirModel)
peTabOpt = setUpCostGradHess(peTabModel, Rodas4P(), 1e-12)
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed)
=#


# Sneyd model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Sneyd_PNAS2002/"
peTabModel = setUpPeTabModel("model_Sneyd_PNAS2002", dirModel)
peTabOpt = setUpCostGradHess(peTabModel, Rodas4P(), 1e-12)
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed)
diff = abs(cost + (319.79177818768756))
if diff > 1e-3
    println("Does not pass ll-test for Sneyd model")
    passTest = false
end


# Weber model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Weber_BMC2015/"
peTabModel = setUpPeTabModel("model_Weber_BMC2015", dirModel)
peTabOpt = setUpCostGradHess(peTabModel, Rodas4P(), 1e-12)
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed)
diff = abs(cost + (-296.2017922646865))
if diff > 1e-3
    println("Does not pass ll-test for Weber model")
    passTest = false
end


# Zheng model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Zheng_PNAS2012/"
peTabModel = setUpPeTabModel("model_Zheng_PNAS2012", dirModel)
peTabOpt = setUpCostGradHess(peTabModel, Rodas4P(), 1e-12)
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed)
diff = abs(cost + (278.33353271001477))
if diff > 1e-3
    println("Does not pass ll-test for Zheng model")
    passTest = false
end


dirModel = pwd() * "/Intermediate/PeTab_models/model_Beer_MolBioSystems2014/"
peTabModel = setUpPeTabModel("model_Beer_MolBioSystems2014", dirModel)
a = 1

# Process PeTab files into type-stable Julia structs 
experimentalConditionsFile, measurementDataFile, parameterDataFile, observablesDataFile = readDataFiles(peTabModel.dirModel, readObs=true)
parameterData = processParameterData(parameterDataFile)
measurementData = processMeasurementData(measurementDataFile, observablesDataFile) 
simulationInfo = getSimulationInfo(measurementDataFile, measurementData)

# Indices for mapping parameter-estimation vector to dynamic, observable and sd parameters correctly when calculating cost
paramEstIndices = getIndicesParam(parameterData, measurementData, peTabModel.odeSystem, experimentalConditionsFile)

setParamToFileValues!(peTabModel.paramMap, peTabModel.stateMap, parameterData)

# The time-span 5e3 is overwritten when performing actual forward simulations 
odeProb = ODEProblem(peTabModel.odeSystem, peTabModel.stateMap, (0.0, 5e3), peTabModel.paramMap, jac=true, sparse=false)
odeProb = remake(odeProb, p = convert.(Float64, odeProb.p), u0 = convert.(Float64, odeProb.u0))

# Functions to map experimental conditions and parameters correctly to the ODE model 
changeToExperimentalCondUse! = (pVec, u0Vec, expID) -> changeExperimentalCond!(pVec, u0Vec, expID, parameterData, experimentalConditionsFile, peTabModel)
changeModelParamUse! = (pVec, u0Vec, paramEst, paramEstNames) -> changeModelParam!(pVec, u0Vec, paramEst, paramEstNames, paramEstIndices, peTabModel)

# Set up function which solves the ODE model for all conditions and stores result 
solveOdeModelAllCondUse! = (solArrayArg, odeProbArg) -> solveOdeModelAllExperimentalCond!(solArrayArg, odeProbArg, changeToExperimentalCondUse!, measurementDataFile, simulationInfo, solver, tol, tol, onlySaveAtTobs=true)

# The starting problem is to correctly identify dynamic parameters between conditions (to build the correct dynParamEst vector).
# Then the second problem to properly allow to ODE solver to realise it must the vector to estimate. I will not have to worry 
# about the transformation (it is handled anyway). 
# i) Fix function looking up for which parameters are dynamic (walk through the experimentalConditionsFile)
# ii) Write change experimentalCond function that can look into dynamic-param vector 
# iii) Integrate this new function into the optimisaiton 



