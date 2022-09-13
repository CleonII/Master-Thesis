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


"""
    testCostGradHess(solver, tol; printRes::Bool=false)

    Compare cost and gradient for Julia PeTab importer against the cost and 
    gradient computed in PyPesto for Boehm model. 
"""
function compareAgainstPyPesto(solver, tol; printRes::Bool=false)

    modelName = "model_Boehm_JProteomeRes2014"
    evalObs = Boehm_JProteomeRes2014
    evalU0 = Boehm_JProteomeRes2014_t0!
    evalSd = Boehm_JProteomeRes2014_sd!

    evalF, evalGradF, evalHessianApproxF, paramVecEstTmp, lowerBounds, upperBounds, idParam = setUpCostFunc(modelName, evalObs, evalU0, evalSd, solver, tol)
    # "Exact" hessian via autodiff (should only be used for smaller models)
    evalH = (hessianMat, paramVec) -> begin hessianMat .= Symmetric(ForwardDiff.hessian(evalF, paramVec)) end

    paramVals = CSV.read(pwd() * "/tests/Boehm/Params.csv", DataFrame)
    paramMat = paramVals[!, Not([:Id, :ratio, :specC17])]
    iUse = [findfirst(x -> x == idParam[i], names(paramMat)) for i in eachindex(idParam)]
    nParam = ncol(paramMat)

    costPython = (CSV.read(pwd() * "/tests/Boehm/Cost.csv", DataFrame))[!, :Cost]
    gradPythonMat = CSV.read(pwd() * "/tests/Boehm/Grad.csv", DataFrame)
    gradPythonMat = gradPythonMat[!, Not([:Id, :ratio, :specC17])]

    for i in 1:nrow(paramMat)
        paramVec = collect(paramMat[i, :])

        costJulia = evalF(paramVec)
        sqDiffCost = (costJulia - costPython[i])^2
        if sqDiffCost > 1e-5
            @printf("sqDiffCost = %.3e\n", sqDiffCost)
            @printf("Does not pass test on cost\n")
            return false
        end

        gradJulia = zeros(nParam); evalGradF(paramVec, gradJulia)
        gradPython = collect(gradPythonMat[i, :])
        sqDiffGrad = sum((gradJulia - gradPython).^2)
        if sqDiffGrad > 1e-5
            @printf("sqDiffGrad = %.3e\n", sqDiffGrad)
            @printf("Does not pass test on gradient\n")
            return false
        end

        if printRes == true
            @printf("sqDiffCost = %.3e\n", sqDiffCost)
            @printf("sqDiffGrad = %.3e\n", sqDiffGrad)
        end
    end

    return true 
end


passTest = compareAgainstPyPesto(Rodas4P(), 1e-12, printRes=false)
if passTest == true
    @printf("Passed test Boehm against PyPesto\n")
else
    @printf("Did not pass test Boehm against PyPesto\n")
end