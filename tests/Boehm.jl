using ModelingToolkit, DifferentialEquations, DataFrames, CSV 
using ForwardDiff
using ReverseDiff
using StatsBase
using Random
using LinearAlgebra
using Distributions
using Printf


# Functions for solving ODE system 
include(joinpath(pwd(), "Additional_functions", "benchmarkSolvers.jl"))

# Ipopt wrapper 
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "CommonParameterEstimationMethods", "CreateIpoptProb.jl"))

# PeTab importer 
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "PeTabImporter.jl"))

# Additional functions for ODE solver 
include(joinpath(pwd(), "Additional_functions", "additional_tools.jl"))

# HyperCube sampling 
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "LatinHyperCubeParameters", "LatinHyperCubeSampledParameters.jl"))

# For converting to SBML 
include(joinpath(pwd(), "Pipeline_SBMLImporter", "main.jl"))

# For building observable 
include(joinpath(pwd(), "Additional_functions", "createObsFile.jl"))


"""
    testCostGradHess(solver, tol; printRes::Bool=false)

    Compare cost and gradient for Julia PeTab importer against the cost and 
    gradient computed in PyPesto for Boehm model. 
"""
function compareAgainstPyPesto(peTabModel::PeTabModel, solver, tol; printRes::Bool=false)

    evalF, evalGradF, evalHessianApproxF, paramVecEstTmp, lowerBounds, upperBounds, idParam = setUpCostFunc(peTabModel, solver, tol)
    # "Exact" hessian via autodiff 
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


peTabModel = setUpPeTabModel("Boehm_JProteomeRes2014", pwd() * "/tests/Boehm/")

passTest = compareAgainstPyPesto(peTabModel, Rodas4P(), 1e-12, printRes=false)
if passTest == true
    @printf("Passed test Boehm against PyPesto\n")
else
    @printf("Did not pass test Boehm against PyPesto\n")
end