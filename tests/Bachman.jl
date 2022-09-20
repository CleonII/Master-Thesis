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
    gradient computed in PyPesto for Bachman model. 
"""
function compareAgainstPyPesto(peTabModel::PeTabModel, solver, tol; printRes::Bool=false)

    evalF, evalGradF, evalHessianApproxF, paramVecEstTmp, lowerBounds, upperBounds, idParam = setUpCostFunc(peTabModel, solver, tol)
    # "Exact" hessian via autodiff 
    evalH = (hessianMat, paramVec) -> begin hessianMat .= Symmetric(ForwardDiff.hessian(evalF, paramVec)) end

    # Parameter values to test gradient at 
    paramVals = CSV.read(pwd() * "/tests/Bachman/Params.csv", DataFrame)
    paramMat = paramVals[!, Not([:Id, :SOCS3RNAEqc, :CISRNAEqc])]

    # PyPesto gradient and cost values 
    costPython = (CSV.read(pwd() * "/tests/Bachman/Cost.csv", DataFrame))[!, :Cost]
    gradPythonMat = CSV.read(pwd() * "/tests/Bachman/Grad.csv", DataFrame)
    gradPythonMat = gradPythonMat[!, Not([:Id, :SOCS3RNAEqc, :CISRNAEqc])]

    # For correct indexing when comparing gradient 
    iUse = [findfirst(x -> x == idParam[i], names(paramMat)) for i in eachindex(idParam)]

    for i in 1:nrow(paramMat)
        
        paramVec = collect(paramMat[i, :])

        # Test cost 
        costJulia = evalF(paramVec)
        sqDiffCost = (costJulia - costPython[i])^2
        if sqDiffCost > 1e-5
            @printf("sqDiffCost = %.3e\n", sqDiffCost)
            @printf("Does not pass test on cost\n")
            return false
        end

        grad = zeros(length(paramVec))
        gradPython = collect(gradPythonMat[i, :])
        evalGradF(paramVec, grad)
        sqDiffGrad = sum((grad - gradPython[iUse]).^2)
        if sqDiffGrad > 1e-5
            @printf("sqDiffGrad = %.3e\n", sqDiffGrad)
            @printf("Does not pass test on grad\n")
            return false
        end

        if printRes == true
            @printf("sqDiffCost = %.3e\n", sqDiffCost)
            @printf("sqDiffGrad = %.3e\n", sqDiffGrad)
        end
    end

    return true
end


peTabModel = setUpPeTabModel("Bachmann_MSB2011", pwd() * "/tests/Bachman/")

passTest = compareAgainstPyPesto(peTabModel, Rodas4P(), 1e-12, printRes=true)
if passTest == true
    @printf("Passed test Boehm against PyPesto\n")
else
    @printf("Did not pass test Boehm against PyPesto\n")
end
