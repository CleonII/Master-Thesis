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

"""
    testCostGradHess(solver, tol; printRes::Bool=false)

    Compare cost and gradient for Julia PeTab importer against the cost and 
    gradient computed in PyPesto for Bachman model. 
"""
function compareAgainstPyPesto(peTabModel::PeTabModel, solver, tol; printRes::Bool=false)

    evalF, evalGradF, evalHessianApproxF, paramVecEstTmp, lowerBounds, upperBounds, idParam = setUpCostGradHess(peTabModel, solver, tol)
    # "Exact" hessian via autodiff 
    evalH = (hessianMat, paramVec) -> begin hessianMat .= Symmetric(ForwardDiff.hessian(evalF, paramVec)) end

    # Parameter values to test gradient at 
    paramVals = CSV.read(pwd() * "/tests/Bachman/Params.csv", DataFrame)
    paramMat = paramVals[!, Not([:Id, :SOCS3RNAEqc, :CISRNAEqc])]

    # PyPesto gradient and cost values 
    costPython = (CSV.read(pwd() * "/tests/Bachman/Cost.csv", DataFrame))[!, :Cost]
    gradPythonMat = CSV.read(pwd() * "/tests/Bachman/Grad.csv", DataFrame)
    gradPythonMat = gradPythonMat[!, Not([:Id, :SOCS3RNAEqc, :CISRNAEqc])]

    # For correct indexing when comparing gradient or when inputing PyPesto vector to Julia 
    iUse = [findfirst(x -> x == idParam[i], names(paramMat)) for i in eachindex(idParam)]

    for i in 1:nrow(paramMat)
        
        paramVec = collect(paramMat[i, :])[iUse]

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
