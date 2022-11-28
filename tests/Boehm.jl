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
using SciMLSensitivity
using Zygote
using BenchmarkTools


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
    gradient computed in PyPesto for Boehm model. 
"""
function compareAgainstPyPesto(peTabModel::PeTabModel, solver, tol; printRes::Bool=false)

    peTabOpt = setUpCostGradHess(peTabModel, solver, tol, sensealg = ForwardDiffSensitivity())

    paramVals = CSV.read(pwd() * "/tests/Boehm/Params.csv", DataFrame)
    paramMat = paramVals[!, Not([:Id, :ratio, :specC17])]
    nParam = ncol(paramMat)

    costPython = (CSV.read(pwd() * "/tests/Boehm/Cost.csv", DataFrame))[!, :Cost]
    gradPythonMat = CSV.read(pwd() * "/tests/Boehm/Grad.csv", DataFrame)
    gradPythonMat = gradPythonMat[!, Not([:Id, :ratio, :specC17])]
    for i in 1:nrow(paramMat)
        paramVec = collect(paramMat[i, :])

        costJulia = peTabOpt.evalF(paramVec)
        sqDiffCost = (costJulia - costPython[i])^2
        if sqDiffCost > 1e-5
            @printf("sqDiffCost = %.3e\n", sqDiffCost)
            @printf("Does not pass test on cost\n")
            return false
        end

        gradJulia = zeros(nParam); peTabOpt.evalGradF(gradJulia, paramVec)
        gradPython = collect(gradPythonMat[i, :])
        sqDiffGrad = sum((gradJulia - gradPython).^2)
        if sqDiffGrad > 1e-5
            @printf("sqDiffGrad = %.3e\n", sqDiffGrad)
            @printf("Does not pass test on gradient\n")
            return false
        end

        costZygote = peTabOpt.evalFZygote(paramVec)
        sqDiffZygote = (costZygote - costPython[i])^2
        if sqDiffZygote > 1e-5
            @printf("sqDiffCostZygote = %.3e\n", sqDiffZygote)
            @printf("Does not pass test on cost for Zygote\n")
            return false
        end

        gradZygote = zeros(nParam)
        peTabOpt.evalGradFZygote(gradZygote, paramVec)
        gradPython = collect(gradPythonMat[i, :])
        sqDiffGradZygote = sum((gradZygote - gradPython).^2)
        if sqDiffGradZygote > 1e-5
            @printf("sqDiffGradZygote = %.3e\n", sqDiffGradZygote)
            @printf("Does not pass test on gradient from Zygote\n")
            return false
        end

        if printRes == true
            @printf("sqDiffCost = %.3e\n", sqDiffCost)
            @printf("sqDiffGrad = %.3e\n", sqDiffGrad)
            @printf("sqDiffCostZygote = %.3e\n", sqDiffZygote)
            @printf("sqDiffGradZygote = %.3e\n", sqDiffGradZygote)
        end
    end

    return true 
end


peTabModel = setUpPeTabModel("Boehm_JProteomeRes2014", pwd() * "/tests/Boehm/", forceBuildJlFile=true)

passTest = compareAgainstPyPesto(peTabModel, Rodas4P(), 1e-12, printRes=true)
if passTest == true
    @printf("Passed test Boehm against PyPesto\n")
else
    @printf("Did not pass test Boehm against PyPesto\n")
end