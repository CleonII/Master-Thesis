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
    gradient computed in PyPesto for Bachman model. 
"""
function compareAgainstPyPesto(solver, tol; printRes::Bool=false)

        modelName = "model_Bachmann_MSB2011"
        evalObs = Bachmann_MSB2011
        evalU0 = Bachmann_MSB2011_t0!
        evalSd = Bachmann_MSB2011_sd!

        # Parameter values to test gradient at 
        paramVals = CSV.read(pwd() * "/tests/Bachman/Params.csv", DataFrame)
        paramMat = paramVals[!, Not([:Id, :SOCS3RNAEqc, :CISRNAEqc])]

        # PyPesto gradient and cost values 
        costPython = (CSV.read(pwd() * "/tests/Bachman/Cost.csv", DataFrame))[!, :Cost]
        gradPythonMat = CSV.read(pwd() * "/tests/Bachman/Grad.csv", DataFrame)
        gradPythonMat = gradPythonMat[!, Not([:Id, :SOCS3RNAEqc, :CISRNAEqc])]


        evalF, evalGradF, evalHessianApproxF, paramVecEstTmp, lowerBounds, upperBounds, idParam = setUpCostFunc(modelName, evalObs, evalU0, evalSd, solver, tol)
        # "Exact" hessian via autodiff (should only be used for smaller models)
        evalH = (hessianMat, paramVec) -> begin hessianMat .= Symmetric(ForwardDiff.hessian(evalF, paramVec)) end

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


passTest = compareAgainstPyPesto(Rodas4P(), 1e-12, printRes=false)
if passTest == true
    @printf("Passed test Boehm against PyPesto\n")
else
    @printf("Did not pass test Boehm against PyPesto\n")
end
