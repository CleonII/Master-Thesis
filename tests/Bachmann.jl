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
using Zygote
using SciMLSensitivity
using Sundials
using YAML
using Test


# Relevant PeTab structs for compuations 
include(joinpath(pwd(), "src", "PeTab_structs.jl"))

# PeTab importer to get cost, grad etc 
include(joinpath(pwd(), "src", "Create_PEtab_model.jl"))

# HyperCube sampling 
include(joinpath(pwd(), "src", "Optimizers", "Lathin_hypercube.jl"))

# For converting to SBML 
include(joinpath(pwd(), "src", "SBML", "SBML_to_ModellingToolkit.jl"))

include(joinpath(pwd(), "tests", "Common.jl"))


"""
    compareAgainstPyPestoBachmann(petabModel::PEtabModel, solver, tol)

    Compare cost and gradient for Julia PeTab importer against the cost and 
    gradient computed in PyPesto for Bachmann model. 
"""
function compareAgainstPyPestoBachmann(petabModel::PEtabModel, solver, tol)

    petabProblem1 = setUpPEtabODEProblem(petabModel, solver, solverAbsTol=tol, solverRelTol=tol, 
                                         sensealgZygote=ForwardDiffSensitivity(), 
                                         odeSolverForwardEquations=CVODE_BDF(), sensealgForwardEquations = ForwardSensitivity(), 
                                         odeSolverAdjoint=solver, solverAdjointAbsTol=tol, solverAdjointRelTol=tol,
                                         sensealgAdjoint=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(true)))

    petabProblem2 = setUpPEtabODEProblem(petabModel, solver, solverAbsTol=tol, solverRelTol=tol, 
                                         sensealgForwardEquations=:AutoDiffForward, odeSolverForwardEquations=solver)


    # Parameter values to test gradient at 
    paramVals = CSV.read(pwd() * "/tests/Bachmann/Params.csv", DataFrame)
    paramMat = paramVals[!, Not([:Id, :SOCS3RNAEqc, :CISRNAEqc])]

    # PyPesto hessian, gradient and cost values 
    costPython = (CSV.read(pwd() * "/tests/Bachmann/Cost.csv", DataFrame))[!, :Cost]
    gradPythonMat = CSV.read(pwd() * "/tests/Bachmann/Grad.csv", DataFrame)
    gradPythonMat = gradPythonMat[!, Not([:Id, :SOCS3RNAEqc, :CISRNAEqc])]
    hessPythonMat = CSV.read(pwd() * "/tests/Bachmann/Hess.csv", DataFrame)
    hessPythonMatCols = names(hessPythonMat)
    hessFilter=findall( x -> occursin("CISRNAEqc", x) || occursin("SOCS3RNAEqc", x), hessPythonMatCols)
    hessPythonMat = hessPythonMat[!, Not(["Id", hessPythonMatCols[hessFilter]...])]
 
    paramEstNames = string.(petabProblem1.θ_estNames)
    # For correct indexing when comparing gradient or when inputing PyPesto vector to Julia 
    iUse = [findfirst(x -> x == paramEstNames[i], names(paramMat)) for i in eachindex(paramEstNames)]
    nParam = ncol(paramMat)

    # For correct indexing when comparing hessian
    namesParamHess = Array{String, 1}(undef,nParam*nParam)
    for k in eachindex(namesParamHess)
        row = 1 + (k-1) % nParam
        col = 1 + (k-1) ÷ nParam
        namesParamHess[k] = paramEstNames[row] * paramEstNames[col]
    end
    iUseHess = [findfirst(x -> x == namesParamHess[i], names(hessPythonMat)) for i in eachindex(namesParamHess)]

    for i in 1:nrow(paramMat)

        p = collect(paramMat[i, iUse])
        referenceCost = costPython[i]
        referenceGradient = collect(gradPythonMat[i, iUse])
        referenceHessian = collect(hessPythonMat[i, iUseHess])

        cost = _testCostGradientOrHessian(petabProblem1, p, cost=true)
        @test cost ≈ referenceCost atol=1e-4
        costZygote = _testCostGradientOrHessian(petabProblem1, p, costZygote=true)
        @test costZygote ≈ referenceCost atol=1e-4
        
        # Test all gradient combinations. Note we test sensitivity equations with and without autodiff 
        gradientAutoDiff = _testCostGradientOrHessian(petabProblem1, p, gradientAutoDiff=true)
        @test norm(gradientAutoDiff - referenceGradient) ≤ 1e-4
        # This currently takes considerble time (hence it does not run for all passes)
        if i == 1
            gradientZygote = _testCostGradientOrHessian(petabProblem1, p, gradientZygote=true)
            @test norm(gradientZygote - referenceGradient) ≤ 1e-4
        end
        gradientAdjoint = _testCostGradientOrHessian(petabProblem1, p, gradientAdjoint=true)
        @test norm(normalize(gradientAdjoint) - normalize((referenceGradient))) ≤ 1e-2
        gradientForwardEquations1 = _testCostGradientOrHessian(petabProblem1, p, gradientForwardEquations=true)
        @test norm(gradientForwardEquations1 - referenceGradient) ≤ 1e-4
        gradientForwardEquations2 = _testCostGradientOrHessian(petabProblem2, p, gradientForwardEquations=true)
        @test norm(gradientForwardEquations2 - referenceGradient) ≤ 1e-4
        
        # Testing "exact" hessian via autodiff 
        hessian = _testCostGradientOrHessian(petabProblem1, p, hessianGN=true)
        @test norm(hessian[:] - referenceHessian) ≤ 1e-2

    end

end


petabModel = readPEtabModel(joinpath(@__DIR__, "Bachmann", "Bachmann_MSB2011.yaml"), forceBuildJuliaFiles=true)
@testset "Against PyPesto : Bachmann" begin 
    compareAgainstPyPestoBachmann(petabModel, Rodas5(), 1e-12)
end
