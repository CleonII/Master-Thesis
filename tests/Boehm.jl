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
    compareAgainstPyPestoBoehm(petabModel::PEtabModel, solver, tol)
    
    Compare cost and gradient for Julia PeTab importer against the cost and 
    gradient computed in PyPesto for Boehm model. 
"""
function compareAgainstPyPestoBoehm(petabModel::PEtabModel, solver, tol)

    petabProblem1 = setUpPEtabODEProblem(petabModel, solver, solverAbsTol=tol, solverRelTol=tol, 
                                         sensealgZygote=ForwardDiffSensitivity(), 
                                         odeSolverForwardEquations=CVODE_BDF(), sensealgForwardEquations = ForwardSensitivity(), 
                                         odeSolverAdjoint=solver, solverAdjointAbsTol=tol, solverAdjointRelTol=tol,
                                         sensealgAdjoint=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(true)))

    petabProblem2 = setUpPEtabODEProblem(petabModel, solver, solverAbsTol=tol, solverRelTol=tol, 
                                         sensealgForwardEquations=:AutoDiffForward, odeSolverForwardEquations=solver)

    paramVals = CSV.read(pwd() * "/tests/Boehm/Params.csv", DataFrame)
    paramMat = paramVals[!, Not([:Id, :ratio, :specC17])]

    costPython = (CSV.read(pwd() * "/tests/Boehm/Cost.csv", DataFrame))[!, :Cost]
    gradPythonMat = CSV.read(pwd() * "/tests/Boehm/Grad.csv", DataFrame)
    gradPythonMat = gradPythonMat[!, Not([:Id, :ratio, :specC17])]
    hessPythonMat = CSV.read(pwd() * "/tests/Boehm/Hess.csv", DataFrame)
    hessPythonMatCols = names(hessPythonMat)
    hessFilter=findall( x -> occursin("ratio", x) || occursin("specC17", x), hessPythonMatCols)
    hessPythonMat = hessPythonMat[!, Not(["Id", hessPythonMatCols[hessFilter]...])]
    
    for i in 1:nrow(paramMat)
        
        p = collect(paramMat[i, :])
        referenceCost = costPython[i]
        referenceGradient = collect(gradPythonMat[i, :])
        referenceHessian = collect(hessPythonMat[i, :])

        cost = _testCostGradientOrHessian(petabProblem1, p, cost=true)
        @test cost ≈ referenceCost atol=1e-4
        costZygote = _testCostGradientOrHessian(petabProblem1, p, costZygote=true)
        @test costZygote ≈ referenceCost atol=1e-4

        # Test all gradient combinations. Note we test sensitivity equations with and without autodiff 
        gradientAutoDiff = _testCostGradientOrHessian(petabProblem1, p, gradientAutoDiff=true)
        @test norm(gradientAutoDiff - referenceGradient) ≤ 1e-4
        gradientZygote = _testCostGradientOrHessian(petabProblem1, p, gradientZygote=true)
        @test norm(gradientZygote - referenceGradient) ≤ 1e-4
        gradientAdjoint = _testCostGradientOrHessian(petabProblem1, p, gradientAdjoint=true)
        @test norm(normalize(gradientAdjoint) - normalize((referenceGradient))) ≤ 1e-2
        gradientForwardEquations1 = _testCostGradientOrHessian(petabProblem1, p, gradientForwardEquations=true)
        @test norm(gradientForwardEquations1 - referenceGradient) ≤ 1e-4
        gradientForwardEquations2 = _testCostGradientOrHessian(petabProblem2, p, gradientForwardEquations=true)
        @test norm(gradientForwardEquations2 - referenceGradient) ≤ 1e-4

        # Testing "exact" hessian via autodiff 
        hessian = _testCostGradientOrHessian(petabProblem1, p, hessianGN=true)
        @test norm(hessian[:] - referenceHessian) ≤ 1e-3
    
    end
end


petabModel = readPEtabModel(joinpath(@__DIR__, "Boehm", "Boehm_JProteomeRes2014.yaml"), forceBuildJuliaFiles=true)
@testset "Against PyPesto : Boehm" begin 
    compareAgainstPyPestoBoehm(petabModel, Rodas4P(), 1e-9)
end
