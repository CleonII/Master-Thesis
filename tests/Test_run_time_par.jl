#= 
    Using the Bachman and Brännmark model this file tests that our parallellisation scales well, 
    and that we compute accurate values.

    Briefly, using Distributed.jl we parallellise over different experimental conditions.
    The Bachman models has 36 different experimental conditions (we solve the ODE 36 times). 
    If our parallellisation works well we should se substantial speedup for the cost, and 
    especially for the gradient, close to the ideal reduction of factor 2 when using 
    two processes. 
    Similarly, the Brännmark model has 8 different experimental conditions, and we can further 
    compute the hessian and an associated hessian approximation and check corresponding accuracy 
    in for each computation.
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
using SciMLSensitivity
using Zygote
using Symbolics
using Sundials
using FiniteDifferences
using BenchmarkTools
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


function testAccuracyDistributed(petabModel::PEtabModel, 
                                 solver, 
                                 tol::Float64;
                                 solverSSRelTol=1e-10, 
                                 solverSSAbsTol=1e-10,
                                 sensealgSS=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(false)),
                                 testHessian::Bool=false)

    println("Model = ", petabModel.modelName)                                 

    petabProblemSeq = setUpPEtabODEProblem(petabModel, solver, solverAbsTol=tol, solverRelTol=tol, 
                                           sensealgForwardEquations=:AutoDiffForward, odeSolverForwardEquations=solver, 
                                           odeSolverAdjoint=solver, solverAdjointAbsTol=tol, solverAdjointRelTol=tol,
                                           sensealgAdjoint=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(true)), 
                                           sensealgAdjointSS=sensealgSS,
                                           solverSSRelTol=solverSSRelTol, solverSSAbsTol=solverSSAbsTol, numberOfprocesses=1)

    petabProblemPar = setUpPEtabODEProblem(petabModel, solver, solverAbsTol=tol, solverRelTol=tol, 
                                           sensealgForwardEquations=:AutoDiffForward, odeSolverForwardEquations=solver, 
                                           odeSolverAdjoint=solver, solverAdjointAbsTol=tol, solverAdjointRelTol=tol,
                                           sensealgAdjoint=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(true)), 
                                           sensealgAdjointSS=sensealgSS,
                                           solverSSRelTol=solverSSRelTol, solverSSAbsTol=solverSSAbsTol, numberOfprocesses=2)       

    θ_est = petabProblemSeq.θ_nominalT               
    costSeq = petabProblemSeq.computeCost(θ_est)                                                 
    costPar = petabProblemPar.computeCost(θ_est)    
    @test costSeq ≈ costPar atol=1e-8
    
    gradientSeqAuto, gradientParAuto = zeros(length(θ_est)), zeros(length(θ_est))
    petabProblemSeq.computeGradientAutoDiff(gradientSeqAuto, θ_est)          
    petabProblemPar.computeGradientAutoDiff(gradientParAuto, θ_est)          
    normDiffAuto = norm(gradientSeqAuto - gradientParAuto)
    @test normDiffAuto ≤ 1e-8

    # To check that the implementation is reasonable the run time is also tested for the gradient 
    bGradSeq = 0.0 
    bGradPar = 0.0 
    for i in 1:20
        bGradSeq += @elapsed petabProblemSeq.computeGradientAutoDiff(gradientSeqAuto, θ_est)          
        sleep(0.1)
        bGradPar += @elapsed petabProblemPar.computeGradientAutoDiff(gradientParAuto, θ_est)          
        sleep(0.1)
    end
    @test bGradPar ≤ bGradSeq * 0.9
    @printf("Ratio run time = %.2f\n", bGradPar / bGradSeq)                                           

    gradientSeqAdjoint, gradientParAdjoint = zeros(length(θ_est)), zeros(length(θ_est))
    petabProblemSeq.computeGradientAdjoint(gradientSeqAdjoint, θ_est)          
    petabProblemPar.computeGradientAdjoint(gradientParAdjoint, θ_est)          
    normDiffAdjoint = norm(gradientSeqAdjoint - gradientParAdjoint)
    @test normDiffAdjoint ≤ 1e-8

    gradientSeqForward, gradientParForward = zeros(length(θ_est)), zeros(length(θ_est))
    petabProblemSeq.computeGradientForwardEquations(gradientSeqForward, θ_est)          
    petabProblemPar.computeGradientForwardEquations(gradientParForward, θ_est)          
    normDiffForward = norm(gradientSeqForward - gradientParForward)
    @test normDiffForward ≤ 1e-8

    if testHessian == true
        hessianAutoDiffSeq, hessianAutoDiffPar = zeros(length(θ_est), length(θ_est)), zeros(length(θ_est), length(θ_est))
        petabProblemSeq.computeHessian(hessianAutoDiffSeq, θ_est)          
        petabProblemPar.computeHessian(hessianAutoDiffPar, θ_est)          
        normDiffHessian = norm(hessianAutoDiffSeq - hessianAutoDiffPar)
        @test normDiffHessian ≤ 1e-8

        hessianBlockDiffSeq, hessianBlockDiffPar = zeros(length(θ_est), length(θ_est)), zeros(length(θ_est), length(θ_est))
        petabProblemSeq.computeHessianGN(hessianBlockDiffSeq, θ_est)          
        petabProblemPar.computeHessianGN(hessianBlockDiffPar, θ_est)          
        normDiffHessianBlock = norm(hessianBlockDiffSeq - hessianBlockDiffPar)
        @test normDiffHessianBlock ≤ 1e-8

        hessianGNSeq, hessianGNPar = zeros(length(θ_est), length(θ_est)), zeros(length(θ_est), length(θ_est))
        petabProblemSeq.computeHessianGN(hessianGNSeq, θ_est)          
        petabProblemPar.computeHessianGN(hessianGNPar, θ_est)          
        normDiffGNHessian = norm(hessianGNSeq - hessianGNPar)
        @test normDiffGNHessian ≤ 1e-8
    end
end



@testset verbose = true "Test distributed" begin

    removeAllProcs()
    pathYML = joinpath(@__DIR__, "..", "Intermediate", "PeTab_models", "model_Bachmann_MSB2011", "Bachmann_MSB2011.yaml")
    petabModel = readPEtabModel(pathYML, verbose=false, forceBuildJuliaFiles=true)
    addprocs(1, exeflags="--project=.")
    @testset "Bachman : Accuracy and run-time distributed" begin
        testAccuracyDistributed(petabModel, Rodas5(), 1e-8)
    end

    removeAllProcs()
    pathYML = joinpath(@__DIR__, "..", "Intermediate", "PeTab_models", "model_Brannmark_JBC2010", "Brannmark_JBC2010.yaml")
    petabModel = readPEtabModel(pathYML, verbose=false, forceBuildJuliaFiles=true)
    addprocs(1, exeflags="--project=.")
    @testset "Brännmark : Accuracy and run-time distributed" begin
        testAccuracyDistributed(petabModel, Rodas5(), 1e-8)
    end
end
