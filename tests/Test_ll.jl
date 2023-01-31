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
using SciMLSensitivity
using Zygote
using Symbolics
using Sundials
using FiniteDifferences
using Test


# Relevant PeTab structs for compuations 
include(joinpath(pwd(), "src", "PeTab_structs.jl"))

# PeTab importer to get cost, grad etc 
include(joinpath(pwd(), "src", "Process_PEtab_files", "Create_cost_grad_hessian.jl"))

# HyperCube sampling 
include(joinpath(pwd(), "src", "Optimizers", "Lathin_hypercube.jl"))

# For converting to SBML 
include(joinpath(pwd(), "src", "SBML", "SBML_to_ModellingToolkit.jl"))


# Used to test cost-value at the nominal parameter value 
function testLogLikelihoodValue(petabModel::PEtabModel, 
                                referenceValue::Float64, 
                                solver; absTol=1e-12, relTol=1e-12, atol=1e-3)

    petabProblem = setUpPEtabODEProblem(petabModel, solver, solverAbsTol=absTol, solverRelTol=relTol)
    cost = petabProblem.computeCost(petabProblem.θ_nominalT)
    costZygote = petabProblem.computeCostZygote(petabProblem.θ_nominalT)
    println("Model : ", petabModel.modelName)
    @test cost ≈ referenceValue atol=atol    
    @test costZygote ≈ referenceValue atol=atol    
end


function testGradientFiniteDifferences(petabModel::PEtabModel, solver, tol::Float64; 
                                       checkAdjoint::Bool=false, 
                                       solverForwardEq=CVODE_BDF(),
                                       checkForwardEquations::Bool=false, 
                                       testTol::Float64=1e-3, 
                                       sensealgSS=SteadyStateAdjoint(), 
                                       sensealgAdjoint=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(true)),
                                       solverSSRelTol=1e-8, solverSSAbsTol=1e-10)

    # Testing the gradient via finite differences 
    petabProblem1 = setUpPEtabODEProblem(petabModel, solver, solverAbsTol=tol, solverRelTol=tol, 
                                         sensealgForwardEquations=:AutoDiffForward, odeSolverForwardEquations=solver, 
                                         odeSolverAdjoint=solver, solverAdjointAbsTol=tol, solverAdjointRelTol=tol,
                                         sensealgAdjoint=sensealgAdjoint, 
                                         solverSSRelTol=solverSSRelTol, solverSSAbsTol=solverSSAbsTol,
                                         sensealgAdjointSS=sensealgSS)
    petabProblem2 = setUpPEtabODEProblem(petabModel, solver, solverAbsTol=tol, solverRelTol=tol, 
                                         solverSSRelTol=solverSSRelTol, solverSSAbsTol=solverSSAbsTol,
                                         sensealgForwardEquations=ForwardSensitivity(), odeSolverForwardEquations=solverForwardEq)                                        
    θ_use = petabProblem1.θ_nominalT

    gradientFinite = FiniteDifferences.grad(central_fdm(5, 1), petabProblem1.computeCost, θ_use)[1]
    gradientForward = zeros(length(θ_use))
    petabProblem1.computeGradientAutoDiff(gradientForward, θ_use)                                       
    @test norm(gradientFinite - gradientForward) ≤ testTol

    if checkForwardEquations == true
        gradientForwardEquations1 = zeros(length(θ_use))
        gradientForwardEquations2 = zeros(length(θ_use))
        petabProblem1.computeGradientForwardEquations(gradientForwardEquations1, θ_use)
        petabProblem2.computeGradientForwardEquations(gradientForwardEquations2, θ_use)
        @test norm(gradientFinite - gradientForwardEquations1) ≤ testTol
        @test norm(gradientFinite - gradientForwardEquations2) ≤ testTol
    end

    if checkAdjoint == true
        gradientAdjoint = zeros(length(θ_use))
        petabProblem1.computeGradientAdjoint(gradientAdjoint, θ_use)
        @test norm(gradientFinite - gradientAdjoint) ≤ testTol
    end
end


@testset "Log likelihood values and gradients for benchmark collection" begin

    # Bachman model 
    dirModel = pwd() * "/Intermediate/PeTab_models/model_Bachmann_MSB2011/"
    petabModel = readPEtabModel("model_Bachmann_MSB2011", dirModel, verbose=false, forceBuildJlFile=true)
    testLogLikelihoodValue(petabModel, -418.40573341425295, Rodas4P())
    testGradientFiniteDifferences(petabModel, Rodas5(), 1e-8)

    # Beer model - Numerically challenging gradient as we have callback rootfinding
    dirModel = pwd() * "/Intermediate/PeTab_models/model_Beer_MolBioSystems2014/"
    petabModel = readPEtabModel("model_Beer_MolBioSystems2014", dirModel, verbose=false, forceBuildJlFile=true)
    testLogLikelihoodValue(petabModel, -58622.9145631413, Rodas4P())
    testGradientFiniteDifferences(petabModel, Rodas4P(), 1e-8, testTol=1e-1)           
    
    # Boehm model 
    dirModel = pwd() * "/Intermediate/PeTab_models/model_Boehm_JProteomeRes2014/"
    petabModel = readPEtabModel("model_Boehm_JProteomeRes2014", dirModel, verbose=false, forceBuildJlFile=true)
    testLogLikelihoodValue(petabModel, 138.22199693517703, Rodas4P())
    testGradientFiniteDifferences(petabModel, Rodas4P(), 1e-8)

    # Brännmark model. Model has pre-equlibration criteria so here we test all gradients. Challenging to compute gradients.
    dirModel = pwd() * "/Intermediate/PeTab_models/model_Brannmark_JBC2010/"
    petabModel = readPEtabModel("model_Brannmark_JBC2010", dirModel, verbose=false, forceBuildJlFile=true)
    testLogLikelihoodValue(petabModel, 141.889113770537, Rodas4P())
    testGradientFiniteDifferences(petabModel, Rodas4P(), 1e-8, checkAdjoint=true, checkForwardEquations=true, testTol=1e-2, sensealgSS=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(false)), sensealgAdjoint=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(false)))

    # Bruno model. Has conditions-specific parameters, hence we test all gradients 
    dirModel = pwd() * "/Intermediate/PeTab_models/model_Bruno_JExpBot2016/"
    petabModel = readPEtabModel("model_Bruno_JExpBot2016", dirModel, verbose=false, forceBuildJlFile=true)
    testLogLikelihoodValue(petabModel, -46.688176988431806, Rodas4P())
    testGradientFiniteDifferences(petabModel, Rodas4P(), 1e-9, checkAdjoint=true, checkForwardEquations=true)

    # Crauste model. The model is numerically challanging and computing a gradient via Finite-differences is not 
    # possible 
    dirModel = pwd() * "/Intermediate/PeTab_models/model_Crauste_CellSystems2017/"
    petabModel = readPEtabModel("model_Crauste_CellSystems2017", dirModel, verbose=false, forceBuildJlFile=true)
    testLogLikelihoodValue(petabModel, 190.96521897435176, Rodas4P(), atol=1e-2)

    # Elowitz model 
    dirModel = pwd() * "/Intermediate/PeTab_models/model_Elowitz_Nature2000/"
    petabModel = readPEtabModel("model_Elowitz_Nature2000", dirModel, verbose=false, forceBuildJlFile=true)
    testLogLikelihoodValue(petabModel, -63.20279991419332, Rodas4P())
    testGradientFiniteDifferences(petabModel, Rodas4P(), 1e-8)

    # Fiedler model 
    dirModel = pwd() * "/Intermediate/PeTab_models/model_Fiedler_BMC2016/"
    petabModel = readPEtabModel("model_Fiedler_BMC2016", dirModel, verbose=false, forceBuildJlFile=true)
    testLogLikelihoodValue(petabModel, -58.58390161681, Rodas4(), absTol=1e-11, relTol=1e-11)
    testGradientFiniteDifferences(petabModel, Rodas4P(), 1e-9)

    # Fujita model. Challangeing to compute accurate gradients  
    dirModel = pwd() * "/Intermediate/PeTab_models/model_Fujita_SciSignal2010/"
    petabModel = readPEtabModel("model_Fujita_SciSignal2010", dirModel, verbose=false, forceBuildJlFile=true)
    testLogLikelihoodValue(petabModel, -53.08377736998929, Rodas4P())
    testGradientFiniteDifferences(petabModel, Rodas4P(), 1e-12, testTol=1e-2)

    # Isensee model. Accurate gradients are computed (but the code takes ages to run with low tolerances)
    dirModel = pwd() * "/Intermediate/PeTab_models/model_Isensee_JCB2018/"
    petabModel = readPEtabModel("model_Isensee_JCB2018", dirModel, verbose=false, forceBuildJlFile=true)
    testLogLikelihoodValue(petabModel, 3949.375966548649-4.45299970460275, Rodas4P(), atol=1e-2)
    testGradientFiniteDifferences(petabModel, Rodas4P(), 1e-8, testTol=1e-2)

    # Lucarelli 
    dirModel = pwd() * "/Intermediate/PeTab_models/model_Lucarelli_CellSystems2018/"
    petabModel = readPEtabModel("model_Lucarelli_CellSystems2018", dirModel, verbose=false, forceBuildJlFile=true)
    testLogLikelihoodValue(petabModel, 1681.6059879426584, Rodas4P())
    testGradientFiniteDifferences(petabModel, Rodas4P(), 1e-9)

    # Schwen model. Model has priors so here we want to test all gradients
    dirModel = pwd() * "/Intermediate/PeTab_models/model_Schwen_PONE2014/"
    petabModel = readPEtabModel("model_Schwen_PONE2014", dirModel, verbose=false, forceBuildJlFile=true)
    testLogLikelihoodValue(petabModel, 943.9992988598723-12.519137073132825, Rodas4P())
    testGradientFiniteDifferences(petabModel, Rodas4P(), 1e-8, checkAdjoint=true, checkForwardEquations=true)

    # Sneyd model 
    dirModel = pwd() * "/Intermediate/PeTab_models/model_Sneyd_PNAS2002/"
    petabModel = readPEtabModel("model_Sneyd_PNAS2002", dirModel, verbose=false, forceBuildJlFile=true)
    testLogLikelihoodValue(petabModel, -319.79177818768756, Rodas4P())
    testGradientFiniteDifferences(petabModel, Rodas4P(), 1e-8)

    # Weber model. Challanging as it sensitivity to steady state tolerances 
    dirModel = pwd() * "/Intermediate/PeTab_models/model_Weber_BMC2015/"
    petabModel = readPEtabModel("model_Weber_BMC2015", dirModel, verbose=false, forceBuildJlFile=true)
    testLogLikelihoodValue(petabModel, 296.2017922646865, Rodas4P())
    testGradientFiniteDifferences(petabModel, Rodas4P(), 1e-12, testTol=1e-2, solverSSRelTol=1e-13, solverSSAbsTol=1e-15)

    # Zheng model 
    dirModel = pwd() * "/Intermediate/PeTab_models/model_Zheng_PNAS2012/"
    petabModel = readPEtabModel("model_Zheng_PNAS2012", dirModel, verbose=false, forceBuildJlFile=true)
    testLogLikelihoodValue(petabModel, -278.33353271001477, Rodas4P())
    testGradientFiniteDifferences(petabModel, Rodas4P(), 1e-8)

end
