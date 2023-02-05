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


@testset "Log likelihood values and gradients for benchmark collection" begin

    # Beer model - Numerically challenging gradient as we have callback rootfinding
    pathYML = joinpath(@__DIR__, "..", "tests", "JuliaImport", "model_Beer_MolBioSystems2014", "Beer_MolBioSystems2014.yaml")
    petabModel = readPEtabModel(pathYML, verbose=false, jlFile=true)     
    testLogLikelihoodValue(petabModel, -58622.9145631413, Rodas4P())
    
    # Boehm model 
    pathYML = joinpath(@__DIR__, "..", "tests", "JuliaImport", "model_Boehm_JProteomeRes2014", "Boehm_JProteomeRes2014.yaml")
    petabModel = readPEtabModel(pathYML, verbose=false, jlFile=true)     
    testLogLikelihoodValue(petabModel, 138.22199693517703, Rodas4P())

    # Brännmark model. Model has pre-equlibration criteria so here we test all gradients. Challenging to compute gradients.
    pathYML = joinpath(@__DIR__, "..", "tests", "JuliaImport", "model_Brannmark_JBC2010", "Brannmark_JBC2010.yaml")
    petabModel = readPEtabModel(pathYML, verbose=false, jlFile=true)
    testLogLikelihoodValue(petabModel, 141.889113770537, Rodas4P())

    # Fujita model. Challangeing to compute accurate gradients  
    pathYML = joinpath(@__DIR__, "..", "tests", "JuliaImport", "model_Fujita_SciSignal2010", "Fujita_SciSignal2010.yaml")
    petabModel = readPEtabModel(pathYML, verbose=false, jlFile=true)
    testLogLikelihoodValue(petabModel, -53.08377736998929, Rodas4P())

    # Isensee model. Accurate gradients are computed (but the code takes ages to run with low tolerances)
    pathYML = joinpath(@__DIR__, "..", "tests", "JuliaImport", "model_Isensee_JCB2018", "Isensee_JCB2018.yaml")
    petabModel = readPEtabModel(pathYML, verbose=false, jlFile=true)
    testLogLikelihoodValue(petabModel, 3949.375966548649-4.45299970460275, Rodas4P(), atol=1e-2)

    # Weber model. Challanging as it sensitivity to steady state tolerances 
    pathYML = joinpath(@__DIR__, "..", "tests", "JuliaImport", "model_Weber_BMC2015", "Weber_BMC2015.yaml")
    petabModel = readPEtabModel(pathYML, verbose=false, jlFile=true)
    testLogLikelihoodValue(petabModel, 296.2017922646865, Rodas4P())

end
