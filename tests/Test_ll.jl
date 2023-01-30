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
function testLogLikelihoodValue(peTabModel::PeTabModel, 
                                referenceValue::Float64, 
                                solver; absTol=1e-12, relTol=1e-12, atol=1e-3)

    petabProblem = setUpPEtabODEProblem(peTabModel, solver, solverAbsTol=absTol, solverRelTol=relTol)
    cost = petabProblem.computeCost(petabProblem.θ_nominalT)
    costZygote = petabProblem.computeCostZygote(petabProblem.θ_nominalT)
    println("Model : ", peTabModel.modelName)
    @test cost ≈ referenceValue atol=atol    
    @test costZygote ≈ referenceValue atol=atol    
end

@testset "Lok likelihood values benchmark collection" begin

# Bachman model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Bachmann_MSB2011/"
peTabModel = setUpPeTabModel("model_Bachmann_MSB2011", dirModel, verbose=false, forceBuildJlFile=true)
testLogLikelihoodValue(peTabModel, -418.40573341425295, Rodas4P())

# Beer model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Beer_MolBioSystems2014/"
peTabModel = setUpPeTabModel("model_Beer_MolBioSystems2014", dirModel, verbose=false, forceBuildJlFile=true)
testLogLikelihoodValue(peTabModel, -58622.9145631413, Rodas4P())

# Boehm model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Boehm_JProteomeRes2014/"
peTabModel = setUpPeTabModel("model_Boehm_JProteomeRes2014", dirModel, verbose=false, forceBuildJlFile=true)
testLogLikelihoodValue(peTabModel, 138.22199693517703, Rodas4P())

# Brännmark model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Brannmark_JBC2010/"
peTabModel = setUpPeTabModel("model_Brannmark_JBC2010", dirModel, verbose=false, forceBuildJlFile=true)
testLogLikelihoodValue(peTabModel, 141.889113770537, Rodas4P())

# Bruno model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Bruno_JExpBot2016/"
peTabModel = setUpPeTabModel("model_Bruno_JExpBot2016", dirModel, verbose=false, forceBuildJlFile=true)
testLogLikelihoodValue(peTabModel, -46.688176988431806, Rodas4P())

# Crauste model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Crauste_CellSystems2017/"
peTabModel = setUpPeTabModel("model_Crauste_CellSystems2017", dirModel, verbose=false, forceBuildJlFile=true)
testLogLikelihoodValue(peTabModel, 190.96521897435176, Rodas4P(), atol=1e-2)

# Elowitz model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Elowitz_Nature2000/"
peTabModel = setUpPeTabModel("model_Elowitz_Nature2000", dirModel, verbose=false, forceBuildJlFile=true)
testLogLikelihoodValue(peTabModel, -63.20279991419332, Rodas4P())

# Fiedler model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Fiedler_BMC2016/"
peTabModel = setUpPeTabModel("model_Fiedler_BMC2016", dirModel, verbose=false, forceBuildJlFile=true)
testLogLikelihoodValue(peTabModel, -58.58390161681, Rodas4(), absTol=1e-11, relTol=1e-11)

# Fujita model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Fujita_SciSignal2010/"
peTabModel = setUpPeTabModel("model_Fujita_SciSignal2010", dirModel, verbose=false, forceBuildJlFile=true)
testLogLikelihoodValue(peTabModel, -53.08377736998929, Rodas4P())

# Isensee model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Isensee_JCB2018/"
peTabModel = setUpPeTabModel("model_Isensee_JCB2018", dirModel, verbose=false, forceBuildJlFile=true)
testLogLikelihoodValue(peTabModel, 3949.375966548649-4.45299970460275, Rodas4P(), atol=1e-2)

# Lucarelli 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Lucarelli_CellSystems2018/"
peTabModel = setUpPeTabModel("model_Lucarelli_CellSystems2018", dirModel, verbose=false, forceBuildJlFile=true)
testLogLikelihoodValue(peTabModel, 1681.6059879426584, Rodas4P())

# Schwen model
dirModel = pwd() * "/Intermediate/PeTab_models/model_Schwen_PONE2014/"
peTabModel = setUpPeTabModel("model_Schwen_PONE2014", dirModel, verbose=false, forceBuildJlFile=true)
testLogLikelihoodValue(peTabModel, 943.9992988598723-12.519137073132825, Rodas4P())

# Sneyd model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Sneyd_PNAS2002/"
peTabModel = setUpPeTabModel("model_Sneyd_PNAS2002", dirModel, verbose=false, forceBuildJlFile=true)
testLogLikelihoodValue(peTabModel, -319.79177818768756, Rodas4P())

# Weber model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Weber_BMC2015/"
peTabModel = setUpPeTabModel("model_Weber_BMC2015", dirModel, verbose=false, forceBuildJlFile=true)
testLogLikelihoodValue(peTabModel, 296.2017922646865, Rodas4P())

# Zheng model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Zheng_PNAS2012/"
peTabModel = setUpPeTabModel("model_Zheng_PNAS2012", dirModel, verbose=false, forceBuildJlFile=true)
testLogLikelihoodValue(peTabModel, -278.33353271001477, Rodas4P())

end