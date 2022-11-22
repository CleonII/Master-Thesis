#= 
    Builds all models.
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



# Relevant PeTab structs for compuations 
include(joinpath(pwd(), "src", "PeTab_structs.jl"))

# PeTab importer to get cost, grad etc 
include(joinpath(pwd(), "src", "PeTab_importer", "Create_cost_grad_hessian.jl"))

# For converting to SBML 
include(joinpath(pwd(), "src", "SBML", "SBML_to_ModellingToolkit.jl"))

# Alla modeller
modelNames = ["model_Alkan_SciSignal2018",
 "model_Bachmann_MSB2011",
 "model_Beer_MolBioSystems2014",
 "model_Bertozzi_PNAS2020",
 "model_Blasi_CellSystems2016",
 "model_Boehm_JProteomeRes2014",
 "model_Borghans_BiophysChem1997",
 "model_Brannmark_JBC2010",
 "model_Bruno_JExpBot2016",
 "model_Crauste_CellSystems2017",
 "model_Elowitz_Nature2000",
 "model_Fiedler_BMC2016",
 "model_Fujita_SciSignal2010",
 "model_Giordano_Nature2020",
 "model_Isensee_JCB2018",
 "model_Laske_PLOSComputBiol2019",
 "model_Lucarelli_CellSystems2018",
 "model_Okuonghae_ChaosSolitonsFractals2020",
 "model_Oliveira_NatCommun2021",
 "model_Perelson_Science1996",
 "model_Rahman_MBS2016",
 "model_Raimundez_PCB2020",
 "model_SalazarCavazos_MBoC2020",
 "model_Schwen_PONE2014",
 "model_Sneyd_PNAS2002",
 "model_Weber_BMC2015",
 "model_Zhao_QuantBiol2020",
 "model_Zheng_PNAS2012"]

# Bara modeller med explicit rules
#modelNames = ["model_Giordano_Nature2020",
# "model_Laske_PLOSComputBiol2019",
# "model_Rahman_MBS2016",
# "model_SalazarCavazos_MBoC2020",
# "model_Zhao_QuantBiol2020"
#]

for modelName in modelNames
    println("Building model " * modelName)
    dirModel = pwd() * "/Intermediate/PeTab_models/" * modelName * "/"

    mdlFile = dirModel * modelName * ".jl"
    if isfile(mdlFile)
        rm(mdlFile)
    end

    ymodFile = dirModel * modelName * "ObsSdU0.jl"
    if isfile(ymodFile)
        rm(ymodFile)
    end

    peTabModel = setUpPeTabModel(modelName, dirModel);
end
