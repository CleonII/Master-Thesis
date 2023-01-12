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
using Sundials


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


passTest = true


# Bachman model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Bachmann_MSB2011/"
peTabModel = setUpPeTabModel("model_Bachmann_MSB2011", dirModel, verbose=false, forceBuildJlFile=true)
peTabOpt = setUpCostGradHess(peTabModel, Rodas4P(), 1e-12)
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed)
costZygote = peTabOpt.evalFZygote(peTabOpt.paramVecTransformed)
diff = cost + (418.40573341425295)
diffZygote = costZygote + (418.40573341425295)
if diff > 1e-3
    println("Does not pass ll-test for Bachman model")
    passTest = false
elseif diffZygote > 1e-3
    println("Does not pass ll-test for Zygote Bachman model")
    passTest = false
end


# Beer model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Beer_MolBioSystems2014/"
peTabModel = setUpPeTabModel("model_Beer_MolBioSystems2014", dirModel, verbose=false, forceBuildJlFile=true)
peTabOpt = setUpCostGradHess(peTabModel, Rodas4P(), 1e-12)
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed)
costZygote = peTabOpt.evalFZygote(peTabOpt.paramVecTransformed)
diff = cost + (58622.9145631413)
diffZygote = costZygote + (58622.9145631413)
if diff > 1e-2
    println("Does not pass ll-test for Beer model")
    passTest = false
elseif diffZygote > 1e-2
    println("Does not pass ll-test for Zygote Beer model")
    passTest = false
end


# Boehm model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Boehm_JProteomeRes2014/"
peTabModel = setUpPeTabModel("model_Boehm_JProteomeRes2014", dirModel, verbose=false, forceBuildJlFile=true)
peTabOpt = setUpCostGradHess(peTabModel, Rodas4P(), 1e-12)
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed)
costZygote = peTabOpt.evalFZygote(peTabOpt.paramVecTransformed)
diff = cost + (-138.22199693517703)
diffZygote = costZygote + (-138.22199693517703)
if diff > 1e-3
    println("Does not pass ll-test for Boehm model")
    passTest = false
elseif diffZygote > 1e-3
    println("Does not pass ll-test for Zygote Boehm model")
    passTest = false
end


# Brännmark model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Brannmark_JBC2010/"
peTabModel = setUpPeTabModel("model_Brannmark_JBC2010", dirModel, verbose=false, forceBuildJlFile=true)
peTabOpt = setUpCostGradHess(peTabModel, Rodas4P(), 1e-12)
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed)
costZygote = peTabOpt.evalFZygote(peTabOpt.paramVecTransformed)
diff = abs(cost + (-141.889113770537))
diffZygote = abs(costZygote + (-141.889113770537))
if diff > 1e-3
    println("Does not pass ll-test for Brännmark model")
    passTest = false
elseif diffZygote > 1e-3
    println("Does not pass ll-test for Zygote Brännmark model")
    passTest = false
end
# Check we get correct gradient for a pre-eq simulation model 
peTabOpt = setUpCostGradHess(peTabModel, Rodas5(), 1e-9, adjSolver=Rodas5(), adjTol=1e-9, sensealg=QuadratureAdjoint(autojacvec=ReverseDiffVJP(false)),
                             sensealgForward = ForwardSensitivity(), solverForward=CVODE_BDF(),
                             adjSensealg=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(false)), 
                             adjSensealgSS=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(false)))
peTabOptAlt = setUpCostGradHess(peTabModel, Rodas5(), 1e-9, adjSolver=Rodas5(), adjTol=1e-9, sensealg=QuadratureAdjoint(autojacvec=ReverseDiffVJP(false)),
                                sensealgForward = :AutoDiffForward, solverForward=Rodas5())
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed)                             
gradFor = zeros(length(peTabOpt.paramVecTransformed))
gradForEqAlt = zeros(length(peTabOpt.paramVecTransformed))
gradAdj = zeros(length(peTabOpt.paramVecTransformed))
gradForEq = zeros(length(peTabOpt.paramVecTransformed))
peTabOpt.evalGradF(gradFor, peTabOpt.paramVecTransformed)
peTabOpt.evalGradFAdjoint(gradAdj, peTabOpt.paramVecTransformed)
peTabOpt.evalGradFForwardEq(gradForEq, peTabOpt.paramVecTransformed)
peTabOptAlt.evalGradFForwardEq(gradForEqAlt, peTabOpt.paramVecTransformed)
if sum((gradFor - gradAdj).^2) > 1e-4
    println("Does not pass prior test for adjoint sensitivity")
    passTest=false
end
if sum((gradFor - gradForEq).^2) > 1e-4
    println("Does not pass prior test for forward sensitiviity equations")
    passTest=false
end
if sum((gradFor - gradForEqAlt).^2) > 1e-4
    println("Does not pass prior test for forward sensitivity equations via autodiff")
    passTest=false
end


# Bruno model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Bruno_JExpBot2016/"
peTabModel = setUpPeTabModel("model_Bruno_JExpBot2016", dirModel, verbose=false, forceBuildJlFile=true)
peTabOpt = setUpCostGradHess(peTabModel, Rodas4P(), 1e-12, sensealg=ForwardDiffSensitivity())
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed)
costZygote = peTabOpt.evalFZygote(peTabOpt.paramVecTransformed)
diff = abs(cost + (46.688176988431806))
diffZygote = abs(costZygote + (46.688176988431806))
if diff > 1e-3
    println("Does not pass ll-test for Bruno model")
    passTest = false
elseif diffZygote > 1e-3
    println("Does not pass ll-test for Zygote Bruno model")
    passTest = false
end
peTabOpt = setUpCostGradHess(peTabModel, Rodas5(), 1e-12, adjSolver=Rodas5(), adjTol=1e-12, sensealg=QuadratureAdjoint(autojacvec=ReverseDiffVJP(false)),
                             sensealgForward = ForwardDiffSensitivity(), solverForward=Rodas4P(autodiff=false),
                             adjSensealg=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(false)))
peTabOptAlt = setUpCostGradHess(peTabModel, Rodas5(), 1e-12, adjSolver=Rodas5(), adjTol=1e-9, sensealg=QuadratureAdjoint(autojacvec=ReverseDiffVJP(false)),
                                sensealgForward = :AutoDiffForward, solverForward=Rodas5())                             
# Check that we can compute gradients for Bruno model 
gradBrunoForward = ForwardDiff.gradient(peTabOpt.evalF, peTabOpt.paramVecTransformed) 
gradBrunoCalc = Calculus.gradient(peTabOpt.evalF, peTabOpt.paramVecTransformed) 
gradImporter = zeros(length(peTabOpt.paramVecTransformed))
gradAdj = zeros(length(peTabOpt.paramVecTransformed))
gradForEq = zeros(length(peTabOpt.paramVecTransformed))
gradForEqAlt = zeros(length(peTabOpt.paramVecTransformed))
peTabOpt.evalGradF(gradImporter, peTabOpt.paramVecTransformed)
gradZygote = Zygote.gradient(peTabOpt.evalFZygote, peTabOpt.paramVecTransformed)[1]
peTabOpt.evalGradFAdjoint(gradAdj, peTabOpt.paramVecTransformed)
peTabOpt.evalGradFForwardEq(gradForEq, peTabOpt.paramVecTransformed)
peTabOptAlt.evalGradFForwardEq(gradForEqAlt, peTabOpt.paramVecTransformed)
if sum((gradImporter - gradBrunoCalc).^2) > 1e-6
    println("Does not pass gradient test for Bruno model")
    passTest = false
end
if sum((gradZygote - gradBrunoCalc).^2) > 1e-6
    println("Does not pass gradient Zygote test for Bruno model")
    passTest = false
end
if sum((gradImporter[1:(end-1)] - gradAdj[1:(end-1)]).^2) > 1e-6
    println("Does not pass gradient Adjoint test for Bruno model")
    passTest = false
end
if sum((gradImporter[1:(end-1)] - gradForEq[1:(end-1)]).^2) > 1e-6
    println("Does not pass gradient ForwardSenseEq test for Bruno model")
    passTest = false
end
if sum((gradImporter[1:(end-1)] - gradForEqAlt[1:(end-1)]).^2) > 1e-6
    println("Does not pass gradient ForwardSenseEq test for Bruno model")
    passTest = false
end


# Crauste model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Crauste_CellSystems2017/"
peTabModel = setUpPeTabModel("model_Crauste_CellSystems2017", dirModel, verbose=false, forceBuildJlFile=true)
peTabOpt = setUpCostGradHess(peTabModel, Rodas4P(), 1e-12)
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed)
costZygote = peTabOpt.evalFZygote(peTabOpt.paramVecTransformed)
diff = abs(cost + (-190.96521897435176))
diffZygote = abs(costZygote + (-190.96521897435176))
if diff > 1e-2
    println("Does not pass ll-test for Crauste model")
    passTest = false
elseif diffZygote > 1e-2
    println("Does not pass ll-test for Zygote Crauste model")
    passTest = false
end


# Elowitz model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Elowitz_Nature2000/"
peTabModel = setUpPeTabModel("model_Elowitz_Nature2000", dirModel, verbose=false, forceBuildJlFile=true)
peTabOpt = setUpCostGradHess(peTabModel, Rodas4P(), 1e-12)
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed)
costZygote = peTabOpt.evalFZygote(peTabOpt.paramVecTransformed)
diff = abs(cost + (63.20279991419332))
diffZygote = abs(costZygote + (63.20279991419332))
if diff > 1e-2
    println("Does not pass ll-test for Elowitz model")
    passTest = false
elseif diffZygote > 1e-2
    println("Does not pass ll-test for Zygote Elowitz model")
    passTest = false
end


# Fiedler model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Fiedler_BMC2016/"
peTabModel = setUpPeTabModel("model_Fiedler_BMC2016", dirModel, verbose=false, forceBuildJlFile=true)
peTabOpt = setUpCostGradHess(peTabModel, Rodas4P(), 1e-9)
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed)
costZygote = peTabOpt.evalFZygote(peTabOpt.paramVecTransformed)
diff = abs(cost + (58.58390161681))
diffZygote = abs(costZygote + (58.58390161681))
if diff > 1e-3
    println("Does not pass ll-test for Fiedler model")
    passTest = false
elseif diffZygote > 1e-3
    println("Does not pass ll-test for Zygote Fiedler model")
    passTest = false
end


# Fujita model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Fujita_SciSignal2010/"
peTabModel = setUpPeTabModel("model_Fujita_SciSignal2010", dirModel, verbose=false, forceBuildJlFile=true)
peTabOpt = setUpCostGradHess(peTabModel, Rodas4P(), 1e-12)
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed)
costZygote = peTabOpt.evalFZygote(peTabOpt.paramVecTransformed)
diff = abs(cost + (53.08377736998929))
diffZygote = abs(costZygote + (53.08377736998929))
if diff > 1e-3
    println("Does not pass ll-test for Fujita model")
    passTest = false
elseif diffZygote > 1e-3
    println("Does not pass ll-test for Zygote Fujita model")
    passTest = false
end


# Isensee model - Extremly strange as it works without priors (must double check)
dirModel = pwd() * "/Intermediate/PeTab_models/model_Isensee_JCB2018/"
peTabModel = setUpPeTabModel("model_Isensee_JCB2018", dirModel, verbose=false, forceBuildJlFile=true)
peTabOpt = setUpCostGradHess(peTabModel, Rodas4P(), 1e-12)
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed) + 4.45299970460275
costZygote = peTabOpt.evalFZygote(peTabOpt.paramVecTransformed) + 4.45299970460275
diff = abs(cost + (-3949.375966548649))
diffZygote = abs(costZygote + (-3949.375966548649))
if diff > 1e-2
    println("Does not pass ll-test for Isensee model")
    passTest = false
elseif diffZygote > 1e-2
    println("Does not pass ll-test for Zygote Isensee model")
    passTest = false
end


# Lucarelli - Breaks code by the same reason as above 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Lucarelli_CellSystems2018/"
peTabModel = setUpPeTabModel("model_Lucarelli_CellSystems2018", dirModel, verbose=false, forceBuildJlFile=true)
peTabOpt = setUpCostGradHess(peTabModel, Rodas4P(), 1e-12)
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed)
costZygote = peTabOpt.evalFZygote(peTabOpt.paramVecTransformed)
diff = abs(cost + (-1681.6059879426584))
diffZygote = abs(costZygote + (-1681.6059879426584))
if diff > 1e-3
    println("Does not pass ll-test for Lucarelli model")
    passTest = false
elseif diffZygote > 1e-3
    println("Does not pass ll-test for Zygote Lucarelli model")
    passTest = false
end


# Schwen model
dirModel = pwd() * "/Intermediate/PeTab_models/model_Schwen_PONE2014/"
peTabModel = setUpPeTabModel("model_Schwen_PONE2014", dirModel, verbose=false, forceBuildJlFile=true)
peTabOpt = setUpCostGradHess(peTabModel, Rodas4P(), 1e-12)
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed) + 12.519137073132825 # We should not include prior value here 
costZygote = peTabOpt.evalFZygote(peTabOpt.paramVecTransformed) + 12.519137073132825
diff = abs(cost + (-943.9992988598723))
diffZygote = abs(costZygote + (-943.9992988598723))
if diff > 1e-3
    println("Does not pass ll-test for Schwen model")
    passTest = false
elseif diffZygote > 1e-3
    println("Does not pass ll-test for Zygote Schwen model")
    passTest = false
end
# For Schwen check that lower level adjoint interface handles priors well 
peTabOpt = setUpCostGradHess(peTabModel, Rodas5(), 1e-9, adjSolver=Rodas5(), adjTol=1e-9, sensealg=QuadratureAdjoint(autojacvec=ReverseDiffVJP(false)),
                             sensealgForward = ForwardDiffSensitivity(), solverForward=Rodas4P(autodiff=false),
                             adjSensealg=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(false)))
peTabOptAlt = setUpCostGradHess(peTabModel, Rodas5(), 1e-9, adjSolver=Rodas5(),
                                sensealgForward = :AutoDiffForward, solverForward=Rodas5())                             
gradFor = zeros(length(peTabOpt.paramVecTransformed))
gradAdj = zeros(length(peTabOpt.paramVecTransformed))
gradForEq = zeros(length(peTabOpt.paramVecTransformed))
gradForEqAlt = zeros(length(peTabOpt.paramVecTransformed))
peTabOpt.evalGradFAdjoint(gradAdj, peTabOpt.paramVecTransformed)
peTabOpt.evalGradFForwardEq(gradForEq, peTabOpt.paramVecTransformed)
peTabOpt.evalGradF(gradFor, peTabOpt.paramVecTransformed)
peTabOptAlt.evalGradFForwardEq(gradForEqAlt, peTabOpt.paramVecTransformed)
if sum((gradFor - gradAdj).^2) > 1e-6
    println("Does not pass prior test for adjoint sensitivity")
    passTest = false
end
if sum((gradFor - gradForEq).^2) > 1e-6
    println("Does not pass prior test for adjoint sensitivity")
    passTest = false
end
if sum((gradFor - gradForEqAlt).^2) > 1e-6
    println("Does not pass gradient ForwardSenseEq test for Schwen model")
    passTest = false
end


# Sneyd model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Sneyd_PNAS2002/"
peTabModel = setUpPeTabModel("model_Sneyd_PNAS2002", dirModel, verbose=false, forceBuildJlFile=true)
peTabOpt = setUpCostGradHess(peTabModel, Rodas4P(), 1e-12)
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed)
costZygote = peTabOpt.evalFZygote(peTabOpt.paramVecTransformed)
diff = abs(cost + (319.79177818768756))
diffZygote = abs(costZygote + (319.79177818768756))
if diff > 1e-3
    println("Does not pass ll-test for Lucarelli model")
    passTest = false
elseif diffZygote > 1e-3
    println("Does not pass ll-test for Zygote Lucarelli model")
    passTest = false
end


# Weber model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Weber_BMC2015/"
peTabModel = setUpPeTabModel("model_Weber_BMC2015", dirModel, verbose=false, forceBuildJlFile=true)
peTabOpt = setUpCostGradHess(peTabModel, Rodas4P(), 1e-12)
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed)
costZygote = peTabOpt.evalFZygote(peTabOpt.paramVecTransformed)
diff = abs(cost + (-296.2017922646865))
diffZygote = abs(costZygote + (-296.2017922646865))
if diff > 1e-3
    println("Does not pass ll-test for Lucarelli model")
    passTest = false
elseif diffZygote > 1e-3
    println("Does not pass ll-test for Zygote Lucarelli model")
    passTest = false
end


# Zheng model 
dirModel = pwd() * "/Intermediate/PeTab_models/model_Zheng_PNAS2012/"
peTabModel = setUpPeTabModel("model_Zheng_PNAS2012", dirModel, verbose=false, forceBuildJlFile=true)
peTabOpt = setUpCostGradHess(peTabModel, Rodas4P(), 1e-12)
cost = peTabOpt.evalF(peTabOpt.paramVecTransformed)
costZygote = peTabOpt.evalFZygote(peTabOpt.paramVecTransformed)
diff = abs(cost + (278.33353271001477))
diffZygote = abs(costZygote + (278.33353271001477))
if diff > 1e-3
    println("Does not pass ll-test for Lucarelli model")
    passTest = false
elseif diffZygote > 1e-3
    println("Does not pass ll-test for Zygote Lucarelli model")
end


if passTest
    println("Passes ll-test")
else
    println("Does not pass ll-test!")
end