#= 
    Check the accruacy of the Julia importer by checking the log-likelihood value against known values for several 
    models containing ifelse-statements and one without.
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


# Beer model 
dirModel = pwd() * "/tests/JuliaImport/model_Beer_MolBioSystems2014/"
peTabModel = setUpPeTabModel("model_Beer_MolBioSystems2014", dirModel, verbose=false, jlFile=true)
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
dirModel = pwd() * "/tests/JuliaImport/model_Boehm_JProteomeRes2014/"
peTabModel = setUpPeTabModel("model_Boehm_JProteomeRes2014", dirModel, verbose=false, jlFile=true)
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
dirModel = pwd() * "/tests/JuliaImport/model_Brannmark_JBC2010/"
peTabModel = setUpPeTabModel("model_Brannmark_JBC2010", dirModel, verbose=false, jlFile=true)
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


# Fujita model 
dirModel = pwd() * "/tests/JuliaImport/model_Fujita_SciSignal2010/"
peTabModel = setUpPeTabModel("model_Fujita_SciSignal2010", dirModel, verbose=false, jlFile=true)
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
dirModel = pwd() * "/tests/JuliaImport/model_Isensee_JCB2018/"
peTabModel = setUpPeTabModel("model_Isensee_JCB2018", dirModel, verbose=false, jlFile=true)
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


# Weber model 
dirModel = pwd() * "/tests/JuliaImport/model_Weber_BMC2015/"
peTabModel = setUpPeTabModel("model_Weber_BMC2015", dirModel, verbose=false, jlFile=true)
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



if passTest
    println("Passes ll-test")
else
    println("Does not pass ll-test!")
end