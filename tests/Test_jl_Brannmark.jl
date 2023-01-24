
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

dirModel = pwd() * "/Intermediate/PeTab_models/model_Brannmark_JL/"
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