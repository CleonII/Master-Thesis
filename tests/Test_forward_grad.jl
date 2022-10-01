#=
    In forward gradient we compute an unbiased gradient estimate g via 
    g = (∇f·v)v, where · is the dot product and v ~ MV(0, I). This script tests 
    i) that g is unbiased
    ii) When running optimization foward-gradient is combined with ADAM, hence 
        this script tests for the Rosenbrock function in 10-dimensions that 
        convergece is reached.
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
using Ipopt
using Optim


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

# Optimizers 
include(joinpath(pwd(), "src", "Optimizers", "Set_up_Ipopt.jl"))
include(joinpath(pwd(), "src", "Optimizers", "Set_up_optim.jl"))
include(joinpath(pwd(), "src", "Optimizers", "Set_up_forward_gradient.jl"))


"""
    testGradForwardGrad(f::Function, p)::Bool

    For forward gradient test if the gradient approximaiton is unbiased by 
    computing a gradient from 100,000 Monte-Carlo samples and then compare 
    it against the gradient obtained from ForwardDiff for function f, 
    and parameter vector p using the squared difference of the normalized 
    gradients.
"""
function testGradForwardGrad(f::Function, p; seed=123, verbose::Bool=false)::Bool

    Random.seed!(seed)

    calcApproxGrad = (pArg) -> calcUnbiasedGrad(pArg, f)
    gradTrue = ForwardDiff.gradient(f, p)
    
    nSamples = 1000000
    gradEst = zeros(length(p))
    for i in 1:nSamples
        gradEst += calcApproxGrad(p)
    end
    gradEst ./= nSamples

    sqDiff = sum((gradEst / norm(gradEst) - gradTrue / norm(gradTrue)).^2)
    if verbose == true
        @printf("sqDiff = %.3e\n", sqDiff)
    end

    if sqDiff < 1e-4
        return true
    else
        @printf("sqDiff = %.3e\n", sqDiff)
        return false
    end
end


# Functions to test again 
function f1(x::AbstractVector)::Real
    return sum(x.^2)
end
function f2(x::AbstractVector)::Real
    return sum(cos.(sin.(x)))
end
function f3(x::AbstractVector)::Real
    return sum(exp.(x.^2))
end
function rosenbrock(x::AbstractVector)
    return( sum( 100*( x[2:end] .- x[1:end-1].^2 ).^2 + ( x[1:end-1] .- 1 ).^2 ) )
end


"""
    testAdamRosenbrock(nParam::Integer)::Bool

    Test the ADAM implementation with forward gradient for the Rosenbrock function 
    with nParam, using 100,000 samples and the zero-vector as start guess. Tests 
    both changing and constant step-length.
"""
function testAdamRosenbrock(nParam::Integer; verbose=false)::Bool
    
    nSamples = 100000
    lB, uB = [-Inf for i in 1:nParam], [Inf for i in 1:nParam]

    # With constant step length 
    adamOpt = createFowardGradientProb(rosenbrock, 1e-4, nParam, 100000, lB, uB)
    minCost, minParams, costVal, paramVal = runAdam(zeros(nParam), adamOpt, verbose=false)         
    if verbose == true
        @printf("minCost = %.3e\n", minCost)
    end
    
    if minCost > 1e-4
        @printf("Constant step-length minCost = %.3e\n", minCost)
        return false
    end

    # With varying step length 
    nSamples = 100000
    stepLengths = log10.(LinRange(exp10(1e-3),exp10(1e-5), nSamples))
    adamOpt = createFowardGradientProb(rosenbrock, stepLengths, nParam, nSamples, lB, uB)
    minCost, minParams, costVal, paramVal = runAdam(zeros(nParam), adamOpt, verbose=false)                                  
    if verbose == true
        @printf("minCost = %.3e\n", minCost)
    end

    if minCost > 1e-4
        @printf("Adaptive step-length minCost = %.3e\n", minCost)
        return false
    end

    return true
end


# Test if gradient estimate is unbiased 
Random.seed!(123)
pTry = [rand(2), rand(4), rand(8)]

passTest = [testGradForwardGrad(f1, pTry[i], verbose=false) for i in eachindex(pTry)]
if all(passTest) == true
    @printf("Passed forward gradient test for f1\n")
else
    @printf("Did not pass forward gradient test for f1\n")
end

passTest = [testGradForwardGrad(f2, pTry[i], verbose=false) for i in eachindex(pTry)]
if all(passTest) == true
    @printf("Passed forward gradient test for f2\n")
else
    @printf("Did not pass forward gradient test for f2\n")
end

passTest = [testGradForwardGrad(f3, pTry[i], verbose=false) for i in eachindex(pTry)]
if all(passTest) == true
    @printf("Passed forward gradient test for f3\n")
else
    @printf("Did not pass forward gradient test for f3\n")
end

passTest = [testGradForwardGrad(rosenbrock, pTry[i], verbose=false) for i in eachindex(pTry)]
if all(passTest) == true
    @printf("Passed forward gradient test for rosenbrock\n")
else
    @printf("Did not pass forward gradient test for rosenbrock\n")
end

passTest = testAdamRosenbrock(10, verbose=false)
if passTest == true
    @printf("Passed forward gradient test with ADAM for rosenbrock\n")
else
    @printf("Did not pass forward gradient test with ADAM for rosenbrock\n")
end

