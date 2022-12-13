#= 
    Using the Bachman and Fideler model this file tests that our parallellisation scales well. 

    Briefly, using Distributed.jl we parallellise over different experimental conditions.
    The Bachman models has 36 different experimental conditions (we solve the ODE 36 times). 
    If our parallellisation works well we should se substantial speedup for the cost, and 
    especially for the gradient, close to the ideal reduction of factor 2 when using 
    two processes. 
    Similarly, the Brännmark model has 8 different experimental conditions, and we can further 
    compute the hessian and an associated hessian approximation.
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
using BenchmarkTools
using Zygote
using Distributed


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


function testRunBenchmark(peTabModel::PeTabModel, solver, tol::Float64, nProcs::Integer; verbose=true, computeH::Bool=false)::Bool

    peTabOptSeq = setUpCostGradHess(peTabModel, solver, tol)
    peTabOptPar = setUpCostGradHess(peTabModel, solver, tol, nProcs=nProcs)

    passTest = true

    # Compute everything first for the sake of pre-compilation to get fair benchmark
    grad1, grad2 = zeros(length(peTabOptSeq.paramVecTransformed)), zeros(length(peTabOptSeq.paramVecTransformed))
    costSeq = peTabOptSeq.evalF(peTabOptSeq.paramVecTransformed)
    costPar = peTabOptPar.evalF(peTabOptSeq.paramVecTransformed)
    peTabOptSeq.evalGradF(grad1, peTabOptSeq.paramVecTransformed)
    peTabOptPar.evalGradF(grad2, peTabOptSeq.paramVecTransformed)

    # Check accuracy of parallell computation before moving on
    if abs(costSeq - costPar) > 1e-5
        println("Does not pass test on accurately computing cost.")
        @printf("Cost sequential = %.e3\n", costSeq)
        @printf("Cost with %d processes = %.e3\n",  costPar)
        return false
    elseif verbose
        @printf("abs-diff cost = %.e3\n",  abs(costSeq - costPar))
    end
    if sum((grad1 - grad2).^2) > 1e-5
        println("Does not pass test on accurately computing gradient")
        @printf("sqDiffGrad = %.e3\n", sum((grad1 - grad2).^2))
        return false
    elseif verbose
        @printf("sqDiffGrad = %.e3\n", sum((grad1 - grad2).^2))
    end
    if computeH == true
        dimP = length(peTabOptSeq.paramVecTransformed)
        hess1, hess2 = zeros(dimP, dimP), zeros(dimP, dimP)

        peTabOptSeq.evalHess(hess1, peTabOptSeq.paramVecTransformed)
        peTabOptPar.evalHess(hess2, peTabOptSeq.paramVecTransformed)
        if sum((hess1 - hess2).^2) > 1e-5
            println("Does not pass test on accurately computing hessian")
            @printf("sqDiffHess = %.e3\n", sum((hess1 - hess2).^2))
            return false
        elseif verbose
            @printf("sqDiffHess = %.e3\n", sum((hess1 - hess2).^2))
        end
           
        hess1, hess2 = zeros(dimP, dimP), zeros(dimP, dimP)
        peTabOptSeq.evalHessApprox(hess1, peTabOptSeq.paramVecTransformed)
        peTabOptPar.evalHessApprox(hess2, peTabOptSeq.paramVecTransformed)
        if sum((hess1 - hess2).^2) > 1e-5
            println("Does not pass test on accurately computing hessian")
            @printf("sqDiffHessApprox = %.e3\n", sum((hess1 - hess2).^2))
            return false
        elseif verbose
            @printf("sqDiffHessApprox = %.e3\n", sum((hess1 - hess2).^2))
        end
    end

    b1, b2 = 0.0, 0.0
    for i in 1:10
        b1Tmp = @elapsed peTabOptSeq.evalF(peTabOptSeq.paramVecTransformed)
        sleep(0.1)
        b2Tmp = @elapsed peTabOptPar.evalF(peTabOptSeq.paramVecTransformed)
        b1 += b1Tmp
        b2 += b2Tmp
        sleep(0.1)
    end
    if b1 < b2
        println("Does not pass test. On cost for model ", peTabModel.modelName, " the run time 
                 using $nProcs processes is slower than using a single process")
        @printf("Run time single proces = %.2e s", b1 / 10)
        @printf("Run time %d processes = %.2e s", nProcs, b2 / 10)
        return false
    elseif verbose
        @printf("For cost\n")
        @printf("Run time single proces = %.2e s\n", b1 / 10)
        @printf("Run time %d processes = %.2e s\n", nProcs, b2 / 10)
        @printf("Ratio = %.2f\n\n", b2 / b1)
    end
    
    b1, b2 = 0.0, 0.0
    for i in 1:10
        b1Tmp = @elapsed peTabOptSeq.evalGradF(grad1, peTabOptSeq.paramVecTransformed)
        sleep(0.1)
        b2Tmp = @elapsed peTabOptPar.evalGradF(grad2, peTabOptSeq.paramVecTransformed)
        b1 += b1Tmp
        b2 += b2Tmp
        sleep(0.1)
    end
    if b1 < b2
        println("Does not pass test. On gradient for model ", peTabModel.modelName, " the run time 
                 using $nProcs processes is slower than using a single process")
        @printf("Run time single proces = %.2e s\n", b1 / 10)
        @printf("Run time %d processes = %.2e s\n", nProcs, b2 / 10)
        return false
    elseif verbose
        @printf("For gradient\n")
        @printf("Run time single proces = %.2e s\n", b1 / 10)
        @printf("Run time %d processes = %.2e s\n", nProcs, b2 / 10)
        @printf("Ratio = %.2f\n\n", b2 / b1)
    end
    
    if computeH == true
        b1, b2 = 0.0, 0.0
        for i in 1:5
            b1Tmp = @elapsed peTabOptSeq.evalHessApprox(hess1, peTabOptSeq.paramVecTransformed)
            sleep(0.1)
            b2Tmp = @elapsed peTabOptSeq.evalHessApprox(hess2, peTabOptSeq.paramVecTransformed)
            b1 += b1Tmp
            b2 += b2Tmp
            sleep(0.1)
        end
        if b1 < b2
            println("Does not pass test. On hessian for model ", peTabModel.modelName, " the run time 
                    using $nProcs processes is slower than using a single process")
            @printf("Run time single proces = %.2e s\n", b1 / 10)
            @printf("Run time %d processes = %.2e s\n", nProcs, b2 / 10)
            return false
        elseif verbose
            @printf("For hessian block approximation\n")
            @printf("Run time single proces = %.2e s\n", b1 / 10)
            @printf("Run time %d processes = %.2e s\n", nProcs, b2 / 10)
            @printf("Ratio = %.2f\n\n", b2 / b1)
        end
    end

    return passTest
end


peTabModel = setUpPeTabModel("Bachmann_MSB2011", pwd() * "/tests/Bachman/")
removeAllProcs()
addprocs(1, exeflags="--project=.")
passTest = testRunBenchmark(peTabModel, Rodas5(), 1e-9, 2)
if passTest == false
    println("Does not pass parallellisation test for Bachman model")
end


dirModel = pwd() * "/Intermediate/PeTab_models/model_Brannmark_JBC2010/"
peTabModel = setUpPeTabModel("model_Brannmark_JBC2010", dirModel)
removeAllProcs()
solver, tol = Rodas5(), 1e-9
addprocs(1, exeflags="--project=.")
passTest = testRunBenchmark(peTabModel, Rodas5(), 1e-9, 2, computeH=true)
if passTest == false
    println("Does not pass parallellisation test for Brännmark model")
end
