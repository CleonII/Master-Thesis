using ModelingToolkit 
using DifferentialEquations
using DataFrames
using CSV 
using ForwardDiff
using ReverseDiff
using StatsBase
using Random
using Printf
using BenchmarkTools
using Zygote
using SciMLSensitivity
using BenchmarkTools
using SparseDiffTools
using LinearAlgebra
using Sundials


# Relevant PeTab structs for computations 
include(joinpath(pwd(), "src", "PeTab_structs.jl"))

# Functions for solving ODE system 
include(joinpath(pwd(), "src", "Solve_ODE_model", "Solve_ode_model.jl"))

# PeTab importer to get cost, grad etc 
include(joinpath(pwd(), "src", "PeTab_importer", "Create_cost_grad_hessian.jl"))

# HyperCube sampling 
include(joinpath(pwd(), "src", "Optimizers", "Lathin_hypercube.jl"))

# For converting to SBML 
include(joinpath(pwd(), "src", "SBML", "SBML_to_ModellingToolkit.jl"))


function getPEtabOpt(gradMethod, sensealg, solverUse, tol)

    if gradMethod == :ForwardSenseEq
        peTabOpt = setUpCostGradHess(peTabModel, solverUse, tol, sensealgForward = sensealg, solverForward=solverUse)
        evalGradF = peTabOpt.evalGradFForwardEq
        return peTabOpt, evalGradF
    elseif gradMethod == :Zygote
        peTabOpt = setUpCostGradHess(peTabModel, solverUse, tol, sensealg = sensealg)
        evalGradF = peTabOpt.evalGradFZygote
        return peTabOpt, evalGradF
    elseif gradMethod == :Adjoint
        peTabOpt = setUpCostGradHess(peTabModel, solverUse, tol, adjSolver=solverUse, adjSensealg=sensealg, adjTol=tol)
        evalGradF = peTabOpt.evalGradFAdjoint
        return peTabOpt, evalGradF
    elseif gradMethod == :ForwardDiff
        peTabOpt = setUpCostGradHess(peTabModel, solverUse, tol)
        evalGradF = peTabOpt.evalGradF
        return peTabOpt, evalGradF
    end
end


function benchmarkCostGrad(peTabModel, modelName::String, gradInfo, solversCheck, pathFileSave, tol; nIter=10, checkHess::Bool=false, checkCost::Bool=false, checkGrad::Bool=false)

    println("Running model $modelName")
    
    for i in eachindex(solversCheck)

        solverUse = solversCheck[i][1]
        solverStr = solversCheck[i][2]
        runTime = Array{Float64, 1}(undef, nIter)

        if checkGrad == true
            what_calc = "Gradient"
            gradMethod, sensealg, methodInfo = gradInfo
            peTabOpt, evalGradF = getPEtabOpt(gradMethod, sensealg, solverUse, tol)

            # Use nominal parameter vector 
            println("Precompiling the code")
            paramVec = peTabOpt.paramVecTransformed
            grad = zeros(length(paramVec))
            # Zygote have problems with SensitivityAdjoint from time-to-time
            if modelName ∈ ["model_Isensee_JCB2018", "model_Brannmark_JBC2010", "model_Weber_BMC2015"] && gradMethod == :Zyogte
                return
            end
            evalGradF(grad, peTabOpt.paramVecTransformed)
            for j in 1:nIter
                bGrad = @elapsed evalGradF(grad, paramVec)
                runTime[j] = bGrad
            end

        elseif checkCost == true
            what_calc = "Cost"
            methodInfo = "Standard"
            peTabOpt = setUpCostGradHess(peTabModel, solverUse, tol)
            paramVec = peTabOpt.paramVecTransformed
            println("Precompiling the code")
            peTabOpt.evalF(paramVec)
            for j in 1:nIter
                bCost = @elapsed cost = peTabOpt.evalF(paramVec) 
                runTime[j] = bCost
            end

        elseif checkHess == true
            what_calc = "Hessian"
            methodInfo = "ForwardDiff"
            peTabOpt = setUpCostGradHess(peTabModel, solverUse, tol, sensealg = sensealg)
            paramVec = peTabOpt.paramVecTransformed
            hess = zeros(length(paramVec), length(paramVec))
            for j in 1:nIter
                bHess = @elapsed peTabOpt.evalHessApprox(hess, paramVec)
                runTime[j] = bHess
            end
        end

        dataSave = DataFrame(Time = runTime, 
                             What_calc=what_calc,
                             Method_info=methodInfo,
                             model = modelName, 
                             tol = tol, 
                             solver = solverStr)

        if isfile(pathFileSave)
            CSV.write(pathFileSave, dataSave, append = true)
        else
            CSV.write(pathFileSave, dataSave)
        end

        GC.gc()
    end
end


if ARGS[1] == "No_pre_eq_models"

    dirSave = pwd() * "/Intermediate/Benchmarks/Cost_grad_hess/"
    pathSave = dirSave * "Cost_grad.csv"
    if !isdir(dirSave)
        mkpath(dirSave)
    end

    modelListTry = ["model_Boehm_JProteomeRes2014", "model_Bachmann_MSB2011", "model_Beer_MolBioSystems2014", 
                    "model_Bruno_JExpBot2016", "model_Crauste_CellSystems2017", 
                    "model_Elowitz_Nature2000", "model_Fiedler_BMC2016", "model_Fujita_SciSignal2010", 
                    "model_Lucarelli_CellSystems2018", "model_Sneyd_PNAS2002"]

    solversCheck = [[Rodas5(), "Rodas5"], 
                    [Rodas5P(), "Rodas5P"], 
                    [QNDF(), "QNDF"]]
    sensealgInfoTot = [[:ForwardDiff, nothing, "ForwardDiff"], 
                       [:Zygote, ForwardDiffSensitivity(), "Zygote_ForwardDiffSensitivity"], 
                       [:Adjoint, InterpolatingAdjoint(autojacvec=ReverseDiffVJP()), "Adj_InterpolatingAdjoint(autojacvec=ReverseDiffVJP())"], 
                       [:Adjoint, QuadratureAdjoint(autojacvec=ReverseDiffVJP()), "Adj_QuadratureAdjoint(autojacvec=ReverseDiffVJP())"], 
                       [:Adjoint, QuadratureAdjoint(autodiff=false, autojacvec=false), "Adj_QuadratureAdjoint(autodiff=false, autojacvec=false)"]]

    solverCheckForEq = [[Rodas5(autodiff=false), "Rodas5"], 
                        [Rodas5P(autodiff=false), "Rodas5P"], 
                        [QNDF(autodiff=false), "QNDF"]]
    sensealgInfoForEq = [[:Adjoint, InterpolatingAdjoint(autodiff=false, autojacvec=false), "Adj_InterpolatingAdjoint(autodiff=false, autojacvec=false)"],
                         [:ForwardSenseEq, ForwardDiffSensitivity(), "ForEq_ForwardDiffSensitivity"], 
                         [:ForwardSenseEq, ForwardSensitivity(), "ForEq_ForwardSensitivity"]]
                    
    for i in eachindex(modelListTry)
        modelName = modelListTry[i]
        dirModel = pwd() * "/Intermediate/PeTab_models/" * modelName * "/"
        peTabModel = setUpPeTabModel(modelName, dirModel)
        tol = 1e-8

        # Check cost 
        benchmarkCostGrad(peTabModel, modelName, nothing, solversCheck, pathSave, tol, checkCost=true)

        # Check Gradient 
        for sensealgInfo in sensealgInfoTot
            benchmarkCostGrad(peTabModel, modelName, sensealgInfo, solversCheck, pathSave, tol, checkGrad=true)
        end

        # Check gradient for sensitivity equations 
        for sensealgInfo in sensealgInfoForEq
            benchmarkCostGrad(peTabModel, modelName, sensealgInfo, solverCheckForEq, pathSave, tol, checkGrad=true)
        end

        # For fun Check CVODE_BDF
        benchmarkCostGrad(peTabModel, modelName, [:ForwardSenseEq, ForwardSensitivity(), "ForEq_ForwardSensitivity"], 
                          [[CVODE_BDF(), "CVODE_BDF"]], pathSave, tol, checkGrad=true)
    end
end


if ARGS[1] == "Large_models"

    dirSave = pwd() * "/Intermediate/Benchmarks/Cost_grad_hess/"
    pathSave = dirSave * "Cost_grad_hess_large.csv"
    if !isdir(dirSave)
        mkpath(dirSave)
    end

    modelListTry = ["model_Chen_MSB2009"]

    solversCheck = [[QNDF(), "QNDF"]]

    for i in eachindex(modelListTry)
        modelName = modelListTry[i]
        dirModel = pwd() * "/Intermediate/PeTab_models/" * modelName * "/"
        peTabModel = setUpPeTabModel(modelName, dirModel)
        # Where we need higher abs- and reltol to solve the ODE 
        tol = 1e-8
        benchmarkCostGrad(peTabModel, modelName, solversCheck, pathSave, tol, checkHess=false, nIter=2, sensealg=QuadratureAdjoint(autojacvec=ReverseDiffVJP()))
    end
end