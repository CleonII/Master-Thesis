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


function benchmarkCostGrad(peTabModel, modelName::String, solversCheck, pathFileSave, tol; nIter=10, checkHess::Bool=false)

    println("Running model $modelName")

    for i in eachindex(solversCheck)

        solverUse = solversCheck[i][1]
        solverStr = solversCheck[i][2]

        peTabOpt = setUpCostGradHess(peTabModel, solverUse, tol, sensealg = ForwardDiffSensitivity())

        # Use nominal parameter vector 
        paramVec = peTabOpt.paramVecTransformed

        runTimeGrad = Array{Float64, 1}(undef, nIter)
        runTimeGradZygote = Array{Float64, 1}(undef, nIter)
        runTimeCost = Array{Float64, 1}(undef, nIter)
        runTimeHess = Array{Float64, 1}(undef, nIter)

        # Run everything a first time to get pre-compilation out of the way 
        println("Precompiling the code")
        grad = zeros(length(paramVec))
        hess = zeros(length(paramVec), length(paramVec))
        peTabOpt.evalF(paramVec)
        peTabOpt.evalGradF(grad, paramVec)
        # Zygote have problems with SensitivityAdjoint from time-to-time
        if modelName ∉ ["model_Isensee_JCB2018", "model_Brannmark_JBC2010", "model_Weber_BMC2015"]
            peTabOpt.evalGradFZygote(grad, paramVec)
        end
        # For these models the Hessian approximation takes too long time (qudratic complexity)
        if checkHess == true
            peTabOpt.evalHessApprox(hess, paramVec)
        end
        println("Done")

        for j in 1:nIter
            bCost = @elapsed cost = peTabOpt.evalF(paramVec) 
            runTimeCost[j] = bCost
        end

        for j in 1:nIter
            bGrad = @elapsed peTabOpt.evalGradF(grad, paramVec) 
            runTimeGrad[j] = bGrad
        end

        for j in 1:nIter
            bGrad = @elapsed peTabOpt.evalGradFZygote(grad, paramVec) 
            runTimeGradZygote[j] = bGrad
        end

        for j in 1:nIter
            # For these models the Hessian approximation takes too long time (qudratic complexity)
            if checkHess == true
                bHess = @elapsed peTabOpt.evalHessApprox(hess, paramVec)
            else
                bHess = Inf
            end
            runTimeHess[j] = bHess
        end

        dataSave = DataFrame(T_cost = runTimeCost, 
                             T_grad = runTimeGrad, 
                             T_grad_zygote = runTimeGradZygote, 
                             T_hess = runTimeHess, 
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

dirSave = pwd() * "/Intermediate/Benchmarks/Cost_grad_hess/"
pathSave = dirSave * "Cost_grad_hess.csv"
if !isdir(dirSave)
    mkpath(dirSave)
end

modelListTry = ["model_Boehm_JProteomeRes2014", "model_Bachmann_MSB2011", "model_Beer_MolBioSystems2014", 
                "model_Bruno_JExpBot2016", "model_Crauste_CellSystems2017", 
                "model_Elowitz_Nature2000", "model_Fiedler_BMC2016", "model_Fujita_SciSignal2010", 
                "model_Lucarelli_CellSystems2018", "model_Sneyd_PNAS2002"]

solversCheck = [[Rodas5(), "Rodas5"], [Rodas5P(), "Rodas5P"], [QNDF(), "QNDF"]]

for i in eachindex(modelListTry)
    modelName = modelListTry[i]
    dirModel = pwd() * "/Intermediate/PeTab_models/" * modelName * "/"
    peTabModel = setUpPeTabModel(modelName, dirModel)
    # Where we need higher abs- and reltol to solve the ODE 
    if modelName ∈ ["model_Boehm_JProteomeRes2014", "model_Beer_MolBioSystems2014", "model_Crauste_CellSystems2017"]
        tol = 1e-8
    else
        tol = 1e-6
    end
    benchmarkCostGrad(peTabModel, modelName, solversCheck, pathSave, tol, checkHess=false)
end
