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
using Zygote
using Symbolics
using Sundials
using YAML


BLAS.set_num_threads(1)


# Relevant PeTab structs for compuations 
include(joinpath(pwd(), "src", "PeTab_structs.jl"))

# PeTab importer to get cost, grad etc 
include(joinpath(pwd(), "src", "Create_PEtab_model.jl"))

# HyperCube sampling 
include(joinpath(pwd(), "src", "Optimizers", "Lathin_hypercube.jl"))

# For converting to SBML 
include(joinpath(pwd(), "src", "SBML", "SBML_to_ModellingToolkit.jl"))

include(joinpath(@__DIR__, "..", "Common.jl"))


function getPEtabProblem(petabModel::PEtabModel, 
                         gradientMethod::Symbol, 
                         sensealg, 
                         odeSolver, 
                         absTol::Float64,
                         relTol::Float64, 
                         sparseJacobian::Bool, 
                         chunkSize::Union{Int64, Nothing})

    if gradientMethod == :ForwardEquations
        petabProblem = setUpPEtabODEProblem(petabModel, odeSolver, solverAbsTol=absTol, solverRelTol=relTol, sensealgForwardEquations = sensealg, odeSolverForwardEquations=odeSolver, sparseJacobian=sparseJacobian, chunkSize=chunkSize)
        computeGradient = petabProblem.computeGradientForwardEquations
        return petabProblem, computeGradient
    elseif gradientMethod == :Zygote
        petabProblem = setUpPEtabODEProblem(petabModel, odeSolver, solverAbsTol=absTol, solverRelTol=relTol, sensealgZygote=sensealg, sparseJacobian=sparseJacobian)
        computeGradient = petabProblem.computeGradientZygote
        return petabProblem, computeGradient
    elseif gradientMethod == :Adjoint
        petabProblem = setUpPEtabODEProblem(petabModel, odeSolver, sensealgAdjoint=sensealg, odeSolverAdjoint=odeSolver, solverAdjointAbsTol=absTol, solverAdjointRelTol=relTol, sparseJacobian=sparseJacobian)
        computeGradient = petabProblem.computeGradientAdjoint
        return petabProblem, computeGradient
    elseif gradientMethod == :ForwardDiff
        petabProblem = setUpPEtabODEProblem(petabModel, odeSolver, solverAbsTol=absTol, solverRelTol=relTol, sparseJacobian=sparseJacobian, chunkSize=chunkSize)
        computeGradient = petabProblem.computeGradientAutoDiff
        return petabProblem, computeGradient
    end
end


function benchmarkCostGrad(petabModel::PEtabModel, 
                           gradientInfo, 
                           solversCheck, 
                           pathFileSave, 
                           absTol::Float64, 
                           relTol::Float64; 
                           checkCost::Bool=false, 
                           checkGradient::Bool=false, 
                           checkHessian::Bool=false, 
                           sparseJacobian::Bool=false, 
                           nParamFixed=nothing, 
                           nRepeat=5, 
                           chunkSize=nothing, 
                           iParameter=0,
                           _θ_est=nothing)

    println("Running model ", petabModel.modelName)
    if isnothing(chunkSize)
        chunkSizeWrite = "Default"
    else
        chunkSizeWrite = chunkSize
    end

    # Get parameter vector 
    if isnothing(_θ_est)
        θ_est = getNominalθ(petabModel)
    else
        θ_est = _θ_est
    end
    
    for i in eachindex(solversCheck)

        odeSolver = solversCheck[i][1]
        odeSolverName = solversCheck[i][2]
        runTime = Vector{Float64}(undef, nRepeat)

        if checkGradient == true
            whatCompute = "Gradient"
            gradientMethod, sensealg, methodInfo = gradientInfo
            petabProblem, computeGradient = getPEtabProblem(petabModel, gradientMethod, sensealg, odeSolver, absTol, relTol, sparseJacobian, chunkSize)

            # Use nominal parameter vector 
            println("Precompiling the code")
            gradient = zeros(length(θ_est))
            # Zygote have problems with Steady-state models 
            if petabModel.modelName ∈ ["model_Isensee_JCB2018", "model_Brannmark_JBC2010", "model_Weber_BMC2015"] && gradientMethod == :Zyogte
                return
            end
            local canEval = true
            try 
                computeGradient(gradient, θ_est)
            catch 
                canEval = false
            end
            if all(gradient .== 1e8) || canEval == false
                runTime .= Inf
            else
                for j in 1:nRepeat
                    bGrad = @elapsed computeGradient(gradient, θ_est)
                    runTime[j] = bGrad 
                end
            end

        elseif checkCost == true
            
            whatCompute = "Cost"
            methodInfo = "Standard"
            petabProblem = setUpPEtabODEProblem(petabModel, odeSolver, solverAbsTol=absTol, solverRelTol=relTol, sparseJacobian=sparseJacobian)
            
            println("Precompiling the code")
            cost = petabProblem.computeCost(θ_est)
            if !isinf(cost)
                for j in 1:nRepeat
                    bCost = @elapsed cost = petabProblem.computeCost(θ_est) 
                    runTime[j] = bCost
                end
            else
                runTime .= Inf
            end

        elseif checkHessian == true
            
            whatCompute = "Hessian"
            methodInfo = "ForwardDiff"
            petabProblem = setUpPEtabODEProblem(petabModel, odeSolver, absTol=absTol, relTol=relTol, sparseJacobian=sparseJacobian)
            hessian = zeros(length(θ_est), length(θ_est))
            cost = petabProblem.computeCost(θ_est)
            if !isinf(cost)
                for j in 1:nIter
                    bHess = @elapsed petabProblem.computeHessian(hessian, θ_est)
                    runTime[j] = bHess
                end
            else
                runTime .= Inf
            end
        end

        writeParamFixed = isnothing(nParamFixed) ? 0 : nParamFixed
        dataSave = DataFrame(Time = runTime, 
                             What_calc=whatCompute,
                             Method_info=methodInfo,
                             Model = petabModel.modelName, 
                             absTol = absTol, 
                             relTol = relTol,
                             N_param_fixed=writeParamFixed,
                             I_parameter=iParameter,
                             chunk_size = chunkSizeWrite,
                             solver = odeSolverName)

        if isfile(pathFileSave)
            CSV.write(pathFileSave, dataSave, append = true)
        else
            CSV.write(pathFileSave, dataSave)
        end

        GC.gc()
    end
end


if ARGS[1] == "No_pre_eq_models"

    dirSave = joinpath(@__DIR__, "..", "..", "Intermediate", "Benchmarks", "Cost_grad_hess")
    pathSave = joinpath(dirSave, "No_preequlibrium.csv")
    if !isdir(dirSave)
        mkpath(dirSave)
    end

    modelList = ["model_Boehm_JProteomeRes2014", "model_Bachmann_MSB2011", "model_Beer_MolBioSystems2014", 
                 "model_Bruno_JExpBot2016", "model_Crauste_CellSystems2017", 
                 "model_Elowitz_Nature2000", "model_Fiedler_BMC2016", "model_Fujita_SciSignal2010", 
                 "model_Lucarelli_CellSystems2018", "model_Sneyd_PNAS2002"]                  

    solversCheck = [[Rodas5(), "Rodas5"], 
                    [Rodas5P(), "Rodas5P"], 
                    [QNDF(), "QNDF"]]
    sensealgsCheck = [[:ForwardDiff, nothing, "ForwardDiff"], 
                      [:ForwardEquations, :AutoDiffForward, "ForEq_AutoDiff"],
                      [:Zygote, ForwardDiffSensitivity(), "Zygote_ForwardDiffSensitivity"], 
                      [:Adjoint, InterpolatingAdjoint(autojacvec=ReverseDiffVJP()), "Adj_InterpolatingAdjoint(autojacvec=ReverseDiffVJP())"], 
                      [:Adjoint, QuadratureAdjoint(autojacvec=ReverseDiffVJP()), "Adj_QuadratureAdjoint(autojacvec=ReverseDiffVJP())"], 
                      [:Adjoint, QuadratureAdjoint(autodiff=false, autojacvec=false), "Adj_QuadratureAdjoint(autodiff=false, autojacvec=false)"]]
                    
    for i in eachindex(modelList)
        
        dirModel = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", modelList[i])
        pathYML = getPathYmlFile(dirModel)
        petabModel = readPEtabModel(pathYML)
        absTol, relTol = 1e-8, 1e-8

        # Check cost 
        benchmarkCostGrad(petabModel, nothing, solversCheck, pathSave, absTol, relTol, checkCost=true, nRepeat=5)

        # Check Gradient 
        for sensealgInfo in sensealgsCheck
            benchmarkCostGrad(petabModel, sensealgInfo, solversCheck, pathSave, absTol, relTol, checkGradient=true, nRepeat=5)
        end

        # For fun Check CVODE_BDF
        benchmarkCostGrad(petabModel, [:ForwardEquations, ForwardSensitivity(), "ForEq_ForwardSensitivity"], 
                          [[CVODE_BDF(), "CVODE_BDF"]], pathSave, absTol, relTol, checkGradient=true, nRepeat=5)
    end
end


if ARGS[1] == "Chen_model"

    dirSave = joinpath(@__DIR__, "..", "..", "Intermediate", "Benchmarks", "Cost_grad_hess")
    pathSave = joinpath(dirSave, "Chen_model.csv")
    if !isdir(dirSave)
        mkpath(dirSave)
    end

    modelList = ["model_Chen_MSB2009"]

    solversCheckCost = [[QNDF(), "QNDF"], 
                        [Rodas5(), "Rodas5"], 
                        [KenCarp4(), "KenCarp4"], 
                        [CVODE_BDF(), "CVODE_BDF"]]
    solversCheckCostS = [[QNDF(), "QNDFS"], 
                         [Rodas5(), "Rodas5S"], 
                         [KenCarp4(), "KenCarp4S"], 
                         [CVODE_BDF(linear_solver=:KLU), "CVODE_BDF_KLU"]]                        

    solversCheckGrad = [[QNDF(), "QNDF"], 
                        [KenCarp4(), "KenCarp4"]]
    solversCheckGradS = [[QNDF(), "QNDFS"], 
                         [KenCarp4(), "KenCarp4S"]]                                                 

    sensealgsTry = [[:Adjoint, InterpolatingAdjoint(autojacvec=ReverseDiffVJP()), "Adj_InterpolatingAdjoint(autojacvec=ReverseDiffVJP())"], 
                    [:Adjoint, QuadratureAdjoint(autojacvec=ReverseDiffVJP()), "Adj_QuadratureAdjoint(autojacvec=ReverseDiffVJP())"]]                          

    dirModel = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", modelList[1])
    pathYML = getPathYmlFile(dirModel)
    petabModel = readPEtabModel(pathYML)

    absTol, relTol = 1e-8, 1e-8
    benchmarkCostGrad(petabModel, nothing, solversCheckCost, pathSave, absTol, relTol, checkCost=true, nRepeat=5)
    benchmarkCostGrad(petabModel, nothing, solversCheckCostS, pathSave, absTol, relTol, checkCost=true, nRepeat=5, sparseJacobian=true)                          

    for sensealgInfo in sensealgsTry                          
        benchmarkCostGrad(petabModel, sensealgInfo, solversCheckGradS, pathSave, absTol, relTol, checkGradient=true, nRepeat=5, sparseJacobian=true)                          
        benchmarkCostGrad(petabModel, sensealgInfo, solversCheckGrad, pathSave, absTol, relTol, checkGradient=true, nRepeat=5)
    end
end


if ARGS[1] == "Fix_parameters"

    if length(ARGS) == 1
        modelList = ["model_Bachmann_MSB2011", "model_Lucarelli_CellSystems2018"]
    elseif ARGS[2] == "Bachman"
        modelList = ["model_Bachmann_MSB2011"]
    elseif ARGS[2] == "Lucarelli"
        modelList = ["model_Lucarelli_CellSystems2018"]
    else
        println("ARGS must either be of length 2, or the second argument must be Bachman or Lucarelli not ", ARGS[2])
        exit(1)
    end

    dirSave = joinpath(@__DIR__, "..", "..", "Intermediate", "Benchmarks", "Cost_grad_hess")
    pathSave = joinpath(dirSave, "Fix_parameters.csv")
    if !isdir(dirSave)
        mkpath(dirSave)
    end

    Random.seed!(123)
    solversCheck = [[QNDF(), "QNDF"]]
    sensealgInfo = [[:ForwardDiff, nothing, "ForwardDiff"]]

    for i in eachindex(modelList)

        dirModel = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", modelList[i])
        pathYML = getPathYmlFile(dirModel)
        petabModel = readPEtabModel(pathYML)
        absTol, relTol = 1e-8, 1e-8
        nDynamicParameters = length(getNominalODEValues(petabModel))

        for nParamFix in 1:(nDynamicParameters-1)
            for j in 1:10
                petabModelFewerParameters = getPEtabModelNparamFixed(petabModel, nParamFix)

                benchmarkCostGrad(petabModelFewerParameters, sensealgInfo[1], solversCheck, pathSave, absTol, relTol, 
                                  checkGradient=true, nParamFixed=nParamFix, nRepeat=10)
                benchmarkCostGrad(petabModelFewerParameters, sensealgInfo[1], solversCheck, pathSave, absTol, relTol, 
                                  checkGradient=true, nParamFixed=nParamFix, nRepeat=10, chunkSize=1)      
                                  
                if isdir(petabModelFewerParameters.dirModel) 
                    rm(petabModelFewerParameters.dirModel, recursive=true)
                end
            end
        end

    end
end


if ARGS[1] == "Test_chunks_random_p"

    if length(ARGS) == 1
        modelList = ["model_Bachmann_MSB2011", "model_Lucarelli_CellSystems2018"]
    elseif ARGS[2] == "Bachman"
        modelList = ["model_Bachmann_MSB2011"]
    elseif ARGS[2] == "Lucarelli"
        modelList = ["model_Lucarelli_CellSystems2018"]
    else
        println("ARGS must either be of length 2, or the second argument must be Bachman or Lucarelli not ", ARGS[2])
        exit(1)
    end

    dirSave = joinpath(@__DIR__, "..", "..", "Intermediate", "Benchmarks", "Cost_grad_hess")
    pathSave = joinpath(dirSave, "Test_chunks.csv")
    if !isdir(dirSave)
        mkpath(dirSave)
    end

    Random.seed!(123)
    solversCheck = [[QNDF(), "QNDF"]]
    sensealgInfo = [[:ForwardDiff, nothing, "ForwardDiff"]]

    for i in eachindex(modelList)

        dirModel = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", modelList[i])
        pathYML = getPathYmlFile(dirModel)
        petabModel = readPEtabModel(pathYML)
        absTol, relTol = 1e-8, 1e-8
        nDynamicParameters = length(getNominalODEValues(petabModel))
        chunkList = 1:nDynamicParameters
        for j in 1:30
            θ_est = getRandomModelParameters(petabModel, Rodas4P(), j, odeSolvers=false)
            for k in eachindex(chunkList)
                benchmarkCostGrad(petabModel, sensealgInfo[1], solversCheck, pathSave, absTol, relTol, 
                                checkGradient=true, nRepeat=10, chunkSize=chunkList[k], iParameter=j, _θ_est=θ_est)
            end
        end
    end
end
