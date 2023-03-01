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
using BenchmarkTools

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
        petabProblem = setUpPEtabODEProblem(petabModel, odeSolver, sensealgAdjoint=sensealg, sensealgAdjointSS=sensealg, odeSolverAdjoint=odeSolver, solverAdjointAbsTol=absTol, solverAdjointRelTol=relTol, sparseJacobian=sparseJacobian)
        computeGradient = petabProblem.computeGradientAdjoint
        return petabProblem, computeGradient
    elseif gradientMethod == :ForwardDiff
        petabProblem = setUpPEtabODEProblem(petabModel, odeSolver, solverAbsTol=absTol, solverRelTol=relTol, sparseJacobian=sparseJacobian, chunkSize=chunkSize)
        computeGradient = petabProblem.computeGradientAutoDiff
        return petabProblem, computeGradient
    end
end


function getRunTime(computeGradient::F, runTime, nRepeat, gradient, θ_est) where F

    computeGradient(gradient, θ_est)
    bGrad =  @benchmark $computeGradient($gradient, $θ_est) samples=(nRepeat+1) seconds=100000 evals=1
    println("bGrad1 = ", bGrad)
    bGrad =  @benchmark $computeGradient($gradient, $θ_est) samples=(nRepeat+1) seconds=100000 evals=1
    println("bGrad2 = ", bGrad)
    runTime .= bGrad.times[2:end] .* 1e-9

end


function benchmarkCostGrad(petabModel::PEtabModel, 
                           gradientInfo, 
                           odeSolver::SciMLAlgorithm,
                           odeSolverName::String, 
                           pathFileSave::String, 
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
                           _θ_est=nothing,
                           shouldSave::Bool=true)

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
    
    runTime = Vector{Float64}(undef, nRepeat)

    if checkGradient == true
        whatCompute = "Gradient"
        gradientMethod, sensealg, methodInfo = gradientInfo
        petabProblem, computeGradient = getPEtabProblem(petabModel, gradientMethod, sensealg, odeSolver, absTol, relTol, sparseJacobian, chunkSize)

        # Use nominal parameter vector 
        gradient = zeros(length(θ_est))
        # Zygote have problems with Steady-state models 
        if petabModel.modelName ∈ ["model_Isensee_JCB2018", "model_Brannmark_JBC2010", "model_Weber_BMC2015", "model_Fiedler_BMC2016"] && gradientMethod == :Zygote
            println("Will not run Zygote")
            return
        end
        # Enzyme does not currently handle callbacks 
        if petabModel.modelName == "model_Isensee_JCB2018" && (methodInfo == "InterpolatingAdjoint(autojacvec=EnzymeVJP())" || methodInfo == "QuadratureAdjoint(autojacvec=EnzymeVJP())")
            return 
        end

        # Callbacks are not handled well at all by the lower level interface 
        if petabModel.modelName ∈ ["model_Fujita_SciSignal2010", "model_Beer_MolBioSystems2014"] && gradientMethod == :Adjoint
            println("Will not run Zygote")
            return
        end
        println("Precompiling the code for method ", methodInfo)
        local canEval = true
        try 
            computeGradient(gradient, θ_est)
        catch 
            canEval = false
        end
        if all(gradient .== 1e8) || canEval == false
            runTime .= Inf
        else
            GC.gc(), GC.gc(), GC.gc()
            sleep(0.2)
            getRunTime(computeGradient, runTime, nRepeat, gradient, θ_est) 
            println("BGrad = ", runTime)
        end

    elseif checkCost == true
            
        whatCompute = "Cost"
        methodInfo = "Standard"
        petabProblem = setUpPEtabODEProblem(petabModel, odeSolver, solverAbsTol=absTol, solverRelTol=relTol, sparseJacobian=sparseJacobian)
            
        println("Precompiling the code")
        cost = petabProblem.computeCost(θ_est)
        if !isinf(cost)
            bCost = @benchmark $petabProblem.computeCost($θ_est) samples=(nRepeat+1) seconds=1e5 evals=1
            runTime .= bCost.times[2:end] .* 1e-9
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
            bHess = @benchmark $petabProblem.computeHessian($hessian, $θ_est) samples=(nRepeat+1) seconds=1e5 evals=1
            runTime .= bHess.times
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

    if isfile(pathFileSave) && shouldSave == true
        CSV.write(pathFileSave, dataSave, append = true)
    elseif shouldSave == true
        CSV.write(pathFileSave, dataSave)
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
                            

    odeSolvers = [Rodas5(), KenCarp4(), QNDF(), AutoVern7(Rodas5())]
    odeSolversName = ["Rodas5", "KenCarp4", "QNDF", "Vern7(Rodas5)"]                 
    sensealgsCheck = [[:ForwardDiff, nothing, "ForwardDiff"], 
                      [:ForwardEquations, :AutoDiffForward, "ForEq_AutoDiff"],
                      [:Zygote, ForwardDiffSensitivity(), "Zygote_ForwardDiffSensitivity"], 
                      [:Adjoint, InterpolatingAdjoint(autojacvec=ReverseDiffVJP()), "Adj_InterpolatingAdjoint(autojacvec=ReverseDiffVJP())"], 
                      [:Adjoint, QuadratureAdjoint(autojacvec=ReverseDiffVJP()), "Adj_QuadratureAdjoint(autojacvec=ReverseDiffVJP())"], 
                      [:Adjoint, QuadratureAdjoint(autodiff=false, autojacvec=false), "Adj_QuadratureAdjoint(autodiff=false, autojacvec=false)"]]

    absTol, relTol = 1e-8, 1e-8                      
    for i in eachindex(modelList)
        dirModel = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", modelList[i])
        pathYML = getPathYmlFile(dirModel)
        petabModel = readPEtabModel(pathYML)
        for j in eachindex(odeSolvers)
        
            # Check cost 
            benchmarkCostGrad(petabModel, nothing, odeSolvers[j], odeSolversName[j], pathSave, absTol, relTol, checkCost=true, nRepeat=10)

            # Check Gradient 
            for sensealgInfo in sensealgsCheck
                benchmarkCostGrad(petabModel, sensealgInfo, odeSolvers[j], odeSolversName[j], pathSave, absTol, relTol, checkGradient=true, nRepeat=10)
            end
        end

        # For fun Check CVODE_BDF
        benchmarkCostGrad(petabModel, [:ForwardEquations, ForwardSensitivity(), "ForEq_ForwardSensitivity"], 
                          CVODE_BDF(), "CVODE_BDF", pathSave, absTol, relTol, checkGradient=true, nRepeat=10)
    end
end


if ARGS[1] == "Test_adjoint"

    dirSave = joinpath(@__DIR__, "..", "..", "Intermediate", "Benchmarks", "Cost_grad_hess")
    pathSave = joinpath(dirSave, "Test_adjoint.csv")
    if !isdir(dirSave)
        mkpath(dirSave)
    end
    modelTest = ARGS[2]


    modelList = ["model_Boehm_JProteomeRes2014", "model_Bachmann_MSB2011",  "model_Lucarelli_CellSystems2018", "model_Isensee_JCB2018", "Smith_BMCSystBiol2013"]
                            
    odeSolvers = [Rodas5P(), QNDF(), CVODE_BDF()]
    odeSolversName = ["Rodas5P", "QNDF", "CVODE_BDF"]                 
    sensealgsCheck = [[:ForwardDiff, nothing, "ForwardDiff"], 
                      [:Adjoint, InterpolatingAdjoint(autojacvec=ReverseDiffVJP(true)), "InterpolatingAdjoint(autojacvec=ReverseDiffVJP(true))"],
                      [:Adjoint, InterpolatingAdjoint(autojacvec=ReverseDiffVJP(false)), "InterpolatingAdjoint(autojacvec=ReverseDiffVJP(false))"],
                      [:Adjoint, InterpolatingAdjoint(autojacvec=EnzymeVJP()), "InterpolatingAdjoint(autojacvec=EnzymeVJP())"], 
                      [:Adjoint, QuadratureAdjoint(abstol=1e-8, reltol=1e-8, autojacvec=ReverseDiffVJP(true)), "QuadratureAdjoint(autojacvec=ReverseDiffVJP(true))"],
                      [:Adjoint, QuadratureAdjoint(abstol=1e-8, reltol=1e-8, autojacvec=ReverseDiffVJP(false)), "QuadratureAdjoint(autojacvec=ReverseDiffVJP(false))"], 
                      [:Adjoint, QuadratureAdjoint(abstol=1e-8, reltol=1e-8, autojacvec=EnzymeVJP()), "QuadratureAdjoint(autojacvec=EnzymeVJP())"]]

    absTol, relTol = 1e-8, 1e-8                      
    dirModel = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", modelTest)
    pathYML = getPathYmlFile(dirModel)
    petabModel = readPEtabModel(pathYML)
    for j in eachindex(odeSolvers)
    
        # Check Gradient 
        for sensealgInfo in sensealgsCheck
            benchmarkCostGrad(petabModel, sensealgInfo, odeSolvers[j], odeSolversName[j], pathSave, absTol, relTol, checkGradient=true, nRepeat=5)
        end
    end

end


if ARGS[1] == "Test_flags"

    dirSave = joinpath(@__DIR__, "..", "..", "Intermediate", "Benchmarks", "Cost_grad_hess")
    pathSave = joinpath(dirSave, "Test_flags" * ARGS[2] * ".csv")
    if !isdir(dirSave)
        mkpath(dirSave)
    end

    modelList = ["model_Boehm_JProteomeRes2014", "model_Bachmann_MSB2011", "model_Beer_MolBioSystems2014", 
                 "model_Bruno_JExpBot2016", "model_Crauste_CellSystems2017", 
                 "model_Elowitz_Nature2000", "model_Fiedler_BMC2016", "model_Fujita_SciSignal2010", 
                 "model_Lucarelli_CellSystems2018", "model_Sneyd_PNAS2002"]
                            

    odeSolvers = [Rodas5(), KenCarp4(), QNDF(), AutoVern7(Rodas5())]
    odeSolversName = ["Rodas5", "KenCarp4", "QNDF", "Vern7(Rodas5)"]                 
    sensealgsCheck = [[:ForwardDiff, nothing, "ForwardDiff"]]

    absTol, relTol = 1e-8, 1e-8                      
    for i in eachindex(modelList)
        dirModel = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", modelList[i])
        pathYML = getPathYmlFile(dirModel)
        petabModel = readPEtabModel(pathYML)
        for j in eachindex(odeSolvers)
        
            # Check cost 
            benchmarkCostGrad(petabModel, nothing, odeSolvers[j], odeSolversName[j], pathSave, absTol, relTol, checkCost=true, nRepeat=10)

            # Check Gradient 
            for sensealgInfo in sensealgsCheck
                benchmarkCostGrad(petabModel, sensealgInfo, odeSolvers[j], odeSolversName[j], pathSave, absTol, relTol, checkGradient=true, nRepeat=10)
                benchmarkCostGrad(petabModel, sensealgInfo, odeSolvers[j], odeSolversName[j], pathSave, absTol, relTol, checkGradient=true, nRepeat=10, chunkSize=1)
            end
        end
    end
end


if ARGS[1] == "Chen_model"

    dirSave = joinpath(@__DIR__, "..", "..", "Intermediate", "Benchmarks", "Cost_grad_hess")
    pathSave = joinpath(dirSave, "Chen_model.csv")
    if !isdir(dirSave)
        mkpath(dirSave)
    end

    modelList = ["model_Chen_MSB2009"]
    odeSolversCost = [QNDF(), Rodas5(), KenCarp4(), CVODE_BDF()]
    odeSolversCostName = ["QNDF", "Rodas5", "KenCarp4", "CVODE_BDF"]                 
    odeSolversCostS = [QNDF(), Rodas5(), KenCarp4(), CVODE_BDF(linear_solver=:KLU)]
    odeSolversCostNameS = ["QNDF_S", "Rodas5_S", "KenCarp4_S", "CVODE_BDF_S"]                 

    odeSolversGradient = [QNDF(), KenCarp4()]
    odeSolversGradientName = ["QNDF", "KenCarp4"]    
    odeSolversGradientS = [QNDF(), KenCarp4()]
    odeSolversGradientSName = ["QNDF_S", "KenCarp4_S"]    

    sensealgsTry = [[:Adjoint, InterpolatingAdjoint(autojacvec=ReverseDiffVJP()), "Adj_InterpolatingAdjoint(autojacvec=ReverseDiffVJP())"], 
                    [:Adjoint, QuadratureAdjoint(autojacvec=ReverseDiffVJP()), "Adj_QuadratureAdjoint(autojacvec=ReverseDiffVJP())"]]                          

    dirModel = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", modelList[1])
    pathYML = getPathYmlFile(dirModel)
    petabModel = readPEtabModel(pathYML)

    absTol, relTol = 1e-8, 1e-8
    for i in eachindex(odeSolversCost)
        benchmarkCostGrad(petabModel, nothing, odeSolversCost[i], odeSolversCostName[i], pathSave, absTol, relTol, checkCost=true, nRepeat=5)
        benchmarkCostGrad(petabModel, nothing, odeSolversCostS[i], odeSolversCostNameS[i], pathSave, absTol, relTol, checkCost=true, nRepeat=5, sparseJacobian=true)                          
    end

    for sensealgInfo in sensealgsTry  
        for i in eachindex(odeSolversGradient)                        
            benchmarkCostGrad(petabModel, sensealgInfo, odeSolversGradient[i], odeSolversGradientName[i], pathSave, absTol, relTol, checkGradient=true, nRepeat=2)                          
            benchmarkCostGrad(petabModel, sensealgInfo, odeSolversGradientS[i], odeSolversGradientSName[i], pathSave, absTol, relTol, checkGradient=true, nRepeat=2, sparseJacobian=true)
        end
    end
end

if ARGS[1] == "Fix_parameters"

    if length(ARGS) == 1
        modelList = ["model_Bachmann_MSB2011", "model_Lucarelli_CellSystems2018", "model_Isensee_JCB2018"]
    elseif ARGS[2] == "Bachman"
        modelList = ["model_Bachmann_MSB2011"]
    elseif ARGS[2] == "Lucarelli"
        modelList = ["model_Lucarelli_CellSystems2018"]
    elseif ARGS[2] == "Isensee"
        modelList = ["model_Isensee_JCB2018"]        
    else
        println("ARGS must either be of length 2, or the second argument must be Bachman or Lucarelli not ", ARGS[2])
        exit(1)
    end

    nParamFix = parse(Int64, ARGS[3])
    seed = parse(Int64, ARGS[4])
    println("Number of parameter fixed = ", nParamFix)

    dirSave = joinpath(@__DIR__, "..", "..", "Intermediate", "Benchmarks", "Cost_grad_hess")
    pathSave = joinpath(dirSave, "Fix_parameters.csv")
    if !isdir(dirSave)
        mkpath(dirSave)
    end

    Random.seed!(seed)
    sensealgInfo = [[:ForwardDiff, nothing, "ForwardDiff"]]
    absTol, relTol = 1e-8, 1e-8

    for i in eachindex(modelList)

        local pathYML = getPathYmlFile(joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", modelList[i]))
        local petabModel = readPEtabModel(pathYML)
        nDynamicParameters = length(getNominalODEValues(petabModel))

        petabModelFewerParameters = getPEtabModelNparamFixed(petabModel, nParamFix)
        println("Number of dynamic parameters = ", length(getNominalODEValues(petabModelFewerParameters)))

        benchmarkCostGrad(petabModelFewerParameters, sensealgInfo[1], QNDF(), "QNDF", pathSave, absTol, relTol, 
                            checkGradient=true, nParamFixed=nParamFix, nRepeat=5, shouldSave=true)
        benchmarkCostGrad(petabModelFewerParameters, sensealgInfo[1], QNDF(), "QNDF", pathSave, absTol, relTol, 
                            checkGradient=true, nParamFixed=nParamFix, nRepeat=5, chunkSize=1, shouldSave=true)      
                                  
        if isdir(petabModelFewerParameters.dirModel) 
            rm(petabModelFewerParameters.dirModel, recursive=true)
            petabModelFewerParameters = nothing
            GC.gc(); GC.gc(); GC.gc()
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
    elseif ARGS[2] == "Isensee"
        modelList = ["model_Isensee_JCB2018"]            
    else
        println("ARGS must either be of length 2, or the second argument must be Bachman or Lucarelli not ", ARGS[2])
        exit(1)
    end

    dirSave = joinpath(@__DIR__, "..", "..", "Intermediate", "Benchmarks", "Cost_grad_hess")
    if !isdir(dirSave)
        mkpath(dirSave)
    end

    if ARGS[3] == "Shuffle" && length(ARGS) == 4 
        shuffleParameters = true
        seed = parse(Int64, ARGS[4])
        pathSave = joinpath(dirSave, "Test_chunks_shuffle_seed_" * ARGS[4] * ".csv")
    else
        shuffleParameters = false
        seed = 123
        Random.seed!(seed)
    end

    solversCheck = [[QNDF(), "QNDF"]]
    sensealgInfo = [[:ForwardDiff, nothing, "ForwardDiff"]]
    absTol, relTol = 1e-8, 1e-8

    for i in eachindex(modelList)

        local pathYML = getPathYmlFile(joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", modelList[i]))
        if shuffleParameters == false
            petabModel = readPEtabModel(pathYML)
        else
            _petabModel = readPEtabModel(pathYML)
            petabModel = getPEtabModelParamPermuted(_petabModel, seed=seed)
        end
        nDynamicParameters = length(getNominalODEValues(petabModel))
        chunkList = 1:nDynamicParameters
        for j in 1:30
            θ_est = getRandomModelParameters(petabModel, Rodas4P(), j, odeSolvers=false)
            for k in eachindex(chunkList)
                benchmarkCostGrad(petabModel, sensealgInfo[1], QNDF(), "QNDF", pathSave, absTol, relTol, 
                                checkGradient=true, nRepeat=5, chunkSize=chunkList[k], iParameter=j, _θ_est=θ_est)
            end
        end
        if shuffleParameters == true
            rm(petabModel.dirModel, recursive=true)
        end
    end
end
