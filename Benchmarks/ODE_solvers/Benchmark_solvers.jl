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
using SciMLSensitivity
using Zygote
using Symbolics
using Sundials
using YAML


# Relevant PeTab structs for compuations 
include(joinpath(pwd(), "src", "PeTab_structs.jl"))

# PeTab importer to get cost, grad etc 
include(joinpath(pwd(), "src", "Create_PEtab_model.jl"))

# HyperCube sampling 
include(joinpath(pwd(), "src", "Optimizers", "Lathin_hypercube.jl"))

include(joinpath(pwd(), "src", "Solve_ODE", "Solve_ode_benchmark.jl"))
include(joinpath(pwd(), "src", "Solve_ODE", "Check_accuracy_ode_solver.jl"))

# For converting to SBML 
include(joinpath(pwd(), "src", "SBML", "SBML_to_ModellingToolkit.jl"))

include(joinpath(@__DIR__, "..", "Common.jl"))


function getSolverInfo(sparseJacobian::Bool, solversCheck)

    lSolver1 = RFLUFactorization()
    lSolver2 = FastLUFactorization()
    lSolver3 = KrylovJL_GMRES()

    solverList = [Vern6(), "Vern6", "nonStiff", "OrdinaryDiffEq", 
                  Vern7(), "Vern7", "nonStiff", "OrdinaryDiffEq", 
                  Vern8(), "Vern8", "nonStiff", "OrdinaryDiffEq", 
                  Vern9(), "Vern9", "nonStiff", "OrdinaryDiffEq", 
                  Tsit5(), "Tsit5", "nonstiff", "OrdinaryDiffEq", 
                  BS3(), "BS3", "nonstiff", "OrdinaryDiffEq", 
                  BS5(), "BS5", "nonstiff", "OrdinaryDiffEq", 
                  DP5(), "DP5", "nonStiff", "OrdinaryDiffEq", 
                  DP8(), "DP8", "nonStiff", "OrdinaryDiffEq", 
                  Feagin14(), "Feagin14", "nonStiff", "OrdinaryDiffEq", 
                  VCABM(), "VCABM", "nonStiff", "OrdinaryDiffEq", 
                  Rosenbrock23(), "Rosenbrock23", "stiff", "OrdinaryDiffEq", 
                  TRBDF2(), "TRBDF2", "stiff", "OrdinaryDiffEq", 
                  Rodas4(), "Rodas4", "stiff", "OrdinaryDiffEq", 
                  Rodas4P(), "Rodas4P", "stiff", "OrdinaryDiffEq", 
                  Rodas5(), "Rodas5", "stiff", "OrdinaryDiffEq", 
                  QNDF(), "QNDF", "stiff", "OrdinaryDiffEq", 
                  FBDF(), "FBDF", "stiff", "OrdinaryDiffEq", 
                  Trapezoid(), "Trapezoid", "stiff", "OrdinaryDiffEq", 
                  KenCarp4(), "KenCarp4", "stiff", "OrdinaryDiffEq", 
                  Kvaerno5(), "Kvaerno5", "stiff", "OrdinaryDiffEq", 
                  RadauIIA3(), "RadauIIA3", "stiff", "OrdinaryDiffEq", 
                  RadauIIA5(), "RadauIIA5", "stiff", "OrdinaryDiffEq", 
                  AutoTsit5(Rosenbrock23()), "Tsit5Rosenbrock23", "composite", "OrdinaryDiffEq", 
                  AutoVern7(Rodas5()), "Vern7Rodas5", "composite", "OrdinaryDiffEq", 
                  AutoVern7(Rodas4P()), "Vern7Rodas4P", "composite", "OrdinaryDiffEq", 
                  AutoVern9(Rodas4P()), "Vern9Rodas4P", "composite", "OrdinaryDiffEq", 
                  CVODE_BDF(), "CVODE_BDF_default", "stiff", "Sundials",
                  CVODE_BDF(linear_solver=:Dense), "CVODE_BDF_Dense", "stiff", "Sundials", 
                  CVODE_BDF(linear_solver=:LapackDense), "CVODE_BDF_LapackDense", "stiff", "Sundials", 
                  CVODE_BDF(linear_solver=:GMRES), "CVODE_BDF_GMRES", "stiff", "Sundials", 
                  CVODE_Adams(linear_solver=:Dense), "CVODE_Adams_Dense", "nonStiff", "Sundials", 
                  CVODE_Adams(linear_solver=:LapackDense), "CVODE_Adams_LapackDense", "nonStiff", "Sundials", 
                  # In Julia v1.8.1 these crash for problems solvers as there is problem with ccall
                  #ARKODE(Sundials.Explicit(), order=4), "ARKODE_Exp4", "nonStiff", "Sundials", 
                  #ARKODE(Sundials.Explicit(), order=8), "ARKODE_Exp8", "nonStiff", "Sundials", 
                  #ARKODE(Sundials.Implicit(), order=3), "ARKODE_Imp3", "stiff", "Sundials", 
                  #ARKODE(Sundials.Implicit(), order=5), "ARKODE_Imp5", "stiff", "Sundials",
                  [:auto], "autoHint", "hint", "OrdinaryDiffEq", 
                  [:nonstiff], "nonstiffHint", "hint", "OrdinaryDiffEq", 
                  [:stiff], "stiffHint", "hint", "OrdinaryDiffEq", 
                  RadauIIA5(linsolve=lSolver1), "RadauIIA5_RFLUF", "stiff", "OrdinaryDiffEq",
                  Rodas5(linsolve=lSolver1), "Rodas5_RFLUF", "stiff", "OrdinaryDiffEq",
                  Rodas4P(linsolve=lSolver1), "Rodas4P_RFLUF", "stiff", "OrdinaryDiffEq",
                  QNDF(linsolve=lSolver1), "QNDF_RFLUF", "stiff", "OrdinaryDiffEq",
                  RadauIIA5(linsolve=lSolver2), "RadauIIA5_FastLU", "stiff", "OrdinaryDiffEq",
                  Rodas5(linsolve=lSolver2), "Rodas5_FastLU", "stiff", "OrdinaryDiffEq",
                  Rodas4P(linsolve=lSolver2), "Rodas4P_FastLU", "stiff", "OrdinaryDiffEq",
                  QNDF(linsolve=lSolver2), "QNDF_FastLU", "stiff", "OrdinaryDiffEq",
                  QNDF(linsolve=lSolver3), "QNDF_GMRES", "stiff", "OrdinaryDiffEq",
                  ImplicitDeuflhardExtrapolation(threading = true), "IDeuflar_Thread", "stiff", "OrdinaryDiffEq",
                  ImplicitHairerWannerExtrapolation(threading = true), "IWanner_Thread", "stiff", "OrdinaryDiffEq",
                  ImplicitEulerExtrapolation(threading = true), "IEuler_Thread", "stiff", "OrdinaryDiffEq",
                  ImplicitDeuflhardExtrapolation(threading = false), "IDeuflar", "stiff", "OrdinaryDiffEq",
                  ImplicitHairerWannerExtrapolation(threading = false), "IWanner", "stiff", "OrdinaryDiffEq",
                  ImplicitEulerExtrapolation(threading = false), "IEuler", "stiff", "OrdinaryDiffEq",
                  ImplicitDeuflhardExtrapolation(threading = false, linsolve=lSolver2), "IDeuflar_FastLU", "stiff", "OrdinaryDiffEq",
                  ImplicitHairerWannerExtrapolation(threading = false, linsolve=lSolver2), "IWanner_FastLU", "stiff", "OrdinaryDiffEq",
                  ImplicitEulerExtrapolation(threading = false, linsolve=lSolver2), "IEuler_FastLU", "stiff", "OrdinaryDiffEq", 
                  ImplicitDeuflhardExtrapolation(threading = false, linsolve=lSolver3), "IDeuflar_GMRES", "stiff", "OrdinaryDiffEq",
                  ImplicitHairerWannerExtrapolation(threading = false, linsolve=lSolver3), "IWanner_GMRES", "stiff", "OrdinaryDiffEq",
                  ImplicitEulerExtrapolation(threading = false, linsolve=lSolver3), "IEuler_GMRES", "stiff", "OrdinaryDiffEq"]                
                
    lSolver1 = KLUFactorization()
    lSolver2 = KrylovJL_GMRES()
    solverListSparse = [RadauIIA5(), "RadauIIA5_S", "stiff", "OrdinaryDiffEq",
                        Rodas5(), "Rodas5_S", "stiff", "OrdinaryDiffEq",
                        Rodas4P(), "Rodas4P_S", "stiff", "OrdinaryDiffEq",
                        QNDF(), "QNDF_S", "stiff", "OrdinaryDiffEq",
                        KenCarp4(), "KenCarp4_S", "stiff", "OrdinaryDiffEq",
                        FBDF(), "FBDF_S", "stiff", "OrdinaryDiffEq",
                        Rosenbrock23(), "Rosenbrock23_S", "stiff", "OrdinaryDiffEq",
                        Rodas4(), "Rodas4_S", "stiff", "OrdinaryDiffEq",
                        TRBDF2(), "TRBDF2_S", "stiff", "OrdinaryDiffEq",
                        RadauIIA5(linsolve=lSolver1), "RadauIIA5_KLU_S", "stiff", "OrdinaryDiffEq",
                        Rodas5(linsolve=lSolver1), "Rodas5_KLU_S", "stiff", "OrdinaryDiffEq",
                        Rodas4P(linsolve=lSolver1), "Rodas4P_KLU_S", "stiff", "OrdinaryDiffEq",
                        QNDF(linsolve=lSolver1), "QNDF_KLU_S", "stiff", "OrdinaryDiffEq",
                        RadauIIA5(linsolve=lSolver2), "RadauIIA5_GMRES_S", "stiff", "OrdinaryDiffEq",
                        Rodas5(linsolve=lSolver2), "Rodas5_GMRES_S", "stiff", "OrdinaryDiffEq",
                        Rodas4P(linsolve=lSolver2), "Rodas4P_GMRES_S", "stiff", "OrdinaryDiffEq",
                        QNDF(linsolve=lSolver2), "QNDF_GMRES_S", "stiff", "OrdinaryDiffEq",
                        CVODE_BDF(linear_solver=:GMRES), "CVODE_BDF_GMRES_S", "stiff", "Sundials",
                        CVODE_BDF(linear_solver=:KLU), "CVODE_BDF_KLU_S", "stiff", "Sundials"]
                  
    nRow = Int(length(solverList) / 4)
    solverList = reshape(solverList, (4, nRow))
    nRow = Int(length(solverListSparse) / 4)
    solverListSparse = reshape(solverListSparse, (4, nRow))
    
    if sparseJacobian == true && solversCheck == "all"
        return solverListSparseMat

    elseif sparseJacobian == true && solversCheck != "all"
        iUse = [findfirst(x -> x == solversCheck[i], solverListSparse[2, :]) for i in eachindex(solversCheck)]
        return solverListSparse[iUse, :]

    elseif sparseJacobian == false && solversCheck == "all"
        return solverList

    else
        iUse = [findfirst(x -> x == solversCheck[i], solverList[2, :]) for i in eachindex(solversCheck)]
        return solverList[iUse, :]
    end         
end


function runBenchmarkOdeSolvers(petabModel::PEtabModel, 
                                pathFileSave::String,
                                sparseJacobian::Bool;
                                solversCheck="all",
                                nTimesRepat::UInt=UInt(3), 
                                tolsCheck=[(1e-6, 1e-6), (1e-9, 1e-9), (1e-12, 1e-12)], 
                                _θ_dynamic=nothing,
                                checkAccuracy::Bool=true)

    println("Working with model ", petabModel.modelName)

    # Process PeTab files into type-stable Julia structs 
    experimentalConditions, measurementsData, parametersData, observablesData = readPEtabFiles(petabModel)
    parameterInfo = processParameters(parametersData) 
    measurementInfo = processMeasurements(measurementsData, observablesData) 
    simulationInfo = processSimulationInfo(petabModel, measurementInfo, absTolSS=1e-10, relTolSS=1e-8)
    θ_indices = computeIndicesθ(parameterInfo, measurementInfo, petabModel.odeSystem, experimentalConditions)
     
    # Set model parameter values to those in the PeTab parameter data ensuring correct value of constant parameters 
    setParamToFileValues!(petabModel.parameterMap, petabModel.stateMap, parameterInfo)
 
    # The time-span 5e3 is overwritten when performing actual forward simulations 
    _odeProblem = ODEProblem(petabModel.odeSystem, petabModel.stateMap, (0.0, 5e3), petabModel.parameterMap, jac=true, sparse=sparseJacobian)
    odeProblem = remake(_odeProblem, p = convert.(Float64, _odeProblem.p), u0 = convert.(Float64, _odeProblem.u0))
    # In case we have provided a random start-guess
    if isnothing(_θ_dynamic)
        θ_dynamic = getFileODEvalues(petabModel)
    else
        θ_dynamic = _θ_dynamic
    end
    # Change to parameters specified by paramVec
    changeODEProblemParameters!(odeProblem.p, odeProblem.u0, θ_dynamic, θ_indices, petabModel)                          

    # High accuracy ODE problem for check solver quality. BigFloats not compatible with sparse solvers 
    _odeProblemAccuarcy = ODEProblem(petabModel.odeSystem, petabModel.stateMap, (0.0, 5e3), petabModel.parameterMap, jac=true, sparse=false)
    odeProblemAccuarcy = remake(_odeProblemAccuarcy, p = convert.(Float64, odeProblem.p), u0 = convert.(Float64, odeProblem.u0))

    # Functions to map experimental conditions and parameters correctly to the ODE model 
    changeExperimentalCondition! = (pODEProblem, u0, conditionId) -> _changeExperimentalCondition!(pODEProblem, u0, conditionId, θ_dynamic, petabModel, θ_indices)
    
    # High accruacy ODE solution. Laske-model does not not solve with BigFloat 
    if checkAccuracy == true
        println("Computing high accuracy model")
        if petabModel.modelName != "model_Laske_PLOSComputBiol2019"
            highAccuracySolutions, statusAccuracy = computeHighAccuracyOdeSolution(odeProblemAccuarcy, changeExperimentalCondition!, simulationInfo, petabModel.computeTStops, absTol=1e-15, relTol=1e-15)
        else
            highAccuracySolutions, tmp = solveODEAllExperimentalConditions(odeProblem, changeExperimentalCondition!, simulationInfo, KenCarp58(), 1e-15, 1e-15, petabModel.computeTStops, nTimePointsSave=100)
            statusAccuracy = true
        end
        println("Done with high accuracy solution and status = $statusAccuracy")
    end

    solverInfo = getSolverInfo(sparseJacobian, solversCheck)

    if checkAccuracy == true && statusAccuracy != true
        # Log failure to disk 
        println("High accuracy solver failed for ", petabModel.modelName)
        open(joinpath(pwd(), "Benchmark", "ODE_solvers", "Log.txt"), "a+") do io
            println(io, "Failed with high accuracy solution for $modelFile")
        end
        return 
    end

    for i in 1:size(solverInfo)[2]

        solver = solverInfo[1, i]
        solverName = solverInfo[2, i]
        solverType = solverInfo[3, i]
        solverLib = solverInfo[4, i]

        println("Trying solver = ", solverName)
        # Crauste crashes as problem is to stiff 
        if !((petabModel.modelName == "model_Crauste_CellSystems2017") && solver == AutoTsit5(Rosenbrock23())) 
            for tol in tolsCheck
                
                absTol, relTol = tol
                runTime = Vector{Float64}(undef, nTimesRepat)

                # Check the accuracy of the ODE solver by comparing with high accuracy solution. In case the squared sum 
                # error cannot be computed the solver crashed and run time is not profiled.
                # If we do not check accuracy, check that we can solve the model using the provided solver.
                local sqDiffSolver = Float64
                _sol = computeAccuracyODESolver(odeProblem, highAccuracySolutions, changeExperimentalCondition!, simulationInfo, solver, absTol, relTol, petabModel.computeTStops)
                if checkAccuracy == true
                    try
                        sqDiffSolver = computeAccuracyODESolver(odeProblem, highAccuracySolutions, changeExperimentalCondition!, simulationInfo, solver, absTol, relTol, petabModel.computeTStops)
                    catch 
                        sqDiffSolver = Inf
                    end
                    println("sqDiffSolver = ", sqDiffSolver)
                    canSolveModel = isinf(sqDiffSolver) ? false : true
                else
                    # Here we need to solver the model for precompiliation purposes 
                    tmp, canSolveModel = solveODEAllExperimentalConditions(odeProblem, changeExperimentalCondition!, simulationInfo, solver, absTol, relTol, petabModel.computeTStops, onlySaveAtObservedTimes=true)
                    sqDiffSolver = 0.0
                end
                    
                if canSolveModel == true
                    for i in 1:nTimesRepat
                        status, runTime = solveODEModelAllConditionsBenchmark(odeProblem, changeExperimentalCondition!, simulationInfo, solver, absTol, relTol, petabModel.computeTStops, onlySaveAtObservedTimes=true, savePreEqTime=true) 
                        runTime[i] = runTime # seconds-precision
                        GC.gc(); GC.gc();GC.gc()
                    end
                else
                    runTime .= NaN
                end
                dataSave = DataFrame(model = petabModel.modelName, 
                                     solver = solverName, 
                                     solverType = solverType,
                                     solverLib = solverLib,
                                     nStates = length(petabModel.stateNames),
                                     nParam = length(petabModel.parameterNames),
                                     reltol = relTol, 
                                     abstol = absTol, 
                                     success = canSolveModel, 
                                     runTime = runTime, 
                                     sqDiff = sqDiffSolver, 
                                     iteration = 1:nTimesRepat)
                if isfile(pathFileSave)
                    CSV.write(pathFileSave, dataSave, append = true)
                else
                    CSV.write(pathFileSave, dataSave)
                end
                GC.gc(); GC.gc();GC.gc()
            end
        else
            for tol in tolsCheck

                absTol, relTol = tol
                dataSave = DataFrame(model = petabModel.modelName, 
                                     solver = solverName, 
                                     solverType = solverType,
                                     solverLib = solverLib,
                                     nStates = length(petabModel.stateNames)-1,
                                     nParam = length(petabModel.paramNames),
                                     reltol = absTol, 
                                     abstol = relTol, 
                                     success = false, 
                                     runTime = NaN,
                                     sqDiff = Inf, 
                                     iteration = 1:nTimesRepat)                
                if isfile(pathFileSave)
                    CSV.write(pathFileSave, dataSave, append = true)
                else
                    CSV.write(pathFileSave, dataSave)
                end
            end
        end
        GC.gc()
    end

    GC.gc()
end


#if ARGS[1] == "Test_all"

    dirSave = joinpath(@__DIR__, "..", "..", "Intermediate", "Benchmarks", "ODE_solvers")
    pathFile = joinpath(dirSave, "All_models.csv")
    pathFileSparse = joinpath(dirSave, "All_models_sparse_jacobian.csv")
    if !isdir(dirSave)
        mkpath(dirSave)
    end

    modelList = ["model_Beer_MolBioSystems2014", "model_Blasi_CellSystems2016", "model_Weber_BMC2015", "model_Schwen_PONE2014", "model_Alkan_SciSignal2018", 
                "model_Bachmann_MSB2011", "model_Bertozzi_PNAS2020", "model_Boehm_JProteomeRes2014", 
                "model_Borghans_BiophysChem1997", "model_Brannmark_JBC2010", "model_Bruno_JExpBot2016", "model_Crauste_CellSystems2017", 
                "model_Elowitz_Nature2000", "model_Fiedler_BMC2016", "model_Fujita_SciSignal2010", "model_Giordano_Nature2020", 
                "model_Isensee_JCB2018", "model_Laske_PLOSComputBiol2019", "model_Lucarelli_CellSystems2018", "model_Okuonghae_ChaosSolitonsFractals2020", 
                "model_Oliveira_NatCommun2021", "model_Perelson_Science1996", "model_Rahman_MBS2016", 
                "model_SalazarCavazos_MBoC2020", "model_Sneyd_PNAS2002", "model_Zhao_QuantBiol2020", "model_Zheng_PNAS2012"]                    

    tolsTry = [(1e-16, 1e-8), (1e-8, 1e-8), (1e-6, 1e-6)]            
    tolsTry = [(1e-6, 1e-6)]            
    for i in eachindex(modelList)
        dirModel = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", modelList[i])
        pathYML = getPathYmlFile(dirModel)
        petabModel = readPEtabModel(pathYML)

        runBenchmarkOdeSolvers(petabModel, pathFile, false, nTimesRepat=UInt(3), tolsCheck=tolsTry)
        GC.gc(); GC.gc();GC.gc()
        runBenchmarkOdeSolvers(petabModel, pathFileSparse, true, nTimesRepat=UInt(3), tolsCheck=tolsTry)
        GC.gc(); GC.gc();GC.gc()
    end
#end


if ARGS[1] == "Large_models"

    dirSave = joinpath(@__DIR__, "..", "..", "Intermediate", "Benchmarks", "ODE_solvers")
    pathSave = joinpath(dirSave, "Large_models.csv") 
    if !isdir(dirSave)
        mkpath(dirSave)
    end

    modelList = ["model_Chen_MSB2009"]
    tolsTry = [(1e-6, 1e-6)]            
    solversCheck = ["KenCarp4", "QNDF", "TRBDF2", "FBDF", "Rodas4P", "CVODE_BDF_default"]

    for i in eachindex(modelList)
        dirModel = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", modelList[i])
        pathYML = getPathYmlFile(dirModel)
        petabModel = readPEtabModel(pathYML)

        runBenchmarkOdeSolvers(petabModel, pathSave, false, nTimesRepat=UInt(3), solversCheck=solversCheck, tolsCheck=tolsTry, checkAccuracy=false)    
    end

    # Now try with sparse Jacobian 
    solversCheckSparse = ["KenCarp4_S", "QNDF_S", "TRBDF2_S", "FBDF_S", "Rodas4P_S", "CVODE_BDF_KLU_S"]
    for i in eachindex(modelList)
        dirModel = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", modelList[i])
        pathYML = getPathYmlFile(dirModel)
        petabModel = readPEtabModel(pathYML)

        runBenchmarkOdeSolvers(petabModel, pathSave, true, nTimesRepat=UInt(3), solversCheck=solversCheckSparse, tolsCheck=tolsTry, checkAccuracy=false)    
    end
end


# Test set of stiff and non-stiff solvers for random parameters values 
if ARGS[1] == "Test_random_parameter"

    dirSave = joinpath(@__DIR__, "..", "..", "Intermediate", "Benchmarks", "ODE_solvers")
    pathSave = joinpath(dirSave, "Random_parameters.csv")
    if !isdir(dirSave)
        mkpath(dirSave)
    end

    modelList = ["model_Perelson_Science1996", "model_Zhao_QuantBiol2020", "model_Crauste_CellSystems2017", "model_Fiedler_BMC2016", 
                 "model_Bruno_JExpBot2016", "model_Okuonghae_ChaosSolitonsFractals2020", "model_Schwen_PONE2014", 
                 "model_Bachmann_MSB2011", "model_Brannmark_JBC2010", "model_Lucarelli_CellSystems2018", 
                 "model_Isensee_JCB2018", "model_Weber_BMC2015"]


    solversCheck = ["Rodas5", "QNDF", "Rodas4P", "CVODE_BDF_default", "Vern7", "Tsit5", "Vern6", "Vern7Rodas4P"]
    tolsTry = [(1e-8, 1e-8)]            
    for i in eachindex(modelList)

        dirModel = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", modelList[i])
        pathYML = getPathYmlFile(dirModel)
        petabModel = readPEtabModel(pathYML)
        
        # 40 random vector 
        for j in 1:40
            θ_dynamic = getRandomModelParameters(petabModel, Rodas4P(), j)
            runBenchmarkOdeSolvers(petabModel, pathSave, false, nTimesRepat=UInt(1), 
                                   solversCheck=solversCheck, tolsCheck=tolsTry, checkAccuracy=false, 
                                   _θ_dynamic=θ_dynamic)    
        end
    end
end