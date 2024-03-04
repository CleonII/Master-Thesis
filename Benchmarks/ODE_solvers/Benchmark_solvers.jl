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

include(joinpath(pwd(), "src", "Solve_ODE", "Solve_ode_benchmark.jl"))
include(joinpath(pwd(), "src", "Solve_ODE", "Check_accuracy_ode_solver.jl"))

# For converting to SBML
include(joinpath(pwd(), "src", "SBML", "SBML_to_ModellingToolkit.jl"))

include(joinpath(@__DIR__, "..", "Common.jl"))


function getSolverInfo(sparseJacobian::Bool, solversCheck)

    lSolver1 = RFLUFactorization()
    lSolver2 = FastLUFactorization()
    lSolver3 = KrylovJL_GMRES()

    solverList = [Vern6(), "Vern6", "nonStiff", "OrdinaryDiffEq", "E_Runge_Kutta",
                  Vern7(), "Vern7", "nonStiff", "OrdinaryDiffEq", "E_Runge_Kutta",
                  Vern8(), "Vern8", "nonStiff", "OrdinaryDiffEq", "E_Runge_Kutta",
                  Vern9(), "Vern9", "nonStiff", "OrdinaryDiffEq", "E_Runge_Kutta",
                  Tsit5(), "Tsit5", "nonstiff", "OrdinaryDiffEq", "E_Runge_Kutta",
                  BS3(), "BS3", "nonstiff", "OrdinaryDiffEq", "E_Runge_Kutta",
                  BS5(), "BS5", "nonstiff", "OrdinaryDiffEq", "E_Runge_Kutta",
                  DP5(), "DP5", "nonStiff", "OrdinaryDiffEq", "E_Runge_Kutta",
                  DP8(), "DP8", "nonStiff", "OrdinaryDiffEq", "E_Runge_Kutta",
                  Feagin14(), "Feagin14", "nonStiff", "OrdinaryDiffEq", "E_Runge_Kutta",
                  VCABM(), "VCABM", "nonStiff", "OrdinaryDiffEq", "Adams_explicit",
                  Rosenbrock23(), "Rosenbrock23", "stiff", "OrdinaryDiffEq", "Rosenbrock",
                  TRBDF2(), "TRBDF2", "stiff", "OrdinaryDiffEq", "SDIRK",
                  Rodas4(), "Rodas4", "stiff", "OrdinaryDiffEq", "Rosenbrock",
                  Rodas4P(), "Rodas4P", "stiff", "OrdinaryDiffEq", "Rosenbrock",
                  Rodas5(), "Rodas5", "stiff", "OrdinaryDiffEq", "Rosenbrock",
                  Rodas5P(), "Rodas5P", "stiff", "OrdinaryDiffEq", "Rosenbrock",
                  QNDF(), "QNDF", "stiff", "OrdinaryDiffEq", "Multistep_bdf",
                  FBDF(), "FBDF", "stiff", "OrdinaryDiffEq", "Multistep_bdf",
                  Trapezoid(), "Trapezoid", "stiff", "OrdinaryDiffEq", "SDIRK",
                  KenCarp4(), "KenCarp4", "stiff", "OrdinaryDiffEq", "SDIRK",
                  KenCarp47(), "KenCarp47", "stiff", "OrdinaryDiffEq", "SDIRK",
                  Kvaerno5(), "Kvaerno5", "stiff", "OrdinaryDiffEq", "SDIRK",
                  RadauIIA3(), "RadauIIA3", "stiff", "OrdinaryDiffEq", "FIRK",
                  RadauIIA5(), "RadauIIA5", "stiff", "OrdinaryDiffEq", "FIRK",
                  AutoTsit5(Rosenbrock23()), "Tsit5Rosenbrock23", "composite", "OrdinaryDiffEq", "composite",
                  AutoVern7(Rodas5P()), "Vern7Rodas5P", "composite", "OrdinaryDiffEq", "composite",
                  AutoVern9(Rodas5P()), "Vern9Rodas5P", "composite", "OrdinaryDiffEq", "composite",
                  AutoVern7(Rodas4P()), "Vern7Rodas4P", "composite", "OrdinaryDiffEq", "composite",
                  AutoVern9(Rodas4P()), "Vern9Rodas4P", "composite", "OrdinaryDiffEq", "composite",
                  CVODE_BDF(), "CVODE_BDF_default", "stiff", "Sundials", "Multistep_bdf",
                  CVODE_BDF(linear_solver=:Dense), "CVODE_BDF_Dense", "stiff", "Sundials", "Multistep_bdf",
                  CVODE_BDF(linear_solver=:LapackDense), "CVODE_BDF_LapackDense", "stiff", "Sundials", "Multistep_bdf",
                  CVODE_BDF(linear_solver=:GMRES), "CVODE_BDF_GMRES", "stiff", "Sundials", "Multistep_bdf",
                  CVODE_Adams(linear_solver=:Dense), "CVODE_Adams_Dense", "nonStiff", "Sundials", "Adams_explicit",
                  CVODE_Adams(linear_solver=:LapackDense), "CVODE_Adams_LapackDense", "nonStiff", "Sundials", "Adams_explicit",
                  [:auto], "autoHint", "hint", "OrdinaryDiffEq", "hint",
                  [:nonstiff], "nonstiffHint", "hint", "OrdinaryDiffEq", "hint" ,
                  [:stiff], "stiffHint", "hint", "OrdinaryDiffEq", "hint",
                  RadauIIA5(linsolve=lSolver1), "RadauIIA5_RFLUF", "stiff", "OrdinaryDiffEq", "FIRK",
                  Rodas5(linsolve=lSolver1), "Rodas5_RFLUF", "stiff", "OrdinaryDiffEq", "Rosenbrock",
                  Rodas5P(linsolve=lSolver1), "Rodas5P_RFLUF", "stiff", "OrdinaryDiffEq", "Rosenbrock",
                  Rodas4P(linsolve=lSolver1), "Rodas4P_RFLUF", "stiff", "OrdinaryDiffEq", "Rosenbrock",
                  QNDF(linsolve=lSolver1), "QNDF_RFLUF", "stiff", "OrdinaryDiffEq", "Multistep_bdf",
                  FBDF(linsolve=lSolver1), "FBDF_RFLUF", "stiff", "OrdinaryDiffEq", "Multistep_bdf",
                  Kvaerno5(linsolve=lSolver1), "Kvaerno5_RFLUF", "stiff", "OrdinaryDiffEq", "SDIRK",
                  KenCarp4(linsolve=lSolver1), "KenCarp4_RFLUF", "stiff", "OrdinaryDiffEq", "SDIRK",
                  KenCarp47(linsolve=lSolver1), "KenCarp47_RFLUF", "stiff", "OrdinaryDiffEq", "SDIRK",
                  RadauIIA5(linsolve=lSolver2), "RadauIIA5_FastLU", "stiff", "OrdinaryDiffEq", "FIRK",
                  Rodas5(linsolve=lSolver2), "Rodas5_FastLU", "stiff", "OrdinaryDiffEq", "Rosenbrock",
                  Rodas5P(linsolve=lSolver2), "Rodas5P_FastLU", "stiff", "OrdinaryDiffEq", "Rosenbrock",
                  Rodas4P(linsolve=lSolver2), "Rodas4P_FastLU", "stiff", "OrdinaryDiffEq", "Rosenbrock",
                  Kvaerno5(linsolve=lSolver2), "Kvaerno5_FastLU", "stiff", "OrdinaryDiffEq", "SDIRK",
                  QNDF(linsolve=lSolver2), "QNDF_FastLU", "stiff", "OrdinaryDiffEq", "Multistep_bdf",
                  FBDF(linsolve=lSolver2), "FBDF_FastLU", "stiff", "OrdinaryDiffEq", "Multistep_bdf",
                  KenCarp4(linsolve=lSolver2), "KenCarp4_FastLU", "stiff", "OrdinaryDiffEq", "SDIRK",
                  KenCarp47(linsolve=lSolver2), "KenCarp47_FastLU", "stiff", "OrdinaryDiffEq", "SDIRK",
                  QNDF(linsolve=lSolver3), "QNDF_GMRES", "stiff", "OrdinaryDiffEq", "Multistep_bdf",
                  ImplicitDeuflhardExtrapolation(threading = false), "IDeuflar", "stiff", "OrdinaryDiffEq", "PIEP",
                  ImplicitHairerWannerExtrapolation(threading = false), "IWanner", "stiff", "OrdinaryDiffEq", "PIEP",
                  ImplicitEulerExtrapolation(threading = false), "IEuler", "stiff", "OrdinaryDiffEq", "PIEP",
                  ImplicitDeuflhardExtrapolation(threading = false, linsolve=lSolver2), "IDeuflar_FastLU", "stiff", "OrdinaryDiffEq", "PIEP",
                  ImplicitHairerWannerExtrapolation(threading = false, linsolve=lSolver2), "IWanner_FastLU", "stiff", "OrdinaryDiffEq", "PIEP",
                  ImplicitEulerExtrapolation(threading = false, linsolve=lSolver2), "IEuler_FastLU", "stiff", "OrdinaryDiffEq", "PIEP",
                  ImplicitDeuflhardExtrapolation(threading = false, linsolve=lSolver3), "IDeuflar_GMRES", "stiff", "OrdinaryDiffEq", "PIEP",
                  ImplicitHairerWannerExtrapolation(threading = false, linsolve=lSolver3), "IWanner_GMRES", "stiff", "OrdinaryDiffEq", "PIEP",
                  ImplicitEulerExtrapolation(threading = false, linsolve=lSolver3), "IEuler_GMRES", "stiff", "OrdinaryDiffEq", "PIEP"]

    lSolver1 = KLUFactorization()
    lSolver2 = KrylovJL_GMRES()
    solverListSparse = [RadauIIA5(), "RadauIIA5_S", "stiff", "OrdinaryDiffEq", "FIRK",
                        Rodas5(), "Rodas5_S", "stiff", "OrdinaryDiffEq", "Rosenbrock",
                        Rodas5P(), "Rodas5P_S", "stiff", "OrdinaryDiffEq", "Rosenbrock",
                        Rodas4P(), "Rodas4P_S", "stiff", "OrdinaryDiffEq", "Rosenbrock",
                        QNDF(), "QNDF_S", "stiff", "OrdinaryDiffEq", "Multistep_bdf",
                        KenCarp4(), "KenCarp4_S", "stiff", "OrdinaryDiffEq", "SDIRK",
                        KenCarp47(), "KenCarp47_S", "stiff", "OrdinaryDiffEq", "SDIRK",
                        Kvaerno5(), "Kvaerno5_S", "stiff", "OrdinaryDiffEq", "SDIRK",
                        FBDF(), "FBDF_S", "stiff", "OrdinaryDiffEq", "Multistep_bdf",
                        Rosenbrock23(), "Rosenbrock23_S", "stiff", "OrdinaryDiffEq", "Rosenbrock",
                        Rodas4(), "Rodas4_S", "stiff", "OrdinaryDiffEq", "Rosenbrock",
                        TRBDF2(), "TRBDF2_S", "stiff", "OrdinaryDiffEq", "SDIRK",
                        RadauIIA5(linsolve=lSolver1), "RadauIIA5_KLU_S", "stiff", "OrdinaryDiffEq", "FIRK",
                        Rodas5(linsolve=lSolver1), "Rodas5_KLU_S", "stiff", "OrdinaryDiffEq", "Rosenbrock",
                        Rodas5P(linsolve=lSolver1), "Rodas5P_KLU_S", "stiff", "OrdinaryDiffEq", "Rosenbrock",
                        Rodas4P(linsolve=lSolver1), "Rodas4P_KLU_S", "stiff", "OrdinaryDiffEq", "Rosenbrock",
                        QNDF(linsolve=lSolver1), "QNDF_KLU_S", "stiff", "OrdinaryDiffEq", "Multistep_bdf",
                        Kvaerno5(linsolve=lSolver1), "Kvaerno5_KLU_S", "stiff", "OrdinaryDiffEq", "SDIRK",
                        KenCarp4(linsolve=lSolver1), "KenCarp4_KLU_S", "stiff", "OrdinaryDiffEq", "SDIRK",
                        KenCarp47(linsolve=lSolver1), "KenCarp47_KLU_S", "stiff", "OrdinaryDiffEq", "SDIRK",
                        FBDF(linsolve=lSolver1), "FBDF_KLU_S", "stiff", "OrdinaryDiffEq", "Multistep_bdf",
                        RadauIIA5(linsolve=lSolver2), "RadauIIA5_GMRES_S", "stiff", "OrdinaryDiffEq", "FIRK",
                        Rodas5(linsolve=lSolver2), "Rodas5_GMRES_S", "stiff", "OrdinaryDiffEq", "Rosenbrock",
                        Rodas5P(linsolve=lSolver2), "Rodas5P_GMRES_S", "stiff", "OrdinaryDiffEq", "Rosenbrock",
                        Rodas4P(linsolve=lSolver2), "Rodas4P_GMRES_S", "stiff", "OrdinaryDiffEq", "Rosenbrock",
                        QNDF(linsolve=lSolver2), "QNDF_GMRES_S", "stiff", "OrdinaryDiffEq", "Multistep_bdf",
                        CVODE_BDF(linear_solver=:GMRES), "CVODE_BDF_GMRES_S", "stiff", "Sundials", "Multistep_bdf",
                        CVODE_BDF(linear_solver=:KLU), "CVODE_BDF_KLU_S", "stiff", "Sundials", "Multistep_bdf",]

    nRow = Int(length(solverList) / 5)
    solverList = reshape(solverList, (5, nRow))
    nRow = Int(length(solverListSparse) / 5)
    solverListSparse = reshape(solverListSparse, (5, nRow))

    if sparseJacobian == true && solversCheck == "all"
        return solverListSparse

    elseif sparseJacobian == true && solversCheck != "all"
        iUse = findall(x -> x ∈ solversCheckSparse, solverListSparse[2, :])
        return solverListSparse[:, iUse]

    elseif sparseJacobian == false && solversCheck == "all"
        return solverList

    else
        iUse = findall(x -> x ∈ solversCheck, solverList[2, :])
        return solverList[:, iUse]
    end
end

function get_highacc_runtime_sqdiff(petabModel, θ_dynamic)

    _θ_dynamic = θ_dynamic[:]
    # Process PeTab files into type-stable Julia structs
    experimentalConditions, measurementsData, parametersData, observablesData = readPEtabFiles(petabModel)
    parameterInfo = processParameters(parametersData)
    measurementInfo = processMeasurements(measurementsData, observablesData)
    simulationInfo = processSimulationInfo(petabModel, measurementInfo, parameterInfo)
    θ_indices = computeIndicesθ(parameterInfo, measurementInfo, petabModel)

    # Set model parameter values to those in the PeTab parameter data ensuring correct value of constant parameters
    setParamToFileValues!(petabModel.parameterMap, petabModel.stateMap, parameterInfo)

    # The time-span 5e3 is overwritten when performing actual forward simulations
    _odeProblem = ODEProblem{true, SciMLBase.FullSpecialize}(petabModel.odeSystem, petabModel.stateMap, (0.0, 5e3), petabModel.parameterMap, jac=true, sparse=false)
    odeProblem = remake(_odeProblem, p = convert.(Float64, _odeProblem.p), u0 = convert.(Float64, _odeProblem.u0))
    # In case we have provided a random start-guess
    θ_dynamic = _θ_dynamic
    # Change to parameters specified by paramVec
    changeODEProblemParameters!(odeProblem.p, odeProblem.u0, θ_dynamic, θ_indices, petabModel)
    # High accuracy ODE problem for check solver quality. BigFloats not compatible with sparse solvers
    _odeProblemAccuarcy = ODEProblem(petabModel.odeSystem, petabModel.stateMap, (0.0, 5e3), petabModel.parameterMap, jac=true, sparse=false)
    odeProblemAccuarcy = remake(_odeProblemAccuarcy, p = convert.(Float64, odeProblem.p), u0 = convert.(Float64, odeProblem.u0))
    # Functions to map experimental conditions and parameters correctly to the ODE model
    changeExperimentalCondition! = (pODEProblem, u0, conditionId) -> _changeExperimentalCondition!(pODEProblem, u0, conditionId, θ_dynamic, petabModel, θ_indices)

    highAccuracySolutions, could_solve_high = solveODEAllExperimentalConditions(odeProblem, changeExperimentalCondition!, simulationInfo, Rodas4P(), 1e-13, 1e-13, petabModel.computeTStops, nTimePointsSave=100)
    if could_solve_high == false
        @info "Failed with Rodas4P trying CVODE_BDF"
        highAccuracySolutions, could_solve_high = solveODEAllExperimentalConditions(odeProblem, changeExperimentalCondition!, simulationInfo, CVODE_BDF(), 1e-13, 1e-13, petabModel.computeTStops, nTimePointsSave=100)
    end
    if could_solve_high == false
        @info "Failed with CVODE_BDF trying KenCarp4"
        highAccuracySolutions, could_solve_high = solveODEAllExperimentalConditions(odeProblem, changeExperimentalCondition!, simulationInfo, KenCarp4(), 1e-13, 1e-13, petabModel.computeTStops, nTimePointsSave=100)
    end
    if could_solve_high == false
        @info "Failed produce high accuracy solution"
    end

    return highAccuracySolutions, could_solve_high
end

function getruntime_sqdiff(petabModel, tolsCheck, solver, θ_dynamic, solver_name, path_save::String, iparameter, highAccuracySolutions, could_solve_high_acc::Bool; nTimesRepat=2)

    _θ_dynamic = θ_dynamic[:]
    # Process PeTab files into type-stable Julia structs
    experimentalConditions, measurementsData, parametersData, observablesData = readPEtabFiles(petabModel)
    parameterInfo = processParameters(parametersData)
    measurementInfo = processMeasurements(measurementsData, observablesData)
    simulationInfo = processSimulationInfo(petabModel, measurementInfo, parameterInfo)
    θ_indices = computeIndicesθ(parameterInfo, measurementInfo, petabModel)

    # Set model parameter values to those in the PeTab parameter data ensuring correct value of constant parameters
    setParamToFileValues!(petabModel.parameterMap, petabModel.stateMap, parameterInfo)

    # The time-span 5e3 is overwritten when performing actual forward simulations
    _odeProblem = ODEProblem{true, SciMLBase.FullSpecialize}(petabModel.odeSystem, petabModel.stateMap, (0.0, 5e3), petabModel.parameterMap, jac=true, sparse=false)
    odeProblem = remake(_odeProblem, p = convert.(Float64, _odeProblem.p), u0 = convert.(Float64, _odeProblem.u0))
    # In case we have provided a random start-guess
    θ_dynamic = _θ_dynamic
    # Change to parameters specified by paramVec
    changeODEProblemParameters!(odeProblem.p, odeProblem.u0, θ_dynamic, θ_indices, petabModel)
    # Functions to map experimental conditions and parameters correctly to the ODE model
    changeExperimentalCondition! = (pODEProblem, u0, conditionId) -> _changeExperimentalCondition!(pODEProblem, u0, conditionId, θ_dynamic, petabModel, θ_indices)

    data_save = DataFrame()
    sqdiff_res = zeros(length(tolsCheck))
    runtime_res = similar(sqdiff_res)
    for (i, tol) in pairs(tolsCheck)

        absTol, relTol = tol, tol
        runTime = Vector{Float64}(undef, nTimesRepat)

        # Check the accuracy of the ODE solver by comparing with high accuracy solution. In case the squared sum
        # error cannot be computed the solver crashed and run time is not profiled
        # If we do not check accuracy, check that we can solve the model using the provided solver.
        sqDiffSolver = Float64
        if could_solve_high_acc == true
            try
                sqDiffSolver = computeAccuracyODESolver(odeProblem, highAccuracySolutions, changeExperimentalCondition!, simulationInfo, solver, absTol, relTol, petabModel.computeTStops)
            catch
                sqDiffSolver = Inf
            end
            @info "sqDiffSolver = $sqDiffSolver"
        end

        if could_solve_high_acc == false
            _, could_solve = solveODEAllExperimentalConditions(odeProblem, changeExperimentalCondition!, simulationInfo, solver, absTol, relTol, petabModel.computeTStops)
            if could_solve == false
                sqDiffSolver = Inf
            # If this else holds could solve local, but not high accuracy
            else
                sqDiffSolver = NaN
            end
        end

        canSolveModel = isinf(sqDiffSolver) ? false : true
        if canSolveModel
            for i in 1:nTimesRepat
                status, _runTime, retcode_ret = solveODEModelAllConditionsBenchmark(odeProblem, changeExperimentalCondition!, simulationInfo, solver, absTol, relTol, petabModel.computeTStops, onlySaveAtObservedTimes=true, savePreEqTime=true)
                runTime[i] = _runTime # seconds-precision
                GC.gc(); GC.gc();GC.gc()
            end
        else
            runTime .= Inf
        end

        sqdiff_res[i] = sqDiffSolver
        runtime_res[i] = mean(runTime)

        _data_save = DataFrame(model = petabModel.modelName,
                              iparameter=iparameter,
                              runtime = runtime_res[i],
                              sqdfiff = sqdiff_res[i],
                              solver = solver_name,
                              tol = tol)
        data_save = vcat(data_save, _data_save)
    end

    CSV.write(path_save, data_save, append=isfile(path_save))

    return runtime_res, sqdiff_res
end

function runBenchmarkOdeSolvers(petabModel::PEtabModel,
                                pathFileSave::String,
                                sparseJacobian::Bool;
                                solversCheck="all",
                                nTimesRepat::UInt=UInt(3),
                                tolsCheck=[(1e-6, 1e-6), (1e-9, 1e-9), (1e-12, 1e-12)],
                                absTolSS::Float64=1e-10,
                                relTolSS::Float64=1e-8,
                                _θ_dynamic=nothing,
                                checkAccuracy::Bool=true,
                                index_param=0)

    println("Working with model ", petabModel.modelName)

    # Process PeTab files into type-stable Julia structs
    experimentalConditions, measurementsData, parametersData, observablesData = readPEtabFiles(petabModel)
    parameterInfo = processParameters(parametersData)
    measurementInfo = processMeasurements(measurementsData, observablesData)
    simulationInfo = processSimulationInfo(petabModel, measurementInfo, parameterInfo, absTolSS=absTolSS, relTolSS=relTolSS)
    θ_indices = computeIndicesθ(parameterInfo, measurementInfo, petabModel)

    # Set model parameter values to those in the PeTab parameter data ensuring correct value of constant parameters
    setParamToFileValues!(petabModel.parameterMap, petabModel.stateMap, parameterInfo)

    # The time-span 5e3 is overwritten when performing actual forward simulations
    _odeProblem = ODEProblem{true, SciMLBase.FullSpecialize}(petabModel.odeSystem, petabModel.stateMap, (0.0, 5e3), petabModel.parameterMap, jac=true, sparse=sparseJacobian)
    odeProblem = remake(_odeProblem, p = convert.(Float64, _odeProblem.p), u0 = convert.(Float64, _odeProblem.u0))
    # In case we have provided a random start-guess
    if isnothing(_θ_dynamic)
        θ_dynamic = getNominalODEValues(petabModel)
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
        solverCategory = solverInfo[5, i]

        println("Trying solver = ", solverName)
        local retcode_ret = "nothing"
        # Crauste crashes as problem is to stiff (or horrible)
        if !((petabModel.modelName == "model_Crauste_CellSystems2017") && solver == AutoTsit5(Rosenbrock23()))
            for tol in tolsCheck

                absTol, relTol = tol
                runTime = Vector{Float64}(undef, nTimesRepat)

                # Check the accuracy of the ODE solver by comparing with high accuracy solution. In case the squared sum
                # error cannot be computed the solver crashed and run time is not profiled.
                # If we do not check accuracy, check that we can solve the model using the provided solver.
                local sqDiffSolver = Float64
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
                    odesols, canSolveModel = solveODEAllExperimentalConditions(odeProblem, changeExperimentalCondition!, simulationInfo, solver, absTol, relTol, petabModel.computeTStops, onlySaveAtObservedTimes=true)
                    sqDiffSolver = 0.0
                    for odesol in values(odesols)
                        if isnothing(odesol)
                            continue
                        end
                        retcode_ret = odesol.retcode
                        if !(retcode_ret == :Success || retcode_ret == :Terminated)
                            break
                        end
                    end

                    # Here we have a fun edge case where the solution errored in the
                    # LinearAlgebra, thus we must catch the error
                    println("retcode = ", retcode_ret)
                    if canSolveModel == false && retcode_ret == "nothing"
                        println("Got here")
                        local _retcode_ret = ""
                        # If it does not error, pre-eq simulation failed
                        try
                            _, _, _retcode_ret = solveODEModelAllConditionsBenchmark(odeProblem, changeExperimentalCondition!, simulationInfo, solver, absTol, relTol, petabModel.computeTStops, onlySaveAtObservedTimes=true, savePreEqTime=true)
                        catch e
                            if e isa BoundsError
                                _retcode_ret = "BoundsError"
                            elseif e isa DomainError
                                _retcode_ret = "DomainError"
                            elseif e isa SingularException
                                _retcode_ret = "SingularException"
                            else
                                rethrow(e)
                            end
                        end
                        retcode_ret = _retcode_ret
                        println("_retcode_ret error = ", _retcode_ret)
                    end
                end
                if canSolveModel == true
                    for i in 1:nTimesRepat
                        status, _runTime, retcode_ret = solveODEModelAllConditionsBenchmark(odeProblem, changeExperimentalCondition!, simulationInfo, solver, absTol, relTol, petabModel.computeTStops, onlySaveAtObservedTimes=true, savePreEqTime=true)
                        runTime[i] = _runTime # seconds-precision
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
                                     absTolSS=absTolSS,
                                     relTolSS=relTolSS,
                                     success = canSolveModel,
                                     runTime = runTime,
                                     sqDiff = sqDiffSolver,
                                     retcode=retcode_ret,
                                     iteration = 1:nTimesRepat,
                                     solverCategory = solverCategory,
                                     index_param = index_param)
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
                                     nParam = length(petabModel.parameterNames),
                                     reltol = absTol,
                                     abstol = relTol,
                                     absTolSS=absTolSS,
                                     relTolSS=relTolSS,
                                     success = false,
                                     runTime = NaN,
                                     sqDiff = Inf,
                                     retcode_ret=retcode_ret,
                                     iteration = 1:nTimesRepat,
                                     solverCategory = solverCategory,
                                     index_param = index_param)
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


# Need to read list with fails
function work_precision_failed_p(petabModel)

    i_fails_df = CSV.read(joinpath(petabModel.dirModel, "Julia_model_files", "I_ode_fail.csv"), DataFrame)
    i_fails = i_fails_df[!, :index]

    @info "length(i_fails) = " * string(length(i_fails))

    # Let us look closer at a work precision diagram where we fail
    tolsCheck = [1e-8, 1e-7, 1e-6, 1e-5]
    path_save = joinpath(@__DIR__, "..", "..", "Intermediate", "Benchmarks", "ODE_solvers", "Work_precision.csv")
    for j in i_fails

        θ_dynamic = getRandomModelParameters(petabModel, Rodas4P(), j)
        high_acc_sol, could_solve = get_highacc_runtime_sqdiff(petabModel, θ_dynamic)

        @info "j = $j"

        if petabModel.modelName == "model_Borghans_BiophysChem1997" && j in [45]
            skipcomposite = true
        else
            skipcomposite = false
        end
        if petabModel.modelName == "model_Crauste_CellSystems2017" && j in [141, 167, 180]
           skipcomposite = true
        else
            skipcomposite = false
        end
        if petabModel.modelName == "model_Weber_BMC2015" && j in [90]
           skip_qndf = true
        else
            skip_qndf = false
        end

        if skip_qndf == false
            runtime_QNDF, sqdiff_QNDF  = getruntime_sqdiff(petabModel, tolsCheck, QNDF(), θ_dynamic, "QNDF", path_save, j, high_acc_sol, could_solve)
        end
        runtime_CVODE, sqdiff_CVODE = getruntime_sqdiff(petabModel, tolsCheck, CVODE_BDF(), θ_dynamic, "CVODE_BDF", path_save, j, high_acc_sol, could_solve)
        if skipcomposite == false
            runtime_auto, sqdiff_auto = getruntime_sqdiff(petabModel, tolsCheck, AutoVern7(Rodas5P()), θ_dynamic, "Vern7(Rodas5P)", path_save, j, high_acc_sol, could_solve)
        end
        runtime_vern7, sqdiff_vern7 = getruntime_sqdiff(petabModel, tolsCheck, Vern7(), θ_dynamic, "Vern7", path_save, j, high_acc_sol, could_solve)
        runtime_rodas5P, sqdiff_rodas5P = getruntime_sqdiff(petabModel, tolsCheck, Rodas5P(), θ_dynamic, "Rodas5P", path_save, j, high_acc_sol, could_solve)
        runtime_tsit5, sqdiff_tsit5 = getruntime_sqdiff(petabModel, tolsCheck, Tsit5(), θ_dynamic, "Tsit5", path_save, j, high_acc_sol, could_solve)
    end
end

if ARGS[1] == "Test_all"

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
    for i in eachindex(modelList)
        dirModel = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", modelList[i])
        pathYML = getPathYmlFile(dirModel)
        petabModel = readPEtabModel(pathYML)

        runBenchmarkOdeSolvers(petabModel, pathFile, false, nTimesRepat=UInt(3), tolsCheck=tolsTry)
        GC.gc(); GC.gc();GC.gc()
        runBenchmarkOdeSolvers(petabModel, pathFileSparse, true, nTimesRepat=UInt(3), tolsCheck=tolsTry)
        GC.gc(); GC.gc();GC.gc()
    end
end


if ARGS[1] == "Large_models_random_p"

    dirSave = joinpath(@__DIR__, "..", "..", "Intermediate", "Benchmarks", "ODE_solvers")
    pathSave = joinpath(dirSave, "Large_models_random_p_new.csv")
    if !isdir(dirSave)
        mkpath(dirSave)
    end
    modelTest = ARGS[2]
    tolsTry = [(1e-8, 1e-8)]

    solversCheck = ["KenCarp47", "KenCarp4", "Kvaerno5", "QNDF", "TRBDF2", "FBDF", "Rodas5P", "CVODE_BDF_default",
                    "KenCarp47_RFLUF", "KenCarp4_RFLUF", "Kvaerno5_RFLUF", "QNDF_RFLUF", "FBDF_RFLUF", "Rodas5P_RFLUF",
                    "KenCarp47_FastLU", "KenCarp4_FastLU", "Kvaerno5_FastLU", "QNDF_FastLU", "FBDF_FastLU", "Rodas5P_FastLU"]
    solversCheckSparse = ["KenCarp47_S", "KenCarp4_S", "Kvaerno5_S", "QNDF_S", "TRBDF2_S", "FBDF_S", "Rodas5P_S", "CVODE_BDF_KLU_S",
                          "KenCarp47_KLU_S", "KenCarp4_KLU_S", "Kvaerno5_KLU_S", "QNDF_KLU_S", "FBDF_KLU_S", "Rodas5P_KLU_S"]

    dirModel = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", modelTest)
    pathYML = getPathYmlFile(dirModel)
    petabModel = readPEtabModel(pathYML)

    # 50 random vector
    for j in 1:50
        if j == 1
            θ_dynamic = getNominalODEValues(petabModel)
        else
            θ_dynamic = getRandomModelParameters(petabModel, Rodas4P(), j)
        end
        runBenchmarkOdeSolvers(petabModel, pathSave, false, nTimesRepat=UInt(3),
                               solversCheck=solversCheck, tolsCheck=tolsTry, checkAccuracy=false,
                               _θ_dynamic=θ_dynamic, index_param=j)
        runBenchmarkOdeSolvers(petabModel, pathSave, true, nTimesRepat=UInt(3),
                               solversCheck=solversCheckSparse, tolsCheck=tolsTry, checkAccuracy=false,
                               _θ_dynamic=θ_dynamic, index_param=j)
    end
end


# Test set of stiff and non-stiff solvers for random parameters values
if ARGS[1] == "Test_random_parameter"

    dirSave = joinpath(@__DIR__, "..", "..", "Intermediate", "Benchmarks", "ODE_solvers")
    pathSave = joinpath(dirSave, "Random_parameters_new.csv")
    if !isdir(dirSave)
        mkpath(dirSave)
    end

    modelList = ["model_Perelson_Science1996", "model_Zhao_QuantBiol2020", "model_Crauste_CellSystems2017",
                 "model_Bertozzi_PNAS2020", "model_Borghans_BiophysChem1997", "model_Giordano_Nature2020",
                 "model_Bruno_JExpBot2016", "model_Okuonghae_ChaosSolitonsFractals2020", "model_Fiedler_BMC2016",
                 "model_Bachmann_MSB2011", "model_Weber_BMC2015", "model_Schwen_PONE2014"]
    solversCheck = ["QNDF", "Rodas5P", "CVODE_BDF_default", "Vern7", "Tsit5", "Vern7Rodas5P"]
    tolsTry = [(1e-8, 1e-8)]
    for i in eachindex(modelList)

        dirModel = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", modelList[i])
        pathYML = getPathYmlFile(dirModel)
        petabModel = readPEtabModel(pathYML)

        # 40 random vector
        for j in 1:200
            θ_dynamic = getRandomModelParameters(petabModel, Rodas4P(), j; nParamCube=200)
            runBenchmarkOdeSolvers(petabModel, pathSave, false, nTimesRepat=UInt(1),
                                   solversCheck=solversCheck, tolsCheck=tolsTry, checkAccuracy=false,
                                   _θ_dynamic=θ_dynamic)
        end
    end
end

# Test set of stiff and non-stiff solvers for random parameters values
if ARGS[1] == "Test_random_parameter_work_prec"

    dirSave = joinpath(@__DIR__, "..", "..", "Intermediate", "Benchmarks", "ODE_solvers")
    pathSave = joinpath(dirSave, "Random_parameters_work_precision.csv")
    if !isdir(dirSave)
        mkpath(dirSave)
    end

    modelList = ["model_Borghans_BiophysChem1997", "model_Perelson_Science1996", "model_Zhao_QuantBiol2020", "model_Crauste_CellSystems2017",
                 "model_Bertozzi_PNAS2020", "model_Giordano_Nature2020",
                 "model_Bruno_JExpBot2016", "model_Okuonghae_ChaosSolitonsFractals2020", "model_Fiedler_BMC2016",
                 "model_Bachmann_MSB2011", "model_Weber_BMC2015", "model_Schwen_PONE2014"]
    for i in eachindex(modelList)


        dirModel = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", modelList[i])
        pathYML = getPathYmlFile(dirModel)
        petabModel = readPEtabModel(pathYML)
        i_fails_path = joinpath(petabModel.dirModel, "Julia_model_files", "I_ode_fail.csv")
        if !isfile(i_fails_path)
            continue
        end
        work_precision_failed_p(petabModel)
    end
end
