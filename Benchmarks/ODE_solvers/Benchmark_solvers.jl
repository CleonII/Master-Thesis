using ModelingToolkit 
using DifferentialEquations
using ODEInterfaceDiffEq
using Sundials
using LSODA
using DataFrames
using CSV 
using ForwardDiff
using ReverseDiff
using StatsBase
using Random
using LinearAlgebra
using Distributions
using Printf
using Plots
using BenchmarkTools


# Relevant PeTab structs for compuations 
include(joinpath(pwd(), "src", "PeTab_structs.jl"))

# Functions for solving ODE system 
include(joinpath(pwd(), "src", "Solve_ODE_model", "Solve_ode_model.jl"))
include(joinpath(pwd(), "src", "Solve_ODE_model", "Check_accuracy_ode_solver.jl"))

# PeTab importer to get cost, grad etc 
include(joinpath(pwd(), "src", "PeTab_importer", "Create_cost_grad_hessian.jl"))

# HyperCube sampling 
include(joinpath(pwd(), "src", "Optimizers", "Lathin_hypercube.jl"))

# For converting to SBML 
include(joinpath(pwd(), "src", "SBML", "SBML_to_ModellingToolkit.jl"))


function getSolverInfo(sparseList::Bool)

    lSolver1 = RFLUFactorization()
    lSolver2 = FastLUFactorization()
    lSolver3 = KrylovJL_GMRES()

    solverList = [Vern6(), "Vern6", "nonStiff", "OrdinaryDiffEq", 
                  Vern7(), "Vern7", "nonStiff", "OrdinaryDiffEq", 
                  Vern8(), "Vern8", "nonStiff", "OrdinaryDiffEq", 
                  Vern9(), "Vern9", "nonStiff", "OrdinaryDiffEq", 
                  Tsit5(), "Tsit5", "nonstiff", "OrdinaryDiffEq", 
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
                  QNDF2(), "QNDF2", "stiff", "OrdinaryDiffEq", 
                  FBDF(), "FDBF", "stiff", "OrdinaryDiffEq", 
                  Trapezoid(), "Trapezoid", "stiff", "OrdinaryDiffEq", 
                  KenCarp4(), "KenCarp4", "stiff", "OrdinaryDiffEq", 
                  Kvaerno5(), "Kvaerno5", "stiff", "OrdinaryDiffEq", 
                  RadauIIA3(), "RadauIIA3", "stiff", "OrdinaryDiffEq", 
                  RadauIIA5(), "RadauIIA5", "stiff", "OrdinaryDiffEq", 
                  AutoTsit5(Rosenbrock23()), "Tsit5Rosenbrock23", "composite", "OrdinaryDiffEq", 
                  AutoVern7(Rodas5()), "Vern7Rodas5", "composite", "OrdinaryDiffEq", 
                  AutoVern7(Rodas4P()), "Vern7Rodas4P", "composite", "OrdinaryDiffEq", 
                  AutoVern9(Rodas4P()), "Vern9Rodas4P", "composite", "OrdinaryDiffEq", 
                  lsoda(), "lsoda", "composite", "LSODA", 
                  CVODE_BDF(linear_solver=:Dense), "CVODE_BDF_Dense", "stiff", "Sundials", 
                  CVODE_BDF(linear_solver=:LapackDense), "CVODE_BDF_LapackDense", "stiff", "Sundials", 
                  CVODE_BDF(linear_solver=:GMRES), "CVODE_BDF_GMRES", "stiff", "Sundials", 
                  CVODE_Adams(linear_solver=:Dense), "CVODE_Adams_Dense", "nonStiff", "Sundials", 
                  CVODE_Adams(linear_solver=:LapackDense), "CVODE_Adams_LapackDense", "nonStiff", "Sundials", 
                  ARKODE(Sundials.Explicit(), order=4), "ARKODE_Exp4", "nonStiff", "Sundials", 
                  ARKODE(Sundials.Explicit(), order=8), "ARKODE_Exp8", "nonStiff", "Sundials", 
                  ARKODE(Sundials.Implicit(), order=3), "ARKODE_Imp3", "stiff", "Sundials", 
                  ARKODE(Sundials.Implicit(), order=5), "ARKODE_Imp5", "stiff", "Sundials",
                  dopri5(), "dopri5", "nonStiff", "ODEInterface", 
                  dop853(), "dop853", "nonStiff", "ODEInterface", 
                  radau(), "radau", "stiff", "ODEInterface", 
                  radau5(), "radau5", "stiff", "ODEInterface", 
                  rodas(), "rodas", "stiff", "ODEInterface", 
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
                  RadauIIA5(linsolve=lSolver3), "RadauIIA5_GMRES", "stiff", "OrdinaryDiffEq",
                  Rodas5(linsolve=lSolver3), "Rodas5_GMRES", "stiff", "OrdinaryDiffEq",
                  Rodas4P(linsolve=lSolver3), "Rodas4P_GMRES", "stiff", "OrdinaryDiffEq",
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
                  
    solverListMat = Array{Any, 2}(undef, (Int(length(solverList)/4), 4))
    for i in 1:size(solverListMat)[1]
        iStart = (i-1)*4 + 1
        iEnd = i*4
        solverListMat[i, :] .= solverList[iStart:iEnd]
    end
    solverListSparseMat = Array{Any, 2}(undef, (Int(length(solverListSparse)/4), 4))
    for i in 1:size(solverListSparseMat)[1]
        iStart = (i-1)*4 + 1
        iEnd = i*4
        solverListSparseMat[i, :] .= solverListSparse[iStart:iEnd]
    end

    if sparseList == true
        return solverListSparseMat
    else
        return solverListMat
    end         
end


function runBenchmarkOdeSolvers(peTabModel::PeTabModel, 
                                pathFileSave::String,
                                sparseLinSolvers::Bool;
                                nTimesRepat::UInt=UInt(3), 
                                tolsCheck::Array{Float64, 1}=[1e-6, 1e-9, 1e-12])
    
    println("Working with model ", peTabModel.modelName)

    # Process PeTab files into type-stable Julia structs 
    experimentalConditionsFile, measurementDataFile, parameterDataFile, observablesDataFile = readDataFiles(peTabModel.dirModel, readObs=true)
    parameterData = processParameterData(parameterDataFile)
    simulationInfo = getSimulationInfo(measurementDataFile)
     
    # Set model parameter values to those in the PeTab parameter data ensuring correct value of constant parameters 
    setParamToFileValues!(peTabModel.paramMap, peTabModel.stateMap, parameterData)
 
    # The time-span 5e3 is overwritten when performing actual forward simulations 
    odeProb = ODEProblem(peTabModel.odeSystem, peTabModel.stateMap, (0.0, 5e3), peTabModel.paramMap, jac=true, sparse=sparseLinSolvers)
    odeProb = remake(odeProb, p = convert.(Float64, odeProb.p), u0 = convert.(Float64, odeProb.u0))
    odeProbHighAcc = ODEProblem(peTabModel.odeSystem, peTabModel.stateMap, (0.0, 5e3), peTabModel.paramMap, jac=true, sparse=false)
    odeProbHighAcc = remake(odeProbHighAcc, p = convert.(Float64, odeProb.p), u0 = convert.(Float64, odeProb.u0))

    # Functions to map experimental conditions and parameters correctly to the ODE model 
    changeToExperimentalCondUse! = (pVec, u0Vec, expID) -> changeExperimentalCond!(pVec, u0Vec, expID, parameterData, experimentalConditionsFile, peTabModel)

    # High accruacy ODE solution 
    println("Computing high accuracy model")
    solArrayHighAcc, statusHighAcc = calcHighAccOdeSolution(odeProbHighAcc, changeToExperimentalCondUse!, measurementDataFile, simulationInfo, tol=1e-14)
    println("Done")

    solverInfoMat = getSolverInfo(sparseLinSolvers)

    if statusHighAcc != true
        # Log failure to disk 
        println("High accuracy solver failed for $modelFile")
        open(joinpath(pwd(), "Benchmakrs", "ODE_solvers", "Log.txt"), "a") do io
            println(io, "Failed with high accuracy solution for $modelFile")
        end
        return 
    else
        println("Done with high accuracy solution")
    end

    for i in 1:size(solverInfoMat)[1]
        solver = solverInfoMat[i, 1]
        solverNames = solverInfoMat[i, 2]
        solverType = solverInfoMat[i, 3]
        solverLib = solverInfoMat[i, 4]

        println("Trying solver = ", solver)
        # Crashes as problem is to stiff 
        if !((peTabModel.modelName == "model_Crauste_CellSystems2017") && alg_solver == AutoTsit5(Rosenbrock23())) 
            for tol in tolsCheck
            
                benchRunTime = Vector{Float64}(undef, nTimesRepat)
                benchMemory = Vector{Float64}(undef, nTimesRepat)
                benchAllocs = Vector{Float64}(undef, nTimesRepat)

                # Check the accuracy of the ODE solver by comparing with high accuracy solution. In case the squared sum 
                # error cannot be computed the solver crashed and run time is not profiled.
                local sqDiffSolver = Float64
                try
                    sqDiffSolver = calcAccuracyOdeSolver(odeProb, solArrayHighAcc, changeToExperimentalCondUse!, measurementDataFile, simulationInfo, solver, tol::Float64)
                catch 
                    sqDiffSolver = Inf
                end
                println("sqDiffSolver = ", sqDiffSolver)
                solverSuccess = isinf(sqDiffSolver) ? false : true
                    
                if solverSuccess
                    for i in 1:nTimesRepat
                        b = @benchmark solveOdeModelAllExperimentalCond($odeProb, $changeToExperimentalCondUse!, $measurementDataFile, $simulationInfo, $solver, $tol) samples=1 evals=1
                        bMin = minimum(b)
                        benchRunTime[i] = bMin.time # microsecond
                        benchMemory[i] = bMin.memory # bytes
                        benchAllocs[i] = bMin.allocs # number of allocations
                    end
                else
                    benchRunTime .= NaN
                    benchMemory .= NaN
                    benchAllocs .= NaN
                end
                dataSave = DataFrame(model = peTabModel.modelName, 
                                     solver = solverNames, 
                                     solverType = solverType,
                                     solverLib = solverLib,
                                     nStates = length(peTabModel.stateNames)-1,
                                     nParam = length(peTabModel.paramNames),
                                     reltol = tol, 
                                     abstol = tol, 
                                     success = solverSuccess, 
                                     runTime = benchRunTime, 
                                     memory = benchMemory, 
                                     allocs = benchAllocs,
                                     sqDiff = sqDiffSolver, 
                                     iteration = 1:nTimesRepat)
                if isfile(pathFileSave)
                    CSV.write(pathFileSave, dataSave, append = true)
                else
                    CSV.write(pathFileSave, dataSave)
                end
            end
        else
            for tol in Tols
                dataSave = DataFrame(model = peTabModel.modelName, 
                                     solver = solverNames, 
                                     solverType = solverType,
                                     solverLib = solverLib,
                                     nStates = length(peTabModel.stateNames)-1,
                                     nParam = length(peTabModel.paramNames),
                                     reltol = tol, 
                                     abstol = tol, 
                                     success = false, 
                                     runTime = NaN,
                                     memory = NaN, 
                                     allocs = NaN,
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

    #=
    for alg_hint in alg_hints
        println("Alg_hint = ", alg_hint[1])
        if ~(modelFile == "model_Crauste_CellSystems2017.jl")
    =#

    GC.gc()
end


dirSave = pwd() * "/Intermediate/Benchmarks/ODE_solvers/"
pathFileSaveNotSparse = dirSave * "Sparse_not_linsolvers.csv"
pathFileSaveSparse = dirSave * "Sparse_linsolvers.csv"
if !isdir(dirSave)
    mkpath(dirSave)
end

modelListTry = ["model_Beer_MolBioSystems2014", "model_Weber_BMC2015", "model_Schwen_PONE2014", "model_Alkan_SciSignal2018", 
                "model_Bachmann_MSB2011", "model_Bertozzi_PNAS2020", "model_Blasi_CellSystems2016", "model_Boehm_JProteomeRes2014", 
                "model_Borghans_BiophysChem1997", "model_Brannmark_JBC2010", "model_Bruno_JExpBot2016", "model_Crauste_CellSystems2017", 
                "model_Elowitz_Nature2000", "model_Fiedler_BMC2016", "model_Fujita_SciSignal2010", "model_Giordano_Nature2020", 
                "model_Isensee_JCB2018", "model_Laske_PLOSComputBiol2019", "model_Lucarelli_CellSystems2018", "model_Okuonghae_ChaosSolitonsFractals2020", 
                "model_Oliveira_NatCommun2021", "model_Perelson_Science1996", "model_Rahman_MBS2016", 
                "model_SalazarCavazos_MBoC2020", "model_Sneyd_PNAS2002", "model_Zhao_QuantBiol2020", "model_Zheng_PNAS2012"]

for i in eachindex(modelListTry)
    modelName = modelListTry[i]
    dirModel = pwd() * "/Intermediate/PeTab_models/" * modelName * "/"
    peTabModel = setUpPeTabModel(modelName, dirModel)
    runBenchmarkOdeSolvers(peTabModel, pathFileSaveNotSparse, false, nTimesRepat=UInt(1))
    runBenchmarkOdeSolvers(peTabModel, pathFileSaveSparse, true, nTimesRepat=UInt(1))

end
