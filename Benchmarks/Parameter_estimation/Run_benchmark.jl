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
using BenchmarkTools


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

# Optimizers 
include(joinpath(pwd(), "src", "Optimizers", "Set_up_Ipopt.jl"))
include(joinpath(pwd(), "src", "Optimizers", "Set_up_optim.jl"))

function writeFile(pathFile::String, 
                   finalCost, 
                   runTime, 
                   retCode, 
                   nIter, 
                   startGuess, 
                   alg::String,
                   solver::String,
                   tol::String)

    # Save after each iteration (do not loose data)
    dataSave = [alg finalCost runTime retCode nIter startGuess solver tol]
    dataSave = DataFrame(dataSave, ["Alg", "Cost", "Run_time", "Ret_code", "N_iter", "Start_guess", "Solver", "tol"])
    shouldAppend = isfile(pathFile) ? true : false
    CSV.write(pathFile, dataSave, append=shouldAppend)

end


function benchmarkParameterEstimation(peTabModel::PeTabModel, 
                                      solver, 
                                      solverStr::String, 
                                      tol::Float64, 
                                      nStartGuess::Integer;
                                      algList=[:IpoptAutoHess, :IpoptBlockAutoDiff, :IpoptLBFGS, :OptimAutoHess, :OptimBlockAutoDiff])

    peTabOpt = setUpCostGradHess(peTabModel, solver, tol)

    pathCube = peTabModel.dirModel * "Cube_benchmark.csv"
    createCube(pathCube, peTabOpt, nStartGuess, seed=123, verbose=true)
    cube = Matrix(CSV.read(pathCube, DataFrame))

    dirRes = pwd() * "/Intermediate/Benchmarks/Parameter_estimation/" * peTabModel.modelName * "/"
    if !isdir(dirRes)
        mkpath(dirRes)
    end
    pathSave = dirRes * "Benchmark_result.csv"

    # Ipopt optimizers 
    ipoptProbHessApprox, iterArrHessApprox = createIpoptProb(peTabOpt, hessianUse=:blockAutoDiff)
    ipoptProbBfgs, iterArrBfgs = createIpoptProb(peTabOpt, hessianUse=:LBFGS)
    ipoptProbAutoHess, iterArrAutoHess = createIpoptProb(peTabOpt, hessianUse=:autoDiff)
    # Optim optimizers 
    optimProbHessApprox = createOptimInteriorNewton(peTabOpt, hessianUse=:blockAutoDiff)
    optimProbAutoHess = createOptimInteriorNewton(peTabOpt, hessianUse=:autoDiff)

    # Make sure to activate allocation of required arrays for Ipopt solvers, and to compile 
    # gradient and required hessian functions to avoid bias in run-times for benchmark. 
    pTmp, gradTmp, hessTmp = cube[1, :], zeros(peTabOpt.nParamEst), zeros(peTabOpt.nParamEst, peTabOpt.nParamEst)
    costTmp = peTabOpt.evalF(pTmp)
    peTabOpt.evalGradF(gradTmp, pTmp)
    if :IpoptAutoHess in algList || :OptimAutoHess in algList
        peTabOpt.evalHess(hessTmp, pTmp)
    end
    if :IpoptBlockAutoDiff in algList || :OptimBlockAutoDiff in algList
        peTabOpt.evalHessApprox(hessTmp, pTmp)
    end

    for i in 1:nStartGuess

        println("Start guess = $i")
        p0 = cube[i, :]

        if :IpoptAutoHess in algList
            ipoptProbAutoHess.x = deepcopy(p0)
            Ipopt.AddIpoptIntOption(ipoptProbAutoHess, "print_level", 0)
            runTime = @elapsed sol_opt = Ipopt.IpoptSolve(ipoptProbAutoHess)
            writeFile(pathSave, ipoptProbAutoHess.obj_val, runTime, ipoptProbAutoHess.status, iterArrAutoHess[1], i, "IpoptAutoHess", solverStr, string(tol))
        end

        if :IpoptBlockAutoDiff in algList
            ipoptProbHessApprox.x = deepcopy(p0)
            Ipopt.AddIpoptIntOption(ipoptProbHessApprox, "print_level", 0)
            runTime = @elapsed sol_opt = Ipopt.IpoptSolve(ipoptProbHessApprox)
            writeFile(pathSave, ipoptProbHessApprox.obj_val, runTime, ipoptProbHessApprox.status, iterArrHessApprox[1], i, "IpoptBlockAutoHess", solverStr, string(tol))
        end

        if :IpoptLBFGS in algList
            ipoptProbBfgs.x = deepcopy(p0)
            Ipopt.AddIpoptIntOption(ipoptProbBfgs, "print_level", 0)
            runTime = @elapsed sol_opt = Ipopt.IpoptSolve(ipoptProbBfgs)
            writeFile(pathSave, ipoptProbBfgs.obj_val, runTime, ipoptProbBfgs.status, iterArrBfgs[1], i, "IpoptLBFGS", solverStr, string(tol))
        end

        if :OptimAutoHess in algList
            res = optimProbAutoHess(p0, showTrace=false)
            writeFile(pathSave, res.minimum, res.time_run, res.f_converged, res.iterations, i, "optimAutoHess", solverStr, string(tol))
        end

        if :OptimBlockAutoDiff in algList
            res = optimProbHessApprox(p0, showTrace=false)
            writeFile(pathSave, res.minimum, res.time_run, res.f_converged, res.iterations, i, "optimBlockAutoHess", solverStr, string(tol))
        end
    end
end

#=
dirModel = pwd() * "/Intermediate/PeTab_models/model_Boehm_JProteomeRes2014/"
peTabModel = setUpPeTabModel("model_Boehm_JProteomeRes2014", dirModel)
benchmarkParameterEstimation(peTabModel, Rodas4P(), "Rodas4", 1e-9, 1000) 
=#

dirModel = pwd() * "/Intermediate/PeTab_models/model_Fiedler_BMC2016/"
peTabModel = setUpPeTabModel("model_Fiedler_BMC2016", dirModel)
benchmarkParameterEstimation(peTabModel, QNDF(), "QNDF", 1e-9, 1000) 
