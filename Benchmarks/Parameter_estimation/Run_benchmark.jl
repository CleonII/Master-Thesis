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
using NLopt
using BenchmarkTools
using Zygote
using SciMLSensitivity


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
include(joinpath(pwd(), "src", "Optimizers", "Set_up_NLopt.jl"))
include(joinpath(pwd(), "src", "Optimizers", "Set_up_fides.jl"))

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
                                      algList=[:IpoptAutoHess, :IpoptBlockAutoDiff, :IpoptLBFGS, :OptimIPNewtonAutoHess, :OptimIPNewtonBlockAutoDiff, :OptimLBFGS, :NLoptLBFGS, :FidesAutoHess, :FidesBlockAutoHess])

    peTabOpt = setUpCostGradHess(peTabModel, solver, tol)

    pathCube = peTabModel.dirModel * "Cube_benchmark.csv"
    createCube(pathCube, peTabOpt, nStartGuess, seed=123, verbose=true)
    cube = Matrix(CSV.read(pathCube, DataFrame))

    dirRes = pwd() * "/Intermediate/Benchmarks/Parameter_estimation/" * peTabModel.modelName * "/"
    if !isdir(dirRes)
        mkpath(dirRes)
    end
    pathSave = dirRes * "Benchmark_result_fides.csv"

    # The termination criteria are set to match Optim which terminates based on;
    # abs(f - f_prev) ≤ f_tol*abs(f), f_tol = 1e-8
    # norm(x - x_prev) ≤ x_tot, x_tol = 0.0
    # norm(grad) ≤ gtol, gtol=1e-6
    # Optim dictates the termination criteria as it has the least flexible termination 
    # criteria. Ipopt terminates on entirely different basis (so hard to compare)

    # Ipopt optimizers 
    ipoptProbHessApprox, iterArrHessApprox = createIpoptProb(peTabOpt, hessianUse=:blockAutoDiff)
    ipoptProbBfgs, iterArrBfgs = createIpoptProb(peTabOpt, hessianUse=:LBFGS)
    ipoptProbAutoHess, iterArrAutoHess = createIpoptProb(peTabOpt, hessianUse=:autoDiff)
    Ipopt.AddIpoptNumOption(ipoptProbHessApprox, "acceptable_tol", 1e-8)
    Ipopt.AddIpoptNumOption(ipoptProbBfgs, "acceptable_tol", 1e-8)
    Ipopt.AddIpoptNumOption(ipoptProbAutoHess, "acceptable_tol", 1e-8)
    
    # Optim optimizers 
    optimProbHessApprox = createOptimProb(peTabOpt, IPNewton(), hessianUse=:blockAutoDiff, 
                                          options=Optim.Options(iterations = 1000, show_trace = false, allow_f_increases=true, 
                                                                successive_f_tol = 3, f_tol=1e-8, g_tol=1e-6, x_tol=0.0))
    optimProbAutoHess = createOptimProb(peTabOpt, IPNewton(), hessianUse=:autoDiff, 
                                        options=Optim.Options(iterations = 1000, show_trace = false, allow_f_increases=true, 
                                                              successive_f_tol = 3, f_tol=1e-8, g_tol=1e-6, x_tol=0.0))
    optimProbLBFGS = createOptimProb(peTabOpt, LBFGS())
    
    # NLopt optimizers 
    NLoptLBFGS = createNLoptProb(peTabOpt, :LD_LBFGS, verbose=false)
    
    # Fides 
    FidesAutoHess = setUpFides(peTabOpt, :autoDiff; verbose=0, 
                               options=py"{'maxiter' : 1000, 'fatol' : 0.0, 'frtol' : 1e-8, 'xtol' : 0.0, 'gatol' : 1e-6, 'grtol' : 0.0}"o)
    FidesAutoHessBlock = setUpFides(peTabOpt, :blockAutoDiff; verbose=0, 
                                    options=py"{'maxiter' : 1000, 'fatol' : 0.0, 'frtol' : 1e-8, 'xtol' : 0.0, 'gatol' : 1e-6, 'grtol' : 0.0}"o)

    # Make sure to activate allocation of required arrays for Ipopt solvers, and to compile 
    # gradient and required hessian functions to avoid bias in run-times for benchmark. 
    pTmp, gradTmp, hessTmp = cube[1, :], zeros(peTabOpt.nParamEst), zeros(peTabOpt.nParamEst, peTabOpt.nParamEst)
    costTmp = peTabOpt.evalF(pTmp)
    peTabOpt.evalGradF(gradTmp, pTmp)
    if :IpoptAutoHess in algList || :OptimIPNewtonAutoHess in algList || :FidesAutoHess in algList
        peTabOpt.evalHess(hessTmp, pTmp)
    end
    if :IpoptBlockAutoDiff in algList || :OptimIPNewtonBlockAutoDiff in algList || :FidesBlockAutoDiff in algList
        peTabOpt.evalHessApprox(hessTmp, pTmp)
    end

    for i in 1:nStartGuess

        p0 = cube[i, :]
        println("Iteration = $i")

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

        if :OptimIPNewtonAutoHess in algList
            res = optimProbAutoHess(p0, showTrace=false)
            writeFile(pathSave, res.minimum, res.time_run, res.f_converged, res.iterations, i, "optimIPNewtonAutoHess", solverStr, string(tol))
        end

        if :OptimIPNewtonBlockAutoDiff in algList
            res = optimProbHessApprox(p0, showTrace=false)
            writeFile(pathSave, res.minimum, res.time_run, res.f_converged, res.iterations, i, "optimIPNewtonBlockAutoHess", solverStr, string(tol))
        end

        if :OptimLBFGS in algList
            res = optimProbLBFGS(p0, showTrace=false)
            writeFile(pathSave, res.minimum, res.time_run, res.f_converged, res.iterations, i, "optimLBFGS", solverStr, string(tol))
        end

        if :NLoptLBFGS in algList
            runTime = @elapsed minF, min, ret = NLopt.optimize(NLoptLBFGS, p0)
            writeFile(pathSave, minF, runTime, string(ret), NLoptLBFGS.numevals, i, "NLoptLBFGS", solverStr, string(tol))
        end

        if :FidesAutoHess in algList
            runTime = @elapsed res, nIter, converged = FidesAutoHess(p0)
            writeFile(pathSave, res[1], runTime, string(converged), nIter, i, "FidesAutoHess", solverStr, string(tol))
        end

        if :FidesBlockAutoHess in algList
            runTime = @elapsed res, nIter, converged = FidesAutoHessBlock(p0)
            writeFile(pathSave, res[1], runTime, string(converged), nIter, i, "FidesAutoHessBlock", solverStr, string(tol))
        end
    end
end


if ARGS[1] == "Fiedler"
    loadFidesFromPython("/home/sebpe/anaconda3/envs/PeTab/bin/python")
    algsTest = [:IpoptAutoHess, :IpoptBlockAutoDiff, :IpoptLBFGS, :OptimIPNewtonAutoHess, :OptimIPNewtonBlockAutoDiff, :OptimLBFGS, :FidesAutoHess, :FidesBlockAutoHess]
    dirModel = pwd() * "/Intermediate/PeTab_models/model_Fiedler_BMC2016/"
    peTabModel = setUpPeTabModel("model_Fiedler_BMC2016", dirModel)
    benchmarkParameterEstimation(peTabModel, QNDF(), "QNDF", 1e-9, 1000, algList=algsTest) 
end

if ARGS[1] == "Boehm"
    loadFidesFromPython("/home/sebpe/anaconda3/envs/PeTab/bin/python")
    algsTest = [:IpoptAutoHess, :IpoptBlockAutoDiff, :IpoptLBFGS, :OptimIPNewtonAutoHess, :OptimIPNewtonBlockAutoDiff, :OptimLBFGS, :FidesAutoHess, :FidesBlockAutoHess]
    dirModel = pwd() * "/Intermediate/PeTab_models/model_Boehm_JProteomeRes2014/"
    peTabModel = setUpPeTabModel("model_Boehm_JProteomeRes2014", dirModel)
    benchmarkParameterEstimation(peTabModel, QNDF(), "QNDF", 1e-9, 1000, algList=algsTest) 
end

if ARGS[1] == "Bachmann"
    dirModel = pwd() * "/Intermediate/PeTab_models/model_Bachmann_MSB2011/"
    peTabModel = setUpPeTabModel("model_Bachmann_MSB2011", dirModel)
    algsTest = [:IpoptLBFGS, :NLoptLBFGS]
    benchmarkParameterEstimation(peTabModel, QNDF(), "QNDF", 1e-6, 1000, algList=algsTest)
end


if ARGS[1] == "Brannmark"
    dirModel = pwd() * "/Intermediate/PeTab_models/model_Brannmark_JBC2010/"
    peTabModel = setUpPeTabModel("model_Brannmark_JBC2010", dirModel)
    algsTest = [:IpoptAutoHess, :IpoptBlockAutoDiff, :IpoptLBFGS, :OptimIPNewtonAutoHess, :OptimIPNewtonBlockAutoDiff, :FidesAutoHess, :FidesBlockAutoHess]
    benchmarkParameterEstimation(peTabModel, Rodas5P(), "Rodas5P", 1e-6, 1000, algList=algsTest)
end


if ARGS[1] == "Fujita"
    dirModel = pwd() * "/Intermediate/PeTab_models/model_Fujita_SciSignal2010/"
    peTabModel = setUpPeTabModel("model_Fujita_SciSignal2010", dirModel)
    algsTest = [:IpoptAutoHess, :IpoptBlockAutoDiff, :IpoptLBFGS, :OptimIPNewtonAutoHess, :OptimIPNewtonBlockAutoDiff, :OptimLBFGS, :FidesAutoHess, :FidesBlockAutoHess]
    benchmarkParameterEstimation(peTabModel, Rodas5(), "Rodas5", 1e-6, 1000, algList=algsTest)
end


if ARGS[1] == "Crauste"
    dirModel = pwd() * "/Intermediate/PeTab_models/model_Crauste_CellSystems2017/"
    peTabModel = setUpPeTabModel("model_Crauste_CellSystems2017", dirModel)
    algsTest = [:IpoptAutoHess, :IpoptLBFGS, :OptimIPNewtonAutoHess, :NLoptLBFGS, :OptimLBFGS]
    benchmarkParameterEstimation(peTabModel, Rodas4P(), "Rodas4P", 1e-8, 1000, algList=algsTest)
end


if ARGS[1] == "Zheng"
    dirModel = pwd() * "/Intermediate/PeTab_models/model_Zheng_PNAS2012/"
    peTabModel = setUpPeTabModel("model_Zheng_PNAS2012", dirModel)
    algsTest = [:IpoptLBFGS, :IpoptBlockAutoDiff, :OptimIPNewtonBlockAutoDiff, :OptimLBFGS, :NLoptLBFGS]
    benchmarkParameterEstimation(peTabModel, QNDF(), "QNDF", 1e-6, 1000, algList=algsTest)
end