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
using Ipopt
using Optim
using NLopt

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
                   absTol::String, 
                   relTol::String)

    # Save after each iteration (do not loose data)
    dataSave = [alg finalCost runTime retCode nIter startGuess solver absTol relTol]
    dataSave = DataFrame(dataSave, ["Alg", "Cost", "Run_time", "Ret_code", "N_iter", "Start_guess", "Solver", "absTol", "relTol"])
    shouldAppend = isfile(pathFile) ? true : false
    CSV.write(pathFile, dataSave, append=shouldAppend)

end


function benchmarkParameterEstimation(petabModel::PEtabModel, 
                                      solver, 
                                      solverStr::String, 
                                      absTol::Float64, 
                                      relTol::Float64,
                                      nStartGuess::Integer;
                                      algList=[:IpoptAutoHess, :IpoptBlockAutoDiff, :IpoptLBFGS, :OptimIPNewtonAutoHess, :OptimIPNewtonBlockAutoDiff, :OptimLBFGS, :NLoptLBFGS, :FidesAutoHess, :FidesBlockAutoHess, :FidesBFGS])

    petabProblem = setUpPEtabODEProblem(petabModel, solver, solverAbsTol=absTol, solverRelTol=relTol)

    pathCube = joinpath(petabModel.dirJulia, "Cube_benchmark.csv")
    createCube(pathCube, petabProblem, nStartGuess, seed=123, verbose=true)
    cube = Matrix(CSV.read(pathCube, DataFrame))

    dirResult = joinpath(@__DIR__, "..", "..", "Intermediate", "Benchmarks", "Parameter_estimation", petabModel.modelName)
    if !isdir(dirResult)
        mkpath(dirResult)
    end
    pathSave = joinpath(dirResult, "Benchmark_result_estimation.csv")

    #=
        The termination criteria are set to match Optim which terminates based on;
        abs(f - f_prev) ≤ f_tol*abs(f), f_tol = 1e-8
        norm(x - x_prev) ≤ x_tot, x_tol = 0.0
        norm(grad) ≤ gtol, gtol=1e-6
        Optim dictates the termination criteria as it has the least flexible termination 
        criteria. Ipopt terminates on entirely different basis (so hard to compare)
    =#

    # Ipopt optimizers 
    ipoptProbHessApprox, iterArrHessApprox = createIpoptProb(petabProblem, hessianUse=:blockAutoDiff)
    ipoptProbBfgs, iterArrBfgs = createIpoptProb(petabProblem, hessianUse=:LBFGS)
    ipoptProbAutoHess, iterArrAutoHess = createIpoptProb(petabProblem, hessianUse=:autoDiff)
    Ipopt.AddIpoptNumOption(ipoptProbHessApprox, "acceptable_tol", 1e-8)
    Ipopt.AddIpoptNumOption(ipoptProbBfgs, "acceptable_tol", 1e-8)
    Ipopt.AddIpoptNumOption(ipoptProbAutoHess, "acceptable_tol", 1e-8)
    
    # Optim optimizers 
    optimProbHessApprox = createOptimProb(petabProblem, IPNewton(), hessianUse=:blockAutoDiff, 
                                          options=Optim.Options(iterations = 1000, show_trace = false, allow_f_increases=true, 
                                                                successive_f_tol = 3, f_tol=1e-8, g_tol=1e-6, x_tol=0.0))
    optimProbAutoHess = createOptimProb(petabProblem, IPNewton(), hessianUse=:autoDiff, 
                                        options=Optim.Options(iterations = 1000, show_trace = false, allow_f_increases=true, 
                                                              successive_f_tol = 3, f_tol=1e-8, g_tol=1e-6, x_tol=0.0))
    optimProbLBFGS = createOptimProb(petabProblem, LBFGS(), 
                                     options=Optim.Options(iterations = 250, 
                                                           show_trace = false, 
                                                           allow_f_increases=true, 
                                                           outer_iterations = 4, 
                                                           successive_f_tol = 3, f_tol=1e-8, g_tol=1e-6, x_tol=0.0))
    
    # NLopt optimizers 
    NLoptLBFGS = createNLoptProb(petabProblem, :LD_LBFGS, verbose=false)
    NLoptLBFGS.ftol_rel = 1e-8
    NLoptLBFGS.xtol_rel = 0.0
    NLoptLBFGS.maxeval = 5000
    
    # Fides 
    FidesAutoHess = setUpFides(petabProblem, :autoDiff; verbose=0, 
                               options=py"{'maxiter' : 1000, 'fatol' : 0.0, 'frtol' : 1e-8, 'xtol' : 0.0, 'gatol' : 1e-6, 'grtol' : 0.0}"o)
    FidesAutoHessBlock = setUpFides(petabProblem, :blockAutoDiff; verbose=0, 
                                    options=py"{'maxiter' : 1000, 'fatol' : 0.0, 'frtol' : 1e-8, 'xtol' : 0.0, 'gatol' : 1e-6, 'grtol' : 0.0}"o)
    FidesBFGS = setUpFides(petabProblem, :None; verbose=0,
                           fidesHessApprox=py"fides.hessian_approximation.BFGS()"o, 
                           options=py"{'maxiter' : 1000, 'fatol' : 0.0, 'frtol' : 1e-8, 'xtol' : 0.0, 'gatol' : 1e-6, 'grtol' : 0.0}"o)

    # Make sure to activate allocation of required arrays for Ipopt solvers, and to compile 
    # gradient and required hessian functions to prevent pre-allocations from affecting run-times 
    # The hessian is only pre-compilied in case of a Hessian based algorithm being used 
    θ_tmp = cube[1, :]
    _gradient = zeros(length(θ_tmp))
    _hessian = zeros(length(θ_tmp), length(θ_tmp))
    _cost = petabProblem.computeCost(θ_tmp)
    petabProblem.computeGradientAutoDiff(_gradient, θ_tmp)
    if :IpoptAutoHess in algList || :OptimIPNewtonAutoHess in algList || :FidesAutoHess in algList
        petabProblem.computeHessian(_hessian, θ_tmp)
    end
    if :IpoptBlockAutoDiff in algList || :OptimIPNewtonBlockAutoDiff in algList || :FidesBlockAutoDiff in algList
        petabProblem.computeHessianBlock(_hessian, θ_tmp)
    end

    for i in 1:nStartGuess

        p0 = cube[i, :]

        if i % 10 == 0 || i == 1
            println("Iteration = $i")
        end

        if :IpoptAutoHess in algList
            ipoptProbAutoHess.x = deepcopy(p0)
            Ipopt.AddIpoptIntOption(ipoptProbAutoHess, "print_level", 0)
            runTime = @elapsed sol_opt = Ipopt.IpoptSolve(ipoptProbAutoHess)
            writeFile(pathSave, ipoptProbAutoHess.obj_val, runTime, ipoptProbAutoHess.status, iterArrAutoHess[1], i, "IpoptAutoHess", solverStr, string(absTol), string(relTol))
        end

        if :IpoptBlockAutoDiff in algList
            ipoptProbHessApprox.x = deepcopy(p0)
            Ipopt.AddIpoptIntOption(ipoptProbHessApprox, "print_level", 0)
            runTime = @elapsed sol_opt = Ipopt.IpoptSolve(ipoptProbHessApprox)
            writeFile(pathSave, ipoptProbHessApprox.obj_val, runTime, ipoptProbHessApprox.status, iterArrHessApprox[1], i, "IpoptBlockAutoHess", solverStr, string(absTol), string(relTol))
        end

        if :IpoptLBFGS in algList
            ipoptProbBfgs.x = deepcopy(p0)
            Ipopt.AddIpoptIntOption(ipoptProbBfgs, "print_level", 0)
            try
                runTime = @elapsed sol_opt = Ipopt.IpoptSolve(ipoptProbBfgs)
                writeFile(pathSave, ipoptProbBfgs.obj_val, runTime, ipoptProbBfgs.status, iterArrBfgs[1], i, "IpoptLBFGS", solverStr, string(absTol), string(relTol))
            catch
                writeFile(pathSave, Inf, Inf, 0, Inf, i, "IpoptLBFGS", solverStr, string(absTol), string(relTol))
            end
        end

        if :OptimIPNewtonAutoHess in algList
            res = optimProbAutoHess(p0, showTrace=false)
            writeFile(pathSave, res.minimum, res.time_run, res.f_converged, res.iterations, i, "optimIPNewtonAutoHess", solverStr, string(absTol), string(relTol))
        end

        if :OptimIPNewtonBlockAutoDiff in algList
            res = optimProbHessApprox(p0, showTrace=false)
            writeFile(pathSave, res.minimum, res.time_run, res.f_converged, res.iterations, i, "optimIPNewtonBlockAutoHess", solverStr, string(absTol), string(relTol))
        end

        if :OptimLBFGS in algList
            try
                res = optimProbLBFGS(p0, showTrace=false)
                println("Final cost value = ", petabProblem.evalF(res.minimizer))
                writeFile(pathSave, res.minimum, res.time_run, res.f_converged, res.iterations, i, "optimLBFGS", solverStr, string(absTol), string(relTol))
            catch
                writeFile(pathSave, Inf, Inf, 0, Inf, i, "optimLBFGS", solverStr, string(absTol), string(relTol))
            end
        end

        if :NLoptLBFGS in algList
            runTime = @elapsed minF, min, ret = NLopt.optimize(NLoptLBFGS, p0)
            writeFile(pathSave, minF, runTime, string(ret), NLoptLBFGS.numevals, i, "NLoptLBFGS", solverStr, string(absTol), string(relTol))
        end

        if :FidesAutoHess in algList
            try
                runTime = @elapsed res, nIter, converged = FidesAutoHess(p0)
                writeFile(pathSave, res[1], runTime, string(converged), nIter, i, "FidesAutoHess", solverStr, string(absTol), string(relTol))
            catch
                writeFile(pathSave, Inf, Inf, 0, Inf, i, "FidesAutoHess", solverStr, string(absTol), string(relTol))
            end
        end

        if :FidesBlockAutoHess in algList
            try
                runTime = @elapsed res, nIter, converged = FidesAutoHessBlock(p0)
                writeFile(pathSave, res[1], runTime, string(converged), nIter, i, "FidesAutoHessBlock", solverStr, string(absTol), string(relTol))
            catch
                writeFile(pathSave, Inf, Inf, 0, Inf, i, "FidesAutoHessBlock", solverStr, string(absTol), string(relTol))
            end
        end

        if :FidesBFGS in algList
            try
                runTime = @elapsed res, nIter, converged = FidesBFGS(p0)
                writeFile(pathSave, res[1], runTime, string(converged), nIter, i, "FidesBFGS", solverStr, string(absTol), string(relTol))
            catch
                writeFile(pathSave, Inf, Inf, 0, Inf, i, "FidesBFGS", solverStr, string(absTol), string(relTol))
            end
        end
    end
end


loadFidesFromPython("/home/sebpe/anaconda3/envs/PeTab/bin/python")


if ARGS[1] == "Fiedler"
    algsTest = [:IpoptAutoHess, :IpoptBlockAutoDiff, :IpoptLBFGS, :OptimIPNewtonAutoHess, :OptimIPNewtonBlockAutoDiff, :FidesAutoHess, :FidesBlockAutoHess]
    pathYML = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", "model_Fiedler_BMC2016", "Fiedler_BMC2016.yaml")
    petabModel = readPEtabModel(pathYML, verbose=true)
    benchmarkParameterEstimation(petabModel, QNDF(), "QNDF", 1e-8, 1e-8, 1000, algList=algsTest) 
end


if ARGS[1] == "Boehm"
    algsTest = [:IpoptAutoHess, :IpoptBlockAutoDiff, :OptimIPNewtonAutoHess, :OptimIPNewtonBlockAutoDiff, :FidesAutoHess, :FidesBlockAutoHess, :FidesBFGS]
    pathYML = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", "model_Boehm_JProteomeRes2014", "Boehm_JProteomeRes2014.yaml")
    petabModel = readPEtabModel(pathYML, verbose=true)     
    benchmarkParameterEstimation(petabModel, QNDF(), "QNDF", 1e-8, 1e-8, 1000, algList=algsTest) 
end


if ARGS[1] == "Elowitz"
    algsTest = [:OptimIPNewtonAutoHess, :OptimIPNewtonBlockAutoDiff, :NLoptLBFGS, :FidesAutoHess, :FidesBlockAutoHess]
    pathYML = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", "model_Elowitz_Nature2000", "Elowitz_Nature2000.yaml")
    petabModel = readPEtabModel(pathYML, verbose=true)
    benchmarkParameterEstimation(petabModel, QNDF(), "QNDF", 1e-8, 1e-8, 1000, algList=algsTest) 
end


if ARGS[1] == "Crauste"
    algsTest = [:OptimLBFGS, :NLoptLBFGS, :OptimIPNewtonAutoHess, :FidesAutoHess]
    pathYML = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", "model_Crauste_CellSystems2017", "Crauste_CellSystems2017.yaml")
    petabModel = readPEtabModel(pathYML, verbose=true)
    benchmarkParameterEstimation(petabModel, AutoVern7(Rodas5()), "Vern7(Rodas5)", 1e-12, 1e-12, 1000, algList=algsTest) 
end


if ARGS[1] == "Weber"
    algsTest = [:IpoptLBFGS, :OptimLBFGS, :OptimIPNewtonAutoHess, :FidesAutoHess]
    pathYML = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", "model_Weber_BMC2015", "Weber_BMC2015.yaml")
    petabModel = readPEtabModel(pathYML, verbose=true)
    benchmarkParameterEstimation(petabModel, Rodas5(), "Rodas5", 1e-8, 1e-8, 1000, algList=algsTest) 
end


if ARGS[1] == "Bachmann"
    pathYML = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", "model_Bachmann_MSB2011", "Bachmann_MSB2011.yaml")
    petabModel = readPEtabModel(pathYML, verbose=true)
    algsTest = [:IpoptLBFGS, :OptimLBFGS, :NLoptLBFGS]
    benchmarkParameterEstimation(petabModel, Rodas4P(), "QNDF", 1e-8, 1e-8, 1000, algList=algsTest)
end


if ARGS[1] == "Brannmark"
    pathYML = joinpath(@__DIR__, "..", "Intermediate", "PeTab_models", "model_Brannmark_JBC2010", "Brannmark_JBC2010.yaml")
    petabModel = readPEtabModel(pathYML, verbose=true)
    algsTest = [:IpoptAutoHess, :IpoptBlockAutoDiff, :IpoptLBFGS, :OptimIPNewtonAutoHess, :OptimIPNewtonBlockAutoDiff, :FidesAutoHess, :FidesBlockAutoHess]
    benchmarkParameterEstimation(petabModel, Rodas5P(), "Rodas5P", 1e-8, 1e-8, 1000, algList=algsTest)
end


if ARGS[1] == "Fujita"
    pathYML = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", "model_Fujita_SciSignal2010", "Fujita_SciSignal2010.yaml")
    petabModel = readPEtabModel(pathYML, verbose=true)
    algsTest = [:IpoptAutoHess, :IpoptBlockAutoDiff, :IpoptLBFGS, :OptimIPNewtonAutoHess, :OptimIPNewtonBlockAutoDiff, :OptimLBFGS, :FidesAutoHess, :FidesBlockAutoHess]
    benchmarkParameterEstimation(petabModel, Rodas5(), "Rodas5", 1e-8, 1e-8, 1000, algList=algsTest)
end


if ARGS[1] == "Zheng"
    pathYML = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", "model_Zheng_PNAS2012", "Zheng_PNAS2012.yaml")
    petabModel = readPEtabModel(pathYML, verbose=true)
    algsTest = [:IpoptLBFGS, :OptimIPNewtonAutoHess, :FidesAutoHess, :NLoptLBFGS]
    benchmarkParameterEstimation(petabModel, QNDF(), "QNDF", 1e-8, 1e-8, 1000, algList=algsTest)
end


if ARGS[1] == "Schwen"
    pathYML = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", "model_Schwen_PONE2014", "Schwen_PONE2014.yaml")
    petabModel = readPEtabModel(pathYML, verbose=true)
    algsTest = [:IpoptBlockAutoDiff, :IpoptLBFGS, :OptimIPNewtonBlockAutoDiff, :OptimLBFGS, :FidesBlockAutoHess]
    benchmarkParameterEstimation(petabModel, QNDF(), "QNDF", 1e-8, 1e-8, 1000, algList=algsTest)
end


if ARGS[1] == "Sneyd"
    pathYML = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", "model_Sneyd_PNAS2002", "Sneyd_PNAS2002.yaml")
    petabModel = readPEtabModel(pathYML, verbose=true)
    algsTest = [:IpoptAutoHess, :IpoptLBFGS, :OptimIPNewtonAutoHess, :OptimLBFGS, :FidesAutoHess]
    benchmarkParameterEstimation(petabModel, QNDF(), "QNDF", 1e-8, 1e-8, 1000, algList=algsTest)
end
