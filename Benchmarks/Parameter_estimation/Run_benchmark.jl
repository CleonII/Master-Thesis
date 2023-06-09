using ModelingToolkit 
using OrdinaryDiffEq
using DiffEqCallbacks
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


function writeFile(dirSave::String,
                   θ_opt::Vector{Float64},
                   parameterNames::Vector{String}, 
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
    pathFile = joinpath(dirSave, "Estimation_statistics.csv")
    dataSave = [alg finalCost runTime retCode nIter startGuess solver absTol relTol]
    dataSave = DataFrame(dataSave, ["Alg", "Cost", "Run_time", "Ret_code", "N_iter", "Start_guess", "Solver", "absTol", "relTol"])
    shouldAppend = isfile(pathFile) ? true : false
    CSV.write(pathFile, dataSave, append=shouldAppend)

    # Save optimal parameter vector
    pathFile = joinpath(dirSave, "Minimizer.csv")
    _dataSave = Matrix{Any}(undef, (1, length(θ_opt)+5))
    _dataSave[:] .= vcat(θ_opt, startGuess, alg, solver, absTol, relTol)
    dataSaveθ = DataFrame(_dataSave, vcat(parameterNames, "Start_guess", "Alg", "Solver", "absTol", "relTol"))
    shouldAppend = isfile(pathFile) ? true : false
    CSV.write(pathFile, dataSaveθ, append=shouldAppend)
    
end


function benchmarkParameterEstimation(petabModel::PEtabModel, 
                                      solver, 
                                      solverStr::String, 
                                      absTol::Float64, 
                                      relTol::Float64,
                                      nStartGuess::Integer;
                                      algList=[:IpoptAutoHess, :IpoptBlockAutoDiff, :IpoptLBFGS, :OptimIPNewtonAutoHess, :OptimIPNewtonBlockAutoDiff, :OptimLBFGS, :NLoptLBFGS, :FidesAutoHess, :FidesBlockAutoHess, :FidesBFGS], 
                                      terminateSSMethod=:Norm, 
                                      solverSSRelTol::Float64=1e-6,
                                      solverSSAbsTol::Float64=1e-6, 
                                      reuseS::Bool=true, 
                                      numberOfprocesses=1,
                                      splitOverConditions::Bool=false)

    petabProblem = setUpPEtabODEProblem(petabModel, solver, solverAbsTol=absTol, solverRelTol=relTol, terminateSSMethod=terminateSSMethod, 
                                        solverSSRelTol=solverSSRelTol, solverSSAbsTol=solverSSAbsTol, 
                                        reuseS=reuseS, sensealgForwardEquations=:AutoDiffForward, odeSolverForwardEquations=solver, 
                                        splitOverConditions=splitOverConditions, numberOfprocesses=numberOfprocesses)
    θ_estNames = string.(petabProblem.θ_estNames)

    pathCube = joinpath(petabModel.dirJulia, "Cube_benchmark.csv")
    createCube(pathCube, petabProblem, nStartGuess, seed=123, verbose=true)
    cube = Matrix(CSV.read(pathCube, DataFrame))

    dirResult = joinpath(@__DIR__, "..", "..", "Intermediate", "Benchmarks", "Parameter_estimation", petabModel.modelName)
    if !isdir(dirResult)
        mkpath(dirResult)
    end

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
    optimProbGN = createOptimProb(petabProblem, IPNewton(), hessianUse=:GaussNewton, 
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
    FidesGN = setUpFides(petabProblem, :GaussNewton; verbose=0, 
                         options=py"{'maxiter' : 1000, 'fatol' : 0.0, 'frtol' : 1e-8, 'xtol' : 0.0, 'gatol' : 1e-6, 'grtol' : 0.0}"o)                                    
    FidesBFGS = setUpFides(petabProblem, :None; verbose=0,
                           fidesHessApprox=py"fides.hessian_approximation.BFGS()"o, 
                           options=py"{'maxiter' : 1000, 'fatol' : 0.0, 'frtol' : 1e-8, 'xtol' : 0.0, 'gatol' : 1e-6, 'grtol' : 0.0}"o)

    # Make sure to activate allocation of required arrays for Ipopt solvers, and to precompile 
    # gradient and required hessian functions to prevent pre-allocations from affecting run-times 
    # The hessian is only pre-compilied in case of a Hessian based algorithm being used 
    θ_tmp = cube[1, :]
    _gradient = zeros(length(θ_tmp))
    _hessian = zeros(length(θ_tmp), length(θ_tmp))
    print("Precompiling cost ... ")
    _cost = petabProblem.computeCost(θ_tmp)
    print("done \n")
    print("Precompiling gradient ... ")
    petabProblem.computeGradientAutoDiff(_gradient, θ_tmp)
    print("done \n")
    GC.gc(); GC.gc(); GC.gc()
    if :IpoptAutoHess in algList || :OptimIPNewtonAutoHess in algList || :FidesAutoHess in algList
        print("Precompiling hessian ... ")
        petabProblem.computeHessian(_hessian, θ_tmp)
        print("done \n")
    end
    if :IpoptBlockAutoDiff in algList || :OptimIPNewtonBlockAutoDiff in algList || :FidesBlockAutoDiff in algList
        petabProblem.computeHessianBlock(_hessian, θ_tmp)
    end
    if :FidesGN in algList || :OptimIPNewtonGN in algList 
        print("Precompiling Gauss Newton ... ")
        petabProblem.computeGradientForwardEquations(_gradient, θ_tmp)
        petabProblem.computeHessianGN(_hessian, θ_tmp)
        print("done \n")
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
            writeFile(dirResult, ipoptProbAutoHess.x, θ_estNames, ipoptProbAutoHess.obj_val, runTime, ipoptProbAutoHess.status, iterArrAutoHess[1], i, "IpoptAutoHess", solverStr, string(absTol), string(relTol))
        end

        if :IpoptBlockAutoDiff in algList
            ipoptProbHessApprox.x = deepcopy(p0)
            Ipopt.AddIpoptIntOption(ipoptProbHessApprox, "print_level", 0)
            runTime = @elapsed sol_opt = Ipopt.IpoptSolve(ipoptProbHessApprox)
            writeFile(dirResult, ipoptProbHessApprox.x, θ_estNames, ipoptProbHessApprox.obj_val, runTime, ipoptProbHessApprox.status, iterArrHessApprox[1], i, "IpoptBlockAutoHess", solverStr, string(absTol), string(relTol))
        end

        if :IpoptLBFGS in algList
            ipoptProbBfgs.x = deepcopy(p0)
            Ipopt.AddIpoptIntOption(ipoptProbBfgs, "print_level", 0)
            try
                runTime = @elapsed sol_opt = Ipopt.IpoptSolve(ipoptProbBfgs)
                writeFile(dirResult, ipoptProbBfgs.x, θ_estNames, ipoptProbBfgs.obj_val, runTime, ipoptProbBfgs.status, iterArrBfgs[1], i, "IpoptLBFGS", solverStr, string(absTol), string(relTol))
            catch
                writeFile(dirResult, ipoptProbBfgs.x, θ_estNames, Inf, Inf, 0, Inf, i, "IpoptLBFGS", solverStr, string(absTol), string(relTol))
            end
        end

        if :OptimIPNewtonAutoHess in algList
            res = optimProbAutoHess(p0, showTrace=false)
            writeFile(dirResult, res.minimizer, θ_estNames, res.minimum, res.time_run, res.f_converged, res.iterations, i, "optimIPNewtonAutoHess", solverStr, string(absTol), string(relTol))
        end

        if :OptimIPNewtonBlockAutoDiff in algList
            res = optimProbHessApprox(p0, showTrace=false)
            writeFile(dirResult, res.minimizer, θ_estNames, res.minimum, res.time_run, res.f_converged, res.iterations, i, "optimIPNewtonBlockAutoHess", solverStr, string(absTol), string(relTol))
        end

        if :OptimIPNewtonGN in algList
            res = optimProbGN(p0, showTrace=false)
            writeFile(dirResult, res.minimizer, θ_estNames, res.minimum, res.time_run, res.f_converged, res.iterations, i, "OptimIPNewtonGN", solverStr, string(absTol), string(relTol))
        end

        if :OptimLBFGS in algList
            try
                res = optimProbLBFGS(p0, showTrace=false)
                println("Final cost value = ", petabProblem.evalF(res.minimizer))
                writeFile(dirResult, res.minimizer, θ_estNames, res.minimum, res.time_run, res.f_converged, res.iterations, i, "optimLBFGS", solverStr, string(absTol), string(relTol))
            catch
                writeFile(dirResult, p0, θ_estNames, Inf, Inf, 0, Inf, i, "optimLBFGS", solverStr, string(absTol), string(relTol))
            end
        end

        if :NLoptLBFGS in algList
            runTime = @elapsed minF, min, ret = NLopt.optimize(NLoptLBFGS, p0)
            writeFile(dirResult, min, θ_estNames, minF, runTime, string(ret), NLoptLBFGS.numevals, i, "NLoptLBFGS", solverStr, string(absTol), string(relTol))
        end

        if :FidesAutoHess in algList
            try
                runTime = @elapsed res, nIter, converged = FidesAutoHess(p0)
                writeFile(dirResult, res[2], θ_estNames, res[1], runTime, string(converged), nIter, i, "FidesAutoHess", solverStr, string(absTol), string(relTol))
            catch
                writeFile(dirResult, p0, θ_estNames, Inf, Inf, 0, Inf, i, "FidesAutoHess", solverStr, string(absTol), string(relTol))
            end
        end

        if :FidesBlockAutoHess in algList
            try
                runTime = @elapsed res, nIter, converged = FidesAutoHessBlock(p0)
                writeFile(dirResult, res[2], θ_estNames, res[1], runTime, string(converged), nIter, i, "FidesAutoHessBlock", solverStr, string(absTol), string(relTol))
            catch
                writeFile(dirResult, p0, θ_estNames, Inf, Inf, 0, Inf, i, "FidesAutoHessBlock", solverStr, string(absTol), string(relTol))
            end
        end

        if :FidesBFGS in algList
            try
                runTime = @elapsed res, nIter, converged = FidesBFGS(p0)
                writeFile(dirResult, res[2], θ_estNames, res[1], runTime, string(converged), nIter, i, "FidesBFGS", solverStr, string(absTol), string(relTol))
            catch
                writeFile(dirResult, p0, θ_estNames, Inf, Inf, 0, Inf, i, "FidesBFGS", solverStr, string(absTol), string(relTol))
            end
        end

        if :FidesGN in algList
            try
                runTime = @elapsed res, nIter, converged = FidesGN(p0)
                writeFile(dirResult, res[2], θ_estNames, res[1], runTime, string(converged), nIter, i, "FidesGN", solverStr, string(absTol), string(relTol))
            catch
                writeFile(dirResult, p0, θ_estNames, Inf, Inf, 0, Inf, i, "FidesGN", solverStr, string(absTol), string(relTol))
            end
        end
    end
end


loadFidesFromPython("/home/sebpe/anaconda3/envs/PeTab/bin/python")


if length(ARGS) < 3
    println("Error : Must be provide at least three command line arguments")
    exit(1)
end


absTol, relTol = 1e-8, 1e-8
modelRun = ARGS[1]
nMultiStarts = parse(Int64, ARGS[2])
optmizersTest = Symbol.(ARGS[3:end])

#=
    With Fides we can reuse the sensitivity matrix when computing the GN hessian approxmiation.
    However, for IPNewton in Optim this does not work well. Hence we run a seperate run for 
    IPNewton with GN in case it is one of the provided algorithms.
=#
iOptimIPNewtonGN = findall(x -> x == :OptimIPNewtonGN, optmizersTest)
iNotOptimIPNewtonGN = findall(x -> x != :OptimIPNewtonGN, optmizersTest)

if ARGS[1] == "Fiedler_BMC2016"
    pathYML = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", "model_Fiedler_BMC2016", "Fiedler_BMC2016.yaml")
    petabModel = readPEtabModel(pathYML, verbose=true)
    benchmarkParameterEstimation(petabModel, QNDF(), "QNDF", absTol, relTol, nMultiStarts, algList=optmizersTest[iNotOptimIPNewtonGN]) 
    benchmarkParameterEstimation(petabModel, QNDF(), "QNDF", absTol, relTol, nMultiStarts, algList=optmizersTest[iOptimIPNewtonGN], reuseS=false) 
end


if ARGS[1] == "Beer_MolBioSystems2014"
    pathYML = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", "model_Beer_MolBioSystems2014", "Beer_MolBioSystems2014.yaml")
    petabModel = readPEtabModel(pathYML, verbose=true)
    benchmarkParameterEstimation(petabModel, Rodas5P(), "Rodas5P", absTol, relTol, nMultiStarts, algList=optmizersTest[iNotOptimIPNewtonGN], splitOverConditions=true) 
    benchmarkParameterEstimation(petabModel, Rodas5P(), "Rodas5P", absTol, relTol, nMultiStarts, algList=optmizersTest[iOptimIPNewtonGN], reuseS=false, splitOverConditions=true) 
end


if ARGS[1] == "Boehm_JProteomeRes2014"
    pathYML = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", "model_Boehm_JProteomeRes2014", "Boehm_JProteomeRes2014.yaml")
    petabModel = readPEtabModel(pathYML, verbose=true)     
    benchmarkParameterEstimation(petabModel, QNDF(), "QNDF", absTol, relTol, nMultiStarts, algList=optmizersTest[iNotOptimIPNewtonGN]) 
    benchmarkParameterEstimation(petabModel, QNDF(), "QNDF", absTol, relTol, nMultiStarts, algList=optmizersTest[iOptimIPNewtonGN], reuseS=false) 
end


if ARGS[1] == "Zheng_PNAS2012"
    pathYML = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", "model_Zheng_PNAS2012", "Zheng_PNAS2012.yaml")
    petabModel = readPEtabModel(pathYML, verbose=true)     
    benchmarkParameterEstimation(petabModel, Rodas5P(), "Rodas5P", absTol, relTol, nMultiStarts, algList=optmizersTest[iNotOptimIPNewtonGN]) 
    benchmarkParameterEstimation(petabModel, Rodas5P(), "Rodas5P", absTol, relTol, nMultiStarts, algList=optmizersTest[iOptimIPNewtonGN], reuseS=false) 
end


if ARGS[1] == "Sneyd_PNAS2002"
    pathYML = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", "model_Sneyd_PNAS2002", "Sneyd_PNAS2002.yaml")
    petabModel = readPEtabModel(pathYML, verbose=true)
    benchmarkParameterEstimation(petabModel, Rodas5P(), "Rodas5P", absTol, relTol, nMultiStarts, algList=optmizersTest[iNotOptimIPNewtonGN]) 
    benchmarkParameterEstimation(petabModel, Rodas5P(), "Rodas5P", absTol, relTol, nMultiStarts, algList=optmizersTest[iOptimIPNewtonGN], reuseS=false) 
end


if ARGS[1] == "Crauste_CellSystems2017"
    pathYML = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", "model_Crauste_CellSystems2017", "Crauste_CellSystems2017.yaml")
    petabModel = readPEtabModel(pathYML, verbose=true)
    benchmarkParameterEstimation(petabModel, Rodas5P(), "Rodas5P", absTol, relTol, nMultiStarts, algList=optmizersTest[iNotOptimIPNewtonGN]) 
    benchmarkParameterEstimation(petabModel, Rodas5P(), "Rodas5P", absTol, relTol, nMultiStarts, algList=optmizersTest[iOptimIPNewtonGN], reuseS=false) 
end


if ARGS[1] == "Lucarelli_CellSystems2018"
    pathYML = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", "model_Lucarelli_CellSystems2018", "Lucarelli_CellSystems2018.yaml")
    petabModel = readPEtabModel(pathYML, verbose=true)
    removeAllProcs()
    addprocs(1, exeflags="--project=.")
    benchmarkParameterEstimation(petabModel, QNDF(), "QNDF", absTol, relTol, nMultiStarts, algList=optmizersTest[iNotOptimIPNewtonGN], numberOfprocesses=2) 
    benchmarkParameterEstimation(petabgModel, QNDF(), "QNDF", absTol, relTol, nMultiStarts, algList=optmizersTest[iOptimIPNewtonGN], reuseS=false, numberOfprocesses=2) 
end


if ARGS[1] == "Weber_BMC2015"
    pathYML = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", "model_Weber_BMC2015", "Weber_BMC2015.yaml")
    petabModel = readPEtabModel(pathYML, verbose=true)
    benchmarkParameterEstimation(petabModel, Rodas5P(), "Rodas5P", absTol, relTol, nMultiStarts, algList=optmizersTest[iNotOptimIPNewtonGN]) 
    benchmarkParameterEstimation(petabModel, Rodas5P(), "Rodas5P", absTol, relTol, nMultiStarts, algList=optmizersTest[iOptimIPNewtonGN], reuseS=false) 
end


if ARGS[1] == "Bachmann_MSB2011"
    pathYML = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", "model_Bachmann_MSB2011", "Bachmann_MSB2011.yaml")
    petabModel = readPEtabModel(pathYML, verbose=true)
    benchmarkParameterEstimation(petabModel, QNDF(), "QNDF", absTol, relTol, nMultiStarts, algList=optmizersTest[iNotOptimIPNewtonGN], reuseS=true) 
    benchmarkParameterEstimation(petabModel, QNDF(), "QNDF", absTol, relTol, nMultiStarts, algList=optmizersTest[iOptimIPNewtonGN], reuseS=false) 
end


if ARGS[1] == "Brannmark_JBC2010"
    pathYML = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", "model_Brannmark_JBC2010", "Brannmark_JBC2010.yaml")
    petabModel = readPEtabModel(pathYML, verbose=true)
    benchmarkParameterEstimation(petabModel, Rodas5P(), "Rodas5P", absTol, relTol, nMultiStarts, algList=optmizersTest[iNotOptimIPNewtonGN], terminateSSMethod=:NewtonNorm, solverSSRelTol=1e-6, solverSSAbsTol=1e-6, reuseS=true)
    benchmarkParameterEstimation(petabModel, Rodas5P(), "Rodas5P", absTol, relTol, nMultiStarts, algList=optmizersTest[iOptimIPNewtonGN], terminateSSMethod=:NewtonNorm, solverSSRelTol=1e-6, solverSSAbsTol=1e-6, reuseS=false)
end


if ARGS[1] == "Fujita_SciSignal2010"
    pathYML = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", "model_Fujita_SciSignal2010", "Fujita_SciSignal2010.yaml")
    petabModel = readPEtabModel(pathYML, verbose=true)
    benchmarkParameterEstimation(petabModel, Rodas5(), "Rodas5", absTol, relTol, nMultiStarts, algList=optmizersTest[iNotOptimIPNewtonGN]) 
    benchmarkParameterEstimation(petabModel, Rodas5(), "Rodas5", absTol, relTol, nMultiStarts, algList=optmizersTest[iOptimIPNewtonGN], reuseS=false) 
end


if ARGS[1] == "Schwen_PONE2014"
    pathYML = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", "model_Schwen_PONE2014", "Schwen_PONE2014.yaml")
    petabModel = readPEtabModel(pathYML, verbose=true)
    benchmarkParameterEstimation(petabModel, Rodas5P(), "Rodas5P", absTol, relTol, nMultiStarts, algList=optmizersTest[iNotOptimIPNewtonGN], reuseS=true)
    benchmarkParameterEstimation(petabModel, Rodas5P(), "Rodas5P", absTol, relTol, nMultiStarts, algList=optmizersTest[iOptimIPNewtonGN], reuseS=false)
end


if ARGS[1] == "Bruno_JExpBot2016"
    pathYML = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", "model_Bruno_JExpBot2016", "Bruno_JExpBot2016.yaml")
    petabModel = readPEtabModel(pathYML, verbose=true)
    benchmarkParameterEstimation(petabModel, Rodas5P(), "Rodas5P", absTol, relTol, nMultiStarts, algList=optmizersTest[iNotOptimIPNewtonGN]) 
    benchmarkParameterEstimation(petabModel, Rodas5P(), "Rodas5P", absTol, relTol, nMultiStarts, algList=optmizersTest[iOptimIPNewtonGN], reuseS=false) 
end


if ARGS[1] == "Elowitz_Nature2000"
    pathYML = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", "model_Elowitz_Nature2000", "Elowitz_Nature2000.yaml")
    petabModel = readPEtabModel(pathYML, verbose=true)
    benchmarkParameterEstimation(petabModel, Rodas5P(), "Rodas5P", absTol, relTol, nMultiStarts, algList=optmizersTest[iNotOptimIPNewtonGN]) 
    benchmarkParameterEstimation(petabModel, Rodas5P(), "Rodas5P", absTol, relTol, nMultiStarts, algList=optmizersTest[iOptimIPNewtonGN], reuseS=false) 
end


if ARGS[1] == "Isensee_JCB2018"
    pathYML = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", "model_Isensee_JCB2018", "Isensee_JCB2018.yaml")
    petabModel = readPEtabModel(pathYML, verbose=true)
    removeAllProcs()
    addprocs(1, exeflags="--project=.")
    benchmarkParameterEstimation(petabModel, QNDF(), "QNDF", absTol, relTol, nMultiStarts, algList=optmizersTest[iNotOptimIPNewtonGN], numberOfprocesses=2, terminateSSMethod=:NewtonNorm, solverSSRelTol=1e-6, solverSSAbsTol=1e-6, reuseS=true) 
    benchmarkParameterEstimation(petabgModel, QNDF(), "QNDF", absTol, relTol, nMultiStarts, algList=optmizersTest[iOptimIPNewtonGN], numberOfprocesses=2, terminateSSMethod=:NewtonNorm, solverSSRelTol=1e-6, solverSSAbsTol=1e-6, reuseS=false) 
end


if ARGS[1] == "Isensee_JCB2018"
    pathYML = joinpath(@__DIR__, "..", "..", "Intermediate", "PeTab_models", "model_Isensee_JCB2018", "Isensee_JCB2018.yaml")
    petabModel = readPEtabModel(pathYML, verbose=true)
    removeAllProcs()
    addprocs(1, exeflags="--project=.")
    benchmarkParameterEstimation(petabModel, QNDF(), "QNDF", absTol, relTol, nMultiStarts, algList=optmizersTest[iNotOptimIPNewtonGN], numberOfprocesses=2, terminateSSMethod=:NewtonNorm, solverSSRelTol=1e-6, solverSSAbsTol=1e-6, reuseS=true) 
    benchmarkParameterEstimation(petabgModel, QNDF(), "QNDF", absTol, relTol, nMultiStarts, algList=optmizersTest[iOptimIPNewtonGN], numberOfprocesses=2, terminateSSMethod=:NewtonNorm, solverSSRelTol=1e-6, solverSSAbsTol=1e-6, reuseS=false) 
end


if ARGS[1] == "Borghans_BiophysChem1997"
    pathYML = pwd() * "/Intermediate/PeTab_models/model_Borghans_BiophysChem1997/Borghans_BiophysChem1997.yaml"
    petabModel = readPEtabModel(pathYML, verbose=true)
    removeAllProcs()
    addprocs(1, exeflags="--project=.")
    benchmarkParameterEstimation(petabModel, Rodas5P(), "Rodas5P", absTol, relTol, nMultiStarts, algList=optmizersTest[iNotOptimIPNewtonGN], reuseS=true) 
    benchmarkParameterEstimation(petabgModel, Rodas5P(), "Rodas5P", absTol, relTol, nMultiStarts, algList=optmizersTest[iOptimIPNewtonGN], reuseS=false) 
end


pathYML = pwd() * "/Intermediate/PeTab_models/model_Borghans_BiophysChem1997/Borghans_BiophysChem1997.yaml"
petabModel = readPEtabModel(pathYML, verbose=true, forceBuildJuliaFiles=true)
petabProblem = setUpPEtabODEProblem(petabModel, Rodas5P(), solverAbsTol=1e-8, solverRelTol=1e-8)
                                        
θ_estNames = string.(petabProblem.θ_estNames)

pathCube = joinpath(petabModel.dirJulia, "Cube_benchmark.csv")
createCube(pathCube, petabProblem, 1000, seed=123, verbose=true)
cube = Matrix(CSV.read(pathCube, DataFrame))

cost = petabProblem.computeCost(petabProblem.θ_nominalT)
gradient = zeros(length(petabProblem.θ_nominalT))
hessian = zeros(length(petabProblem.θ_nominalT), length(petabProblem.θ_nominalT))
@elapsed petabProblem.computeGradientAutoDiff(gradient, petabProblem.θ_nominalT)
@elapsed petabProblem.computeHessian(hessian, petabProblem.θ_nominalT)

