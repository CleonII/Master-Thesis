using ModelingToolkit, DifferentialEquations, BenchmarkTools, DataFrames, CSV, LSODA, Sundials, ODEInterface, ODEInterfaceDiffEq
using ForwardDiff
using ReverseDiff
using Plots
using StatsBase
using Random
using LinearAlgebra
using Calculus
using Ipopt
using Optim
using Distributions
#plotly()

# TODO : Fix input order of hessian functions 
# TODO : Fix ultra slow calc-log lik. Can be done by precomputing all indices 

include(joinpath(pwd(), "Pipeline_ModelSolver", "BigFloatODEProblem.jl"))
include(joinpath(pwd(), "Additional_functions", "additional_tools.jl"))
include(joinpath(pwd(), "Additional_functions", "Solver_info.jl"))
include(joinpath(pwd(), "Additional_functions", "benchmarkSolvers.jl"))

# Ipopt wrapper 
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "CommonParameterEstimationMethods", "CreateIpoptProb.jl"))

# PeTab importer 
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "PeTabImporter.jl"))

# Observation function, needs to be generalised 
include(pwd() * "/Pipeline_ModelParameterEstimation/Data/model_Boehm_JProteomeRes2014/Boehm_JProteomeRes2014Obs.jl")
include(pwd() * "/Pipeline_ModelParameterEstimation/Data/model_Bachmann_MSB2011/Bachmann_MSB2011Obs.jl")

allModelFunctionVector = includeAllModels(getModelFiles(joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels")), 
        joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels"))
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "LatinHyperCubeParameters", "LatinHyperCubeSampledParameters.jl"))

# Will be easier to write general code when a struct has been implemented for holding a model. 
function runBoehm(solver, tol, fileNameSave; useAutoH::Bool=false)

        modelName = "model_Boehm_JProteomeRes2014"
        evalObs = Boehm_JProteomeRes2014
        evalU0 = Boehm_JProteomeRes2014_t0!
        evalSd = Boehm_JProteomeRes2014_sd!

        evalF, evalGradF, evalHessianApproxF, paramVecEstTmp, lowerBounds, upperBounds, idParam = setUpCostFunc(modelName, evalObs, evalU0, evalSd, solver, tol)
        # "Exact" hessian via autodiff (should only be used for smaller models)
        evalH = (hessianMat, paramVec) -> begin hessianMat .= Symmetric(ForwardDiff.hessian(evalF, paramVec)) end

        # Create a Lathin-hypercube (generating start guesses)
        dirSave = pwd() * "/Intermediate/model_Boehm_JProteomeRes2014/"
        fileSaveCube = dirSave * "CubeOpt.csv"
        println("Creating cube")
        createCube(500, lowerBounds, upperBounds, fileSaveCube, evalF)
        println("Done creating cube")
        cube = Matrix(CSV.read(fileSaveCube, DataFrame))

        fileSave = dirSave * fileNameSave
        runOptimizer(fileSave, evalF, evalGradF, evalHessianApproxF, evalH, lowerBounds, upperBounds, cube, useAutoH)
end



function writeFile(fileSave, cost, runTime, retCode, nIter, startGuess, alg::String)
    # Save after each iteration (do not loose data)
    dataSave = [alg cost runTime retCode nIter startGuess]
    dataSave = DataFrame(dataSave, ["Alg", "Cost", "Run_time", "Ret_code", "N_iter", "Start_guess"])
    shouldAppend = isfile(fileSave) ? true : false
    CSV.write(fileSave, dataSave, append=shouldAppend)

end


function runOptimizer(fileSave, evalF, evalGradF, evalHessianApproxF, evalH, lowerBounds, upperBounds, cube::Array{Float64, 2}, useExactH::Bool)

    println("Running benchmark, saving results in $fileSave")

    # Three different ipopt problems 
    nParam = length(lowerBounds)
    ipoptProbHessApprox, iterArrHessApprox = createIpoptProbNew(evalF, evalGradF, evalHessianApproxF, lowerBounds, upperBounds, emptyH=false) 
    ipoptProbBfgs, iterArrBfgs = createIpoptProbNew(evalF, evalGradF, evalHessianApproxF, lowerBounds, upperBounds, emptyH=true) 
    if useExactH == true
        ipoptProbAutoHess, iterArrAutoHess = createIpoptProbNew(evalF, evalGradF, evalH, lowerBounds, upperBounds, emptyH=false) 
    else
        ipoptProbAutoHess, iterArrAutoHess = nothing, nothing
    end

    # Two different Optim problems 
    optimProbHessApprox = createOptimProb(evalF, evalGradF, evalHessianApproxF, lowerBounds, upperBounds, showTrace=false)
    if useExactH == true
        optimProbAutoHess = createOptimProb(evalF, evalGradF, evalH, lowerBounds, upperBounds, showTrace=false)
    else    
        optimProbAutoHess = nothing
    end

    # Evaluate gradient + hessian (allow precompilation to get fair timings for benchmark)
    println("Evaluting gradient and hessians first time")
    ipoptProbHessApprox.eval_f(cube[1, :])
    ipoptProbHessApprox.eval_grad_f(cube[1, :], zeros(nParam))
    evalHessianApproxF(zeros(nParam, nParam), cube[1, :])
    if useExactH == true
        evalH(zeros(nParam, nParam), cube[1, :])
    end

    
    nEvals = size(cube)[1]
    println("Starting benchmark")
    for i in 1:nEvals

        println("I = $i of $nEvals")
        p0 = cube[i, :] # Sample from hypercube 
        
        # Ipopt with hessian approximation
        ipoptProbHessApprox.x = deepcopy(p0)
        benchRunTime = @elapsed  sol_opt = Ipopt.IpoptSolve(ipoptProbHessApprox)
        writeFile(fileSave, ipoptProbHessApprox.obj_val, benchRunTime, ipoptProbHessApprox.status, iterArrHessApprox[1], i, "ipoptHessApprox")

        # Ipopt BFGS 
        ipoptProbBfgs.x = deepcopy(p0)
        benchRunTime = @elapsed  sol_opt = Ipopt.IpoptSolve(ipoptProbBfgs)
        writeFile(fileSave, ipoptProbBfgs.obj_val, benchRunTime, ipoptProbBfgs.status, iterArrBfgs[1], i, "ipoptBfgs")

        # In case of exact hessian approximation 
        if useExactH == true
            ipoptProbAutoHess.x = deepcopy(p0)
            benchRunTime = @elapsed  sol_opt = Ipopt.IpoptSolve(ipoptProbAutoHess)
            writeFile(fileSave, ipoptProbAutoHess.obj_val, benchRunTime, ipoptProbAutoHess.status, iterArrAutoHess[1], i, "ipoptHessAuto")
        end

        # Optim with approximate hessian 
        res = optimProbHessApprox(p0)
        writeFile(fileSave, res.minimum, res.time_run, res.f_converged, res.iterations, i, "optimHessApprox")

        # Optim with Hessian 
        if useExactH == true
            res = optimProbAutoHess(p0)
            writeFile(fileSave, res.minimum, res.time_run, res.f_converged, res.iterations, i, "optimHessAuto")
        end        
    end
end



solver = QNDF()
tol = 1e-9
fileName = "Bohem_opt_qndf.csv"
runBoehm(solver, tol, fileName, useAutoH=true)

#=


# For Bachman model 
modelName = "model_Bachmann_MSB2011"
evalObs = Bachmann_MSB2011
evalU0 = Bachmann_MSB2011_t0!
evalSd = Bachmann_MSB2011_sd!

solver = QNDF()
tol = 1e-8

evalF, evalGradF, evalHessianApproxF, paramVecEstTmp, lowerBounds, upperBounds, idParam = setUpCostFunc(modelName, evalObs, evalU0, evalSd, solver, tol)
ipoptProbHessApprox, iterArrHessApprox = createIpoptProbNew(evalF, evalGradF, evalHessianApproxF, lowerBounds, upperBounds, emptyH=true) 
p0 = [rand(Uniform(lowerBounds[i], upperBounds[i])) for i in eachindex(lowerBounds)]
costStart = ipoptProbHessApprox.eval_f(p0)
sol_opt = Ipopt.IpoptSolve(ipoptProbHessApprox)
a = 1




#=x0
println("Starting compute gradients")
gradUse = zeros(length(paramVecEstTmp))
b = @elapsed evalGradF(paramVecEstTmp, gradUse)
=#

#hessianUse = zeros(length(paramVecEstTmp), length(paramVecEstTmp))
#b = @elapsed evalHessianApproxF(hessianUse, paramVecEstTmp)
#4534.0076376876295
=#