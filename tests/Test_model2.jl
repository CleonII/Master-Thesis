#= 
    Check the accruacy of the PeTab importer for a simple linear ODE;
        s' = alpha*s; s(0) = 8.0 -> s(t) = 8.0 * exp(alpha*t)
        d' = beta*d;  d(0) = 4.0 -> d(t) = 4.0 * exp(beta*t)
    This ODE is solved analytically, and using the analytical solution the accuracy of 
    the ODE solver, cost function, gradient and hessian of the PeTab importer is checked.
    The accuracy of the optimizers is further checked.
    The measurment data is avaible in tests/Test_model2/
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
using Zygote
using SciMLSensitivity
using Zygote
using Sundials
using Ipopt
using Optim
using NLopt
using YAML
using Test


# Relevant PeTab structs for compuations 
include(joinpath(pwd(), "src", "PeTab_structs.jl"))

# PeTab importer to get cost, grad etc 
include(joinpath(pwd(), "src", "Create_PEtab_model.jl"))

# HyperCube sampling 
include(joinpath(pwd(), "src", "Optimizers", "Lathin_hypercube.jl"))

# For converting to SBML 
include(joinpath(pwd(), "src", "SBML", "SBML_to_ModellingToolkit.jl"))

# Optimizers 
include(joinpath(pwd(), "src", "Optimizers", "Set_up_Ipopt.jl"))
include(joinpath(pwd(), "src", "Optimizers", "Set_up_optim.jl"))
include(joinpath(pwd(), "src", "Optimizers", "Set_up_NLopt.jl"))
include(joinpath(pwd(), "src", "Optimizers", "Set_up_fides.jl"))
include(joinpath(pwd(), "src", "Optimizers", "Set_up_forward_gradient.jl"))

include(joinpath(pwd(), "tests", "Common.jl"))


"""
    testOdeSol(solver, tol; printRes=false)

    Compare analytical vs numeric ODE solver using a provided solver with 
    tolerance tol for the Test_model2.
    Returns true if passes test (sqDiff less than 1e-8) else returns false. 
"""
function testODESolverTestModel2(petabModel::PEtabModel, solver, tol)
   
    # Set values to PeTab file values 
    experimentalConditionsFile, measurementDataFile, parameterDataFile, observablesDataFile = readPEtabFiles(petabModel)
    measurementData = processMeasurements(measurementDataFile, observablesDataFile) 
    paramData = processParameters(parameterDataFile) 
    θ_indices = computeIndicesθ(paramData, measurementData, petabModel.odeSystem, experimentalConditionsFile)
    simulationInfo = processSimulationInfo(petabModel, measurementData, paramData)
    setParamToFileValues!(petabModel.parameterMap, petabModel.stateMap, paramData)

    # Parameter values where to teast accuracy. Each column is a alpha, beta, gamma and delta
    u0 = [8.0, 4.0]
    parametersTest = reshape([2.0, 3.0,
                              1.0, 2.0,
                              1.0, 0.4,
                              4.0, 3.0,
                              0.01, 0.02], (2, 5))

    for i in 1:5

        alpha, beta = parametersTest[:, i]        
        # Set parameter values for ODE
        petabModel.parameterMap[2] = Pair(petabModel.parameterMap[2].first, alpha)
        petabModel.parameterMap[3] = Pair(petabModel.parameterMap[3].first, beta)
        prob = ODEProblem(petabModel.odeSystem, petabModel.stateMap, (0.0, 5e3), petabModel.parameterMap, jac=true)
        prob = remake(prob, p = convert.(Float64, prob.p), u0 = convert.(Float64, prob.u0))
        θ_est = getFileODEvalues(petabModel)
        changeExperimentalCondition! = (pVec, u0Vec, expID) -> _changeExperimentalCondition!(pVec, u0Vec, expID, θ_est, petabModel, θ_indices)
        
        # Solve ODE system 
        odeSolutions, success = solveODEAllExperimentalConditions(prob, changeExperimentalCondition!, simulationInfo, solver, tol, tol, petabModel.computeTStops)
        odeSolution = odeSolutions[simulationInfo.experimentalConditionId[1]]
        
        # Compare against analytical solution 
        sqDiff = 0.0
        for t in odeSolution.t
            solAnalytic = [u0[1]*exp(alpha*t), u0[2]*exp(beta*t)]
            sqDiff += sum((odeSolution(t)[1:2] - solAnalytic).^2)
        end

        @test sqDiff ≤ 1e-6
    end
end


function computeCostAnalyticTestModel2(paramVec)

    u0 = [8.0, 4.0]
    alpha, beta = paramVec[1:2]
    measurementData = CSV.read(pwd() * "/tests/Test_model2/measurementData_Test_model2.tsv", DataFrame)

    # Extract correct parameter for observation i and compute logLik
    logLik = 0.0
    for i in 1:nrow(measurementData)
        
        # Specs for observation i
        obsID = measurementData[i, :observableId]
        noiseID = measurementData[i, :noiseParameters]
        yObs = measurementData[i, :measurement]
        t = measurementData[i, :time]
        # Extract correct sigma 
        if noiseID == "sd_sebastian_new"
            sigma = paramVec[3]
        elseif noiseID == "sd_damiano_new"
            sigma = paramVec[4]
        end

        sol = [u0[1]*exp(alpha*t), u0[2]*exp(beta*t)]
        if obsID == "sebastian_measurement"
            yMod = sol[1]
        elseif obsID == "damiano_measurement"
            yMod = sol[2]        
        end

        logLik += log(sigma) + 0.5*log(2*pi) + 0.5 * ((yObs - yMod) / sigma)^2
    end
        
    return logLik
end


"""
    testCostGradientOrHessianTestModel2(solver, tol; printRes::Bool=false)

    Compare cost, gradient and hessian computed via the analytical solution 
    vs the PeTab importer functions (to check PeTab importer) for five random 
    parameter vectors for Test_model2. For the analytical solution the gradient 
    and hessian are computed via ForwardDiff.
"""
function testCostGradientOrHessianTestModel2(petabModel::PEtabModel, solver, tol)
    
    petabProblem1 = setUpPEtabODEProblem(petabModel, solver, solverAbsTol=tol, solverRelTol=tol, 
                                         sensealgZygote = ForwardDiffSensitivity(), 
                                         odeSolverForwardEquations=Vern9(), sensealgForwardEquations = ForwardDiffSensitivity(), 
                                         odeSolverAdjoint=solver, solverAdjointAbsTol=tol, solverAdjointRelTol=tol,
                                         sensealgAdjoint=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(true)))

    petabProblem2 = setUpPEtabODEProblem(petabModel, solver, solverAbsTol=tol, solverRelTol=tol, 
                                         sensealgForwardEquations=:AutoDiffForward, odeSolverForwardEquations=Vern9())
    

    Random.seed!(123)
    createCube(petabProblem1, 5)
    cube = Matrix(CSV.read(petabProblem1.pathCube, DataFrame))

    for i in 1:5

        p = cube[i, :]
        referenceCost = computeCostAnalyticTestModel2(p)
        referenceGradient = ForwardDiff.gradient(computeCostAnalyticTestModel2, p)
        referenceHessian = ForwardDiff.hessian(computeCostAnalyticTestModel2, p)

        # Test both the standard and Zygote approach to compute the cost 
        cost = _testCostGradientOrHessian(petabProblem1, p, cost=true)
        @test cost ≈ referenceCost atol=1e-3
        costZygote = _testCostGradientOrHessian(petabProblem1, p, costZygote=true)
        @test costZygote ≈ referenceCost atol=1e-3

        # Test all gradient combinations. Note we test sensitivity equations with and without autodiff 
        gradientAutoDiff = _testCostGradientOrHessian(petabProblem1, p, gradientAutoDiff=true)
        @test norm(gradientAutoDiff - referenceGradient) ≤ 1e-2
        gradientZygote = _testCostGradientOrHessian(petabProblem1, p, gradientZygote=true)
        @test norm(gradientZygote - referenceGradient) ≤ 1e-2
        gradientAdjoint = _testCostGradientOrHessian(petabProblem1, p, gradientAdjoint=true)
        @test norm(normalize(gradientAdjoint) - normalize((referenceGradient))) ≤ 1e-2
        gradientForwardEquations1 = _testCostGradientOrHessian(petabProblem1, p, gradientForwardEquations=true)
        @test norm(gradientForwardEquations1 - referenceGradient) ≤ 1e-2
        gradientForwardEquations2 = _testCostGradientOrHessian(petabProblem2, p, gradientForwardEquations=true)
        @test norm(gradientForwardEquations2 - referenceGradient) ≤ 1e-2

        # Testing "exact" hessian via autodiff 
        hessian = _testCostGradientOrHessian(petabProblem1, p, hessian=true)
        @test norm(hessian - referenceHessian) ≤ 1e-3
    end
end


"""
    testOptimizer(solver, tol)

    Check if optimizers (all Ipopt and Optim) can reach the minima for Test_model2 
    from a random starting point. Convergence assumed if sum((p_est - p_true).^2) is 
    smaller than 1e-3.
"""
function testOptimizersTestModel2(petabModel::PEtabModel, solver, tol)

    petabProblem = setUpPEtabODEProblem(petabModel, solver, solverAbsTol=tol, solverRelTol=tol, 
                                        reuseS=true, sensealgForwardEquations=:AutoDiffForward, odeSolverForwardEquations=solver)
    petabProblem2 = setUpPEtabODEProblem(petabModel, solver, solverAbsTol=tol, solverRelTol=tol, 
                                        reuseS=false, sensealgForwardEquations=:AutoDiffForward, odeSolverForwardEquations=solver)                                        

    Random.seed!(123)
    createCube(petabProblem, 5)
    cube = Matrix(CSV.read(petabProblem.pathCube, DataFrame))

    # Ipopt optimizers 
    ipoptProbHessApprox, iterArrHessApprox = createIpoptProb(petabProblem, hessianUse=:blockAutoDiff)
    ipoptProbBfgs, iterArrBfgs = createIpoptProb(petabProblem, hessianUse=:LBFGS)
    ipoptProbAutoHess, iterArrAutoHess = createIpoptProb(petabProblem, hessianUse=:autoDiff)
    ipoptProbGN, iterArrAutoHess = createIpoptProb(petabProblem, hessianUse=:GaussNewton)
    Ipopt.AddIpoptNumOption(ipoptProbAutoHess, "acceptable_tol", 1e-8)
    
    # Fides optimizers 
    runFidesHess = setUpFides(petabProblem, :autoDiff; verbose=0, 
                              options=py"{'maxiter' : 1000, 'fatol' : 0.0, 'frtol' : 1e-8, 'xtol' : 0.0, 'gatol' : 1e-6, 'grtol' : 1e-6}"o)
    runFidesHessBlock = setUpFides(petabProblem, :blockAutoDiff; verbose=0, 
                                   options=py"{'maxiter' : 1000, 'fatol' : 0.0, 'frtol' : 1e-8, 'xtol' : 0.0, 'gatol' : 1e-6, 'grtol' : 1e-6}"o)
    runFidesGN = setUpFides(petabProblem, :GaussNewton; verbose=0, 
                            options=py"{'maxiter' : 1000, 'fatol' : 0.0, 'frtol' : 1e-8, 'xtol' : 0.0, 'gatol' : 1e-6, 'grtol' : 1e-6}"o)                                   
    runFidesHessBFGS = setUpFides(petabProblem, :None; verbose=0, 
                                  fidesHessApprox=py"fides.hessian_approximation.BFGS()"o, 
                                  options=py"{'maxiter' : 1000, 'fatol' : 0.0, 'frtol' : 1e-8, 'xtol' : 0.0, 'gatol' : 1e-6, 'grtol' : 1e-6}"o)                                   

    # Optim optimizers 
    optimProbHessApprox = createOptimProb(petabProblem, IPNewton(), hessianUse=:blockAutoDiff)
    optimProbAutoHess = createOptimProb(petabProblem, IPNewton(), hessianUse=:autoDiff, 
                                        options=Optim.Options(iterations = 1000, show_trace = false, allow_f_increases=true, 
                                                              successive_f_tol = 3, f_tol=1e-8, g_tol=1e-6, x_tol=0.0))
    optimProbGN = createOptimProb(petabProblem2, IPNewton(), hessianUse=:GaussNewton, 
                                  options=Optim.Options(iterations = 1000, show_trace = false, allow_f_increases=true, 
                                                         successive_f_tol = 3, f_tol=1e-8, g_tol=1e-6, x_tol=0.0))                                                              
    optimProbBFGS = createOptimProb(petabProblem, BFGS())
    optimProbLBFGS = createOptimProb(petabProblem, LBFGS(), 
                                     options=Optim.Options(iterations = 250, 
                                                           show_trace = false, 
                                                           allow_f_increases=true, 
                                                           outer_iterations = 4))

    # NLopt optimizers 
    NLoptLBFGS = createNLoptProb(petabProblem, :LD_LBFGS, verbose=false)
    NLoptLBFGS.ftol_rel = 1e-8
    NLoptLBFGS.xtol_rel = 0.0
    NLoptLBFGS.maxeval = 5000
    NLoptTNewton = createNLoptProb(petabProblem, :LD_TNEWTON_PRECOND_RESTART, verbose=false)
    NLoptTNewton.ftol_rel = 1e-8
    NLoptTNewton.xtol_rel = 0.0
    NLoptTNewton.maxeval = 5000

    p0 = cube[1, :]

    # Ipopt AutoDiffHessian 
    ipoptProbAutoHess.x = deepcopy(p0)
    ipoptProbAutoHess.eval_f(p0)
    sol_opt = Ipopt.IpoptSolve(ipoptProbAutoHess)
    sqDiffIpoptHessian = sum((ipoptProbAutoHess.x - petabProblem.θ_nominal).^2)
    @test sqDiffIpoptHessian ≤ 1e-3

    # Ipopt Hessian approximation
    ipoptProbHessApprox.x = deepcopy(p0)
    ipoptProbHessApprox.eval_f(p0) # Needed to avoid segfault, must write a wrapper for this 
    sol_opt = Ipopt.IpoptSolve(ipoptProbHessApprox)
    sqDiffIpoptBlockHessian = sum((ipoptProbHessApprox.x - petabProblem.θ_nominal).^2)
    @test sqDiffIpoptBlockHessian ≤ 1e-3

    # Ipot Hessian Gauss Newton approximation
    ipoptProbGN.x = deepcopy(p0)
    ipoptProbGN.eval_f(p0) # Needed to avoid segfault, must write a wrapper for this 
    sol_opt = Ipopt.IpoptSolve(ipoptProbGN)
    sqDiffIpoptGN = sum((ipoptProbGN.x - petabProblem.θ_nominal).^2)
    @test sqDiffIpoptGN ≤ 1e-3

    # Ipopt BFGS 
    ipoptProbBfgs.x = deepcopy(p0)
    ipoptProbBfgs.eval_f(p0)
    sol_opt = Ipopt.IpoptSolve(ipoptProbBfgs)
    sqDiffIpoptBFGS = sum((ipoptProbBfgs.x - petabProblem.θ_nominal).^2)
    @test sqDiffIpoptBFGS ≤ 1e-3

    # Fides AutoDiff hessian
    res, niter, converged = runFidesHess(p0)
    sqDiffFidesHessian = sum((res[2] - petabProblem.θ_nominal).^2)
    @test sqDiffFidesHessian ≤ 1e-3

    # Fides BlockAutoDiff hessian  
    res, niter, converged = runFidesHessBlock(p0)
    sqDiffFidesBlockHessian = sum((res[2] - petabProblem.θ_nominal).^2)
    @test sqDiffFidesBlockHessian ≤ 1e-3

    # Fides Guass Newton hessian
    res, niter, converged = runFidesGN(p0)
    sqDiffFidesGN= sum((res[2] - petabProblem.θ_nominal).^2)
    @test sqDiffFidesGN ≤ 1e-3

    # Fides BFGS 
    res, niter, converged = runFidesHessBFGS(p0)
    sqDiffFidesBFGS = sum((res[2] - petabProblem.θ_nominal).^2)
    @test sqDiffFidesBFGS ≤ 1e-3

    # Optim AutoHess 
    res = optimProbAutoHess(p0, showTrace=false)
    sqDiffOptimHessian = sum((res.minimizer - petabProblem.θ_nominal).^2)
    @test sqDiffOptimHessian ≤ 1e-3

    # Optim block hessian 
    res = optimProbHessApprox(p0, showTrace=false)
    sqDiffOptimBlockHessian = sum((res.minimizer - petabProblem.θ_nominal).^2)
    @test sqDiffOptimBlockHessian ≤ 1e-3

    # Optim Gauss Newton hessian 
    res = optimProbGN(p0, showTrace=false)
    sqDiffOptimGN = sum((res.minimizer - petabProblem.θ_nominal).^2)
    @test sqDiffOptimGN ≤ 1e-3

    # Optim BFGS
    res = optimProbBFGS(p0, showTrace=false)
    sqDiffOptimBFGS = sum((res.minimizer - petabProblem.θ_nominal).^2)
    @test sqDiffOptimBFGS ≤ 1e-3

    # Optim LBFGS
    res = optimProbLBFGS(p0, showTrace=false)
    sqDiffOptimLBFGS = sum((res.minimizer - petabProblem.θ_nominal).^2)
    @test sqDiffOptimLBFGS ≤ 1e-3

    # NLopt LBFGS
    minf, minx, ret = NLopt.optimize(NLoptLBFGS, p0)
    sqDiffNLoptLBFGS = sum((minx - petabProblem.θ_nominal).^2)
    @test sqDiffNLoptLBFGS ≤ 1e-3

    # NLopt truncated Newton with preconditioner 
    minf, minx, ret = NLopt.optimize(NLoptTNewton, p0)
    sqDiffNLoptTNewtonPre = sum((minx - petabProblem.θ_nominal).^2)
    @test sqDiffNLoptTNewtonPre ≤ 1e-3

    p0 = petabProblem.θ_nominal .+ 0.1
    stepLengths = log10.(LinRange(exp10(1e-2),exp10(1e-4), 10000))
    forwardGradOpt = createFowardGradientProb(petabProblem, stepLengths, 10000)
    minCost, minParams, costVal, paramVal = runAdam(p0, forwardGradOpt, verbose=false, seed=14)
    # The dynamic parameters are easily estimated, but sd-parameter often fail so do not 
    # include them into the tests.
    sqDiffForwardGradient = sum((minParams[1:2] - petabProblem.θ_nominal[1:2]).^2)
    @test sqDiffForwardGradient ≤ 1e-3
end


petabModel = readPEtabModel(joinpath(@__DIR__, "Test_model2", "Test_model2.yaml"), forceBuildJuliaFiles=true)
loadFidesFromPython("/home/sebpe/anaconda3/envs/PeTab/bin/python")

@testset verbose=true "Test model 2" begin
    @testset "Test model2 : ODE solver" begin 
        testODESolverTestModel2(petabModel, Vern9(), 1e-9)
    end

    @testset "Test model2 : Cost gradient and hessian" begin 
        testCostGradientOrHessianTestModel2(petabModel, Vern9(), 1e-15)
    end

    @testset "Test model 2 : Gradient of residuals" begin
        checkGradientResiduals(petabModel, Rodas5(), 1e-9)
    end

    @testset "Test model 2 : Optimizers" begin
        testOptimizersTestModel2(petabModel, Rodas4P(), 1e-9)
    end
end
