# Check the accruacy of the PeTab importer for a simple linear ODE;
#   s' = alpha*s; s(0) = 8.0 -> s(t) = 8.0 * exp(alpha*t)
#   d' = beta*d;  d(0) = 4.0 -> d(t) = 4.0 * exp(beta*t)
# This ODE is solved analytically, and using the analytical solution the accuracy of 
# the ODE solver, cost function, gradient and hessian of the PeTab importer is checked.
# The accuracy of the optimizers is further checked.
# The measurment data is avaible in tests/Test_model2/


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
using LineSearches
using SciMLSensitivity
using Zygote
using Sundials


# Relevant PeTab structs for compuations 
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
include(joinpath(pwd(), "src", "Optimizers", "Set_up_forward_gradient.jl"))


"""
    testOdeSol(solver, tol; printRes=false)

    Compare analytical vs numeric ODE solver using a provided solver with 
    tolerance tol for the Test_model2.
    Returns true if passes test (sqDiff less than 1e-8) else returns false. 
"""
function testOdeSol(peTabModel::PeTabModel, solver, tol; printRes=false)
   
    # Set values to PeTab file values 
    experimentalConditionsFile, measurementDataFile, parameterDataFile, observablesDataFile = readDataFiles(peTabModel.dirModel, readObs=true)
    measurementData = processMeasurementData(measurementDataFile, observablesDataFile) 
    paramData = processParameterData(parameterDataFile) 
    setParamToFileValues!(peTabModel.paramMap, peTabModel.stateMap, paramData)
    
    # Extract experimental conditions for simulations 
    simulationInfo = getSimulationInfo(peTabModel, measurementDataFile, measurementData)

    # Parameter values where to teast accuracy. Each column is a alpha, beta, gamma and delta
    u0 = [8.0, 4.0]
    paramMat = reshape([2.0, 3.0,
                        1.0, 2.0,
                        1.0, 0.4,
                        4.0, 3.0,
                        0.01, 0.02], (2, 5))

    for i in 1:5

        alpha, beta = paramMat[:, i]        
        # Set parameter values for ODE
        peTabModel.paramMap[2] = Pair(peTabModel.paramMap[2].first, alpha)
        peTabModel.paramMap[3] = Pair(peTabModel.paramMap[3].first, beta)
        prob = ODEProblem(peTabModel.odeSystem, peTabModel.stateMap, (0.0, 5e3), peTabModel.paramMap, jac=true)
        prob = remake(prob, p = convert.(Float64, prob.p), u0 = convert.(Float64, prob.u0))
        changeToExperimentalCondUse! = (pVec, u0Vec, expID) -> changeExperimentalCond!(pVec, u0Vec, expID, paramData, experimentalConditionsFile, peTabModel)
        
        # Solve ODE system 
        solArray, success = solveOdeModelAllExperimentalCond(prob, changeToExperimentalCondUse!, simulationInfo, solver, tol, tol, peTabModel.getTStops)
        solNumeric = solArray[1]
        
        # Compare against analytical solution 
        sqDiff = 0.0
        for t in solNumeric.t
            solAnalytic = [u0[1]*exp(alpha*t), u0[2]*exp(beta*t)]
            sqDiff += sum((solNumeric(t)[1:2] - solAnalytic).^2)
        end

        if sqDiff > 1e-8
            @printf("sqDiff = %.3e\n", sqDiff)
            @printf("Does not pass test\n")
            return false
        end

        if printRes == true
            @printf("sqDiff = %.3e\n", sqDiff)
        end
    end

    return true
end


function calcCostAnalytic(paramVec)

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
    testCostGradHess(solver, tol; printRes::Bool=false)

    Compare cost, gradient and hessian computed via the analytical solution 
    vs the PeTab importer functions (to check PeTab importer) for five random 
    parameter vectors for Test_model2. For the analytical solution the gradient 
    and hessian are computed via ForwardDiff.
"""
function testCostGradHess(peTabModel::PeTabModel, solver, tol; printRes::Bool=false)
    
    peTabOpt = setUpCostGradHess(peTabModel, solver, tol, sensealg = ForwardDiffSensitivity(), 
                                 sensealgForward = ForwardSensitivity(), solverForward=solver,
                                 adjSolver=solver, adjTol=tol, adjSensealg=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(true)))
    peTabOptAlt = setUpCostGradHess(peTabModel, solver, tol, 
                                    sensealgForward = :AutoDiffForward, solverForward=Vern9())                                 

    Random.seed!(123)
    createCube(peTabOpt, 5)
    cube = Matrix(CSV.read(peTabOpt.pathCube, DataFrame))
    nParam = size(cube)[2]

    for i in 1:5

        paramVec = cube[i, :]

        # Evaluate cost 
        costPeTab = peTabOpt.evalF(paramVec)
        costAnalytic = calcCostAnalytic(paramVec)
        sqDiffCost = (costPeTab - costAnalytic)^2
        if sqDiffCost > 1e-6
            @printf("sqDiffCost = %.3e\n", sqDiffCost)
            @printf("Does not pass test on cost\n")
            return false
        end

        # Evaluate gradient 
        gradAnalytic = ForwardDiff.gradient(calcCostAnalytic, paramVec)
        gradNumeric = zeros(nParam); peTabOpt.evalGradF(gradNumeric, paramVec)
        sqDiffGrad = sum((gradAnalytic - gradNumeric).^2)
        if sqDiffGrad > 1e-4
            @printf("sqDiffGrad = %.3e\n", sqDiffGrad)
            @printf("Does not pass test on gradient\n")
            return false
        end

        # Evalute Zygote cost function 
        costPeTab = peTabOpt.evalFZygote(paramVec)
        sqDiffCostZygote = (costPeTab - costAnalytic)^2
        if sqDiffCostZygote > 1e-6
            @printf("sqDiffCost = %.3e\n", sqDiffCost)
            @printf("Does not pass test on cost for Zygote cost function\n")
            return false
        end

        # Evalute gradient obtained via Zygote and sensealg 
        gradZygoteSensealg = Zygote.gradient(peTabOpt.evalFZygote, paramVec)[1]
        sqDiffGradZygote = sum((gradAnalytic - gradZygoteSensealg).^2)
        if sqDiffGradZygote > 1e-6
            @printf("sqDiffGrad = %.3e\n", sqDiffGrad)
            @printf("Does not pass test on gradient with Zygote\n")
            return false
        end

        # Forward sensitivity equations gradient
        gradForwardEq = zeros(nParam)
        peTabOpt.evalGradFForwardEq(gradForwardEq, paramVec)
        sqDiffGradForwardEq = sum((gradForwardEq - gradAnalytic).^2)
        if sqDiffGradForwardEq > 1e-5
            @printf("sqDiffGradForwardEq = %.3e\n", sqDiffGradForwardEq)
            @printf("Does not pass test on gradient from Forward sensitivity equations\n")
            return false
        end

        # Forward sensitivity equations using autodiff 
        gradForwardEqAuto = zeros(nParam)
        peTabOptAlt.evalGradFForwardEq(gradForwardEqAuto, paramVec)
        sqDiffGradForwardEqAuto = sum((gradForwardEq - gradAnalytic).^2)
        if sqDiffGradForwardEqAuto > 1e-5
            @printf("sqDiffGradForwardEqAuto = %.3e\n", sqDiffGradForwardEqAuto)
            @printf("Does not pass test on gradient from Forward sensitivity equations\n")
            return false
        end

        # Evaluate lower level adjoint sensitivity interfance gradient 
        gradAdj = zeros(nParam)
        peTabOpt.evalGradFAdjoint(gradAdj, paramVec)
        sqDiffGradAdjoint1 = sum((gradAdj ./ norm(gradAdj) - gradAnalytic ./ norm(gradAnalytic)).^2)
        if sqDiffGradAdjoint1 > 1e-4
            @printf("sqDiffGradAdjointOpt1 = %.3e\n", sqDiffGradAdjoint1)
            @printf("Does not pass test on adjoint gradient gradient\n")
            return false
        end

        # Evaluate hessian 
        hessAnalytic = ForwardDiff.hessian(calcCostAnalytic, paramVec)
        hessNumeric = zeros(nParam, nParam); peTabOpt.evalHess(hessNumeric, paramVec)
        sqDiffHess = sum((hessAnalytic - hessNumeric).^2)
        if sqDiffHess > 1e-3
            @printf("sqDiffHess = %.3e\n", sqDiffHess)
            @printf("Does not pass test on hessian\n")
            return false
        end
                
        if printRes == true
            @printf("sqDiffCost = %.3e\n", sqDiffCost)
            @printf("sqDiffGrad = %.3e\n", sqDiffGrad)
            @printf("sqDiffHess = %.3e\n", sqDiffHess)
            @printf("sqDiffCostZygote = %.3e\n", sqDiffCostZygote)
            @printf("sqDiffGradZygote = %.3e\n", sqDiffGradZygote)
            @printf("sqDiffGradForwardEq = %.3e\n", sqDiffGradForwardEq)
            @printf("sqDiffGradForwardEqAuto = %.3e\n", sqDiffGradForwardEqAuto)
            @printf("sqDiffGradAdjointOpt1 = %.3e\n", sqDiffGradAdjoint1)
        end
    end

    return true
end


"""
    testOptimizer(solver, tol)

    Check if optimizers (all Ipopt and Optim) can reach the minima for Test_model2 
    from a random starting point. Convergence assumed if sum((p_est - p_true).^2) is 
    smaller than 1e-3.
"""
function testOptimizer(peTabModel::PeTabModel, solver, tol)

    peTabOpt = setUpCostGradHess(peTabModel, solver, tol)

    Random.seed!(123)
    createCube(peTabOpt, 5)
    cube = Matrix(CSV.read(peTabOpt.pathCube, DataFrame))

    # Ipopt optimizers 
    ipoptProbHessApprox, iterArrHessApprox = createIpoptProb(peTabOpt, hessianUse=:blockAutoDiff)
    ipoptProbBfgs, iterArrBfgs = createIpoptProb(peTabOpt, hessianUse=:LBFGS)
    ipoptProbAutoHess, iterArrAutoHess = createIpoptProb(peTabOpt, hessianUse=:autoDiff)
    Ipopt.AddIpoptNumOption(ipoptProbAutoHess, "acceptable_tol", 1e-8)
    
    # Fides optimizers 
    runFidesHess = setUpFides(peTabOpt, :autoDiff; verbose=0, 
                              options=py"{'maxiter' : 1000, 'fatol' : 0.0, 'frtol' : 1e-8, 'xtol' : 0.0, 'gatol' : 1e-6, 'grtol' : 1e-6}"o)
    runFidesHessBlock = setUpFides(peTabOpt, :blockAutoDiff; verbose=0, 
                                   options=py"{'maxiter' : 1000, 'fatol' : 0.0, 'frtol' : 1e-8, 'xtol' : 0.0, 'gatol' : 1e-6, 'grtol' : 1e-6}"o)

    # Optim optimizers 
    optimProbHessApprox = createOptimProb(peTabOpt, IPNewton(), hessianUse=:blockAutoDiff)
    optimProbAutoHess = createOptimProb(peTabOpt, IPNewton(), hessianUse=:autoDiff, 
                                        options=Optim.Options(iterations = 1000, show_trace = false, allow_f_increases=true, 
                                                              successive_f_tol = 3, f_tol=1e-8, g_tol=1e-6, x_tol=0.0))
    optimProbBFGS = createOptimProb(peTabOpt, BFGS())
    optimProbLBFGS = createOptimProb(peTabOpt, LBFGS(), 
                                     options=Optim.Options(iterations = 250, 
                                                           show_trace = false, 
                                                           allow_f_increases=true, 
                                                           outer_iterations = 4))

    # NLopt optimizers 
    NLoptLBFGS = createNLoptProb(peTabOpt, :LD_LBFGS, verbose=false)
    NLoptLBFGS.ftol_rel = 1e-8
    NLoptLBFGS.xtol_rel = 0.0
    NLoptLBFGS.maxeval = 5000
    NLoptTNewton = createNLoptProb(peTabOpt, :LD_TNEWTON_PRECOND_RESTART, verbose=false)
    NLoptTNewton.ftol_rel = 1e-8
    NLoptTNewton.xtol_rel = 0.0
    NLoptTNewton.maxeval = 5000

    p0 = cube[1, :]
    # Ipopt Hessian approximation
    ipoptProbHessApprox.x = deepcopy(p0)
    ipoptProbHessApprox.eval_f(p0) # Needed to avoid segfault, must write a wrapper for this 
    sol_opt = Ipopt.IpoptSolve(ipoptProbHessApprox)
    sqDiff = sum((ipoptProbHessApprox.x - peTabOpt.paramVecNotTransformed).^2)
    if sqDiff > 1e-3
        @printf("sqDiffHessApprox = %.3e\n", sqDiff)
        @printf("Failed on optimization for Ipopt hessian approximation\n")
        return false
    else
        @printf("Passed test for Ipopt hess approx\n")
    end
    
    # Ipopt AutoHessian 
    ipoptProbAutoHess.x = deepcopy(p0)
    ipoptProbAutoHess.eval_f(p0)
    sol_opt = Ipopt.IpoptSolve(ipoptProbAutoHess)
    sqDiff = sum((ipoptProbAutoHess.x - peTabOpt.paramVecNotTransformed).^2)
    if sqDiff > 1e-3
        @printf("sqDiffHessApprox = %.3e\n", sqDiff)
        @printf("Failed on optimization for Ipopt auto hessian\n")
        return false
    else
        @printf("Passed test for Ipopt hess autodiff\n")
    end

    # Ipopt BFGS 
    ipoptProbBfgs.x = deepcopy(p0)
    ipoptProbBfgs.eval_f(p0)
    sol_opt = Ipopt.IpoptSolve(ipoptProbBfgs)
    sqDiff = sum((ipoptProbBfgs.x - peTabOpt.paramVecNotTransformed).^2)
    if sqDiff > 1e-3
        @printf("sqDiffHessApprox = %.3e\n", sqDiff)
        @printf("Failed on optimization for Ipopt BFGS\n")
        return false
    else
        @printf("Passed test for Ipopt LBFGS\n")
    end

    # Fides autoDiff hessian 
    res, niter, converged = runFidesHess(p0)
    sqDiff = sum((res[2] - peTabOpt.paramVecNotTransformed).^2)
    if sqDiff > 1e-3
        @printf("sqDiffHessApprox = %.3e\n", sqDiff)
        @printf("Failed on optimization for Fides hessian autodiff\n")
        return false
    else
        @printf("Passed test for Fides hessian autodiff\n")
    end

    # Fides block auto diff hessian 
    res, niter, converged = runFidesHessBlock(p0)
    sqDiff = sum((res[2] - peTabOpt.paramVecNotTransformed).^2)
    if sqDiff > 1e-3
        @printf("sqDiffHessApprox = %.3e\n", sqDiff)
        @printf("Failed on optimization for Fides hessian block autodiff\n")
        return false
    else
        @printf("Passed test for Fides hessian block autodiff\n")
    end

    # Optim HessApprox 
    res = optimProbHessApprox(p0, showTrace=false)
    sqDiff = sum((res.minimizer - peTabOpt.paramVecNotTransformed).^2)
    if sqDiff > 1e-3
        @printf("sqDiffHessApprox = %.3e\n", sqDiff)
        @printf("Failed on optimization for HessApprox\n")
        return false
    else
        @printf("Passed test for Optim interior point newton hessian approximation\n")
    end

    # Optim AutoHess 
    res = optimProbAutoHess(p0, showTrace=false)
    sqDiff = sum((res.minimizer - peTabOpt.paramVecNotTransformed).^2)
    if sqDiff > 1e-3
        @printf("sqDiffHessApprox = %.3e\n", sqDiff)
        @printf("Failed on optimization for Optim AutoHess\n")
        return false
    else
        @printf("Passed test for Optim interior point newton hessian autoDiff\n")
    end

    # Optim BFGS
    res = optimProbBFGS(p0, showTrace=false)
    sqDiff = sum((res.minimizer - peTabOpt.paramVecNotTransformed).^2)
    if sqDiff > 1.0
        @printf("sqDiffBFGS = %.3e\n", sqDiff)
        @printf("Failed on optimization for Optim BFGS\n")
        return false
    else
        @printf("Passed test for Optim BFGS\n")
    end

    # Optim BFGS
    res = optimProbLBFGS(p0, showTrace=false)
    sqDiff = sum((res.minimizer - peTabOpt.paramVecNotTransformed).^2)
    if sqDiff > 1.0
        @printf("sqDiffLBFGS = %.3e\n", sqDiff)
        @printf("Failed on optimization for Optim LBFGS\n")
        return false
    else
        @printf("Passed test for Optim LBFGS\n")
    end

    # NLopt LBFGS
    minf, minx, ret = NLopt.optimize(NLoptLBFGS, p0)
    sqDiff = sum((minx - peTabOpt.paramVecNotTransformed).^2)
    if sqDiff > 1.0
        @printf("sqDiffLBFGS = %.3e\n", sqDiff)
        @printf("Failed on optimization for NLopt LBFGS\n")
        return false
    else
        @printf("Passed test for NLopt LBFGS\n")
    end

    # NLopt truncated Newton with preconditioner 
    minf, minx, ret = NLopt.optimize(NLoptTNewton, p0)
    sqDiff = sum((minx - peTabOpt.paramVecNotTransformed).^2)
    if sqDiff > 1e-3
        @printf("sqDiffLBFGS = %.3e\n", sqDiff)
        @printf("Failed on optimization for NLopt Newton with preconditioner \n")
        return false
    else
        @printf("Passed test for NLopt Newton with preconditioner \n")
    end

    p0 = peTabOpt.paramVecNotTransformed .+ 0.1
    stepLengths = log10.(LinRange(exp10(1e-2),exp10(1e-4), 10000))
    forwardGradOpt = createFowardGradientProb(peTabOpt, stepLengths, 10000)
    minCost, minParams, costVal, paramVal = runAdam(p0, forwardGradOpt, verbose=false, seed=14)
    # The dynamic parameters are easily estimated, but sd-parameter often fail so do not 
    # include them into the tests.
    sqDiff = sum((minParams[1:2] - peTabOpt.paramVecNotTransformed[1:2]).^2)
    if sqDiff > 1e-3
        @printf("sqDiffForwardGradient = %.3e\n", sqDiff)
        @printf("Failed optimization for forward gradient\n")
        return false
    else
        @printf("Passed test for optimization forward gradient\n")
    end

    return true
end


peTabModel = setUpPeTabModel("Test_model2", pwd() * "/tests/Test_model2/", forceBuildJlFile=true)

passTest = testOdeSol(peTabModel, Vern9(), 1e-9, printRes=true)
if passTest == true
    @printf("Passed test for ODE solution\n")
else
    @printf("Did not pass test for ODE solution\n")
end

passTest = testCostGradHess(peTabModel, Vern9(), 1e-15, printRes=true)
if passTest == true
    @printf("Passed test for cost, gradient and hessian\n")
else
    @printf("Did not pass test for cost, gradient and hessian\n")
end

# Will have to change parameters 
loadFidesFromPython("/home/sebpe/anaconda3/envs/PeTab/bin/python")
passTest = testOptimizer(peTabModel, Rodas4P(), 1e-9)
if passTest == true
    @printf("Passed test for checking optimizers\n")
else
    @printf("Did not pass test for checking optimizers\n")
end
