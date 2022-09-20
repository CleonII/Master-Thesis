# Check the accruacy of the PeTab importer for a simple linear ODE;
#   s' = alpha*s; s(0) = 8.0 -> s(t) = 8.0 * exp(alpha*t)
#   d' = beta*d;  d(0) = 4.0 -> d(t) = 4.0 * exp(beta*t)
# This ODE is solved analytically, and using the analytical solution the accuracy of 
# the ODE solver, cost function, gradient and hessian of the PeTab importer is checked.
# The accuracy of the optimizers is further checked.


using ModelingToolkit, DifferentialEquations, DataFrames, CSV 
using ForwardDiff
using ReverseDiff
using StatsBase
using Random
using LinearAlgebra
using Distributions
using Printf
using Ipopt
using Optim


# Functions for solving ODE system 
include(joinpath(pwd(), "Additional_functions", "benchmarkSolvers.jl"))

# Ipopt wrapper 
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "CommonParameterEstimationMethods", "CreateIpoptProb.jl"))

# PeTab importer 
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "PeTabImporter.jl"))

# Additional functions for ODE solver 
include(joinpath(pwd(), "Additional_functions", "additional_tools.jl"))

# HyperCube sampling 
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "LatinHyperCubeParameters", "LatinHyperCubeSampledParameters.jl"))

# For converting to SBML 
include(joinpath(pwd(), "Pipeline_SBMLImporter", "main.jl"))

# For building observable 
include(joinpath(pwd(), "Additional_functions", "createObsFile.jl"))


"""
    testOdeSol(solver, tol; printRes=false)

    Compare analytical vs numeric ODE solver using a provided solver with 
    tolerance tol for the Test_model2.
    Returns true if passes test (sqDiff less than 1e-8) else returns false. 
"""
function testOdeSol(peTabModel::PeTabModel, solver, tol; printRes=false)
   
    # Set values to PeTab file values 
    experimentalConditionsFile, measurementDataFile, parameterBoundsFile = readDataFiles(peTabModel.dirModel)
    paramData = processParameterData(parameterBoundsFile) 
    setParamToParamFileVal!(peTabModel.paramMap, peTabModel.stateMap, paramData)
    
    # Extract experimental conditions for simulations 
    simulationInfo = getSimulationInfo(measurementDataFile, peTabModel.odeSystem)

    # Parameter values where to teast accuracy. Each column is a alpha, beta, gamma and delta
    u0 = [8.0, 4.0]
    paramMat = reshape([2.0, 3.0,
                        1.0, 2.0,
                        1.0, 0.4,
                        4.0, 3.0,
                        0.01, 0.02], (2, 5))

    for i in 1:5
        i = 1
        alpha, beta = paramMat[:, i]        
        # Set parameter values for ODE
        peTabModel.paramMap[2] = Pair(peTabModel.paramMap[2].first, alpha)
        peTabModel.paramMap[3] = Pair(peTabModel.paramMap[3].first, beta)
        prob = ODEProblem(peTabModel.odeSystem, peTabModel.stateMap, (0.0, 5e3), peTabModel.paramMap, jac=true)
        prob = remake(prob, p = convert.(Float64, prob.p), u0 = convert.(Float64, prob.u0))
        changeToExperimentalCondUse! = (pVec, u0Vec, expID) -> changeExperimentalCond!(pVec, u0Vec, expID, paramData, experimentalConditionsFile, peTabModel)
        
        # Solve ODE system 
        solArray, success = solveOdeModelAllCond(prob, changeToExperimentalCondUse!, measurementDataFile, simulationInfo, solver, tol)
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
    
    evalF, evalGradF, evalHessianApproxF, paramVecEstTmp, lowerBounds, upperBounds, idParam = setUpCostFunc(peTabModel, solver, tol)
    # "Exact" hessian via autodiff 
    evalH = (hessianMat, paramVec) -> begin hessianMat .= Symmetric(ForwardDiff.hessian(evalF, paramVec)) end

    Random.seed!(123)
    fileSaveCube = pwd() * "/tests/Test_model2/CubeTestModel2.csv"
    createCube(5, lowerBounds, upperBounds, fileSaveCube, evalF)
    cube = Matrix(CSV.read(fileSaveCube, DataFrame))
    nParam = size(cube)[2]

    for i in 1:5
        i = 1
        paramVec = cube[i, :]

        # Evaluate cost 
        costPeTab = evalF(paramVec)
        costAnalytic = calcCostAnalytic(paramVec)
        sqDiffCost = (costPeTab - costAnalytic)^2
        if sqDiffCost > 1e-6
            @printf("sqDiffCost = %.3e\n", sqDiffCost)
            @printf("Does not pass test on cost\n")
            return false
        end

        # Evaluate gradient 
        gradAnalytic = ForwardDiff.gradient(calcCostAnalytic, paramVec)
        gradNumeric = zeros(nParam); evalGradF(paramVec, gradNumeric)
        sqDiffGrad = sum((gradAnalytic - gradNumeric).^2)
        if sqDiffGrad > 1e-4
            @printf("sqDiffGrad = %.3e\n", sqDiffGrad)
            @printf("Does not pass test on gradient\n")
            return false
        end

        # Evaluate hessian 
        hessAnalytic = ForwardDiff.hessian(calcCostAnalytic, paramVec)
        hessNumeric = zeros(nParam, nParam); evalH(hessNumeric, paramVec)
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

    evalF, evalGradF, evalHessianApproxF, paramVecEstTmp, lowerBounds, upperBounds, idParam = setUpCostFunc(peTabModel, solver, tol)
    # "Exact" hessian via autodiff 
    evalH = (hessianMat, paramVec) -> begin hessianMat .= Symmetric(ForwardDiff.hessian(evalF, paramVec)) end

    Random.seed!(123)
    fileSaveCube = pwd() * "/tests/Test_model2/CubeTestModel2.csv"
    createCube(5, lowerBounds, upperBounds, fileSaveCube, evalF)
    cube = Matrix(CSV.read(fileSaveCube, DataFrame))
    nParam = size(cube)[2]

    # Ipopt optimizers 
    ipoptProbHessApprox, iterArrHessApprox = createIpoptProbNew(evalF, evalGradF, evalHessianApproxF, lowerBounds, upperBounds, emptyH=false) 
    ipoptProbBfgs, iterArrBfgs = createIpoptProbNew(evalF, evalGradF, evalHessianApproxF, lowerBounds, upperBounds, emptyH=true) 
    ipoptProbAutoHess, iterArrAutoHess = createIpoptProbNew(evalF, evalGradF, evalH, lowerBounds, upperBounds, emptyH=false) 

    # Optim optimizers 
    optimProbHessApprox = createOptimProb(evalF, evalGradF, evalHessianApproxF, lowerBounds, upperBounds, showTrace=false)
    optimProbAutoHess = createOptimProb(evalF, evalGradF, evalH, lowerBounds, upperBounds, showTrace=false)

    p0 = cube[1, :]
    # Ipopt Hessian approximation
    ipoptProbHessApprox.x = deepcopy(p0)
    sol_opt = Ipopt.IpoptSolve(ipoptProbHessApprox)
    sqDiff = sum((ipoptProbHessApprox.x - paramVecEstTmp).^2)
    if sqDiff > 1e-3
        @printf("sqDiffHessApprox = %.3e\n", sqDiff)
        @printf("Failed on optimization for Ipopt hessian approximation\n")
        return false
    end
    
    # Ipopt AutoHessian 
    ipoptProbAutoHess.x = deepcopy(p0)
    sol_opt = Ipopt.IpoptSolve(ipoptProbAutoHess)
    sqDiff = sum((ipoptProbAutoHess.x - paramVecEstTmp).^2)
    if sqDiff > 1e-3
        @printf("sqDiffHessApprox = %.3e\n", sqDiff)
        @printf("Failed on optimization for Ipopt auto hessian\n")
        return false
    end

    # Ipopt BFGS 
    ipoptProbBfgs.x = deepcopy(p0)
    sol_opt = Ipopt.IpoptSolve(ipoptProbBfgs)
    sqDiff = sum((ipoptProbBfgs.x - paramVecEstTmp).^2)
    if sqDiff > 1e-3
        @printf("sqDiffHessApprox = %.3e\n", sqDiff)
        @printf("Failed on optimization for Ipopt BFGS\n")
        return false
    end

    # Optim AutoHess 
    res = optimProbAutoHess(p0)
    sqDiff = sum((res.minimizer - paramVecEstTmp).^2)
    if sqDiff > 1e-3
        @printf("sqDiffHessApprox = %.3e\n", sqDiff)
        @printf("Failed on optimization for Optim AutoHess\n")
        return false
    end

    # Optim HessApprox 
    res = optimProbHessApprox(p0)
    sqDiff = sum((res.minimizer - paramVecEstTmp).^2)
    if sqDiff > 1e-3
        @printf("sqDiffHessApprox = %.3e\n", sqDiff)
        @printf("Failed on optimization for HessApprox\n")
        return false
    end

    return true
end


peTabModel = setUpPeTabModel("Test_model2", pwd() * "/tests/Test_model2/")

passTest = testOdeSol(peTabModel, Vern9(), 1e-9, printRes=false)
if passTest == true
    @printf("Passed test for ODE solution\n")
else
    @printf("Did not pass test for ODE solution\n")
end

passTest = testCostGradHess(peTabModel, Vern9(), 1e-15, printRes=false)
if passTest == true
    @printf("Passed test for cost, gradient and hessian\n")
else
    @printf("Did not pass test for cost, gradient and hessian\n")
end

passTest = testOptimizer(peTabModel, Vern9(), 1e-12)
if passTest == true
    @printf("Passed test for checking optimizers\n")
else
    @printf("Did not pass test for checking optimizers\n")
end
