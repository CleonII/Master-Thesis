# Check the accruacy of the PeTab importer for a simple linear ODE;
#   x' = Ax; x(0) = u0 and A = [alpha beta ; gamma delta].
# This ODE is solved analytically, and using the analytical solution the accuracy of 
# the ODE solver, cost function, gradient and hessian of the PeTab importer is checked.
# The measurment data is avaible in tests/Test_model1/


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


"""
    solveOde2x2Lin(t, u0, alpha, beta, gamma, delta)

    Solve ODE system x' = Ax with A = [alpha beta ; gamma delta] a time point 
    t given initial values x(0)=u0. 
"""
function solveOde2x2Lin(t, u0, alpha, beta, gamma, delta)

    # Calculate eigenvectors of A (PQ-formula)
    m = (alpha + delta) / 2
    p = (alpha*delta) - (beta*gamma)
    eVal1 = m + sqrt(m^2 - p)
    eVal2 = m - sqrt(m^2 - p)

    # Eigenvectors of A in matrix 
    e1, e3 = [ beta, -(alpha-eVal1)]
    e2, e4 = [ beta, -(alpha-eVal2)]
    eVecMat = [e1 e2; e3 e4]

    # Solve for constants c1 and c2 
    c1, c2 = eVecMat \ u0

    # ODE solution at time t
    return c1*exp(eVal1*t)*eVecMat[:, 1] + c2*exp(eVal2*t)*eVecMat[:, 2]
end


"""
    testOdeSol(solver, tol; printRes=false)

    Compare analytical vs numeric ODE solver using a provided solver with 
    tolerance tol for the Test_model1.
    Returns true if passes test (sqDiff less than 1e-8) else returns false. 
"""
function testOdeSol(peTabModel::PeTabModel, solver, tol; printRes=false)

    # Set values to PeTab file values 
    experimentalConditionsFile, measurementDataFile, parameterDataFile, observablesDataFile = readDataFiles(peTabModel.dirModel, readObs=true)
    measurementData = processMeasurementData(measurementDataFile, observablesDataFile) 
    paramData = processParameterData(parameterDataFile) 
    setParamToFileValues!(peTabModel.paramMap, peTabModel.stateMap, paramData)
    
    # Extract experimental conditions for simulations 
    simulationInfo = getSimulationInfo(measurementDataFile, measurementData)

    # Parameter values where to teast accuracy. Each column is a alpha, beta, gamma and delta
    u0 = [8.0, 4.0]
    paramMat = reshape([2.0, 3.0, 3.0, 5.0, 
                        1.0, 2.0, 3.0, 3.0,
                        1.0, 0.4, 0.3, 0.5,
                        4.0, 3.0, 2.0, 1.0,
                        0.01, 0.02, 0.03, 0.04], (4, 5))

    for i in 1:5

        alpha, beta, gamma, delta = paramMat[:, i]
        
        # Set parameter values for ODE
        peTabModel.paramMap[2] = Pair(peTabModel.paramMap[2].first, alpha)
        peTabModel.paramMap[4] = Pair(peTabModel.paramMap[4].first, gamma)
        peTabModel.paramMap[5] = Pair(peTabModel.paramMap[5].first, delta)
        peTabModel.paramMap[6] = Pair(peTabModel.paramMap[6].first, beta)
        prob = ODEProblem(peTabModel.odeSystem, peTabModel.stateMap, (0.0, 5e3), peTabModel.paramMap, jac=true)
        prob = remake(prob, p = convert.(Float64, prob.p), u0 = convert.(Float64, prob.u0))
        changeToExperimentalCondUse! = (pVec, u0Vec, expID) -> changeExperimentalCond!(pVec, u0Vec, expID, paramData, experimentalConditionsFile, peTabModel)
        
        # Solve ODE system 
        solArray, success = solveOdeModelAllExperimentalCond(prob, changeToExperimentalCondUse!, measurementDataFile, simulationInfo, solver, tol)
        solNumeric = solArray[1]
        
        # Compare against analytical solution 
        sqDiff = 0.0
        for t in solNumeric.t
            solAnalytic = solveOde2x2Lin(t, u0, alpha*0.5, beta, gamma, delta)
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
    alpha, beta, gamma, delta = paramVec[1:4]
    measurementData = CSV.read(pwd() * "/tests/Test_model1/measurementData_Test_model1.tsv", DataFrame)

    # Extract correct parameter for observation i and compute logLik
    logLik = 0.0
    for i in 1:nrow(measurementData)

        # Specs for observation i
        obsID = measurementData[i, :observableId]
        noiseID = measurementData[i, :noiseParameters]
        yObs = measurementData[i, :measurement]
        t = measurementData[i, :time]
        expId = measurementData[i, :simulationConditionId]

        # Extract correct sigma 
        if noiseID == "sd_sebastian_new"
            sigma = paramVec[5]
        elseif noiseID == "sd_damiano_new"
            sigma = paramVec[6]
        end

        # Consider experimental conditions 
        if expId == "model1_data1"
            alpha_use = alpha * 0.5
        else
            alpha_use = alpha * 1.0
        end

        sol = solveOde2x2Lin(t, u0, alpha_use, beta, gamma, delta)
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
    parameter vectors for Test_model1. For the analytical solution the gradient 
    and hessian are computed via ForwardDiff.
"""
function testCostGradHess(peTabModel::PeTabModel, solver, tol; printRes::Bool=false)

    peTabOpt = setUpCostGradHess(peTabModel, solver, tol)

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
        if sqDiffGrad > 1e-6
            @printf("sqDiffGrad = %.3e\n", sqDiffGrad)
            @printf("Does not pass test on gradient\n")
            return false
        end

        # Evaluate hessian 
        hessAnalytic = ForwardDiff.hessian(calcCostAnalytic, paramVec)
        hessNumeric = zeros(nParam, nParam); peTabOpt.evalHess(hessNumeric, paramVec)
        sqDiffHess = sum((hessAnalytic - hessNumeric).^2)
        if sqDiffHess > 1e-4
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


peTabModel = setUpPeTabModel("Test_model1", pwd() * "/tests/Test_model1/", forceBuildJlFile=true)

passTest = testOdeSol(peTabModel, Vern9(), 1e-9, printRes=false)
if passTest == true
    @printf("Passed test for ODE solution\n")
else
    @printf("Did not pass test for ODE solution\n")
end

passTest = testCostGradHess(peTabModel, Vern9(), 1e-12, printRes=true)
if passTest == true
    @printf("Passed test for cost, gradient and hessian\n")
else
    @printf("Did not pass test for cost, gradient and hessian\n")
end
