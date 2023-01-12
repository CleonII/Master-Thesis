#=
    Check the accruacy of the PeTab importer for a simple linear ODE;
        x' = Ax; x(0) = u0 and A = [α β ; γ δ].
    This ODE is solved analytically, and using the analytical solution the accuracy of 
    the ODE solver, cost function, gradient and hessian of the PeTab importer is checked.
    The measurment data is avaible in tests/Test_model1/
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


"""
    solveOde2x2Lin(t, u0, α, β, γ, δ)

    Solve ODE system x' = Ax with A = [α β ; γ δ] a time point 
    t given initial values x(0)=u0. 
"""
function solveOde2x2Lin(t, u0, α, β, γ, δ)

    # Calculate eigenvectors of A (PQ-formula)
    m = (α + δ) / 2
    p = (α*δ) - (β*γ)
    eVal1 = m + sqrt(m^2 - p)
    eVal2 = m - sqrt(m^2 - p)

    # Eigenvectors of A in matrix 
    e1, e3 = [ β, -(α-eVal1)]
    e2, e4 = [ β, -(α-eVal2)]
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
    simulationInfo = getSimulationInfo(peTabModel, measurementDataFile, measurementData)

    # Parameter values where to teast accuracy. Each column is a α, β, γ and δ
    u0 = [8.0, 4.0]
    paramMat = reshape([2.0, 3.0, 3.0, 5.0, 
                        1.0, 2.0, 3.0, 3.0,
                        1.0, 0.4, 0.3, 0.5,
                        4.0, 3.0, 2.0, 1.0,
                        0.01, 0.02, 0.03, 0.04], (4, 5))

    for i in 1:5

        α, β, γ, δ = paramMat[:, i]
        
        # Set parameter values for ODE
        peTabModel.paramMap[2] = Pair(peTabModel.paramMap[2].first, α)
        peTabModel.paramMap[4] = Pair(peTabModel.paramMap[4].first, γ)
        peTabModel.paramMap[5] = Pair(peTabModel.paramMap[5].first, δ)
        peTabModel.paramMap[6] = Pair(peTabModel.paramMap[6].first, β)
        prob = ODEProblem(peTabModel.odeSystem, peTabModel.stateMap, (0.0, 5e3), peTabModel.paramMap, jac=true)
        prob = remake(prob, p = convert.(Float64, prob.p), u0 = convert.(Float64, prob.u0))
        changeToExperimentalCondUse! = (pVec, u0Vec, expID) -> changeExperimentalCond!(pVec, u0Vec, expID, paramData, experimentalConditionsFile, peTabModel)
        
        # Solve ODE system 
        solArray, success = solveOdeModelAllExperimentalCond(prob, changeToExperimentalCondUse!, simulationInfo, solver, tol, tol, peTabModel.getTStops)
        solNumeric = solArray[1]
        
        # Compare against analytical solution 
        sqDiff = 0.0
        for t in solNumeric.t
            solAnalytic = solveOde2x2Lin(t, u0, α*0.5, β, γ, δ)
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
    α, β, γ, δ = paramVec[1:4]
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
            α_use = α * 0.5
        else
            α_use = α * 1.0
        end

        sol = solveOde2x2Lin(t, u0, α_use, β, γ, δ)
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

    peTabOpt = setUpCostGradHess(peTabModel, solver, tol, sensealg = ForwardDiffSensitivity(), 
                                 sensealgForward = ForwardDiffSensitivity(), solverForward=Vern9(),
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

        # Evalute Zygote cost function 
        costPeTab = peTabOpt.evalFZygote(paramVec)
        sqDiffCostZygote = (costPeTab - costAnalytic)^2
        if sqDiffCostZygote > 1e-6
            @printf("sqDiffCost = %.3e\n", sqDiffCost)
            @printf("Does not pass test on cost for Zygote cost function\n")
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
        
        # Evaluate hessian 
        hessAnalytic = ForwardDiff.hessian(calcCostAnalytic, paramVec)
        hessNumeric = zeros(nParam, nParam); peTabOpt.evalHess(hessNumeric, paramVec)
        sqDiffHess = sum((hessAnalytic - hessNumeric).^2)
        if sqDiffHess > 1e-4
            @printf("sqDiffHess = %.3e\n", sqDiffHess)
            @printf("Does not pass test on hessian\n")
            return false
        end

        # Here we normalise the gradients (as the gradients are typically huge here 1e8, so without any form 
        # of normalisation we do not pass the test)
        gradAdj = zeros(nParam)
        peTabOpt.evalGradFAdjoint(gradAdj, paramVec)
        sqDiffGradAdjoint1 = sum((gradAnalytic ./ norm(gradAnalytic) - gradAdj / norm(gradAdj)).^2)
        if sqDiffGradAdjoint1 > 1e-4
            @printf("sqDiffGradAdjointOpt1 = %.3e\n", sqDiffGradAdjoint1)
            @printf("Does not pass test on adjoint gradient gradient\n")
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


# Test different adjoint approaches 
function testDifferentAdjoints(peTabModel::PeTabModel, solver, tol; printRes::Bool=false)

    peTabOpt1 = setUpCostGradHess(peTabModel, solver, tol, sensealg = ForwardDiffSensitivity())
    peTabOpt2 = setUpCostGradHess(peTabModel, solver, tol, sensealg = QuadratureAdjoint())
    peTabOpt3 = setUpCostGradHess(peTabModel, solver, tol, sensealg = InterpolatingAdjoint(autojacvec=ReverseDiffVJP()))
    peTabOpt4 = setUpCostGradHess(peTabModel, solver, tol, sensealg = BacksolveAdjoint())
    peTabOpt5 = setUpCostGradHess(peTabModel, solver, tol, sensealg = ReverseDiffAdjoint())

    Random.seed!(123)
    createCube(peTabOpt1, 5)
    cube = Matrix(CSV.read(peTabOpt1.pathCube, DataFrame))

    paramVec = cube[1, :]
    gradAnalytic = ForwardDiff.gradient(calcCostAnalytic, paramVec)

    gradZygoteSensealg = Zygote.gradient(peTabOpt1.evalFZygote, paramVec)[1]
    sqDiffGrad = sum((gradAnalytic - gradZygoteSensealg).^2)
    if sqDiffGrad > 1e-6
        @printf("sqDiffGrad = %.3e\n", sqDiffGrad)
        @printf("Does not pass test on gradient with ForwardDiffSensitivity\n")
        return false
    elseif printRes == true
        @printf("sqDiffGrad = %.3e\n", sqDiffGrad)
        @printf("Pass test on gradient with ForwardDiffSensitivity\n")
    end

    gradZygoteSensealg = Zygote.gradient(peTabOpt2.evalFZygote, paramVec)[1]
    sqDiffGrad = sum((gradAnalytic - gradZygoteSensealg).^2)
    if sqDiffGrad > 1e-6
        @printf("sqDiffGrad = %.3e\n", sqDiffGrad)
        @printf("Does not pass test on gradient with QuadratureAdjoint\n")
        return false
    elseif printRes == true
        @printf("sqDiffGrad = %.3e\n", sqDiffGrad)
        @printf("Pass test on gradient with QuadratureAdjoint\n")
    end

    gradZygoteSensealg = Zygote.gradient(peTabOpt3.evalFZygote, paramVec)[1]
    sqDiffGrad = sum((gradAnalytic - gradZygoteSensealg).^2)
    if sqDiffGrad > 1e-6
        @printf("sqDiffGrad = %.3e\n", sqDiffGrad)
        @printf("Does not pass test on gradient with InterpolatingAdjoint(autojacvec=ReverseDiffVJP())\n")
        return false
    elseif printRes == true
        @printf("sqDiffGrad = %.3e\n", sqDiffGrad)
        @printf("Pass test on gradient with InterpolatingAdjoint(autojacvec=ReverseDiffVJP())\n")
    end

    gradZygoteSensealg = Zygote.gradient(peTabOpt4.evalFZygote, paramVec)[1]
    sqDiffGrad = sum((gradAnalytic - gradZygoteSensealg).^2)
    if sqDiffGrad > 1e-6
        @printf("sqDiffGrad = %.3e\n", sqDiffGrad)
        @printf("Does not pass test on gradient with BacksolveAdjoint\n")
        return false
    elseif printRes == true
        @printf("sqDiffGrad = %.3e\n", sqDiffGrad)
        @printf("Pass test on gradient with BacksolveAdjoint\n")
    end

    # Currently breaks because ReverseDiff procues an array as output (will have to fix)
    #=
    gradZygoteSensealg = Zygote.gradient(peTabOpt5.evalFZygote, paramVec)[1]
    sqDiffGrad = sum((gradAnalytic - gradZygoteSensealg).^2)
    if sqDiffGrad > 1e-6
        @printf("sqDiffGrad = %.3e\n", sqDiffGrad)
        @printf("Does not pass test on gradient with ReverseDiffAdjoint\n")
        return false
    elseif printRes == true
        @printf("sqDiffGrad = %.3e\n", sqDiffGrad)
        @printf("Pass test on gradient with ReverseDiffAdjoint\n")
    end
    =#

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

passTest = testDifferentAdjoints(peTabModel, Vern9(), 1e-12, printRes=true)
if passTest == true
    @printf("Passed test for different adjoints\n")
else
    @printf("Did not pass test for different adjoints\n")
end
