# Check the accruacy of the PeTab importer for a simple linear ODE;
#   x' = a - bx + cy; x(0) = 0 
#   y' = bx - cy - dy;  y(0) = 0
# where the model has a pre-equilibrium condition. That is he simulated data for 
# this ODE model is generated by starting from the steady state;
# x* = a / b + ( a * c ) / ( b * d )
# y* = a / d
# and when computing the cost in the PeTab importer the model is first simualted 
# to a steady state, and then the mian simulation matched against data is 
# performed.
# This test compares the ODE-solution, cost, gradient and hessian when 
# i) solving the ODE using the SS-equations as initial condition, and ii) when 
# first simulating the model to the steady state.
# Accruacy of both the hessian and gradient are strongly dependent on the tolerances 
# used in the TerminateSteadyState used to terminate steady state simulation.


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

# PeTab importer to get cost, grad etc 
include(joinpath(pwd(), "src", "Process_PEtab_files", "Create_cost_grad_hessian.jl"))

# HyperCube sampling 
include(joinpath(pwd(), "src", "Optimizers", "Lathin_hypercube.jl"))

# For converting to SBML 
include(joinpath(pwd(), "src", "SBML", "SBML_to_ModellingToolkit.jl"))


function getSolAlgebraicSS(peTabModel::PeTabModel, solver, tol::Float64, a::T1, b::T1, c::T1, d::T1) where T1<:Real

    # ODE solution with algebraically computed initial values (instead of ss pre-simulation)
    odeProb = ODEProblem(peTabModel.odeSystem, peTabModel.stateMap, (0.0, 9.7), peTabModel.paramMap, jac=true)
    odeProb = remake(odeProb, p = convert.(eltype(a), odeProb.p), u0 = convert.(eltype(a), odeProb.u0))
    solArray = Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1}(undef, 2)

    # Set model parameter values to ensure initial steady state 
    odeProb.p[5], odeProb.p[3], odeProb.p[1], odeProb.p[6] = a, b, c, d
    odeProb.u0[1] = a / b + ( a * c ) / ( b * d ) # x0
    odeProb.u0[2] = a / d # y0

    odeProb.p[4] = 2.0 # a_scale
    solArray[1] = solve(odeProb, solver, abstol=tol, reltol=tol)
    odeProb.p[4] = 0.5 # a_scale
    solArray[2] = solve(odeProb, solver, abstol=tol, reltol=tol)

    return solArray
end


function calcCostAlgebraic(paramVec, peTabModel, solver, tol)
    
    a, b, c, d = paramVec

    experimentalConditionsFile, measurementDataFile, parameterDataFile, observablesDataFile = readDataFiles(peTabModel.dirModel, readObs=true)
    measurementData = processMeasurementData(measurementDataFile, observablesDataFile) 

    solArrayAlg = getSolAlgebraicSS(peTabModel, solver, tol, a, b, c, d)
    logLik = 0.0
    for i in eachindex(measurementData.tObs)
        yObs = measurementData.yObsNotTransformed[i]
        t = measurementData.tObs[i]
        if measurementData.conditionId[i] == "defaultdouble"
            yMod = solArrayAlg[1](t)[1]
        else
            yMod = solArrayAlg[2](t)[2]
        end
        sigma = 0.04
        logLik += log(sigma) + 0.5*log(2*pi) + 0.5 * ((yObs - yMod) / sigma)^2
    end
    
    return logLik
end


function testOdeSol(peTabModel::PeTabModel, solver, tol; printRes=false)
   
    # Set values to PeTab file values 
    experimentalConditionsFile, measurementDataFile, parameterDataFile, observablesDataFile = readDataFiles(peTabModel.dirModel, readObs=true)
    measurementData = processMeasurementData(measurementDataFile, observablesDataFile) 
    paramData = processParameterData(parameterDataFile) 
    setParamToFileValues!(peTabModel.paramMap, peTabModel.stateMap, paramData)
    
    # Extract experimental conditions for simulations 
    simulationInfo = getSimulationInfo(peTabModel, measurementDataFile, measurementData, absTolSS=1e-12, relTolSS=1e-10)

    # Parameter values where to teast accuracy. Each column is a alpha, beta, gamma and delta
    # a, b, c, d
    paramMat = reshape([1.0, 2.0, 3.0, 4.0,
                        0.1, 0.2, 0.3, 0.4,
                        4.0, 3.0, 2.0, 1.0,
                        1.0, 1.0, 1.0, 1.0,
                        2.5, 7.0, 3.0, 3.0,], (4, 5))

    for i in 1:5    
        a, b, c, d = paramMat[:, i]        
        # Set parameter values for ODE
        peTabModel.paramMap[1] = Pair(peTabModel.paramMap[1].first, c)
        peTabModel.paramMap[3] = Pair(peTabModel.paramMap[3].first, b)
        peTabModel.paramMap[5] = Pair(peTabModel.paramMap[5].first, a)
        peTabModel.paramMap[6] = Pair(peTabModel.paramMap[6].first, d)

        prob = ODEProblem(peTabModel.odeSystem, peTabModel.stateMap, (0.0, 9.7), peTabModel.paramMap, jac=true)
        prob = remake(prob, p = convert.(Float64, prob.p), u0 = convert.(Float64, prob.u0))
        changeToExperimentalCondUse! = (pVec, u0Vec, expID) -> changeExperimentalCond!(pVec, u0Vec, expID, paramData, experimentalConditionsFile, peTabModel)
        
        # Solve ODE system with steady state simulation 
        solArray, success = solveOdeModelAllExperimentalCond(prob, changeToExperimentalCondUse!, simulationInfo, solver, tol, tol, peTabModel.getTStops)

        # Solve ODE system with algebraic intial values 
        solArrayAlg = getSolAlgebraicSS(peTabModel, solver, tol, a, b, c, d)
        
        # Compare against analytical solution 
        sqDiff = 0.0
        for i in eachindex(solArray)
            solNum = solArray[i]
            solAlg = solArrayAlg[i]
            sqDiff += sum((Array(solNum) - Array(solAlg(solNum.t))).^2)
        end

        if sqDiff > 1e-6
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


function testCostGradHess(peTabModel::PeTabModel, solver, tol; printRes::Bool=false)

    solver = Rodas4P(autodiff=false)
    tol = 1e-12
    peTabOpt = setUpCostGradHess(peTabModel, solver, tol, absTolSS=1e-12, relTolSS=1e-10, sensealg=ForwardDiffSensitivity(), 
                                 sensealgForward = ForwardSensitivity(), solverForward=CVODE_BDF(),
                                 adjTol=1e-12, adjSensealg=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(false)))
    peTabOptAdj = setUpCostGradHess(peTabModel, solver, tol, absTolSS=1e-12, relTolSS=1e-10, sensealg=ForwardDiffSensitivity(), 
                                    adjTol=1e-12, adjSensealg=adjSensealg=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(false)),
                                    adjSensealgSS=QuadratureAdjoint(autojacvec=ReverseDiffVJP(false)))   
    peTabOptAlt = setUpCostGradHess(peTabModel, solver, tol, 
                                    sensealgForward = :AutoDiffForward, solverForward=solver)                                                                                                   
    calcCostAlg = (paramVec) -> calcCostAlgebraic(paramVec, peTabModel, solver, tol)

    Random.seed!(123)
    createCube(peTabOpt, 5)
    cube = Matrix(CSV.read(peTabOpt.pathCube, DataFrame))

    for i in 1:5

        println("i = $i")
        paramVec = cube[i, :]

        # Evaluate cost 
        costPeTab = peTabOpt.evalF(paramVec)
        costAlg = calcCostAlg(paramVec)
        sqDiffCost = (costPeTab - costAlg)^2
        if sqDiffCost > 1e-4
            @printf("sqDiffCost = %.3e\n", sqDiffCost)
            @printf("Does not pass test on cost\n")
            return false
        end

        # Evaluate cost for Zygote 
        costZygote = peTabOpt.evalFZygote(paramVec)
        sqDiffCostZygote = (costZygote - costAlg)^2
        if sqDiffCostZygote > 1e-4
            @printf("sqDiffCostZygote = %.3e\n", sqDiffCostZygote)
            @printf("Does not pass test on cost\n")
            return false
        end

        # Evaluate gradient 
        gradAlg = ForwardDiff.gradient(calcCostAlg, paramVec)
        grad = zeros(4)
        peTabOpt.evalGradF(grad, paramVec)
        sqDiffGrad = sum((grad - gradAlg).^2)
        if sqDiffGrad > 1e-4
            @printf("sqDiffGrad = %.3e\n", sqDiffGrad)
            @printf("Does not pass test on gradient\n")
            return false
        end

        # Evaluate gradient via Zygote 
        gradZygote = Zygote.gradient(peTabOpt.evalFZygote, paramVec)[1]
        sqDiffGradZygote = sum((gradZygote - gradAlg).^2)
        if sqDiffGradZygote > 1e-4
            @printf("sqDiffGradZygote = %.3e\n", sqDiffGradZygote)
            @printf("Does not pass test on gradient\n")
            return false
        end

        # Forward sensitivity equations gradient
        gradForwardEq = zeros(4)
        peTabOpt.evalGradFForwardEq(gradForwardEq, paramVec)
        sqDiffGradForwardEq = sum((gradForwardEq - gradAlg).^2)
        if sqDiffGradForwardEq > 1e-5
            @printf("sqDiffGradForwardEq = %.3e\n", sqDiffGradForwardEq)
            @printf("Does not pass test on gradient from Forward sensitivity equations\n")
            return false
        end

        # Forward sensitivity equations using autodiff 
        gradForwardEqAuto = zeros(4)
        peTabOptAlt.evalGradFForwardEq(gradForwardEqAuto, paramVec)
        sqDiffGradForwardEqAuto = sum((gradForwardEq - gradAlg).^2)
        if sqDiffGradForwardEqAuto > 1e-5
            @printf("sqDiffGradForwardEqAuto = %.3e\n", sqDiffGradForwardEqAuto)
            @printf("Does not pass test on gradient from Forward sensitivity equations\n")
            return false
        end

        # Evaluate lower level adjoint sensitivity interfance gradient 
        gradAdj = zeros(4)
        peTabOpt.evalGradFAdjoint(gradAdj, paramVec)
        sqDiffGradAdjoint1 = sum((gradAdj - gradAlg).^2)
        if sqDiffGradAdjoint1 > 1e-4
            @printf("sqDiffGradAdjointOpt1 = %.3e\n", sqDiffGradAdjoint1)
            @printf("Does not pass test on adjoint gradient gradient\n")
            return false
        end

        # Evaluate lower level adjoint sensitivity interfance gradient 
        gradAdj = zeros(4)
        peTabOptAdj.evalGradFAdjoint(gradAdj, paramVec)
        sqDiffGradAdjoint2 = sum((gradAdj - gradAlg).^2)
        if sqDiffGradAdjoint2 > 1e-4
            @printf("sqDiffGradAdjointOpt2 = %.3e\n", sqDiffGradAdjoint2)
            @printf("Does not pass test on adjoint without steady state adjoint approach\n")
            return false
        end

        # Evaluate hessian 
        hessAlg = ForwardDiff.hessian(calcCostAlg, paramVec)
        hess = zeros(4, 4); peTabOpt.evalHess(hess, paramVec)
        sqDiffHess = sum((hessAlg - hess).^2)
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
            @printf("sqDiffGradAdjointOpt1 = %.3e\n", sqDiffGradAdjoint1)
            @printf("sqDiffGradAdjointOpt2 = %.3e\n", sqDiffGradAdjoint2)
            @printf("sqDiffGradForwardEq = %.3e\n", sqDiffGradForwardEq)
            @printf("sqDiffGradForwardEqAuto = %.3e\n", sqDiffGradForwardEqAuto)
        end
    end

    return true 
end


peTabModel = setUpPeTabModel("Test_model3", pwd() * "/tests/Test_model3/", forceBuildJlFile=true)
passTest = testOdeSol(peTabModel, Rodas4P(), 1e-12, printRes=true)
if passTest == true
    @printf("Passed test for ODE solution\n")
else
    @printf("Did not pass test for ODE solution\n")
end

passTest = testCostGradHess(peTabModel, Rodas4P(), 1e-12, printRes=true)
if passTest == true
    @printf("Passed test for cost, gradient and hessian\n")
else
    @printf("Did not pass test for cost, gradient and hessian\n")
end

passTest = checkGradientResiduals(peTabModel, Rodas5(), 1e-9)
if passTest == true
    @printf("Passed test for gradient for residuals\n")
else
    @printf("Did not pass test for gradient for residuals\n")
end