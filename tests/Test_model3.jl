#=
    Check the accruacy of the PeTab importer for a simple linear ODE;
    x' = a - bx + cy; x(0) = 0 
    y' = bx - cy - dy;  y(0) = 0
    where the model has a pre-equilibrium condition. That is he simulated data for 
    this ODE model is generated by starting from the steady state;
    x* = a / b + ( a * c ) / ( b * d )
    y* = a / d
    and when computing the cost in the PeTab importer the model is first simualted 
    to a steady state, and then the mian simulation matched against data is 
    performed.
    This test compares the ODE-solution, cost, gradient and hessian when 
    i) solving the ODE using the SS-equations as initial condition, and ii) when 
    first simulating the model to the steady state.
    Accruacy of both the hessian and gradient are strongly dependent on the tolerances 
    used in the TerminateSteadyState callback.
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

# PeTab importer to get cost, grad etc 
include(joinpath(pwd(), "src", "Process_PEtab_files", "Create_cost_grad_hessian.jl"))

# HyperCube sampling 
include(joinpath(pwd(), "src", "Optimizers", "Lathin_hypercube.jl"))

# For converting to SBML 
include(joinpath(pwd(), "src", "SBML", "SBML_to_ModellingToolkit.jl"))

include(joinpath(pwd(), "tests", "Common.jl"))


function getSolAlgebraicSS(petabModel::PEtabModel, solver, tol::Float64, a::T1, b::T1, c::T1, d::T1) where T1<:Real

    # ODE solution with algebraically computed initial values (instead of ss pre-simulation)
    odeProb = ODEProblem(petabModel.odeSystem, petabModel.stateMap, (0.0, 9.7), petabModel.parameterMap, jac=true)
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


function computeCostAlgebraic(paramVec, petabModel, solver, tol)
    
    a, b, c, d = paramVec

    experimentalConditionsFile, measurementDataFile, parameterDataFile, observablesDataFile = readPEtabFiles(petabModel.dirModel, readObservables=true)
    measurementData = processMeasurements(measurementDataFile, observablesDataFile) 

    solArrayAlg = getSolAlgebraicSS(petabModel, solver, tol, a, b, c, d)
    logLik = 0.0
    for i in eachindex(measurementData.time)
        yObs = measurementData.measurement[i]
        t = measurementData.time[i]
        if measurementData.simulationConditionId[i] == :double
            yMod = solArrayAlg[1](t)[1]
        else
            yMod = solArrayAlg[2](t)[2]
        end
        sigma = 0.04
        logLik += log(sigma) + 0.5*log(2*pi) + 0.5 * ((yObs - yMod) / sigma)^2
    end
    
    return logLik
end


function testODESolverTestModel3(petabModel::PEtabModel, solver, tol)
   
    # Set values to PeTab file values 
    experimentalConditionsFile, measurementDataFile, parameterDataFile, observablesDataFile = readPEtabFiles(petabModel.dirModel, readObservables=true)
    measurementData = processMeasurements(measurementDataFile, observablesDataFile) 
    paramData = processParameters(parameterDataFile) 
    setParamToFileValues!(petabModel.parameterMap, petabModel.stateMap, paramData)
    θ_indices = computeIndicesθ(paramData, measurementData, petabModel.odeSystem, experimentalConditionsFile)
    
    # Extract experimental conditions for simulations 
    simulationInfo = processSimulationInfo(petabModel, measurementData, absTolSS=1e-12, relTolSS=1e-10)

    # Parameter values where to teast accuracy. Each column is a alpha, beta, gamma and delta
    # a, b, c, d
    parametersTest = reshape([1.0, 2.0, 3.0, 4.0,
                              0.1, 0.2, 0.3, 0.4,
                              4.0, 3.0, 2.0, 1.0,
                              1.0, 1.0, 1.0, 1.0,
                              2.5, 7.0, 3.0, 3.0,], (4, 5))

    for i in 1:5    
        a, b, c, d = parametersTest[:, i]        
        # Set parameter values for ODE
        petabModel.parameterMap[1] = Pair(petabModel.parameterMap[1].first, c)
        petabModel.parameterMap[3] = Pair(petabModel.parameterMap[3].first, b)
        petabModel.parameterMap[5] = Pair(petabModel.parameterMap[5].first, a)
        petabModel.parameterMap[6] = Pair(petabModel.parameterMap[6].first, d)

        prob = ODEProblem(petabModel.odeSystem, petabModel.stateMap, (0.0, 9.7), petabModel.parameterMap, jac=true)
        prob = remake(prob, p = convert.(Float64, prob.p), u0 = convert.(Float64, prob.u0))
        θ_est = getFileODEvalues(petabModel)
        changeExperimentalCondition! = (pVec, u0Vec, expID) -> _changeExperimentalCondition!(pVec, u0Vec, expID, θ_est, petabModel, θ_indices)
        
        # Solve ODE system 
        odeSolutions, success = solveODEAllExperimentalConditions(prob, changeExperimentalCondition!, simulationInfo, solver, tol, tol, petabModel.computeTStops)
        # Solve ODE system with algebraic intial values 
        algebraicODESolutions = getSolAlgebraicSS(petabModel, solver, tol, a, b, c, d)
        
        # Compare against analytical solution 
        sqDiff = 0.0
        for i in eachindex(simulationInfo.experimentalConditionId)
            solNum = odeSolutions[simulationInfo.experimentalConditionId[i]]
            solAlg = algebraicODESolutions[i]
            sqDiff += sum((Array(solNum) - Array(solAlg(solNum.t))).^2)
        end

        @test sqDiff ≤ 1e-6
    end
end


function testCostGradientOrHessianTestModel3(petabModel::PEtabModel, solver, tol)

    petabProblem1 = setUpPEtabODEProblem(petabModel, solver, solverAbsTol=tol, solverRelTol=tol, 
                                         sensealgZygote = ForwardDiffSensitivity(), 
                                         solverSSRelTol=1e-10, solverSSAbsTol=1e-12,
                                         odeSolverForwardEquations=CVODE_BDF(), sensealgForwardEquations = ForwardSensitivity(), 
                                         odeSolverAdjoint=solver, solverAdjointAbsTol=tol, solverAdjointRelTol=tol,
                                         sensealgAdjoint=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(true)))

    petabProblem2 = setUpPEtabODEProblem(petabModel, solver, solverAbsTol=tol, solverRelTol=tol, 
                                         odeSolverAdjoint=solver, solverAdjointAbsTol=tol, solverAdjointRelTol=tol,
                                         solverSSRelTol=1e-10, solverSSAbsTol=1e-12,
                                         sensealgAdjoint=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(true)), 
                                         sensealgAdjointSS=QuadratureAdjoint(autojacvec=ReverseDiffVJP(false)), 
                                         sensealgForwardEquations=:AutoDiffForward, odeSolverForwardEquations=solver)                                         

    _computeCostAlgebraic = (pArg) -> computeCostAlgebraic(pArg, petabModel, solver, tol)

    Random.seed!(123)
    createCube(petabProblem1, 5)
    cube = Matrix(CSV.read(petabProblem1.pathCube, DataFrame))

    for i in 1:5

        p = cube[i, :]
        
        referenceCost = _computeCostAlgebraic(p)
        referenceGradient = ForwardDiff.gradient(_computeCostAlgebraic, p)
        referenceHessian = ForwardDiff.hessian(_computeCostAlgebraic, p)

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

    return true 
end


petabModel = readPEtabModel("Test_model3", pwd() * "/tests/Test_model3/", forceBuildJlFile=true)

@testset verbose=true "Test model3" begin
    @testset "Test model3 : ODE solver" begin 
        testODESolverTestModel3(petabModel, Rodas4P(), 1e-12)
    end

    @testset "Test model3 : Cost gradient and hessian" begin 
        testCostGradientOrHessianTestModel3(petabModel, Rodas4P(autodiff=false), 1e-12)
    end

    @testset "Test model 3 : Gradient of residuals" begin
        checkGradientResiduals(petabModel, Rodas5(), 1e-9)
    end
end
