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
using Test
using YAML


# Relevant PeTab structs for compuations 
include(joinpath(pwd(), "src", "PeTab_structs.jl"))

# PeTab importer to get cost, grad etc 
include(joinpath(pwd(), "src", "Create_PEtab_model.jl"))

# HyperCube sampling 
include(joinpath(pwd(), "src", "Optimizers", "Lathin_hypercube.jl"))

# For converting to SBML 
include(joinpath(pwd(), "src", "SBML", "SBML_to_ModellingToolkit.jl"))

include(joinpath(pwd(), "tests", "Common.jl"))


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


function computeCostAnalytic(paramVec)

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


function testODESolverTestModel1(petabModel::PEtabModel, solver, tol)

    # Set values to PeTab file values 
    experimentalConditionsFile, measurementDataFile, parameterDataFile, observablesDataFile = readPEtabFiles(petabModel)
    measurementData = processMeasurements(measurementDataFile, observablesDataFile) 
    paramData = processParameters(parameterDataFile) 
    setParamToFileValues!(petabModel.parameterMap, petabModel.stateMap, paramData)
    θ_indices = computeIndicesθ(paramData, measurementData, petabModel.odeSystem, experimentalConditionsFile)
    
    # Extract experimental conditions for simulations 
    simulationInfo = processSimulationInfo(petabModel, measurementData, paramData)
    
    # Parameter values where to test accuracy. Each column is a α, β, γ and δ
    u0 = [8.0, 4.0]
    parametersTest = reshape([2.0, 3.0, 3.0, 5.0, 
                              1.0, 2.0, 3.0, 3.0,
                              1.0, 0.4, 0.3, 0.5,
                              4.0, 3.0, 2.0, 1.0,
                              0.01, 0.02, 0.03, 0.04], (4, 5))

    for i in 1:5

        α, β, γ, δ = parametersTest[:, i]
        
        # Set parameter values for the ODE-model 
        petabModel.parameterMap[2] = Pair(petabModel.parameterMap[2].first, α)
        petabModel.parameterMap[4] = Pair(petabModel.parameterMap[4].first, γ)
        petabModel.parameterMap[5] = Pair(petabModel.parameterMap[5].first, δ)
        petabModel.parameterMap[6] = Pair(petabModel.parameterMap[6].first, β)
        prob = ODEProblem(petabModel.odeSystem, petabModel.stateMap, (0.0, 5e3), petabModel.parameterMap, jac=true)
        prob = remake(prob, p = convert.(Float64, prob.p), u0 = convert.(Float64, prob.u0))
        
        # Solve the ODE-model using the PEtab solver 
        θ_est = getFileODEvalues(petabModel)
        changeExperimentalCondition! = (pVec, u0Vec, expID) -> _changeExperimentalCondition!(pVec, u0Vec, expID, θ_est, petabModel, θ_indices)
        odeSolutions, success = solveODEAllExperimentalConditions(prob, changeExperimentalCondition!, simulationInfo, solver, tol, tol, petabModel.computeTStops)
        odeSolution = odeSolutions[simulationInfo.experimentalConditionId[1]]
        
        # Compare against analytical solution 
        sqDiff = 0.0
        for t in odeSolution.t
            solutionAnalytic = solveOde2x2Lin(t, u0, α*0.5, β, γ, δ)
            sqDiff += sum((odeSolution(t)[1:2] - solutionAnalytic).^2)
        end

        @test sqDiff ≤ 1e-8
    end
end


function testCostGradientOrHessianTestModel1(petabModel::PEtabModel, solver, tol)

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
        referenceCost = computeCostAnalytic(p)
        referenceGradient = ForwardDiff.gradient(computeCostAnalytic, p)
        referenceHessian = ForwardDiff.hessian(computeCostAnalytic, p)

        # Test both the standard and Zygote approach to compute the cost 
        cost = _testCostGradientOrHessian(petabProblem1, p, cost=true)
        @test cost ≈ referenceCost atol=1e-4
        costZygote = _testCostGradientOrHessian(petabProblem1, p, costZygote=true)
        @test costZygote ≈ referenceCost atol=1e-4

        # Test all gradient combinations. Note we test sensitivity equations with and without autodiff 
        gradientAutoDiff = _testCostGradientOrHessian(petabProblem1, p, gradientAutoDiff=true)
        @test norm(gradientAutoDiff - referenceGradient) ≤ 1e-4
        gradientZygote = _testCostGradientOrHessian(petabProblem1, p, gradientZygote=true)
        @test norm(gradientZygote - referenceGradient) ≤ 1e-4
        gradientAdjoint = _testCostGradientOrHessian(petabProblem1, p, gradientAdjoint=true)
        @test norm(normalize(gradientAdjoint) - normalize((referenceGradient))) ≤ 1e-4
        gradientForwardEquations1 = _testCostGradientOrHessian(petabProblem1, p, gradientForwardEquations=true)
        @test norm(gradientForwardEquations1 - referenceGradient) ≤ 1e-4
        gradientForwardEquations2 = _testCostGradientOrHessian(petabProblem2, p, gradientForwardEquations=true)
        @test norm(gradientForwardEquations2 - referenceGradient) ≤ 1e-4

        # Testing "exact" hessian via autodiff 
        hessian = _testCostGradientOrHessian(petabProblem1, p, hessian=true)
        @test norm(hessian - referenceHessian) ≤ 1e-3
    end
end


# Test different adjoint approaches 
function testZYgoteAdjointsTestModel1(petabModel::PEtabModel, solver, tol)

    petabProblem1 = setUpPEtabODEProblem(petabModel, solver, solverAbsTol=tol, solverRelTol=tol, sensealgZygote = ForwardDiffSensitivity())
    petabProblem2 = setUpPEtabODEProblem(petabModel, solver, solverAbsTol=tol, solverRelTol=tol, sensealgZygote = QuadratureAdjoint())
    petabProblem3 = setUpPEtabODEProblem(petabModel, solver, solverAbsTol=tol, solverRelTol=tol, sensealgZygote = InterpolatingAdjoint(autojacvec=ReverseDiffVJP()))
    petabProblem4 = setUpPEtabODEProblem(petabModel, solver, solverAbsTol=tol, solverRelTol=tol, sensealgZygote = BacksolveAdjoint())
    petabProblem5 = setUpPEtabODEProblem(petabModel, solver, solverAbsTol=tol, solverRelTol=tol, sensealgZygote = ReverseDiffAdjoint())

    Random.seed!(123)
    createCube(petabProblem1, 5)
    cube = Matrix(CSV.read(petabProblem1.pathCube, DataFrame))

    p = cube[1, :]
    referenceGradient = ForwardDiff.gradient(computeCostAnalytic, p)

    gradient1 = _testCostGradientOrHessian(petabProblem1, p, gradientZygote=true)
    @test norm(gradient1 - referenceGradient) ≤ 1e-4
    gradient2 = _testCostGradientOrHessian(petabProblem2, p, gradientZygote=true)
    @test norm(gradient2 - referenceGradient) ≤ 1e-3
    gradient3 = _testCostGradientOrHessian(petabProblem3, p, gradientZygote=true)
    @test norm(gradient3 - referenceGradient) ≤ 1e-3
    gradient4 = _testCostGradientOrHessian(petabProblem4, p, gradientZygote=true)
    @test norm(gradient4 - referenceGradient) ≤ 1e-3
end


petabModel = readPEtabModel(joinpath(@__DIR__, "Test_model1", "Test_model1.yaml"), forceBuildJuliaFiles=true)

@testset verbose = true "Test model 1" begin
    @testset "Test model1 : ODE solver" begin 
        testODESolverTestModel1(petabModel, Vern9(), 1e-9)
    end

    @testset "Test model1 : Cost gradient and hessian" begin 
        testCostGradientOrHessianTestModel1(petabModel, Vern9(), 1e-12)
    end

    @testset "Test model 1 : Zygote adjoints" begin
        testZYgoteAdjointsTestModel1(petabModel, Vern9(), 1e-12)
    end

    @testset "Test model 1 : Gradient of residuals" begin
        checkGradientResiduals(petabModel, Rodas5(), 1e-9)
    end
end
