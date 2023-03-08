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
using YAML
using FiniteDifferences


# Relevant PeTab structs for compuations 
include(joinpath(pwd(), "src", "PeTab_structs.jl"))

# PeTab importer to get cost, grad etc 
include(joinpath(pwd(), "src", "Create_PEtab_model.jl"))

# For converting to SBML 
include(joinpath(pwd(), "src", "SBML", "SBML_to_ModellingToolkit.jl"))

include(joinpath(pwd(), "tests", "Common.jl"))


petabModel = readPEtabModel(joinpath(@__DIR__, "Bachmann", "Bachmann_MSB2011.yaml"), forceBuildJuliaFiles=false)

solver, tol = Rodas5P(), 1e-9
petabProblem1 = setUpPEtabODEProblem(petabModel, solver, solverAbsTol=tol, solverRelTol=tol, 
                                     odeSolverAdjoint=solver, solverAdjointAbsTol=tol, solverAdjointRelTol=tol,
                                     sensealgAdjoint=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(true)))

# File with random parameter values 
paramVals = CSV.read(pwd() * "/tests/Bachmann/Params.csv", DataFrame)
paramMat = paramVals[!, Not([:Id, :SOCS3RNAEqc, :CISRNAEqc])]

# PyPesto hessian, gradient and cost values for the random parameter vectors 
costPython = (CSV.read(pwd() * "/tests/Bachmann/Cost.csv", DataFrame))[!, :Cost]
gradPythonMat = CSV.read(pwd() * "/tests/Bachmann/Grad.csv", DataFrame)
gradPythonMat = gradPythonMat[!, Not([:Id, :SOCS3RNAEqc, :CISRNAEqc])]
 
paramEstNames = string.(petabProblem1.θ_estNames)
# For correct indexing when comparing gradient (we have different ordering compared to PyPesto)
iUse = [findfirst(x -> x == paramEstNames[i], names(paramMat)) for i in eachindex(paramEstNames)]
nParam = ncol(paramMat)

i = 3
p = collect(paramMat[i, iUse])
referenceCost = costPython[i] # Cost from PyPesto
referenceGradient = collect(gradPythonMat[i, iUse]) # Cost from 

cost = petabProblem1.computeCost(p)
@printf("Cost python = %.3e and cost Julia = %.3e\n", referenceCost, cost)
@printf("abs(Difference cost) = %.3e\n", abs(referenceCost - cost))

# Gradient via Forward-mode automatic differentitation 
gradientForwardDiff = zeros(length(p))        
petabProblem1.computeGradientAutoDiff(gradientForwardDiff, p)
@printf("norm(gradientPyPesto - gradientForwardDiffJulia) = %.3e\n", norm(gradientForwardDiff - referenceGradient))

# Gradient via adjoint sensitivity analysis
gradientAdjoint = zeros(length(p))
petabProblem1.computeGradientAdjoint(gradientAdjoint, p)
@printf("norm(gradientPyPesto - gradientAdjointJulia) = %.3e\n", norm(gradientAdjoint - referenceGradient))
        
# Gradient from FiniteDifferences 
gradientFinite = FiniteDifferences.grad(central_fdm(5, 1), petabProblem1.computeCost, p)[1]
@printf("norm(gradientPyPesto - gradientFiniteDifferences) = %.3e\n", norm(gradientFinite - referenceGradient))

println("Gradient Adjoint = ", gradientAdjoint)
println("Gradient PyPesto = ", referenceGradient)