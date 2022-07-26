using ModelingToolkit, DifferentialEquations, BenchmarkTools, DataFrames, CSV, LSODA, Sundials, ODEInterface, ODEInterfaceDiffEq
using ForwardDiff
using Plots
using StatsBase
using Random
using LinearAlgebra
using Calculus
using Ipopt
using Distributions
#plotly()

# TODO : Fix observebleTransformation
# TODO : Model-solver why Beer fails 

include(joinpath(pwd(), "Pipeline_ModelSolver", "BigFloatODEProblem.jl"))
include(joinpath(pwd(), "Additional_functions", "additional_tools.jl"))
include(joinpath(pwd(), "Additional_functions", "Solver_info.jl"))
include(joinpath(pwd(), "Additional_functions", "benchmarkSolvers.jl"))

# Ipopt wrapper 
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "CommonParameterEstimationMethods", "CreateIpoptProb.jl"))

# PeTab importer 
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "PeTabImporter.jl"))

# Observation function, needs to be generalised 
include("/home/sebpe/Dropbox/PhD/Projects/Master-Thesis/Pipeline_ModelParameterEstimation/Data/model_Boehm_JProteomeRes2014/Boehm_JProteomeRes2014Obs.jl")
include("/home/sebpe/Dropbox/PhD/Projects/Master-Thesis/Pipeline_ModelParameterEstimation/Data/model_Bachmann_MSB2011/Bachmann_MSB2011Obs.jl")

allModelFunctionVector = includeAllModels(getModelFiles(joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels")), 
        joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels"))

modelName = "model_Boehm_JProteomeRes2014"
solver = QNDF()
tol = 1e-9

evalObs = Boehm_JProteomeRes2014
evalU0 = Boehm_JProteomeRes2014_t0!
evalSd = Boehm_JProteomeRes2014_sd!

evalF, evalGradF, evalHessianApproxF, paramVecEstTmp, lowerBounds, upperBounds, idParam = setUpCostFunc(modelName, evalObs, evalU0, evalSd, solver, tol)


cost = evalF(paramVecEstTmp)

gradUseF = (vec, grad) -> grad .= ForwardDiff.gradient(evalF, vec)
Random.seed!(123)
p0 = [rand(Uniform(lowerBounds[i], upperBounds[i])) for i in eachindex(lowerBounds)]

println("Starting compute gradients")
grad1, grad2 = zeros(length(paramVecEstTmp)), zeros(length(paramVecEstTmp))
evalGradF(p0, grad1)
gradUseF(p0, grad2)

println(exp10.(p0))

p0 .= paramVecEstTmp .+ 0.1

evalH = (hessianMat, paramVec) -> begin hessianMat .= Symmetric(ForwardDiff.hessian(evalF, paramVecEstTmp)) end
hess = zeros(9, 9)
hessApprox = zeros(9, 9)
evalH(hess, p0)
evalHessianApproxF(hessApprox, p0)

#=
ipoptProb, iterVec = createIpoptProbNew(evalF, evalGradF, evalH, lowerBounds, upperBounds)


Ipopt.AddIpoptStrOption(ipoptProb, "hessian_approximation", "exact")
Ipopt.AddIpoptIntOption(ipoptProb, "print_level", 5)
Ipopt.AddIpoptIntOption(ipoptProb, "max_iter", 500)
ipoptProb.x .= deepcopy(p0)
cost = ipoptProb.eval_f(p0)
sol_opt = Ipopt.IpoptSolve(ipoptProb)

=#

# For Bachman model 
#=
modelName = "model_Bachmann_MSB2011"
solver = QNDF()
tol = 1e-9

evalObs = Bachmann_MSB2011
evalU0 = Bachmann_MSB2011_t0!
evalSd = Bachmann_MSB2011_sd!

evalF, evalGradF, evalHessianApproxF, paramVecEstTmp, lowerBounds, upperBounds, idParam = setUpCostFunc(modelName, evalObs, evalU0, evalSd, solver, tol)

println("Starting compute gradients")
gradUse = zeros(length(paramVecEstTmp))
evalGradF(paramVecEstTmp, gradUse)

hessianUse = zeros(length(paramVecEstTmp), length(paramVecEstTmp))
evalHessianApproxF(hessianUse, paramVecEstTmp)
=#