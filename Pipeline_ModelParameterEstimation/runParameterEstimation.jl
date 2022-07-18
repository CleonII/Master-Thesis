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

allModelFunctionVector = includeAllModels(getModelFiles(joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels")), 
        joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels"))



modelName = "model_Boehm_JProteomeRes2014"
solver = Rodas4P()
tol = 1e-9

evalF, evalGradF, evalHessianApproxF, paramVecEstTmp, lb, ub, stateMap = setUpCostFunc(modelName, solver, tol)
evalH = (hessianMat, paramVec) -> begin hessianMat .= Symmetric(ForwardDiff.hessian(evalF, paramVecEstTmp)) end

ipoptProb, iterVec = createIpoptProbNew(evalF, evalGradF, evalH, lb, ub)

Random.seed!(123)
p0 = [rand(Uniform(lb[i], ub[i])) for i in eachindex(lb)]
Ipopt.AddIpoptStrOption(ipoptProb, "hessian_approximation", "exact")
Ipopt.AddIpoptIntOption(ipoptProb, "print_level", 5)
Ipopt.AddIpoptIntOption(ipoptProb, "max_iter", 1000)
ipoptProb.x .= deepcopy(p0)
sol_opt = Ipopt.IpoptSolve(ipoptProb)
