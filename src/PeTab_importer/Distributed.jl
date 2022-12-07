using Distributed
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
using SciMLSensitivity
using BenchmarkTools
using Zygote


include(joinpath(pwd(), "src", "PeTab_structs.jl"))
include(joinpath(pwd(), "src", "Solve_ODE_model", "Solve_ode_model.jl"))
include(joinpath(pwd(), "src", "PeTab_importer", "Create_cost_grad_hessian.jl"))
include(joinpath(pwd(), "src", "PeTab_importer", "Common.jl"))
include(joinpath(pwd(), "src", "PeTab_importer", "Map_parameters.jl"))
include(joinpath(pwd(), "src", "PeTab_importer", "Process_PeTab_files.jl"))
include(joinpath(pwd(), "src", "Common.jl"))

# Helper function to remove all processes 
function removeAllProcs()
    if length(procs()) > 1
        procsAvailble = procs()[2:end]
        rmprocs(procsAvailble)
    end
    return 
end


# When using MultiProcessing we need to ensure that appropriate pacakges are loaded 
# for each process 
function loadPackages()
    println("Loading required packages for each process")
    @eval @everywhere begin 
                        macro LoadLib()
                            quote
                                using DifferentialEquations
                                using ModelingToolkit
                                using DataFrames
                                using LinearAlgebra
                                using ForwardDiff
                                using ForwardDiff
                                using ReverseDiff
                                using Zygote
                                using Printf
                            end
                        end
                    end
    @eval @everywhere @LoadLib()
    println("Done")
end
 

# Load the relevant structus and packages in order to evaluate the cost 
function loadFunctionsAndStructs()
    println("Loading required functions and structs for each process")
    @eval @everywhere begin 
                        macro LoadFuncStruct()
                            quote
                                include(joinpath(pwd(), "src", "PeTab_structs.jl"))
                                include(joinpath(pwd(), "src", "Solve_ODE_model", "Solve_ode_model.jl"))
                                include(joinpath(pwd(), "src", "PeTab_importer", "Create_cost_grad_hessian.jl"))
                                include(joinpath(pwd(), "src", "PeTab_importer", "Common.jl"))
                                include(joinpath(pwd(), "src", "PeTab_importer", "Map_parameters.jl"))
                                include(joinpath(pwd(), "src", "PeTab_importer", "Process_PeTab_files.jl"))
                                include(joinpath(pwd(), "src", "Common.jl"))
                            end
                        end
                    end
    @eval @everywhere @LoadFuncStruct()
    println("Done")
end
 

function loadYmodSdU0(peTabModel::PeTabModel)
    println("Loading yMod, SD and u0 functions")
    pathObsSdU0 = peTabModel.dirModel * peTabModel.modelName * "ObsSdU0.jl"
    @eval @everywhere include($pathObsSdU0)
    println("Done")
end



dirModel = pwd() * "/Intermediate/PeTab_models/model_Fiedler_BMC2016/"
peTabModel = setUpPeTabModel("model_Fiedler_BMC2016", dirModel, verbose=false)
solver, tol = Rodas4P(), 1e-9

# Process PeTab files into type-stable Julia structs 
experimentalConditionsFile, measurementDataFile, parameterDataFile, observablesDataFile = readDataFiles(peTabModel.dirModel, readObs=true)
parameterData = processParameterData(parameterDataFile)
measurementData = processMeasurementData(measurementDataFile, observablesDataFile) 
simulationInfo = getSimulationInfo(measurementDataFile, measurementData, absTolSS=1e-8, relTolSS=1e-8)

# Indices for mapping parameter-estimation vector to dynamic, observable and sd parameters correctly when calculating cost
paramEstIndices = getIndicesParam(parameterData, measurementData, peTabModel.odeSystem, experimentalConditionsFile)
    
# Set up potential prior for the parameters to estimate 
priorInfo = getPriorInfo(paramEstIndices::ParameterIndices, parameterDataFile::DataFrame)::PriorInfo

# Set model parameter values to those in the PeTab parameter data ensuring correct value of constant parameters 
setParamToFileValues!(peTabModel.paramMap, peTabModel.stateMap, parameterData)

# The time-span 5e3 is overwritten when performing actual forward simulations 
odeProb = ODEProblem(peTabModel.odeSystem, peTabModel.stateMap, (0.0, 5e3), peTabModel.paramMap, jac=true, sparse=false)
odeProb = remake(odeProb, p = convert.(Float64, odeProb.p), u0 = convert.(Float64, odeProb.u0))

removeAllProcs()
addprocs(1, exeflags="--project=.")

nProcs = nprocs()

# Make functions, structs and packages aware for each process 
loadPackages()
loadFunctionsAndStructs()
loadYmodSdU0(peTabModel)

@everywhere function runProcess(jobs, results) 
    
    # Import actual ODE model
    odeProb::ODEProblem = take!(jobs)[1]
    put!(results, tuple(:Done))

    # Import structs needed to compute the cost, gradient, and hessian
    peTabModel::PeTabModel = take!(jobs)[1]
    put!(results, tuple(:Done))
    parameterData::ParamData = take!(jobs)[1]
    put!(results, tuple(:Done))
    measurementData::MeasurementData = take!(jobs)[1]
    put!(results, tuple(:Done))
    simulationInfo::SimulationInfo = take!(jobs)[1]
    put!(results, tuple(:Done))
    paramEstIndices::ParameterIndices = take!(jobs)[1]
    put!(results, tuple(:Done))
    priorInfo::PriorInfo = take!(jobs)[1]
    put!(results, tuple(:Done))

    println("Done loading structs for ", myid())

    solver, tol::Float64 = take!(jobs)
    put!(results, tuple(:Done))

    expIDs::Array{String, 1} = take!(jobs)[1]

    # Set up cost, gradient, and hessian functions 
    changeToExperimentalCondUse! = (pVec, u0Vec, expID, dynParamEst) -> changeExperimentalCondEst!(pVec, u0Vec, expID, dynParamEst, peTabModel, paramEstIndices)
    changeModelParamUse! = (pVec, u0Vec, paramEst) -> changeModelParam!(pVec, u0Vec, paramEst, paramEstIndices, peTabModel)

    # Set up function which solves the ODE model for all conditions and stores result 
    solveOdeModelAllCondUse! = (solArrayArg, odeProbArg, dynParamEst, expIDSolveArg) -> solveOdeModelAllExperimentalCond!(solArrayArg, odeProbArg, dynParamEst, changeToExperimentalCondUse!, simulationInfo, solver, tol, tol, onlySaveAtTobs=true, expIDSolve=expIDSolveArg)
    evalF = (paramVecEst) -> calcCost(paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, priorInfo, expIDSolve=expIDs)
    evalGradF = (grad, paramVecEst) -> calcGradCost!(grad, paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, priorInfo, expIDSolve=expIDs)
    evalHessApprox = (hessianMat, paramVecEst) -> calcHessianApprox!(hessianMat, paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, priorInfo, expIDSolve=expIDs)
    
    _evalHess = (paramVecEst) -> calcCost(paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, priorInfo, calcHessian=true, expIDSolve=expIDs)
    evalHess = (hessianMat, paramVec) ->    begin 
                                                computeH = true
                                                @inbounds for i in eachindex(simulationInfo.solArray)
                                                    if (expIDs[1] == "all" || simulationInfo.conditionIdSol[i] ∈ expIDs) 
                                                        if simulationInfo.solArray[i].retcode != :Success
                                                            computeH = false
                                                        end
                                                    end
                                                end
                                                if computeH
                                                    try 
                                                        hessianMat .= Symmetric(ForwardDiff.hessian(_evalHess, paramVec))
                                                    catch
                                                        hessianMat .= 0.0
                                                    end
                                                else
                                                    hessianMat .= 0.0
                                                end
                                            end

    namesParamEst = paramEstIndices.namesParamEst                                         
    gradVec = zeros(length(namesParamEst))
    hessMat = zeros(length(namesParamEst), length(namesParamEst))                                         
    println("Done setting up cost, grad, and hessian for process ", myid())    
    put!(results, tuple(:Done))  

    while true 
        paramVec::Vector{Float64}, task::Symbol = take!(jobs)
        if task == :Cost
            println("Starting to compute cost")
            cost = evalF(paramVec)
            println("Done computing cost with cost = $cost")
            put!(results, tuple(:Done, cost))
        end
        
        if task == :Gradient
            println("Starting to compute gradient")
            evalGradF(gradVec, paramVec)
            println("Done computing gradient")
            put!(results, tuple(:Done, gradVec))
        end

        if task == :HessianApprox
            println("Starting to compute hessian approximation")
            evalHessApprox(hessMat, paramVec)
            println("Done computing hessian approximation")
            put!(results, tuple(:Done, hessMat))
        end

        if task == :Hessian
            println("Starting to compute hessian")
            evalHess(hessMat, paramVec)
            println("Done computing hessian")
            put!(results, tuple(:Done, hessMat))
        end
    end
end


function sendPEtabStruct(structSend::Union{PeTabModel, ParamData, MeasurementData, SimulationInfo, ParameterIndices, PriorInfo}, 
                         job::RemoteChannel, 
                         result::RemoteChannel, 
                         strSend::String, 
                         pID::Integer)

    @async put!(job, tuple(deepcopy(structSend)))   
    status = take!(result)[1]                     
    if status != :Done
        println("Error : Could not send $strSend to process ", pID)
    end
end

# Divide the experimental conditions between the number of processes and set up channels for 
# communicating with each process
idsEachProcess = collect(Iterators.partition(simulationInfo.conditionIdSol, Int(round(length(simulationInfo.conditionIdSol) /nProcs))))
jobs = [RemoteChannel(()->Channel{Tuple}(1)) for i in 1:nProcs]
results = [RemoteChannel(()->Channel{Tuple}(1)) for i in 1:nProcs]

# Send ODE-problem, and simultaneously launch the processes 
for i in 1:nProcs
    @async put!(jobs[i], tuple(deepcopy(odeProb))) 
    remote_do(runProcess, procs()[i], jobs[i], results[i])
    status = take!(results[i])[1]
    if status != :Done
        println("Error : Could not send ODE problem to proces ", procs()[i])
    end
end

# Send required PEtab structs to processes 
for i in 1:nProcs
    sendPEtabStruct(peTabModel, jobs[i], results[i], "PEtab model", procs()[i])
    sendPEtabStruct(parameterData, jobs[i], results[i], "Parameter data", procs()[i])
    sendPEtabStruct(measurementData, jobs[i], results[i], "Measurement data", procs()[i])
    sendPEtabStruct(simulationInfo, jobs[i], results[i], "Simulation info", procs()[i])
    sendPEtabStruct(paramEstIndices, jobs[i], results[i], "Parameter indices", procs()[i])
    sendPEtabStruct(priorInfo, jobs[i], results[i], "Prior info", procs()[i])
end

# Send solver and solver tolerance 
for i in 1:nProcs
    @async put!(jobs[i], tuple(solver, tol)) 
    status = take!(results[i])[1]
    if status != :Done
        println("Error : Could not send ExpIds problem to proces ", procs()[i])
    end
end


# Send experimental ID:s process each process works with 
for i in 1:nProcs
    @async put!(jobs[i], tuple(idsEachProcess[i])) 
    status = take!(results[i])[1]
    if status != :Done
        println("Error : Could not send ExpIds problem to proces ", procs()[i])
    end
end


# Compute the cost in paralell the cost 
namesParamEst = paramEstIndices.namesParamEst                                         
paramVecNominal = [parameterData.paramVal[findfirst(x -> x == namesParamEst[i], parameterData.parameterID)] for i in eachindex(namesParamEst)]
pVec = transformParamVec(paramVecNominal, namesParamEst, parameterData, revTransform=true)
@elapsed begin
costTot = 0.0
for i in 1:nProcs
    @async put!(jobs[i], tuple(pVec, :Cost)) 
end
for i in 1:nProcs
    status, cost = take!(results[i])
    if status != :Done
        println("Error : Could not send ODE problem to proces ", procs()[i])
    end
    costTot += cost
end
end

grad = zeros(length(pVec))
for i in 1:nProcs
    @async put!(jobs[i], tuple(pVec, :Gradient)) 
end
for i in 1:nProcs
    status, gradPart = take!(results[i])
    if status != :Done
        println("Error : Could not send ODE problem to proces ", procs()[i])
    end
    grad .+= gradPart
end

@elapsed begin
hess = zeros(length(pVec), length(pVec))
for i in 1:nProcs
    @async put!(jobs[i], tuple(pVec, :HessianApprox)) 
end
for i in 1:nProcs
    status, hessPart = take!(results[i])
    if status != :Done
        println("Error : Could not send ODE problem to proces ", procs()[i])
    end
    hess .+= hessPart
end
end

hess = zeros(length(pVec), length(pVec))
for i in 1:nProcs
    @async put!(jobs[i], tuple(pVec, :Hessian)) 
end
for i in 1:nProcs
    status, hessPart = take!(results[i])
    if status != :Done
        println("Error : Could not send ODE problem to proces ", procs()[i])
    end
    hess .+= hessPart
end



changeToExperimentalCondUse! = (pVec, u0Vec, expID, dynParamEst) -> changeExperimentalCondEst!(pVec, u0Vec, expID, dynParamEst, peTabModel, paramEstIndices)
changeModelParamUse! = (pVec, u0Vec, paramEst) -> changeModelParam!(pVec, u0Vec, paramEst, paramEstIndices, peTabModel)

expIDs::Array{String, 1} = idsEachProcess[1]


# Set up function which solves the ODE model for all conditions and stores result 
solveOdeModelAllCondUse! = (solArrayArg, odeProbArg, dynParamEst, expIDSolveArg) -> solveOdeModelAllExperimentalCond!(solArrayArg, odeProbArg, dynParamEst, changeToExperimentalCondUse!, simulationInfo, solver, tol, tol, onlySaveAtTobs=true, expIDSolve=expIDSolveArg)
evalF = (paramVecEst) -> calcCost(paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, priorInfo, expIDSolve=expIDs)
evalGradF = (grad, paramVecEst) -> calcGradCost!(grad, paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, priorInfo, expIDSolve=expIDs)
evalHessApprox = (hessianMat, paramVecEst) -> calcHessianApprox!(hessianMat, paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, priorInfo, expIDSolve=expIDs)
    
_evalHess = (paramVecEst) -> calcCost(paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, priorInfo, calcHessian=true, expIDSolve=expIDs)
evalHess = (hessianMat, paramVec) -> begin 
                                        computeH = true
                                        @inbounds for i in eachindex(simulationInfo.solArray)
                                            if (expIDs[1] == "all" || simulationInfo.conditionIdSol[i] ∈ expIDs) 
                                                if simulationInfo.solArray[i].retcode != :Success
                                                    computeH = false
                                                end
                                            end
                                        end
                                        if computeH
                                            try 
                                                hessianMat .= Symmetric(ForwardDiff.hessian(_evalHess, paramVec))
                                            catch
                                                hessianMat .= 0.0
                                            end
                                        else
                                            hessianMat .= 0.0
                                        end
                                        end
evalF(pVec)
evalHess(hess, pVec)                                        