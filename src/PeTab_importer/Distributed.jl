using Distributed


function setUpPEtabOptDistributed(peTabModel::PeTabModel,
                                  solver::SciMLAlgorithm, 
                                  tol::Float64,
                                  adjSolver::SciMLAlgorithm,
                                  adjSensealg::SciMLSensitivity.AbstractAdjointSensitivityAlgorithm,
                                  adjSensealgSS::SciMLSensitivity.AbstractAdjointSensitivityAlgorithm,
                                  adjTol::Float64,
                                  forwardSolver::SciMLAlgorithm, 
                                  sensealgForward::Union{Symbol, SciMLSensitivity.AbstractForwardSensitivityAlgorithm},
                                  parameterData::ParamData,
                                  measurementData::MeasurementData, 
                                  simulationInfo::SimulationInfo, 
                                  paramEstIndices::ParameterIndices, 
                                  pirorInfo::PriorInfo,
                                  odeProb::ODEProblem)

    println("Setting up cost, grad, and hessian to be computed on several processes using Distributed.jl")

    # Make functions, structs and packages aware for each process 
    loadPackages()
    loadFunctionsAndStructs()
    loadYmodSdU0(peTabModel)
    
    nProcs = nprocs()                                  
    if nProcs > length(simulationInfo.conditionIdSol)
        println("Warning - There are less experimental conditions than processes. Hence some processes will run empty")
        println("Number of processes = $nProcs, number of experimental conditions = ", length(simulationInfo.conditionIdSol))
    end
    idsEachProcess = collect(Iterators.partition(simulationInfo.conditionIdSol, Int(round(length(simulationInfo.conditionIdSol) /nProcs))))

    # Divide the experimental conditions between the number of processes and set up channels for 
    # communicating with each process
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
        sendPEtabStruct(pirorInfo, jobs[i], results[i], "Prior info", procs()[i])
    end

    # Send solver and solver tolerance 
    for i in 1:nProcs
        @async put!(jobs[i], tuple(solver, tol)) 
        status = take!(results[i])[1]
        if status != :Done
            println("Error : Could not send solver and tolerance problem to proces ", procs()[i])
        end
    end

    # Send adjoint solver, adjont tolerance and adjoint solver tolerance
    for i in 1:nProcs
        @async put!(jobs[i], tuple(adjSolver, adjTol, adjSensealg, adjSensealgSS)) 
        status = take!(results[i])[1]
        if status != :Done
            println("Error : Could not send adjSolver, adjTolerance and adjSensealg to proces ", procs()[i])
        end
    end

    # Send forward sensitivity equations solver and associated 
    for i in 1:nProcs
        @async put!(jobs[i], tuple(forwardSolver, sensealgForward)) 
        status = take!(results[i])[1]
        if status != :Done
            println("Error : Could not send adjSolver, adjTolerance and adjSensealg to proces ", procs()[i])
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

    evalF = (pVec) ->           begin
                                    costTot::Float64 = 0.0
                                    @inbounds for i in nProcs:-1:1
                                        @async put!(jobs[i], tuple(pVec, :Cost)) 
                                    end
                                    @inbounds for i in nProcs:-1:1
                                        status::Symbol, cost::Float64 = take!(results[i])
                                        if status != :Done
                                            println("Error : Could not send ODE problem to proces ", procs()[i])
                                        end
                                        costTot += cost
                                    end
                                    return costTot
                                end

    evalGradF = (grad, pVec) -> begin
                                    grad .= 0.0
                                    @inbounds for i in nProcs:-1:1
                                        @async put!(jobs[i], tuple(pVec, :Gradient)) 
                                    end
                                    @inbounds for i in nProcs:-1:1
                                        status::Symbol, gradPart::Vector{Float64} = take!(results[i])
                                        if status != :Done
                                            println("Error : Could not send ODE problem to proces ", procs()[i])
                                        end
                                        grad .+= gradPart
                                    end
                                end

    evalGradFAdjoint = (grad, pVec) -> begin
                                    grad .= 0.0
                                    @inbounds for i in nProcs:-1:1
                                        @async put!(jobs[i], tuple(pVec, :AdjGradient)) 
                                    end
                                    @inbounds for i in nProcs:-1:1
                                        status::Symbol, gradPart::Vector{Float64} = take!(results[i])
                                        if status != :Done
                                            println("Error : Could not send ODE problem to proces ", procs()[i])
                                        end
                                        grad .+= gradPart
                                    end
                                end       
                                
    evalGradFForwardEq = (grad, pVec) -> begin
                                    grad .= 0.0
                                    @inbounds for i in nProcs:-1:1
                                        @async put!(jobs[i], tuple(pVec, :ForwardSenseEqGradient)) 
                                    end
                                    @inbounds for i in nProcs:-1:1
                                        status::Symbol, gradPart::Vector{Float64} = take!(results[i])
                                        if status != :Done
                                            println("Error : Could not send ODE problem to proces ", procs()[i])
                                        end
                                        grad .+= gradPart
                                    end
                                end                                           

    evalHessA = (hess, pVec) -> begin
                                    hess .= 0.0
                                    @inbounds for i in nProcs:-1:1
                                        @async put!(jobs[i], tuple(pVec, :HessianApprox)) 
                                    end
                                    @inbounds for i in nProcs:-1:1
                                        status::Symbol, hessPart::Matrix{Float64} = take!(results[i])
                                        if status != :Done
                                            println("Error : Could not send ODE problem to proces ", procs()[i])
                                        end
                                        hess .+= hessPart
                                    end        
                                end

    evalHess = (hess, pVec) ->  begin
                                    hess .= 0.0
                                    @inbounds for i in nProcs:-1:1
                                        @async put!(jobs[i], tuple(pVec, :Hessian)) 
                                    end
                                    @inbounds for i in nProcs:-1:1
                                        status::Symbol, hessPart::Matrix{Float64} = take!(results[i])
                                        if status != :Done
                                            println("Error : Could not send ODE problem to proces ", procs()[i])
                                        end
                                        hess .+= hessPart
                                    end
                                end

    return tuple(evalF, evalGradF, evalGradFForwardEq, evalGradFAdjoint, evalHess, evalHessA)
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


function removeAllProcs()
    if length(procs()) > 1
        procsAvailble = procs()[2:end]
        rmprocs(procsAvailble)
    end
    return 
end


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
                                using SciMLSensitivity
                            end
                        end
                    end
    @eval @everywhere @LoadLib()
    println("Done")
end
 

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
                                include(joinpath(pwd(), "src", "PeTab_importer", "Distributed_run.jl"))
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
    pathDObsSdU0 = peTabModel.dirModel * peTabModel.modelName * "DObsSdU0.jl"
    pathCallbacks = peTabModel.dirModel * peTabModel.modelName * "Callbacks_time_piecewise.jl"
    @eval @everywhere include($pathObsSdU0)
    @eval @everywhere include($pathDObsSdU0)
    @eval @everywhere include($pathCallbacks)
    println("Done")
end

