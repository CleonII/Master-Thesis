function runProcess(jobs, results) 
    
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
                                                    if simulationInfo.conditionIdSol[i] ∈ expIDs
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
            cost = evalF(paramVec)
            put!(results, tuple(:Done, cost))
        end
        
        if task == :Gradient
            evalGradF(gradVec, paramVec)
            put!(results, tuple(:Done, gradVec))
        end

        if task == :HessianApprox
            evalHessApprox(hessMat, paramVec)
            put!(results, tuple(:Done, hessMat))
        end

        if task == :Hessian
            evalHess(hessMat, paramVec)
            put!(results, tuple(:Done, hessMat))
        end
    end
end