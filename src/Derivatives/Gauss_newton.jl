function computeJacobianResidualsDynamicθ!(jacobian::Union{Matrix{Float64}, SubArray},
                                           θ_dynamic::Vector{Float64},
                                           θ_sd::Vector{Float64},
                                           θ_observable::Vector{Float64},
                                           θ_nonDynamic::Vector{Float64},
                                           S::Matrix{Float64},
                                           peTabModel::PeTabModel,
                                           odeProblem::ODEProblem,
                                           simulationInfo::SimulationInfo,
                                           θ_indices::ParameterIndices,
                                           measurementData ::MeasurementData, 
                                           parameterInfo::ParameterInfo, 
                                           changeODEProblemParameters!::Function,
                                           solveOdeModelAllConditions!::Function;
                                           expIDSolve::Array{String, 1} = ["all"])

    θ_dynamicT = transformθ(θ_dynamic, θ_indices.θ_dynamicNames, parameterInfo)
    θ_sdT = transformθ(θ_sd, θ_indices.θ_sdNames, parameterInfo)
    θ_observableT = transformθ(θ_observable, θ_indices.θ_observableNames, parameterInfo)
    θ_nonDynamicT = transformθ(θ_nonDynamic, θ_indices.θ_nonDynamicNames, parameterInfo)

    # Solve the expanded ODE system for the sensitivites
    success = solveForSensitivites(S, odeProblem, simulationInfo, peTabModel, :AutoDiff, θ_dynamicT, 
                                   solveOdeModelAllConditions!, changeODEProblemParameters!, expIDSolve)
    if success != true
        println("Failed to solve sensitivity equations")
        jacobian .= 1e8
        return
    end

    jacobian .= 0.0
    # Compute the gradient by looping through all experimental conditions.
    for conditionID in keys(measurementData.iPerConditionId)

        if expIDSolve[1] != "all" && conditionID ∉ expIDSolve
            continue
        end

        whichForwardODESolution = findfirst(x -> x == conditionID, simulationInfo.conditionIdSol)
        sol = simulationInfo.solArrayGrad[whichForwardODESolution]
        positionInSolArray = simulationInfo.posInSolArray[conditionID]

        # In case the model is simulated first to a steady state we need to keep track of the post-equlibrium experimental 
        # condition Id to identify parameters specific to an experimental condition.
        postEqulibriumId = simulationInfo.simulateSS == true ? simulationInfo.postEqIdSol[whichForwardODESolution] : conditionID

        # If we have a callback it needs to be properly handled 
        computeJacobianResidualsExpCond!(jacobian, sol, S, θ_dynamicT, θ_sdT, θ_observableT, θ_nonDynamicT,
                                         conditionID, postEqulibriumId, positionInSolArray, peTabModel, θ_indices, 
                                         measurementData, parameterInfo)
    end
end


function computeJacobianResidualsExpCond!(jacobian::AbstractMatrix,
                                          sol::ODESolution,
                                          S::Matrix{Float64},
                                          θ_dynamic::Vector{Float64},
                                          θ_sd::Vector{Float64}, 
                                          θ_observable::Vector{Float64}, 
                                          θ_nonDynamic::Vector{Float64},
                                          conditionID::String,
                                          postEqulibriumId::String,
                                          positionInSolArray::UnitRange{Int64},
                                          peTabModel::PeTabModel,
                                          θ_indices::ParameterIndices,
                                          measurementData::MeasurementData, 
                                          parameterInfo::ParameterInfo)

    iPerTimePoint = measurementData.iGroupedTObs[conditionID]
    # Pre allcoate vectors needed for computations 
    ∂h∂u, ∂σ∂u, ∂h∂p, ∂σ∂p = allocateObservableFunctionDerivatives(sol, peTabModel) 
    
    # To compute 
    compute∂G∂u = (out, u, p, t, i, it) -> begin compute∂G∂_(out, u, p, t, i, it, 
                                                         measurementData, parameterInfo, 
                                                         θ_indices, peTabModel, 
                                                         θ_dynamic, θ_sd, θ_observable, θ_nonDynamic, 
                                                         ∂h∂u, ∂σ∂u, compute∂G∂U=true, 
                                                         computeResiduals=true)
                                            end
    compute∂G∂p = (out, u, p, t, i, it) -> begin compute∂G∂_(out, u, p, t, i, it, 
                                                         measurementData, parameterInfo, 
                                                         θ_indices, peTabModel, 
                                                         θ_dynamic, θ_sd, θ_observable, θ_nonDynamic, 
                                                         ∂h∂p, ∂σ∂p, compute∂G∂U=false, 
                                                         computeResiduals=true)
                                            end
    
    # Extract relevant parameters for the experimental conditions 
    mapConditionId = θ_indices.mapsConiditionId[Symbol(postEqulibriumId)]                                 
    iθ_experimentalCondition = vcat(θ_indices.mapODEProblem.iθDynamic, mapConditionId.iθDynamic)                                                                     

    # Loop through solution and extract sensitivites                                                 
    tSaveAt = measurementData.tVecSave[conditionID]
    nModelStates = length(peTabModel.stateNames)
    p = dualToFloat.(sol.prob.p)
    ∂G∂p = zeros(Float64, length(sol.prob.p))    
    ∂G∂u = zeros(Float64, nModelStates)
    for i in eachindex(tSaveAt)     
        u = dualToFloat.(sol[:, i])
        t = tSaveAt[i]
        iStart, iEnd = (positionInSolArray[i]-1)*nModelStates+1, (positionInSolArray[i]-1)*nModelStates + nModelStates
        _S = @view S[iStart:iEnd, iθ_experimentalCondition]
        for indexMeasurementData in iPerTimePoint[i]
            compute∂G∂u(∂G∂u, u, p, t, 1, [[indexMeasurementData]])
            compute∂G∂p(∂G∂p, u, p, t, 1, [[indexMeasurementData]])
            _gradient = _S'*∂G∂u 
            
            # Thus far have have computed dY/dθ for the residuals, but for parameters on the log-scale we want dY/dθ_log. 
            # We can adjust via; dY/dθ_log = log(10) * θ * dY/dθ
            adjustGradientTransformedParameters!((@views jacobian[iθ_experimentalCondition, indexMeasurementData]), _gradient, 
                                                 ∂G∂p, θ_dynamic, θ_indices, parameterInfo, postEqulibriumId, autoDiffSensitivites=true)
                             
        end
    end
end


# To compute the gradient for non-dynamic parameters 
function computeResidualsNotSolveODE(θ_dynamic::Vector{Float64},
                                     θ_sd::AbstractVector, 
                                     θ_observable::AbstractVector, 
                                     θ_nonDynamic::AbstractVector,
                                     peTabModel::PeTabModel,
                                     simulationInfo::SimulationInfo,
                                     θ_indices::ParameterIndices,
                                     measurementData::MeasurementData,
                                     parameterInfo::ParameterInfo; 
                                     expIDSolve::Array{String, 1} = ["all"])::AbstractVector 

    residuals = zeros(eltype(θ_sd), length(measurementData.tObs))

    # To be able to use ReverseDiff sdParamEstUse and obsParamEstUse cannot be overwritten. 
    # Hence new vectors have to be created.
    θ_dynamicT = transformθ(θ_dynamic, θ_indices.θ_dynamicNames, parameterInfo)
    θ_sdT = transformθ(θ_sd, θ_indices.θ_sdNames, parameterInfo)
    θ_observableT = transformθ(θ_observable, θ_indices.θ_observableNames, parameterInfo)
    θ_nonDynamicT = transformθ(θ_nonDynamic, θ_indices.θ_nonDynamicNames, parameterInfo)
    
    # Compute residuals per experimental conditions
    odeSolutionArray = simulationInfo.solArrayGrad
    for conditionID in keys(measurementData.iPerConditionId)

        if expIDSolve[1] != "all" && conditionID ∉ expIDSolve
            continue
        end

        odeSolution = odeSolutionArray[findfirst(x -> x == conditionID, simulationInfo.conditionIdSol)]   
        sucess = computeResidualsExpCond!(residuals, odeSolution, θ_dynamicT, θ_sdT, θ_observableT, θ_nonDynamicT,
                                          peTabModel, conditionID, θ_indices, measurementData, parameterInfo)
        if sucess == false
            residuals .= Inf
            break
        end
    end

    return residuals
end


# For an experimental condition compute residuals 
function computeResidualsExpCond!(residuals::AbstractVector, 
                                  odeSolution::Union{ODESolution, OrdinaryDiffEq.ODECompositeSolution},
                                  θ_dynamic::Vector{Float64},
                                  θ_sd::AbstractVector, 
                                  θ_observable::AbstractVector, 
                                  θ_nonDynamic::AbstractVector,
                                  peTabModel::PeTabModel,
                                  conditionID::String,
                                  θ_indices::ParameterIndices,
                                  measurementData::MeasurementData,
                                  parameterInfo::ParameterInfo)::Bool

    if !(odeSolution.retcode == :Success || odeSolution.retcode == :Terminated)
        return false
    end

    # Compute yMod and sd for all observations having id conditionID 
    nModelStates = length(peTabModel.stateNames)
    for i in measurementData.iPerConditionId[conditionID]
        
        t = measurementData.tObs[i]
        u = dualToFloat.(odeSolution[1:nModelStates, measurementData.iTObs[i]])
        hTransformed = computehTransformed(u, t, θ_dynamic, θ_observable, θ_nonDynamic, peTabModel, i, measurementData, θ_indices, parameterInfo)
        σ = computeσ(u, t, θ_dynamic, θ_sd, θ_nonDynamic, peTabModel, i, measurementData, θ_indices, parameterInfo)
            
        # By default a positive ODE solution is not enforced (even though the user can provide it as option).
        # In case with transformations on the data the code can crash, hence Inf is returned in case the 
        # model data transformation can not be perfomred. 
        if isinf(hTransformed)
            return false
        end

        if measurementData.transformData[i] == :lin
            residuals[i] = (hTransformed - measurementData.yObsNotTransformed[i]) / σ
        elseif measurementData.transformData[i] == :log10
            residuals[i] = (hTransformed - measurementData.yObsTransformed[i]) / σ 
        else
            println("Transformation ", measurementData.transformData[i], "not yet supported.")
            return false
        end   
    end
    return true
end
