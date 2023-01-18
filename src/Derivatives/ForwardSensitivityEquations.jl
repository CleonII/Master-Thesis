#=
    Functions specific to gradient compuations via forward sensitivity equations. Notice that we can solve the 
    forward system either via i) solving the expanded ODE-system or ii) by using AutoDiff to obtain the sensitivites, 
    which efficiently are the jacobian of the ODESolution 
=#


function computeGradientForwardEqDynParam!(gradient::Vector{Float64},
                                           θ_dynamic::Vector{Float64},
                                           θ_sd::Vector{Float64},
                                           θ_observable::Vector{Float64},
                                           θ_nonDynamic::Vector{Float64},
                                           S::Matrix{Float64},
                                           peTabModel::PeTabModel,
                                           sensealg::Union{Symbol, SciMLSensitivity.AbstractForwardSensitivityAlgorithm},
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
    success = solveForSensitivites(S, odeProblem, simulationInfo, peTabModel, sensealg, θ_dynamicT, 
                                   solveOdeModelAllConditions!, changeODEProblemParameters!, expIDSolve)
    if success != true
        println("Failed to solve sensitivity equations")
        gradient .= 1e8
        return
    end

    gradient .= 0.0
    # Compute the gradient per experimental conditions. TODO : Find better solution for positionInSolArray 
    for conditionID in keys(measurementData.iPerConditionId)

        if expIDSolve[1] != "all" && conditionID ∉ expIDSolve
            continue
        end

        sol = simulationInfo.solArrayGrad[findfirst(x -> x == conditionID, simulationInfo.conditionIdSol)]
        positionInSolArray = simulationInfo.posInSolArray[conditionID]

        # In case the model is simulated first to a steady state we need to keep track of the post-equlibrium experimental 
        # condition Id to identify parameters specific to an experimental condition.
        postEqulibriumId = simulationInfo.simulateSS == true ? simulationInfo.postEqIdSol[whichForwardSol] : conditionID
        
        # If we have a callback it needs to be properly handled 
        computeGradientForwardExpCond(gradient, sol, S, sensealg, θ_dynamicT, θ_sdT, θ_observableT, θ_nonDynamicT,
                                       conditionID, postEqulibriumId, positionInSolArray, peTabModel, θ_indices, 
                                       measurementData, parameterInfo)
    end
end


function solveForSensitivites(S::Matrix{Float64},
                              odeProblem::ODEProblem, 
                              simulationInfo::SimulationInfo,
                              peTabModel::PeTabModel,
                              sensealg::SciMLSensitivity.AbstractForwardSensitivityAlgorithm,
                              θ_dynamic::AbstractVector, 
                              changeODEProblemParameters!::Function,
                              solveOdeModelAllConditions!::Function,
                              expIDSolve::Vector{String})

    nModelStates = length(peTabModel.stateNames)
    _odeProblem = remake(odeProblem, p = convert.(eltype(θ_dynamic), odeProblem.p), u0 = convert.(eltype(θ_dynamic), odeProblem.u0))
    changeODEProblemParameters!(_odeProblem.p, (@view _odeProblem.u0[1:nModelStates]), θ_dynamic)
    success = solveOdeModelAllConditions!(simulationInfo.solArrayGrad, _odeProblem, θ_dynamic, expIDSolve)
    return success
end
function solveForSensitivites(S::Matrix{Float64},
                              odeProblem::ODEProblem, 
                              simulationInfo::SimulationInfo,
                              peTabModel::PeTabModel,
                              sensealg::Symbol,
                              θ_dynamic::AbstractVector, 
                              changeODEProblemParameters!::Function,
                              solveOdeModelAllConditions!::Function,
                              expIDSolve::Vector{String})

    _odeProblem = remake(odeProblem, p = convert.(eltype(θ_dynamic), odeProblem.p), u0 = convert.(eltype(θ_dynamic), odeProblem.u0))
    success = solveOdeModelAllCondUse!(simulationInfo.solArrayGrad, S, _odeProblem, θ_dynamic, expIDSolve)

    return success
end


function computeGradientForwardExpCond!(gradient::Vector{Float64},
                                        sol::ODESolution,
                                        S::Matrix{Float64},
                                        sensealg::SciMLSensitivity.AbstractForwardSensitivityAlgorithm,
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
    compute∂G∂u = (out, u, p, t, i) -> begin compute∂G∂x(out, u, p, t, i, iPerTimePoint, 
                                                         measurementData, parameterInfo, 
                                                         θ_indices, peTabModel, 
                                                         θ_dynamic, θ_sd, θ_observable, θ_nonDynamic, 
                                                         ∂h∂u, ∂σ∂u, compute∂G∂U=true)
                                            end
    compute∂G∂p = (out, u, p, t, i) -> begin compute∂G∂x(out, u, p, t, i, iPerTimePoint, 
                                                         measurementData, parameterInfo, 
                                                         θ_indices, peTabModel, 
                                                         θ_dynamic, θ_sd, θ_observable, θ_nonDynamic, 
                                                         ∂h∂p, ∂σ∂p, compute∂G∂U=false)
                                        end
      
    # Loop through solution and extract sensitivites                                                 
    p = sol.prob.p
    tSaveAt = measurementData.tVecSave[conditionID]
    ∂G∂p, ∂G∂p_ = zeros(Float64, length(p)), zeros(Float64, length(p)) 
    ∂G∂u = zeros(Float64, length(peTabModel.stateNames))
    _gradient = zeros(Float64, length(p))
    for i in eachindex(tSaveAt)     
        u, _S = extract_local_sensitivities(sol, i, true)
        compute∂G∂u(∂G∂u, u, p, tSaveAt[i], i)
        compute∂G∂p(∂G∂p_, u, p, tSaveAt[i], i)
        _gradient .+= _S'*dgDuOut 
        ∂G∂p .+= ∂G∂p_
    end

    # Thus far have have computed dY/dθ, but for parameters on the log-scale we want dY/dθ_log. We can adjust via;
    # dY/dθ_log = log(10) * θ * dY/dθ
    adjustGradientToTransformedParameters!(gradient, _gradient, ∂G∂p, θ_dynamic, θ_indices, parameterInfo, 
                                           postEqulibriumId, autoDiffSensitivites=false)

end
function computeGradientForwardExpCond!(gradient::Vector{Float64},
                                        sol::ODESolution,
                                        S::Matrix{Float64},
                                        sensealg::Symbol,
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
    compute∂G∂u = (out, u, p, t, i) -> begin compute∂G∂x(out, u, p, t, i, iPerTimePoint, 
                                                         measurementData, parameterInfo, 
                                                         θ_indices, peTabModel, 
                                                         θ_dynamic, θ_sd, θ_observable, θ_nonDynamic, 
                                                         ∂h∂u, ∂σ∂u, compute∂G∂U=true)
                                            end
    compute∂G∂p = (out, u, p, t, i) -> begin compute∂G∂x(out, u, p, t, i, iPerTimePoint, 
                                                         measurementData, parameterInfo, 
                                                         θ_indices, peTabModel, 
                                                         θ_dynamic, θ_sd, θ_observable, θ_nonDynamic, 
                                                         ∂h∂p, ∂σ∂p, compute∂G∂U=false)
                                        end
      
    # Extract which parameters we compute gradient for in this specific experimental condition 
    whichExpMap = findfirst(x -> x == postEqId, [θ_indices.mapExpCond[i].condID for i in eachindex(θ_indices.mapExpCond)])
    iθ_experimentalCondition = vcat(θ_indices.mapDynParEst.iDynParamInVecEst, θ_indices.mapExpCond[whichExpMap])                                                                     
    
    # Loop through solution and extract sensitivites                                                 
    p = dualToFloat.(sol.prob.p)
    nModelStates = length(peTabModel.stateNames)
    tSaveAt = measurementData.tVecSave[conditionID]
    ∂G∂p, ∂G∂p_ = zeros(Float64, length(p)), zeros(Float64, length(p)) 
    ∂G∂u = zeros(Float64, nModelStates)
    _gradient = zeros(Float64, length(θ_indices.iθ_dynamic))
    for i in eachindex(tSaveAt)     
        u = dualToFloat.(sol[:, i])
        compute∂G∂u(∂G∂u, u, p, tSaveAt[i], i)
        compute∂G∂p(∂G∂p_, u, p, tSaveAt[i], i)
        # We need to extract the correct indices from the big sensitivity matrix (row is observation at specific time
        # point). Overall, positions are precomputed in positionInSolArray
        iStart, iEnd = (positionInSolArray[i]-1)*nModelStates+1, (positionInSolArray[i]-1)*nModelStates + nModelStates
        _S = @view S[iStart:iEnd, iθ_experimentalCondition]
        @views _gradient[iθ_experimentalCondition] .+= _S'*∂G∂u 
        ∂G∂p .+= ∂G∂p_
    end

    # Thus far have have computed dY/dθ, but for parameters on the log-scale we want dY/dθ_log. We can adjust via;
    # dY/dθ_log = log(10) * θ * dY/dθ
    adjustGradientToTransformedParameters!(gradient, _gradient, ∂G∂p, θ_dynamic, θ_indices, parameterInfo, 
                                           postEqulibriumId, autoDiffSensitivites=true)
end
