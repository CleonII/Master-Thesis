
# Compute cost (likelihood) for Zygote, which means only using out-of-place functions 

function computeCostZygote(θ_est,
                           odeProblem::ODEProblem,  
                           peTabModel::PeTabModel,
                           simulationInfo::SimulationInfo,
                           θ_indices::ParameterIndices,
                           measurementData::MeasurementData,
                           parameterInfo::ParameterInfo,
                           changeODEProblemParameters::Function,
                           solveOdeModelAllConditions::Function,
                           priorInfo::PriorInfo)

    θ_dynamic, θ_observable, θ_sd, θ_nonDynamic = splitParameterVector(θ_est, θ_indices)                     

    cost = _computeCostZygote(θ_dynamic, θ_sd, θ_observable, θ_nonDynamic, odeProblem,
                              peTabModel, simulationInfo, θ_indices, measurementData,
                              parameterInfo, changeODEProblemParameters, solveOdeModelAllConditions, 
                              computeGradientDynamicθ=false)

    if priorInfo.hasPriors == true
        θ_estT = transformθ(θ_est, θ_indices.namesParamEst, parameterInfo)
        cost += evalPriors(θ_estT, θ_est, θ_indices.namesParamEst, θ_indices, priorInfo)
    end                                  

    return cost                          
end


# Computes the likelihood in such a in a Zygote compatible way, which mainly means that no arrays are mutated.
function _computeCostZygote(θ_dynamic,
                            θ_sd,
                            θ_observable,
                            θ_nonDynamic,
                            odeProblem::ODEProblem,  
                            peTabModel::PeTabModel,
                            simulationInfo::SimulationInfo,
                            θ_indices::ParameterIndices,
                            measurementData::MeasurementData,
                            paramterInfo::ParameterInfo,
                            changeODEProblemParameters::Function,
                            solveOdeModelAllConditions::Function;
                            computeGradientDynamicθ::Bool=false)::Real

    θ_dynamicT = transformθ(θ_dynamic, θ_indices.θ_dynamicNames, parameterInfo)
    θ_sdT = transformθ(θ_sd, θ_indices.θ_sdNames, parameterInfo)
    θ_observableT = transformθ(θ_observable, θ_indices.θ_observableNames, parameterInfo)
    θ_nonDynamicT = transformθ(θ_nonDynamic, θ_indices.θ_nonDynamicNames, parameterInfo)
                                                            
    _p, _u0 = changeODEProblemParameters(odeProblem.p, θ_dynamicT)
    _odeProblem = remake(odeProblem, p = convert.(eltype(θ_dynamic), _p), u0 = convert.(eltype(θ_dynamic), _u0))
    
    # Compute yMod and sd-val by looping through all experimental conditons. At the end 
    # update the likelihood 
    cost = convert(eltype(θ_dynamic), 0.0)
    for conditionID in keys(measurementData.iPerConditionId)
        
        tMax = simulationInfo.tMaxForwardSim[findfirst(x -> x == conditionID, simulationInfo.conditionIdSol)]
        odeSolution, success = solveOdeModelAllConditions(_odeProblem, conditionID, θ_dynamicT, tMax)
        if success != true
            return Inf
        end

        cost += _computeCost(odeSolution, θ_dynamicT, θ_sdT, θ_observableT, θ_nonDynamicT, peTabModel, 
                             conditionID, θ_indices, measurementData, paramterInfo, 
                             computeGradientNotSolveAdjoint=computeGradientDynamicθ)

        if isinf(cost)
            return cost
        end
    end

    return cost
end