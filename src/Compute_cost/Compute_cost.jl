function computeCost(θ_est::AbstractVector,
                     odeProblem::ODEProblem,  
                     peTabModel::PeTabModel,
                     simulationInfo::SimulationInfo,
                     θ_indices::ParameterIndices,
                     measurementData::MeasurementData,
                     parameterInfo::ParameterInfo,
                     changeODEProblemParameters!::Function,
                     solveOdeModelAllConditions!::Function,
                     priorInfo::PriorInfo;
                     expIDSolve::Array{String, 1} = ["all"],
                     computeHessian::Bool=false, 
                     computeResiduals::Bool=false)::Real

    θ_dynamic, θ_observable, θ_sd, θ_nonDynamic = splitParameterVector(θ_est, θ_indices) 

    cost = computeCostSolveODE(θ_dynamic, θ_sd, θ_observable, θ_nonDynamic, odeProblem, peTabModel, simulationInfo, 
                               θ_indices, measurementData, parameterInfo, changeODEProblemParameters!, solveOdeModelAllConditions!, 
                               computeHessian=computeHessian, computeResiduals=computeResiduals, expIDSolve=expIDSolve)

    if priorInfo.hasPriors == true && computeHessian == false  
        θ_estT = transformθ(θ_est, θ_indices.namesParamEst, parameterInfo)
        cost += evalPriors(θ_estT, θ_est, θ_indices.namesParamEst, θ_indices, priorInfo)
    end

    return cost
end


function computeCostSolveODE(θ_dynamic::AbstractVector,
                             θ_sd::AbstractVector,
                             θ_observable::AbstractVector,
                             θ_nonDynamic::AbstractVector,
                             odeProblem::ODEProblem,
                             peTabModel::PeTabModel,
                             simulationInfo::SimulationInfo,
                             θ_indices::ParameterIndices,
                             measurementData::MeasurementData,
                             parameterInfo::ParameterInfo,
                             changeODEProblemParameters!::Function,
                             solveOdeModelAllConditions!::Function;
                             computeHessian::Bool=false, 
                             computeGradientDynamicθ::Bool=false, 
                             computeResiduals::Bool=false,
                             expIDSolve::Array{String, 1} = ["all"])::Real

    θ_dynamicT = transformθ(θ_dynamic, θ_indices.θ_dynamicNames, parameterInfo)
    θ_sdT = transformθ(θ_sd, θ_indices.θ_sdNames, parameterInfo)
    θ_observableT = transformθ(θ_observable, θ_indices.θ_observableNames, parameterInfo)
    θ_nonDynamicT = transformθ(θ_nonDynamic, θ_indices.θ_nonDynamicNames, parameterInfo)

    _odeProblem = remake(odeProblem, p = convert.(eltype(θ_dynamicT), odeProblem.p), u0 = convert.(eltype(θ_dynamicT), odeProblem.u0))
    changeODEProblemParameters!(_odeProblem.p, _odeProblem.u0, θ_dynamicT)
    
    # If computing hessian or gradient store ODE solution in arrary with dual numbers, else use 
    # solution array with floats
    if computeHessian == true || computeGradientDynamicθ == true
        success = solveOdeModelAllConditions!(simulationInfo.solArrayGrad, _odeProblem, θ_dynamicT, expIDSolve)
    else
        success = solveOdeModelAllConditions!(simulationInfo.solArray, _odeProblem, θ_dynamicT, expIDSolve)
    end
    if success != true
        println("Failed to solve ODE model")
        return Inf
    end

    cost = _computeCost(θ_dynamicT, θ_sdT, θ_observableT, θ_nonDynamicT, peTabModel, simulationInfo, θ_indices, measurementData, 
                        parameterInfo, expIDSolve, 
                        computeHessian=computeHessian, 
                        computeGradientDynamicθ=computeGradientDynamicθ, 
                        computeResiduals=computeResiduals)

    return cost
end


function computeCostNotSolveODE(θ_dynamic::Vector{Float64},
                                θ_sd::AbstractVector,
                                θ_observable::AbstractVector,
                                θ_nonDynamic::AbstractVector,
                                peTabModel::PeTabModel,
                                simulationInfo::SimulationInfo,
                                θ_indices::ParameterIndices,
                                measurementData::MeasurementData,
                                parameterInfo::ParameterInfo;
                                computeGradientNotSolveAutoDiff::Bool=false,
                                computeGradientNotSolveAdjoint::Bool=false,
                                computeGradientNotSolveForward::Bool=false, 
                                expIDSolve::Array{String, 1} = ["all"])::Real 

    # To be able to use ReverseDiff sdParamEstUse and obsParamEstUse cannot be overwritten. 
    # Hence new vectors have to be created. Minimal overhead.
    θ_dynamicT = transformθ(θ_dynamic, θ_indices.θ_dynamicNames, parameterInfo)
    θ_sdT = transformθ(θ_sd, θ_indices.θ_sdNames, parameterInfo)
    θ_observableT = transformθ(θ_observable, θ_indices.θ_observableNames, parameterInfo)
    θ_nonDynamicT = transformθ(θ_nonDynamic, θ_indices.θ_nonDynamicNames, parameterInfo)

    cost = _computeCost(θ_dynamicT, θ_sdT, θ_observableT, θ_nonDynamicT, peTabModel, simulationInfo, θ_indices, 
                        measurementData, parameterInfo, expIDSolve, 
                        computeGradientNotSolveAutoDiff=computeGradientNotSolveAutoDiff, 
                        computeGradientNotSolveAdjoint=computeGradientNotSolveAdjoint, 
                        computeGradientNotSolveForward=computeGradientNotSolveForward)

    return cost
end


function _computeCost(θ_dynamic::AbstractVector,
                      θ_sd::AbstractVector, 
                      θ_observable::AbstractVector, 
                      θ_nonDynamic::AbstractVector,
                      peTabModel::PeTabModel,
                      simulationInfo::SimulationInfo,
                      θ_indices::ParameterIndices,
                      measurementData::MeasurementData, 
                      parameterInfo::ParameterInfo,
                      expIDSolve::Array{String, 1} = ["all"];
                      computeHessian::Bool=false, 
                      computeGradientDynamicθ::Bool=false, 
                      computeResiduals::Bool=false,
                      computeGradientNotSolveAdjoint::Bool=false, 
                      computeGradientNotSolveForward::Bool=false, 
                      computeGradientNotSolveAutoDiff::Bool=false)::Real 

    if computeHessian == true || computeGradientDynamicθ == true || computeGradientNotSolveAdjoint == true || computeGradientNotSolveForward == true || computeGradientNotSolveAutoDiff == true 
        odeSolutionArray = simulationInfo.solArrayGrad
    else
        odeSolutionArray = simulationInfo.solArray
    end

    cost = 0.0
    for conditionID in keys(measurementData.iPerConditionId)

        if expIDSolve[1] != "all" && conditionID ∉ expIDSolve
            continue
        end

        # Extract the ODE-solution for specific condition ID
        odeSolution = odeSolutionArray[findfirst(x -> x == conditionID, simulationInfo.conditionIdSol)]   
        cost += computeCostExpCond(odeSolution, θ_dynamic, θ_sd, θ_observable, θ_nonDynamic, peTabModel, 
                                   conditionID, θ_indices, measurementData, parameterInfo, 
                                   computeResiduals=computeResiduals,
                                   computeGradientNotSolveAdjoint=computeGradientNotSolveAdjoint, 
                                   computeGradientNotSolveForward=computeGradientNotSolveForward, 
                                   computeGradientNotSolveAutoDiff=computeGradientNotSolveAutoDiff)

        if isinf(cost)
            return Inf
        end
    end

    return cost
end


function computeCostExpCond(odeSolution::Union{ODESolution, OrdinaryDiffEq.ODECompositeSolution},
                            θ_dynamic::AbstractVector,
                            θ_sd::AbstractVector, 
                            θ_observable::AbstractVector, 
                            θ_nonDynamic::AbstractVector,
                            peTabModel::PeTabModel,
                            conditionID::String,
                            θ_indices::ParameterIndices,
                            measurementData::MeasurementData,
                            parameterInfo::ParameterInfo; 
                            computeResiduals::Bool=false,
                            computeGradientNotSolveAdjoint::Bool=false, 
                            computeGradientNotSolveForward::Bool=false, 
                            computeGradientNotSolveAutoDiff::Bool=false)::Real

    if !(odeSolution.retcode == :Success || odeSolution.retcode == :Terminated)
        return Inf
    end

    # Compute yMod and sd for all observations having id conditionID 
    cost = 0.0
    for i in measurementData.iPerConditionId[conditionID]
        
        t = measurementData.tObs[i]

        # In these cases we only save the ODE at observed time-points and we do not want 
        # to extract Dual ODE solution 
        if computeGradientNotSolveForward == true || computeGradientNotSolveAutoDiff == true
            nModelStates = length(peTabModel.stateNames)
            u = dualToFloat.(odeSolution[1:nModelStates, measurementData.iTObs[i]]) 
        # For adjoint sensitivity analysis we have a dense-ode solution 
        elseif computeGradientNotSolveAdjoint == true
            # In case we only have sol.t = 0.0 (or similar) interpolation does not work
            u = length(odeSolution.t) > 1 ? odeSolution(t) : odeSolution[1]
        # When we want to extract dual number from the ODE solution 
        else
            u = odeSolution[:, measurementData.iTObs[i]] 
        end

        hTransformed = computehTransformed(u, t, θ_dynamic, θ_observable, θ_nonDynamic, peTabModel, i, measurementData, θ_indices, parameterInfo)
        σ = computeσ(u, t, θ_dynamic, θ_sd, θ_nonDynamic, peTabModel, i, measurementData, θ_indices, parameterInfo)
            
        # By default a positive ODE solution is not enforced (even though the user can provide it as option).
        # In case with transformations on the data the code can crash, hence Inf is returned in case the 
        # model data transformation can not be perfomred. 
        if isinf(hTransformed)
            return Inf
        end

        # Update log-likelihood. In case of guass newton approximation we are only interested in the residuals, and here 
        # we allow the residuals to be computed to test the gauss-newton implementation 
        if computeResiduals == false
            if measurementData.transformData[i] == :lin
                cost += log(σ) + 0.5*log(2*pi) + 0.5*((hTransformed - measurementData.yObsNotTransformed[i]) / σ)^2
            elseif measurementData.transformData[i] == :log10
                cost += log(σ) + 0.5*log(2*pi) + log(log(10)) + log(10)*measurementData.yObsTransformed[i] + 0.5*((hTransformed - measurementData.yObsTransformed[i]) / σ)^2
            else
                println("Transformation ", measurementData.transformData[i], "not yet supported.")
                return Inf
            end   
        else
            if measurementData.transformData[i] == :lin
                cost += ((hTransformed - measurementData.yObsNotTransformed[i]) / σ)
            elseif measurementData.transformData[i] == :log10
                cost += ((hTransformed - measurementData.yObsTransformed[i]) / σ)
            else
                println("Transformation ", measurementData.transformData[i], "not yet supported.")
                return Inf
            end   
        end
    end

    return cost
end