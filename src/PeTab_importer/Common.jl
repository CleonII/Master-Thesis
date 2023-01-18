"""
    isNumber(x::String)::Bool

    Check if a string x is a number (Float).
"""
function isNumber(x::AbstractString)::Bool
    re1 = r"^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)$" # Picks up scientific notation
    re2 = r"^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)$"
    return (occursin(re1, x) || occursin(re2, x))
end
"""
    isNumber(x::SubString{String})::Bool
"""
function isNumber(x::SubString{String})::Bool
    re1 = r"^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)$" # Picks up scientific notation
    re2 = r"^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)$"
    return (occursin(re1, x) || occursin(re2, x))
end


"""
    transformYobsOrYmodArr!(vals, transformationArr::Array{Symbol, 1})

    Transform the Yobs or Ymod arrays (vals) in place using for each value 
    in vals the transformation specifed in transformationArr.

    Currently :lin, :log10 and :log transforamtions are supported, see
    `setUpCostFunc`.
"""
function transformYobsOrYmodArr!(vals, transformationArr::Array{Symbol, 1})
    for i in eachindex(vals)
        vals[i] = transformObsOrData(vals[i], transformationArr[i])
    end
end


"""
    transformYobsOrYmodArr!(vals, transformationArr::Array{Symbol, 1})

    Transform val using either :lin (identify), :log10 and :log transforamtions.
"""
function transformObsOrData(val, transform::Symbol)
    if transform == :lin
        return val
    elseif transform == :log10
        return val > 0 ? log10(val) : Inf
    elseif transform == :log
        return val > 0 ? log(val) : Inf
    else
        println("Error : $transform is not an allowed transformation")
        println("Only :lin, :log10 and :log are supported.")
    end
end


function splitParameterVector(θ_est::AbstractVector, 
                              θ_indices::ParameterIndices)::Tuple{AbstractVector, AbstractVector, AbstractVector, AbstractVector} 

    θ_dynamic = θ_est[θ_indices.iθ_dynamic]
    θ_observable = θ_est[θ_indices.iθ_observable]
    θ_sd = θ_est[θ_indices.iθ_sd]
    θ_nonDynamic = θ_est[θ_indices.iθ_sd]

    return θ_dynamic, θ_observable, θ_sd, θ_nonDynamic
end


function couldSolveODEModel(simulationInfo::SimulationInfo, expIDSolve::Vector{String})::Bool
    @inbounds for i in eachindex(simulationInfo.solArrayGrad)
        if expIDSolve[1] == "all" || simulationInfo.conditionIdSol[i] ∈ expIDSolve
            if simulationInfo.solArrayGrad[i].retcode != :Success
                return false
            end
        end
    end
    return true
end


function getIndicesParametersNotInODESystem(θ_indices::ParameterIndices)::Tuple

    θ_observableNames = θ_indices.θ_observableNames
    θ_sdNames = θ_indices.θ_sdNames
    θ_nonDynamicNames = θ_indices.θ_nonDynamicNames
    iθ_notOdeSystemNames = θ_indices.iθ_notOdeSystemNames

    iθ_sd = [findfirst(x -> x == θ_sdNames[i], iθ_notOdeSystemNames) for i in eachindex(θ_sdNames)]
    iθ_observable = [findfirst(x -> x == θ_observableNames[i],  iθ_notOdeSystemNames) for i in eachindex(θ_observableNames)]
    iθ_nonDynamic = [findfirst(x -> x == θ_nonDynamicNames[i],  iθ_notOdeSystemNames) for i in eachindex(θ_nonDynamicNames)]
    iθ_notOdeSystem = θ_indices.iθ_notOdeSystem

    return iθ_sd, iθ_observable, iθ_nonDynamic, iθ_notOdeSystem
end


# Allocate derivates needed when computing ∂G∂u and ∂G∂p
function allocateObservableFunctionDerivatives(sol::ODESolution, peTabModel::PeTabModel)

    nModelStates = length(peTabModel.stateNames)
    ∂h∂u = zeros(Float64, nModelStates)
    ∂σ∂u = zeros(Float64, nModelStates)
    ∂h∂p = zeros(Float64, length(sol.prob.p))
    ∂σ∂p = zeros(Float64, length(sol.prob.p))
    return ∂h∂u, ∂σ∂u, ∂h∂p, ∂σ∂p
end


function adjustGradientTransformedParameters!(gradient::AbstractVector, 
                                              _gradient::AbstractVector, 
                                              ∂G∂p::AbstractVector, 
                                              θ_dynamic::Vector{Float64},
                                              θ_indices::ParameterIndices,
                                              parameterInfo::ParameterInfo, 
                                              postEqulibriumId::String; 
                                              autoDiffSensitivites::Bool=false)

    # In case we compute the sensitivtes via automatic differentation the parameters in _gradient=S'*∂G∂u will appear in the 
    # same order as they appear in θ_est. In case we do not compute sensitivtes via autodiff, or do adjoint sensitity analysis, 
    # the parameters in _gradient=S'∂G∂u appear in the same order as in odeProblem.p.
    if autoDiffSensitivites == true
        i_gradient1 = θ_indices.mapDynParEst.iDynParamInVecEst
        i_gradient2 = expMap.iDynEstVec
    else
        i_gradient1 = θ_indices.mapDynParEst.iDynParamInSys
        i_gradient2 = expMap.iOdeProbDynParam
    end
    
    # Transform gradient parameter that for each experimental condition appear in the ODE system                                                             
    gradient[θ_indices.mapDynParEst.iDynParamInVecEst] .+= transformGradient(_gradient[i_gradient] .+ ∂G∂p[θ_indices.mapDynParEst.iDynParamInSys], 
                                                                             θ_dynamic[θ_indices.mapDynParEst.iDynParamInVecEst], 
                                                                             θ_indices.θ_dynamicNames[θ_indices.mapDynParEst.iDynParamInVecEst], 
                                                                             parameterInfo)
    
    # Transform gradient for parameters which are specific to certain experimental conditions. 
    whichExpMap = findfirst(x -> x == postEqulibriumId, [θ_indices.mapExpCond[i].condID for i in eachindex(paramIndices.mapExpCond)])
    expMap = θ_indices.mapExpCond[whichExpMap]                                          
    gradient[expMap.iDynEstVec] .+= transformParamVecGrad(_gradient[i_gradient2] .+ ∂G∂p[expMap.iOdeProbDynParam], 
                                                          θ_dynamic[expMap.iDynEstVec], 
                                                          θ_indices.θ_dynamicNames[expMap.iDynEstVec], 
                                                          parameterInfo)                                   

end