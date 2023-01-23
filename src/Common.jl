# Functions used by both the ODE-solvers and PeTab importer. 


"""
    setParamToFileValues!(paramMap, stateMap, paramData::ParamData)

    Function that sets the parameter and state values in paramMap and stateMap 
    to those in the PeTab parameters file. 

    Used when setting up the PeTab cost function, and when solving the ODE-system 
    for the values in the parameters-file. 
"""
function setParamToFileValues!(paramMap, stateMap, paramData::ParameterInfo)

    parameterNames = paramData.parameterID
    parameterNamesStr = string.([paramMap[i].first for i in eachindex(paramMap)])
    stateNamesStr = replace.(string.([stateMap[i].first for i in eachindex(stateMap)]), "(t)" => "")
    for i in eachindex(parameterNames)
        
        parameterName = parameterNames[i]
        valChangeTo = paramData.paramVal[i]
        
        # Check for value to change to in parameter file 
        i_param = findfirst(x -> x == parameterName, parameterNamesStr)
        i_state = findfirst(x -> x == parameterName, stateNamesStr)

        if !isnothing(i_param)
            paramMap[i_param] = Pair(paramMap[i_param].first, valChangeTo) 
        elseif !isnothing(i_state)
            stateMap[i_state] = Pair(stateMap[i_state].first, valChangeTo)
        end
    end
end


function splitParameterVector(θ_est::AbstractVector, 
                              θ_indices::ParameterIndices)::Tuple{AbstractVector, AbstractVector, AbstractVector, AbstractVector} 

    θ_dynamic = θ_est[θ_indices.iθ_dynamic]
    θ_observable = θ_est[θ_indices.iθ_observable]
    θ_sd = θ_est[θ_indices.iθ_sd]
    θ_nonDynamic = θ_est[θ_indices.iθ_nonDynamic]

    return θ_dynamic, θ_observable, θ_sd, θ_nonDynamic
end


function computeσ(u::AbstractVector,
                  t::Float64,
                  θ_dynamic::AbstractVector,
                  θ_sd::AbstractVector, 
                  θ_nonDynamic::AbstractVector,
                  peTabModel::PeTabModel,
                  iMeasurement::Int64, 
                  measurementData::MeasurementData,
                  θ_indices::ParameterIndices, 
                  parameterInfo::ParameterInfo)::Real

    # Compute associated SD-value or extract said number if it is known 
    if typeof(measurementData.sdParams[iMeasurement]) <: AbstractFloat
        σ = measurementData.sdParams[iMeasurement]
    else
        mapSdParam = θ_indices.mapArraySdParam[θ_indices.indexSdParamMap[iMeasurement]]
        σ = peTabModel.evalSd!(u, t, θ_sd, θ_dynamic, θ_nonDynamic, parameterInfo, measurementData.observebleID[iMeasurement], mapSdParam)
    end

    return σ
end


# Compute observation function h
function computehTransformed(u::AbstractVector,
                             t::Float64,
                             θ_dynamic::AbstractVector,
                             θ_observable::AbstractVector, 
                             θ_nonDynamic::AbstractVector,
                             peTabModel::PeTabModel,
                             iMeasurement::Int64, 
                             measurementData::MeasurementData,
                             θ_indices::ParameterIndices, 
                             parameterInfo::ParameterInfo)::Real

    mapObsParam = θ_indices.mapArrayObsParam[θ_indices.indexObsParamMap[iMeasurement]]
    h = peTabModel.evalYmod(u, t, θ_dynamic, θ_observable, θ_nonDynamic, parameterInfo, measurementData.observebleID[iMeasurement], mapObsParam) 
    # Transform yMod is necessary
    hTransformed = transformObsOrData(h, measurementData.transformData[iMeasurement])

    return hTransformed
end


# Transform parameter from log10 scale to normal scale, or reverse transform
function transformθ!(θ::AbstractVector, 
                     θ_names::Vector{String}, 
                     parameterInfo::ParameterInfo; 
                     reverseTransform::Bool=false) 
    
    @inbounds for i in eachindex(θ)
        iParam = findfirst(x -> x == θ_names[i], parameterInfo.parameterID)
        if isnothing(iParam)
            println("Warning : Could not find parameter ID for ", θ_names[i])
        end
        if parameterInfo.logScale[iParam] == true && reverseTransform == false
            θ[i] = exp10(θ[i])
        elseif parameterInfo.logScale[iParam] == true && reverseTransform == true
            θ[i] = log10(θ[i])
        end
    end
end


# Transform parameter from log10 scale to normal scale, or reverse transform
function transformθ(θ::AbstractVector, 
                    θ_names::Vector{String}, 
                    parameterInfo::ParameterInfo; 
                    reverseTransform::Bool=false)::AbstractVector
    
    iθ = [findfirst(x -> x == θ_names[i], parameterInfo.parameterID) for i in eachindex(θ_names)]
    shouldTransform = [parameterInfo.logScale[i] == true ? true : false for i in iθ]
    shouldNotTransform = .!shouldTransform

    if reverseTransform == false
        out = exp10.(θ) .* shouldTransform .+ θ .* shouldNotTransform
    else
        out = log10.(θ) .* shouldTransform .+ θ .* shouldNotTransform
    end
    return out
end


function changeODEProblemParameters!(pODEProblem::AbstractVector, 
                                     u0::AbstractVector,
                                     θ::AbstractVector,
                                     θ_indices::ParameterIndices,
                                     peTabModel::PeTabModel)

    mapDynParam = θ_indices.mapDynParEst
    pODEProblem[mapDynParam.iDynParamInSys] .= θ[mapDynParam.iDynParamInVecEst]
    peTabModel.evalU0!(u0, pODEProblem) 
    
    return nothing
end


function changeODEProblemParameters(pODEProblem::AbstractVector, 
                                    θ::AbstractVector,
                                    θ_indices::ParameterIndices,
                                    peTabModel::PeTabModel)

    # Helper function to not-inplace map parameters 
    function mapParamToEst(j::Integer, mapDynParam::MapDynParEst)
        whichIndex = findfirst(x -> x == j, mapDynParam.iDynParamInSys)
        return mapDynParam.iDynParamInVecEst[whichIndex]
    end

    mapDynParam = θ_indices.mapDynParEst
    outpODEProblem = [i ∈ mapDynParam.iDynParamInSys ? θ[mapParamToEst(i, mapDynParam)] : pODEProblem[i] for i in eachindex(pODEProblem)]
    outu0Ret = peTabModel.evalU0(outpODEProblem) 
    
    return outpODEProblem, outu0Ret
end


"""
    dualToFloat(x::ForwardDiff.Dual)::Real
    
    Via recursion convert a Dual to a Float.
"""
function dualToFloat(x::ForwardDiff.Dual)::Real
    return dualToFloat(x.value)
end
"""
    dualToFloat(x::AbstractFloat)::AbstractFloat
"""
function dualToFloat(x::AbstractFloat)::AbstractFloat
    return x
end


# Compute prior contribution to log-likelihood 
function computeGradientPrior!(gradient::AbstractVector, 
                               θ::AbstractVector, 
                               θ_indices::ParameterIndices, 
                               priorInfo::PriorInfo, 
                               parameterInfo::ParameterInfo)
            

    _evalPriors = (θ_est) -> begin
                                θ_estT = transformθ(θ_est, θ_indices.namesParamEst, parameterInfo)
                                return evalPriors(θ_estT, θ_est, θ_indices.namesParamEst, θ_indices, priorInfo)
                            end
    gradient .+= ForwardDiff.gradient(_evalPriors, θ)                                
end


# Compute prior contribution to log-likelihood 
function computeHessianPrior!(gradient::AbstractVector, 
                              θ::AbstractVector, 
                              θ_indices::ParameterIndices, 
                              priorInfo::PriorInfo, 
                              parameterInfo::ParameterInfo)

    _evalPriors = (θ_est) -> begin
                                θ_estT =  transformθ(θ_est, θ_indices.namesParamEst, parameterInfo)
                                return evalPriors(θ_estT, θ_est, θ_indices.namesParamEst, θ_indices, priorInfo)
                            end
    gradient .+= ForwardDiff.hessian(_evalPriors, θ)                                
end