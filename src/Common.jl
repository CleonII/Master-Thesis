# Functions used by both the ODE-solvers and PeTab importer. 


"""
    setParamToFileValues!(paramMap, stateMap, paramData::ParamData)

    Function that sets the parameter and state values in paramMap and stateMap 
    to those in the PeTab parameters file. 

    Used when setting up the PeTab cost function, and when solving the ODE-system 
    for the values in the parameters-file. 
"""
function setParamToFileValues!(paramMap, stateMap, paramData::ParamData)

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


function computeσ(u::AbstractVector,
                  t::Float64,
                  θ_dynamic::Vector{Float64},
                  θ_sd::Vector{Float64}, 
                  θ_nonDynamic::Vector{Float64},
                  peTabModel::PeTabModel,
                  iMeasurement::Int64, 
                  θ_indices::ParameterIndices)::Real

    # Compute associated SD-value or extract said number if it is known 
    if typeof(measurementData.sdParams[iMeasurement]) <: AbstractFloat
        σ = measurementData.sdParams[iMeasurement]
    else
        mapSdParam = θ_indices.mapArraySdParam[θ_indices.indexSdParamMap[iMeasurement]]
        σ = peTabModel.evalSd!(u, t, θ_sd, θ_dynamic, θ_nonDynamic, parameterData, measurementData.observebleID[iMeasurement], mapSdParam)
    end

    return σ
end


# Compute observation function h
function computehTransformed(u::AbstractVector,
                             t::Float64,
                             θ_dynamic::Vector{Float64},
                             θ_observable::Vector{Float64}, 
                             θ_nonDynamic::Vector{Float64},
                             peTabModel::PeTabModel,
                             iMeasurement::Int64, 
                             θ_indices::ParameterIndices, 
                             parameterInfo::ParameterInfo)::Real

    mapObsParam = θ_indices.mapArrayObsParam[θ_indices.indexObsParamMap[iMeasurement]]
    h = peTabModel.evalYmod(u, t, θ_dynamic, θ_observable, θ_nonDynamic, parameterInfo, measurementData.observebleID[iMeasurement], mapObsParam) 
    # Transform yMod is necessary
    hTransformed = transformObsOrData(h, measurementData.transformData[iMeasurement])

    return hTransformed
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
