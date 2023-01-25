# Functions used by both the ODE-solvers and PeTab importer. 


"""
    setParamToFileValues!(paramMap, stateMap, paramData::ParamData)

    Function that sets the parameter and state values in paramMap and stateMap 
    to those in the PeTab parameters file. 

    Used when setting up the PeTab cost function, and when solving the ODE-system 
    for the values in the parameters-file. 
"""
function setParamToFileValues!(paramMap, stateMap, paramData::ParametersInfo)

    parameterNames = paramData.parameterId
    parameterNamesStr = string.([paramMap[i].first for i in eachindex(paramMap)])
    stateNamesStr = replace.(string.([stateMap[i].first for i in eachindex(stateMap)]), "(t)" => "")
    for i in eachindex(parameterNames)
        
        parameterName = parameterNames[i]
        valChangeTo = paramData.nominalValue[i]
        
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
                  measurementInfo::MeasurementsInfo,
                  θ_indices::ParameterIndices, 
                  parameterInfo::ParametersInfo)::Real

    # Compute associated SD-value or extract said number if it is known 
    mapθ_sd = θ_indices.mapθ_sd[iMeasurement]
    if mapθ_sd.isSingleConstant == true
        σ = mapθ_sd.constantValues[1]
    else
        σ = peTabModel.evalSd!(u, t, θ_sd, θ_dynamic, θ_nonDynamic, parameterInfo, measurementInfo.observableId[iMeasurement], mapθ_sd)
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
                             measurementInfo::MeasurementsInfo,
                             θ_indices::ParameterIndices, 
                             parameterInfo::ParametersInfo)::Real

    mapθ_observable = θ_indices.mapθ_observable[iMeasurement]
    h = peTabModel.evalYmod(u, t, θ_dynamic, θ_observable, θ_nonDynamic, parameterInfo, measurementInfo.observableId[iMeasurement], mapθ_observable) 
    # Transform yMod is necessary
    hTransformed = transformObsOrData(h, measurementInfo.measurementTransformation[iMeasurement])

    return hTransformed
end


# Function to extract observable or noise parameters when computing h or σ
function getObsOrSdParam(θ::AbstractVector, parameterMap::θObsOrSdParameterMap)

    # Helper function to map SD or obs-parameters in non-mutating way 
    function map1Tmp(iVal)
        whichI = sum(parameterMap.shouldEstimate[1:iVal])
        return parameterMap.indexInθ[whichI]
    end
    function map2Tmp(iVal)
        whichI = sum(.!parameterMap.shouldEstimate[1:iVal])
        return whichI
    end

    # In case of no SD/observable parameter exit function
    if parameterMap.nParameters == 0
        return
    end

    # In case of single-value return do not have to return an array and think about type
    if parameterMap.nParameters == 1 
        if parameterMap.shouldEstimate[1] == true
            return θ[parameterMap.indexInθ][1]
        else
            return parameterMap.constantValues[1]
        end
    end

    nParametersToEstimte = sum(parameterMap.shouldEstimate)
    if nParametersToEstimte == parameterMap.nParameters
        return θ[parameterMap.indexInθ]

    elseif nParametersToEstimte == 0
        return parameterMap.constantValues        

    # Computaionally most demanding case. Here a subset of the parameters 
    # are to be estimated. This code must be non-mutating to support Zygote which 
    # negatively affects performance 
    elseif nParametersToEstimte > 0
        _values = [parameterMap.shouldEstimate[i] == true ? θ[map1Tmp(i)] : 0.0 for i in 1:parameterMap.nParameters]
        values = [parameterMap.shouldEstimate[i] == false ? parameterMap.constantValues[map2Tmp(i)] : _values[i] for i in 1:parameterMap.nParameters]
        return values
    end
end


# Transform parameter from log10 scale to normal scale, or reverse transform
function transformθ!(θ::AbstractVector, 
                     θ_names::Vector{Symbol}, 
                     parameterInfo::ParametersInfo; 
                     reverseTransform::Bool=false) 
    
    @inbounds for i in eachindex(θ)
        iParam = findfirst(x -> x == θ_names[i], parameterInfo.parameterId)
        if isnothing(iParam)
            println("Warning : Could not find parameter ID for ", θ_names[i])
        end
        if parameterInfo.parameterScale[iParam] == :lin
            θ[i] = θ[i]
        elseif parameterInfo.parameterScale[iParam] == :log10
            θ[i] = reverseTransform == true ? log10(θ[i]) : exp10(θ[i]) 
        elseif parameterInfo.parameterScale[iParam] == :log
            θ[i] = reverseTransform == true ? log(θ[i]) : exp(θ[i]) 
        end
    end
end


# Transform parameter from log10 scale to normal scale, or reverse transform
function transformθ(θ::AbstractVector, 
                    θ_names::Vector{Symbol}, 
                    parameterInfo::ParametersInfo; 
                    reverseTransform::Bool=false)::AbstractVector
    
    iθ = [findfirst(x -> x == θ_names[i], parameterInfo.parameterId) for i in eachindex(θ_names)]
    shouldTransform = [parameterInfo.parameterScale[i] == :log10 ? true : false for i in iθ]
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

    mapODEProblem = θ_indices.mapODEProblem
    pODEProblem[mapODEProblem.iODEProblemθDynamic] .= θ[mapODEProblem.iθDynamic]
    peTabModel.evalU0!(u0, pODEProblem) 
    
    return nothing
end


function changeODEProblemParameters(pODEProblem::AbstractVector, 
                                    θ::AbstractVector,
                                    θ_indices::ParameterIndices,
                                    peTabModel::PeTabModel)

    # Helper function to not-inplace map parameters 
    function mapParamToEst(j::Integer, mapDynParam::MapODEProblem)
        whichIndex = findfirst(x -> x == j, mapDynParam.iODEProblemθDynamic)
        return mapODEProblem.iθDynamic[whichIndex]
    end

    mapODEProblem = θ_indices.mapODEProblem
    outpODEProblem = [i ∈ mapODEProblem.iODEProblemθDynamic ? θ[mapParamToEst(i, mapODEProblem)] : pODEProblem[i] for i in eachindex(pODEProblem)]
    outu0 = peTabModel.evalU0(outpODEProblem) 
    
    return outpODEProblem, outu0
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


function getFileODEvalues(peTabModel::PeTabModel)
  
    # Change model parameters 
    experimentalConditionsFile, measurementDataFile, parameterDataFile, observablesDataFile = readPEtabFiles(peTabModel.dirModel, readObservables=true)
    parameterInfo = processParameters(parameterDataFile)
    measurementInfo = processMeasurements(measurementDataFile, observablesDataFile) 
    θ_indices = computeIndicesθ(parameterInfo, measurementInfo, peTabModel.odeSystem, experimentalConditionsFile)

    θ_estNames = θ_indices.θ_estNames
    θ_est = parameterInfo.nominalValue[findall(x -> x ∈ θ_estNames, parameterInfo.parameterId)]

    return θ_est[θ_indices.iθ_dynamic]
end