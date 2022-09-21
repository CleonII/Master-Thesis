# Functions for efficently mapping parameters. Specifcially, mapping the optimizsing 
# parameter vector to dynamic, observable and sd parameter vectors, and for mapping 
# these three vectors efficently when computing the observables, u0 and sd-values.


"""
    getIndicesParam(paramData::ParamData, measurementData::MeasurementData)::ParameterIndices

    For a PeTab-model creates index vectors for extracting the dynamic-, observable- and sd-parameters from 
    the  optimizsing parameter vector, and extract names for all parameters and stores the information in 
    the ParameterIndices-struct.

    TODO: Add index for mapping observables

    See also: [`ParameterIndices`]
"""
function getIndicesParam(paramData::ParamData, measurementData::MeasurementData)::ParameterIndices

    namesObsParam = getIdEst(measurementData.obsParam, paramData)
    isObsParam = [paramData.parameterID[i] in namesObsParam for i in eachindex(paramData.parameterID)]

    namesSdParam = getIdEst(measurementData.sdParams, paramData)
    isSdParam = [paramData.parameterID[i] in namesSdParam for i in eachindex(paramData.parameterID)]

    isDynamicParam = (paramData.shouldEst .&& .!isSdParam .&& .!isObsParam)
    namesParamDyn = paramData.parameterID[isDynamicParam]

    namesParamEst::Array{String, 1} = string.(vcat(string.(namesParamDyn), string.(namesSdParam), string.(namesObsParam)))
    
    # Index vector for the dynamic and sd parameters as UInt32 vectors 
    iDynPar::Array{UInt32, 1} = convert(Array{UInt32, 1}, collect(1:length(namesParamDyn)))
    iSdPar::Array{UInt32, 1} = convert(Array{UInt32, 1}, collect((length(namesParamDyn)+1):(length(namesParamDyn) + length(namesSdParam))))
    iObsPar::Array{UInt32, 1} = convert(Array{UInt32, 1}, collect((length(namesParamDyn) + length(namesSdParam) + 1):(length(namesParamDyn) + length(namesSdParam) + length(namesObsParam))))
    
    paramIndicies = ParameterIndices(iDynPar, 
                                     iObsPar, 
                                     iSdPar, 
                                     string.(namesParamDyn), 
                                     string.(namesObsParam), 
                                     string.(namesSdParam),
                                     namesParamEst)

    return paramIndicies
end


# Function for extracting Observable parameter when computing the cost. Will be heavily altered as slow.
function getObsOrSdParam(vecParam,
                         paramData::ParamData, 
                         measurementData::MeasurementData, 
                         observableId::String,
                         simulationConditionId::String,
                         t::Float64;
                         getObsPar::Bool=true)

    if getObsPar == true
        idVec = getIdEst(measurementData.obsParam, paramData)
    else
        idVec = getIdEst(measurementData.sdParams, paramData)
    end

    nObs = length(measurementData.tObs)
    # Acquire which observation is being handled (for really large parameter files smarter data structs might be needed)
    iUse = findfirst(x -> measurementData.observebleID[x] == observableId && measurementData.conditionId[x] == simulationConditionId && measurementData.tObs[x] == t, 1:nObs)
    if isnothing(iUse)
        println("Warning : Cannot identify an observation in getObsParam ")
    end
    if getObsPar == true && isempty(measurementData.obsParam[iUse])
        return nothing
    elseif getObsPar == false && isempty(measurementData.sdParams[iUse])
        return nothing
    end

    if getObsPar == true
        paramsRet = split(measurementData.obsParam[iUse], ';')
    else
        paramsRet = split(measurementData.sdParams[iUse], ';')
    end

    output = Array{eltype(vecParam), 1}(undef, length(paramsRet))
    for i in eachindex(paramsRet)
        # Hard coded constant number 
        if isNumber(paramsRet[i])
            output[i] = parse(Float64, paramsRet[i])

        # Parameter to be estimated 
        elseif paramsRet[i] in idVec
            output[i] = vecParam[findfirst(x -> x == paramsRet[i], idVec)]
        
        # Constant parameter 
        elseif paramsRet[i] in paramData.parameterID
            output[i] = paramData.paramVal[findfirst(x -> x == paramsRet[i], paramData.parameterID)]
        else
            println("Warning : cannot find matching for ", paramsRet[i])
        end
    end

    if length(output) == 1
        return output[1]
    else
        return output
    end
end


# Helper function for extracting Observable parameter when computing the cost. Will be heavily altered as slow.
function getIdEst(idsInStr::Array{String, 1}, paramData::ParamData)
    idVec = String[]

    for i in eachindex(idsInStr)
        if isempty(idsInStr[i])
            continue
        else
            idsInStrSplit = split(idsInStr[i], ';')
            for idInStr in idsInStrSplit

                # Disregard Id if parameters should not be estimated, or 
                iParam = findfirst(x -> x == idInStr, paramData.parameterID)
                if isNumber(idInStr)
                    continue
                elseif isnothing(iParam)                    
                    println("Warning : param $idInStr could not be found in parameter file")
                elseif idInStr in idVec
                    continue
                elseif paramData.shouldEst[iParam] == false
                    continue
                else
                    idVec = vcat(idVec, idInStr)
                end
            end
        end
    end

    return idVec
end