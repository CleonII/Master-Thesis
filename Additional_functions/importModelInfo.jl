

struct ModelData
    observableLogTransformation::Vector{Bool}
    optParameterIndices::Vector{Int64}
    inputParameterIndices::Vector{Int64}
    initVariableIndices::Vector{Int64}
    observableVariableIndices::Vector{Int64}
    parameterInObservableIndices::Vector{Int64}
    parameterInU0Indices::Vector{Int64}
end

function ModelData(new_sys, prob, observables, experimentalConditions, initVariableNames, observableVariableNames, parameterInU0Names, parameterInObservableNames)
    observableLogTransformation = [observables[i, :observableTransformation] == "log10" for i in eachindex(observables[!, :observableTransformation])]

    variableNames = string.(states(new_sys))
    numUsedParameters = length(prob.p)
    problemParameterSymbols = parameters(new_sys)
    inputParameterSymbols = Symbol.(names(experimentalConditions)[2:end])
    tmp = [pPS.name for pPS in problemParameterSymbols]
    inputParameterSymbols = inputParameterSymbols[[any(inpPar .== tmp) for inpPar in inputParameterSymbols]]
    inputParameterIndices = [findfirst([pPS.name == inputParameterSymbols[i] for pPS in problemParameterSymbols]) for i in 1:length(inputParameterSymbols)]
    optParameterIndices = collect(1:numUsedParameters)[[pPS.name ∉ inputParameterSymbols for pPS in problemParameterSymbols]]
    parameterNames = string.(problemParameterSymbols)

    initVariableNames = initVariableNames .* "(t)"
    initVariableIndices = [findfirst(initVar .== variableNames) for initVar in initVariableNames]

    observableVariableNames = observableVariableNames .* "(t)"
    observableVariableIndices = [findfirst(observableVar .== variableNames) for observableVar in observableVariableNames]

    parameterInU0Indices = [findfirst(parIU0N .== parameterNames) for parIU0N in parameterInU0Names]

    parameterInObservableIndices = [findfirst(parION .== parameterNames) for parION in parameterInObservableNames]

    modelData = ModelData(observableLogTransformation, optParameterIndices, inputParameterIndices, 
        initVariableIndices, observableVariableIndices, parameterInObservableIndices, parameterInU0Indices)

    return modelData, inputParameterSymbols
end


struct ExperimentalData
    measurementForCondObs::Array{Vector{Float64}, 2}
    timeStepsForCond::Vector{Vector{Float64}}
    observedAtIndexForCondObs::Array{Vector{Int64}, 2}
    numConditions::Int64
    numObservables::Int64
    numDataForCondObs::Array{Int64, 2}
    inputParameterValuesForCond::Vector{Vector{Float64}}
    observedObservableForCond::Vector{Vector{Int64}}
    observablesTimeIndexIndicesForCond::Vector{Array{Vector{Int64}, 2}}
end

function ExperimentalData(observables, experimentalConditions, measurementData, inputParameterSymbols)
    observableNames = observables[!, 1]

    numConditions = length(experimentalConditions[!,1])
    numObservables = length(observableNames) 

    ## Setting up for different conditions
    inputParameterValuesForCond = Vector{Vector{Float64}}(undef, numConditions)

    timeStepsForCond = Vector{Vector{Float64}}(undef, numConditions)
    numTimeStepsForCond = Vector{Int64}(undef, numConditions)
    observedObservableForCond = Vector{Vector{Int64}}(undef, numConditions)

    observedAtForCond = Vector{Vector{Vector{Float64}}}(undef, numConditions)
    observedAtIndexForCondObs = Array{Vector{Int64}, 2}(undef, numConditions, numObservables)
    measurementForCondObs = Array{Vector{Float64}, 2}(undef, numConditions, numObservables)
    numDataForCondObs = Array{Int64, 2}(undef, numConditions, numObservables)
    observablesTimeIndexIndicesForCond = Vector{Array{Vector{Int64}, 2}}(undef, numConditions)

    for (iCond, condId) in enumerate(experimentalConditions[!,1])
        inputParameterValuesForCond[iCond] = [experimentalConditions[iCond, inPar] for inPar in inputParameterSymbols]

        relevantMeasurementData = measurementData[measurementData[:,3] .== condId, :]
        observableIDs = relevantMeasurementData[:,1]
        timeStepsForCond[iCond] = sort(unique(relevantMeasurementData[:,5]))
        numTimeStepsForCond[iCond] = length(timeStepsForCond[iCond])
        isObserved = [observableNames[i] in unique(observableIDs) for i=1:numObservables]
        observedObservableForCond[iCond] = collect(1:numObservables)[isObserved]

        observedAtForCond[iCond] = Vector{Vector{Float64}}(undef, numObservables)
        for (iObs, obsId) = enumerate(observableNames)
            observedAtForCond[iCond][iObs] = relevantMeasurementData[observableIDs .== obsId, 5]
        end
        for iObs = 1:numObservables
            observedAtIndexForCondObs[iCond, iObs] = indexin(observedAtForCond[iCond][iObs], timeStepsForCond[iCond])
        end
        for (iObs, obsId) = enumerate(observableNames)
            measurementForCondObs[iCond, iObs] = relevantMeasurementData[observableIDs .== obsId, 4]
            numDataForCondObs[iCond, iObs] = length(measurementForCondObs[iCond, iObs])
        end
        observablesTimeIndexIndicesForCond[iCond] = Array{Vector{Int64}, 2}(undef, numObservables, numTimeStepsForCond[iCond])
        for i = 1:numObservables
            for j = 1:numTimeStepsForCond[iCond]
                observablesTimeIndexIndicesForCond[iCond][i,j] = findall(x->x==j, observedAtIndexForCondObs[iCond, i])
            end
        end
    end

    experimentalData = ExperimentalData(measurementForCondObs, timeStepsForCond, observedAtIndexForCondObs, 
            numConditions, numObservables, numDataForCondObs, inputParameterValuesForCond, 
            observedObservableForCond, observablesTimeIndexIndicesForCond)

    return experimentalData
end


struct ModelParameters
    scaleIndices::Vector{Int64}
    offsetIndices::Vector{Int64}
    varianceIndices::Vector{Int64}
    parameterIndices::Vector{Int64}
    #
    scaleMap::Array{Int64, 2}
    offsetMap::Array{Int64, 2}
    varianceMap::Array{Int64, 2}
    #
    numScale::Int64
    numOffset::Int64
    numVariance::Int64
    numOptParameters::Int64
    numUsedParameters::Int64
    numAllParameters::Int64
    numVariables::Int64
    #
    scaleNames::Vector{String}
    offsetNames::Vector{String}
    varianceNames::Vector{String}
    optParameterNames::Vector{String}
    #
    scaleVector::Vector{Float64}
    offsetVector::Vector{Float64}
    varianceVector::Vector{Float64}
    dynamicParametersVector::Vector{Float64}
    u0Vector::Vector{Float64}
    allParameters::Vector{Float64}
end


function ModelParameters(new_sys, prob, parameterBounds, experimentalConditions, measurementData, observables, experimentalData; 
    scaleDeterminer = "scale", offsetDeterminer = "offset", varianceDeterminer = "sd_")
    numConditions = experimentalData.numConditions
    numObservables = experimentalData.numObservables
    observableNames = observables[!, 1]

    ## Creating a map from a condition, observable pair to a scale ##
    scaleNames = parameterBounds[findall(occursin.(scaleDeterminer, parameterBounds[!, 1])), 1]
    numScale = length(findall(occursin.(scaleDeterminer, parameterBounds[!, 1])))
    # scale[numScale+1] = 1
    scaleMap = Array{Int64, 2}(undef, numConditions, numObservables) # Have as sparse?
    for (iCond, Cond) in enumerate(experimentalConditions[!, 1])
        for (iObs, Obs) in enumerate(observableNames)
            measurementIndex = findfirst((measurementData[!, 1] .== Obs) .& (measurementData[!, 3] .== Cond))
            if measurementIndex === nothing
                scaleMap[iCond, iObs] = numScale+1
            else
                obsPar = measurementData[measurementIndex, 6]
                if obsPar === missing
                    scaleMap[iCond, iObs] = numScale+1
                else
                    if occursin(';', obsPar)
                        obsPars = split(obsPar, ';')
                        for op in obsPars
                            if occursin(scaleDeterminer, op)
                                obsPar = op 
                            end
                        end
                    end
                    if occursin(scaleDeterminer, obsPar)
                        scaleIndex = findfirst(scaleNames .== obsPar)
                        scaleMap[iCond, iObs] = scaleIndex
                    else
                        scaleMap[iCond, iObs] = numScale+1
                    end
                end
            end
        end
    end

    ## Creating a map from a condition, observable pair to an offset ##
    offsetNames = parameterBounds[findall(occursin.(offsetDeterminer, parameterBounds[!, 1])), 1]
    numOffset = length(findall(occursin.(offsetDeterminer, parameterBounds[!, 1])))
    # offset[numOffset+1] = 0, offset[numOffset+2] = 1
    offsetMap = Array{Int64, 2}(undef, numConditions, numObservables)
    for (iCond, Cond) in enumerate(experimentalConditions[!, 1])
        for (iObs, Obs) in enumerate(observableNames)
            measurementIndex = findfirst((measurementData[!, 1] .== Obs) .& (measurementData[!, 3] .== Cond))
            if measurementIndex === nothing
                offsetMap[iCond, iObs] = numOffset+1
            else
                obsPar = measurementData[measurementIndex, 6]
                if obsPar === missing
                    if observables[iObs, :observableFormula][end-2:end] == "+ 1"
                        offsetMap[iCond, iObs] = numOffset+2
                    else
                        offsetMap[iCond, iObs] = numOffset+1
                    end
                else
                    if occursin(';', obsPar)
                        obsPar = split(obsPar, ';')[1]
                        if occursin(offsetDeterminer, obsPar)
                            offsetIndex = findfirst(offsetNames .== obsPar)
                            offsetMap[iCond, iObs] = offsetIndex
                        elseif obsPar == "0"
                            offsetMap[iCond, iObs] = numOffset+1
                        else
                            offsetMap[iCond, iObs] = numOffset+2
                        end
                    else
                        if occursin(offsetDeterminer, obsPar)
                            offsetIndex = findfirst(offsetNames .== obsPar)
                            offsetMap[iCond, iObs] = offsetIndex
                        elseif observables[iObs, :observableFormula][end-2:end] == "+ 1"
                            if obsPar == "0"
                                offsetMap[iCond, iObs] = numOffset+1
                            else
                                offsetMap[iCond, iObs] = numOffset+2
                            end
                        else
                            offsetMap[iCond, iObs] = numOffset+1
                        end
                    end
                end
            end
        end
    end

    ## Creating a map from a condition, observable pair to a variance ##
    varianceNames = parameterBounds[findall(occursin.(varianceDeterminer, parameterBounds[!, 1])), 1]
    numVariance = length(findall(occursin.(varianceDeterminer, parameterBounds[!, 1])))
    varianceMap = Array{Int64, 2}(undef, numConditions, numObservables)
    for (iCond, Cond) in enumerate(experimentalConditions[!, 1])
        for (iObs, Obs) in enumerate(observableNames)
            measurementIndex = findfirst((measurementData[!, 1] .== Obs) .& (measurementData[!, 3] .== Cond))
            if measurementIndex === nothing
                varianceMap[iCond, iObs] = numVariance+1
            else
                obsPar = measurementData[measurementIndex, 7]
                if obsPar === missing
                    varianceMap[iCond, iObs] = numVariance+1
                else
                    if occursin(';', obsPar)
                        obsPars = split(obsPar, ';')
                        for op in obsPars
                            if occursin(varianceDeterminer, op)
                                obsPar = op 
                            end
                        end
                    end
                    if occursin(varianceDeterminer, obsPar)
                        varianceIndex = findfirst(varianceNames .== obsPar)
                        varianceMap[iCond, iObs] = varianceIndex
                    else
                        varianceMap[iCond, iObs] = numVariance+1
                    end
                end
            end
        end
    end

    numUsedParameters = length(prob.p)
    numVariables = length(prob.u0)
    problemParameterSymbols = parameters(new_sys)
    inputParameterSymbols = Symbol.(names(experimentalConditions)[2:end])
    tmp = [pPS.name for pPS in problemParameterSymbols]
    inputParameterSymbols = inputParameterSymbols[[any(inpPar .== tmp) for inpPar in inputParameterSymbols]]

    optParameterIndices = collect(1:numUsedParameters)[[pPS.name ∉ inputParameterSymbols for pPS in problemParameterSymbols]]
    numOptParameters = length(optParameterIndices)
    parameterNames = string.(problemParameterSymbols)
    optParameterNames = parameterNames[optParameterIndices]

    scaleIndices = collect(1:numScale)
    offsetIndices = maximum(scaleIndices) .+ collect(1:numOffset)
    varianceIndices = maximum(offsetIndices) .+ collect(1:numVariance)
    parameterIndices = maximum(varianceIndices) .+ collect(1:numOptParameters)

    scaleVector = Vector{Float64}(undef, numScale + 1)
    scaleVector[numScale + 1] = 1.0
    offsetVector = Vector{Float64}(undef, numOffset + 2)
    offsetVector[numOffset + 1] = 0.0
    offsetVector[numOffset + 2] = 1.0
    varianceVector = Vector{Float64}(undef, numVariance + 1)
    varianceVector[numVariance + 1] = 1.0
    dynamicParametersVector = Vector{Float64}(undef, numUsedParameters)
    u0Vector = prob.u0
    numAllParameters = numScale + numOffset + numVariance + numOptParameters
    allParameters = Vector{Float64}(undef, numAllParameters)

    modelParameters = ModelParameters(scaleIndices, offsetIndices, varianceIndices, parameterIndices, 
        scaleMap, offsetMap, varianceMap, numScale, numOffset, numVariance, numOptParameters, numUsedParameters, numAllParameters, numVariables, 
        scaleNames, offsetNames, varianceNames, optParameterNames, 
        scaleVector, offsetVector, varianceVector, dynamicParametersVector, u0Vector, allParameters)

    return modelParameters
end


struct DualModelParameters
    dualScaleVector::Vector{ForwardDiff.Dual}
    dualOffsetVector::Vector{ForwardDiff.Dual}
    dualVarianceVector::Vector{ForwardDiff.Dual}
    dualDynamicParametersVector::Vector{ForwardDiff.Dual}
    dualU0Vector::Vector{ForwardDiff.Dual}
end

function DualModelParameters(prob, parameterBounds; 
    scaleDeterminer = "scale", offsetDeterminer = "offset", varianceDeterminer = "sd_")

    numScale = length(findall(occursin.(scaleDeterminer, parameterBounds[!, 1])))
    numOffset = length(findall(occursin.(offsetDeterminer, parameterBounds[!, 1])))
    numVariance = length(findall(occursin.(varianceDeterminer, parameterBounds[!, 1])))
    numUsedParameters = length(prob.p)

    dualScaleVector = Vector{ForwardDiff.Dual}(undef, numScale + 1)
    dualOffsetVector = Vector{ForwardDiff.Dual}(undef, numOffset + 2)
    dualVarianceVector = Vector{ForwardDiff.Dual}(undef, numVariance + 1)
    dualDynamicParameters = Vector{ForwardDiff.Dual}(undef, numUsedParameters)
    dualU0Vector = Vector{ForwardDiff.Dual}(undef, length(prob.u0))

    dualModelParameters = DualModelParameters(dualScaleVector, dualOffsetVector, dualVarianceVector, dualDynamicParameters, dualU0Vector)

    return dualModelParameters
end


struct ParameterSpace
    scaleLowerBound::Vector{Float64}
    scaleUpperBound::Vector{Float64}
    offsetLowerBound::Vector{Float64}
    offsetUpperBound::Vector{Float64}
    varianceLowerBound::Vector{Float64}
    varianceUpperBound::Vector{Float64}
    parameterLowerBound::Vector{Float64}
    parameterUpperBound::Vector{Float64}
    doLogSearch::Vector{Int64}
end

function ParameterSpace(modelParameters, parameterBounds; 
    scaleDeterminer = "scale", offsetDeterminer = "offset", varianceDeterminer = "sd_")
    numScale = modelParameters.numScale
    numOffset = modelParameters.numOffset
    numVariance = modelParameters.numVariance
    numOptParameters = modelParameters.numOptParameters

    ## Setting up lower and upper bounds
    slb = Vector{Float64}(undef, numScale)
    sub = Vector{Float64}(undef, numScale)
    isLogScale = Vector{Bool}(undef, numScale) 
    olb = Vector{Float64}(undef, numOffset)
    oub = Vector{Float64}(undef, numOffset)
    isLogOffset = Vector{Bool}(undef, numOffset) 
    vlb = Vector{Float64}(undef, numVariance)
    vub = Vector{Float64}(undef, numVariance)
    isLogVariance = Vector{Bool}(undef, numVariance) 
    plb = Vector{Float64}(undef, numOptParameters)
    pub = Vector{Float64}(undef, numOptParameters) 
    isLogParameter = Vector{Bool}(undef, numOptParameters) 

    scaleNames = modelParameters.scaleNames
    offsetNames = modelParameters.offsetNames
    varianceNames = modelParameters.varianceNames
    optParameterNames = modelParameters.optParameterNames

    for (i, parId) in enumerate(parameterBounds[!, 1])
        if occursin(scaleDeterminer, parId)
            sIndex = findfirst(lowercase(parId) .== lowercase.(scaleNames))
            if parameterBounds[i, 3] == "log10"
                isLogScale[sIndex] = true
                slb[sIndex] = log10(parameterBounds[i, 4])
                sub[sIndex] = log10(parameterBounds[i, 5])
            else
                isLogScale[sIndex] = false
                slb[sIndex] = parameterBounds[i, 4]
                sub[sIndex] = parameterBounds[i, 5]
            end
        elseif occursin(offsetDeterminer, parId)
            oIndex = findfirst(lowercase(parId) .== lowercase.(offsetNames))
            if parameterBounds[i, 3] == "log10"
                isLogOffset[oIndex] = true
                olb[oIndex] = log10(parameterBounds[i, 4])
                oub[oIndex] = log10(parameterBounds[i, 5])
            else
                isLogOffset[oIndex] = false
                olb[oIndex] = parameterBounds[i, 4]
                oub[oIndex] = parameterBounds[i, 5]
            end
        elseif occursin(varianceDeterminer, parId)
            vIndex = findfirst(lowercase(parId) .== lowercase.(varianceNames))
            if parameterBounds[i, 3] == "log10"
                isLogVariance[vIndex] = true
                vlb[vIndex] = log10(parameterBounds[i, 4])
                vub[vIndex] = log10(parameterBounds[i, 5])
            else
                isLogVariance[vIndex] = false
                vlb[vIndex] = parameterBounds[i, 4]
                vub[vIndex] = parameterBounds[i, 5]
            end
        else
            pIndex = findfirst(lowercase(parId) .== lowercase.(optParameterNames))
            if parameterBounds[i, 3] == "log10"
                isLogParameter[pIndex] = true
                plb[pIndex] = log10(parameterBounds[i, 4])
                pub[pIndex] = log10(parameterBounds[i, 5])
            else
                isLogParameter[pIndex] = false
                plb[pIndex] = parameterBounds[i, 4]
                pub[pIndex] = parameterBounds[i, 5]
            end
        end
    end

    scaleIndices = modelParameters.scaleIndices
    offsetIndices = modelParameters.offsetIndices
    varianceIndices = modelParameters.varianceIndices
    parameterIndices = modelParameters.parameterIndices

    doLogSearchScale = scaleIndices[isLogScale]
    doLogSearchOffset = offsetIndices[isLogOffset]
    doLogSearchVariance = varianceIndices[isLogVariance]
    doLogSearchParameter = parameterIndices[isLogParameter]
    doLogSearch = vcat(doLogSearchScale, doLogSearchOffset, doLogSearchVariance, doLogSearchParameter)

    parameterSpace = ParameterSpace(slb, sub, olb, oub, vlb, vub, plb, pub, doLogSearch)

    lowerBounds = vcat(slb, olb, vlb, plb)
    upperBounds = vcat(sub, oub, vub, pub)
    numAllStartParameters = length(lowerBounds)

    return parameterSpace, numAllStartParameters, lowerBounds, upperBounds
end


struct ModelOutput{T1<:Union{ForwardDiff.Dual, Float64}}
    sols::Vector{ODESolution}
    uForCond::Vector{Array{Float64, 2}}
    dpForCond::Vector{Vector{Array{Float64, 2}}}
    h_barForCondObs::Array{Vector{T1}, 2}
    h_hatForCondObs::Array{Vector{T1}, 2}
    costForCondObs::Array{T1, 2}
    scaleGrad::Vector{Float64}
    offsetGrad::Vector{Float64}
    varianceGrad::Vector{Float64}
    ∂g∂p::Vector{Float64}
    dynParGrad::Vector{Float64}
    gradMatrix::Array{Float64, 2}
    allParametersGrad::Vector{Float64}
end

function ModelOutput(usedType, experimentalData, modelParameters)
    numConditions = experimentalData.numConditions
    numObservables = experimentalData.numObservables
    measurementForCondObs = experimentalData.measurementForCondObs

    sols = Vector{ODESolution}(undef, numConditions)
    uForCond = Vector{Array{Float64, 2}}(undef, numConditions)
    dpForCond = Vector{Vector{Array{Float64, 2}}}(undef, numConditions)

    h_barForCondObs = convert(Array{Vector{usedType}, 2}, copy(measurementForCondObs))
    h_hatForCondObs = convert(Array{Vector{usedType}, 2}, copy(measurementForCondObs))
    costForCondObs = zeros(usedType, numConditions, numObservables)

    scaleGrad = Vector{Float64}(undef, length(modelParameters.scaleVector))
    offsetGrad = Vector{Float64}(undef, length(modelParameters.offsetVector))
    varianceGrad = Vector{Float64}(undef, length(modelParameters.varianceVector))
    ∂g∂p = Vector{Float64}(undef, length(modelParameters.dynamicParametersVector))
    dynamicParameterGrad = Vector{Float64}(undef, length(modelParameters.dynamicParametersVector))
    gradMatrix = Array{Float64, 2}(undef, numConditions, length(modelParameters.allParameters))
    allParametersGrad = Vector{Float64}(undef, modelParameters.numAllParameters)

    modelOutput = ModelOutput{usedType}(sols, uForCond, dpForCond, h_barForCondObs, h_hatForCondObs, costForCondObs, 
            scaleGrad, offsetGrad, varianceGrad, ∂g∂p, dynamicParameterGrad, gradMatrix, allParametersGrad)

    return modelOutput
end


struct FilesAndPaths
    modelName::String
    modelPath::String
    modelFile::String
    methodPath::String
    writePath::String
    writeFile::String
end

