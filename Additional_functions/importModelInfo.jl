

struct ModelData{
        T1 <: Vector{Bool},
        T2 <: Vector{Int64},
        T3 <: Vector{Int64},
        T4 <: Vector{Int64},
        T5 <: Vector{Int64},
        T6 <: Vector{Int64},
        T7 <: Vector{Int64},
        T8 <: Vector{Symbol}
        }
    observableLogTransformation::T1
    optParameterIndices::T2
    inputParameterIndices::T3
    initVariableIndices::T4
    observableVariableIndices::T5
    parameterInObservableIndices::T6
    parameterInU0Indices::T7
    inputParameterSymbols::T8
end

function createModelData(new_sys::ODESystem, 
                         prob::ODEProblem, 
                         observables::DataFrame, 
                         experimentalConditions::DataFrame, 
                         initVariableNames::Vector{String}, 
                         observableVariableNames::Vector{String}, 
                         parameterInU0Names::Vector{String}, 
                         parameterInObservableNames::Vector{String})::ModelData
                         observableLogTransformation::Vector{Bool} = [observables[i, :observableTransformation] == "log10" for i in eachindex(observables[!, :observableTransformation])]

    variableNames::Vector{String} = string.(states(new_sys))
    numUsedParameters::Int64 = length(prob.p)
    tmpvec::Vector{Any} = parameters(new_sys)
    problemParameterSymbols::Vector{Symbol} = [tmp.name for tmp in tmpvec]
    inputParameterSymbols_proto::Vector{Symbol} = Symbol.(names(experimentalConditions))
    inputParameterSymbols::Vector{Symbol} = inputParameterSymbols_proto[[any(inpPar .== problemParameterSymbols) for inpPar in inputParameterSymbols_proto]]
    inputParameterIndices::Vector{Int64} = [findfirst([pPS == inputParameterSymbols[i] for pPS in problemParameterSymbols]) for i in 1:length(inputParameterSymbols)]
    optParameterIndices::Vector{Int64} = collect(1:numUsedParameters)[[pPS ∉ inputParameterSymbols for pPS in problemParameterSymbols]]
    parameterNames::Vector{String} = string.(problemParameterSymbols)

    initVariableNames .*= "(t)"
    initVariableIndices::Vector{Int64} = [findfirst(initVar .== variableNames) for initVar in initVariableNames]

    observableVariableNames .*= "(t)"
    observableVariableIndices::Vector{Int64} = [findfirst(observableVar .== variableNames) for observableVar in observableVariableNames]

    parameterInU0Indices::Vector{Int64} = [findfirst(parIU0N .== parameterNames) for parIU0N in parameterInU0Names]

    parameterInObservableIndices::Vector{Int64} = [findfirst(parION .== parameterNames) for parION in parameterInObservableNames]

    modelData = ModelData(observableLogTransformation, optParameterIndices, inputParameterIndices, 
        initVariableIndices, observableVariableIndices, parameterInObservableIndices, parameterInU0Indices, inputParameterSymbols)

    return modelData
end


struct ExperimentalData{
        T1 <: Array{Vector{Float64}, 2},
        T2 <: Vector{Vector{Float64}},
        T3 <: Array{Vector{Int64}, 2},
        T4 <: Int64,
        T5 <: Array{Int64, 2},
        T6 <: Vector{Vector{Float64}},
        T7 <: Vector{Vector{Int64}},
        T8 <: Vector{Array{Vector{Int64}, 2}}
        }
    measurementForCondObs::T1
    timeStepsForCond::T2
    observedAtIndexForCondObs::T3
    numConditions::T4
    numObservables::T4
    numDataForCondObs::T5
    inputParameterValuesForCond::T6
    observedObservableForCond::T7
    observablesTimeIndexIndicesForCond::T8
end

function createExperimentalData(observables::DataFrame, 
                                experimentalConditions::DataFrame, 
                                measurementData::DataFrame, 
                                modelData::ModelData)::ExperimentalData

    observableNames::Vector{String} = observables[!, 1]
    inputParameterSymbols = modelData.inputParameterSymbols
    observableLogTransformation = modelData.observableLogTransformation

    numConditions::Int64 = length(experimentalConditions[!,1])
    numObservables::Int64 = length(observableNames) 

    ## Setting up for different conditions
    inputParameterValuesForCond::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, numConditions)

    timeStepsForCond::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, numConditions)
    numTimeStepsForCond::Vector{Int64} = Vector{Int64}(undef, numConditions)
    observedObservableForCond::Vector{Vector{Int64}} = Vector{Vector{Int64}}(undef, numConditions)

    observedAtForCond::Vector{Vector{Vector{Float64}}} = Vector{Vector{Vector{Float64}}}(undef, numConditions)
    observedAtIndexForCondObs::Array{Vector{Int64}, 2} = Array{Vector{Int64}, 2}(undef, numConditions, numObservables)
    measurementForCondObs::Array{Vector{Float64}, 2} = Array{Vector{Float64}, 2}(undef, numConditions, numObservables)
    numDataForCondObs::Array{Int64, 2} = Array{Int64, 2}(undef, numConditions, numObservables)
    observablesTimeIndexIndicesForCond::Vector{Array{Vector{Int64}, 2}} = Vector{Array{Vector{Int64}, 2}}(undef, numConditions)

    for (iCond, condId::String) in enumerate(experimentalConditions[!,1])
        inputParameterValuesForCond[iCond] = [experimentalConditions[iCond, inPar] for inPar in inputParameterSymbols]

        relevantMeasurementData::DataFrame = measurementData[measurementData[:,3] .== condId, :]
        observableIDs::Vector{String} = relevantMeasurementData[:,1]
        timeStepsForCond[iCond] = sort(unique(relevantMeasurementData[:,5]))
        numTimeStepsForCond[iCond] = length(timeStepsForCond[iCond])
        isObserved::Vector{Bool} = [observableNames[i] in unique(observableIDs) for i=1:numObservables]
        observedObservableForCond[iCond] = collect(1:numObservables)[isObserved]

        observedAtForCond[iCond] = Vector{Vector{Float64}}(undef, numObservables)
        for (iObs, obsId) = enumerate(observableNames)
            observedAtForCond[iCond][iObs] = relevantMeasurementData[observableIDs .== obsId, 5]
        end
        for iObs = 1:numObservables
            observedAtIndexForCondObs[iCond, iObs] = indexin(observedAtForCond[iCond][iObs], timeStepsForCond[iCond])
        end
        for (iObs, obsId) = enumerate(observableNames)
            if observableLogTransformation[iObs]
                measurementForCondObs[iCond, iObs] = log10.(relevantMeasurementData[observableIDs .== obsId, 4])
            else
                measurementForCondObs[iCond, iObs] = relevantMeasurementData[observableIDs .== obsId, 4]
            end
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


struct ModelParameters{
        T1 <: Vector{Int64},
        T2 <: Vector{Int64},
        T3 <: Vector{Int64},
        T4 <: Vector{Int64},
        T5 <: Array{Int64, 2},
        T6 <: Array{Vector{Int64}, 2},
        T7 <: Int64,
        T8 <: Vector{String},
        T9 <: Vector{String},
        T10 <: Vector{String},
        T11 <: Vector{String},
        T12 <: Vector{Float64},
        T13 <: Vector{Float64},
        T14 <: Vector{Float64},
        T15 <: Vector{Float64},
        T16 <: Vector{Float64},
        T17 <: Vector{Float64}
        }
    scaleIndices::T1
    offsetIndices::T2
    varianceIndices::T3
    parameterIndices::T4
    #
    scaleMap::T5
    offsetMap::T5
    varianceMap::T6
    #
    numScale::T7
    numOffset::T7
    numVariance::T7
    numOptParameters::T7
    numUsedParameters::T7
    numAllParameters::T7
    numVariables::T7
    #
    scaleNames::T8
    offsetNames::T9
    varianceNames::T10
    optParameterNames::T11
    #
    scaleVector::T12
    offsetVector::T13
    varianceVector::T14
    dynamicParametersVector::T15
    u0Vector::T16
    allParameters::T17
end


function createModelParameters(new_sys::ODESystem, 
                               prob::ODEProblem, 
                               parameterBounds::DataFrame, 
                               experimentalConditions::DataFrame, 
                               measurementData::DataFrame, 
                               observables::DataFrame, 
                               experimentalData::ExperimentalData; 
                               scaleDeterminer::String = "scale", 
                               offsetDeterminer::String = "offset", 
                               varianceDeterminer::String = "sd_")::ModelParameters

    numConditions = experimentalData.numConditions
    numObservables = experimentalData.numObservables
    observableNames::Vector{String} = observables[!, 1]

    ## Creating a map from a condition, observable pair to a scale ##
    scaleNames::Vector{String} = parameterBounds[findall(occursin.(scaleDeterminer, parameterBounds[!, 1])), 1]
    numScale::Int64 = length(scaleNames)
    # scale[numScale+1] = 1
    scaleMap::Array{Int64, 2} = Array{Int64, 2}(undef, numConditions, numObservables) # Have as sparse?
    for (iCond, cond) in enumerate(experimentalConditions[!, 1])
        for (iObs, obs) in enumerate(observableNames)
            measurementIndex = findfirst((measurementData[!, 1] .== obs) .& (measurementData[!, 3] .== cond))
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
                        scaleIndex::Int64 = findfirst(scaleNames .== obsPar)
                        scaleMap[iCond, iObs] = scaleIndex
                    else
                        scaleMap[iCond, iObs] = numScale+1
                    end
                end
            end
        end
    end

    ## Creating a map from a condition, observable pair to an offset ##
    offsetNames::Vector{String} = parameterBounds[findall(occursin.(offsetDeterminer, parameterBounds[!, 1])), 1]
    numOffset::Int64 = length(offsetNames)
    # offset[numOffset+1] = 0, offset[numOffset+2] = 1
    offsetMap::Array{Int64, 2} = Array{Int64, 2}(undef, numConditions, numObservables)
    offsetIndex::Int64 = 0
    for (iCond, cond) in enumerate(experimentalConditions[!, 1])
        for (iObs, obs) in enumerate(observableNames)
            measurementIndex = findfirst((measurementData[!, 1] .== obs) .& (measurementData[!, 3] .== cond))
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
    varianceNames::Vector{String} = parameterBounds[findall(occursin.(varianceDeterminer, parameterBounds[!, 1])), 1]
    numVariance::Int64 = length(varianceNames)
    allKnownVariances::Vector{Float64} = [1.0]
    numKnownVariance::Int64 = 1
    if numVariance == 0 && typeof(measurementData[:, :noiseParameters]) <: Vector{Float64}
        allKnownVariances = vcat(unique(measurementData[:, :noiseParameters]), allKnownVariances)
        numKnownVariance = length(allKnownVariances)
    end
    varianceMap::Array{Vector{Int64}, 2} = Array{Vector{Int64}, 2}(undef, numConditions, numObservables)
    for (iCond, cond) in enumerate(experimentalConditions[!, 1])
        for (iObs, obs) in enumerate(observableNames)
            measurementIndex = findfirst((measurementData[!, 1] .== obs) .& (measurementData[!, 3] .== cond))
            if measurementIndex === nothing
                varianceMap[iCond, iObs] = [numVariance + numKnownVariance]
            else
                obsPar = measurementData[measurementIndex, :noiseParameters]
                if obsPar === missing
                    varianceMap[iCond, iObs] = [numVariance + numKnownVariance]
                else
                    if occursin(';', obsPar)
                        obsPars = split(obsPar, ';')
                        tmpvec = Vector{Int64}(undef, length(obsPars))
                        for (i, obsPar) in enumerate(obsPars)
                            if typeof(obsPar) <: Float64
                                knownVarianceIndex = findfirst(allKnownVariances .== obsPar)::Int64
                                tmpvec[i] = numVariance + knownVarianceIndex
                            else
                                if occursin(varianceDeterminer, obsPar)
                                    varianceIndex = findfirst(varianceNames .== obsPar)::Int64
                                    tmpvec[i] = varianceIndex
                                else
                                    tmpvec[i] = numVariance + numKnownVariance
                                end
                            end
                        end
                        varianceMap[iCond, iObs] = copy(tmpvec)
                    else
                        if typeof(obsPar) <: Float64
                            knownVarianceIndex::Int64 = findfirst(allKnownVariances .== obsPar)
                            varianceMap[iCond, iObs] = [numVariance + knownVarianceIndex]
                        else
                            if occursin(varianceDeterminer, obsPar)
                                varianceIndex::Int64 = findfirst(varianceNames .== obsPar)
                                varianceMap[iCond, iObs] = [varianceIndex]
                            else
                                varianceMap[iCond, iObs] = [numVariance + numKnownVariance]
                            end
                        end
                    end
                end
            end
        end
    end

    numUsedParameters::Int64 = length(prob.p)
    numVariables::Int64 = length(prob.u0)
    tmpvec = parameters(new_sys)
    problemParameterSymbols::Vector{Symbol} = [tmp.name for tmp in tmpvec]
    inputParameterSymbols_proto::Vector{Symbol} = Symbol.(names(experimentalConditions))

    inputParameterSymbols::Vector{Symbol} = inputParameterSymbols_proto[[any(inpPar .== problemParameterSymbols) for inpPar in inputParameterSymbols_proto]]

    optParameterIndices::Vector{Int64} = collect(1:numUsedParameters)[[pPS ∉ inputParameterSymbols for pPS in problemParameterSymbols]]
    numOptParameters::Int64 = length(optParameterIndices)
    parameterNames::Vector{String} = string.(problemParameterSymbols)
    optParameterNames::Vector{String} = parameterNames[optParameterIndices]

    scaleIndices::Vector{Int64} = collect(1:numScale)
    offsetIndices::Vector{Int64} = numScale .+ collect(1:numOffset)
    varianceIndices::Vector{Int64} = numScale + numOffset .+ collect(1:numVariance)
    parameterIndices::Vector{Int64} = numScale + numOffset + numVariance .+ collect(1:numOptParameters)

    scaleVector::Vector{Float64} = Vector{Float64}(undef, numScale + 1)
    scaleVector[numScale + 1] = 1.0
    offsetVector::Vector{Float64} = Vector{Float64}(undef, numOffset + 2)
    offsetVector[numOffset + 1] = 0.0
    offsetVector[numOffset + 2] = 1.0
    varianceVector::Vector{Float64} = Vector{Float64}(undef, numVariance + numKnownVariance)
    varianceVector[numVariance + 1:end] = allKnownVariances
    dynamicParametersVector::Vector{Float64} = Vector{Float64}(undef, numUsedParameters)
    u0Vector::Vector{Float64} = prob.u0
    numAllParameters::Int64 = numScale + numOffset + numVariance + numOptParameters
    allParameters::Vector{Float64} = Vector{Float64}(undef, numAllParameters)

    modelParameters = ModelParameters(scaleIndices, offsetIndices, varianceIndices, parameterIndices, 
        scaleMap, offsetMap, varianceMap, numScale, numOffset, numVariance, numOptParameters, numUsedParameters, numAllParameters, numVariables, 
        scaleNames, offsetNames, varianceNames, optParameterNames, 
        scaleVector, offsetVector, varianceVector, dynamicParametersVector, u0Vector, allParameters)

    return modelParameters
end


struct DualModelParameters{
        T1 <: Vector{ForwardDiff.Dual},
        T2 <: Vector{ForwardDiff.Dual},
        T3 <: Vector{ForwardDiff.Dual},
        T4 <: Vector{ForwardDiff.Dual},
        T5 <: Vector{ForwardDiff.Dual}
        }
    dualScaleVector::T1
    dualOffsetVector::T2
    dualVarianceVector::T3
    dualDynamicParametersVector::T4
    dualU0Vector::T5
end

function createDualModelParameters(modelParameters::ModelParameters)::DualModelParameters

    dualScaleVector::Vector{ForwardDiff.Dual} = Vector{ForwardDiff.Dual}(undef, length(modelParameters.scaleVector))
    dualOffsetVector::Vector{ForwardDiff.Dual} = Vector{ForwardDiff.Dual}(undef, length(modelParameters.offsetVector))
    dualVarianceVector::Vector{ForwardDiff.Dual} = Vector{ForwardDiff.Dual}(undef, length(modelParameters.varianceVector))
    dualDynamicParameters::Vector{ForwardDiff.Dual} = Vector{ForwardDiff.Dual}(undef, length(modelParameters.dynamicParametersVector))
    dualU0Vector::Vector{ForwardDiff.Dual} = Vector{ForwardDiff.Dual}(undef, length(modelParameters.u0Vector))

    dualModelParameters = DualModelParameters(dualScaleVector, dualOffsetVector, dualVarianceVector, dualDynamicParameters, dualU0Vector)

    return dualModelParameters
end


struct ParameterSpace{
        T1 <: Vector{Float64},
        T2 <: Vector{Float64},
        T3 <: Vector{Float64},
        T4 <: Vector{Float64},
        T5 <: Vector{Int64}
        }
    scaleLowerBound::T1
    scaleUpperBound::T1
    offsetLowerBound::T2
    offsetUpperBound::T2
    varianceLowerBound::T3
    varianceUpperBound::T3
    parameterLowerBound::T4
    parameterUpperBound::T4
    doLogSearch::T5
end

function createParameterSpace(modelParameters::ModelParameters, 
                              parameterBounds::DataFrame; 
                              scaleDeterminer::String = "scale", 
                              offsetDeterminer::String = "offset", 
                              varianceDeterminer::String = "sd_")::Tuple{ParameterSpace, Int64, Vector{Float64}, Vector{Float64}}
    
                              numScale = modelParameters.numScale
    numOffset = modelParameters.numOffset
    numVariance = modelParameters.numVariance
    numOptParameters = modelParameters.numOptParameters

    ## Setting up lower and upper bounds
    slb::Vector{Float64} = Vector{Float64}(undef, numScale)
    sub::Vector{Float64} = Vector{Float64}(undef, numScale)
    isLogScale::Vector{Bool} = Vector{Bool}(undef, numScale) 
    olb::Vector{Float64} = Vector{Float64}(undef, numOffset)
    oub::Vector{Float64} = Vector{Float64}(undef, numOffset)
    isLogOffset::Vector{Bool} = Vector{Bool}(undef, numOffset) 
    vlb::Vector{Float64} = Vector{Float64}(undef, numVariance)
    vub::Vector{Float64} = Vector{Float64}(undef, numVariance)
    isLogVariance::Vector{Bool} = Vector{Bool}(undef, numVariance) 
    plb::Vector{Float64} = Vector{Float64}(undef, numOptParameters)
    pub::Vector{Float64} = Vector{Float64}(undef, numOptParameters) 
    isLogParameter::Vector{Bool} = Vector{Bool}(undef, numOptParameters) 

    scaleNames = modelParameters.scaleNames
    offsetNames = modelParameters.offsetNames
    varianceNames = modelParameters.varianceNames
    optParameterNames = modelParameters.optParameterNames

    for (i, parId) in enumerate(parameterBounds[!, 1])
        if occursin(scaleDeterminer, parId)
            sIndex::Int64 = findfirst(lowercase(parId) .== lowercase.(scaleNames))
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
            oIndex::Int64 = findfirst(lowercase(parId) .== lowercase.(offsetNames))
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
            vIndex::Int64 = findfirst(lowercase(parId) .== lowercase.(varianceNames))
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
            pIndex::Int64 = findfirst(lowercase(parId) .== lowercase.(optParameterNames))
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

    doLogSearchScale::Vector{Int64} = scaleIndices[isLogScale]
    doLogSearchOffset::Vector{Int64} = offsetIndices[isLogOffset]
    doLogSearchVariance::Vector{Int64} = varianceIndices[isLogVariance]
    doLogSearchParameter::Vector{Int64} = parameterIndices[isLogParameter]
    doLogSearch::Vector{Int64} = vcat(doLogSearchScale, doLogSearchOffset, doLogSearchVariance, doLogSearchParameter)

    parameterSpace = ParameterSpace(slb, sub, olb, oub, vlb, vub, plb, pub, doLogSearch)

    lowerBounds::Vector{Float64} = vcat(slb, olb, vlb, plb)
    upperBounds::Vector{Float64} = vcat(sub, oub, vub, pub)
    numAllStartParameters::Int64 = length(lowerBounds)

    return parameterSpace, numAllStartParameters, lowerBounds, upperBounds
end


struct ModelOutput{
        T0 <: Union{ODESolution, OrdinaryDiffEq.ODECompositeSolution},
        T1 <: Vector{T0},
        T2 <: Vector{Array{Float64, 2}},
        T3 <: Vector{Vector{Array{Float64, 2}}},
        T4 <: Union{ForwardDiff.Dual, Float64},
        T5 <: Array{Vector{T4}, 2},
        T6 <: Array{T4, 2},
        T7 <: Vector{Float64},
        T8 <: Vector{Float64},
        T9 <: Vector{Float64},
        T10 <: Vector{Float64},
        T11 <: Array{Float64, 2},
        T12 <: Vector{Float64}
        }
    sols::T1
    uForCond::T2
    dpForCond::T3
    h_barForCondObs::T5
    h_hatForCondObs::T5
    costForCondObs::T6
    scaleGrad::T7
    offsetGrad::T8
    varianceGrad::T9
    ∂g∂p::T10
    dynParGrad::T10
    gradMatrix::T11
    allParametersGrad::T12
end

function createModelOutput(usedODESolType, 
                           usedType, 
                           experimentalData::ExperimentalData, 
                           modelParameters::ModelParameters)
                           
    numConditions = experimentalData.numConditions
    numObservables = experimentalData.numObservables
    measurementForCondObs = experimentalData.measurementForCondObs

    sols = Vector{usedODESolType}(undef, numConditions)
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

    modelOutput = ModelOutput(sols, uForCond, dpForCond, h_barForCondObs, h_hatForCondObs, costForCondObs, 
            scaleGrad, offsetGrad, varianceGrad, ∂g∂p, dynamicParameterGrad, gradMatrix, allParametersGrad)

    return modelOutput
end


struct FilesAndPaths{
        T1 <: String
        }
    modelName::T1
    modelPath::T1
    modelFile::T1
    methodPath::T1
    writePath::T1
    writeFile::T1
end