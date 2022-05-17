


function calcScaledObservable_forwAD_proto(scaleVector::Vector{T1}, offsetVector::Vector{T1}, h_barFCO::Array{Vector{T2}, 2}, h_hatFCO::Array{Vector{T3}, 2}, modelParameters::ModelParameters, 
        modelData::ModelData, experimentalData::ExperimentalData, iCond::Int64)::Nothing where {T1 <: Union{ForwardDiff.Dual, Float64}, T2 <: Union{ForwardDiff.Dual, Float64}, T3 <: Union{ForwardDiff.Dual, Float64}}

    scaleMap = modelParameters.scaleMap
    offsetMap = modelParameters.offsetMap

    observedOFC = experimentalData.observedObservableForCond

    observableLT = modelData.observableLogTransformation

    for iObs = observedOFC[iCond]
        if observableLT[iObs]
            h_hatFCO[iCond, iObs] = log10.(scaleVector[scaleMap[iCond, iObs]] * h_barFCO[iCond, iObs] .+ offsetVector[offsetMap[iCond, iObs]])
        else
            h_hatFCO[iCond, iObs] = scaleVector[scaleMap[iCond, iObs]] * h_barFCO[iCond, iObs] .+ offsetVector[offsetMap[iCond, iObs]]
        end
    end

    nothing
end

function calcCost_forwAD_proto(varianceVector::Vector{T1}, h_hatFCO::Array{Vector{T2}, 2}, costFCO::Array{T3, 2}, modelParameters::ModelParameters, 
        experimentalData::ExperimentalData, iCond::Int64)::T3 where {T1 <: Union{ForwardDiff.Dual, Float64}, T2 <: Union{ForwardDiff.Dual, Float64}, T3 <: Union{ForwardDiff.Dual, Float64}}

    varianceMap = modelParameters.varianceMap

    numDataFCO = experimentalData.numDataForCondObs

    measurementFCO = experimentalData.measurementForCondObs
    observedOFC = experimentalData.observedObservableForCond

    for iObs in observedOFC[iCond]
        if length(varianceMap[iCond, iObs]) > 1
            variance = sum(sqrt.(varianceVector[varianceMap[iCond, iObs]]))^2
        else
            variance = varianceVector[varianceMap[iCond, iObs]][1]
        end
        costFCO[iCond, iObs] = log(2*pi*variance) * numDataFCO[iCond, iObs] + (dot(measurementFCO[iCond, iObs], measurementFCO[iCond, iObs]) - 
            2*dot(measurementFCO[iCond, iObs], h_hatFCO[iCond, iObs]) + 
            dot(h_hatFCO[iCond, iObs], h_hatFCO[iCond, iObs])) / (variance)
    end

    return sum(costFCO[iCond, observedOFC[iCond]])
end

function allConditionsCost_forwAD_float_proto(modelParameters::ModelParameters, experimentalData::ExperimentalData, modelData::ModelData, 
        solveODESystem::Function, calcUnscaledObservable::Function, calcScaledObservable::Function, calcCost::Function)::Float64

    dynParVector = modelParameters.dynamicParametersVector
    inputPVFC = experimentalData.inputParameterValuesForCond
    inputPI = modelData.inputParameterIndices

    cost = 0.0
    for (iCond, inputParameterValues) in enumerate(inputPVFC)
        dynParVector[inputPI] = inputParameterValues
        solveODESystem(iCond)
        calcUnscaledObservable(iCond)
        calcScaledObservable(iCond)
        cost += calcCost(iCond)
    end

    return cost
end

function allConditionsCost_forwAD_dual_proto(modelParameters::ModelParameters, dualModelParameters::DualModelParameters, experimentalData::ExperimentalData, 
        modelData::ModelData, solveODESystem_dual::Function, calcUnscaledObservable_dual::Function, 
        calcScaledObservable_dual_float::Function, calcCost_dual_float::Function, dynPar::Vector{T1})::T1 where T1 <: ForwardDiff.Dual

    dynParVector = modelParameters.dynamicParametersVector
    dualDynParVector = dualModelParameters.dualDynamicParametersVector
    dualDynParVector .= convert.(eltype(dynPar), dynParVector)::Vector{eltype(dynPar)}
    dualDynParVector[modelData.optParameterIndices] .= dynPar

    u0Vector = modelParameters.u0Vector
    dualU0Vector = convert(Vector{eltype(dynPar)}, u0Vector)::Vector{eltype(dynPar)}

    inputPVFC = experimentalData.inputParameterValuesForCond
    inputPI = modelData.inputParameterIndices

    cost = 0.0
    for (iCond, inputParameterValues) in enumerate(inputPVFC)
        dualDynParVector[inputPI] = convert.(eltype(dynPar), inputParameterValues)::Vector{eltype(dynPar)}
        solveODESystem_dual(dualU0Vector, iCond)
        calcUnscaledObservable_dual(iCond)
        calcScaledObservable_dual_float(iCond)
        cost += calcCost_dual_float(iCond)
    end

    return cost
end

function allConditionsCost_noODE_forwAD_dual_proto(experimentalData::ExperimentalData, updateNonDynamicDualParameterVectors::Function, 
            calcScaledObservable_dual::Function, calcCost_dual::Function, nonDymPar::Vector{T1})::T1 where T1 <: ForwardDiff.Dual

    updateNonDynamicDualParameterVectors(nonDymPar)
    
    cost = 0.0
    for iCond in 1:experimentalData.numConditions
        calcScaledObservable_dual(iCond)
        cost += calcCost_dual(iCond)
    end

    return cost
end


function f_cost_forwAD_proto(updateAllParameterVectors::Function, allConditionsCost_float::Function, parameterSpace::ParameterSpace, modelParameters::ModelParameters, p...)::Float64
    doLogSearch = parameterSpace.doLogSearch
    allPar = modelParameters.allParameters
    allPar .= p
    view(allPar, doLogSearch) .= exp10.(view(allPar, doLogSearch))

    updateAllParameterVectors()

    cost = allConditionsCost_float()

    return cost
end

function f_grad_forwAD_proto(grad, results::Vector{<:DiffResults.MutableDiffResult}, updateAllParameterVectors::Function, allConditionsCost_noODE_dual::Function, allConditionsCost_dual::Function, cfg::GradientConfig, parameterSpace::ParameterSpace, 
        modelParameters::ModelParameters, modelOutput::ModelOutput, p...)

    doLogSearch = parameterSpace.doLogSearch
    allPar = modelParameters.allParameters
    allPar .= p
    view(allPar, doLogSearch) .= exp10.(view(allPar, doLogSearch))

    updateAllParameterVectors()

    scaleIndices = modelParameters.scaleIndices
    offsetIndices = modelParameters.offsetIndices
    varianceIndices = modelParameters.varianceIndices
    nonDynPar = view(allPar, vcat(scaleIndices, offsetIndices, varianceIndices))

    ForwardDiff.gradient!(results[1], allConditionsCost_noODE_dual, nonDynPar)

    dynParIndices = modelParameters.parameterIndices
    dynPar = view(allPar, dynParIndices)

    ForwardDiff.gradient!(results[2], allConditionsCost_dual, dynPar, cfg)

    allParametersGrad = modelOutput.allParametersGrad
    allParametersGrad[:] = vcat(DiffResults.gradient(results[1]), DiffResults.gradient(results[2]))
    view(allParametersGrad, doLogSearch) .*= view(allPar, doLogSearch) * log(10)

    grad[:] = allParametersGrad

    nothing
end

