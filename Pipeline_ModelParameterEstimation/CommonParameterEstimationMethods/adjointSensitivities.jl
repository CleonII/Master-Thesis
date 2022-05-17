

function g_scaledObservationFunctions_AdjSens_proto(type::Type, h_bar::Vector{Vector{T1}}, scaleVector::Vector{T2}, offsetVector::Vector{T3}, 
        iCond::Int64, modelParameters::ModelParameters, experimentalData::ExperimentalData, modelData::ModelData)::Vector{Vector{type}} where 
        {T1 <: Union{ForwardDiff.Dual, Float64}, T2 <: Union{ForwardDiff.Dual, Float64}, T3 <: Union{ForwardDiff.Dual, Float64}}
    
        scaleMap = modelParameters.scaleMap
    offsetMap = modelParameters.offsetMap
    observedOFC = experimentalData.observedObservableForCond

    observableLT = modelData.observableLogTransformation
    h_hat = Vector{Vector{type}}(undef, experimentalData.numObservables)
    for iObs = observedOFC[iCond]
        if observableLT[iObs]
            h_hat[iObs] = log10.(scaleVector[scaleMap[iCond, iObs]] * h_bar[iObs] .+ offsetVector[offsetMap[iCond, iObs]])
        else
            h_hat[iObs] = scaleVector[scaleMap[iCond, iObs]] * h_bar[iObs] .+ offsetVector[offsetMap[iCond, iObs]]
        end
    end

    return h_hat
end

function g_cost_AdjSens_proto(type::Type, h_hat::Vector{Vector{T1}}, varianceVector::Vector{T2}, i::Int64, iCond::Int64, modelParameters::ModelParameters, 
        experimentalData::ExperimentalData)::type where {T1 <: Union{ForwardDiff.Dual, Float64}, T2 <: Union{ForwardDiff.Dual, Float64}}
    varianceMap = modelParameters.varianceMap
    measurementFCO = experimentalData.measurementForCondObs
    observablesTIIFC = experimentalData.observablesTimeIndexIndicesForCond[iCond]
    observedOFC = experimentalData.observedObservableForCond

    costFCO = Vector{type}(undef, experimentalData.numObservables)
    for iObs in observedOFC[iCond]
        if length(varianceMap[iCond, iObs]) > 1
            variance = sum(sqrt.(varianceVector[varianceMap[iCond, iObs]]))^2
        else
            variance = varianceVector[varianceMap[iCond, iObs]][1]
        end
        costFCO[iObs] = log(2*pi*variance) * length(h_hat[iObs]) + 
                (dot(measurementFCO[iCond, iObs][observablesTIIFC[iObs, i]], measurementFCO[iCond, iObs][observablesTIIFC[iObs, i]]) - 
                2*dot(measurementFCO[iCond, iObs][observablesTIIFC[iObs, i]], h_hat[iObs]) + 
                dot(h_hat[iObs], h_hat[iObs])) / (variance)
    end

    return sum(costFCO[observedOFC[iCond]])
end



function g_AdjSens_proto(u::Vector{T1}, dynPar::Vector{T2}, scale::Vector{T3}, offset::Vector{T4}, variance::Vector{T5}, i::Int64, iCond::Int64, 
        g_unscaledObservables::Function, g_scaledObservationFunctions::Function, g_cost::Function; 
        type = get_type([u, dynPar, scale, offset, variance]))::type where {T1 <: Union{ForwardDiff.Dual, Float64}, T2 <: Union{ForwardDiff.Dual, Float64}, 
        T3 <: Union{ForwardDiff.Dual, Float64}, T4 <: Union{ForwardDiff.Dual, Float64}, T5 <: Union{ForwardDiff.Dual, Float64}, T6 <: Union{ForwardDiff.Dual, Float64}}

    h_bar = g_unscaledObservables(type, u, dynPar, i, iCond)

    h_hat = g_scaledObservationFunctions(type, h_bar, scale, offset, iCond)

    cost = g_cost(type, h_hat, variance, i, iCond)

    return cost
end

function G_AdjSens_proto(solveODESystem::Function, g::Function, iCond::Int64, modelParameters::ModelParameters, experimentalData::ExperimentalData, 
        modelOutput::ModelOutput)::Float64
    
    solveODESystem(iCond)
    
    dynPar = modelParameters.dynamicParametersVector
    scale = modelParameters.scaleVector
    offset = modelParameters.offsetVector
    variance = modelParameters.varianceVector
    costAt = (u, i) -> g(u, dynPar, scale, offset, variance, i, iCond, type = Float64)
    sol = modelOutput.sols[iCond]

    cost = 0.0
    ts = experimentalData.timeStepsForCond[iCond]
    for (i, t) in enumerate(ts)
        cost += costAt(sol(t), i)
    end

    return cost
end

function G_specifiedDynPar_AdjSens_proto(specifiedDynPar::Vector{T1}, specifiedDynParIndices::Vector{Int64}, g::Function, solveODESystem::Function, 
        iCond::Int64, modelParameters::ModelParameters, dualModelParameters::DualModelParameters, experimentalData::ExperimentalData, 
        modelOutput::ModelOutput)::T1 where T1 <: ForwardDiff.Dual
    
    dynParVector = modelParameters.dynamicParametersVector
    dualDynParVector = dualModelParameters.dualDynamicParametersVector
    dualDynParVector .= convert.(eltype(specifiedDynPar), dynParVector)
    dualDynParVector[specifiedDynParIndices] .= specifiedDynPar

    u0Vector = modelParameters.u0Vector
    dualU0Vector = convert(Vector{eltype(specifiedDynPar)}, u0Vector)

    solveODESystem(dualDynParVector, dualU0Vector, iCond)

    scale = modelParameters.scaleVector
    offset = modelParameters.offsetVector
    variance = modelParameters.varianceVector

    costAt = (u, i) -> g(u, dualDynParVector, scale, offset, variance, i, iCond, type = eltype(specifiedDynPar))
    sol = modelOutput.sols[iCond]

    cost = 0.0
    ts = experimentalData.timeStepsForCond[iCond]
    for (i, t) in enumerate(ts)
        cost += costAt(sol(t), i)
    end

    return cost
end



function allConditionsCost_AdjSens_proto(parameterSpace::ParameterSpace, modelParameters::ModelParameters, experimentalData::ExperimentalData, 
        modelData::ModelData, updateAllParameterVectors::Function, G::Function, p...)

    doLogSearch = parameterSpace.doLogSearch
    allParameters = modelParameters.allParameters
    allParameters .= p 
    view(allParameters, doLogSearch) .= exp10.(view(allParameters, doLogSearch))

    updateAllParameterVectors()

    dynParVector = modelParameters.dynamicParametersVector
    inputPVFC = experimentalData.inputParameterValuesForCond
    inputPI = modelData.inputParameterIndices

    cost = 0.0
    for (iCond, inputParameterValues) in enumerate(inputPVFC)
        dynParVector[inputPI] = inputParameterValues
        cost += G(iCond)
    end

    return cost
end

function dg_AdjSens_proto!(out, u, p, t, i, gu) 

    guAt = (u) -> gu(u, i)
    out .= ForwardDiff.gradient(guAt, u)

    nothing
end

function calcCostGrad_AdjSens_proto(g::Function, dg!::Function, G_specifiedDynPar::Function, iCond::Int64, modelParameters::ModelParameters, 
        modelData::ModelData, experimentalData::ExperimentalData, modelOutput::ModelOutput, senseAlg, solver)::Nothing

    timeSteps = experimentalData.timeStepsForCond[iCond]
    dynParGrad = modelOutput.dynParGrad

    sol = modelOutput.sols[iCond]

    ~, dynParGrad[:] = adjoint_sensitivities(sol, solver, dg!, timeSteps, 
            sensealg = senseAlg, reltol = 1e-9, abstol = 1e-9)

    # when a parameter is included in the observation function the gradient is incorrect, has to correct with adding dgdp 

    # gradient for scales, offsets, variances and dgdp 
    dynPar = modelParameters.dynamicParametersVector
    ∂g∂p = modelOutput.∂g∂p
    ∂g∂p .= 0.0
    scale = modelParameters.scaleVector
    scaleGrad = modelOutput.scaleGrad
    scaleGrad .= 0.0
    offset = modelParameters.offsetVector
    offsetGrad = modelOutput.offsetGrad
    offsetGrad .= 0.0
    variance = modelParameters.varianceVector
    varianceGrad = modelOutput.varianceGrad
    varianceGrad .= 0.0

    gdynPar = (u, dynPar, i) -> g(u, dynPar, scale, offset, variance, i, iCond)
    gscale = (u, scale, i) -> g(u, dynPar, scale, offset, variance, i, iCond)
    goffset = (u, offset, i) -> g(u, dynPar, scale, offset, variance, i, iCond)
    gvariance = (u, variance, i) -> g(u, dynPar, scale, offset, variance, i, iCond)
    for (i, t) in enumerate(timeSteps)
        gdynParAt = (dynPar) -> gdynPar(sol(t), dynPar, i)
        ∂g∂p += ForwardDiff.gradient(gdynParAt, dynPar)

        gscaleAt = (scale) -> gscale(sol(t), scale, i)
        scaleGrad += ForwardDiff.gradient(gscaleAt, scale)

        goffsetAt = (offset) -> goffset(sol(t), offset, i)
        offsetGrad += ForwardDiff.gradient(goffsetAt, offset)

        gvarianceAt = (variance) -> gvariance(sol(t), variance, i)
        varianceGrad += ForwardDiff.gradient(gvarianceAt, variance)
    end

    # For parameters that influence the starting concentrations 

    specifiedDynParIndices = modelData.parameterInU0Indices
    specifiedDynPar = dynPar[specifiedDynParIndices]

    G_specDynPar = (specifiedDynPar) -> G_specifiedDynPar(specifiedDynPar, iCond)
    correctSpecDynParGrad = ForwardDiff.gradient(G_specDynPar, specifiedDynPar)

    dynParGrad *= -1 # adjoint sensitivities seems to return the negative gradient
    dynParGrad += ∂g∂p
    dynParGrad[specifiedDynParIndices] = correctSpecDynParGrad

    # Insert into common gradient
    gradMatrix = modelOutput.gradMatrix

    scaleIndices = modelParameters.scaleIndices
    gradMatrix[iCond, scaleIndices] = view(scaleGrad, 1:modelParameters.numScale)
    offsetIndices = modelParameters.offsetIndices
    gradMatrix[iCond, offsetIndices] = view(offsetGrad, 1:modelParameters.numOffset)
    varianceIndices = modelParameters.varianceIndices
    gradMatrix[iCond, varianceIndices] = view(varianceGrad, 1:modelParameters.numVariance)

    parameterIndices = modelParameters.parameterIndices
    optParameterIndices = modelData.optParameterIndices
    gradMatrix[iCond, parameterIndices] = view(dynParGrad, optParameterIndices)

    nothing
end


function allConditionsCostGrad_AdjSens_proto(parameterSpace::ParameterSpace, modelParameters::ModelParameters, modelData::ModelData, modelOutput::ModelOutput, 
        experimentalData::ExperimentalData, updateAllParameterVectors::Function, calcCostGrad::Function, g::Function, grad, p...)

    doLogSearch = parameterSpace.doLogSearch
    allParameters = modelParameters.allParameters
    allParameters .= p 
    view(allParameters, doLogSearch) .= exp10.(view(allParameters, doLogSearch))

    updateAllParameterVectors()

    dynPar = modelParameters.dynamicParametersVector
    scale = modelParameters.scaleVector
    offset = modelParameters.offsetVector
    variance = modelParameters.varianceVector
    inputPVFC = experimentalData.inputParameterValuesForCond
    inputPI = modelData.inputParameterIndices
    gradMatrix = modelOutput.gradMatrix

    for (iCond, inputParameterValues) in enumerate(inputPVFC)
        dynPar[inputPI] = inputParameterValues
        gu = (u, i) -> g(u, dynPar, scale, offset, variance, i, iCond)
        dg! = (out, u, p, t, i) -> dg_AdjSens_proto!(out, u, p, t, i, gu)
        calcCostGrad(dg!, iCond)
    end

    replace!(gradMatrix, NaN=>0.0)

    grad[:] = vec(sum(gradMatrix, dims=1))
    view(grad, doLogSearch) .*= view(allParameters, doLogSearch) * log(10)

    nothing
end


