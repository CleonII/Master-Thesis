

function solveODESystem_forwGrad_proto_model_Boehm_JProteomeRes2014(prob::ODEProblem, solver, dynParVector::Vector{T1}, u0Vector::Vector{T2}, modelData::ModelData, 
    modelOutput::ModelOutput, iCond::Int64)::Bool where {T1 <: ForwardDiff.Dual, T2 <: ForwardDiff.Dual}

    initVariable = modelData.initVariableIndices
    parameterInU0Indices = modelData.parameterInU0Indices
    u0Vector[initVariable[1]] = 207.6 * dynParVector[parameterInU0Indices[1]]
    u0Vector[initVariable[2]] = 207.6 - 207.6 * dynParVector[parameterInU0Indices[1]]

    _prob = remake(prob, u0 = u0Vector, p = dynParVector)

    try
        modelOutput.sols[iCond] = OrdinaryDiffEq.solve(_prob, solver, reltol=1e-9, abstol=1e-9)
        if modelOutput.sols[iCond].t[end] !== prob.tspan[2]
                return true
        else
                return false
        end
    catch err
        return true
    end

end

function calcUnscaledObservable_forwGrad_proto_model_Boehm_JProteomeRes2014(dynParVector::Vector{T1}, sols::Vector{T2}, h_barFCO::Array{Vector{T1}, 2}, modelData::ModelData, 
    experimentalData::ExperimentalData, iCond::Int64)::Nothing where {T1 <: ForwardDiff.Dual, T2 <: Union{ODESolution, OrdinaryDiffEq.ODECompositeSolution}}
    
    timeSteps = experimentalData.timeStepsForCond[iCond]
    sol = sols[iCond](timeSteps)

    observedAIFCO = experimentalData.observedAtIndexForCondObs
    optParIOI = modelData.parameterInObservableIndices

    h_barFCO[iCond, 1] = ( (100 * sol[Symbol("pApB(t)")] + 200 * sol[Symbol("pApA(t)")] * dynParVector[optParIOI[1]] ) ./ 
            ( sol[Symbol("pApB(t)")] + sol[Symbol("STAT5A(t)")] * dynParVector[optParIOI[1]] + 2 * sol[Symbol("pApA(t)")] * dynParVector[optParIOI[1]]) )[observedAIFCO[iCond, 1]]
    h_barFCO[iCond, 2] = ( -1 * ( 100 * sol[Symbol("pApB(t)")] - 200 * sol[Symbol("pBpB(t)")] * (dynParVector[optParIOI[1]] - 1) ) ./ 
            ( sol[Symbol("STAT5B(t)")] * (dynParVector[optParIOI[1]] - 1) - sol[Symbol("pApB(t)")] + 2 * sol[Symbol("pBpB(t)")] * (dynParVector[optParIOI[1]] - 1) ) )[observedAIFCO[iCond, 2]]
    h_barFCO[iCond, 3] = ( ( 100 * sol[Symbol("pApB(t)")] + 100 * sol[Symbol("STAT5A(t)")] * dynParVector[optParIOI[1]] + 200 * sol[Symbol("pApA(t)")] * dynParVector[optParIOI[1]] ) ./ 
            (2 * sol[Symbol("pApB(t)")] + sol[Symbol("STAT5A(t)")] * dynParVector[optParIOI[1]] + 2 * sol[Symbol("pApA(t)")] * dynParVector[optParIOI[1]] - 
            sol[Symbol("STAT5B(t)")] * (dynParVector[optParIOI[1]] - 1) - 2 * sol[Symbol("pBpB(t)")] * (dynParVector[optParIOI[1]] - 1)) )[observedAIFCO[iCond, 3]]

    nothing
end



function forwardGradient_model_Boehm_JProteomeRes2014(modelFunction::Function, iStartPar::Int64, n_it::Int64, b2::Float64, stepRange::Vector{<:AbstractFloat}, solver, 
    timeEnd::Float64, experimentalConditions::DataFrame, measurementData::DataFrame, observables::DataFrame, parameterBounds::DataFrame)::Tuple{Function, Adam, Vector{Int64}}

    sys, initialSpeciesValues, trueParameterValues = modelFunction()
    new_sys = ode_order_lowering(sys)
    u0 = initialSpeciesValues
    pars = trueParameterValues 
    tspan = (0.0, timeEnd)
    prob = ODEProblem(new_sys,u0,tspan,pars,jac=true)

    initVariableNames = ["STAT5A", "STAT5B"]

    observableVariableNames = ["pApB", "pApA", "STAT5A", "pBpB", "STAT5B"]

    parameterInU0Names = ["ratio"]

    parameterInObservableNames = ["specC17"]

    # Initialize structs

    modelData = createModelData(new_sys, prob, observables, experimentalConditions, 
            initVariableNames, observableVariableNames, parameterInU0Names, parameterInObservableNames)

    experimentalData = createExperimentalData(observables, experimentalConditions, measurementData, modelData)

    modelParameters = createModelParameters(new_sys, prob, parameterBounds, experimentalConditions, measurementData, observables, experimentalData)

    dualModelParameters = createDualModelParameters(modelParameters)

    parameterSpace, numAllStartParameters, lowerBounds, upperBounds = createParameterSpace(modelParameters, parameterBounds)

    if typeof(solver) <: CompositeAlgorithm
        modelOutput_dual = createModelOutput(OrdinaryDiffEq.ODECompositeSolution, ForwardDiff.Dual, experimentalData, modelParameters)
    else
        modelOutput_dual = createModelOutput(ODESolution, ForwardDiff.Dual, experimentalData, modelParameters)
    end

    # Getting start values

    allStartParameters = getSamples(numAllStartParameters, lowerBounds, upperBounds, iStartPar)

    # Initialize functions

    keepInBounds! = (theta) -> keepInBounds_forwGrad_proto!(theta, parameterSpace, modelParameters)

    solveODESystem_dual = (u0Vector, iCond) -> solveODESystem_forwGrad_proto_model_Boehm_JProteomeRes2014(prob, solver, dualModelParameters.dualDynamicParametersVector, u0Vector, 
            modelData, modelOutput_dual, iCond)

    calcUnscaledObservable_dual = (iCond) -> calcUnscaledObservable_forwGrad_proto_model_Boehm_JProteomeRes2014(dualModelParameters.dualDynamicParametersVector, modelOutput_dual.sols, 
            modelOutput_dual.h_barForCondObs, modelData, experimentalData, iCond)

    calcScaledObservable_dual = (iCond) -> calcScaledObservable_forwGrad_proto(dualModelParameters.dualScaleVector, dualModelParameters.dualOffsetVector, 
            modelOutput_dual.h_barForCondObs, modelOutput_dual.h_hatForCondObs, modelParameters, modelData, experimentalData, iCond)

    calcCost_dual = (iCond) -> calcCost_forwGrad_proto(dualModelParameters.dualVarianceVector, modelOutput_dual.h_hatForCondObs, modelOutput_dual.costForCondObs, 
            modelParameters, experimentalData, iCond)

    updateAllDualParameterVectors = (p) -> updateAllDualParameterVectors_proto(modelParameters, dualModelParameters, modelData, p)

    allConditionsCost_dual = (p) -> allConditionsCost_forwGrad_proto(modelParameters, dualModelParameters, experimentalData, modelData,
            updateAllDualParameterVectors, solveODESystem_dual, calcUnscaledObservable_dual, calcScaledObservable_dual, calcCost_dual, p)

    result = DiffResults.DiffResult(0.0, 0.0)

    grad_g = (p) -> unbiasedGrad_forwGrad_proto(p, allConditionsCost_dual, parameterSpace, modelParameters, result)

    stepLengths = exp10.(range(stepRange[1], stepRange[2], n_it))
    adam_opt = createAdam(allStartParameters, grad_g, stepLengths; b2 = b2)
    step = () -> step_forwGrad_proto!(adam_opt, keepInBounds!)

    return step, adam_opt, parameterSpace.doLogSearch
end
