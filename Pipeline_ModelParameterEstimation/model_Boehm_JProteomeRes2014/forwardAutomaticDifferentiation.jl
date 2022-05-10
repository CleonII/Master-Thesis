

function solveODESystem_forwAD_proto(prob, dynParVector, u0Vector, modelData, modelOutput, iCond)

    initVariable = modelData.initVariableIndices
    parameterInU0Indices = modelData.parameterInU0Indices
    u0Vector[initVariable[1]] = 207.6 * dynParVector[parameterInU0Indices[1]]
    u0Vector[initVariable[2]] = 207.6 - 207.6 * dynParVector[parameterInU0Indices[1]]

    _prob = remake(prob, u0 = u0Vector, p = dynParVector)
    modelOutput.sols[iCond] = OrdinaryDiffEq.solve(_prob, Rodas4P(), reltol=1e-9, abstol=1e-9)

    nothing
end


function calcUnscaledObservable_forwAD_proto(dynParVector, modelData, experimentalData, modelOutput, iCond)

    timeSteps = experimentalData.timeStepsForCond[iCond]
    sol = modelOutput.sols[iCond](timeSteps)

    h_barFCO = modelOutput.h_barForCondObs
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

function calcScaledObservable_forwAD_proto(scaleVector, offsetVector, modelParameters, modelData, modelOutput, experimentalData, iCond)

    scaleMap = modelParameters.scaleMap
    offsetMap = modelParameters.offsetMap

    observedOFC = experimentalData.observedObservableForCond
    
    h_barFCO = modelOutput.h_barForCondObs
    h_hatFCO = modelOutput.h_hatForCondObs
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

function calcCost_forwAD_proto(varianceVector, modelParameters, modelOutput, experimentalData, iCond)

    varianceMap = modelParameters.varianceMap

    h_hatFCO = modelOutput.h_hatForCondObs
    numDataFCO = experimentalData.numDataForCondObs
    
    costFCO = modelOutput.costForCondObs
    measurementFCO = experimentalData.measurementForCondObs
    observedOFC = experimentalData.observedObservableForCond

    for iObs in observedOFC[iCond]
        costFCO[iCond, iObs] = log(2*pi*varianceVector[varianceMap[iCond, iObs]]) * numDataFCO[iCond, iObs] + (dot(measurementFCO[iCond, iObs], measurementFCO[iCond, iObs]) - 
            2*dot(measurementFCO[iCond, iObs], h_hatFCO[iCond, iObs]) + 
            dot(h_hatFCO[iCond, iObs], h_hatFCO[iCond, iObs])) / (varianceVector[varianceMap[iCond, iObs]])
    end

    return sum(costFCO[iCond, observedOFC[iCond]])
end

function allConditionsCost_forwAD_float_proto(modelParameters, experimentalData, modelData, updateAllParameterVectors, solveODESystem, calcUnscaledObservable, calcScaledObservable, calcCost)

    updateAllParameterVectors()

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

function allConditionsCost_forwAD_dual_proto(modelParameters, dualModelParameters, experimentalData, modelData, updateAllDualParameterVectors, solveODESystem, calcUnscaledObservable, calcScaledObservable, calcCost, p)

    updateAllDualParameterVectors(p)

    u0Vector = modelParameters.u0Vector
    dualU0Vector = convert(Vector{eltype(p)}, u0Vector)

    dualDynParVector = dualModelParameters.dualDynamicParametersVector
    inputPVFC = experimentalData.inputParameterValuesForCond
    inputPI = modelData.inputParameterIndices

    cost = 0.0
    for (iCond, inputParameterValues) in enumerate(inputPVFC)
        dualDynParVector[inputPI] = convert.(eltype(p), inputParameterValues)
        solveODESystem(dualU0Vector, iCond)
        calcUnscaledObservable(iCond)
        calcScaledObservable(iCond)
        cost += calcCost(iCond)
    end

    return cost
end


function f_cost_forwAD_proto(allConditionsCost::Function, parameterSpace, modelParameters, p...)::Float64
    doLogSearch = parameterSpace.doLogSearch
    allPar = modelParameters.allParameters
    allPar .= p
    view(allPar, doLogSearch) .= exp10.(view(allPar, doLogSearch))

    cost = allConditionsCost()

    return cost
end

function f_grad_forwAD_proto(grad, result::DiffResults.MutableDiffResult, allConditionsCost::Function, cfg::GradientConfig, parameterSpace, modelParameters, modelOutput, p...)
    doLogSearch = parameterSpace.doLogSearch
    allPar = modelParameters.allParameters
    allPar .= p
    view(allPar, doLogSearch) .= exp10.(view(allPar, doLogSearch))

    ForwardDiff.gradient!(result, allConditionsCost, allPar, cfg)

    allParametersGrad = modelOutput.allParametersGrad
    allParametersGrad[:] = DiffResults.gradient(result)
    view(allParametersGrad, doLogSearch) .*= view(allPar, doLogSearch) * log(10)

    grad[:] = allParametersGrad

    nothing
end






function forwardAutomaticDifferentiation(modelFunction, iStartPar, chunkSize, optAlg,
        timeEnd, experimentalConditions, measurementData, observables, parameterBounds)

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

    modelData = ModelData(new_sys, prob, observables, experimentalConditions, 
            initVariableNames, observableVariableNames, parameterInU0Names, parameterInObservableNames)

    experimentalData = ExperimentalData(observables, experimentalConditions, measurementData, modelData)

    modelParameters = ModelParameters(new_sys, prob, parameterBounds, experimentalConditions, measurementData, observables, experimentalData)

    dualModelParameters = DualModelParameters(modelParameters)

    parameterSpace, numAllStartParameters, lowerBounds, upperBounds = ParameterSpace(modelParameters, parameterBounds)

    modelOutput_float = ModelOutput(Float64, experimentalData, modelParameters)

    modelOutput_dual = ModelOutput(ForwardDiff.Dual, experimentalData, modelParameters)

    # Getting start values

    allStartParameters = getSamples(numAllStartParameters, lowerBounds, upperBounds, iStartPar)

    # Initialize functions

    solveODESystem_float = (iCond) -> solveODESystem_forwAD_proto(prob, modelParameters.dynamicParametersVector, modelParameters.u0Vector, modelData, modelOutput_float, iCond)

    solveODESystem_dual = (dualU0Vector, iCond) -> solveODESystem_forwAD_proto(prob, dualModelParameters.dualDynamicParametersVector, dualU0Vector, modelData, modelOutput_dual, iCond)

    calcUnscaledObservable_float = (iCond) -> calcUnscaledObservable_forwAD_proto(modelParameters.dynamicParametersVector, modelData, experimentalData, modelOutput_float, iCond)

    calcUnscaledObservable_dual = (iCond) -> calcUnscaledObservable_forwAD_proto(dualModelParameters.dualDynamicParametersVector, modelData, experimentalData, modelOutput_dual, iCond)

    calcScaledObservable_float = (iCond) -> calcScaledObservable_forwAD_proto(modelParameters.scaleVector, modelParameters.offsetVector, modelParameters, modelData, modelOutput_float, experimentalData, iCond)

    calcScaledObservable_dual = (iCond) -> calcScaledObservable_forwAD_proto(dualModelParameters.dualScaleVector, dualModelParameters.dualOffsetVector, modelParameters, modelData, modelOutput_dual, experimentalData, iCond)

    calcCost_float = (iCond) -> calcCost_forwAD_proto(modelParameters.varianceVector, modelParameters, modelOutput_float, experimentalData, iCond)

    calcCost_dual = (iCond) -> calcCost_forwAD_proto(dualModelParameters.dualVarianceVector, modelParameters, modelOutput_dual, experimentalData, iCond)

    updateAllParameterVectors = () -> updateAllParameterVectors_proto(modelParameters, modelData)

    updateAllDualParameterVectors = (p) -> updateAllDualParameterVectors_proto(modelParameters, dualModelParameters, modelData, p)

    allConditionsCost_float = () -> allConditionsCost_forwAD_float_proto(modelParameters, experimentalData, modelData, 
            updateAllParameterVectors, solveODESystem_float, calcUnscaledObservable_float, calcScaledObservable_float, calcCost_float)

    allConditionsCost_dual = (p) -> allConditionsCost_forwAD_dual_proto(modelParameters, dualModelParameters, experimentalData, modelData,
            updateAllDualParameterVectors, solveODESystem_dual, calcUnscaledObservable_dual, calcScaledObservable_dual, calcCost_dual, p)

    cfg = GradientConfig(allConditionsCost_dual, allStartParameters, Chunk{chunkSize}())

    result = DiffResults.GradientResult(allStartParameters::Vector{Float64})

    f = (p_tuple...) -> f_cost_forwAD_proto(allConditionsCost_float, parameterSpace, modelParameters, p_tuple...)
    f_grad = (grad, p_tuple...) -> f_grad_forwAD_proto(grad, result, allConditionsCost_dual, cfg, parameterSpace, modelParameters, modelOutput_dual, p_tuple...)

    if optAlg == :Ipopt
        model = Model(Ipopt.Optimizer)
        set_optimizer_attribute(model, "print_level", 0)
        set_optimizer_attribute(model, "max_iter", 500)
        set_optimizer_attribute(model, "hessian_approximation", "limited-memory")
        set_optimizer_attribute(model, "tol", 1e-6)
        set_optimizer_attribute(model, "acceptable_tol", 1e-4)
        
    else
        model = Model(NLopt.Optimizer)
        #model.moi_backend.optimizer.model.options
        set_optimizer_attribute(model, "algorithm", optAlg)
        set_optimizer_attribute(model, "maxeval", 500)
        set_optimizer_attribute(model, "ftol_rel", 1e-6)
        set_optimizer_attribute(model, "xtol_rel", 1e-4)
    end

    #JuMP.register(model::Model, s::Symbol, dimension::Integer, f::Function, ∇f::Function, ∇²f::Function)
    register(model, :f, numAllStartParameters, f, f_grad)
    
    @variable(model, lowerBounds[i] <= p[i = 1:numAllStartParameters] <= upperBounds[i]) # fix lower and upper bounds
    for i in 1:numAllStartParameters
        set_start_value(p[i], allStartParameters[i])
    end
    @NLobjective(model, Min, f(p...))
    
    return model, p, parameterSpace.doLogSearch
    
end
