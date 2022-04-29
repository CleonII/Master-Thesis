

function solveODESystem_forwAD_proto(prob, modelParameters, dualModelParameters, modelData, modelOutput, iCond, p)
    parameterIndices = modelParameters.parameterIndices
    dynPars = view(p, parameterIndices)
    dynParVector = modelParameters.dynamicParametersVector
    dualDynParVector = dualModelParameters.dualDynamicParametersVector
    dualDynParVector .= convert.(eltype(dynPars), dynParVector)
    dualDynParVector[modelData.optParameterIndices] .= dynPars

    dualU0 = convert(Vector{eltype(p)}, prob.u0)
    initVariable = modelData.initVariableIndices
    parameterInU0Indices = modelData.parameterInU0Indices
    dualU0[initVariable[1]] = dualDynParVector[parameterInU0Indices[1]]
    dualU0[initVariable[2]] = dualDynParVector[parameterInU0Indices[2]] * (dualDynParVector[parameterInU0Indices[3]] * dualDynParVector[parameterInU0Indices[4]] + 1)
    dualU0[initVariable[3]] = dualDynParVector[parameterInU0Indices[5]]
    dualU0[initVariable[4]] = dualDynParVector[parameterInU0Indices[6]]
    dualU0[initVariable[5]] = dualDynParVector[parameterInU0Indices[7]] * dualDynParVector[parameterInU0Indices[8]] * dualDynParVector[parameterInU0Indices[9]]
    dualU0[initVariable[6]] = dualDynParVector[parameterInU0Indices[10]] * dualDynParVector[parameterInU0Indices[11]] * dualDynParVector[parameterInU0Indices[12]]

    _prob = remake(prob, u0 = dualU0, p = dualDynParVector)

    modelOutput.sols[iCond] = OrdinaryDiffEq.solve(_prob, Rodas5(), reltol=1e-9, abstol=1e-9)
end


function calcUnscaledObservable_forwAD_proto(modelParameters, dualModelParameters, modelData, experimentalData, modelOutput, iCond, p)
    parameterIndices = modelParameters.parameterIndices
    dynPars = view(p, parameterIndices)
    dynParVector = modelParameters.dynamicParametersVector
    dualDynParVector = dualModelParameters.dualDynamicParametersVector
    dualDynParVector .= convert.(eltype(dynPars), dynParVector)
    dualDynParVector[modelData.optParameterIndices] .= dynPars

    timeSteps = experimentalData.timeStepsForCond[iCond]
    sol = modelOutput.sols[iCond](timeSteps)

    h_barFCO = modelOutput.h_barForCondObs
    observedAIFCO = experimentalData.observedAtIndexForCondObs

    optParIOI = modelData.parameterInObservableIndices

    h_barFCO[iCond, 1] = ( sol[Symbol("CISRNA(t)")] ./ dualDynParVector[optParIOI[1]] )[observedAIFCO[iCond, 1]]
    h_barFCO[iCond, 2] = ( sol[Symbol("CISRNA(t)")] ./ dualDynParVector[optParIOI[1]] )[observedAIFCO[iCond, 2]]
    h_barFCO[iCond, 3] = ( sol[Symbol("CISRNA(t)")] ./ dualDynParVector[optParIOI[1]] )[observedAIFCO[iCond, 3]]
    h_barFCO[iCond, 4] = ( sol[Symbol("CIS(t)")] )[observedAIFCO[iCond, 4]]
    h_barFCO[iCond, 5] = ( sol[Symbol("CIS(t)")] ./ dualDynParVector[optParIOI[2]] )[observedAIFCO[iCond, 5]]
    h_barFCO[iCond, 6] = ( sol[Symbol("CIS(t)")] ./ dualDynParVector[optParIOI[2]] )[observedAIFCO[iCond, 6]]
    h_barFCO[iCond, 7] = ( sol[Symbol("CIS(t)")] ./ dualDynParVector[optParIOI[2]] )[observedAIFCO[iCond, 7]]
    h_barFCO[iCond, 8] = ( sol[Symbol("SHP1(t)")] + sol[Symbol("SHP1Act(t)")] )[observedAIFCO[iCond, 8]]
    h_barFCO[iCond, 9] = ( sol[Symbol("SOCS3RNA(t)")] ./ dualDynParVector[optParIOI[3]] )[observedAIFCO[iCond, 9]]
    h_barFCO[iCond, 10] = ( sol[Symbol("SOCS3RNA(t)")] ./ dualDynParVector[optParIOI[3]] )[observedAIFCO[iCond, 10]]
    h_barFCO[iCond, 11] = ( sol[Symbol("SOCS3RNA(t)")] ./ dualDynParVector[optParIOI[3]] )[observedAIFCO[iCond, 11]]
    h_barFCO[iCond, 12] = ( sol[Symbol("SOCS3(t)")] )[observedAIFCO[iCond, 12]]
    h_barFCO[iCond, 13] = ( sol[Symbol("SOCS3(t)")] ./ dualDynParVector[optParIOI[4]] )[observedAIFCO[iCond, 13]]
    h_barFCO[iCond, 14] = ( sol[Symbol("STAT5(t)")] )[observedAIFCO[iCond, 14]]
    h_barFCO[iCond, 15] = ( 16 * (sol[Symbol("p12EpoRpJAK2(t)")] + sol[Symbol("p1EpoRpJAK2(t)")] + 
            sol[Symbol("p2EpoRpJAK2(t)")]) / dualDynParVector[optParIOI[5]] )[observedAIFCO[iCond, 15]]
    h_barFCO[iCond, 16] = ( 2 * (sol[Symbol("EpoRpJAK2(t)")] + sol[Symbol("p12EpoRpJAK2(t)")] + 
            sol[Symbol("p1EpoRpJAK2(t)")] + sol[Symbol("p2EpoRpJAK2(t)")]) / dualDynParVector[optParIOI[5]] )[observedAIFCO[iCond, 16]]
    h_barFCO[iCond, 17] = ( (100 * sol[Symbol("pSTAT5(t)")]) ./ (sol[Symbol("STAT5(t)")] + sol[Symbol("pSTAT5(t)")]) )[observedAIFCO[iCond, 17]]
    h_barFCO[iCond, 18] = ( sol[Symbol("pSTAT5(t)")] / dualDynParVector[optParIOI[6]] )[observedAIFCO[iCond, 18]]
    h_barFCO[iCond, 19] = ( (sol[Symbol("SHP1(t)")] + sol[Symbol("SHP1Act(t)")]) / dualDynParVector[optParIOI[7]] )[observedAIFCO[iCond, 19]]
    h_barFCO[iCond, 20] = ( (sol[Symbol("STAT5(t)")] + sol[Symbol("pSTAT5(t)")]) / dualDynParVector[optParIOI[6]] )[observedAIFCO[iCond, 20]]

    nothing
end

function calcScaledObservable_forwAD_proto(modelParameters, dualModelParameters, modelData, modelOutput, experimentalData, p)
    scaleIndices = modelParameters.scaleIndices
    scale = view(p, scaleIndices)
    scaleVector = modelParameters.scaleVector
    dualScaleVector = dualModelParameters.dualScaleVector
    dualScaleVector .= convert.(eltype(scale), scaleVector)
    dualScaleVector[1:modelParameters.numScale] .= scale
    scaleMap = modelParameters.scaleMap

    offsetIndices = modelParameters.offsetIndices
    offset = view(p, offsetIndices)
    offsetVector = modelParameters.offsetVector
    dualOffsetVector = dualModelParameters.dualOffsetVector
    dualOffsetVector .= convert.(eltype(offset), offsetVector)
    dualOffsetVector[1:modelParameters.numOffset] .= offset
    offsetMap = modelParameters.offsetMap

    numConditions = experimentalData.numConditions
    observedOFC = experimentalData.observedObservableForCond
    
    h_barFCO = modelOutput.h_barForCondObs
    h_hatFCO = modelOutput.h_hatForCondObs
    observableLT = modelData.observableLogTransformation

    for iCond in 1:numConditions
        for iObs = observedOFC[iCond]
            if observableLT[iObs]
                h_hatFCO[iCond, iObs] = log10.(dualScaleVector[scaleMap[iCond, iObs]] * h_barFCO[iCond, iObs] .+ dualOffsetVector[offsetMap[iCond, iObs]])
            else
                h_hatFCO[iCond, iObs] = dualScaleVector[scaleMap[iCond, iObs]] * h_barFCO[iCond, iObs] .+ dualOffsetVector[offsetMap[iCond, iObs]]
            end
        end
    end

    nothing
end

function calcCost_forwAD_proto(modelParameters, dualModelParameters, modelOutput, experimentalData, p)
    varianceIndices = modelParameters.varianceIndices
    variance = view(p, varianceIndices)
    dualVarianceVector = dualModelParameters.dualVarianceVector
    dualVarianceVector[1:modelParameters.numVariance] .= variance
    varianceMap = modelParameters.varianceMap

    h_hatFCO = modelOutput.h_hatForCondObs
    numDataFCO = experimentalData.numDataForCondObs
    
    costFCO = modelOutput.costForCondObs
    measurementFCO = experimentalData.measurementForCondObs
    observedOFC = experimentalData.observedObservableForCond

    cost = 0.0
    for iCond in 1:experimentalData.numConditions
        for iObs in observedOFC[iCond]
            costFCO[iCond, iObs] = log(2*pi*dualVarianceVector[varianceMap[iCond, iObs]]) * numDataFCO[iCond, iObs] + (dot(measurementFCO[iCond, iObs], measurementFCO[iCond, iObs]) - 
                2*dot(measurementFCO[iCond, iObs], h_hatFCO[iCond, iObs]) + 
                dot(h_hatFCO[iCond, iObs], h_hatFCO[iCond, iObs])) / (2 * dualVarianceVector[varianceMap[iCond, iObs]])
        end

        cost += sum(costFCO[iCond, observedOFC[iCond]])
    end

    return cost
end

function allConditionsCost_forwAD_proto(modelParameters, experimentalData, modelData, solveODESystem, calcUnscaledObservable, calcScaledObservable, calcCost, p)
    dynParVector = modelParameters.dynamicParametersVector
    inputPVFC = experimentalData.inputParameterValuesForCond
    inputPI = modelData.inputParameterIndices

    for (iCond, inputParameterValues) in enumerate(inputPVFC)
        dynParVector[inputPI] = inputParameterValues
        solveODESystem(iCond, p)
        calcUnscaledObservable(iCond, p)
    end

    calcScaledObservable(p)

    cost = calcCost(p)

    return cost
end



function f_cost_forwAD_proto(result::DiffResults.MutableDiffResult, f::Function, cfg::GradientConfig, parameterSpace, modelParameters, modelOutput, p...)::Float64
    doLogSearch = parameterSpace.doLogSearch
    allPar = modelParameters.allParameters
    allPar .= p
    view(allPar, doLogSearch) .= exp10.(view(allPar, doLogSearch))

    ForwardDiff.gradient!(result, f, allPar, cfg)
    allParametersGrad = modelOutput.allParametersGrad
    allParametersGrad[:] = DiffResults.gradient(result)
    view(allParametersGrad, doLogSearch) .*= view(allPar, doLogSearch) * log(10)
    
    println(DiffResults.value(result), "\n")

    return DiffResults.value(result)
end

function f_grad_forwAD_proto(grad, modelOutput, p...)
    grad[:] = modelOutput.allParametersGrad

    nothing
end






function forwardAutomaticDifferentiation_proto(iStartPar, chunkSize, optAlg,
        filesAndPaths, timeEnd, experimentalConditions, measurementData, observables, parameterBounds)

    include(joinpath(filesAndPaths.modelPath, filesAndPaths.modelFile))
    new_sys = ode_order_lowering(sys)
    u0 = initialSpeciesValues
    p = trueParameterValues 
    c = trueConstantsValues
    pars = [p;c]
    tspan = (0.0, timeEnd)
    prob = ODEProblem(new_sys,u0,tspan,pars,jac=true)

    initVariableNames = ["EpoRJAK2_CIS", "SHP1", "STAT5", "EpoRJAK2", "SOCS3", "CIS"]

    observableVariableNames = ["CISRNA", "CIS", "SHP1", "SHP1Act", "SOCS3RNA", "SOCS3", "STAT5", 
            "p12EpoRpJAK2", "p1EpoRpJAK2", "p2EpoRpJAK2", "EpoRpJAK2", "pSTAT5"]

    parameterInU0Names = ["init_EpoRJAK2_CIS", "init_SHP1", "init_SHP1_multiplier", "SHP1ProOE", "init_STAT5", 
            "init_EpoRJAK2", "init_SOCS3_multiplier", "SOCS3EqcOE", "SOCS3Eqc", "init_CIS_multiplier", "CISEqc", "CISEqcOE"]

    parameterInObservableNames = ["CISRNAEqc", "CISEqc", "SOCS3RNAEqc", "SOCS3Eqc", "init_EpoRJAK2", "init_STAT5", "init_SHP1"]

    # Initialize structs

    modelData, inputParameterSymbols = ModelData(new_sys, prob, observables, experimentalConditions, 
            initVariableNames, observableVariableNames, parameterInU0Names, parameterInObservableNames)

    experimentalData = ExperimentalData(observables, experimentalConditions, measurementData, inputParameterSymbols)

    modelParameters = ModelParameters(new_sys, prob, parameterBounds, experimentalConditions, measurementData, observables, experimentalData)

    dualModelParameters = DualModelParameters(prob, parameterBounds)

    parameterSpace, numAllStartParameters, lowerBounds, upperBounds = ParameterSpace(modelParameters, parameterBounds)

    modelOutput = ModelOutput(ForwardDiff.Dual, experimentalData, modelParameters)

    # Getting start values

    allStartParameters = getSamples(numAllStartParameters, lowerBounds, upperBounds, iStartPar)

    # Initialize functions

    solveODESystem = (iCond, p) -> solveODESystem_forwAD_proto(prob, modelParameters, dualModelParameters, modelData, modelOutput, iCond, p)

    calcUnscaledObservable = (iCond, p) -> calcUnscaledObservable_forwAD_proto(modelParameters, dualModelParameters, modelData, experimentalData, modelOutput, iCond, p)

    calcScaledObservable = (p) -> calcScaledObservable_forwAD_proto(modelParameters, dualModelParameters, modelData, modelOutput, experimentalData, p)

    calcCost = (p) -> calcCost_forwAD_proto(modelParameters, dualModelParameters, modelOutput, experimentalData, p)

    allConditionsCost = (p) -> allConditionsCost_forwAD_proto(modelParameters, experimentalData, modelData,
            solveODESystem, calcUnscaledObservable, calcScaledObservable, calcCost, p)

    cfg = GradientConfig(allConditionsCost, allStartParameters, Chunk{chunkSize}())

    result = DiffResults.GradientResult(allStartParameters::Vector{Float64})

    f = (p_tuple...) -> f_cost_forwAD_proto(result, allConditionsCost, cfg, parameterSpace, modelParameters, modelOutput, p_tuple...)
    f_grad = (grad, p_tuple...) -> f_grad_forwAD_proto(grad, modelOutput, p_tuple...)

    model = Model(NLopt.Optimizer)
    #JuMP.register(model::Model, s::Symbol, dimension::Integer, f::Function, ∇f::Function, ∇²f::Function)
    register(model, :f, numAllStartParameters, f, f_grad)

    set_optimizer_attribute(model, "algorithm", optAlg)
    set_optimizer_attribute(model, "maxeval", 1000)
    set_optimizer_attribute(model, "ftol_rel", 1e-8)
    set_optimizer_attribute(model, "xtol_rel", 1e-6)
    
    @variable(model, lowerBounds[i] <= p[i = 1:numAllStartParameters] <= upperBounds[i]) # fix lower and upper bounds
    for i in 1:numAllStartParameters
        set_start_value(p[i], allStartParameters[i])
    end
    @NLobjective(model, Min, f(p...))
    
    return model, p, parameterSpace.doLogSearch
    
end
