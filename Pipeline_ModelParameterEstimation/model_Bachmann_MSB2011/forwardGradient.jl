

mutable struct Adam
    theta::AbstractArray{Float64} # Parameter array
    grad::Function                # Gradient function
    loss::Float64
    m::AbstractArray{Float64}     # First moment
    v::AbstractArray{Float64}     # Second moment
    b1::Float64                   # Exp. decay first moment
    b2::Float64                   # Exp. decay second moment
    a::Vector{Float64}                    # list of (decaying) step sizes
    eps::Float64                  # Epsilon for stability 
    t::Int                        # Time step (iteration)
    β::Float64                    # Reduction factor (scales the step length)
    r::Float64                    # Reduction multiplier (shortens the step length at failed step)
    c::Float64                    # Increas multiplyer (Increases the step length at successful step, up to 1.0)
    fail::Int64                   # Number of infeasible points tested in a row
end
  
# Outer constructor
function Adam(theta::AbstractArray{Float64}, grad::Function, a::Vector{Float64}; b2=0.999::Float64)
    loss = 0.0
    m   = zeros(size(theta))
    v   = zeros(size(theta))
    b1  = 0.9
    a   = a
    eps = 1e-8
    t   = 0
    β   = 1.0
    r   = 0.2
    c   = 1.3
    fail = 0
    Adam(theta, grad, loss, m, v, b1, b2, a, eps, t, β, r, c, fail)
end

function keepInBounds_forwGrad_proto!(theta, parameterSpace, modelParameters)
    slb = parameterSpace.scaleLowerBound
    sub = parameterSpace.scaleUpperBound
    scaleIndices = modelParameters.scaleIndices
    scale = view(theta, scaleIndices)
    scale .= min.(max.(scale, slb), sub)
    olb = parameterSpace.offsetLowerBound
    oub = parameterSpace.offsetUpperBound
    offsetIndices = modelParameters.offsetIndices
    offset = view(theta, offsetIndices)
    offset .= min.(max.(offset, olb), oub)
    vlb = parameterSpace.varianceLowerBound
    vub = parameterSpace.varianceUpperBound
    varianceIndices = modelParameters.varianceIndices
    variance  = view(theta, varianceIndices)
    variance .= min.(max.(variance, vlb), vub)
    plb = parameterSpace.parameterLowerBound
    pub = parameterSpace.parameterUpperBound
    parameterIndices = modelParameters.parameterIndices
    pars = view(theta, parameterIndices)
    pars .= min.(max.(pars, plb), pub)
end

# Step function with optional keyword arguments for the data passed to grad()
function step_forwGrad_proto!(opt::Adam, keepInBounds!)
    opt.t += 1
    gt, loss    = opt.grad(opt.theta)

    if loss == Inf
        opt.fail += 1
        opt.β *= opt.r 
    else
        opt.β = min(opt.β * opt.c, 1.0)
        opt.fail = 0
        opt.loss = loss
        opt.m = opt.b1 .* opt.m + (1 - opt.b1) .* gt
        opt.v = opt.b2 .* opt.v + (1 - opt.b2) .* gt .^ 2
        mhat = opt.m ./ (1 - opt.b1^opt.t)
        vhat = opt.v ./ (1 - opt.b2^opt.t)
        opt.theta -= opt.β * opt.a[opt.t] * (mhat ./ (sqrt.(vhat) .+ opt.eps))

        keepInBounds!(opt.theta)

    end
end

function solveODESystem_forwGrad_proto(prob, modelParameters, dualModelParameters, modelData, modelOutput, iCond, p)
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

    try
        modelOutput.sols[iCond] = OrdinaryDiffEq.solve(_prob, Rodas5(), reltol=1e-9, abstol=1e-9)
        return false
    catch err
        return true
    end

end

function calcUnscaledObservable_forwGrad_proto(modelParameters, dualModelParameters, modelData, experimentalData, modelOutput, iCond, p)
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

function calcScaledObservable_forwGrad_proto(modelParameters, dualModelParameters, modelData, modelOutput, experimentalData, p)
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

function calcCost_forwGrad_proto(modelParameters, dualModelParameters, modelOutput, experimentalData, p)
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

function allConditionsCost_forwGrad_proto(modelParameters, experimentalData, modelData, solveODESystem, calcUnscaledObservable, calcScaledObservable, calcCost, p)
    dynParVector = modelParameters.dynamicParametersVector
    inputPVFC = experimentalData.inputParameterValuesForCond
    inputPI = modelData.inputParameterIndices

    for (iCond, inputParameterValues) in enumerate(inputPVFC)
        dynParVector[inputPI] = inputParameterValues
        fail = solveODESystem(iCond, p)
        if fail
            # If infeasible parameters, return dual (Inf, NaN)
            return ForwardDiff.Dual{ForwardDiff.Tag{eltype(p), Float64},Float64, 1}(Inf, ForwardDiff.Partials{1, Float64}((NaN,)))
        end
        calcUnscaledObservable(iCond, p)
    end

    calcScaledObservable(p)

    cost = calcCost(p)

    return cost
end


# F - function on p 
function unbiasedGrad_forwGrad_proto(p::T1, f::Function, parameterSpace, modelParameters, result)::Tuple{Array{Float64, 1}, Float64} where T1<:Array{<:AbstractFloat, 1}
    doLogSearch = parameterSpace.doLogSearch
    allParameters = modelParameters.allParameters
    allParameters .= p 
    view(allParameters, doLogSearch) .= exp10.(view(allParameters, doLogSearch))

    #updateAllParameterVectors()

    v::Array{Float64, 1} = randn(length(p)) # v ~ MV(0, I)

    g = r -> f(allParameters + r*v)
    result = ForwardDiff.derivative!(result, g, 0.0)
    uGrad = DiffResults.gradient(result) * v
    view(uGrad, doLogSearch) .*= view(allParameters, doLogSearch) * log(10)
    loss = DiffResults.value(result)
    return uGrad, loss
end





function forwardGradient_proto(iStartPar, n_it, b2, stepRange, 
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

    usedType = ForwardDiff.Dual
    modelOutput = ModelOutput(usedType, experimentalData, modelParameters)

    # Getting start values

    allStartParameters = getSamples(numAllStartParameters, lowerBounds, upperBounds, iStartPar)

    # Initialize functions

    keepInBounds! = (theta) -> keepInBounds_forwGrad_proto!(theta, parameterSpace, modelParameters)

    solveODESystem = (iCond, p) -> solveODESystem_forwGrad_proto(prob, modelParameters, dualModelParameters, modelData, modelOutput, iCond, p)

    calcUnscaledObservable = (iCond, p) -> calcUnscaledObservable_forwGrad_proto(modelParameters, dualModelParameters, modelData, experimentalData, modelOutput, iCond, p)

    calcScaledObservable = (p) -> calcScaledObservable_forwGrad_proto(modelParameters, dualModelParameters, modelData, modelOutput, experimentalData, p)

    calcCost = (p) -> calcCost_forwGrad_proto(modelParameters, dualModelParameters, modelOutput, experimentalData, p)

    allConditionsCost = (p) -> allConditionsCost_forwGrad_proto(modelParameters, experimentalData, modelData,
            solveODESystem, calcUnscaledObservable, calcScaledObservable, calcCost, p)

    result = DiffResults.DiffResult(0.0, 0.0)

    grad_g = (p) -> unbiasedGrad_forwGrad_proto(p, allConditionsCost, parameterSpace, modelParameters, result)

    stepLengths = exp10.(range(stepRange[1], stepRange[1], n_it))
    adam_opt = Adam(allStartParameters, grad_g, stepLengths; b2 = b2)
    step = () -> step_forwGrad_proto!(adam_opt, keepInBounds!)
    
    return step, adam_opt, parameterSpace.doLogSearch
end



