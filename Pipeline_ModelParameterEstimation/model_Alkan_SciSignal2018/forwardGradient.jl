

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

function solveODESystem_forwGrad_proto(prob, dynParVector, u0Vector, modelData, modelOutput, iCond)

    initVariable = modelData.initVariableIndices
    parameterInU0Indices = modelData.parameterInU0Indices
    u0Vector[initVariable[1]] = -dynParVector[parameterInU0Indices[1]] * dynParVector[parameterInU0Indices[2]] * (dynParVector[parameterInU0Indices[3]] - 1)
    u0Vector[initVariable[2]] = dynParVector[parameterInU0Indices[1]] * dynParVector[parameterInU0Indices[3]]
    u0Vector[initVariable[3]] = dynParVector[parameterInU0Indices[1]] * (dynParVector[parameterInU0Indices[2]] - 1) * (dynParVector[parameterInU0Indices[3]] - 1)

    _prob = remake(prob, u0 = u0Vector, p = dynParVector)

    try
        modelOutput.sols[iCond] = OrdinaryDiffEq.solve(_prob, Rodas5(), reltol=1e-9, abstol=1e-9)
        if minimum(modelOutput.sols[iCond][:,:]).value < 0.0
            return true
        else
            return false
        end
    catch err
        return true
    end

end

function calcUnscaledObservable_forwGrad_proto(dynParVector, modelData, experimentalData, modelOutput, iCond)
    
    timeSteps = experimentalData.timeStepsForCond[iCond]
    sol = modelOutput.sols[iCond](timeSteps)

    h_barFCO = modelOutput.h_barForCondObs
    observedAIFCO = experimentalData.observedAtIndexForCondObs

    h_barFCO[iCond, 1] = ( sol[Symbol("Cells_Apo4(t)")] + sol[Symbol("Cells_Apo_ReOx(t)")] )[observedAIFCO[iCond, 1]]
    h_barFCO[iCond, 2] = ( sol[Symbol("Cells(t)")] + sol[Symbol("Cells_Cycle_S(t)")] + sol[Symbol("Cells_Cycle_G2(t)")] + 
            sol[Symbol("Cells_SSBDamage_S(t)")] + sol[Symbol("Cells_DSBDamage_S(t)")] + sol[Symbol("Cells_DSBDamage_G2(t)")] + 
            sol[Symbol("Cells_Apo(t)")] + sol[Symbol("Cells_Apo1(t)")] + sol[Symbol("Cells_Apo2(t)")] + 
            sol[Symbol("Cells_Apo3(t)")] + sol[Symbol("Cells_Apo4(t)")] +  + sol[Symbol("Cells_Apo_ReOx(t)")] )[observedAIFCO[iCond, 2]]
    h_barFCO[iCond, 3] = ( sol[Symbol("Wip1_mRNA_G2(t)")] + sol[Symbol("Wip1_mRNA_S(t)")] )[observedAIFCO[iCond, 3]]
    h_barFCO[iCond, 4] = ( sol[Symbol("p21_mRNA_S(t)")] + sol[Symbol("p21_mRNA_G2(t)")] )[observedAIFCO[iCond, 4]]
    h_barFCO[iCond, 5] = ( sol[Symbol("pATM_G2(t)")] + sol[Symbol("pATM_S(t)")] )[observedAIFCO[iCond, 5]]
    h_barFCO[iCond, 6] = ( sol[Symbol("pChk1_G2(t)")] + sol[Symbol("pChk1_S(t)")] )[observedAIFCO[iCond, 6]]
    h_barFCO[iCond, 7] = ( sol[Symbol("pChk2_G2(t)")] + sol[Symbol("pChk2_S(t)")] )[observedAIFCO[iCond, 7]]
    h_barFCO[iCond, 8] = ( sol[Symbol("pDNAPK_G2(t)")] + sol[Symbol("pDNAPK_S(t)")] )[observedAIFCO[iCond, 8]]
    h_barFCO[iCond, 9] = ( sol[Symbol("pp53_S(t)")] + sol[Symbol("pp53_G2(t)")] )[observedAIFCO[iCond, 9]]
    h_barFCO[iCond, 10] = ( sol[Symbol("p21_S(t)")] + sol[Symbol("p21_G2(t)")] )[observedAIFCO[iCond, 10]]
    h_barFCO[iCond, 11] = ( sol[Symbol("pp53_S(t)")] + sol[Symbol("pp53_G2(t)")] )[observedAIFCO[iCond, 11]]
    h_barFCO[iCond, 12] = ( sol[Symbol("yH2AX_S(t)")] + sol[Symbol("yH2AX_G2(t)")] )[observedAIFCO[iCond, 12]]

    nothing
end

function calcScaledObservable_forwGrad_proto(scaleVector, offsetVector, modelParameters, modelData, modelOutput, experimentalData, iCond)

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

function calcCost_forwGrad_proto(varianceVector, modelParameters, modelOutput, experimentalData, iCond)

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

function allConditionsCost_forwGrad_proto(modelParameters, dualModelParameters, experimentalData, modelData, updateAllDualParameterVectors, solveODESystem, calcUnscaledObservable, calcScaledObservable, calcCost, p)

    updateAllDualParameterVectors(p)

    u0Vector = modelParameters.u0Vector
    dualU0Vector = convert(Vector{eltype(p)}, u0Vector)

    dualDynParVector = dualModelParameters.dualDynamicParametersVector
    inputPVFC = experimentalData.inputParameterValuesForCond
    inputPI = modelData.inputParameterIndices

    cost = 0.0
    for (iCond, inputParameterValues) in enumerate(inputPVFC)
        dualDynParVector[inputPI] = convert.(eltype(p), inputParameterValues)
        fail = solveODESystem(dualU0Vector, iCond)
        if fail
            # If infeasible parameters, return dual (Inf, NaN)
            println("infeasible step")
            return convert(eltype(p), Inf)
        end
        calcUnscaledObservable(iCond)
        calcScaledObservable(iCond)
        cost += calcCost(iCond)
    end

    return cost
end


# F - function on p 
function unbiasedGrad_forwGrad_proto(p::T1, f::Function, parameterSpace, modelParameters, result)::Tuple{Array{Float64, 1}, Float64} where T1<:Array{<:AbstractFloat, 1}
    doLogSearch = parameterSpace.doLogSearch
    allParameters = modelParameters.allParameters
    allParameters .= p 
    view(allParameters, doLogSearch) .= exp10.(view(allParameters, doLogSearch))

    v::Array{Float64, 1} = randn(length(p)) # v ~ MV(0, I)

    g = r -> f(allParameters + r*v)
    result = ForwardDiff.derivative!(result, g, 0.0)
    uGrad = DiffResults.gradient(result) * v
    view(uGrad, doLogSearch) .*= view(allParameters, doLogSearch) * log(10)
    loss = DiffResults.value(result)

    return uGrad, loss
end



function forwardGradient(modelFunction, iStartPar, n_it, b2, stepRange, 
        timeEnd, experimentalConditions, measurementData, observables, parameterBounds)

    sys, initialSpeciesValues, trueParameterValues = modelFunction()
    new_sys = ode_order_lowering(sys)
    u0 = initialSpeciesValues
    pars = trueParameterValues 
    tspan = (0.0, timeEnd)
    prob = ODEProblem(new_sys,u0,tspan,pars,jac=true)

    initVariableNames = ["Cells_Cycle_G2", "Cells_Cycle_S", "Cells"]

    observableVariableNames = ["Cells_Apo4", "Cells_Apo_ReOx", "Cells", "Cells_Cycle_S", "Cells_Cycle_G2", "Cells_SSBDamage_S", 
            "Cells_DSBDamage_S", "Cells_DSBDamage_G2", "Cells_Apo", "Cells_Apo1", "Cells_Apo2", "Cells_Apo3",
            "Wip1_mRNA_G2", "Wip1_mRNA_S", "p21_mRNA_S", "p21_mRNA_G2", "pATM_G2", "pATM_S", "pChk1_G2", "pChk1_S", "pChk2_S", "pChk2_G2", 
            "pDNAPK_G2", "pDNAPK_S", "pp53_S", "pp53_G2", "p21_S", "p21_G2", "yH2AX_S", "yH2AX_G2"]

    parameterInU0Names = ["init_Cells", "init_Cells_Cycle_G2_rel", "init_Cells_Cycle_S_rel"]

    parameterInObservableNames = []

    # Initialize structs

    modelData = ModelData(new_sys, prob, observables, experimentalConditions, 
            initVariableNames, observableVariableNames, parameterInU0Names, parameterInObservableNames)

    experimentalData = ExperimentalData(observables, experimentalConditions, measurementData, modelData)

    modelParameters = ModelParameters(new_sys, prob, parameterBounds, experimentalConditions, measurementData, observables, experimentalData)

    dualModelParameters = DualModelParameters(modelParameters)

    parameterSpace, numAllStartParameters, lowerBounds, upperBounds = ParameterSpace(modelParameters, parameterBounds)

    modelOutput_dual = ModelOutput(ForwardDiff.Dual, experimentalData, modelParameters)

    # Getting start values

    allStartParameters = getSamples(numAllStartParameters, lowerBounds, upperBounds, iStartPar)

    # Initialize functions

    keepInBounds! = (theta) -> keepInBounds_forwGrad_proto!(theta, parameterSpace, modelParameters)

    solveODESystem_dual = (u0Vector, iCond) -> solveODESystem_forwGrad_proto(prob, dualModelParameters.dualDynamicParametersVector, u0Vector, modelData, modelOutput_dual, iCond)

    calcUnscaledObservable_dual = (iCond) -> calcUnscaledObservable_forwGrad_proto(dualModelParameters.dualDynamicParametersVector, modelData, experimentalData, modelOutput_dual, iCond)

    calcScaledObservable_dual = (iCond) -> calcScaledObservable_forwGrad_proto(dualModelParameters.dualScaleVector, dualModelParameters.dualOffsetVector, modelParameters, modelData, modelOutput_dual, experimentalData, iCond)

    calcCost_dual = (iCond) -> calcCost_forwGrad_proto(dualModelParameters.dualVarianceVector, modelParameters, modelOutput_dual, experimentalData, iCond)

    updateAllDualParameterVectors = (p) -> updateAllDualParameterVectors_proto(modelParameters, dualModelParameters, modelData, p)

    allConditionsCost_dual = (p) -> allConditionsCost_forwGrad_proto(modelParameters, dualModelParameters, experimentalData, modelData,
            updateAllDualParameterVectors, solveODESystem_dual, calcUnscaledObservable_dual, calcScaledObservable_dual, calcCost_dual, p)

    result = DiffResults.DiffResult(0.0, 0.0)

    grad_g = (p) -> unbiasedGrad_forwGrad_proto(p, allConditionsCost_dual, parameterSpace, modelParameters, result)

    stepLengths = exp10.(range(stepRange[1], stepRange[1], n_it))
    adam_opt = Adam(allStartParameters, grad_g, stepLengths; b2 = b2)
    step = () -> step_forwGrad_proto!(adam_opt, keepInBounds!)

    return step, adam_opt, parameterSpace.doLogSearch
end