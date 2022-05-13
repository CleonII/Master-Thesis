

mutable struct Adam{
        T1 <: AbstractArray{Float64},
        T2 <: Function,
        T3 <: Float64,
        T4 <: Vector{Float64},
        T5 <: Int64
        }
    theta::T1 # Parameter array
    grad::T2                # Gradient function
    loss::T3
    m::T1     # First moment
    v::T1     # Second moment
    b1::T3                   # Exp. decay first moment
    b2::T3                   # Exp. decay second moment
    a::T4                    # list of (decaying) step sizes
    eps::T3                  # Epsilon for stability 
    t::T5                        # Time step (iteration)
    β::T3                    # Reduction factor (scales the step length)
    r::T3                    # Reduction multiplier (shortens the step length at failed step)
    c::T3                    # Increas multiplyer (Increases the step length at successful step, up to 1.0)
    fail::T5                   # Number of infeasible points tested in a row
end
  
# Outer constructor
function createAdam(theta::AbstractArray{Float64}, grad::Function, a::Vector{Float64}; b2=0.999::Float64)
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

function keepInBounds_forwGrad_proto!(theta::Vector{Float64}, parameterSpace::ParameterSpace, modelParameters::ModelParameters)
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
function step_forwGrad_proto!(opt::Adam, keepInBounds!::Function)
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

    nothing
end

function solveODESystem_forwGrad_proto(prob::ODEProblem, solver, dynParVector::Vector{T1}, u0Vector::Vector{T2}, modelData::ModelData, 
        modelOutput::ModelOutput, iCond::Int64)::Bool where {T1 <: ForwardDiff.Dual, T2 <: ForwardDiff.Dual}

    initVariable = modelData.initVariableIndices
    parameterInU0Indices = modelData.parameterInU0Indices
    u0Vector[initVariable[1]] = dynParVector[parameterInU0Indices[1]]
    u0Vector[initVariable[2]] = dynParVector[parameterInU0Indices[2]] * (dynParVector[parameterInU0Indices[3]] * dynParVector[parameterInU0Indices[4]] + 1)
    u0Vector[initVariable[3]] = dynParVector[parameterInU0Indices[5]]
    u0Vector[initVariable[4]] = dynParVector[parameterInU0Indices[6]]
    u0Vector[initVariable[5]] = dynParVector[parameterInU0Indices[7]] * dynParVector[parameterInU0Indices[8]] * dynParVector[parameterInU0Indices[9]]
    u0Vector[initVariable[6]] = dynParVector[parameterInU0Indices[10]] * dynParVector[parameterInU0Indices[11]] * dynParVector[parameterInU0Indices[12]]

    _prob = remake(prob, u0 = u0Vector, p = dynParVector)

    try
        modelOutput.sols[iCond] = OrdinaryDiffEq.solve(_prob, solver, reltol=1e-9, abstol=1e-9)
        if minimum(modelOutput.sols[iCond][:,:]).value < 0.0
            return true
        else
            return false
        end
    catch err
        return true
    end

end

function calcUnscaledObservable_forwGrad_proto(dynParVector::Vector{T1}, sols::Vector{T2}, h_barFCO::Array{Vector{T1}, 2}, modelData::ModelData, 
        experimentalData::ExperimentalData, iCond::Int64)::Nothing where {T1 <: ForwardDiff.Dual, T2 <: Union{ODESolution, OrdinaryDiffEq.ODECompositeSolution}}
    
    timeSteps = experimentalData.timeStepsForCond[iCond]
    sol = sols[iCond](timeSteps)

    observedAIFCO = experimentalData.observedAtIndexForCondObs

    optParIOI = modelData.parameterInObservableIndices

    h_barFCO[iCond, 1] = ( sol[Symbol("CISRNA(t)")] ./ dynParVector[optParIOI[1]] )[observedAIFCO[iCond, 1]]
    h_barFCO[iCond, 2] = ( sol[Symbol("CISRNA(t)")] ./ dynParVector[optParIOI[1]] )[observedAIFCO[iCond, 2]]
    h_barFCO[iCond, 3] = ( sol[Symbol("CISRNA(t)")] ./ dynParVector[optParIOI[1]] )[observedAIFCO[iCond, 3]]
    h_barFCO[iCond, 4] = ( sol[Symbol("CIS(t)")] )[observedAIFCO[iCond, 4]]
    h_barFCO[iCond, 5] = ( sol[Symbol("CIS(t)")] ./ dynParVector[optParIOI[2]] )[observedAIFCO[iCond, 5]]
    h_barFCO[iCond, 6] = ( sol[Symbol("CIS(t)")] ./ dynParVector[optParIOI[2]] )[observedAIFCO[iCond, 6]]
    h_barFCO[iCond, 7] = ( sol[Symbol("CIS(t)")] ./ dynParVector[optParIOI[2]] )[observedAIFCO[iCond, 7]]
    h_barFCO[iCond, 8] = ( sol[Symbol("SHP1(t)")] + sol[Symbol("SHP1Act(t)")] )[observedAIFCO[iCond, 8]]
    h_barFCO[iCond, 9] = ( sol[Symbol("SOCS3RNA(t)")] ./ dynParVector[optParIOI[3]] )[observedAIFCO[iCond, 9]]
    h_barFCO[iCond, 10] = ( sol[Symbol("SOCS3RNA(t)")] ./ dynParVector[optParIOI[3]] )[observedAIFCO[iCond, 10]]
    h_barFCO[iCond, 11] = ( sol[Symbol("SOCS3RNA(t)")] ./ dynParVector[optParIOI[3]] )[observedAIFCO[iCond, 11]]
    h_barFCO[iCond, 12] = ( sol[Symbol("SOCS3(t)")] )[observedAIFCO[iCond, 12]]
    h_barFCO[iCond, 13] = ( sol[Symbol("SOCS3(t)")] ./ dynParVector[optParIOI[4]] )[observedAIFCO[iCond, 13]]
    h_barFCO[iCond, 14] = ( sol[Symbol("STAT5(t)")] )[observedAIFCO[iCond, 14]]
    h_barFCO[iCond, 15] = ( 16 * (sol[Symbol("p12EpoRpJAK2(t)")] + sol[Symbol("p1EpoRpJAK2(t)")] + 
            sol[Symbol("p2EpoRpJAK2(t)")]) / dynParVector[optParIOI[5]] )[observedAIFCO[iCond, 15]]
    h_barFCO[iCond, 16] = ( 2 * (sol[Symbol("EpoRpJAK2(t)")] + sol[Symbol("p12EpoRpJAK2(t)")] + 
            sol[Symbol("p1EpoRpJAK2(t)")] + sol[Symbol("p2EpoRpJAK2(t)")]) / dynParVector[optParIOI[5]] )[observedAIFCO[iCond, 16]]
    h_barFCO[iCond, 17] = ( (100 * sol[Symbol("pSTAT5(t)")]) ./ (sol[Symbol("STAT5(t)")] + sol[Symbol("pSTAT5(t)")]) )[observedAIFCO[iCond, 17]]
    h_barFCO[iCond, 18] = ( sol[Symbol("pSTAT5(t)")] / dynParVector[optParIOI[6]] )[observedAIFCO[iCond, 18]]
    h_barFCO[iCond, 19] = ( (sol[Symbol("SHP1(t)")] + sol[Symbol("SHP1Act(t)")]) / dynParVector[optParIOI[7]] )[observedAIFCO[iCond, 19]]
    h_barFCO[iCond, 20] = ( (sol[Symbol("STAT5(t)")] + sol[Symbol("pSTAT5(t)")]) / dynParVector[optParIOI[6]] )[observedAIFCO[iCond, 20]]

    nothing
end

function calcScaledObservable_forwGrad_proto(scaleVector::Vector{T1}, offsetVector::Vector{T1}, h_barFCO::Array{Vector{T1}, 2}, h_hatFCO::Array{Vector{T1}, 2}, 
        modelParameters::ModelParameters, modelData::ModelData, experimentalData::ExperimentalData, iCond::Int64)::Nothing where {T1 <: ForwardDiff.Dual}

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

function calcCost_forwGrad_proto(varianceVector::Vector{T1}, h_hatFCO::Array{Vector{T1}, 2}, costFCO::Array{T1, 2}, modelParameters::ModelParameters, 
        experimentalData::ExperimentalData, iCond::Int64)::T1 where T1 <: ForwardDiff.Dual

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

function allConditionsCost_forwGrad_proto(modelParameters::ModelParameters, dualModelParameters::DualModelParameters, experimentalData::ExperimentalData, 
        modelData::ModelData, updateAllDualParameterVectors::Function, solveODESystem::Function, calcUnscaledObservable::Function, calcScaledObservable::Function, 
        calcCost::Function, p::Vector{T1})::T1 where T1 <: ForwardDiff.Dual

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
            # If infeasible parameters, return dual (Inf, 0.0)
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
function unbiasedGrad_forwGrad_proto(p::T1, f::Function, parameterSpace::ParameterSpace, modelParameters::ModelParameters, 
        result)::Tuple{Array{Float64, 1}, Float64} where T1<:Array{<:AbstractFloat, 1}

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



function forwardGradient(modelFunction::Function, iStartPar::Int64, n_it::Int64, b2::Float64, stepRange::Vector{<:AbstractFloat}, solver, 
        timeEnd::Float64, experimentalConditions::DataFrame, measurementData::DataFrame, observables::DataFrame, parameterBounds::DataFrame)::Tuple{Function, Adam, Vector{Int64}}

    sys, initialSpeciesValues, trueParameterValues = modelFunction()
    new_sys = ode_order_lowering(sys)
    u0 = initialSpeciesValues
    pars = trueParameterValues 
    tspan = (0.0, timeEnd)
    prob = ODEProblem(new_sys,u0,tspan,pars,jac=true)

    initVariableNames = ["EpoRJAK2_CIS", "SHP1", "STAT5", "EpoRJAK2", "SOCS3", "CIS"]

    observableVariableNames = ["CISRNA", "CIS", "SHP1", "SHP1Act", "SOCS3RNA", "SOCS3", "STAT5", 
            "p12EpoRpJAK2", "p1EpoRpJAK2", "p2EpoRpJAK2", "EpoRpJAK2", "pSTAT5"]

    parameterInU0Names = ["init_EpoRJAK2_CIS", "init_SHP1", "init_SHP1_multiplier", "SHP1ProOE", "init_STAT5", 
            "init_EpoRJAK2", "init_SOCS3_multiplier", "SOCS3EqcOE", "SOCS3Eqc", "init_CIS_multiplier", "CISEqc", "CISEqcOE"]

    parameterInObservableNames = ["CISRNAEqc", "CISEqc", "SOCS3RNAEqc", "SOCS3Eqc", "init_EpoRJAK2", "init_STAT5", "init_SHP1"]

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

    solveODESystem_dual = (u0Vector, iCond) -> solveODESystem_forwGrad_proto(prob, solver, dualModelParameters.dualDynamicParametersVector, u0Vector, 
            modelData, modelOutput_dual, iCond)

    calcUnscaledObservable_dual = (iCond) -> calcUnscaledObservable_forwGrad_proto(dualModelParameters.dualDynamicParametersVector, modelOutput_dual.sols, 
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

    stepLengths = exp10.(range(stepRange[1], stepRange[1], n_it))
    adam_opt = createAdam(allStartParameters, grad_g, stepLengths; b2 = b2)
    step = () -> step_forwGrad_proto!(adam_opt, keepInBounds!)

    return step, adam_opt, parameterSpace.doLogSearch
end
