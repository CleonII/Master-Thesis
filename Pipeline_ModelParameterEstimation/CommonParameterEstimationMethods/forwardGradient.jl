

mutable struct Adam{
        T1 <: AbstractArray{Float64},
        T2 <: Function,
        T3 <: Float64,
        T4 <: Vector{Float64},
        T5 <: Int64
        }
    theta::T1                # Parameter array
    theta_old::T1            # Old parameter array
    grad::T2                 # Gradient function
    loss::T3                 # value of the cost function
    m::T1                    # First moment
    v::T1                    # Second moment
    b1::T3                   # Exp. decay first moment
    b2::T3                   # Exp. decay second moment
    a::T4                    # list of (decaying) step sizes
    eps::T3                  # Epsilon for stability 
    t::T5                    # Time step (iteration)
    β::T3                    # Reduction factor (scales the step length)
    r::T3                    # Reduction multiplier (shortens the step length at failed step)
    c::T3                    # Increas multiplyer (Increases the step length at successful step, up to 1.0)
    fail::T5                 # Number of infeasible points tested in a row
end
  

# Outer constructor
function createAdam(theta::AbstractArray{Float64}, grad::Function, a::Vector{Float64}; b2=0.999::Float64)
    theta_old = zeros(size(theta))
    loss = 0.0
    b1  = 0.9
    m   = zeros(size(theta))
    v   = zeros(size(theta))
    a   = a
    eps = 1e-8
    t   = 0
    β   = 1.0
    r   = 0.2
    c   = 1.3
    fail = 0
    Adam(theta, theta_old, grad, loss, m, v, b1, b2, a, eps, t, β, r, c, fail)
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
        if opt.t == 1
            opt.fail = Inf
        else
            opt.fail += 1
            opt.β *= opt.r 
            opt.theta .= opt.theta_old
        end
    else
        if all(opt.theta_old .!== opt.theta)
            opt.β = min(opt.β * opt.c, 1.0)
            opt.fail = 0
            opt.loss = loss
            opt.theta_old .= opt.theta
        end
        
        opt.m = opt.b1 .* opt.m + (1 - opt.b1) .* gt
        opt.v = opt.b2 .* opt.v + (1 - opt.b2) .* gt .^ 2
        mhat = opt.m ./ (1 - opt.b1^opt.t)
        vhat = opt.v ./ (1 - opt.b2^opt.t)
        opt.theta -= opt.β * opt.a[opt.t] * (mhat ./ (sqrt.(vhat) .+ opt.eps))

        keepInBounds!(opt.theta)

    end

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
