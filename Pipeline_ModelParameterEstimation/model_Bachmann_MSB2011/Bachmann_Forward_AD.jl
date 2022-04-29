using ModelingToolkit, DifferentialEquations, BenchmarkTools, DataFrames, CSV
using JuMP, NLopt, LinearAlgebra, DiffEqSensitivity, ForwardDiff
using ModelingToolkit: varmap_to_vars
using ForwardDiff: GradientConfig, Chunk

println("Done loading modules")

include(joinpath(pwd(), "Additional_functions", "additional_tools.jl"))
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "LatinHyperCubeParameters", "LatinHyperCubeSampledParameters.jl"))

struct ParameterSpace
    scaleLowerBound::Vector{Float64}
    scaleUpperBound::Vector{Float64}
    offsetLowerBound::Vector{Float64}
    offsetUpperBound::Vector{Float64}
    varianceLowerBound::Vector{Float64}
    varianceUpperBound::Vector{Float64}
    parameterLowerBound::Vector{Float64}
    parameterUpperBound::Vector{Float64}
    doLogSearch::Vector{Int64}
end

struct ExperimentalData
    measurementForCondObs::Array{Vector{Float64}, 2}
    timeStepsForCond::Vector{Vector{Float64}}
    observedAtIndexForCondObs::Array{Vector{Int64}, 2}
    numConditions::Int64
    numObservables::Int64
    numDataForCondObs::Array{Int64, 2}
    inputParameterValuesForCond::Vector{Vector{Float64}}
    observedObservableForCond::Vector{Vector{Int64}}
end

struct ModelData
    observableLogTransformation::Vector{Bool}
    optParameterIndices::Vector{Int64}
    inputParameterIndices::Vector{Int64}
    initVariableIndices::Vector{Int64}
    observableVariableIndices::Vector{Int64}
    parameterInObservableIndices::Vector{Int64}
    parameterInU0Indices::Vector{Int64}
end

struct ModelParameters
    scaleIndices::Vector{Int64}
    offsetIndices::Vector{Int64}
    varianceIndices::Vector{Int64}
    parameterIndices::Vector{Int64}
    #
    scaleMap::Array{Int64, 2}
    offsetMap::Array{Int64, 2}
    varianceMap::Array{Int64, 2}
    #
    numScale::Int64
    numOffset::Int64
    numVariance::Int64
    numOptParameters::Int64
    numUsedParameters::Int64
    #
    scaleNames::Vector{String}
    offsetNames::Vector{String}
    varianceNames::Vector{String}
    optParameterNames::Vector{String}
    #
    dualScaleVector::Vector{ForwardDiff.Dual}
    dualOffsetVector::Vector{ForwardDiff.Dual}
    dualVarianceVector::Vector{ForwardDiff.Dual}
    dualDynamicParametersVector::Vector{ForwardDiff.Dual}
    scaleVector::Vector{Float64}
    offsetVector::Vector{Float64}
    varianceVector::Vector{Float64}
    dynamicParametersVector::Vector{Float64}
    allParameters::Vector{Float64}
end

struct ModelOutput
    h_barForCondObs::Array{Vector{ForwardDiff.Dual}, 2}
    h_hatForCondObs::Array{Vector{ForwardDiff.Dual}, 2}
    costForCondObs::Array{ForwardDiff.Dual, 2}
    grad::Vector{Float64}
end



struct FilesAndPaths
    modelFile::String
    modelPath::String
    writefile_solver::String
    writefile_sensealg::String
end

mutable struct Adam
    theta::AbstractArray{Float64} # Parameter array
    grad::Function                # Gradient function
    loss::Float64
    m::AbstractArray{Float64}     # First moment
    v::AbstractArray{Float64}     # Second moment
    b1::Float64                   # Exp. decay first moment
    b2::Float64                   # Exp. decay second moment
    a::Float64                    # Step size
    eps::Float64                  # Epsilon for stability
    t::Int                        # Time step (iteration)
end
  
# Outer constructor
function Adam(theta::AbstractArray{Float64}, grad::Function; a::Float64=0.005)
    loss = 0
    m   = zeros(size(theta))
    v   = zeros(size(theta))
    b1  = 0.9
    b2  = 0.999
    a   = a
    eps = 1e-8
    t   = 0
    Adam(theta, grad, loss, m, v, b1, b2, a, eps, t)
end

# Step function with optional keyword arguments for the data passed to grad()
function step!(opt::Adam, parameterSpace, modelParameters)
    opt.t += 1
    #println(opt.theta)
    gt, opt.loss    = opt.grad(opt.theta)
    #println(opt.theta)
    opt.m = opt.b1 .* opt.m + (1 - opt.b1) .* gt
    opt.v = opt.b2 .* opt.v + (1 - opt.b2) .* gt .^ 2
    mhat = opt.m ./ (1 - opt.b1^opt.t)
    vhat = opt.v ./ (1 - opt.b2^opt.t)
    opt.theta -= opt.a .* (mhat ./ (sqrt.(vhat) .+ opt.eps))

    slb = parameterSpace.scaleLowerBound
    sub = parameterSpace.scaleUpperBound
    scaleIndices = modelParameters.scaleIndices
    scale = view(opt.theta, scaleIndices)
    scale .= min.(max.(scale, slb), sub)
    olb = parameterSpace.offsetLowerBound
    oub = parameterSpace.offsetUpperBound
    offsetIndices = modelParameters.offsetIndices
    offset = view(opt.theta, offsetIndices)
    offset .= min.(max.(offset, olb), oub)
    vlb = parameterSpace.varianceLowerBound
    vub = parameterSpace.varianceUpperBound
    varianceIndices = modelParameters.varianceIndices
    variance  = view(opt.theta, varianceIndices)
    variance .= min.(max.(variance, vlb), vub)
    plb = parameterSpace.parameterLowerBound
    pub = parameterSpace.parameterUpperBound
    parameterIndices = modelParameters.parameterIndices
    pars = view(opt.theta, parameterIndices)
    pars .= min.(max.(pars, plb), pub)
end

function calcUnscaledObservable_proto(prob, modelParameters, modelData, experimentalData, modelOutput, iCond, p)
    parameterIndices = modelParameters.parameterIndices
    dynPars = view(p, parameterIndices)
    dynParVector = modelParameters.dynamicParametersVector
    dualDynParVector = modelParameters.dualDynamicParametersVector
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

    sol = OrdinaryDiffEq.solve(_prob, Rodas5(), reltol=1e-9, abstol=1e-9, saveat = experimentalData.timeStepsForCond[iCond])

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

function calcScaledObservable_proto(modelParameters, modelData, modelOutput, experimentalData, p)
    scaleIndices = modelParameters.scaleIndices
    scale = view(p, scaleIndices)
    scaleVector = modelParameters.scaleVector
    dualScaleVector = modelParameters.dualScaleVector
    dualScaleVector .= convert.(eltype(scale), scaleVector)
    dualScaleVector[1:modelParameters.numScale] .= scale
    scaleMap = modelParameters.scaleMap

    offsetIndices = modelParameters.offsetIndices
    offset = view(p, offsetIndices)
    offsetVector = modelParameters.offsetVector
    dualOffsetVector = modelParameters.dualOffsetVector
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

function calcCost_proto(modelParameters, modelOutput, experimentalData, p)
    varianceIndices = modelParameters.varianceIndices
    variance = view(p, varianceIndices)
    dualVarianceVector = modelParameters.dualVarianceVector
    dualVarianceVector .= variance
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

function allConditionsCost_proto(modelParameters, experimentalData, modelData, calcUnscaledObservable, calcScaledObservable, calcCost, p)
    dynParVector = modelParameters.dynamicParametersVector
    inputPVFC = experimentalData.inputParameterValuesForCond
    inputPI = modelData.inputParameterIndices

    for (iCond, inputParameterValues) in enumerate(inputPVFC)
        dynParVector[inputPI] = inputParameterValues
        calcUnscaledObservable(iCond, p)
    end

    calcScaledObservable(p)

    cost = calcCost(p)

    return cost
end



function f_cost_proto(result::DiffResults.MutableDiffResult, f::Function, cfg::GradientConfig, parameterSpace, modelParameters, modelOutput, p...)::Float64
    doLogSearch = parameterSpace.doLogSearch
    allPar = modelParameters.allParameters
    allPar .= p
    view(allPar, doLogSearch) .= exp10.(view(allPar, doLogSearch))

    ForwardDiff.gradient!(result, f, allPar, cfg)
    modelOutput.grad[:] = DiffResults.gradient(result)
    
    println(DiffResults.value(result), "\n")
    return DiffResults.value(result)
end

function f_grad_proto(grad, modelOutput, p...)
    grad[:] = modelOutput.grad
    nothing
end





function ParameterSpace(modelParameters, parameterBounds; 
        scaleDeterminer = "scale", offsetDeterminer = "offset", varianceDeterminer = "sd_")
    numScale = modelParameters.numScale
    numOffset = modelParameters.numOffset
    numVariance = modelParameters.numVariance
    numOptParameters = modelParameters.numOptParameters

    ## Setting up lower and upper bounds
    slb = Vector{Float64}(undef, numScale)
    sub = Vector{Float64}(undef, numScale)
    isLogScale = Vector{Bool}(undef, numScale) 
    olb = Vector{Float64}(undef, numOffset)
    oub = Vector{Float64}(undef, numOffset)
    isLogOffset = Vector{Bool}(undef, numOffset) 
    vlb = Vector{Float64}(undef, numVariance)
    vub = Vector{Float64}(undef, numVariance)
    isLogVariance = Vector{Bool}(undef, numVariance) 
    plb = Vector{Float64}(undef, numOptParameters)
    pub = Vector{Float64}(undef, numOptParameters) 
    isLogParameter = Vector{Bool}(undef, numOptParameters) 

    scaleNames = modelParameters.scaleNames
    offsetNames = modelParameters.offsetNames
    varianceNames = modelParameters.varianceNames
    optParameterNames = modelParameters.optParameterNames

    for (i, parId) in enumerate(parameterBounds[!, 1])
        if occursin(scaleDeterminer, parId)
            sIndex = findfirst(lowercase(parId) .== lowercase.(scaleNames))
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
            oIndex = findfirst(lowercase(parId) .== lowercase.(offsetNames))
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
            vIndex = findfirst(lowercase(parId) .== lowercase.(varianceNames))
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
            pIndex = findfirst(lowercase(parId) .== lowercase.(optParameterNames))
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

    doLogSearchScale = scaleIndices[isLogScale]
    doLogSearchOffset = offsetIndices[isLogOffset]
    doLogSearchVariance = varianceIndices[isLogVariance]
    doLogSearchParameter = parameterIndices[isLogParameter]
    doLogSearch = vcat(doLogSearchScale, doLogSearchOffset, doLogSearchVariance, doLogSearchParameter)

    parameterSpace = ParameterSpace(slb, sub, olb, oub, vlb, vub, plb, pub, doLogSearch)

    lowerBounds = vcat(slb, olb, vlb, plb)
    upperBounds = vcat(sub, oub, vub, pub)
    numAllStartParameters = length(lowerBounds)

    return parameterSpace, numAllStartParameters, lowerBounds, upperBounds
end

function ModelParameters(new_sys, prob, parameterBounds, experimentalConditions, measurementData, observables, experimentalData; 
    scaleDeterminer = "scale", offsetDeterminer = "offset", varianceDeterminer = "sd_")
    numConditions = experimentalData.numConditions
    numObservables = experimentalData.numObservables
    observableNames = observables[!, 1]

    ## Creating a map from a condition, observable pair to a scale ##
    scaleNames = parameterBounds[findall(occursin.(scaleDeterminer, parameterBounds[!, 1])), 1]
    numScale = length(findall(occursin.(scaleDeterminer, parameterBounds[!, 1])))
    # scale[numScale+1] = 1
    scaleMap = Array{Int64, 2}(undef, numConditions, numObservables) # Have as sparse?
    for (iCond, Cond) in enumerate(experimentalConditions[!, 1])
        for (iObs, Obs) in enumerate(observableNames)
            measurementIndex = findfirst((measurementData[!, 1] .== Obs) .& (measurementData[!, 3] .== Cond))
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
                        scaleIndex = findfirst(scaleNames .== obsPar)
                        scaleMap[iCond, iObs] = scaleIndex
                    else
                        scaleMap[iCond, iObs] = numScale+1
                    end
                end
            end
        end
    end

    ## Creating a map from a condition, observable pair to an offset ##
    offsetNames = parameterBounds[findall(occursin.(offsetDeterminer, parameterBounds[!, 1])), 1]
    numOffset = length(findall(occursin.(offsetDeterminer, parameterBounds[!, 1])))
    # offset[numOffset+1] = 0, offset[numOffset+2] = 1
    offsetMap = Array{Int64, 2}(undef, numConditions, numObservables)
    for (iCond, Cond) in enumerate(experimentalConditions[!, 1])
        for (iObs, Obs) in enumerate(observableNames)
            measurementIndex = findfirst((measurementData[!, 1] .== Obs) .& (measurementData[!, 3] .== Cond))
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
    varianceNames = parameterBounds[findall(occursin.(varianceDeterminer, parameterBounds[!, 1])), 1]
    numVariance = length(findall(occursin.(varianceDeterminer, parameterBounds[!, 1])))
    varianceMap = Array{Int64, 2}(undef, numConditions, numObservables)
    for (iCond, Cond) in enumerate(experimentalConditions[!, 1])
        for (iObs, Obs) in enumerate(observableNames)
            measurementIndex = findfirst((measurementData[!, 1] .== Obs) .& (measurementData[!, 3] .== Cond))
            if measurementIndex === nothing
                varianceMap[iCond, iObs] = numVariance+1
            else
                obsPar = measurementData[measurementIndex, 7]
                if obsPar === missing
                    varianceMap[iCond, iObs] = numVariance+1
                else
                    if occursin(';', obsPar)
                        obsPars = split(obsPar, ';')
                        for op in obsPars
                            if occursin(varianceDeterminer, op)
                                obsPar = op 
                            end
                        end
                    end
                    if occursin(varianceDeterminer, obsPar)
                        varianceIndex = findfirst(varianceNames .== obsPar)
                        varianceMap[iCond, iObs] = varianceIndex
                    else
                        varianceMap[iCond, iObs] = numVariance+1
                    end
                end
            end
        end
    end

    numUsedParameters = length(prob.p)
    problemParameterSymbols = parameters(new_sys)
    inputParameterSymbols = Symbol.(names(experimentalConditions)[2:end])
    tmp = [pPS.name for pPS in problemParameterSymbols]
    inputParameterSymbols = inputParameterSymbols[[any(inpPar .== tmp) for inpPar in inputParameterSymbols]]

    optParameterIndices = collect(1:numUsedParameters)[[pPS.name ∉ inputParameterSymbols for pPS in problemParameterSymbols]]
    numOptParameters = length(optParameterIndices)
    parameterNames = string.(problemParameterSymbols)
    optParameterNames = parameterNames[optParameterIndices]

    scaleIndices = collect(1:numScale)
    offsetIndices = maximum(scaleIndices) .+ collect(1:numOffset)
    varianceIndices = maximum(offsetIndices) .+ collect(1:numVariance)
    parameterIndices = maximum(varianceIndices) .+ collect(1:numOptParameters)

    dualScaleVector = Vector{ForwardDiff.Dual}(undef, numScale + 1)
    dualOffsetVector = Vector{ForwardDiff.Dual}(undef, numOffset + 2)
    dualVarianceVector = Vector{ForwardDiff.Dual}(undef, numVariance)
    dualDynamicParameters = Vector{ForwardDiff.Dual}(undef, numUsedParameters)

    scaleVector = Vector{Float64}(undef, numScale + 1)
    scaleVector[numScale + 1] = 1.0
    offsetVector = Vector{Float64}(undef, numOffset + 2)
    offsetVector[numOffset + 1] = 0.0
    offsetVector[numOffset + 2] = 1.0
    varianceVector = Vector{Float64}(undef, numVariance)
    dynamicParametersVector = Vector{Float64}(undef, numUsedParameters)
    allParameters = Vector{Float64}(undef, numScale + numOffset + numVariance + numOptParameters)

    modelParameters = ModelParameters(scaleIndices, offsetIndices, varianceIndices, parameterIndices, 
        scaleMap, offsetMap, varianceMap, numScale, numOffset, numVariance, numOptParameters, numUsedParameters, 
        scaleNames, offsetNames, varianceNames, optParameterNames, dualScaleVector, dualOffsetVector, dualVarianceVector, dualDynamicParameters, 
        scaleVector, offsetVector, varianceVector, dynamicParametersVector, allParameters)

    return modelParameters
end

function ModelData(new_sys, prob, observables, experimentalConditions)
    observableLogTransformation = [observables[i, :observableTransformation] == "log10" for i in eachindex(observables[!, :observableTransformation])]

    variableNames = string.(states(new_sys))
    numUsedParameters = length(prob.p)
    problemParameterSymbols = parameters(new_sys)
    inputParameterSymbols = Symbol.(names(experimentalConditions)[2:end])
    tmp = [pPS.name for pPS in problemParameterSymbols]
    inputParameterSymbols = inputParameterSymbols[[any(inpPar .== tmp) for inpPar in inputParameterSymbols]]
    inputParameterIndices = [findfirst([pPS.name == inputParameterSymbols[i] for pPS in problemParameterSymbols]) for i in 1:length(inputParameterSymbols)]
    optParameterIndices = collect(1:numUsedParameters)[[pPS.name ∉ inputParameterSymbols for pPS in problemParameterSymbols]]
    parameterNames = string.(problemParameterSymbols)

    initVariableNames = ["EpoRJAK2_CIS", "SHP1", "STAT5", "EpoRJAK2", "SOCS3", "CIS"]
    initVariableNames = initVariableNames .* "(t)"
    initVariableIndices = [findfirst(initVar .== variableNames) for initVar in initVariableNames]

    observableVariableNames = ["CISRNA", "CIS", "SHP1", "SHP1Act", "SOCS3RNA", "SOCS3", "STAT5", 
            "p12EpoRpJAK2", "p1EpoRpJAK2", "p2EpoRpJAK2", "EpoRpJAK2", "pSTAT5"]
    observableVariableNames = observableVariableNames .* "(t)"
    observableVariableIndices = [findfirst(observableVar .== variableNames) for observableVar in observableVariableNames]

    parameterInU0Names = ["init_EpoRJAK2_CIS", "init_SHP1", "init_SHP1_multiplier", "SHP1ProOE", "init_STAT5", 
            "init_EpoRJAK2", "init_SOCS3_multiplier", "SOCS3EqcOE", "SOCS3Eqc", "init_CIS_multiplier", "CISEqc", "CISEqcOE"]
    parameterInU0Indices = [findfirst(parIU0N .== parameterNames) for parIU0N in parameterInU0Names]

    parameterInObservableNames = ["CISRNAEqc", "CISEqc", "SOCS3RNAEqc", "SOCS3Eqc", "init_EpoRJAK2", "init_STAT5", "init_SHP1"]
    parameterInObservableIndices = [findfirst(parION .== parameterNames) for parION in parameterInObservableNames]

    modelData = ModelData(observableLogTransformation, optParameterIndices, inputParameterIndices, 
        initVariableIndices, observableVariableIndices, parameterInObservableIndices, parameterInU0Indices)

    return modelData, inputParameterSymbols
end

function ExperimentalData(observables, experimentalConditions, measurementData, inputParameterSymbols)
    observableNames = observables[!, 1]

    numConditions = length(experimentalConditions[!,1])
    numObservables = length(observableNames) 

    ## Setting up for different conditions
    inputParameterValuesForCond = Vector{Vector{Float64}}(undef, numConditions)

    timeStepsForCond = Vector{Vector{Float64}}(undef, numConditions)
    numTimeStepsForCond = Vector{Int64}(undef, numConditions)
    observedObservableForCond = Vector{Vector{Int64}}(undef, numConditions)

    observedAtForCond = Vector{Vector{Vector{Float64}}}(undef, numConditions)
    observedAtIndexForCondObs = Array{Vector{Int64}, 2}(undef, numConditions, numObservables)
    measurementForCondObs = Array{Vector{Float64}, 2}(undef, numConditions, numObservables)
    numDataForCondObs = Array{Int64, 2}(undef, numConditions, numObservables)

    for (iCond, condId) in enumerate(experimentalConditions[!,1])
        inputParameterValuesForCond[iCond] = [experimentalConditions[iCond, inPar] for inPar in inputParameterSymbols]

        relevantMeasurementData = measurementData[measurementData[:,3] .== condId, :]
        observableIDs = relevantMeasurementData[:,1]
        timeStepsForCond[iCond] = sort(unique(relevantMeasurementData[:,5]))
        numTimeStepsForCond[iCond] = length(timeStepsForCond[iCond])
        isObserved = [observableNames[i] in unique(observableIDs) for i=1:numObservables]
        observedObservableForCond[iCond] = collect(1:numObservables)[isObserved]

        observedAtForCond[iCond] = Vector{Vector{Float64}}(undef, numObservables)
        for (iObs, obsId) = enumerate(observableNames)
            observedAtForCond[iCond][iObs] = relevantMeasurementData[observableIDs .== obsId, 5]
        end
        for iObs = 1:numObservables
            observedAtIndexForCondObs[iCond, iObs] = indexin(observedAtForCond[iCond][iObs], timeStepsForCond[iCond])
        end
        for (iObs, obsId) = enumerate(observableNames)
            measurementForCondObs[iCond, iObs] = relevantMeasurementData[observableIDs .== obsId, 4]
            numDataForCondObs[iCond, iObs] = length(measurementForCondObs[iCond, iObs])
        end
    end

    experimentalData = ExperimentalData(measurementForCondObs, timeStepsForCond, observedAtIndexForCondObs, numConditions, numObservables, numDataForCondObs, inputParameterValuesForCond, observedObservableForCond)

    return experimentalData
end

function ModelOutput(usedType, experimentalData, modelParameters)
    numConditions = experimentalData.numConditions
    numObservables = experimentalData.numObservables

    h_barForCondObs = Array{Vector{usedType}, 2}(undef, numConditions, numObservables)
    h_hatForCondObs = Array{Vector{usedType}, 2}(undef, numConditions, numObservables)
    costForCondObs = zeros(usedType, numConditions, numObservables)
    grad = Vector{Float64}(undef, length(modelParameters.allParameters))

    modelOutput = ModelOutput(h_barForCondObs, h_hatForCondObs, costForCondObs, grad)

    return modelOutput
end





function GradCalc_forwardDiff(filesAndPaths, timeEnd, experimentalConditions, measurementData, observables, parameterBounds)
    include(filesAndPaths.modelPath)
    new_sys = ode_order_lowering(sys)
    u0 = initialSpeciesValues
    p = trueParameterValues 
    c = trueConstantsValues
    pars = [p;c]
    tspan = (0.0, timeEnd)
    prob = ODEProblem(new_sys,u0,tspan,pars,jac=true)

    # Initialize structs

    modelData, inputParameterSymbols = ModelData(new_sys, prob, observables, experimentalConditions)

    experimentalData = ExperimentalData(observables, experimentalConditions, measurementData, inputParameterSymbols)

    modelParameters = ModelParameters(new_sys, prob, parameterBounds, experimentalConditions, measurementData, observables, experimentalData)

    parameterSpace, numAllStartParameters, lowerBounds, upperBounds = ParameterSpace(modelParameters, parameterBounds)

    usedType = ForwardDiff.Dual
    modelOutput = ModelOutput(usedType, experimentalData, modelParameters)

    # Initialize functions

    calcUnscaledObservable = (iCond, p) -> calcUnscaledObservable_proto(prob, modelParameters, modelData, experimentalData, modelOutput, iCond, p)

    calcScaledObservable = (p) -> calcScaledObservable_proto(modelParameters, modelData, modelOutput, experimentalData, p)

    calcCost = (p) -> calcCost_proto(modelParameters, modelOutput, experimentalData, p)

    allConditionsCost = (p) -> allConditionsCost_proto(modelParameters, experimentalData, modelData,
            calcUnscaledObservable, calcScaledObservable, calcCost, p)

    allStartParameters = getSamples(numAllStartParameters, lowerBounds, upperBounds, 1)

    chunkSize = 23
    cfg = GradientConfig(allConditionsCost, allStartParameters, Chunk{chunkSize}())

    result = DiffResults.GradientResult(allStartParameters::Vector{Float64})

    f = (p_tuple...) -> f_cost_proto(result, allConditionsCost, cfg, parameterSpace, modelParameters, modelOutput, p_tuple...)
    f_grad = (grad, p_tuple...) -> f_grad_proto(grad, modelOutput, p_tuple...)

    
    model = Model(NLopt.Optimizer)
    #JuMP.register(model::Model, s::Symbol, dimension::Integer, f::Function, ∇f::Function, ∇²f::Function)
    register(model, :f, numAllStartParameters, f, f_grad)
    set_optimizer_attribute(model, "algorithm", :LD_MMA)
    @variable(model, lowerBounds[i] <= p[i = 1:numAllStartParameters] <= upperBounds[i]) # fix lower and upper bounds
    for i in 1:numAllStartParameters
        set_start_value(p[i], allStartParameters[i])
    end
    @NLobjective(model, Min, f(p...))
    JuMP.optimize!(model)
    @show termination_status(model)
    @show primal_status(model)
    println("Done!")
    p_opt = [value(p[i]) for i=1:numAllStartParameters]
    cost_opt = objective_value(model)
    println("Final cost: ", cost_opt)

    optimizedScales = exp10.(p_opt[scaleIndices])
    optimizedOffsets = exp10.(p_opt[offsetIndices])
    optimizedVariances = exp10.(p_opt[varianceIndices])
    optimizedParameters = exp10.(p_opt[parameterIndices])
    println(optimizedParameters)
    
end


















function main()
    
    readModelPath = joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels")
    modelFile = "model_Bachmann_MSB2011.jl" #"model_Alkan_SciSignal2018.jl"
    modelPath = joinpath(readModelPath, modelFile)
    writePath = joinpath(pwd(), "Pipeline_ModelParameterEstimation", "IntermediaryResults")
    fixDirectories(writePath)
    writefile_solver = joinpath(writePath, "benchmark_solver_" * string(getNumberOfFiles(writePath) + 1) * ".csv")
    writefile_sensealg = joinpath(writePath, "benchmark_sensealg_" * string(getNumberOfFiles(writePath) + 1) * ".csv")

    filesAndPaths = FilesAndPaths(modelFile, modelPath, writefile_solver, writefile_sensealg)

    readDataPath = joinpath(pwd(), "Pipeline_ModelParameterEstimation", "Data", modelFile[1:end-3])
    dataEnding = modelFile[7:end-3] * ".tsv"
    experimentalConditions = CSV.read(joinpath(readDataPath, "experimentalCondition_" * dataEnding), DataFrame)
    measurementData = CSV.read(joinpath(readDataPath, "measurementData_" * dataEnding), DataFrame)
    timeEnd = maximum(measurementData[!, 5])
    observables = CSV.read(joinpath(readDataPath, "observables_" * dataEnding), DataFrame)
    parameterBounds = CSV.read(joinpath(readDataPath, "parameters_" * dataEnding), DataFrame)

    println("Starting optimisation")
    GradCalc_forwardDiff(filesAndPaths, timeEnd, experimentalConditions, measurementData, observables, parameterBounds)

end

main()