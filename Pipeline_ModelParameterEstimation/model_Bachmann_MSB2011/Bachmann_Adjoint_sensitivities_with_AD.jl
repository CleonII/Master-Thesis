using ModelingToolkit, DifferentialEquations, BenchmarkTools, DataFrames, CSV
using JuMP, NLopt, LinearAlgebra, DiffEqSensitivity, ForwardDiff
using ModelingToolkit: varmap_to_vars

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
    numTimeStepsForCondObs::Array{Int64, 2}
    inputParameterValuesForCond::Vector{Vector{Float64}}
    observedObservableForCond::Vector{Vector{Int64}}
    observablesTimeIndexIndicesForCond::Vector{Array{Vector{Int64}, 2}}
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
    scaleVector::Vector{Float64}
    offsetVector::Vector{Float64}
    varianceVector::Vector{Float64}
    dynamicParametersVector::Vector{Float64}
    allParameters::Vector{Float64}
end

struct ModelOutput
    sols::Vector{ODESolution}
    h_barForCondObs::Array{Vector{Float64}, 2}
    h_hatForCondObs::Array{Vector{Float64}, 2}
    costForCondObs::Array{Float64, 2}
    scaleGrad::Vector{Float64}
    offsetGrad::Vector{Float64}
    varianceGrad::Vector{Float64}
    tmp_grad::Vector{Float64}
    gradMatrix::Array{Float64, 2}
end

struct FilesAndPaths
    modelFile::String
    modelPath::String
    writefile_solver::String
    writefile_sensealg::String
end




function solveODESystem_proto(prob, modelParameters, modelData, modelOutput, iCond)
    dynParVector = modelParameters.dynamicParametersVector

    U0Vector = prob.u0
    initVariable = modelData.initVariableIndices
    parameterInU0Indices = modelData.parameterInU0Indices
    U0Vector[initVariable[1]] = dynParVector[parameterInU0Indices[1]]
    U0Vector[initVariable[2]] = dynParVector[parameterInU0Indices[2]] * (dynParVector[parameterInU0Indices[3]] * dynParVector[parameterInU0Indices[4]] + 1)
    U0Vector[initVariable[3]] = dynParVector[parameterInU0Indices[5]]
    U0Vector[initVariable[4]] = dynParVector[parameterInU0Indices[6]]
    U0Vector[initVariable[5]] = dynParVector[parameterInU0Indices[7]] * dynParVector[parameterInU0Indices[8]] * dynParVector[parameterInU0Indices[9]]
    U0Vector[initVariable[6]] = dynParVector[parameterInU0Indices[10]] * dynParVector[parameterInU0Indices[11]] * dynParVector[parameterInU0Indices[12]]

    _prob = remake(prob, u0 = U0Vector, p = dynParVector)
    modelOutput.sols[iCond] = OrdinaryDiffEq.solve(_prob, Rodas5(), reltol=1e-12, abstol=1e-12)
    
    nothing
end



function g_unscaledObservables_proto(type, u, dynPar, i, iCond, modelData, experimentalData)

    optParIOI = modelData.parameterInObservableIndices
    observableVI = modelData.observableVariableIndices
    observablesTIIFC = experimentalData.observablesTimeIndexIndicesForCond[iCond]
    
    h_bar = Vector{Vector{type}}(undef, experimentalData.numObservables)
    
    h_bar[1] = [ u[observableVI[1]] / dynPar[optParIOI[1]] ][ones(Int64, length(observablesTIIFC[1, i]))]             # log10
    h_bar[2] = [ u[observableVI[1]] / dynPar[optParIOI[1]] ][ones(Int64, length(observablesTIIFC[2, i]))]                                            # log10
    h_bar[3] = [ u[observableVI[1]] / dynPar[optParIOI[1]] ][ones(Int64, length(observablesTIIFC[3, i]))]                                            # log10
    h_bar[4] = [ u[observableVI[2]] ][ones(Int64, length(observablesTIIFC[4, i]))]                                                                             # log10
    h_bar[5] = [ u[observableVI[2]] / dynPar[optParIOI[2]] ][ones(Int64, length(observablesTIIFC[5, i]))]                                               # log10
    h_bar[6] = [ u[observableVI[2]] / dynPar[optParIOI[2]] ][ones(Int64, length(observablesTIIFC[6, i]))]                                               # log10
    h_bar[7] = [ u[observableVI[2]] / dynPar[optParIOI[2]] ][ones(Int64, length(observablesTIIFC[7, i]))]                                               # log10
    h_bar[8] = [ u[observableVI[3]] + dynPar[observableVI[4]] ][ones(Int64, length(observablesTIIFC[8, i]))]                                                # log10
    h_bar[9] = [ u[observableVI[5]] / dynPar[optParIOI[3]] ][ones(Int64, length(observablesTIIFC[9, i]))]                                          # log10
    h_bar[10] = [ u[observableVI[5]] / dynPar[optParIOI[3]] ][ones(Int64, length(observablesTIIFC[10, i]))]                                        # log10
    h_bar[11] = [ u[observableVI[5]] / dynPar[optParIOI[3]] ][ones(Int64, length(observablesTIIFC[11, i]))]                                        # log10
    h_bar[12] = [ u[observableVI[6]] ][ones(Int64, length(observablesTIIFC[12, i]))]                                                                         # log10
    h_bar[13] = [ u[observableVI[6]] / dynPar[optParIOI[4]] ][ones(Int64, length(observablesTIIFC[13, i]))]                                           # log10
    h_bar[14] = [ u[observableVI[7]] ][ones(Int64, length(observablesTIIFC[14, i]))]                                                                         # log10
    h_bar[15] = [ 16 * (u[observableVI[8]] + u[observableVI[9]] + 
            u[observableVI[10]]) / dynPar[optParIOI[5]] ][ones(Int64, length(observablesTIIFC[15, i]))]                                                     # log10
    h_bar[16] = [ 2 * (u[observableVI[11]] + u[observableVI[8]] + 
            u[observableVI[9]] + u[observableVI[10]]) / dynPar[optParIOI[5]] ][ones(Int64, length(observablesTIIFC[16, i]))]                     # log10
    h_bar[17] = [ (100 * u[observableVI[12]]) / (u[observableVI[7]] + u[observableVI[12]]) ][ones(Int64, length(observablesTIIFC[17, i]))]        # lin
    h_bar[18] = [ u[observableVI[12]] / dynPar[optParIOI[6]] ][ones(Int64, length(observablesTIIFC[18, i]))]                                           # log10
    h_bar[19] = [ (u[observableVI[3]] + u[observableVI[4]]) / dynPar[optParIOI[7]] ][ones(Int64, length(observablesTIIFC[19, i]))]               # log10
    h_bar[20] = [ (u[observableVI[7]] + u[observableVI[12]]) / dynPar[optParIOI[6]] ][ones(Int64, length(observablesTIIFC[20, i]))]               # log10

    return h_bar
end

function g_scaledObservationFunctions_proto(type, h_bar, scaleVector, offsetVector, iCond, modelParameters, experimentalData, modelData)
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

function g_cost_proto(type, h_hat, varianceVector, i, iCond, modelParameters, experimentalData)
    varianceMap = modelParameters.varianceMap
    measurementFCO = experimentalData.measurementForCondObs
    observablesTIIFC = experimentalData.observablesTimeIndexIndicesForCond[iCond]
    observedOFC = experimentalData.observedObservableForCond

    costFCO = Vector{type}(undef, experimentalData.numObservables)
    for iObs in observedOFC[iCond]
        costFCO[iObs] = log(2*pi*varianceVector[varianceMap[iCond, iObs]]) * length(h_hat[iObs]) + 
                (dot(measurementFCO[iCond, iObs][observablesTIIFC[iObs, i]], measurementFCO[iCond, iObs][observablesTIIFC[iObs, i]]) - 
                2*dot(measurementFCO[iCond, iObs][observablesTIIFC[iObs, i]], h_hat[iObs]) + 
                dot(h_hat[iObs], h_hat[iObs])) / (2 * varianceVector[varianceMap[iCond, iObs]])
    end
    return sum(costFCO[observedOFC[iCond]])
end

function get_type(lists)
    type = Float64
    if eltype(lists) <: AbstractVector
        for list in lists
            if eltype(list) !== Float64
                type = eltype(list)
                break
            end
        end
    else 
        type = eltype(lists)
    end
    return type
end

function g_proto(u, dynPar, scale, offset, variance, i, iCond, g_unscaledObservables, g_scaledObservationFunctions, g_cost; type = get_type([u, dynPar, scale, offset, variance]))

    h_bar = g_unscaledObservables(type, u, dynPar, i, iCond)

    h_hat = g_scaledObservationFunctions(type, h_bar, scale, offset, iCond)

    cost = g_cost(type, h_hat, variance, i, iCond)

    return cost
end

function G_proto(solveODESystem, g, iCond, modelParameters, experimentalData, modelOutput)
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

function updateAllParameterVectors_proto(modelParameters, modelData)
    allParameters = modelParameters.allParameters
    parameterIndices = modelParameters.parameterIndices
    pars = view(allParameters, parameterIndices)
    dynPar = modelParameters.dynamicParametersVector
    dynPar[modelData.optParameterIndices] .= pars

    scaleIndices = modelParameters.scaleIndices
    scale = view(allParameters, scaleIndices)
    scaleVector = modelParameters.scaleVector
    scaleVector[1:modelParameters.numScale] .= scale

    offsetIndices = modelParameters.offsetIndices
    offset = view(allParameters, offsetIndices)
    offsetVector = modelParameters.offsetVector
    offsetVector[1:modelParameters.numOffset] .= offset

    varianceIndices = modelParameters.varianceIndices
    variance = view(allParameters, varianceIndices)
    varianceVector = modelParameters.varianceVector
    varianceVector[1:modelParameters.numVariance] .= variance
end

function allConditionsCost_proto(parameterSpace, modelParameters, experimentalData, modelData, 
    updateAllParameterVectors, G, p...)

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

    println("Cost: ", cost)

    return cost
end

function dg_proto!(out, u, p, t, i, gu) 

    guAt = (u) -> gu(u, i)
    out .= ForwardDiff.gradient(guAt, u)

    nothing
end

function calcCostGrad_proto(g, dg!, iCond, modelParameters, modelData, experimentalData, modelOutput)

    timeSteps = experimentalData.timeStepsForCond[iCond]
    tmp_grad = modelOutput.tmp_grad

    sol = modelOutput.sols[iCond]

    ~, tmp_grad[:] = adjoint_sensitivities(sol, Rodas4P(), dg!, timeSteps, 
            sensalg = QuadratureAdjoint(autojacvec = ReverseDiffVJP(true)), reltol = 1e-12, abstol = 1e-12)

    # when a parameter is included in the observation function the gradient is incorrect, has to correct with adding dgdp 


    # gradient for scales, offsets and variances 
    dynPar = modelParameters.dynamicParametersVector
    dgddynPar = zeros(length(dynPar))
    scale = modelParameters.scaleVector
    dgdscale = zeros(length(scale))
    offset = modelParameters.offsetVector
    dgdoffset = zeros(length(offset))
    variance = modelParameters.varianceVector
    dgdvariance = zeros(length(variance))

    gdynPar = (u, dynPar, i) -> g(u, dynPar, scale, offset, variance, i, iCond)
    gscale = (u, scale, i) -> g(u, dynPar, scale, offset, variance, i, iCond)
    goffset = (u, offset, i) -> g(u, dynPar, scale, offset, variance, i, iCond)
    gvariance = (u, variance, i) -> g(u, dynPar, scale, offset, variance, i, iCond)
    for (i, t) in enumerate(timeSteps)
        gdynParAt = (dynPar) -> gdynPar(sol(t), dynPar, i)
        dgddynPar += ForwardDiff.gradient(gdynParAt, dynPar)

        gscaleAt = (scale) -> gscale(sol(t), scale, i)
        dgdscale += ForwardDiff.gradient(gscaleAt, scale)

        goffsetAt = (offset) -> goffset(sol(t), offset, i)
        dgdoffset += ForwardDiff.gradient(goffsetAt, offset)

        gvarianceAt = (variance) -> gvariance(sol(t), variance, i)
        dgdvariance += ForwardDiff.gradient(gvarianceAt, variance)
    end


    # Insert into common gradient
    gradMatrix = modelOutput.gradMatrix

    scaleIndices = modelParameters.scaleIndices
    gradMatrix[iCond, scaleIndices] = view(dgdscale, 1:modelParameters.numScale)
    offsetIndices = modelParameters.offsetIndices
    gradMatrix[iCond, offsetIndices] = view(dgdoffset, 1:modelParameters.numOffset)
    varianceIndices = modelParameters.varianceIndices
    gradMatrix[iCond, varianceIndices] = view(dgdvariance, 1:modelParameters.numVariance)

    parameterIndices = modelParameters.parameterIndices
    optParameterIndices = modelData.optParameterIndices
    gradMatrix[iCond, parameterIndices] = (tmp_grad[optParameterIndices] .* -1) + dgddynPar[optParameterIndices]

    nothing
end


function allConditionsCostGrad_proto(parameterSpace, modelParameters, modelData, modelOutput, experimentalData, 
    updateAllParameterVectors, calcCostGrad, g, grad, p...)
    println("In f_prime()")

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
        dg! = (out, u, p, t, i) -> dg_proto!(out, u, p, t, i, gu)
        calcCostGrad(dg!, iCond)
    end

    replace!(gradMatrix, NaN=>0.0)

    grad[:] = vec(sum(gradMatrix, dims=1))# / (norm(vec(sum(gradMatrix, dims=1))) + eps(Float64))
    println("Grad")
    println(grad, "\n")

    nothing
end
















#=
function calcUnscaledObservable_proto(prob, modelParameters, modelData, experimentalData, modelOutput, iCond)
    dynParVector = modelParameters.dynamicParametersVector

    U0Vector = prob.u0
    initVariable = modelData.initVariableIndices
    parameterInU0Indices = modelData.parameterInU0Indices
    U0Vector[initVariable[1]] = dynParVector[parameterInU0Indices[1]]
    U0Vector[initVariable[2]] = dynParVector[parameterInU0Indices[2]] * (dynParVector[parameterInU0Indices[3]] * dynParVector[parameterInU0Indices[4]] + 1)
    U0Vector[initVariable[3]] = dynParVector[parameterInU0Indices[5]]
    U0Vector[initVariable[4]] = dynParVector[parameterInU0Indices[6]]
    U0Vector[initVariable[5]] = dynParVector[parameterInU0Indices[7]] * dynParVector[parameterInU0Indices[8]] * dynParVector[parameterInU0Indices[9]]
    U0Vector[initVariable[6]] = dynParVector[parameterInU0Indices[10]] * dynParVector[parameterInU0Indices[11]] * dynParVector[parameterInU0Indices[12]]

    _prob = remake(prob, u0 = U0Vector, p = dynParVector)
    sol = OrdinaryDiffEq.solve(_prob, Rodas5(), reltol=1e-9, abstol=1e-9, saveat = experimentalData.timeStepsForCond[iCond])

    h_barFCO = modelOutput.h_barForCondObs
    observedAIFCO = experimentalData.observedAtIndexForCondObs

    optParIOI = modelData.parameterInObservableIndices

    h_barFCO[iCond, 1] = ( sol[Symbol("CISRNA(t)")] ./ dynParVector[optParIOI[1]] )[observedAIFCO[iCond, 1]]                                            # log10
    h_barFCO[iCond, 2] = ( sol[Symbol("CISRNA(t)")] ./ dynParVector[optParIOI[1]] )[observedAIFCO[iCond, 2]]                                            # log10
    h_barFCO[iCond, 3] = ( sol[Symbol("CISRNA(t)")] ./ dynParVector[optParIOI[1]] )[observedAIFCO[iCond, 3]]                                            # log10
    h_barFCO[iCond, 4] = ( sol[Symbol("CIS(t)")] )[observedAIFCO[iCond, 4]]                                                                             # log10
    h_barFCO[iCond, 5] = ( sol[Symbol("CIS(t)")] ./ dynParVector[optParIOI[2]] )[observedAIFCO[iCond, 5]]                                               # log10
    h_barFCO[iCond, 6] = ( sol[Symbol("CIS(t)")] ./ dynParVector[optParIOI[2]] )[observedAIFCO[iCond, 6]]                                               # log10
    h_barFCO[iCond, 7] = ( sol[Symbol("CIS(t)")] ./ dynParVector[optParIOI[2]] )[observedAIFCO[iCond, 7]]                                               # log10
    h_barFCO[iCond, 8] = ( sol[Symbol("SHP1(t)")] + sol[Symbol("SHP1Act(t)")] )[observedAIFCO[iCond, 8]]                                                # log10
    h_barFCO[iCond, 9] = ( sol[Symbol("SOCS3RNA(t)")] ./ dynParVector[optParIOI[3]] )[observedAIFCO[iCond, 9]]                                          # log10
    h_barFCO[iCond, 10] = ( sol[Symbol("SOCS3RNA(t)")] ./ dynParVector[optParIOI[3]] )[observedAIFCO[iCond, 10]]                                        # log10
    h_barFCO[iCond, 11] = ( sol[Symbol("SOCS3RNA(t)")] ./ dynParVector[optParIOI[3]] )[observedAIFCO[iCond, 11]]                                        # log10
    h_barFCO[iCond, 12] = ( sol[Symbol("SOCS3(t)")] )[observedAIFCO[iCond, 12]]                                                                         # log10
    h_barFCO[iCond, 13] = ( sol[Symbol("SOCS3(t)")] ./ dynParVector[optParIOI[4]] )[observedAIFCO[iCond, 13]]                                           # log10
    h_barFCO[iCond, 14] = ( sol[Symbol("STAT5(t)")] )[observedAIFCO[iCond, 14]]                                                                         # log10
    h_barFCO[iCond, 15] = ( 16 * (sol[Symbol("p12EpoRpJAK2(t)")] + sol[Symbol("p1EpoRpJAK2(t)")] + 
            sol[Symbol("p2EpoRpJAK2(t)")]) / dynParVector[optParIOI[5]] )[observedAIFCO[iCond, 15]]                                                     # log10
    h_barFCO[iCond, 16] = ( 2 * (sol[Symbol("EpoRpJAK2(t)")] + sol[Symbol("p12EpoRpJAK2(t)")] + 
            sol[Symbol("p1EpoRpJAK2(t)")] + sol[Symbol("p2EpoRpJAK2(t)")]) / dynParVector[optParIOI[5]] )[observedAIFCO[iCond, 16]]                     # log10
    h_barFCO[iCond, 17] = ( (100 * sol[Symbol("pSTAT5(t)")]) ./ (sol[Symbol("STAT5(t)")] + sol[Symbol("pSTAT5(t)")]) )[observedAIFCO[iCond, 17]]        # lin
    h_barFCO[iCond, 18] = ( sol[Symbol("pSTAT5(t)")] / dynParVector[optParIOI[6]] )[observedAIFCO[iCond, 18]]                                           # log10
    h_barFCO[iCond, 19] = ( (sol[Symbol("SHP1(t)")] + sol[Symbol("SHP1Act(t)")]) / dynParVector[optParIOI[7]] )[observedAIFCO[iCond, 19]]               # log10
    h_barFCO[iCond, 20] = ( (sol[Symbol("STAT5(t)")] + sol[Symbol("pSTAT5(t)")]) / dynParVector[optParIOI[6]] )[observedAIFCO[iCond, 20]]               # log10

    nothing
end

function calcScaledObservable_proto(modelParameters, modelData, modelOutput, experimentalData, iCond)
    
    scaleVector = modelParameters.scaleVector
    scaleMap = modelParameters.scaleMap

    offsetVector = modelParameters.offsetVector
    offsetMap = modelParameters.offsetMap

    observedOFC = experimentalData.observedObservableForCond
    
    h_barFCO = modelOutput.h_barForCondObs
    h_hatFCO = modelOutput.h_hatForCondObs
    linearH_hatFCO = modelOutput.linearH_hatForCondObs
    observableLT = modelData.observableLogTransformation

    for iObs = observedOFC[iCond]
        linearH_hatFCO[iCond, iObs] = scaleVector[scaleMap[iCond, iObs]] * h_barFCO[iCond, iObs] .+ offsetVector[offsetMap[iCond, iObs]]
        if observableLT[iObs]
            h_hatFCO[iCond, iObs] = log10.(scaleVector[scaleMap[iCond, iObs]] * h_barFCO[iCond, iObs] .+ offsetVector[offsetMap[iCond, iObs]])
        else
            h_hatFCO[iCond, iObs] = scaleVector[scaleMap[iCond, iObs]] * h_barFCO[iCond, iObs] .+ offsetVector[offsetMap[iCond, iObs]]
        end
    end

    nothing
end

function calcCost_proto(modelParameters, modelOutput, experimentalData, iCond)
    varianceVector = modelParameters.varianceVector
    varianceMap = modelParameters.varianceMap

    h_hatFCO = modelOutput.h_hatForCondObs
    
    costFCO = modelOutput.costForCondObs
    measurementFCO = experimentalData.measurementForCondObs
    observedOFC = experimentalData.observedObservableForCond

    for iObs in observedOFC[iCond]
        costFCO[iCond, iObs] = log(2*pi*varianceVector[varianceMap[iCond, iObs]]) + (dot(measurementFCO[iCond, iObs], measurementFCO[iCond, iObs]) - 
                2*dot(measurementFCO[iCond, iObs], h_hatFCO[iCond, iObs]) + 
                dot(h_hatFCO[iCond, iObs], h_hatFCO[iCond, iObs])) / (2 * varianceVector[varianceMap[iCond, iObs]])
    end

    return sum(costFCO[iCond, observedOFC[iCond]])
end

function updateAllParameterVectors_proto(modelParameters, modelData)
    allParameters = modelParameters.allParameters
    parameterIndices = modelParameters.parameterIndices
    pars = view(allParameters, parameterIndices)
    dynPar = modelParameters.dynamicParametersVector
    dynPar[modelData.optParameterIndices] .= pars

    scaleIndices = modelParameters.scaleIndices
    scale = view(allParameters, scaleIndices)
    scaleVector = modelParameters.scaleVector
    scaleVector[1:modelParameters.numScale] .= scale

    offsetIndices = modelParameters.offsetIndices
    offset = view(allParameters, offsetIndices)
    offsetVector = modelParameters.offsetVector
    offsetVector[1:modelParameters.numOffset] .= offset

    varianceIndices = modelParameters.varianceIndices
    variance = view(allParameters, varianceIndices)
    varianceVector = modelParameters.varianceVector
    varianceVector[1:modelParameters.numVariance] .= variance
end

function allConditionsCost_proto(parameterSpace, modelParameters, experimentalData, modelData, 
            updateAllParameterVectors, calcUnscaledObservable, calcScaledObservable, calcCost, p...)
    println("In f()")
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
        calcUnscaledObservable(iCond)
        calcScaledObservable(iCond)
        cost += calcCost(iCond)
    end


    println("Cost: ", cost)

    return cost
end


# the total cost, G, is the integral/sum of the costs at each time step, g
# here we have the derivative, dg, of the cost at each time step
# Note: adjoint_sensitivities seems to use the negative gradient due to a sign error
# Each out[i] corresponds with the gradient of the cost with regards to state i at time t
function dg_proto!(out,u,p,t,i, g, iCond) 

    gu = (u) -> g(u, p, t, i, iCond)
    out .= ForwardDiff.gradient(gu, u)

    gp = (p) -> g(u, p, t, i, iCond)
    dgdp += ForwardDiff.gradient(gp, p) 

    nothing
end


function calcCostGrad_proto(prob, modelParameters, modelData, experimentalData, modelOutput, dg!, iCond)
    dynParVector = modelParameters.dynamicParametersVector

    U0Vector = prob.u0
    initVariable = modelData.initVariableIndices
    parameterInU0Indices = modelData.parameterInU0Indices
    U0Vector[initVariable[1]] = dynParVector[parameterInU0Indices[1]]
    U0Vector[initVariable[2]] = dynParVector[parameterInU0Indices[2]] * (dynParVector[parameterInU0Indices[3]] * dynParVector[parameterInU0Indices[4]] + 1)
    U0Vector[initVariable[3]] = dynParVector[parameterInU0Indices[5]]
    U0Vector[initVariable[4]] = dynParVector[parameterInU0Indices[6]]
    U0Vector[initVariable[5]] = dynParVector[parameterInU0Indices[7]] * dynParVector[parameterInU0Indices[8]] * dynParVector[parameterInU0Indices[9]]
    U0Vector[initVariable[6]] = dynParVector[parameterInU0Indices[10]] * dynParVector[parameterInU0Indices[11]] * dynParVector[parameterInU0Indices[12]]

    _prob = remake(prob, u0 = U0Vector, p = dynParVector)

    timeSteps = experimentalData.timeStepsForCond[iCond]
    tmp_grad = modelOutput.tmp_grad

    sol = OrdinaryDiffEq.solve(_prob, Rodas5(), reltol=1e-9, abstol=1e-9, 
            sensealg=InterpolatingAdjoint(checkpointing=true))

    ~, tmp_grad[:] = adjoint_sensitivities(sol, Rodas4P(), dg!, timeSteps, 
            sensalg = QuadratureAdjoint(autojacvec = ReverseDiffVJP(true)))

    # when a parameter is included in the observation function the gradient is incorrect, has to correct with adding dgdp 


    # gradient for scales, offsets and variances 
    scaleGrad = modelOutput.scaleGrad
    scaleGrad .= 0.0
    offsetGrad = modelOutput.offsetGrad
    offsetGrad .= 0.0
    varianceGrad = modelOutput.varianceGrad
    varianceGrad .= 0.0

    scaleMap = modelParameters.scaleMap
    offsetMap = modelParameters.offsetMap
    varianceVector = modelParameters.varianceVector
    varianceMap = modelParameters.varianceMap

    h_hatFCO = modelOutput.h_hatForCondObs
    linearH_hatFCO = modelOutput.linearH_hatForCondObs
    h_barFCO = modelOutput.h_barForCondObs
    measurementFCO = experimentalData.measurementForCondObs
    observedO = experimentalData.observedObservableForCond[iCond]
    observableLT = modelData.observableLogTransformation
    numTimeStepsFor = experimentalData.numTimeStepsForCondObs


    for iObs in observedO
        tmp = h_hatFCO[iCond, iObs] .- measurementFCO[iCond, iObs]
        if observableLT[iObs]
            scaleGrad[scaleMap[iCond, iObs]] += dot(h_barFCO[iCond, iObs] ./ linearH_hatFCO[iCond, iObs], tmp) / (varianceVector[varianceMap[iCond, iObs]] * log(10))
            offsetGrad[offsetMap[iCond, iObs]] += dot(1/linearH_hatFCO[iCond, iObs], tmp) / (varianceVector[varianceMap[iCond, iObs]] * log(10))
            varianceGrad[varianceMap[iCond, iObs]] += (0.5/varianceVector[varianceMap[iCond, iObs]]) * (numTimeStepsFor[iCond, iObs] - 
                    (1/varianceVector[varianceMap[iCond, iObs]]) * dot(tmp, tmp))
        else
            scaleGrad[scaleMap[iCond, iObs]] += dot(h_barFCO[iCond, iObs], tmp) / varianceVector[varianceMap[iCond, iObs]]
            offsetGrad[offsetMap[iCond, iObs]] += sum(tmp) / varianceVector[varianceMap[iCond, iObs]]
            varianceGrad[varianceMap[iCond, iObs]] += (0.5/varianceVector[varianceMap[iCond, iObs]]) * (numTimeStepsFor[iCond, iObs] - 
                    (1/varianceVector[varianceMap[iCond, iObs]]) * dot(tmp, tmp))
        end
    end

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
    gradMatrix[iCond, parameterIndices] = tmp_grad[optParameterIndices] .* -1

    nothing
end

function allConditionsCostGrad_proto(parameterSpace, modelParameters, modelData, modelOutput, experimentalData, 
            updateAllParameterVectors, calcCostGrad, grad, p...)
    println("In f_prime()")

    doLogSearch = parameterSpace.doLogSearch
    allParameters = modelParameters.allParameters
    allParameters .= p 
    view(allParameters, doLogSearch) .= exp10.(view(allParameters, doLogSearch))

    updateAllParameterVectors()

    dynParVector = modelParameters.dynamicParametersVector
    inputPVFC = experimentalData.inputParameterValuesForCond
    inputPI = modelData.inputParameterIndices
    gradMatrix = modelOutput.gradMatrix

    for (iCond, inputParameterValues) in enumerate(inputPVFC)
        dynParVector[inputPI] = inputParameterValues
        dg! = (out, u, p, t, i) -> dg_proto!(out, u, p, t, i, parameterSpace, modelParameters, modelData, modelOutput, experimentalData, iCond)
        calcCostGrad(dg!, iCond)
    end

    replace!(gradMatrix, NaN=>0.0)
    
    grad[:] = vec(sum(gradMatrix, dims=1))# / (norm(vec(sum(gradMatrix, dims=1))) + eps(Float64))
    println("Grad")
    println(grad, "\n")

    nothing
end
=#





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

    scaleVector = Vector{Float64}(undef, numScale + 1)
    scaleVector[numScale + 1] = 1.0
    offsetVector = Vector{Float64}(undef, numOffset + 2)
    offsetVector[numOffset + 1] = 0.0
    offsetVector[numOffset + 2] = 1.0
    varianceVector = Vector{Float64}(undef, numVariance + 1)
    varianceVector[numVariance + 1] = 1.0
    dynamicParametersVector = Vector{Float64}(undef, numUsedParameters)
    allParameters = Vector{Float64}(undef, numScale + numOffset + numVariance + numOptParameters)

    modelParameters = ModelParameters(scaleIndices, offsetIndices, varianceIndices, parameterIndices, 
        scaleMap, offsetMap, varianceMap, numScale, numOffset, numVariance, numOptParameters, numUsedParameters, 
        scaleNames, offsetNames, varianceNames, optParameterNames, 
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
    numTimeStepsForCondObs = Array{Int64}(undef, numConditions, numObservables)
    observedObservableForCond = Vector{Vector{Int64}}(undef, numConditions)

    observedAtForCond = Vector{Vector{Vector{Float64}}}(undef, numConditions)
    observedAtIndexForCondObs = Array{Vector{Int64}, 2}(undef, numConditions, numObservables)
    measurementForCondObs = Array{Vector{Float64}, 2}(undef, numConditions, numObservables)
    observablesTimeIndexIndicesForCond = Vector{Array{Vector{Int64}, 2}}(undef, numConditions)

    for (iCond, condId) in enumerate(experimentalConditions[!,1])
        inputParameterValuesForCond[iCond] = [experimentalConditions[iCond, inPar] for inPar in inputParameterSymbols]

        relevantMeasurementData = measurementData[measurementData[:,3] .== condId, :]
        observableIDs = relevantMeasurementData[:,1]
        timeStepsForCond[iCond] = sort(unique(relevantMeasurementData[:,5]))
        numTimeStepsForCond[iCond] = length(timeStepsForCond[iCond])
        isObserved = [observableNames[i] in unique(observableIDs) for i=1:numObservables]
        observedObservableForCond[iCond] = collect(1:numObservables)[isObserved]

        observedAtForCond[iCond] = Vector{Vector{Float64}}(undef, numObservables)
        for (i, obsId) = enumerate(observableNames)
            observedAtForCond[iCond][i] = relevantMeasurementData[observableIDs .== obsId, 5]
            numTimeStepsForCondObs[iCond, i] = length(observedAtForCond[iCond][i])
        end
        for i = 1:numObservables
            observedAtIndexForCondObs[iCond, i] = indexin(observedAtForCond[iCond][i], timeStepsForCond[iCond])
        end
        for (i, obsId) = enumerate(observableNames)
            measurementForCondObs[iCond, i] = relevantMeasurementData[observableIDs .== obsId, 4]
        end
        observablesTimeIndexIndicesForCond[iCond] = Array{Vector{Int64}, 2}(undef, numObservables, numTimeStepsForCond[iCond])
        for i = 1:numObservables
            for j = 1:numTimeStepsForCond[iCond]
                observablesTimeIndexIndicesForCond[iCond][i,j] = findall(x->x==j, observedAtIndexForCondObs[iCond, i])
            end
        end
    end

    experimentalData = ExperimentalData(measurementForCondObs, timeStepsForCond, observedAtIndexForCondObs, numConditions, 
            numObservables, numTimeStepsForCondObs, inputParameterValuesForCond, observedObservableForCond, observablesTimeIndexIndicesForCond)

    return experimentalData
end

function ModelOutput(usedType, experimentalData, modelParameters)
    numConditions = experimentalData.numConditions
    numObservables = experimentalData.numObservables
    measurementForCondObs = experimentalData.measurementForCondObs

    sols = Vector{ODESolution}(undef, numConditions)

    h_barForCondObs = convert(Array{Vector{usedType}, 2}, copy(measurementForCondObs))
    #Array{Vector{usedType}, 2}(undef, numConditions, numObservables)
    h_hatForCondObs = convert(Array{Vector{usedType}, 2}, copy(measurementForCondObs))
    costForCondObs = zeros(usedType, numConditions, numObservables)

    scaleGrad = Vector{Float64}(undef, length(modelParameters.scaleVector))
    offsetGrad = Vector{Float64}(undef, length(modelParameters.offsetVector))
    varianceGrad = Vector{Float64}(undef, length(modelParameters.varianceVector))
    dynamicParameterGrad = Vector{Float64}(undef, length(modelParameters.dynamicParametersVector))
    gradMatrix = Array{Float64, 2}(undef, numConditions, length(modelParameters.allParameters))

    modelOutput = ModelOutput(sols, h_barForCondObs, h_hatForCondObs, costForCondObs, scaleGrad, offsetGrad, varianceGrad, dynamicParameterGrad, gradMatrix)

    return modelOutput
end








function GradCalc_adjoint(filesAndPaths, timeEnd, experimentalConditions, measurementData, observables, parameterBounds)

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

    usedType = Float64
    modelOutput = ModelOutput(usedType, experimentalData, modelParameters)

    allStartParameters = getSamples(numAllStartParameters, lowerBounds, upperBounds, 1)

    # Initialize functions

    solveODESystem = (iCond) -> solveODESystem_proto(prob, modelParameters, modelData, modelOutput, iCond)

    g_unscaledObservables = (type, u, dynPar, i, iCond) -> g_unscaledObservables_proto(type, u, dynPar, i, iCond, modelData, experimentalData)

    g_scaledObservationFunctions = (type, h_bar, scaleVector, offsetVector, iCond) -> g_scaledObservationFunctions_proto(type, h_bar, scaleVector, offsetVector, iCond, modelParameters, experimentalData, modelData)

    g_cost = (type, h_hat, varianceVector, i, iCond) -> g_cost_proto(type, h_hat, varianceVector, i, iCond, modelParameters, experimentalData)

    g = (u, dynPar, scale, offset, variance, i, iCond; type = get_type([u, dynPar, scale, offset, variance])) -> g_proto(u, dynPar, scale, offset, variance, i, iCond, g_unscaledObservables, g_scaledObservationFunctions, g_cost; type = get_type([u, dynPar, scale, offset, variance]))

    G = (iCond) -> G_proto(solveODESystem, g, iCond, modelParameters, experimentalData, modelOutput)

    updateAllParameterVectors = () -> updateAllParameterVectors_proto(modelParameters, modelData)

    allConditionsCost = (p...) -> allConditionsCost_proto(parameterSpace, modelParameters, experimentalData, modelData, updateAllParameterVectors, G, p...)

    calcCostGrad = (dg!, iCond) -> calcCostGrad_proto(g, dg!, iCond, modelParameters, modelData, experimentalData, modelOutput)

    allConditionsCostGrad = (grad, p...) -> allConditionsCostGrad_proto(parameterSpace, modelParameters, modelData, modelOutput, experimentalData, 
    updateAllParameterVectors, calcCostGrad, g, grad, p...)




    #=
    calcUnscaledObservable = (iCond) -> calcUnscaledObservable_proto(prob, modelParameters, modelData, experimentalData, modelOutput, iCond)

    calcScaledObservable = (iCond) -> calcScaledObservable_proto(modelParameters, modelData, modelOutput, experimentalData, iCond)

    calcCost = (iCond) -> calcCost_proto(modelParameters, modelOutput, experimentalData, iCond)

    updateAllParameterVectors = () -> updateAllParameterVectors_proto(modelParameters, modelData)

    allConditionsCost = (p_tuple...) -> allConditionsCost_proto(parameterSpace, modelParameters, experimentalData, modelData, 
            updateAllParameterVectors, calcUnscaledObservable, calcScaledObservable, calcCost, p_tuple...)

    calcCostGrad = (dg!, iCond) -> calcCostGrad_proto(prob, modelParameters, modelData, experimentalData, modelOutput, dg!, iCond)

    allConditionsCostGrad = (grad, p_tuple...) -> allConditionsCostGrad_proto(parameterSpace, modelParameters, modelData, modelOutput, experimentalData, 
            updateAllParameterVectors, calcCostGrad, grad, p_tuple...)
    =#
    
    # test to add dg_0 dp
    
    #=
    cost = allConditionsCost(allStartParameters...)
    grad = Vector{Float64}(undef, numAllStartParameters)
    println("Cost: ", cost)
    allConditionsCostGrad(grad, allStartParameters...)
    println("Grad: ", grad)

    
    corrGrad = [-0.04021212913073861, -0.00027585198339496564, -4.884570903749361e-8, -5.124695936522999e-8, -3.840135151862373e-8, 3.9594942638311417e-16, 0.0044802381015997235, 2.4894798981731915e-14, -7.343705803174255e-14, -5.554622248177194e-13, -1.6578523867101942, -2.747962778043845e-5, -1.4342632477552844e-5, -3.425147799191447e-5, -3.6364458350250004e-8, -7.331175003006753e-7, 1.872764989572625e-5, -0.003735194572929876, 1.8422889138514828e-8, 2.38219097786897e-9, 4.5577015477591606e-8, 5.35382101515278e-9, 1.7569386641870226e-8, -0.00029369708070514736, -7.898328093305923e-5, 3.7313108443756425e-8, -9.932159057931086e-7, -1.3464983281842777e-5, 1.0143017765358342e-8, 1.421759138063622e-8, -0.0005868755764216341, 1.1771009720309698e-8, 1.7235046710579977e-9, 5.378407710609003e-10, -0.019108478958611327, -1.6227506754055986e-6, -719206.0479240953, 3.0454199264363384e-7, -5.434781215721821e-5, -3.486837933691301e-5, -482082.85690103454, -6727.395038761593, -130327.96359959114, 0.0001434012002170814, 0.00038384533078328527, 0.0009265076038541756, -30.29506511331806, -0.8938988334802855, -8.868310144605134, -7.1269707164188505, 0.0021917135347774216, -50186.661002802524, 0.4914585375205777, 0.07717135714862139, 0.5008262014771981, 0.10966560393014929, 0.11517909505138836, -2780.678500229983, -5534.423397431363, 0.5898581386710189, -92.57912784974596, -1668.1287728345274, 0.4308668927955358, 0.9323999833337377, -17086.186771010096, 0.6129562666031785, 0.40259091542548187, 0.06429112399970456, -4.1941685198729937e6, -373.51471896791566, -3.737850807820179, 41.91636647030552, -40196.9044119608, -11567.420819833222, -413693.857680959, 0.3208063311751015, -3017.6523743971875, -67.69635577781581, -294.83538492659324, 0.040585687154953164, 0.002604367033410381, 1.0311047418818973, -3.432764705340477e9, -4.238573899343332e8, -1.6086552461628643, -3331.2217561259977, -264.22933639162477, -1.4528681729763915e-18, -0.00814814952400725, -0.00015083065633797558, 4.2322187236844267e-7, 1.0763019688583152e-8, -4.905059065169967e-13, 114.99516255285361, -154.64010349043667, -131.32089631937237, -1.4228164300644703e-9, 5.669699407447015e-23, 3.9645971956524596e-22, -0.0646578263221752, -26.312475106815263, -2.983187872663642e-10, 5.844154681256194e-6, 2521.3503638135794, -2.8485791313453995e7, 1.0806041334887841e-12, -0.29440884155185726, -344425.2215469204, -184.2669024812995, 1.9915529518868966e-16, -5.470461495897123, 0.00010987442270134895, 0.00011467384126148284, -91014.49851909376, -0.010773468856076739]
    scaleIndices = modelParameters.scaleIndices
    offsetIndices = modelParameters.offsetIndices
    varianceIndices = modelParameters.varianceIndices
    parameterIndices = modelParameters.parameterIndices
    scaleGrad = grad[scaleIndices]
    scaleDiff = sum((scaleGrad .- corrGrad[scaleIndices]).^2)
    offsetGrad = grad[offsetIndices]
    offsetDiff = sum((offsetGrad .- corrGrad[offsetIndices]).^2)
    varianceGrad = grad[varianceIndices]
    varianceDiff = sum((varianceGrad .- corrGrad[varianceIndices]).^2)
    dynParGrad = grad[parameterIndices]
    dynParDiff = sum((dynParGrad .- corrGrad[parameterIndices]).^2)

    println("")
    println(scaleGrad, "\n")
    println(offsetGrad, "\n")
    println(varianceGrad, "\n")
    println(dynParGrad, "\n")
    println(scaleDiff, "\n", offsetDiff, "\n", varianceDiff, "\n", dynParDiff)
    =#

    #= Result
    Cost: 6.8094589076663265e6
    [-0.040215717384824115, -0.00027588171720311684, -4.884511958365966e-8, -5.124697517966012e-8, -3.840084530817628e-8, 3.9661675943746355e-16, 0.0044802381015997235, 2.4903914963718664e-14, -7.45107036384727e-14, -5.556572952747958e-13, -1.6578523867101942, -2.747963410327273e-5, -1.434263485866402e-5, -3.425148723370016e-5, -3.6363519326192874e-8, -7.331193447228313e-7, 1.872764989602129e-5, -0.003735223988700809, 1.842294492312467e-8, 2.3822003134824323e-9, 4.5577827430778716e-8, 5.3539923944404814e-9, 1.756952125227467e-8, -0.00029369823297981324, -7.898479459335733e-5, 3.73133355032059e-8, -9.93215905790229e-7, -1.3464983281862482e-5, 1.0143017765353918e-8, 1.4217591380636184e-8, -0.0005868755764205766, 1.177100972024535e-8, 1.7235046710568887e-9, 5.378407710608245e-10, -0.019108478976694793, -1.6227506754554586e-6, -719206.0479209339, 3.0454199268305825e-7, -5.434781216147312e-5, -3.4868379339027596e-5, -482082.8569011356, -6727.395038762998, -130327.96359959248]

    [0.0001434012002170814, 0.00038384533078328527, 0.0009265076038541755, -30.295065113317943, -0.893898833480176, -8.868310144710787, -7.126970716154661, 0.0021917135347774216, -50186.660999771186, 0.4914585375205778, 0.07717135714862139, 0.5008262014730508, 0.10966560393014926, 0.11517909505138794, -2780.678500199318, -5534.423396553522, 0.5898581386711746, -92.57912784974597, -1668.128772834527, 0.4308668927955358, 0.9323999833337377, -17086.186771010096, 0.6129562666031787, 0.40259091542548187, 0.06429112399970456, -4.1941685198729914e6, -373.5147189679157, -3.737850807820179, 41.91636647030553, -40196.9044119608, -11567.420819833224]

    [-413704.29918497166, 0.24303118863208248, -3141.056380987115, -66.38609866622598, -294.60883060635564, 0.04058568715495316, 0.002604366085771056, 1.4454627162219853, -3.4327647053404775e9, -4.2386307979601926e8, -1.6718980070056293, -3331.2217561259977]

    [0.5367225585749565, -1.4528925017114376e-18, -0.00814881459999649, -0.00015086729501357212, 4.2463044987726773e-7, 1.0766809693449906e-8, 3437.801904236584, 114.99518176479549, -154.6471915342424, -131.61392912755306, -1.4237807543577703e-9, 1.1042880813768036e-21, 3.964666282756976e-22, -0.06466581443241999, -26.314731954975766, -2.983628569154093e-10, 5.8450196776385934e-6, 2272.6370057796366, -2.8487092480833277e7, -6.460859714240395, -0.33530605764728555, -345427.76523977745, 213.74897485666904, -1.1153470637574956e-5, -5.4712426677535415, 9.704611048345861e-6, 1.4630987436878326e-5, -91220.7966619139, -9.170933357136406e-15]

    2.288182355434283e-11
    9.960563542425413e-12
    2.2317877196999846e9
    1.4849589130970541e7
    =#

    #= Reference
    Cost: 6.80957880225908e6
    [-0.04021212913073861, -0.00027585198339496564, -4.884570903749361e-8, -5.124695936522999e-8, -3.840135151862373e-8, 3.9594942638311417e-16, 0.0044802381015997235, 2.4894798981731915e-14, -7.343705803174255e-14, -5.554622248177194e-13, -1.6578523867101942, -2.747962778043845e-5, -1.4342632477552844e-5, -3.425147799191447e-5, -3.6364458350250004e-8, -7.331175003006753e-7, 1.872764989572625e-5, -0.003735194572929876, 1.8422889138514828e-8, 2.38219097786897e-9, 4.5577015477591606e-8, 5.35382101515278e-9, 1.7569386641870226e-8, -0.00029369708070514736, -7.898328093305923e-5, 3.7313108443756425e-8, -9.932159057931086e-7, -1.3464983281842777e-5, 1.0143017765358342e-8, 1.421759138063622e-8, -0.0005868755764216341, 1.1771009720309698e-8, 1.7235046710579977e-9, 5.378407710609003e-10, -0.019108478958611327, -1.6227506754055986e-6, -719206.0479240953, 3.0454199264363384e-7, -5.434781215721821e-5, -3.486837933691301e-5, -482082.85690103454, -6727.395038761593, -130327.96359959114]

    [0.0001434012002170814, 0.00038384533078328527, 0.0009265076038541756, -30.29506511331806, -0.8938988334802855, -8.868310144605134, -7.1269707164188505, 0.0021917135347774216, -50186.661002802524, 0.4914585375205777, 0.07717135714862139, 0.5008262014771981, 0.10966560393014929, 0.11517909505138836, -2780.678500229983, -5534.423397431363, 0.5898581386710189, -92.57912784974596, -1668.1287728345274, 0.4308668927955358, 0.9323999833337377, -17086.186771010096, 0.6129562666031785, 0.40259091542548187, 0.06429112399970456, -4.1941685198729937e6, -373.51471896791566, -3.737850807820179, 41.91636647030552, -40196.9044119608, -11567.420819833222]

    [-413693.857680959, 0.3208063311751015, -3017.6523743971875, -67.69635577781581, -294.83538492659324, 0.040585687154953164, 0.002604367033410381, 1.0311047418818973, -3.432764705340477e9, -4.238573899343332e8, -1.6086552461628643, -3331.2217561259977]

    [-264.22933639162477, -1.4528681729763915e-18, -0.00814814952400725, -0.00015083065633797558, 4.2322187236844267e-7, 1.0763019688583152e-8, -4.905059065169967e-13, 114.99516255285361, -154.64010349043667, -131.32089631937237, -1.4228164300644703e-9, 5.669699407447015e-23, 3.9645971956524596e-22, -0.0646578263221752, -26.312475106815263, -2.983187872663642e-10, 5.844154681256194e-6, 2521.3503638135794, -2.8485791313453995e7, 1.0806041334887841e-12, -0.29440884155185726, -344425.2215469204, -184.2669024812995, 1.9915529518868966e-16, -5.470461495897123, 0.00010987442270134895, 0.00011467384126148284, -91014.49851909376, -0.010773468856076739]
    =#


    model = Model(NLopt.Optimizer)
    #JuMP.register(model::Model, s::Symbol, dimension::Integer, f::Function, ∇f::Function, ∇²f::Function)
    register(model, :f, numAllStartParameters, allConditionsCost, allConditionsCostGrad)
    set_optimizer_attribute(model, "algorithm", :LD_MMA)
    
    @variable(model, lowerBounds[i] <= p[i = 1:numAllStartParameters] <= upperBounds[i])
    for i in 1:numAllStartParameters
        set_start_value(p[i], allStartParameters[i])
    end
    @NLobjective(model, Min, f(p...))
    JuMP.optimize!(model)
    #@show termination_status(model)
    #@show primal_status(model)
    p_opt = [value(p[i]) for i=1:numAllStartParameters]
    doLogSearch = parameterSpace.doLogSearch
    view(p_opt, doLogSearch) .= exp10.(view(p_opt, doLogSearch))
    println("Done!")
    cost_opt = objective_value(model)
    println(cost_opt)

    scaleIndices = modelParameters.scaleIndices
    offsetIndices = modelParameters.offsetIndices
    varianceIndices = modelParameters.varianceIndices
    parameterIndices = modelParameters.parameterIndices
    optimizedScales = p_opt[scaleIndices]
    optimizedOffsets = p_opt[offsetIndices]
    optimizedVariances = p_opt[varianceIndices]
    optimizedParameters = p_opt[parameterIndices]
    println(optimizedParameters)
    
    #return p_opt, cost_opt
    
    
    nothing
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
    GradCalc_adjoint(filesAndPaths, timeEnd, experimentalConditions, measurementData, observables, parameterBounds)

end

main()