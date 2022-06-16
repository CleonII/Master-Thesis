

function solveODESystem_forwAD_proto_model_Boehm_JProteomeRes2014(prob::ODEProblem, solver, dynParVector::Vector{T1}, u0Vector::Vector{T2}, modelData::ModelData, 
    modelOutput::ModelOutput, iCond::Int64)::Nothing where {T1 <: Union{ForwardDiff.Dual, Float64}, T2 <: Union{ForwardDiff.Dual, Float64}}

    initVariable = modelData.initVariableIndices
    parameterInU0Indices = modelData.parameterInU0Indices
    u0Vector[initVariable[1]] = 207.6 * dynParVector[parameterInU0Indices[1]]
    u0Vector[initVariable[2]] = 207.6 - 207.6 * dynParVector[parameterInU0Indices[1]]

    _prob = remake(prob, u0 = u0Vector, p = dynParVector)
    modelOutput.sols[iCond] = OrdinaryDiffEq.solve(_prob, solver, reltol=1e-9, abstol=1e-9)

    nothing
end


function calcUnscaledObservable_forwAD_proto_model_Boehm_JProteomeRes2014(dynParVector::Vector{T1}, sols::Vector{T2}, h_barFCO::Array{Vector{T1}, 2}, modelData::ModelData, experimentalData::ExperimentalData, 
    iCond::Int64)::Nothing where {T1 <: Union{ForwardDiff.Dual, Float64}, T2 <: Union{ODESolution, OrdinaryDiffEq.ODECompositeSolution}}

    timeSteps = experimentalData.timeStepsForCond[iCond]
    sol = sols[iCond](timeSteps)

    observedAIFCO = experimentalData.observedAtIndexForCondObs
    optParIOI = modelData.parameterInObservableIndices
    
    termAdd::Float64 = 0.0
    for odeSol in sols[iCond]
        if odeSol.retcode != :Success 
            termAdd = Inf
        else
            termAdd = 0.0
        end
    end

    h_barFCO[iCond, 1] = ( (100 * sol[Symbol("pApB(t)")] + 200 * sol[Symbol("pApA(t)")] * dynParVector[optParIOI[1]] ) ./ 
            ( sol[Symbol("pApB(t)")] + sol[Symbol("STAT5A(t)")] * dynParVector[optParIOI[1]] + 2 * sol[Symbol("pApA(t)")] * dynParVector[optParIOI[1]]) )[observedAIFCO[iCond, 1]] .+ termAdd
    h_barFCO[iCond, 2] = ( -1 * ( 100 * sol[Symbol("pApB(t)")] - 200 * sol[Symbol("pBpB(t)")] * (dynParVector[optParIOI[1]] - 1) ) ./ 
            ( sol[Symbol("STAT5B(t)")] * (dynParVector[optParIOI[1]] - 1) - sol[Symbol("pApB(t)")] + 2 * sol[Symbol("pBpB(t)")] * (dynParVector[optParIOI[1]] - 1) ) )[observedAIFCO[iCond, 2]] .+ termAdd
    h_barFCO[iCond, 3] = ( ( 100 * sol[Symbol("pApB(t)")] + 100 * sol[Symbol("STAT5A(t)")] * dynParVector[optParIOI[1]] + 200 * sol[Symbol("pApA(t)")] * dynParVector[optParIOI[1]] ) ./ 
            (2 * sol[Symbol("pApB(t)")] + sol[Symbol("STAT5A(t)")] * dynParVector[optParIOI[1]] + 2 * sol[Symbol("pApA(t)")] * dynParVector[optParIOI[1]] - 
            sol[Symbol("STAT5B(t)")] * (dynParVector[optParIOI[1]] - 1) - 2 * sol[Symbol("pBpB(t)")] * (dynParVector[optParIOI[1]] - 1)) )[observedAIFCO[iCond, 3]] .+ termAdd

    nothing
end

function calcScaledObservable_forwAD_proto_model_Boehm_JProteomeRes2014(scaleVector::Vector{T1}, offsetVector::Vector{T1}, h_barFCO::Array{Vector{T1}, 2}, h_hatFCO::Array{Vector{T1}, 2}, modelParameters::ModelParameters, 
        modelData::ModelData, experimentalData::ExperimentalData, iCond::Int64)::Nothing where T1 <: Union{ForwardDiff.Dual, Float64}

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

function calcCost_forwAD_proto_model_Boehm_JProteomeRes2014(varianceVector::Vector{T1}, h_hatFCO::Array{Vector{T1}, 2}, costFCO::Array{T1, 2}, modelParameters::ModelParameters,
        experimentalData::ExperimentalData, iCond::Int64)::T1 where T1 <: Union{ForwardDiff.Dual, Float64}

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
        costFCO[iCond, iObs] = log(2*pi*variance) * numDataFCO[iCond, iObs] + 
            (dot(measurementFCO[iCond, iObs], measurementFCO[iCond, iObs]) - 
            2*dot(measurementFCO[iCond, iObs], h_hatFCO[iCond, iObs]) + 
            dot(h_hatFCO[iCond, iObs], h_hatFCO[iCond, iObs])) / variance
    end

    return sum(costFCO[iCond, observedOFC[iCond]])
end

function allConditionsCost_forwAD_float_proto_model_Boehm_JProteomeRes2014(modelParameters::ModelParameters, experimentalData::ExperimentalData, modelData::ModelData, 
        solveODESystem::Function, calcUnscaledObservable::Function, calcScaledObservable::Function, calcCost::Function)::Float64

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

function allConditionsCost_forwAD_dual_proto_model_Boehm_JProteomeRes2014(modelParameters::ModelParameters, dualModelParameters::DualModelParameters, experimentalData::ExperimentalData, 
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
        solveODESystem(dualU0Vector, iCond)
        calcUnscaledObservable(iCond)
        calcScaledObservable(iCond)
        cost += calcCost(iCond)
    end

    return cost
end


function f_cost_forwAD_proto_model_Boehm_JProteomeRes2014(updateAllParameterVectors::Function, allConditionsCost_float::Function, 
        parameterSpace::ParameterSpace, modelParameters::ModelParameters, p...)::Float64
    
    doLogSearch = parameterSpace.doLogSearch
    allPar = modelParameters.allParameters
    allPar .= p
    view(allPar, doLogSearch) .= exp10.(view(allPar, doLogSearch))

    updateAllParameterVectors()

    cost = allConditionsCost_float()

    return cost
end

function f_grad_forwAD_proto_model_Boehm_JProteomeRes2014(grad, result::DiffResults.MutableDiffResult, allConditionsCost_dual::Function, cfg::GradientConfig, parameterSpace::ParameterSpace, 
        modelParameters::ModelParameters, modelOutput::ModelOutput, p...)::Nothing
    
    doLogSearch = parameterSpace.doLogSearch
    allPar = modelParameters.allParameters
    allPar .= p
    view(allPar, doLogSearch) .= exp10.(view(allPar, doLogSearch))

    ForwardDiff.gradient!(result, allConditionsCost_dual, allPar, cfg)

    allParametersGrad = modelOutput.allParametersGrad
    allParametersGrad[:] = DiffResults.gradient(result)
    view(allParametersGrad, doLogSearch) .*= view(allPar, doLogSearch) * log(10)

    grad[:] = allParametersGrad

    nothing
end

function f_hessian_forwAD_proto_model_Boehm_JProteomeRes2014(hess, result::DiffResults.MutableDiffResult, allConditionsCost_dual::Function, cfg::ForwardDiff.HessianConfig, parameterSpace::ParameterSpace, 
        modelParameters::ModelParameters, modelOutput::ModelOutput, p...)::Nothing

    doLogSearch = parameterSpace.doLogSearch
    allPar = modelParameters.allParameters
    allPar .= p
    view(allPar, doLogSearch) .= exp10.(view(allPar, doLogSearch))

    ForwardDiff.hessian!(result, allConditionsCost_dual, allPar)

    grad = modelOutput.allParametersGrad
    grad .= DiffResults.gradient(result) 
    hess .= DiffResults.hessian(result)

    for i in eachindex(allPar)
        if i in doLogSearch
            for j in eachindex(allPar)
                if i == j
                    hess[i, j] = log(10)^2 * allPar[i] * ( grad[i] + allPar[i] * hess[i, j] )
                elseif j in doLogSearch
                    hess[i, j] = log(10)^2 * allPar[i] * allPar[j] * hess[i, j]
                else
                    hess[i, j] = log(10) * allPar[i] * hess[i, j]
                end
            end
        else
            for j in doLogSearch
                hess[i, j] = log(10) * allPar[j] * hess[i, j]
            end
        end
    end

    nothing
end


function forwardAutomaticDifferentiation_hessian_proto_model_Boehm_JProteomeRes2014(modelFunction::Function, solver,
        experimentalConditions::DataFrame, measurementData::DataFrame, observables::DataFrame, parameterBounds::DataFrame)

    sys, initialSpeciesValues, trueParameterValues = modelFunction()
    new_sys = ode_order_lowering(sys)
    u0 = initialSpeciesValues
    pars = trueParameterValues 
    tspan = (0.0, 100.0)
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
        modelOutput_float = createModelOutput(OrdinaryDiffEq.ODECompositeSolution, Float64, experimentalData, modelParameters)
        modelOutput_dual = createModelOutput(OrdinaryDiffEq.ODECompositeSolution, ForwardDiff.Dual, experimentalData, modelParameters)
    else
        modelOutput_float = createModelOutput(ODESolution, Float64, experimentalData, modelParameters)
        modelOutput_dual = createModelOutput(ODESolution, ForwardDiff.Dual, experimentalData, modelParameters)
    end

    # Getting start values
    allStartParameters = zeros(Float64, numAllStartParameters)

    # Initialize functions
    # General
    updateAllParameterVectors = () -> updateAllParameterVectors_proto(modelParameters, modelData)
    updateAllDualParameterVectors = (p) -> updateAllDualParameterVectors_proto(modelParameters, dualModelParameters, modelData, p)

    # For cost calc
    solveODESystem_float = (iCond) -> solveODESystem_forwAD_proto_model_Boehm_JProteomeRes2014(prob, solver, modelParameters.dynamicParametersVector, modelParameters.u0Vector, modelData, modelOutput_float, iCond)
    calcUnscaledObservable_float = (iCond) -> calcUnscaledObservable_forwAD_proto_model_Boehm_JProteomeRes2014(modelParameters.dynamicParametersVector, modelOutput_float.sols, modelOutput_float.h_barForCondObs, modelData, experimentalData, iCond)
    calcScaledObservable_float = (iCond) -> calcScaledObservable_forwAD_proto_model_Boehm_JProteomeRes2014(modelParameters.scaleVector, modelParameters.offsetVector, modelOutput_float.h_barForCondObs, modelOutput_float.h_hatForCondObs, modelParameters, modelData, experimentalData, iCond)
    calcCost_float = (iCond) -> calcCost_forwAD_proto_model_Boehm_JProteomeRes2014(modelParameters.varianceVector, modelOutput_float.h_hatForCondObs, modelOutput_float.costForCondObs, modelParameters, experimentalData, iCond)
    allConditionsCost_float = () -> allConditionsCost_forwAD_float_proto_model_Boehm_JProteomeRes2014(modelParameters, experimentalData, modelData, 
            solveODESystem_float, calcUnscaledObservable_float, calcScaledObservable_float, calcCost_float)
    f = (p_tuple...) -> f_cost_forwAD_proto_model_Boehm_JProteomeRes2014(updateAllParameterVectors, allConditionsCost_float, parameterSpace, modelParameters, p_tuple...)

    # For gradient
    solveODESystem_dual = (dualU0Vector, iCond) -> solveODESystem_forwAD_proto_model_Boehm_JProteomeRes2014(prob, solver, dualModelParameters.dualDynamicParametersVector, dualU0Vector, modelData, modelOutput_dual, iCond)
    calcUnscaledObservable_dual = (iCond) -> calcUnscaledObservable_forwAD_proto_model_Boehm_JProteomeRes2014(dualModelParameters.dualDynamicParametersVector, modelOutput_dual.sols, modelOutput_dual.h_barForCondObs, modelData, experimentalData, iCond)
    calcScaledObservable_dual = (iCond) -> calcScaledObservable_forwAD_proto_model_Boehm_JProteomeRes2014(dualModelParameters.dualScaleVector, dualModelParameters.dualOffsetVector, modelOutput_dual.h_barForCondObs, modelOutput_dual.h_hatForCondObs, modelParameters, modelData, experimentalData, iCond)
    calcCost_dual = (iCond) -> calcCost_forwAD_proto_model_Boehm_JProteomeRes2014(dualModelParameters.dualVarianceVector, modelOutput_dual.h_hatForCondObs, modelOutput_dual.costForCondObs, modelParameters, experimentalData, iCond)
    allConditionsCost_dual = (p) -> allConditionsCost_forwAD_dual_proto_model_Boehm_JProteomeRes2014(modelParameters, dualModelParameters, experimentalData, modelData,
            updateAllDualParameterVectors, solveODESystem_dual, calcUnscaledObservable_dual, calcScaledObservable_dual, calcCost_dual, p)
    cfg_grad = GradientConfig(allConditionsCost_dual, allStartParameters, Chunk{numAllStartParameters}())
    result_grad = DiffResults.GradientResult(allStartParameters::Vector{Float64})
    f_grad = (grad, p_tuple...) -> f_grad_forwAD_proto_model_Boehm_JProteomeRes2014(grad, result_grad, allConditionsCost_dual, cfg_grad, parameterSpace, modelParameters, modelOutput_dual, p_tuple...)

    result_hess = DiffResults.HessianResult(allStartParameters::Vector{Float64})
    cfg_hess = ForwardDiff.HessianConfig(allConditionsCost_dual, result_hess, allStartParameters, Chunk{numAllStartParameters}())
    
    f_hess = (hess, p...) -> f_hessian_forwAD_proto_model_Boehm_JProteomeRes2014(hess, result_hess, allConditionsCost_dual, cfg_hess, parameterSpace, 
            modelParameters, modelOutput_dual, p...)

    return f, f_grad, f_hess, lowerBounds, upperBounds
end


# cost = 4.84858095419511e6
# grad = [-1.1126285781895082e7, -37423.589107106935, -221.234057192666, -8.005741701395443e-10, -4.557051303245403e-11, -4.669236946001506e-12, 8.366385806615513e-11, -2.0928875208997618e-15, -3.4398312033081646e-13, -11.055986833432552, -25.917687434590214]