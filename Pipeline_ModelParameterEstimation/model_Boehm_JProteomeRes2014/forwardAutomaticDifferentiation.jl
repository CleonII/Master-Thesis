

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

    h_barFCO[iCond, 1] = ( (100 * sol[Symbol("pApB(t)")] + 200 * sol[Symbol("pApA(t)")] * dynParVector[optParIOI[1]] ) ./ 
            ( sol[Symbol("pApB(t)")] + sol[Symbol("STAT5A(t)")] * dynParVector[optParIOI[1]] + 2 * sol[Symbol("pApA(t)")] * dynParVector[optParIOI[1]]) )[observedAIFCO[iCond, 1]]
    h_barFCO[iCond, 2] = ( -1 * ( 100 * sol[Symbol("pApB(t)")] - 200 * sol[Symbol("pBpB(t)")] * (dynParVector[optParIOI[1]] - 1) ) ./ 
            ( sol[Symbol("STAT5B(t)")] * (dynParVector[optParIOI[1]] - 1) - sol[Symbol("pApB(t)")] + 2 * sol[Symbol("pBpB(t)")] * (dynParVector[optParIOI[1]] - 1) ) )[observedAIFCO[iCond, 2]]
    h_barFCO[iCond, 3] = ( ( 100 * sol[Symbol("pApB(t)")] + 100 * sol[Symbol("STAT5A(t)")] * dynParVector[optParIOI[1]] + 200 * sol[Symbol("pApA(t)")] * dynParVector[optParIOI[1]] ) ./ 
            (2 * sol[Symbol("pApB(t)")] + sol[Symbol("STAT5A(t)")] * dynParVector[optParIOI[1]] + 2 * sol[Symbol("pApA(t)")] * dynParVector[optParIOI[1]] - 
            sol[Symbol("STAT5B(t)")] * (dynParVector[optParIOI[1]] - 1) - 2 * sol[Symbol("pBpB(t)")] * (dynParVector[optParIOI[1]] - 1)) )[observedAIFCO[iCond, 3]]

    nothing
end



function forwardAutomaticDifferentiation_proto_model_Boehm_JProteomeRes2014(modelFunction::Function, iStartPar::Int64, optAlg::Symbol, solver,
        timeEnd::AbstractFloat, experimentalConditions::DataFrame, measurementData::DataFrame, observables::DataFrame, parameterBounds::DataFrame)

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
        modelOutput_float = createModelOutput(OrdinaryDiffEq.ODECompositeSolution, Float64, experimentalData, modelParameters)
        modelOutput_dual = createModelOutput(OrdinaryDiffEq.ODECompositeSolution, ForwardDiff.Dual, experimentalData, modelParameters)
    else
        modelOutput_float = createModelOutput(ODESolution, Float64, experimentalData, modelParameters)
        modelOutput_dual = createModelOutput(ODESolution, ForwardDiff.Dual, experimentalData, modelParameters)
    end

    # Getting start values

    allStartParameters = getSamples(numAllStartParameters, lowerBounds, upperBounds, iStartPar)

    ### Initialize functions
    # General
    updateAllParameterVectors = () -> updateAllParameterVectors_proto(modelParameters, modelData)
    updateAllDualParameterVectors = (p) -> updateAllDualParameterVectors_proto(modelParameters, dualModelParameters, modelData, p)
    updateNonDynamicDualParameterVectors = (p) -> updateNonDynamicDualParameterVectors_proto(modelParameters, dualModelParameters, p)

    # For cost calc
    solveODESystem_float = (iCond) -> solveODESystem_forwAD_proto_model_Boehm_JProteomeRes2014(prob, solver, modelParameters.dynamicParametersVector, modelParameters.u0Vector, modelData, modelOutput_float, iCond)
    calcUnscaledObservable_float = (iCond) -> calcUnscaledObservable_forwAD_proto_model_Boehm_JProteomeRes2014(modelParameters.dynamicParametersVector, modelOutput_float.sols, modelOutput_float.h_barForCondObs, modelData, experimentalData, iCond)
    calcScaledObservable_float = (iCond) -> calcScaledObservable_forwAD_proto(modelParameters.scaleVector, modelParameters.offsetVector, modelOutput_float.h_barForCondObs, modelOutput_float.h_hatForCondObs, modelParameters, modelData, experimentalData, iCond)
    calcCost_float = (iCond) -> calcCost_forwAD_proto(modelParameters.varianceVector, modelOutput_float.h_hatForCondObs, modelOutput_float.costForCondObs, modelParameters, experimentalData, iCond)
    allConditionsCost_float = () -> allConditionsCost_forwAD_float_proto(modelParameters, experimentalData, modelData, 
            solveODESystem_float, calcUnscaledObservable_float, calcScaledObservable_float, calcCost_float)
    f = (p_tuple...) -> f_cost_forwAD_proto(updateAllParameterVectors, allConditionsCost_float, parameterSpace, modelParameters, p_tuple...)

    # For dynamic parameter grad
    solveODESystem_dual = (dualU0Vector, iCond) -> solveODESystem_forwAD_proto_model_Boehm_JProteomeRes2014(prob, solver, dualModelParameters.dualDynamicParametersVector, dualU0Vector, modelData, modelOutput_dual, iCond)
    calcUnscaledObservable_dual = (iCond) -> calcUnscaledObservable_forwAD_proto_model_Boehm_JProteomeRes2014(dualModelParameters.dualDynamicParametersVector, modelOutput_dual.sols, modelOutput_dual.h_barForCondObs, modelData, experimentalData, iCond)
    calcScaledObservable_dual_float = (iCond) -> calcScaledObservable_forwAD_proto(modelParameters.scaleVector, modelParameters.offsetVector, modelOutput_dual.h_barForCondObs, modelOutput_dual.h_hatForCondObs, modelParameters, modelData, experimentalData, iCond)
    calcCost_dual_float = (iCond) -> calcCost_forwAD_proto(modelParameters.varianceVector, modelOutput_dual.h_hatForCondObs, modelOutput_dual.costForCondObs, modelParameters, experimentalData, iCond)
    allConditionsCost_dual = (p) -> allConditionsCost_forwAD_dual_proto(modelParameters, dualModelParameters, experimentalData, modelData,
            solveODESystem_dual, calcUnscaledObservable_dual, calcScaledObservable_dual_float, calcCost_dual_float, p)

    # for non-dynamic parameter grad
    calcScaledObservable_float_dual = (iCond) -> calcScaledObservable_forwAD_proto(dualModelParameters.dualScaleVector, dualModelParameters.dualOffsetVector, modelOutput_float.h_barForCondObs, modelOutput_dual.h_hatForCondObs, modelParameters, modelData, experimentalData, iCond)
    calcCost_dual = (iCond) -> calcCost_forwAD_proto(dualModelParameters.dualVarianceVector, modelOutput_dual.h_hatForCondObs, modelOutput_dual.costForCondObs, modelParameters, experimentalData, iCond)
    allConditionsCost_noODE_dual = (p) -> allConditionsCost_noODE_forwAD_dual_proto(experimentalData, updateNonDynamicDualParameterVectors, 
            calcScaledObservable_float_dual, calcCost_dual, p)
    

    dynPar = allStartParameters[modelParameters.parameterIndices]
    noDynPar = allStartParameters[vcat(modelParameters.scaleIndices, modelParameters.offsetIndices, modelParameters.varianceIndices)]
    result1 = DiffResults.GradientResult(noDynPar::Vector{Float64})
    cfg = GradientConfig(allConditionsCost_dual, dynPar, Chunk{length(dynPar)}())
    result2 = DiffResults.GradientResult(dynPar::Vector{Float64})
    
    f_grad = (grad, p_tuple...) -> f_grad_forwAD_proto(grad, [result1, result2], updateAllParameterVectors, allConditionsCost_noODE_dual, allConditionsCost_dual, cfg, parameterSpace, modelParameters, modelOutput_dual, p_tuple...)


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
