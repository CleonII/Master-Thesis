

function solveODESystem_forwAD_proto_model_Bachmann_MSB2011(prob::ODEProblem, solver, dynParVector::Vector{T1}, u0Vector::Vector{T2}, modelData::ModelData, 
        modelOutput::ModelOutput, iCond::Int64)::Nothing where {T1 <: Union{ForwardDiff.Dual, Float64}, T2 <: Union{ForwardDiff.Dual, Float64}}

    initVariable = modelData.initVariableIndices
    parameterInU0Indices = modelData.parameterInU0Indices
    u0Vector[initVariable[1]] = dynParVector[parameterInU0Indices[1]]
    u0Vector[initVariable[2]] = dynParVector[parameterInU0Indices[2]] * (dynParVector[parameterInU0Indices[3]] * dynParVector[parameterInU0Indices[4]] + 1)
    u0Vector[initVariable[3]] = dynParVector[parameterInU0Indices[5]]
    u0Vector[initVariable[4]] = dynParVector[parameterInU0Indices[6]]
    u0Vector[initVariable[5]] = dynParVector[parameterInU0Indices[7]] * dynParVector[parameterInU0Indices[8]] * dynParVector[parameterInU0Indices[9]]
    u0Vector[initVariable[6]] = dynParVector[parameterInU0Indices[10]] * dynParVector[parameterInU0Indices[11]] * dynParVector[parameterInU0Indices[12]]

    _prob = remake(prob, u0 = u0Vector, p = dynParVector)
    modelOutput.sols[iCond] = OrdinaryDiffEq.solve(_prob, solver, reltol=1e-9, abstol=1e-9)

    nothing
end


function calcUnscaledObservable_forwAD_proto_model_Bachmann_MSB2011(dynParVector::Vector{T1}, sols::Vector{T2}, h_barFCO::Array{Vector{T1}, 2}, modelData::ModelData, experimentalData::ExperimentalData, 
        iCond::Int64)::Nothing where {T1 <: Union{ForwardDiff.Dual, Float64}, T2 <: Union{ODESolution, OrdinaryDiffEq.ODECompositeSolution}}

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


function forwardAutomaticDifferentiation_model_Bachmann_MSB2011(modelFunction::Function, iStartPar::Int64, optAlg::Symbol, solver,
        timeEnd::AbstractFloat, experimentalConditions::DataFrame, measurementData::DataFrame, observables::DataFrame, parameterBounds::DataFrame)

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

    ### Initialize structs

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

    ### Getting start values

    allStartParameters = getSamples(numAllStartParameters, lowerBounds, upperBounds, iStartPar)

    ### Initialize functions
    # General
    updateAllParameterVectors = () -> updateAllParameterVectors_proto(modelParameters, modelData)
    updateAllDualParameterVectors = (p) -> updateAllDualParameterVectors_proto(modelParameters, dualModelParameters, modelData, p)
    updateNonDynamicDualParameterVectors = (p) -> updateNonDynamicDualParameterVectors_proto(modelParameters, dualModelParameters, p)

    # For cost calc
    solveODESystem_float = (iCond) -> solveODESystem_forwAD_proto_model_Bachmann_MSB2011(prob, solver, modelParameters.dynamicParametersVector, modelParameters.u0Vector, modelData, modelOutput_float, iCond)
    calcUnscaledObservable_float = (iCond) -> calcUnscaledObservable_forwAD_proto_model_Bachmann_MSB2011(modelParameters.dynamicParametersVector, modelOutput_float.sols, modelOutput_float.h_barForCondObs, modelData, experimentalData, iCond)
    calcScaledObservable_float = (iCond) -> calcScaledObservable_forwAD_proto(modelParameters.scaleVector, modelParameters.offsetVector, modelOutput_float.h_barForCondObs, modelOutput_float.h_hatForCondObs, modelParameters, modelData, experimentalData, iCond)
    calcCost_float = (iCond) -> calcCost_forwAD_proto(modelParameters.varianceVector, modelOutput_float.h_hatForCondObs, modelOutput_float.costForCondObs, modelParameters, experimentalData, iCond)
    allConditionsCost_float = () -> allConditionsCost_forwAD_float_proto(modelParameters, experimentalData, modelData, 
            solveODESystem_float, calcUnscaledObservable_float, calcScaledObservable_float, calcCost_float)
    f = (p_tuple...) -> f_cost_forwAD_proto(updateAllParameterVectors, allConditionsCost_float, parameterSpace, modelParameters, p_tuple...)

    # For dynamic parameter grad
    solveODESystem_dual = (dualU0Vector, iCond) -> solveODESystem_forwAD_proto_model_Bachmann_MSB2011(prob, solver, dualModelParameters.dualDynamicParametersVector, dualU0Vector, modelData, modelOutput_dual, iCond)
    calcUnscaledObservable_dual = (iCond) -> calcUnscaledObservable_forwAD_proto_model_Bachmann_MSB2011(dualModelParameters.dualDynamicParametersVector, modelOutput_dual.sols, modelOutput_dual.h_barForCondObs, modelData, experimentalData, iCond)
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