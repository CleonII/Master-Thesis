

function solveODESystem_forwAD_proto_model_Alkan_SciSignal2018(prob::ODEProblem, solver, dynParVector::Vector{T1}, u0Vector::Vector{T2}, modelData::ModelData, 
    modelOutput::ModelOutput, iCond::Int64)::Nothing where {T1 <: Union{ForwardDiff.Dual, Float64}, T2 <: Union{ForwardDiff.Dual, Float64}}

    initVariable = modelData.initVariableIndices
    parameterInU0Indices = modelData.parameterInU0Indices
    u0Vector[initVariable[1]] = -dynParVector[parameterInU0Indices[1]] * dynParVector[parameterInU0Indices[2]] * (dynParVector[parameterInU0Indices[3]] - 1)
    u0Vector[initVariable[2]] = dynParVector[parameterInU0Indices[1]] * dynParVector[parameterInU0Indices[3]]
    u0Vector[initVariable[3]] = dynParVector[parameterInU0Indices[1]] * (dynParVector[parameterInU0Indices[2]] - 1) * (dynParVector[parameterInU0Indices[3]] - 1)
    
    _prob = remake(prob, u0 = u0Vector, p = dynParVector)
    modelOutput.sols[iCond] = OrdinaryDiffEq.solve(_prob, solver, reltol=1e-9, abstol=1e-9)
end


function calcUnscaledObservable_forwAD_proto_model_Alkan_SciSignal2018(dynParVector::Vector{T1}, sols::Vector{T2}, h_barFCO::Array{Vector{T1}, 2}, modelData::ModelData, experimentalData::ExperimentalData, 
    iCond::Int64)::Nothing where {T1 <: Union{ForwardDiff.Dual, Float64}, T2 <: Union{ODESolution, OrdinaryDiffEq.ODECompositeSolution}}

    timeSteps = experimentalData.timeStepsForCond[iCond]
    sol = sols[iCond](timeSteps)

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



function forwardAutomaticDifferentiation_model_Alkan_SciSignal2018(modelFunction::Function, iStartPar::Int64, optAlg::Symbol, solver,
    timeEnd::AbstractFloat, experimentalConditions::DataFrame, measurementData::DataFrame, observables::DataFrame, parameterBounds::DataFrame)

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

    parameterInObservableNames = String[]

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
    solveODESystem_float = (iCond) -> solveODESystem_forwAD_proto_model_Alkan_SciSignal2018(prob, solver, modelParameters.dynamicParametersVector, modelParameters.u0Vector, modelData, modelOutput_float, iCond)
    calcUnscaledObservable_float = (iCond) -> calcUnscaledObservable_forwAD_proto_model_Alkan_SciSignal2018(modelParameters.dynamicParametersVector, modelOutput_float.sols, modelOutput_float.h_barForCondObs, modelData, experimentalData, iCond)
    calcScaledObservable_float = (iCond) -> calcScaledObservable_forwAD_proto(modelParameters.scaleVector, modelParameters.offsetVector, modelOutput_float.h_barForCondObs, modelOutput_float.h_hatForCondObs, modelParameters, modelData, experimentalData, iCond)
    calcCost_float = (iCond) -> calcCost_forwAD_proto(modelParameters.varianceVector, modelOutput_float.h_hatForCondObs, modelOutput_float.costForCondObs, modelParameters, experimentalData, iCond)
    allConditionsCost_float = () -> allConditionsCost_forwAD_float_proto(modelParameters, experimentalData, modelData, 
            solveODESystem_float, calcUnscaledObservable_float, calcScaledObservable_float, calcCost_float)
    f = (p_tuple...) -> f_cost_forwAD_proto(updateAllParameterVectors, allConditionsCost_float, parameterSpace, modelParameters, p_tuple...)

    # For dynamic parameter grad
    solveODESystem_dual = (dualU0Vector, iCond) -> solveODESystem_forwAD_proto_model_Alkan_SciSignal2018(prob, solver, dualModelParameters.dualDynamicParametersVector, dualU0Vector, modelData, modelOutput_dual, iCond)
    calcUnscaledObservable_dual = (iCond) -> calcUnscaledObservable_forwAD_proto_model_Alkan_SciSignal2018(dualModelParameters.dualDynamicParametersVector, modelOutput_dual.sols, modelOutput_dual.h_barForCondObs, modelData, experimentalData, iCond)
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