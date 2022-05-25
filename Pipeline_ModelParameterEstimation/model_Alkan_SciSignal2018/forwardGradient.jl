

function solveODESystem_forwGrad_proto_model_Alkan_SciSignal2018(prob::ODEProblem, solver, dynParVector::Vector{T1}, u0Vector::Vector{T2}, modelData::ModelData, 
    modelOutput::ModelOutput, iCond::Int64)::Bool where {T1 <: ForwardDiff.Dual, T2 <: ForwardDiff.Dual}

    initVariable = modelData.initVariableIndices
    parameterInU0Indices = modelData.parameterInU0Indices
    u0Vector[initVariable[1]] = -dynParVector[parameterInU0Indices[1]] * dynParVector[parameterInU0Indices[2]] * (dynParVector[parameterInU0Indices[3]] - 1)
    u0Vector[initVariable[2]] = dynParVector[parameterInU0Indices[1]] * dynParVector[parameterInU0Indices[3]]
    u0Vector[initVariable[3]] = dynParVector[parameterInU0Indices[1]] * (dynParVector[parameterInU0Indices[2]] - 1) * (dynParVector[parameterInU0Indices[3]] - 1)

    _prob = remake(prob, u0 = u0Vector, p = dynParVector)

    try
        modelOutput.sols[iCond] = OrdinaryDiffEq.solve(_prob, solver, reltol=1e-9, abstol=1e-9)
        if minimum(modelOutput.sols[iCond][:,:]).value < 0.0 || modelOutput.sols[iCond].t[end] !== prob.tspan[2]
            return true
        else
            return false
        end
    catch err
        return true
    end

end

function calcUnscaledObservable_forwGrad_proto_model_Alkan_SciSignal2018(dynParVector::Vector{T1}, sols::Vector{T2}, h_barFCO::Array{Vector{T1}, 2}, modelData::ModelData, 
    experimentalData::ExperimentalData, iCond::Int64)::Nothing where {T1 <: ForwardDiff.Dual, T2 <: Union{ODESolution, OrdinaryDiffEq.ODECompositeSolution}}
    
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



function forwardGradient_model_Alkan_SciSignal2018(modelFunction::Function, iStartPar::Int64, n_it::Int64, b2::Float64, stepRange::Vector{<:AbstractFloat}, solver, 
    timeEnd::Float64, experimentalConditions::DataFrame, measurementData::DataFrame, observables::DataFrame, parameterBounds::DataFrame)::Tuple{Function, Adam, Vector{Int64}}

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

    solveODESystem_dual = (u0Vector, iCond) -> solveODESystem_forwGrad_proto_model_Alkan_SciSignal2018(prob, solver, dualModelParameters.dualDynamicParametersVector, u0Vector, 
            modelData, modelOutput_dual, iCond)

    calcUnscaledObservable_dual = (iCond) -> calcUnscaledObservable_forwGrad_proto_model_Alkan_SciSignal2018(dualModelParameters.dualDynamicParametersVector, modelOutput_dual.sols, 
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

    stepLengths = exp10.(range(stepRange[1], stepRange[2], n_it))
    adam_opt = createAdam(allStartParameters, grad_g, stepLengths; b2 = b2)
    step = () -> step_forwGrad_proto!(adam_opt, keepInBounds!)

    return step, adam_opt, parameterSpace.doLogSearch
end
