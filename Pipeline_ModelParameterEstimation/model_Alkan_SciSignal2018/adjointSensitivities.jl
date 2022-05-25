

function solveODESystem_AdjSens_proto_model_Alkan_SciSignal2018(prob::ODEProblem, solver, dynParVector::Vector{T1}, u0Vector::Vector{T2}, modelData::ModelData, 
    modelOutput::ModelOutput, iCond::Int64)::Nothing where {T1 <: Union{ForwardDiff.Dual, Float64}, T2 <: Union{ForwardDiff.Dual, Float64}}

    initVariable = modelData.initVariableIndices
    parameterInU0Indices = modelData.parameterInU0Indices
    u0Vector[initVariable[1]] = -dynParVector[parameterInU0Indices[1]] * dynParVector[parameterInU0Indices[2]] * (dynParVector[parameterInU0Indices[3]] - 1)
    u0Vector[initVariable[2]] = dynParVector[parameterInU0Indices[1]] * dynParVector[parameterInU0Indices[3]]
    u0Vector[initVariable[3]] = dynParVector[parameterInU0Indices[1]] * (dynParVector[parameterInU0Indices[2]] - 1) * (dynParVector[parameterInU0Indices[3]] - 1)

    _prob = remake(prob, u0 = u0Vector, p = dynParVector)
    modelOutput.sols[iCond] = OrdinaryDiffEq.solve(_prob, solver, reltol=1e-9, abstol=1e-9)
    
    nothing
end


function g_unscaledObservables_AdjSens_proto_model_Alkan_SciSignal2018(type::Type, u::Vector{T1}, dynPar::Vector{T2}, i::Int64, iCond::Int64, modelData::ModelData, 
    experimentalData::ExperimentalData)::Vector{Vector{type}} where {T1 <: Union{ForwardDiff.Dual, Float64}, 
    T2 <: Union{ForwardDiff.Dual, Float64}}

    observableVI = modelData.observableVariableIndices
    observablesTIIFC = experimentalData.observablesTimeIndexIndicesForCond[iCond]
    
    h_bar = Vector{Vector{type}}(undef, experimentalData.numObservables)

    h_bar[1] = [ u[observableVI[1]] + u[observableVI[2]] ][ones(Int64, length(observablesTIIFC[1, i]))]  
    h_bar[2] = [ u[observableVI[3]] + u[observableVI[4]] + u[observableVI[5]] + u[observableVI[6]] + u[observableVI[7]] + 
            u[observableVI[8]] + u[observableVI[9]] + u[observableVI[10]] + u[observableVI[11]] + u[observableVI[12]] + 
            u[observableVI[1]] + u[observableVI[2]] ][ones(Int64, length(observablesTIIFC[2, i]))]
    h_bar[3] = [ u[observableVI[13]] + u[observableVI[14]] ][ones(Int64, length(observablesTIIFC[3, i]))]                                 
    h_bar[4] = [ u[observableVI[15]] + u[observableVI[16]] ][ones(Int64, length(observablesTIIFC[4, i]))]                                                                  
    h_bar[5] = [ u[observableVI[17]] + u[observableVI[18]] ][ones(Int64, length(observablesTIIFC[5, i]))]                                    
    h_bar[6] = [ u[observableVI[19]] + u[observableVI[20]] ][ones(Int64, length(observablesTIIFC[6, i]))]                                    
    h_bar[7] = [ u[observableVI[21]] + u[observableVI[22]] ][ones(Int64, length(observablesTIIFC[7, i]))]                                    
    h_bar[8] = [ u[observableVI[23]] + u[observableVI[24]] ][ones(Int64, length(observablesTIIFC[8, i]))]                                     
    h_bar[9] = [ u[observableVI[25]] + u[observableVI[26]] ][ones(Int64, length(observablesTIIFC[9, i]))]                               
    h_bar[10] = [ u[observableVI[27]] + u[observableVI[28]] ][ones(Int64, length(observablesTIIFC[10, i]))]                             
    h_bar[11] = [ u[observableVI[25]] + u[observableVI[26]] ][ones(Int64, length(observablesTIIFC[11, i]))]             
    h_bar[12] = [ u[observableVI[29]] + u[observableVI[30]] ][ones(Int64, length(observablesTIIFC[12, i]))]                                             

    return h_bar
end



function adjointSensitivities_model_Alkan_SciSignal2018(modelFunction::Function, iStartPar::Int64, senseAlg, optAlg::Symbol, solver, solver_adjoint, 
    timeEnd::AbstractFloat, experimentalConditions::DataFrame, measurementData::DataFrame, 
    observables::DataFrame, parameterBounds::DataFrame)

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
        modelOutput_float = createModelOutput(OrdinaryDiffEq.ODECompositeSolution, Float64, experimentalData, modelParameters)
        modelOutput_dual = createModelOutput(OrdinaryDiffEq.ODECompositeSolution, ForwardDiff.Dual, experimentalData, modelParameters)
    else
        modelOutput_float = createModelOutput(ODESolution, Float64, experimentalData, modelParameters)
        modelOutput_dual = createModelOutput(ODESolution, ForwardDiff.Dual, experimentalData, modelParameters)
    end

    ### Getting start values

    allStartParameters = getSamples(numAllStartParameters, lowerBounds, upperBounds, iStartPar)

    ### Initialize functions

    solveODESystem_float = (iCond) -> solveODESystem_AdjSens_proto_model_Alkan_SciSignal2018(prob, solver, modelParameters.dynamicParametersVector, modelParameters.u0Vector, modelData, modelOutput_float, iCond)

    solveODESystem_dual = (dualPars, dualU0Vector, iCond) -> solveODESystem_AdjSens_proto_model_Alkan_SciSignal2018(prob, solver, dualPars, dualU0Vector, modelData, modelOutput_dual, iCond)

    g_unscaledObservables = (type, u, dynPar, i, iCond) -> g_unscaledObservables_AdjSens_proto_model_Alkan_SciSignal2018(type, u, dynPar, i, iCond, modelData, experimentalData)

    g_scaledObservationFunctions = (type, h_bar, scaleVector, offsetVector, iCond) -> g_scaledObservationFunctions_AdjSens_proto(type, h_bar, scaleVector, offsetVector, iCond, modelParameters, experimentalData, modelData)

    g_cost = (type, h_hat, varianceVector, i, iCond) -> g_cost_AdjSens_proto(type, h_hat, varianceVector, i, iCond, modelParameters, experimentalData)

    g = (u, dynPar, scale, offset, variance, i, iCond; type = get_type([u, dynPar, scale, offset, variance])) -> g_AdjSens_proto(u, dynPar, scale, offset, variance, i, iCond, g_unscaledObservables, g_scaledObservationFunctions, g_cost; type)

    G = (iCond) -> G_AdjSens_proto(solveODESystem_float, g, iCond, modelParameters, experimentalData, modelOutput_float)

    G_specifiedDynPar = (specifiedDynPar, iCond) -> G_specifiedDynPar_AdjSens_proto(specifiedDynPar, modelData.parameterInU0Indices, g, solveODESystem_dual, iCond, modelParameters, dualModelParameters, experimentalData, modelOutput_dual)

    updateAllParameterVectors = () -> updateAllParameterVectors_proto(modelParameters, modelData)

    allConditionsCost = (p...) -> allConditionsCost_AdjSens_proto(parameterSpace, modelParameters, experimentalData, modelData, updateAllParameterVectors, G, p...)

    calcCostGrad = (dg!, iCond) -> calcCostGrad_AdjSens_proto(g, dg!, G_specifiedDynPar, iCond, modelParameters, modelData, experimentalData, modelOutput_float, senseAlg, solver_adjoint)

    allConditionsCostGrad = (grad, p...) -> allConditionsCostGrad_AdjSens_proto(parameterSpace, modelParameters, modelData, modelOutput_float, experimentalData, 
            updateAllParameterVectors, calcCostGrad, g, grad, p...)

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
    register(model, :f, numAllStartParameters, allConditionsCost, allConditionsCostGrad)
    
    @variable(model, lowerBounds[i] <= p[i = 1:numAllStartParameters] <= upperBounds[i])
    for i in 1:numAllStartParameters
        set_start_value(p[i], allStartParameters[i])
    end
    @NLobjective(model, Min, f(p...))
    
    return model, p, parameterSpace.doLogSearch
    
end
