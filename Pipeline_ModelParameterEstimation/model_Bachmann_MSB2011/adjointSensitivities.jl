

function solveODESystem_AdjSens_proto_model_Bachmann_MSB2011(prob::ODEProblem, solver, dynParVector::Vector{T1}, U0Vector::Vector{T2}, modelData::ModelData, 
        modelOutput::ModelOutput, iCond::Int64)::Nothing where {T1 <: Union{ForwardDiff.Dual, Float64}, T2 <: Union{ForwardDiff.Dual, Float64}}

    initVariable = modelData.initVariableIndices
    parameterInU0Indices = modelData.parameterInU0Indices
    U0Vector[initVariable[1]] = dynParVector[parameterInU0Indices[1]]
    U0Vector[initVariable[2]] = dynParVector[parameterInU0Indices[2]] * (dynParVector[parameterInU0Indices[3]] * dynParVector[parameterInU0Indices[4]] + 1)
    U0Vector[initVariable[3]] = dynParVector[parameterInU0Indices[5]]
    U0Vector[initVariable[4]] = dynParVector[parameterInU0Indices[6]]
    U0Vector[initVariable[5]] = dynParVector[parameterInU0Indices[7]] * dynParVector[parameterInU0Indices[8]] * dynParVector[parameterInU0Indices[9]]
    U0Vector[initVariable[6]] = dynParVector[parameterInU0Indices[10]] * dynParVector[parameterInU0Indices[11]] * dynParVector[parameterInU0Indices[12]]

    _prob = remake(prob, u0 = U0Vector, p = dynParVector)
    modelOutput.sols[iCond] = OrdinaryDiffEq.solve(_prob, solver, reltol=1e-9, abstol=1e-9)
    
    nothing
end


function g_unscaledObservables_AdjSens_proto_model_Bachmann_MSB2011(type::Type, u::Vector{T1}, dynPar::Vector{T2}, i::Int64, iCond::Int64, modelData::ModelData, 
        experimentalData::ExperimentalData)::Vector{Vector{type}} where {T1 <: Union{ForwardDiff.Dual, Float64}, 
        T2 <: Union{ForwardDiff.Dual, Float64}}

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
    h_bar[8] = [ u[observableVI[3]] + u[observableVI[4]] ][ones(Int64, length(observablesTIIFC[8, i]))]                                                # log10
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



function adjointSensitivities_model_Bachmann_MSB2011(modelFunction::Function, iStartPar::Int64, senseAlg, optAlg::Symbol, solver, solver_adjoint, 
        timeEnd::AbstractFloat, experimentalConditions::DataFrame, measurementData::DataFrame, 
        observables::DataFrame, parameterBounds::DataFrame)

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

    solveODESystem_float = (iCond) -> solveODESystem_AdjSens_proto_model_Bachmann_MSB2011(prob, solver, modelParameters.dynamicParametersVector, modelParameters.u0Vector, modelData, modelOutput_float, iCond)

    solveODESystem_dual = (dualPars, dualU0Vector, iCond) -> solveODESystem_AdjSens_proto_model_Bachmann_MSB2011(prob, solver, dualPars, dualU0Vector, modelData, modelOutput_dual, iCond)

    g_unscaledObservables = (type, u, dynPar, i, iCond) -> g_unscaledObservables_AdjSens_proto_model_Bachmann_MSB2011(type, u, dynPar, i, iCond, modelData, experimentalData)

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
