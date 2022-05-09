

function solveODESystem_AdjSens_proto(prob, dynParVector, U0Vector, modelData, modelOutput, iCond)

    initVariable = modelData.initVariableIndices
    parameterInU0Indices = modelData.parameterInU0Indices
    U0Vector[initVariable[1]] = dynParVector[parameterInU0Indices[1]]
    U0Vector[initVariable[2]] = dynParVector[parameterInU0Indices[2]] * (dynParVector[parameterInU0Indices[3]] * dynParVector[parameterInU0Indices[4]] + 1)
    U0Vector[initVariable[3]] = dynParVector[parameterInU0Indices[5]]
    U0Vector[initVariable[4]] = dynParVector[parameterInU0Indices[6]]
    U0Vector[initVariable[5]] = dynParVector[parameterInU0Indices[7]] * dynParVector[parameterInU0Indices[8]] * dynParVector[parameterInU0Indices[9]]
    U0Vector[initVariable[6]] = dynParVector[parameterInU0Indices[10]] * dynParVector[parameterInU0Indices[11]] * dynParVector[parameterInU0Indices[12]]

    _prob = remake(prob, u0 = U0Vector, p = dynParVector)
    modelOutput.sols[iCond] = OrdinaryDiffEq.solve(_prob, Rodas5(), reltol=1e-9, abstol=1e-9)
    
    nothing
end


function g_unscaledObservables_AdjSens_proto(type, u, dynPar, i, iCond, modelData, experimentalData)

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

function g_scaledObservationFunctions_AdjSens_proto(type, h_bar, scaleVector, offsetVector, iCond, modelParameters, experimentalData, modelData)
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

function g_cost_AdjSens_proto(type, h_hat, varianceVector, i, iCond, modelParameters, experimentalData)
    varianceMap = modelParameters.varianceMap
    measurementFCO = experimentalData.measurementForCondObs
    observablesTIIFC = experimentalData.observablesTimeIndexIndicesForCond[iCond]
    observedOFC = experimentalData.observedObservableForCond

    costFCO = Vector{type}(undef, experimentalData.numObservables)
    for iObs in observedOFC[iCond]
        costFCO[iObs] = log(2*pi*varianceVector[varianceMap[iCond, iObs]]) * length(h_hat[iObs]) + 
                (dot(measurementFCO[iCond, iObs][observablesTIIFC[iObs, i]], measurementFCO[iCond, iObs][observablesTIIFC[iObs, i]]) - 
                2*dot(measurementFCO[iCond, iObs][observablesTIIFC[iObs, i]], h_hat[iObs]) + 
                dot(h_hat[iObs], h_hat[iObs])) / (varianceVector[varianceMap[iCond, iObs]])
    end

    return sum(costFCO[observedOFC[iCond]])
end



function g_AdjSens_proto(u, dynPar, scale, offset, variance, i, iCond, g_unscaledObservables, g_scaledObservationFunctions, g_cost; type = get_type([u, dynPar, scale, offset, variance]))

    h_bar = g_unscaledObservables(type, u, dynPar, i, iCond)

    h_hat = g_scaledObservationFunctions(type, h_bar, scale, offset, iCond)

    cost = g_cost(type, h_hat, variance, i, iCond)

    return cost
end

function G_AdjSens_proto(solveODESystem, g, iCond, modelParameters, experimentalData, modelOutput)
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

function G_specifiedDynPar_AdjSens_proto(specifiedDynPar, specifiedDynParIndices, g, solveODESystem, iCond, modelParameters, dualModelParameters, experimentalData, modelOutput)
    
    dynParVector = modelParameters.dynamicParametersVector
    dualDynParVector = dualModelParameters.dualDynamicParametersVector
    dualDynParVector .= convert.(eltype(specifiedDynPar), dynParVector)
    dualDynParVector[specifiedDynParIndices] .= specifiedDynPar

    u0Vector = modelParameters.u0Vector
    dualU0Vector = convert(Vector{eltype(specifiedDynPar)}, u0Vector)

    solveODESystem(dualDynParVector, dualU0Vector, iCond)

    scale = modelParameters.scaleVector
    offset = modelParameters.offsetVector
    variance = modelParameters.varianceVector

    costAt = (u, i) -> g(u, dualDynParVector, scale, offset, variance, i, iCond, type = eltype(specifiedDynPar))
    sol = modelOutput.sols[iCond]

    cost = 0.0
    ts = experimentalData.timeStepsForCond[iCond]
    for (i, t) in enumerate(ts)
        cost += costAt(sol(t), i)
    end

    return cost

end



function allConditionsCost_AdjSens_proto(parameterSpace, modelParameters, experimentalData, modelData, 
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

    return cost
end

function dg_AdjSens_proto!(out, u, p, t, i, gu) 

    guAt = (u) -> gu(u, i)
    out .= ForwardDiff.gradient(guAt, u)

    nothing
end

function calcCostGrad_AdjSens_proto(g, dg!, G_specifiedDynPar, iCond, modelParameters, modelData, experimentalData, modelOutput, senseAlg)

    timeSteps = experimentalData.timeStepsForCond[iCond]
    dynParGrad = modelOutput.dynParGrad

    sol = modelOutput.sols[iCond]

    ~, dynParGrad[:] = adjoint_sensitivities(sol, Rodas4P(), dg!, timeSteps, 
        sensealg = senseAlg, reltol = 1e-9, abstol = 1e-9)

    # when a parameter is included in the observation function the gradient is incorrect, has to correct with adding dgdp 

    # gradient for scales, offsets, variances and dgdp 
    dynPar = modelParameters.dynamicParametersVector
    ∂g∂p = modelOutput.∂g∂p
    ∂g∂p .= 0.0
    scale = modelParameters.scaleVector
    scaleGrad = modelOutput.scaleGrad
    scaleGrad .= 0.0
    offset = modelParameters.offsetVector
    offsetGrad = modelOutput.offsetGrad
    offsetGrad .= 0.0
    variance = modelParameters.varianceVector
    varianceGrad = modelOutput.varianceGrad
    varianceGrad .= 0.0

    gdynPar = (u, dynPar, i) -> g(u, dynPar, scale, offset, variance, i, iCond)
    gscale = (u, scale, i) -> g(u, dynPar, scale, offset, variance, i, iCond)
    goffset = (u, offset, i) -> g(u, dynPar, scale, offset, variance, i, iCond)
    gvariance = (u, variance, i) -> g(u, dynPar, scale, offset, variance, i, iCond)
    for (i, t) in enumerate(timeSteps)
        gdynParAt = (dynPar) -> gdynPar(sol(t), dynPar, i)
        ∂g∂p += ForwardDiff.gradient(gdynParAt, dynPar)

        gscaleAt = (scale) -> gscale(sol(t), scale, i)
        scaleGrad += ForwardDiff.gradient(gscaleAt, scale)

        goffsetAt = (offset) -> goffset(sol(t), offset, i)
        offsetGrad += ForwardDiff.gradient(goffsetAt, offset)

        gvarianceAt = (variance) -> gvariance(sol(t), variance, i)
        varianceGrad += ForwardDiff.gradient(gvarianceAt, variance)
    end

    # For parameters that influence the starting concentrations 

    specifiedDynParIndices = modelData.parameterInU0Indices
    specifiedDynPar = dynPar[specifiedDynParIndices]

    G_specDynPar = (specifiedDynPar) -> G_specifiedDynPar(specifiedDynPar, iCond)
    correctSpecDynParGrad = ForwardDiff.gradient(G_specDynPar, specifiedDynPar)

    dynParGrad *= -1 # adjoint sensitivities seems to return the negative gradient
    dynParGrad += ∂g∂p
    dynParGrad[specifiedDynParIndices] = correctSpecDynParGrad

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
    gradMatrix[iCond, parameterIndices] = view(dynParGrad, optParameterIndices)

    nothing
end


function allConditionsCostGrad_AdjSens_proto(parameterSpace, modelParameters, modelData, modelOutput, experimentalData, 
    updateAllParameterVectors, calcCostGrad, g, grad, p...)

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
        dg! = (out, u, p, t, i) -> dg_AdjSens_proto!(out, u, p, t, i, gu)
        calcCostGrad(dg!, iCond)
    end

    replace!(gradMatrix, NaN=>0.0)

    grad[:] = vec(sum(gradMatrix, dims=1))
    view(grad, doLogSearch) .*= view(allParameters, doLogSearch) * log(10)

    nothing
end




function adjointSensitivities(modelFunction, iStartPar, senseAlg, optAlg, 
        timeEnd, experimentalConditions, measurementData, observables, parameterBounds)

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

    # Initialize structs

    modelData = ModelData(new_sys, prob, observables, experimentalConditions, 
            initVariableNames, observableVariableNames, parameterInU0Names, parameterInObservableNames)

    experimentalData = ExperimentalData(observables, experimentalConditions, measurementData, modelData)

    modelParameters = ModelParameters(new_sys, prob, parameterBounds, experimentalConditions, measurementData, observables, experimentalData)

    dualModelParameters = DualModelParameters(modelParameters)

    parameterSpace, numAllStartParameters, lowerBounds, upperBounds = ParameterSpace(modelParameters, parameterBounds)

    modelOutput = ModelOutput(Float64, experimentalData, modelParameters)

    modelOutput_dual = ModelOutput(ForwardDiff.Dual, experimentalData, modelParameters)

    # Getting start values

    allStartParameters = getSamples(numAllStartParameters, lowerBounds, upperBounds, iStartPar)

    # Initialize functions

    solveODESystem = (iCond) -> solveODESystem_AdjSens_proto(prob, modelParameters.dynamicParametersVector, modelParameters.u0Vector, modelData, modelOutput, iCond)

    solveODESystem_dual = (dualPars, dualU0Vector, iCond) -> solveODESystem_AdjSens_proto(prob, dualPars, dualU0Vector, modelData, modelOutput_dual, iCond)

    g_unscaledObservables = (type, u, dynPar, i, iCond) -> g_unscaledObservables_AdjSens_proto(type, u, dynPar, i, iCond, modelData, experimentalData)

    g_scaledObservationFunctions = (type, h_bar, scaleVector, offsetVector, iCond) -> g_scaledObservationFunctions_AdjSens_proto(type, h_bar, scaleVector, offsetVector, iCond, modelParameters, experimentalData, modelData)

    g_cost = (type, h_hat, varianceVector, i, iCond) -> g_cost_AdjSens_proto(type, h_hat, varianceVector, i, iCond, modelParameters, experimentalData)

    g = (u, dynPar, scale, offset, variance, i, iCond; type = get_type([u, dynPar, scale, offset, variance])) -> g_AdjSens_proto(u, dynPar, scale, offset, variance, i, iCond, g_unscaledObservables, g_scaledObservationFunctions, g_cost; type)

    G = (iCond) -> G_AdjSens_proto(solveODESystem, g, iCond, modelParameters, experimentalData, modelOutput)

    G_specifiedDynPar = (specifiedDynPar, iCond) -> G_specifiedDynPar_AdjSens_proto(specifiedDynPar, modelData.parameterInU0Indices, g, solveODESystem_dual, iCond, modelParameters, dualModelParameters, experimentalData, modelOutput_dual)

    updateAllParameterVectors = () -> updateAllParameterVectors_proto(modelParameters, modelData)

    allConditionsCost = (p...) -> allConditionsCost_AdjSens_proto(parameterSpace, modelParameters, experimentalData, modelData, updateAllParameterVectors, G, p...)

    calcCostGrad = (dg!, iCond) -> calcCostGrad_AdjSens_proto(g, dg!, G_specifiedDynPar, iCond, modelParameters, modelData, experimentalData, modelOutput, senseAlg)

    allConditionsCostGrad = (grad, p...) -> allConditionsCostGrad_AdjSens_proto(parameterSpace, modelParameters, modelData, modelOutput, experimentalData, 
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
