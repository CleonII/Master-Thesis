include(joinpath(pwd(), "src", "PeTab_importer", "Common.jl"))
include(joinpath(pwd(), "src", "PeTab_importer", "Map_parameters.jl"))
include(joinpath(pwd(), "src", "PeTab_importer", "Create_obs_u0_sd_functions.jl"))
include(joinpath(pwd(), "src", "PeTab_importer", "Process_PeTab_files.jl"))
include(joinpath(pwd(), "src", "Common.jl"))

"""
    setUpCostGradHess(peTabModel::PeTabModel, solver, tol::Float64)

    For a PeTab-model set up functions for computing i) the likelihood, ii) likelhood gradient, 
    and iii) likelhood Hessian block approximation. The functions are stored in PeTabOpt-struct 
    that can be used as input to the optimizers. 

    Currently the gradient for dynamic parameters (part of ODE-system) is computed via ForwardDiff, 
    and ReverseDiff is used for observable and sd parameters. The hessian approximation assumes the 
    interaction betweeen dynamic and (observable, sd) parameters is zero.
"""
function setUpCostGradHess(peTabModel::PeTabModel, solver, tol::Float64; sparseJac::Bool=false)::PeTabOpt

    # Process PeTab files into type-stable Julia structs 
    experimentalConditionsFile, measurementDataFile, parameterDataFile, observablesDataFile = readDataFiles(peTabModel.dirModel, readObs=true)
    parameterData = processParameterData(parameterDataFile)
    measurementData = processMeasurementData(measurementDataFile, observablesDataFile) 
    simulationInfo = getSimulationInfo(measurementDataFile)

    # Indices for mapping parameter-estimation vector to dynamic, observable and sd parameters correctly when calculating cost
    paramEstIndices = getIndicesParam(parameterData, measurementData)
    
    # Set model parameter values to those in the PeTab parameter data ensuring correct value of constant parameters 
    setParamToFileValues!(peTabModel.paramMap, peTabModel.stateMap, parameterData)

    # The time-span 5e3 is overwritten when performing actual forward simulations 
    odeProb = ODEProblem(peTabModel.odeSystem, peTabModel.stateMap, (0.0, 5e3), peTabModel.paramMap, jac=true, sparse=sparseJac)
    odeProb = remake(odeProb, p = convert.(Float64, odeProb.p), u0 = convert.(Float64, odeProb.u0))

    # Functions to map experimental conditions and parameters correctly to the ODE model 
    changeToExperimentalCondUse! = (pVec, u0Vec, expID) -> changeExperimentalCond!(pVec, u0Vec, expID, parameterData, experimentalConditionsFile, peTabModel)
    changeModelParamUse! = (pVec, u0Vec, paramEst, paramEstNames) -> changeModelParam!(pVec, u0Vec, paramEst, paramEstNames, peTabModel)

    # Set up function which solves the ODE model for all conditions and stores result 
    solveOdeModelAllCondUse! = (solArrayArg, odeProbArg) -> solveOdeModelAllExperimentalCond!(solArrayArg, odeProbArg, changeToExperimentalCondUse!, measurementDataFile, simulationInfo, solver, tol)
    
    evalF = (paramVecEst) -> calcCost(paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!)
    evalGradF = (grad, paramVecEst) -> calcGradCost!(grad, paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!)
    evalHessApprox = (hessianMat, paramVecEst) -> calcHessianApprox!(hessianMat, paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!)
    evalHess = (hessianMat, paramVec) -> begin hessianMat .= Symmetric(ForwardDiff.hessian(evalF, paramVec)) end
    
    # TODO: Move this to another function 
    # Lower and upper bounds for parameters to estimate 
    namesParamEst = paramEstIndices.namesParamEst
    lowerBounds = [parameterData.lowerBounds[findfirst(x -> x == namesParamEst[i], parameterData.parameterID)] for i in eachindex(namesParamEst)] 
    upperBounds = [parameterData.upperBounds[findfirst(x -> x == namesParamEst[i], parameterData.parameterID)] for i in eachindex(namesParamEst)] 
    # Parameter with nominal values in PeTab file 
    paramVecNominal = [parameterData.paramVal[findfirst(x -> x == namesParamEst[i], parameterData.parameterID)] for i in eachindex(namesParamEst)]

    # Transform upper and lower bounds if the case 
    transformParamVec!(lowerBounds, namesParamEst, parameterData, revTransform=true)
    transformParamVec!(upperBounds, namesParamEst, parameterData, revTransform=true)
    paramVecNominalTransformed = transformParamVec(paramVecNominal, namesParamEst, parameterData, revTransform=true)

    peTabOpt = PeTabOpt(evalF, 
                        evalGradF, 
                        evalHess,
                        evalHessApprox, 
                        length(namesParamEst), 
                        namesParamEst, 
                        paramVecNominal, 
                        paramVecNominalTransformed, 
                        lowerBounds, 
                        upperBounds, 
                        peTabModel.dirModel * "Cube" * peTabModel.modelName * ".csv",
                        peTabModel)
    return peTabOpt
end


"""
    calcCost(paramVecEst, 
             odeProb::ODEProblem,  
             peTabModel::PeTabModel,
             simulationInfo::SimulationInfo,
             paramIndices::ParameterIndices,
             measurementData::MeasurementData,
             parameterData::ParamData,
             changeModelParamUse!::Function,
             solveOdeModelAllCondUse!::Function)

    For a PeTab model compute the cost (likelhood) for a parameter vector 
    paramVecEst. With respect to paramVecEst (all other inputs fixed) 
    the function is compatible with ForwardDiff. 

    To compute the cost an ODE-problem, peTabModel, ODE simulation info, 
    indices to map parameter from paramVecEst, measurement data, parameter 
    data (e.g constant parameters), function to map parameters correctly to 
    ODE-model, and a function to solve the ODE model are required. These 
    are all set up correctly by the `setUpCostGradHess` function. 

    See also: [`setUpCostGradHess`]
"""
function calcCost(paramVecEst::T1,
                  odeProb::ODEProblem,  
                  peTabModel::PeTabModel,
                  simulationInfo::SimulationInfo,
                  paramIndices::ParameterIndices,
                  measurementData::MeasurementData,
                  parameterData::ParamData,
                  changeModelParamUse!::Function,
                  solveOdeModelAllCondUse!::Function) where T1<:Vector{<:Real}
                                                            

    # Correctly map paramVecEst to dynmaic, observable and sd param. The vectors, 
    # e.g dynamicParamEst, are distinct copies so transforming them will not change 
    # paramVecEst.
    dynamicParamEst = paramVecEst[paramIndices.iDynParam]
    obsParEst = paramVecEst[paramIndices.iObsParam]
    sdParamEst = paramVecEst[paramIndices.iSdParam]

    logLik = calcLogLikSolveODE(dynamicParamEst, sdParamEst, obsParEst, odeProb, peTabModel, simulationInfo, paramIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!)

    return logLik
end


"""
    calcGradCost!(grad::T1, 
                  paramVecEst, 
                  odeProb::ODEProblem,  
                  peTabModel::PeTabModel,
                  simulationInfo::SimulationInfo,
                  paramIndices::ParameterIndices,
                  measurementData::MeasurementData,
                  parameterData::ParamData, 
                  changeModelParamUse!::Function, 
                  solveOdeModelAllCondUse!::Function) where T1<:Array{<:AbstractFloat, 1}

    For a PeTab model compute inplace the gradient  of the cost (likelhood) for 
    a parameter vector paramVecEst. 

    Currently the gradient for dynamic parameters (part of ODE-system) is computed via ForwardDiff, 
    and ReverseDiff is used for observable and sd parameters. The input arguements are the same 
    as for `calcCost`, and everything is setup by `setUpCostGradHess` function.

    See also: [`setUpCostGradHess`]
"""
function calcGradCost!(grad::T1, 
                       paramVecEst::T2, 
                       odeProb::ODEProblem,  
                       peTabModel::PeTabModel,
                       simulationInfo::SimulationInfo,
                       paramIndices::ParameterIndices,
                       measurementData::MeasurementData,
                       parameterData::ParamData, 
                       changeModelParamUse!::Function, 
                       solveOdeModelAllCondUse!::Function) where {T1<:Array{<:AbstractFloat, 1}, 
                                                                  T2<:Vector{<:Real}}
    
    # Split input into observeble and dynamic parameters 
    dynamicParamEst = paramVecEst[paramIndices.iDynParam]
    obsParEst = paramVecEst[paramIndices.iObsParam]
    sdParamEst = paramVecEst[paramIndices.iSdParam]

    # Calculate gradient seperately for dynamic and non dynamic parameter. 

    # I have tried to decrease run time here with chunking without success (deafult value performs best). Might be 
    # worth to look into a parellisation over the chunks (as for larger models each call takes relatively long time). 
    # Also parellisation of the chunks should be faster than paralellisation over experimental condtions.
    calcCostDyn = (x) -> calcLogLikSolveODE(x, sdParamEst, obsParEst, odeProb, peTabModel, simulationInfo, paramIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, computeGradOrHess=true)
    grad[paramIndices.iDynParam] .= ForwardDiff.gradient(calcCostDyn, dynamicParamEst)::Vector{Float64}

    # Compute gradient for none dynamic parameters 
    noneDynPar = vcat(sdParamEst, obsParEst)
    iSdUse, iObsUse = 1:length(sdParamEst), (length(sdParamEst)+1):(length(sdParamEst) + length(obsParEst))
    # An ODE sol without dual numbers is needed here so ODE-system is resolved TODO: Make this the cost solution
    dynParamSolve = deepcopy(dynamicParamEst)
    transformParamVec!(dynParamSolve, paramIndices.namesDynParam, parameterData)
    odeProbUse = remake(odeProb, p = convert.(eltype(dynParamSolve), odeProb.p), u0 = convert.(eltype(dynParamSolve), odeProb.u0))
    changeModelParamUse!(odeProbUse.p, odeProbUse.u0, dynParamSolve, paramIndices.namesDynParam)
    solveOdeModelAllCondUse!(simulationInfo.solArray, odeProbUse)

    calcCostNonDyn = (x) -> calcLogLikNotSolveODE(dynamicParamEst, x[iSdUse], x[iObsUse], peTabModel, simulationInfo, paramIndices, measurementData, parameterData)
    # TODO : Make choice of gradient availble 
    @views ReverseDiff.gradient!(grad[vcat(paramIndices.iSdParam, paramIndices.iObsParam)], calcCostNonDyn, noneDynPar)
end


"""
    calcHessianApprox!(hessian::T1, 
                       paramVecEst, 
                       odeProb::ODEProblem,  
                       peTabModel::PeTabModel,
                       simulationInfo::SimulationInfo,
                       paramIndices::ParameterIndices,
                       measurementData::MeasurementData,
                       parameterData::ParamData, 
                       changeModelParamUse!::Function,
                       solveOdeModelAllCondUse!::Function) where T1<:Array{<:AbstractFloat, 2}

    For a PeTab model compute inplace hessian approximation of the cost (likelhood) for 
    a parameter vector paramVecEst. 

    The hessian approximation assumes the interaction betweeen dynamic and (observable, sd) parameters is zero. 
    The input arguements are the same as for `calcCost`, and everything is setup by `setUpCostGradHess` function.

    See also: [`setUpCostGradHess`]
"""
function calcHessianApprox!(hessian::T1, 
                            paramVecEst::T2,
                            odeProb::ODEProblem,  
                            peTabModel::PeTabModel,
                            simulationInfo::SimulationInfo,
                            paramIndices::ParameterIndices,
                            measurementData::MeasurementData,
                            parameterData::ParamData, 
                            changeModelParamUse!::Function,
                            solveOdeModelAllCondUse!::Function) where {T1<:Array{<:AbstractFloat, 2}, 
                                                                       T2<:Vector{<:Real}}

    # Avoid incorrect non-zero values 
    hessian .= 0.0

    dynamicParamEst = paramVecEst[paramIndices.iDynParam]
    obsParEst = paramVecEst[paramIndices.iObsParam]
    sdParamEst = paramVecEst[paramIndices.iSdParam]

    # Calculate gradient seperately for dynamic and non dynamic parameter. 
    calcCostDyn = (x) -> calcLogLikSolveODE(x, sdParamEst, obsParEst, odeProb, peTabModel, simulationInfo, paramIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, computeGradOrHess=true)
    @views hessian[paramIndices.iDynParam, paramIndices.iDynParam] .= ForwardDiff.hessian(calcCostDyn, dynamicParamEst)

    # Compute gradient for none dynamic parameters 
    noneDynPar = vcat(sdParamEst, obsParEst)
    iSdUse, iObsUse = 1:length(sdParamEst), (length(sdParamEst)+1):(length(sdParamEst) + length(obsParEst))
    # An ODE sol without dual numbers is needed here so ODE-system is resolved TODO: Make this the cost solution
    dynParamSolve = deepcopy(dynamicParamEst)
    transformParamVec!(dynParamSolve, paramIndices.namesDynParam, parameterData)
    odeProbUse = remake(odeProb, p = convert.(eltype(dynParamSolve), odeProb.p), u0 = convert.(eltype(dynParamSolve), odeProb.u0))
    changeModelParamUse!(odeProbUse.p, odeProbUse.u0, dynParamSolve, paramIndices.namesDynParam)
    solveOdeModelAllCondUse!(simulationInfo.solArray, odeProbUse)

    calcCostNonDyn = (x) -> calcLogLikNotSolveODE(dynamicParamEst, x[iSdUse], x[iObsUse], peTabModel, simulationInfo, paramIndices, measurementData, parameterData)
    @views hessian[vcat(paramIndices.iSdParam, paramIndices.iObsParam), vcat(paramIndices.iSdParam, paramIndices.iObsParam)] .= ForwardDiff.hessian(calcCostNonDyn, noneDynPar)
end


"""
    calcLogLikSolveODE(dynamicParamEst, 
                       sdParamEst, 
                       obsParEst, 
                       odeProb::ODEProblem,
                       peTabModel::PeTabModel,
                       simulationInfo::SimulationInfo,
                       measurementData ::MeasurementData, 
                       parameterData::ParamData, 
                       changeModelParamUse!::Function,
                       solveOdeModelAllCondUse!::Function, 
                       computeGradOrHess::Bool=false)

    Helper function computing the likelhood by solving the ODE system for all 
    PeTab-specifed experimental conditions using the dynamic-parameters, 
    sd-parameters and observable parameters. 

    When computing the cost and gradient/hessian for dynamic parameters the ODE 
    system must be solved before getting the likelhood. Besides the different 
    parameter vector the input arguements are the same as for `calcCost`, and 
    everything is setup by `setUpCostGradHess` function.

    See also: [`calcCost`, `setUpCostGradHess`]
"""
function calcLogLikSolveODE(dynamicParamEst::T1, 
                            sdParamEst,
                            obsParEst,
                            odeProb::ODEProblem,
                            peTabModel::PeTabModel,
                            simulationInfo::SimulationInfo,
                            paramIndices::ParameterIndices,
                            measurementData ::MeasurementData, 
                            parameterData::ParamData, 
                            changeModelParamUse!::Function,
                            solveOdeModelAllCondUse!::Function;
                            computeGradOrHess::Bool=false)::Real where T1<:Vector{<:Real}

    # If computing hessian or gradient the parameter vectors must not be overwritten, 
    # because here they are not copied as in the compute calcCost function (see first rows).
    if computeGradOrHess == true
        dynamicParamEst = transformParamVec(dynamicParamEst, paramIndices.namesDynParam, parameterData)
        sdParamEst = transformParamVec(sdParamEst, paramIndices.namesSdParam, parameterData)
        obsParEst = transformParamVec(obsParEst, paramIndices.namesObsParam, parameterData)

    # When computing the cost 
    else computeGradOrHess == false
        transformParamVec!(dynamicParamEst, paramIndices.namesDynParam, parameterData)
        transformParamVec!(sdParamEst, paramIndices.namesSdParam, parameterData)
        transformParamVec!(obsParEst, paramIndices.namesObsParam, parameterData)
    end

    odeProbUse = remake(odeProb, p = convert.(eltype(dynamicParamEst), odeProb.p), u0 = convert.(eltype(dynamicParamEst), odeProb.u0))
    changeModelParamUse!(odeProbUse.p, odeProbUse.u0, dynamicParamEst, paramIndices.namesDynParam)
    
    success = solveOdeModelAllCondUse!(simulationInfo.solArray, odeProbUse)
    if success != true
        return Inf
    end

    logLik = calcLogLik(dynamicParamEst, sdParamEst, obsParEst, peTabModel, simulationInfo, paramIndices, measurementData, parameterData)

    return logLik
end


"""
    calcLogLikNotSolveODE(dynamicParamEst, 
                          sdParamEst,
                          obsParamEst,
                          peTabModel::PeTabModel,
                          simulationInfo::SimulationInfo,
                          paramIndices::ParameterIndices,
                          measurementData::MeasurementData,
                          parameterData::ParamData)    

    Helper function computing the likelhood by given  an already existing ODE-solution stored 
    in simulationInfo using the dynamic-parameters, sd-parameters and observable parameters. 

    When computing the cost and gradient/hessian for the sd- and observable-parameters 
    only a solved ODE-system is needed (no need to resolve). This greatly reduces run-time 
    since a lot of dual numbers do not have to be propegated through the ODE solver.
    Besides the different parameter vector the input arguements are the same as for `calcCost`, 
    and everything is setup by `setUpCostGradHess` function.

    See also: [`calcCost`, `setUpCostGradHess`]
"""
function calcLogLikNotSolveODE(dynamicParamEst::T1, 
                               sdParamEst,
                               obsParamEst,
                               peTabModel::PeTabModel,
                               simulationInfo::SimulationInfo,
                               paramIndices::ParameterIndices,
                               measurementData::MeasurementData,
                               parameterData::ParamData)::Real where T1<:Vector{<:Real}

    # To be able to use ReverseDiff sdParamEstUse and obsParamEstUse cannot be overwritten. 
    # Hence new vectors have to be created. Minimal overhead.
    dynamicParamEstUse = transformParamVec(dynamicParamEst, paramIndices.namesDynParam, parameterData)
    sdParamEstUse = transformParamVec(sdParamEst, paramIndices.namesSdParam, parameterData)
    obsParamEstUse = transformParamVec(obsParamEst, paramIndices.namesObsParam, parameterData)

    logLik = calcLogLik(dynamicParamEstUse, sdParamEstUse, obsParamEstUse, peTabModel, simulationInfo, paramIndices, measurementData, parameterData)
    
    return logLik
end


"""
    calcLogLikNotSolveODE(dynamicParamEst, 
                          sdParamEst,
                          obsParamEst,
                          peTabModel::PeTabModel,
                          simulationInfo::SimulationInfo,
                          paramIndices::ParameterIndices,
                          measurementData::MeasurementData,
                          parameterData::ParamData)    

    Helper function computing the likelhood by given after solving the ODE-system using
    using the dynamic-parameters, sd-parameters and observable parameters. 

    Currently for Gaussian data log10 and non-transformed data is accepted.

    See also: [`calcCost`, `setUpCostGradHess`]
"""
function calcLogLik(dynamicParamEst::T1,
                    sdParamEst, 
                    obsPar, 
                    peTabModel::PeTabModel,
                    simulationInfo::SimulationInfo,
                    paramIndices::ParameterIndices,
                    measurementData::MeasurementData, 
                    parameterData::ParamData)::Real where T1<:Vector{<:Real}
                                                      
    # Ensure correct type of yMod-array (dual only if obs-par or ODE solution are dual)
    if !isempty(obsPar) && !(typeof(obsPar[1]) <: AbstractFloat)
        yMod = Array{eltype(obsPar), 1}(undef, length(measurementData.yObsNotTransformed))
    else
        yMod = Array{eltype(dynamicParamEst), 1}(undef, length(measurementData.yObsNotTransformed))
    end

    # Compute yMod value by looping through all observables 
    @inbounds for i in eachindex(yMod)
        whichForwardSol = findfirst(x -> x == measurementData.conditionId[i], simulationInfo.conditionIdSol)
        t = measurementData.tObs[i]
        odeSol = simulationInfo.solArray[whichForwardSol](t)
        mapObsParam = paramIndices.mapArrayObsParam[paramIndices.indexObsParamMap[i]]
        yMod[i] = peTabModel.evalYmod(odeSol, t, dynamicParamEst, obsPar, parameterData, measurementData.observebleID[i], mapObsParam) 
    end

    # Transform model output if required
    transformYobsOrYmodArr!(yMod, measurementData.transformData)

    # Calculate (or extract) standard deviation for each observation 
    sdVal = Array{eltype(sdParamEst), 1}(undef, length(measurementData.yObsNotTransformed))
    @inbounds for i in eachindex(sdVal)
        # Known SD-value 
        if isNumber(measurementData.sdParams[i])
            sdVal[i] = parse(Float64, measurementData.sdParams[i])
            continue
        end
        # Compute SD in case not known 
        whichForwardSol = findfirst(x -> x == measurementData.conditionId[i], simulationInfo.conditionIdSol)
        t = measurementData.tObs[i]
        odeSol = simulationInfo.solArray[whichForwardSol](t)
        mapSdParam = paramIndices.mapArraySdParam[paramIndices.indexSdParamMap[i]]
        sdVal[i] = peTabModel.evalSd!(odeSol, t, sdParamEst, dynamicParamEst, parameterData, measurementData.observebleID[i], mapSdParam)
    end

    logLik = 0.0
    @inbounds for i in eachindex(sdVal)
        if measurementData.transformData[i] == :lin
            logLik += log(sdVal[i]) + 0.5*log(2*pi) + 0.5*((yMod[i] - measurementData.yObsNotTransformed[i]) / sdVal[i])^2
        elseif measurementData.transformData[i] == :log10
            logLik += log(sdVal[i]) + 0.5*log(2*pi) + log(log(10)) + log(exp10(measurementData.yObsTransformed[i])) + 0.5*( ( log(exp10(yMod[i])) - log(exp10(measurementData.yObsTransformed[i])) ) / (log(10)*sdVal[i]))^2
            #logLik += log(sdVal[i]) + 0.5*log(2*pi) + log(log(10)) + 0.5* ((yMod[i] - measurementData.yObsTransformed[i]) / sdVal[i])^2
        else
            println("Transformation ", measurementData.transformData[i], "not yet supported.")
        end
    end

    return logLik
end


"""
    changeModelParam!(paramVecOdeModel, 
                      stateVecOdeModel,
                      paramVecEst,
                      paramEstNames::Array{String, 1},
                      peTabModel::PeTabModel)

    Change the ODE parameter vector (paramVecOdeModel) and initial value vector (stateVecOdeModel)
    values to the values in parameter vector used for parameter estimation paramVecEst. 
    Basically, map the parameter-estiamtion vector to the ODE model.

    The function can handle that paramVecEst is a Float64 vector or a vector of Duals for the 
    gradient calculations. This function is used when computing the cost, and everything 
    is set up by `setUpCostGradHess`. 
     
    See also: [`setUpCostGradHess`]
"""
function changeModelParam!(paramVecOdeModel, 
                           stateVecOdeModel,
                           paramVecEst,
                           paramEstNames::Array{String, 1},
                           peTabModel::PeTabModel)

    # Allow the code to propegate dual numbers for gradients 
    paramMapUse = convert.(Pair{Num, eltype(paramVecOdeModel)}, peTabModel.paramMap)
    # TODO: Precompute this step 
    parameterNamesStr = string.([paramMapUse[i].first for i in eachindex(paramMapUse)])

    # Change parameters (and states) to current iterations in parameter estimation vector 
    for i in eachindex(paramEstNames)
        
        paramChangeName = paramEstNames[i]
        valChangeTo = paramVecEst[i]
        # Propegate dual numbers correctly. In the reference paramMap want to have floats (not duals)
        valChangeToFloat = typeof(valChangeTo) <: Union{ForwardDiff.Dual, ForwardDiff.Dual{<:ForwardDiff.Dual}} ? valChangeTo.value : valChangeTo
        
        i_param = findfirst(x -> x == paramChangeName, parameterNamesStr)
        if !isnothing(i_param)
            paramMapUse[i_param] = Pair(paramMapUse[i_param].first, valChangeTo) 
            # Update reference paramMap 
            if typeof(valChangeToFloat) <: AbstractFloat
                peTabModel.paramMap[i_param] = Pair(peTabModel.paramMap[i_param].first, valChangeToFloat) 
            end
        else
            println("Error : Simulation parameter to change not found for experimental condition $expID")
        end
    end

    # Use ModellingToolkit and u0 function to correctly map parameters to ODE-system 
    newVal = ModelingToolkit.varmap_to_vars(paramMapUse, peTabModel.paramNames)
    paramVecOdeModel .= newVal
    peTabModel.evalU0!(stateVecOdeModel, paramVecOdeModel) 
    
    return nothing
end


"""
    transformParamVec!(paramVec, namesParam::Array{String, 1}, paramData::ParamData; revTransform::Bool=false)

    Helper function which transforms in-place a parameter vector with parameters specied in namesParam according to the 
    transformation for said parameter specifid in paramData.shouldTransform. In case revTransform is true 
    performs the inverse parameter transformation (e.g exp10 instead of log10)
"""
function transformParamVec!(paramVec, 
                            namesParam::Array{String, 1}, 
                            paramData::ParamData; 
                            revTransform::Bool=false) 
    
    @inbounds for i in eachindex(paramVec)
        iParam = findfirst(x -> x == namesParam[i], paramData.parameterID)
        if isnothing(iParam)
            println("Warning : Could not find paramID for $namesParam")
        end
        if paramData.logScale[iParam] == true && revTransform == false
            paramVec[i] = exp10(paramVec[i])
        elseif paramData.logScale[iParam] == true && revTransform == true
            paramVec[i] = log10(paramVec[i])
        end
    end
end


"""
    transformParamVec!(paramVec, namesParam::Array{String, 1}, paramData::ParamData; revTransform::Bool=false)

    Helper function which returns a transformed parameter vector with parameters specied in namesParam according to the 
    transformation for said parameter specifid in paramData.shouldTransform. In case revTransform is true 
    performs the inverse parameter transformation (e.g exp10 instead of log10)
"""
function transformParamVec(paramVec, 
                           namesParam::Array{String, 1}, 
                           paramData::ParamData; 
                           revTransform::Bool=false)
    
    shouldTransform = zeros(Float64, length(paramVec))
    shouldNotTransform = ones(Float64, length(paramVec))
    @inbounds for i in eachindex(paramVec)
        iParam = findfirst(x -> x == namesParam[i], paramData.parameterID)
        if isnothing(iParam)
            println("Warning : Could not find paramID for $namesParam")
        end
        if paramData.logScale[iParam]
            shouldTransform[i] = 1.0
            shouldNotTransform[i] = 0.0
        end
    end

    if revTransform == false
        return exp10.(paramVec) .* shouldTransform .+ paramVec .* shouldNotTransform
    else
        return log10.(paramVec) .* shouldTransform .+ paramVec .* shouldNotTransform
    end
end