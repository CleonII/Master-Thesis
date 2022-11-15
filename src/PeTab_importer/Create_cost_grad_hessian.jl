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
function setUpCostGradHess(peTabModel::PeTabModel, 
                           solver, 
                           tol::Float64; 
                           sparseJac::Bool=false, 
                           absTolSS::Float64=1e-8, 
                           relTolSS::Float64=1e-6)::PeTabOpt

    # Process PeTab files into type-stable Julia structs 
    experimentalConditionsFile, measurementDataFile, parameterDataFile, observablesDataFile = readDataFiles(peTabModel.dirModel, readObs=true)
    parameterData = processParameterData(parameterDataFile)
    measurementData = processMeasurementData(measurementDataFile, observablesDataFile) 
    simulationInfo = getSimulationInfo(measurementDataFile, measurementData, absTolSS=absTolSS, relTolSS=relTolSS)

    # Indices for mapping parameter-estimation vector to dynamic, observable and sd parameters correctly when calculating cost
    paramEstIndices = getIndicesParam(parameterData, measurementData)
    
    # Set model parameter values to those in the PeTab parameter data ensuring correct value of constant parameters 
    setParamToFileValues!(peTabModel.paramMap, peTabModel.stateMap, parameterData)

    # The time-span 5e3 is overwritten when performing actual forward simulations 
    odeProb = ODEProblem(peTabModel.odeSystem, peTabModel.stateMap, (0.0, 5e3), peTabModel.paramMap, jac=true, sparse=sparseJac)
    odeProb = remake(odeProb, p = convert.(Float64, odeProb.p), u0 = convert.(Float64, odeProb.u0))

    # Functions to map experimental conditions and parameters correctly to the ODE model 
    changeToExperimentalCondUse! = (pVec, u0Vec, expID) -> changeExperimentalCond!(pVec, u0Vec, expID, parameterData, experimentalConditionsFile, peTabModel)
    changeModelParamUse! = (pVec, u0Vec, paramEst, paramEstNames) -> changeModelParam!(pVec, u0Vec, paramEst, paramEstNames, paramEstIndices, peTabModel)

    # Set up function which solves the ODE model for all conditions and stores result 
    solveOdeModelAllCondUse! = (solArrayArg, odeProbArg) -> solveOdeModelAllExperimentalCond!(solArrayArg, odeProbArg, changeToExperimentalCondUse!, measurementDataFile, simulationInfo, solver, tol, onlySaveAtTobs=true)
    
    evalF = (paramVecEst) -> calcCost(paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!)
    evalGradF = (grad, paramVecEst) -> calcGradCost!(grad, paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!)
    evalHessApprox = (hessianMat, paramVecEst) -> calcHessianApprox!(hessianMat, paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!)
    # This is subtle. When computing the hessian via autodiff it is important that the ODE-solution arrary with dual 
    # numbers is used, else dual numbers will be present when computing the cost which will crash the code when taking 
    # the gradient of non-dynamic parameters in optim. 
    _evalHess = (paramVecEst) -> calcCost(paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, calcHessian=true)
    evalHess = (hessianMat, paramVec) -> begin hessianMat .= Symmetric(ForwardDiff.hessian(_evalHess, paramVec)) end
    
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
             solveOdeModelAllCondUse!::Function;
             calcHessian::Bool=false)

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
function calcCost(paramVecEst,
                  odeProb::ODEProblem,  
                  peTabModel::PeTabModel,
                  simulationInfo::SimulationInfo,
                  paramIndices::ParameterIndices,
                  measurementData::MeasurementData,
                  parameterData::ParamData,
                  changeModelParamUse!::Function,
                  solveOdeModelAllCondUse!::Function;
                  calcHessian::Bool=false)
                                                            

    # Correctly map paramVecEst to dynmaic, observable and sd param. The vectors, 
    # e.g dynamicParamEst, are distinct copies so transforming them will not change 
    # paramVecEst.
    dynamicParamEst = paramVecEst[paramIndices.iDynParam]
    obsParEst = paramVecEst[paramIndices.iObsParam]
    sdParamEst = paramVecEst[paramIndices.iSdParam]

    logLik = calcLogLikSolveODE(dynamicParamEst, sdParamEst, obsParEst, odeProb, peTabModel, simulationInfo, paramIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, calcHessDynParam=calcHessian)

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
    namesSdParam = paramIndices.namesSdParam
    namesObsParam = paramIndices.namesObsParam
    namesSdObsParam = paramIndices.namesSdObsParam

    # Calculate gradient seperately for dynamic and non dynamic parameter. 

    # I have tried to decrease run time here with chunking without success (deafult value performs best). Might be 
    # worth to look into a parellisation over the chunks (as for larger models each call takes relatively long time). 
    # Also parellisation of the chunks should be faster than paralellisation over experimental condtions.
    calcCostDyn = (x) -> calcLogLikSolveODE(x, sdParamEst, obsParEst, odeProb, peTabModel, simulationInfo, paramIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, calcGradDynParam=true)
    grad[paramIndices.iDynParam] .= ForwardDiff.gradient(calcCostDyn, dynamicParamEst)::Vector{Float64}

    # Here it is crucial to account for that obs- and sd parameter can be overlapping. Thus, a name-map 
    # of both Sd and Obs param is used to account for this 
    noneDynPar = paramVecEst[paramIndices.iSdObsPar]
    iSdUse = [findfirst(x -> x == namesSdParam[i], namesSdObsParam) for i in eachindex(namesSdParam)]
    iObsUse = [findfirst(x -> x == namesObsParam[i],  namesSdObsParam) for i in eachindex(namesObsParam)]

    # TODO : Make choice of gradient availble 
    calcCostNonDyn = (x) -> calcLogLikNotSolveODE(dynamicParamEst, x[iSdUse], x[iObsUse], peTabModel, simulationInfo, paramIndices, measurementData, parameterData, calcGradObsSdParam=true)
    @views ReverseDiff.gradient!(grad[paramIndices.iSdObsPar], calcCostNonDyn, noneDynPar)
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
    namesSdParam = paramIndices.namesSdParam
    namesObsParam = paramIndices.namesObsParam
    namesSdObsParam = paramIndices.namesSdObsParam

    # Calculate gradient seperately for dynamic and non dynamic parameter. 
    calcCostDyn = (x) -> calcLogLikSolveODE(x, sdParamEst, obsParEst, odeProb, peTabModel, simulationInfo, paramIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, calcHessDynParam=true)
    hessian[paramIndices.iDynParam, paramIndices.iDynParam] .= ForwardDiff.hessian(calcCostDyn, dynamicParamEst)::Matrix{Float64}

    # An ODE sol without dual numbers is needed here so ODE-system is resolved TODO: Make this the cost solution
    #=
    dynParamSolve = deepcopy(dynamicParamEst)
    transformParamVec!(dynParamSolve, paramIndices.namesDynParam, parameterData)
    odeProbUse = remake(odeProb, p = convert.(eltype(dynParamSolve), odeProb.p), u0 = convert.(eltype(dynParamSolve), odeProb.u0))
    changeModelParamUse!(odeProbUse.p, odeProbUse.u0, dynParamSolve, paramIndices.namesDynParam)
    solveOdeModelAllCondUse!(simulationInfo.solArray, odeProbUse)
    =#

    # Here it is crucial to account for that obs- and sd parameter can be overlapping. Thus, a name-map 
    # of both Sd and Obs param is used to account for this 
    noneDynPar::Array{Float64, 1} = paramVecEst[paramIndices.iSdObsPar]
    iSdUse = [findfirst(x -> x == namesSdParam[i], namesSdObsParam) for i in eachindex(namesSdParam)]
    iObsUse = [findfirst(x -> x == namesObsParam[i],  namesSdObsParam) for i in eachindex(namesObsParam)]
    
    # Compute hessian for none dynamic parameters 
    calcCostNonDyn = (x) -> calcLogLikNotSolveODE(dynamicParamEst, x[iSdUse], x[iObsUse], peTabModel, simulationInfo, paramIndices, measurementData, parameterData)
    @views ReverseDiff.hessian!(hessian[paramIndices.iSdObsPar, paramIndices.iSdObsPar], calcCostNonDyn, noneDynPar)

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
function calcLogLikSolveODE(dynamicParamEst,
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
                            calcHessDynParam::Bool=false, 
                            calcGradDynParam::Bool=false)::Real

    dynamicParamEstUse = transformParamVec(dynamicParamEst, paramIndices.namesDynParam, parameterData)
    sdParamEstUse = transformParamVec(sdParamEst, paramIndices.namesSdParam, parameterData)
    obsParEstUse = transformParamVec(obsParEst, paramIndices.namesObsParam, parameterData)

    odeProbUse = remake(odeProb, p = convert.(eltype(dynamicParamEstUse), odeProb.p), u0 = convert.(eltype(dynamicParamEstUse), odeProb.u0))
    changeModelParamUse!(odeProbUse.p, odeProbUse.u0, dynamicParamEstUse, paramIndices.namesDynParam)
    
    # If computing hessian or gradient store ODE solution in arrary with dual numbers, else use 
    # solution array with floats
    if calcHessDynParam == true || calcGradDynParam == true
        success = solveOdeModelAllCondUse!(simulationInfo.solArrayGrad, odeProbUse)
    else
        success = solveOdeModelAllCondUse!(simulationInfo.solArray, odeProbUse)
    end
    if success != true
        println("Failed to solve ODE model")
    end

    logLik = calcLogLik(dynamicParamEstUse, sdParamEstUse, obsParEstUse, peTabModel, simulationInfo, paramIndices, measurementData, parameterData, calcHessDynParam=calcHessDynParam, calcGradDynParam=calcGradDynParam)

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
                               parameterData::ParamData;
                               calcGradObsSdParam::Bool=false)::Real where T1<:Vector{<:Real}

    # To be able to use ReverseDiff sdParamEstUse and obsParamEstUse cannot be overwritten. 
    # Hence new vectors have to be created. Minimal overhead.
    dynamicParamEstUse = transformParamVec(dynamicParamEst, paramIndices.namesDynParam, parameterData)
    sdParamEstUse = transformParamVec(sdParamEst, paramIndices.namesSdParam, parameterData)
    obsParamEstUse = transformParamVec(obsParamEst, paramIndices.namesObsParam, parameterData)

    logLik = calcLogLik(dynamicParamEstUse, sdParamEstUse, obsParamEstUse, peTabModel, simulationInfo, paramIndices, measurementData, parameterData, calcGradObsSdParam=calcGradObsSdParam)
    
    return logLik
end


"""
    calcLogLik(dynamicParamEst::T1,
                    sdParamEst, 
                    obsPar, 
                    peTabModel::PeTabModel,
                    simulationInfo::SimulationInfo,
                    paramIndices::ParameterIndices,
                    measurementData::MeasurementData, 
                    parameterData::ParamData;
                    gradHessDynParam::Bool=false)::Real where T1<:Vector{<:Real}  

    Helper function computing the likelhood by given after solving the ODE-system using
    using the dynamic-parameters, sd-parameters and observable parameters. 

    Currently for Gaussian data log10 and non-transformed data is accepted.

    See also: [`calcCost`, `setUpCostGradHess`]
"""
function calcLogLik(dynamicParamEst,
                    sdParamEst, 
                    obsPar, 
                    peTabModel::PeTabModel,
                    simulationInfo::SimulationInfo,
                    paramIndices::ParameterIndices,
                    measurementData::MeasurementData, 
                    parameterData::ParamData;
                    calcHessDynParam::Bool=false, 
                    calcGradDynParam::Bool=false, 
                    calcGradObsSdParam::Bool=false)::Real 

    if calcHessDynParam == true || calcGradDynParam == true || calcGradObsSdParam == true
        odeSolArray = simulationInfo.solArrayGrad
    else
        odeSolArray = simulationInfo.solArray
    end

    # Compute yMod and sd-val by looping through all experimental conditons. At the end 
    # update the likelihood 
    logLik = 0.0
    for conditionID in keys(measurementData.iPerConditionId)
        # Extract the ODE-solution for specific condition ID
        whichForwardSol = findfirst(x -> x == conditionID, simulationInfo.conditionIdSol)
        odeSol = odeSolArray[whichForwardSol]   
        logLik += calcLogLikExpCond(odeSol, dynamicParamEst, sdParamEst, obsPar, 
                                    peTabModel, conditionID, paramIndices,
                                    measurementData, parameterData, calcGradObsSdParam)

        if isinf(logLik)
            return Inf
        end
    end

    return logLik
end


function calcLogLikExpCond(odeSol::ODESolution,
                           dynamicParamEst,
                           sdParamEst, 
                           obsPar, 
                           peTabModel::PeTabModel,
                           conditionID::String,
                           paramIndices::ParameterIndices,
                           measurementData::MeasurementData,
                           parameterData::ParamData, 
                           calcGradObsSdParam::Bool)::Real

    # Compute yMod and sd for all observations having id conditionID 
    logLik = 0.0
    for i in measurementData.iPerConditionId[conditionID]
        # Compute Y-mod value 
        t = measurementData.tObs[i]
        if calcGradObsSdParam == true
            odeSolAtT = dualVecToFloatVec(odeSol[:, measurementData.iTObs[i]])
        else
            odeSolAtT = odeSol[:, measurementData.iTObs[i]]
        end
        mapObsParam = paramIndices.mapArrayObsParam[paramIndices.indexObsParamMap[i]]
        yMod = peTabModel.evalYmod(odeSolAtT, t, dynamicParamEst, obsPar, parameterData, measurementData.observebleID[i], mapObsParam) 

        # Compute associated SD-value or extract said number if it is known 
        if isNumber(measurementData.sdParams[i])
            sdVal = parse(Float64, measurementData.sdParams[i])
        else
            mapSdParam = paramIndices.mapArraySdParam[paramIndices.indexSdParamMap[i]]
            sdVal = peTabModel.evalSd!(odeSolAtT, t, sdParamEst, dynamicParamEst, parameterData, measurementData.observebleID[i], mapSdParam)
        end

        # Transform yMod is necessary
        yModTrans = transformObsOrData(yMod, measurementData.transformData[i])
            
        # By default a positive ODE solution is not enforced (even though the user can provide it as option).
        # In case with transformations on the data the code can crash, hence Inf is returned in case the 
        # model data transformation can not be perfomred. 
        if isinf(yModTrans)
            return Inf
        end

        # Update log-likelihood 
        if measurementData.transformData[i] == :lin
            logLik += log(sdVal) + 0.5*log(2*pi) + 0.5*((yModTrans - measurementData.yObsNotTransformed[i]) / sdVal)^2
        elseif measurementData.transformData[i] == :log10
            logLik += log(sdVal) + 0.5*log(2*pi) + log(log(10)) + log(exp10(measurementData.yObsTransformed[i])) + 0.5*( ( log(exp10(yModTrans)) - log(exp10(measurementData.yObsTransformed[i])) ) / (log(10)*sdVal))^2
        else
            println("Transformation ", measurementData.transformData[i], "not yet supported.")
            return Inf
        end   
    end

    return logLik
end


"""
    changeModelParam!(paramVecOdeModel, 
                      stateVecOdeModel,
                      paramVecEst,
                      paramEstNames::Array{String, 1},
                      paramIndices::ParameterIndices,
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
                           paramIndices::ParameterIndices,
                           peTabModel::PeTabModel)

    # Allow the code to propegate dual numbers for gradients 
    paramMapUse = convert.(Pair{Num, eltype(paramVecOdeModel)}, peTabModel.paramMap)
    # TODO: Precompute this step 
    parameterNamesStr = string.([paramMapUse[i].first for i in eachindex(paramMapUse)])

    # Keep correct indices in ParamMap (for downstream computations when changing experimental conditions)
    # TODO : Check if this can be removed
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
    paramVecOdeModel[paramIndices.iMapDynParam] .= paramVecEst
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


function dualVecToFloatVec(dualVec::T1)::Vector{Float64} where T1<:Vector{<:ForwardDiff.Dual}
    return [dualVec[i].value for i in eachindex(dualVec)]
end