function checkForPeTabFile(fileSearchFor::String, dirModel::String)

    filesDirModel = readdir(dirModel)
    iUse = findall(x -> occursin(fileSearchFor, x), filesDirModel)
    if length(iUse) > 1 
        @printf("Error : More than 1 file starting with %s in %s\n", fileSearchFor, filesDirModel)
    end
    if length(iUse) == 0
        @printf("Error : No file starting with %s in %s\n", fileSearchFor, filesDirModel)
    end

    return dirModel * filesDirModel[iUse[1]]
end


"""
    setUpPeTabModel(modelName::String, dirModel::String)::PeTabModel

    Given a model directory (dirModel) containing the PeTab files and a 
    xml-file on format modelName.xml will return a PeTabModel struct holding 
    paths to PeTab files, ode-system in ModellingToolkit format, functions for 
    evaluating yMod, u0 and standard deviations, and a parameter and state maps 
    for how parameters and states are mapped in the ModellingToolkit ODE system
    along with state and parameter names.

    dirModel must contain a SBML file named modelName.xml, and files starting with 
    measurementData, experimentalCondition, parameter, and observables (tsv-files).
    The latter files must be unique (e.g only one file starting with measurementData)

    TODO : Example  
"""
function setUpPeTabModel(modelName::String, dirModel::String)::PeTabModel

    # Sanity check user input 
    modelFileXml = dirModel * modelName * ".xml"
    modelFileJl = dirModel * modelName * ".jl"
    if !isdir(dirModel)
        @printf("Model directory %s does not exist\n", dirModel)
    end
    if !isfile(modelFileXml)
        @printf("Model directory does not contain xml-file with name %s\n", modelName * "xml")
    end
    # If Julia model file does exists build it 
    if !isfile(modelFileJl)
        @printf("Julia model file does not exist - will build it\n")
        XmlToModellingToolkit(modelFileXml, modelName, dirModel)
    else
        @printf("Julia model file exists at %s - will not rebuild it\n", modelFileJl)
    end

    # Extract ODE-system and mapping of maps of how to map parameters to states and model parmaeters 
    include(modelFileJl)
    expr = Expr(:call, Symbol("getODEModel_" * modelName))
    odeSys, stateMap, paramMap = eval(expr)
    odeSysUse = ode_order_lowering(odeSys)
    # Parameter and state names for ODE-system 
    parameterNames = parameters(odeSysUse)
    stateNames = states(odeSysUse)

    # Sanity check for presence of all PeTab-files 
    pathMeasurementData = checkForPeTabFile("measurementData", dirModel)
    pathExperimentalCond = checkForPeTabFile("experimentalCondition", dirModel)
    pathParameters = checkForPeTabFile("parameters", dirModel)
    pathObservables = checkForPeTabFile("observables", dirModel)

    # Build functions for observables, sd and u0 if does not exist and include
    pathObsSdU0 = dirModel * modelName * "ObsSdU0.jl"
    if !isfile(pathObsSdU0)
        @printf("File for yMod, U0 and Sd does not exist - building it\n")
        createFilesModel(modelName, dirModel, odeSysUse, paramMap, stateMap)
    else
        @printf("File for yMod, U0 and Sd does exist - will not rebuild it\n")
    end
    include(pathObsSdU0)

    peTabModel = PeTabModel(modelName,
                            evalYmod,
                            evalU0!,
                            evalSd!,
                            odeSysUse,
                            paramMap,
                            stateMap,
                            parameterNames, 
                            stateNames,
                            dirModel,
                            pathMeasurementData,
                            pathExperimentalCond,
                            pathObservables, 
                            pathParameters)

    return peTabModel
end


# Transform a parameterVec with corresponding parameterId based on the reference paramData. 
# Currently supports log-scale but extending to other transforms is straighforward. 
function transformVector!(parameterVec, parameterId::Array{String, 1}, paramData::ParamData; revTransform::Bool=false)
    for i in eachindex(parameterVec)
        iParam = findfirst(x -> x == parameterId[i], paramData.parameterID)
        if isnothing(iParam)
            println("Warning : Could not find paramID for $parameterId")
        end
        if paramData.logScale[iParam] == true && revTransform == false
            parameterVec[i] = exp10(parameterVec[i])
        elseif paramData.logScale[iParam] == true && revTransform == true
            parameterVec[i] = log10(parameterVec[i])
        end
    end
end
function transformVector(parameterVec, parameterId::Array{String, 1}, paramData::ParamData; revTransform::Bool=false)
    
    shouldTransform = zeros(Float64, length(parameterVec))
    shouldNotTransform = ones(Float64, length(parameterVec))
    for i in eachindex(parameterVec)
        iParam = findfirst(x -> x == parameterId[i], paramData.parameterID)
        if isnothing(iParam)
            println("Warning : Could not find paramID for $parameterId")
        end
        if paramData.logScale[iParam]
            shouldTransform[i] = 1.0
            shouldNotTransform[i] = 0.0
        end
    end

    if revTransform == false
        return exp10.(parameterVec) .* shouldTransform .+ parameterVec .* shouldNotTransform
    else
        return log10.(parameterVec) .* shouldTransform .+ parameterVec .* shouldNotTransform
    end
end


# Set up cost function, gradient function, and approximate hessian function. 
# An exact hessian up to numerical precision is obtained by using ForwardDiff
# on the cost function (as the code can propagate dual numbers)
function setUpCostFunc(peTabModel::PeTabModel, solver, tol::Float64)

    # Read data on experimental conditions and parameter values and process into Julia structs 
    experimentalConditionsFile, measurementDataFile, parameterDataFile, observablesDataFile = readDataFiles(peTabModel.dirModel, readObs=true)
    parameterData = processParameterData(parameterDataFile)
    measurementData = processObsData(measurementDataFile, observablesDataFile) 

    # Experimental ID:s and array for storing forward simulations 
    simulationInfo = getSimulationInfo(measurementDataFile, peTabModel.odeSystem)
    
    # Set model parameters to those in the parameter PeTab-file 
    setParamToParamFileVal!(peTabModel.paramMap, peTabModel.stateMap, parameterData)

    # Build ODE-prob from ODE-system. Note 5e3 time-span is overwritten in the forward simulations 
    odeProb = ODEProblem(peTabModel.odeSystem, peTabModel.stateMap, (0.0, 5e3), peTabModel.paramMap, jac=true)
    odeProb = remake(odeProb, p = convert.(Float64, odeProb.p), u0 = convert.(Float64, odeProb.u0))

    # Functions to map experimental conditions and parameter to estimate correctly to the ODE-model 
    changeToExperimentalCondUse! = (pVec, u0Vec, expID) -> changeExperimentalCond!(pVec, u0Vec, expID, parameterData, experimentalConditionsFile, peTabModel)
    changeModelParamUse! = (pVec, u0Vec, paramEst, paramEstNames) -> changeModelParam!(pVec, u0Vec, paramEst, paramEstNames, parameterData, peTabModel)

    # Set up function which solves the ODE model for all conditions and stores result 
    solveOdeModelAllCondUse! = (solArrayArg, odeProbArg) -> solveOdeModelAllCond!(solArrayArg, odeProbArg, changeToExperimentalCondUse!, measurementDataFile, simulationInfo, solver, tol)

    # Indices for mapping parameter-estimation vector to dynamic, observable and sd parameters correctly when calculating cost
    paramEstIndices = getIndicesParam(parameterData, measurementData)

    # Lower and upper bounds for parameters to estimate 
    namesParamEst = paramEstIndices.namesParamEst
    lowerBounds = [parameterData.lowerBounds[findfirst(x -> x == namesParamEst[i], parameterData.parameterID)] for i in eachindex(namesParamEst)] 
    upperBounds = [parameterData.upperBounds[findfirst(x -> x == namesParamEst[i], parameterData.parameterID)] for i in eachindex(namesParamEst)] 
    
    # Set up cost funciton and derivative functions 
    evalF = (paramVecEst) -> calcCost(paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!)
    evalGradF = (paramVecEst, grad) -> calcGradCost!(grad, paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!)
    evalHessianApproxF = (hessian, paramVecEst) -> calcHessianApprox!(hessian, paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!)

    # Parameter with nominal values in PeTab file 
    paramVecNominal = [parameterData.paramVal[findfirst(x -> x == namesParamEst[i], parameterData.parameterID)] for i in eachindex(namesParamEst)]
    
    # Transform upper and lower bounds if the case 
    transformVector!(lowerBounds, namesParamEst, parameterData, revTransform=true)
    transformVector!(upperBounds, namesParamEst, parameterData, revTransform=true)
    transformVector!(paramVecNominal, namesParamEst, parameterData, revTransform=true)

    return evalF, evalGradF, evalHessianApproxF, paramVecNominal, lowerBounds, upperBounds, namesParamEst
end


# Calculate the likelhood for a given model by solving the ODE system for all experimental 
# conditions. Exact hessian can be obtained by using ForwardDiff on this function. 
function calcCost(paramVecEst, 
                  odeProb::ODEProblem,  
                  peTabModel::PeTabModel,
                  simulationInfo::SimulationInfo,
                  paramIndices::ParameterIndices,
                  measurementData::ObservedData,
                  parameterData::ParamData,
                  changeModelParamUse!::Function,
                  solveOdeModelAllCondUse!::Function)

    # Split input into observeble and dynamic parameters 
    dynamicParamEst = paramVecEst[paramIndices.iDynParam]
    obsParEst = paramVecEst[paramIndices.iObsParam]
    sdParamEst = paramVecEst[paramIndices.iSdParam]

    # Correctly transform parameters 
    transformVector!(dynamicParamEst, paramIndices.namesDynParam, parameterData)
    transformVector!(sdParamEst, paramIndices.namesSdParam, parameterData)
    transformVector!(obsParEst, paramIndices.namesObsParam, parameterData)

    # Convert ODE-problem to propegate correct numbers (e.g dual)
    odeProbUse = remake(odeProb, p = convert.(eltype(paramVecEst), odeProb.p), u0 = convert.(eltype(paramVecEst), odeProb.u0))
    # Correctly map parameters to ODE-system
    changeModelParamUse!(odeProbUse.p, odeProbUse.u0, dynamicParamEst, paramIndices.namesDynParam)

    logLik = calcLogLikSolveODE(dynamicParamEst, sdParamEst, obsParEst, odeProbUse, peTabModel, simulationInfo, measurementData, parameterData, solveOdeModelAllCondUse!)

    return logLik
end


# Calculate the log-likelihood given the ODE solution with correct transformation for the model output. 
# Need to add support for scale and offset parameters. 
function calcLogLik(dynamicParamEst, 
                    sdParamEst, 
                    obsPar, 
                    peTabModel::PeTabModel,
                    simulationInfo::SimulationInfo,
                    measurementData::ObservedData, 
                    parameterData::ParamData)

    # Ensure correct type of yMod-array (dual only if obs-par or ODE solution are dual)
    if !isempty(obsPar) && !(typeof(obsPar[1]) <: AbstractFloat)
        yMod = Array{eltype(obsPar), 1}(undef, length(measurementData.yObs))
    else
        yMod = Array{eltype(dynamicParamEst), 1}(undef, length(measurementData.yObs))
    end
    # Compute yMod value by looping through all observables 
    @inbounds for i in eachindex(yMod)
        whichForwardSol = findfirst(x -> x == measurementData.conditionId[i], simulationInfo.conditionIdSol)
        t = measurementData.tObs[i]
        odeSol = simulationInfo.solArray[whichForwardSol](t)
        yMod[i] = peTabModel.evalYmod(odeSol, t, dynamicParamEst, obsPar, parameterData, measurementData, measurementData.observebleDd[i], measurementData.conditionId[i]) 
    end

    # Transform model output if required
    transformObsOrDataArr!(yMod, measurementData.transformData)

    # Calculate (or extract) standard deviation for each observation 
    sdVal = Array{eltype(sdParamEst), 1}(undef, length(measurementData.yObs))
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
        sdVal[i] = peTabModel.evalSd!(odeSol, t, sdParamEst, dynamicParamEst, parameterData, measurementData, measurementData.observebleDd[i], measurementData.conditionId[i])
    end

    logLik = 0.0
    @inbounds for i in eachindex(sdVal)
        if measurementData.transformData[i] == :lin
            logLik += log(sdVal[i]) + 0.5*log(2*pi) + 0.5*((yMod[i] - measurementData.yObs[i]) / sdVal[i])^2
        elseif measurementData.transformData[i] == :log10
            logLik += log(sdVal[i]) + 0.5*log(2*pi) + log(log(10)) + log(exp10(measurementData.yObs[i])) + 0.5*( ( log(exp10(yMod[i])) - log(exp10(measurementData.yObs[i])) ) / (log(10)*sdVal[i]))^2
        else
            println("Transformation ", measurementData.transformData[i], "not yet supported.")
        end
    end

    return logLik
end


# Calculate the log-likelihood by solving the ODE-problem for all experimental conditions. 
function calcLogLikSolveODE(dynamicParamEst, 
                            sdParamEst, 
                            obsParEst, 
                            odeProb::ODEProblem,
                            peTabModel::PeTabModel,
                            simulationInfo::SimulationInfo,
                            measurementData ::ObservedData, 
                            parameterData::ParamData, 
                            solveOdeModelAllCondUse!::Function)


    # In place solving (store ODE-solution in solArray)
    success = solveOdeModelAllCondUse!(simulationInfo.solArray, odeProb)
    if success != true
        return Inf
    end

    logLik = calcLogLik(dynamicParamEst, sdParamEst, obsParEst, peTabModel, simulationInfo, measurementData, parameterData)

    return logLik
end


# Calculate the gradient of the likelihood in two steps. First calculate for the parameters which are 
# a part of the ODE-system. This is expansive as dual numbers have to be propagated through the ODE system. 
# Secondly, calculate the gradient for the observeble parameters given the ODE solution, thus avoiding 
# the step of propegating dual numbers through the ODE solver. 
function calcGradCost!(grad::T1, 
                       paramVecEst, 
                       odeProb::ODEProblem,  
                       peTabModel::PeTabModel,
                       simulationInfo::SimulationInfo,
                       paramIndices::ParameterIndices,
                       measurementData::ObservedData,
                       parameterData::ParamData, 
                       changeModelParamUse!::Function, 
                       solveOdeModelAllCondUse!::Function) where T1<:Array{<:AbstractFloat, 1}
    
    # Split input into observeble and dynamic parameters 
    dynamicParamEst = paramVecEst[paramIndices.iDynParam]
    obsParEst = paramVecEst[paramIndices.iObsParam]
    sdParamEst = paramVecEst[paramIndices.iSdParam]

    # Calculate gradient seperately for dynamic and non dynamic parameter. 
    calcCostDyn = (x) -> calcLogLikSolveOdeGrad(x, sdParamEst, obsParEst, odeProb, peTabModel, simulationInfo, paramIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!)
    grad[paramIndices.iDynParam] .= ForwardDiff.gradient(calcCostDyn, dynamicParamEst)

    # Compute gradient for none dynamic parameters 
    noneDynPar = vcat(sdParamEst, obsParEst)
    iSdUse, iObsUse = 1:length(sdParamEst), (length(sdParamEst)+1):(length(sdParamEst) + length(obsParEst))
    # An ODE sol without dual numbers is needed here so ODE-system is resolved TODO: Make this the cost solution
    dynParamSolve = deepcopy(dynamicParamEst)
    transformVector!(dynParamSolve, paramIndices.namesDynParam, parameterData)
    odeProbUse = remake(odeProb, p = convert.(eltype(dynParamSolve), odeProb.p), u0 = convert.(eltype(dynParamSolve), odeProb.u0))
    changeModelParamUse!(odeProbUse.p, odeProbUse.u0, dynParamSolve, paramIndices.namesDynParam)
    solveOdeModelAllCondUse!(simulationInfo.solArray, odeProbUse)

    calcCostNonDyn = (x) -> calcLogLikNotSolveODE(dynamicParamEst, x[iSdUse], x[iObsUse], peTabModel, simulationInfo, paramIndices, measurementData, parameterData)
    # TODO : Make choice of gradient availble 
    @views ReverseDiff.gradient!(grad[vcat(paramIndices.iSdParam, paramIndices.iObsParam)], calcCostNonDyn, noneDynPar)

end


# Calculate the log-likelihood by solving the ODE system. This function is used for gradient calculations, 
# and this function is used to obtain the gradient for the dynamic parameters. 
function calcLogLikSolveOdeGrad(dynamicParamEst, 
                                sdParamEst, 
                                obsPar, 
                                odeProb::ODEProblem,
                                peTabModel::PeTabModel, 
                                simulationInfo::SimulationInfo, 
                                paramIndices::ParameterIndices, 
                                measurementData::ObservedData, 
                                parameterData::ParamData, 
                                changeModelParamUse!::Function, 
                                solveOdeModelAllCondUse!::Function)

    # Prevent transformation from overwriting other parts of calculations 
    dynamicParamEstUse = dynamicParamEst[:]
    sdParamEstUse = sdParamEst[:]
    obsParamEstUse = obsPar[:]

    # Correctly transform parameters within AD part of code 
    transformVector!(dynamicParamEstUse, paramIndices.namesDynParam, parameterData)
    transformVector!(sdParamEstUse, paramIndices.namesSdParam, parameterData)
    transformVector!(obsParamEstUse, paramIndices.namesObsParam, parameterData)

    # Function solving ODE and calculating likelihood. 
    odeProbUse = remake(odeProb, p = convert.(eltype(dynamicParamEstUse), odeProb.p), u0 = convert.(eltype(dynamicParamEstUse), odeProb.u0))
    # Map dynamic parameters to ODE-problem via varmap n modelling-toolkit 
    changeModelParamUse!(odeProbUse.p, odeProbUse.u0, dynamicParamEstUse, paramIndices.namesDynParam)
    
    logLik = calcLogLikSolveODE(dynamicParamEstUse, sdParamEstUse, obsParamEstUse, odeProbUse, peTabModel, simulationInfo, measurementData, parameterData, solveOdeModelAllCondUse!)

    return logLik
end


# Calculate the log-likelihood given the solution of the ODE system. This function is used for gradient calculations, 
# and this function is used to obtain the gradient for the observable parameters. 
function calcLogLikNotSolveODE(dynamicParamEst, 
                               sdParamEst,
                               obsParamEst,
                               peTabModel::PeTabModel,
                               simulationInfo::SimulationInfo,
                               paramIndices::ParameterIndices,
                               measurementData::ObservedData,
                               parameterData::ParamData)    

    # Avoid having transformation overwriting other part of code 
    dynamicParamEstUse = dynamicParamEst[:]    

    # To be able to use ReverseDiff sdParamEstUse and obsParamEstUse cannot be overwritten. 
    # Hence new vectors have to be created. Minimal overhead.
    transformVector!(dynamicParamEstUse, paramIndices.namesDynParam, parameterData)
    sdParamEstUse = transformVector(sdParamEst, paramIndices.namesSdParam, parameterData)
    obsParamEstUse = transformVector(obsParamEst, paramIndices.namesObsParam, parameterData)

    logLik = calcLogLik(dynamicParamEstUse, sdParamEstUse, obsParamEstUse, peTabModel, simulationInfo, measurementData, parameterData)
    
    return logLik
end


# Calculate an approximation of the hessian where the interaction terms between observeble and 
# dynamic parameters is assumed to be zero. 
function calcHessianApprox!(hessian::T1, 
                            paramVecEst, 
                            odeProb::ODEProblem,  
                            peTabModel::PeTabModel,
                            simulationInfo::SimulationInfo,
                            paramIndices::ParameterIndices,
                            measurementData::ObservedData,
                            parameterData::ParamData, 
                            changeModelParamUse!::Function,
                            solveOdeModelAllCondUse!::Function) where T1<:Array{<:AbstractFloat, 2}

    # Avoid incorrect non-zero values 
    hessian .= 0.0

    dynamicParamEst = paramVecEst[paramIndices.iDynParam]
    obsParEst = paramVecEst[paramIndices.iObsParam]
    sdParamEst = paramVecEst[paramIndices.iSdParam]

    # Calculate gradient seperately for dynamic and non dynamic parameter. 
    calcCostDyn = (x) -> calcLogLikSolveOdeGrad(x, sdParamEst, obsParEst, odeProb, peTabModel, simulationInfo, paramIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!)
    @views hessian[paramIndices.iDynParam, paramIndices.iDynParam] .= ForwardDiff.hessian(calcCostDyn, dynamicParamEst)

    # Compute gradient for none dynamic parameters 
    noneDynPar = vcat(sdParamEst, obsParEst)
    iSdUse, iObsUse = 1:length(sdParamEst), (length(sdParamEst)+1):(length(sdParamEst) + length(obsParEst))
    # An ODE sol without dual numbers is needed here so ODE-system is resolved TODO: Make this the cost solution
    dynParamSolve = deepcopy(dynamicParamEst)
    transformVector!(dynParamSolve, paramIndices.namesDynParam, parameterData)
    odeProbUse = remake(odeProb, p = convert.(eltype(dynParamSolve), odeProb.p), u0 = convert.(eltype(dynParamSolve), odeProb.u0))
    changeModelParamUse!(odeProbUse.p, odeProbUse.u0, dynParamSolve, paramIndices.namesDynParam)
    solveOdeModelAllCondUse!(simulationInfo.solArray, odeProbUse)

    calcCostNonDyn = (x) -> calcLogLikNotSolveODE(dynamicParamEst, x[iSdUse], x[iObsUse], peTabModel, simulationInfo, paramIndices, measurementData, parameterData)
    @views hessian[vcat(paramIndices.iSdParam, paramIndices.iObsParam), vcat(paramIndices.iSdParam, paramIndices.iObsParam)] .= ForwardDiff.hessian(calcCostNonDyn, noneDynPar)
end


