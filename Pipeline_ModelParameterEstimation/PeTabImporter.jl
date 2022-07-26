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


# TODO : Add support for obsPar (parameters that aren't noise parameter but present in observation function)
# TODO : Add support for scale and offset parameters 


# Set up cost function, gradient function, and approximate hessian function. 
# An exact hessian up to numerical precision is obtained by using ForwardDiff
# on the cost function (as the code can propagate dual numbers)
function setUpCostFunc(modelName, evalObs, evalU0, evalSd, solver, tol::Float64)

    # Set up files where the actual model files are stored 
    modelFile = modelName * ".jl"
    readPath = joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels")
    writePath = joinpath(pwd(), "Pipeline_ModelSolver", "IntermediaryResults")
    fixDirectories(writePath)

    # Extract relevant model file (ode-system)
    allModelFiles = getModelFiles(readPath)    
    usedModelFunctionVector = allModelFunctionVector[[allModelFile .== modelFile for allModelFile in allModelFiles]][1]
    
    # Create an ODE-system and run via ModellingToolkit
    sys, initialSpeciesValues, trueParameterValues = usedModelFunctionVector()
    new_sys = ode_order_lowering(sys)

    # Read data on experimental conditions and parameter values 
    experimentalConditions, measurementData, parameterBounds, observableData = readDataFiles(modelName, readObs=true)
    paramData = processParameterData(parameterBounds) # Model data in convient structure 
    firstExpIds, shiftExpIds, simulateSS, parameterNames, stateNames = getSimulationInfo(measurementData, sys)

    # Read file with observed values and put into struct 
    obsData = processObsData(measurementData, observableData)
    
    # Set up a stateMap and paramMap to map parameter correctly to the ODE model via ModellingToolkit
    stateMap = initialSpeciesValues
    paramMap = trueParameterValues     

    # Ensure correct mapping of constant parameters from the Parameters-PeTab file 
    setParamToParamFileVal!(paramMap, stateMap, paramData)

    # Build actual ODE-problem used for further calculations 
    prob = ODEProblem(new_sys, stateMap, (0.0, 5e3), paramMap, jac=true)
    prob = remake(prob, p = convert.(Float64, prob.p), u0 = convert.(Float64, prob.u0))

    # Set up functions allowing mapping to different experiment condtions and of the parameter vector 
    changeToCondUse! = (pVec, u0Vec, expID) -> changeToCondEst!(pVec, u0Vec, expID, paramData, experimentalConditions, parameterNames, paramMap, evalU0)
    changeModelParamUse! = (pVec, u0Vec, paramEstVec, paramIdChange) -> changeModelParam!(pVec, u0Vec, paramEstVec, paramIdChange, paramData, parameterNames, paramMap, evalU0)

    # Set up function which solves the ODE model for all conditions and stores result 
    solveOdeModelAllCondUse! = (solArrayArg, conditionIdArg, probArg) -> solveOdeModelAllCond!(solArrayArg, conditionIdArg, probArg, changeToCondUse!, simulateSS, measurementData, firstExpIds, shiftExpIds, tol, solver)

    # Set up to bookeeping of indices between observed and noise parameters in the large input vector 
    iDynPar, iSdPar, iObsPar, idParamDyn, idSdParam, idObsParam = getIndicesParam(paramData, obsData)
    idParamDyn, idSdParam, idObsParam = string.(idParamDyn), string.(idSdParam), string.(idObsParam)

    # Lower and upper bounds for parameters to estimate 
    idParam = string.(vcat(idParamDyn, idSdParam, idObsParam))
    lowerBounds = [paramData.lowerBounds[findfirst(x -> x == idParam[i], paramData.parameterID)] for i in eachindex(idParam)] 
    upperBounds = [paramData.upperBounds[findfirst(x -> x == idParam[i], paramData.parameterID)] for i in eachindex(idParam)] 
    
    # Arrays to store ODE-solutions and experimental ID solution corresponds to 
    if simulateSS == true
        nForwardSol = Int(sum([length(shiftExpIds[i]) for i in eachindex(shiftExpIds)]))
    else
        nForwardSol = Int64(length(firstExpIds))
    end
    solArray = Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1}(undef, nForwardSol)
    conditionIdSol = Array{String, 1}(undef, nForwardSol)

    # Set up cost funciton and derivative functions 
    evalF = (paramVecEst) -> calcCost(paramVecEst, prob, solArray, conditionIdSol, changeModelParamUse!, solveOdeModelAllCondUse!, obsData, iDynPar, iSdPar, iObsPar, idParamDyn, paramData, idSdParam, idObsParam, evalObs, evalSd)
    evalGradF = (paramVecEst, grad) -> calcGradCost!(grad, paramVecEst, prob, solArray, conditionIdSol, changeModelParamUse!, solveOdeModelAllCondUse!, obsData, iDynPar, iSdPar, iObsPar, idParamDyn, paramData, idSdParam, idObsParam, evalObs, evalSd)
    evalHessianApproxF = (hessian, paramVecEst) -> calcHessianApprox!(hessian, paramVecEst, prob, solArray, conditionIdSol, changeModelParamUse!, solveOdeModelAllCondUse!, obsData, iDynPar, iSdPar, iObsPar, idParamDyn, paramData, idSdParam, idObsParam, evalObs, evalSd)

    paramVecTmp = [paramData.paramVal[findfirst(x -> x == idParam[i], paramData.parameterID)] for i in eachindex(idParam)]
    
    # Transform upper and lower bounds if the case 
    transformVector!(lowerBounds, idParam, paramData, revTransform=true)
    transformVector!(upperBounds, idParam, paramData, revTransform=true)
    transformVector!(paramVecTmp, idParam, paramData, revTransform=true)

    return evalF, evalGradF, evalHessianApproxF, paramVecTmp, lowerBounds, upperBounds, idParam
end


# Calculate the likelhood for a given model by solving the ODE system for all experimental 
# conditions. Exact hessian can be obtained by using ForwardDiff on this function. 
function calcCost(paramVecEst, 
                  prob::ODEProblem,  
                  solArray::Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1},
                  conditionIdSol::Array{String, 1},
                  changeModelParamUse!::Function,
                  solveOdeModelAllCondUse!::Function,
                  obsData::ObservedData,
                  iDynPar, 
                  iSdPar, 
                  iObsPar, 
                  idParamDyn, 
                  paramData::ParamData, 
                  idParamSd, 
                  idParamObs, 
                  evalObs::Function, 
                  evalSd::Function)

    # Split input into observeble and dynamic parameters 
    dynamicParamEst = paramVecEst[iDynPar]
    obsParEst = paramVecEst[iObsPar]
    sdParamEst = paramVecEst[iSdPar]

    # Correctly transform parameters 
    transformVector!(dynamicParamEst, idParamDyn, paramData)
    transformVector!(sdParamEst, idParamSd, paramData)
    transformVector!(obsParEst, idParamObs, paramData)

    # Function solving ODE and calculating likelihood. 
    probUse = remake(prob, p = convert.(eltype(paramVecEst), prob.p), u0 = convert.(eltype(paramVecEst), prob.u0))
    # Map dynamic parameters to ODE-problem via varmap n modelling-toolkit 
    changeModelParamUse!(probUse.p, probUse.u0, dynamicParamEst, idParamDyn)

    logLik = calcLogLikSolveODE(dynamicParamEst, sdParamEst, obsParEst, probUse, solArray, conditionIdSol, obsData, solveOdeModelAllCondUse!, paramData, idParamSd, evalObs, evalSd)

    return logLik
end


# Calculate likelhood value when solving ODE-system 
function calcCostSolveODE(paramVecEst, 
                          prob::ODEProblem,  
                          solArray::Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1},
                          conditionIdSol::Array{String, 1},
                          changeModelParamUse!::Function,
                          solveOdeModelAllCondUse!::Function,
                          obsData::ObservedData,
                          iDynPar, 
                          iSdPar, 
                          idParamDyn, 
                          paramData::ParamData, 
                          idParamSd, 
                          evalObs::Function, 
                          evalSd::Function)

    # Split input into observeble and dynamic parameters 
    dynamicParamEst = paramVecEst[iDynPar]
    obsPar = []
    sdParamEst = paramVecEst[iSdPar]

    # Correctly transform parameters 
    transformVector!(dynamicParamEst, idParamDyn, paramData)
    transformVector!(sdParamEst, idParamSd, paramData)

    # Function solving ODE and calculating likelihood. 
    probUse = remake(prob, p = convert.(eltype(paramVecEst), prob.p), u0 = convert.(eltype(paramVecEst), prob.u0))
    # Map dynamic parameters to ODE-problem via varmap n modelling-toolkit 
    changeModelParamUse!(probUse.p, probUse.u0, dynamicParamEst, idParamDyn)

    logLik = calcLogLikSolveODE(dynamicParamEst, sdParamEst, obsPar, probUse, solArray, conditionIdSol, obsData, solveOdeModelAllCondUse!, paramData, idParamSd, evalObs, evalSd)

    return logLik
end


# Calculate the log-likelihood given the ODE solution with correct transformation for the model output. 
# Need to add support for scale and offset parameters. 
function calcLogLik(dynamicParamEst, 
                    sdParamEst, 
                    obsPar, 
                    solArray, 
                    conditionIdSol, 
                    obsData::ObservedData, 
                    paramData::ParamData, 
                    idParamSd, 
                    evalObs::Function, 
                    evalSd::Function; 
                    calcObsGrad::Bool=false)

    # Setup yMod and fill each location 
    nStates = length(solArray[1].prob.u0)   
    if !isempty(obsPar) && typeof(obsPar[1]) <: ForwardDiff.Dual
        yMod = Array{eltype(obsPar), 1}(undef, length(obsData.yObs))
    else
        yMod = Array{eltype(dynamicParamEst), 1}(undef, length(obsData.yObs))
    end
    for i in eachindex(yMod)
        whichForwardSol = findfirst(x -> x == obsData.conditionId[i], conditionIdSol)
        t = obsData.tObs[i]

        # Allow the code to correctly propagate dual numbers for gradients
        if calcObsGrad == true
            simVal = [solArray[whichForwardSol](t)[i].value for i in 1:nStates]
            if typeof(simVal[1]) <: ForwardDiff.Dual
                simVal = [solArray[whichForwardSol](t)[i].value for i in 1:nStates]
            end
        else
            simVal = solArray[whichForwardSol](t)
        end
        yMod[i] = evalObs(simVal, t, dynamicParamEst, obsPar, paramData, obsData, obsData.observebleDd[i], obsData.conditionId[i]) 
    end
    # Transform model output if required
    transformObsOrDataArr!(yMod, obsData.transformData)

    # Calculate log-lik - I need to map correctly 
    sdVal = Array{eltype(sdParamEst), 1}(undef, length(obsData.yObs))
    for i in eachindex(sdVal)
        if isNumber(obsData.sdParams[i])
            sdVal[i] = parse(Float64, obsData.sdParams[i])
            continue
        end

        whichForwardSol = findfirst(x -> x == obsData.conditionId[i], conditionIdSol)
        t = obsData.tObs[i]
        simVal = solArray[whichForwardSol](t)
        sdVal[i] = evalSd(simVal, t, sdParamEst, dynamicParamEst, paramData, obsData, obsData.observebleDd[i], obsData.conditionId[i])
    end

    logLik = 0.0
    @inbounds for i in eachindex(sdVal)
        if obsData.transformData[i] == :lin
            logLik += log(sdVal[i]) + 0.5*log(2*pi) + 0.5*((yMod[i] - obsData.yObs[i]) / sdVal[i])^2
        elseif obsData.transformData[i] == :log10
            logLik += log(sdVal[i] * exp10(yMod[i]) * log(10)) + 0.5*log(2*pi) + 0.5*((yMod[i] - obsData.yObs[i]) / sdVal[i])^2
        else
            println("Transformation ", obsData.transformData[i], "not yet supported.")
        end
    end

    return logLik
end


# Calculate the log-likelihood by solving the ODE-problem for all experimental conditions. 
function calcLogLikSolveODE(dynamicParamEst, 
                            sdParamEst, 
                            obsPar, 
                            prob::ODEProblem,
                            solArray::Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1},
                            conditionIdSol::Array{String, 1}, 
                            obsData::ObservedData, 
                            solveOdeModelAllCondUse!::Function, 
                            paramData::ParamData, 
                            idParamSd, 
                            evalObs::Function, 
                            evalSd::Function)

    # In place solving (store ODE-solution in solArray)
    success = solveOdeModelAllCondUse!(solArray, conditionIdSol, prob)
    if success != true
        return Inf
    end

    logLik = calcLogLik(dynamicParamEst, sdParamEst, obsPar, solArray, conditionIdSol, obsData, paramData, idParamSd, evalObs, evalSd)

    return logLik
end


# Calculate likelhood value when solving ODE-system 
function calcCostSolveODE(paramVecEst, 
                          prob::ODEProblem,  
                          solArray::Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1},
                          conditionIdSol::Array{String, 1},
                          changeModelParamUse!::Function,
                          solveOdeModelAllCondUse!::Function,
                          obsData::ObservedData,
                          iDynPar, 
                          iSdPar, 
                          iObsPar,
                          idParamDyn, 
                          paramData::ParamData, 
                          idParamSd, 
                          idParamObs,
                          evalObs::Function, 
                          evalSd::Function)

    # Split input into observeble and dynamic parameters 
    dynamicParamEst = paramVecEst[iDynPar]
    obsParEst = [iObsPar]
    sdParamEst = paramVecEst[iSdPar]

    # Correctly transform parameters 
    transformVector!(dynamicParamEst, idParamDyn, paramData)
    transformVector!(sdParamEst, idParamSd, paramData)
    transformVector!(obsParEst, idParamObs, paramData)

    # Function solving ODE and calculating likelihood. 
    probUse = remake(prob, p = convert.(eltype(paramVecEst), prob.p), u0 = convert.(eltype(paramVecEst), prob.u0))
    # Map dynamic parameters to ODE-problem via varmap n modelling-toolkit 
    changeModelParamUse!(probUse.p, probUse.u0, dynamicParamEst, idParamDyn)

    logLik = calcLogLikSolveODE(dynamicParamEst, sdParamEst, obsParEst, probUse, solArray, conditionIdSol, obsData, solveOdeModelAllCondUse!, paramData, idParamSd, evalObs, evalSd)

    return logLik
end

# Calculate the gradient of the likelihood in two steps. First calculate for the parameters which are 
# a part of the ODE-system. This is expansive as dual numbers have to be propagated through the ODE system. 
# Secondly, calculate the gradient for the observeble parameters given the ODE solution, thus avoiding 
# the step of propegating dual numbers through the ODE solver. 
function calcGradCost!(grad, 
                       paramVecEst, 
                       prob::ODEProblem,  
                       solArray::Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1},
                       conditionIdSol::Array{String, 1},
                       changeModelParamUse!::Function,
                       solveOdeModelAllCondUse!::Function,
                       obsData::ObservedData,
                       iDynPar, 
                       iSdPar, 
                       iObsPar,
                       idParamDyn, 
                       paramData::ParamData, 
                       idParamSd::Array{String, 1}, 
                       idParamObs,
                       evalObs::Function, 
                       evalSd::Function)

    # Split input into observeble and dynamic parameters 
    dynamicParamEst = paramVecEst[iDynPar]
    obsPar = paramVecEst[iObsPar]
    sdParamEst = paramVecEst[iSdPar]

    nonDynPar = vcat(sdParamEst, obsPar)
    iSdUse = 1:length(sdParamEst)
    iObsUse = (length(sdParamEst)+1):(length(sdParamEst) + length(obsPar))

    calcCostDyn = (x) -> calcLogLikSolveOdeGrad(x, sdParamEst, obsPar, prob, solArray, conditionIdSol, obsData, solveOdeModelAllCondUse!, paramData, changeModelParamUse!, idParamDyn, idParamSd, idParamObs, evalObs, evalSd)
    grad[iDynPar] .= ForwardDiff.gradient(calcCostDyn, dynamicParamEst)

    # Correctly transform parameters beofe generating the soluation arrary 
    dynParamSolve = deepcopy(dynamicParamEst)
    transformVector!(dynParamSolve, idParamDyn, paramData)
    probUse = remake(prob, p = convert.(eltype(dynParamSolve), prob.p), u0 = convert.(eltype(dynParamSolve), prob.u0))
    changeModelParamUse!(probUse.p, probUse.u0, dynParamSolve, idParamDyn)
    solveOdeModelAllCondUse!(solArray, conditionIdSol, probUse)

    calcCostNonDyn = (x) -> calcLogLikNotSolveODE(dynamicParamEst, x[iSdUse], x[iObsUse], solArray, conditionIdSol, obsData, paramData, idParamSd, idParamDyn, idParamObs, evalObs, evalSd, calcObsGrad=false)
    grad[vcat(iSdPar, iObsPar)] .= ForwardDiff.gradient(calcCostNonDyn, nonDynPar)
end


# Calculate the log-likelihood by solving the ODE system. This function is used for gradient calculations, 
# and this function is used to obtain the gradient for the dynamic parameters. 
function calcLogLikSolveOdeGrad(dynamicParamEst, 
                                sdParamEst, 
                                obsPar, 
                                prob::ODEProblem,
                                solArray::Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1},
                                conditionIdSol::Array{String, 1}, 
                                obsData::ObservedData, 
                                solveOdeModelAllCondUse!::Function, 
                                paramData::ParamData, 
                                changeModelParamUse!::Function, 
                                idParamDyn::Array{String, 1}, 
                                idParamSd::Array{String, 1}, 
                                idParamObs::Array{String, 1}, 
                                evalObs::Function, 
                                evalSd::Function)

    # Avoid having transformation overwriting other part of code 
    dynamicParamEstUse = dynamicParamEst[:]
    sdParamEstUse = sdParamEst[:]
    obsParamEstUse = obsPar[:]

    # Correctly transform parameters within AD part of code 
    transformVector!(dynamicParamEstUse, idParamDyn, paramData)
    transformVector!(sdParamEstUse, idParamSd, paramData)
    transformVector!(obsParamEstUse, idParamObs, paramData)

    # Function solving ODE and calculating likelihood. 
    probUse = remake(prob, p = convert.(eltype(dynamicParamEstUse), prob.p), u0 = convert.(eltype(dynamicParamEstUse), prob.u0))
    # Map dynamic parameters to ODE-problem via varmap n modelling-toolkit 
    changeModelParamUse!(probUse.p, probUse.u0, dynamicParamEstUse, idParamDyn)
    logLik = calcLogLikSolveODE(dynamicParamEstUse, sdParamEstUse, obsParamEstUse, probUse, solArray, conditionIdSol, obsData, solveOdeModelAllCondUse!, paramData, idParamSd, evalObs, evalSd)

    return logLik
end


# Calculate the log-likelihood given the solution of the ODE system. This function is used for gradient calculations, 
# and this function is used to obtain the gradient for the observable parameters. 
function calcLogLikNotSolveODE(dynamicParamEst, 
                               sdParamEst,
                               obsPar,
                               solArray::Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1},
                               conditionIdSol::Array{String, 1}, 
                               obsData::ObservedData, 
                               paramData::ParamData, 
                               idParamSd::Array{String, 1}, 
                               idParamDyn::Array{String, 1}, 
                               idParamObs::Array{String, 1},
                               evalObs::Function, 
                               evalSd::Function;
                               calcObsGrad=true)

    # Avoid having transformation overwriting other part of code 
    dynamicParamEstUse = dynamicParamEst[:]
    sdParamEstUse = sdParamEst[:]
    obsParamEstUse = obsPar[:]

    # Correctly transform parameters within AD part of code 
    transformVector!(dynamicParamEstUse, idParamDyn, paramData)
    transformVector!(sdParamEstUse, idParamSd, paramData)
    transformVector!(obsParamEstUse, idParamObs, paramData)

    logLik = calcLogLik(dynamicParamEstUse, sdParamEstUse, obsParamEstUse, solArray, conditionIdSol, obsData, paramData, idParamSd, evalObs, evalSd, calcObsGrad=calcObsGrad)

    return logLik
end


# Calculate an approximation of the hessian where the interaction terms between observeble and 
# dynamic parameters is assumed to be zero. 
function calcHessianApprox!(hessian, 
                            paramVecEst, 
                            prob::ODEProblem,  
                            solArray::Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1},
                            conditionIdSol::Array{String, 1},
                            changeModelParamUse!::Function,
                            solveOdeModelAllCondUse!::Function,
                            obsData::ObservedData,
                            iDynPar, 
                            iSdPar, 
                            iObsPar,
                            idParamDyn, 
                            paramData::ParamData, 
                            idParamSd::Array{String, 1}, 
                            idParamObs, 
                            evalObs, 
                            evalSd)

    # Avoid incorrect non-zero values 
    hessian .= 0.0

    dynamicParamEst = paramVecEst[iDynPar]
    obsPar = paramVecEst[iObsPar]
    sdParamEst = paramVecEst[iSdPar]

    nonDynPar = vcat(sdParamEst, obsPar)
    iSdUse = 1:length(sdParamEst)
    iObsUse = (length(sdParamEst)+1):(length(sdParamEst) + length(obsPar))

    calcCostDyn = (x) -> calcLogLikSolveOdeGrad(x, sdParamEst, obsPar, prob, solArray, conditionIdSol, obsData, solveOdeModelAllCondUse!, paramData, changeModelParamUse!, idParamDyn, idParamSd, idParamObs, evalObs, evalSd)
    hessian[iDynPar, iDynPar] .= ForwardDiff.hessian(calcCostDyn, dynamicParamEst)
    
    # To get types correct a reasonable solArray is needed, could in future create specific solArray from cost function evaluation and another for gradient 
    dynParamSolve = deepcopy(dynamicParamEst)
    transformVector!(dynParamSolve, idParamDyn, paramData)
    probUse = remake(prob, p = convert.(eltype(dynParamSolve), prob.p), u0 = convert.(eltype(dynParamSolve), prob.u0))
    changeModelParamUse!(probUse.p, probUse.u0, dynParamSolve, idParamDyn)
    solveOdeModelAllCondUse!(solArray, conditionIdSol, probUse)

    calcCostNonDyn = (x) -> calcLogLikNotSolveODE(dynamicParamEst, x[iSdUse], x[iObsUse], solArray, conditionIdSol, obsData, paramData, idParamSd, idParamDyn, idParamObs, evalObs, evalSd, calcObsGrad=false)
    hessian[vcat(iSdPar, iObsPar), vcat(iSdPar, iObsPar)] .= ForwardDiff.hessian(calcCostNonDyn, nonDynPar)

end


function testOptimizerHessian(modelName, dirSave, solver, nEvals)


    # Set up directory where to save results 
    if !isdir(dirSave)
        mkpath(dirSave)
    end

    # Set up paths for data files 
    modelPath = joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels")
    modelNameShort = modelName[7:end]
    modelFile = modelName * ".jl" 
    methodPath = joinpath(pwd(), "Pipeline_ModelParameterEstimation", modelName)
    # Reading data for model
    readDataPath = joinpath(pwd(), "Pipeline_ModelParameterEstimation", "Data", modelName)
    dataEnding = modelNameShort * ".tsv"
    experimentalConditions = CSV.read(joinpath(readDataPath, "experimentalCondition_" * dataEnding), DataFrame)
    measurementData = CSV.read(joinpath(readDataPath, "measurementData_" * dataEnding), DataFrame)
    observables = CSV.read(joinpath(readDataPath, "observables_" * dataEnding), DataFrame)
    parameterBounds = CSV.read(joinpath(readDataPath, "parameters_" * dataEnding), DataFrame)
    allModelFiles = getModelFiles(modelPath)
    usedModelFunction = allModelFunctionVector[[allModelFile in [modelFile] for allModelFile in allModelFiles]][1]

    # Generate a cost function, inplace gradient function, in place hessian and lower and upper bounds 
    fPre, fGradPre, fHessPre, lowerBounds, upperBounds = forwardAutomaticDifferentiation_hessian_proto_model_Boehm_JProteomeRes2014(usedModelFunction, solver, experimentalConditions, measurementData, observables, parameterBounds)
    nParam = length(lowerBounds)

    # Set up functions into a usable format for hupercube and get start-guesses
    fCube = (x) -> fPre(x...)
    fileSaveCube = dirSave * "Cube.csv"
    createCube(1000, lowerBounds, upperBounds, fileSaveCube, fCube)
    # Read cube 
    cube = Matrix(CSV.read(fileSaveCube, DataFrame))

    # Ipopt and Optim problem 
    ipoptProb, iterArr = createIpoptProb(fPre, fGradPre, fHessPre, lowerBounds, upperBounds)
    Ipopt.AddIpoptIntOption(ipoptProb, "print_level", 0)
    Ipopt.AddIpoptIntOption(ipoptProb, "max_iter", 1000)
    dfOpt, dfcOpt = createOptimProb(fPre, fGradPre, fHessPre, lowerBounds, upperBounds)

    # Evaluate gradient + hessian (allow precompilation)
    dfOpt.f(cube[1, :])
    dfOpt.df(zeros(nParam), cube[1, :])
    dfOpt.h(zeros(nParam, nParam), cube[1, :])

    # Where to save results 
    # Optim 
    costOptim = zeros(Float64, nEvals)
    nIterOptim = zeros(Int, nEvals)
    retCodeOptim = Array{Any, 1}(undef, nEvals)
    runTimeOptim = Array{Any, 1}(undef, nEvals)
    # Ipopt 
    costIpopt = zeros(Float64, nEvals)
    nIterIpopt = zeros(Int, nEvals)
    retCodeIpopt = Array{Any, 1}(undef, nEvals)
    runTimeIpopt = Array{Any, 1}(undef, nEvals)
    # Ipopt no hessian 
    costIpoptBFGS = zeros(Float64, nEvals)
    nIterIpoptBFGS = zeros(Int, nEvals)
    retCodeIpoptBFGS = Array{Any, 1}(undef, nEvals)
    runTimeIpoptBFGS = Array{Any, 1}(undef, nEvals)

    fileSave = dirSave * "OptRes.csv"
    for i in 1:nEvals

        println("I = $i of $nEvals")

        p0 = cube[i, :] # Sample from hypercube 
        # Ipopt with hessian 
        Ipopt.AddIpoptStrOption(ipoptProb, "hessian_approximation", "exact")
        ipoptProb.x = deepcopy(p0)
        benchRunTime = @elapsed  sol_opt = Ipopt.IpoptSolve(ipoptProb)
        # Save relevant values 
        costIpopt[i] = ipoptProb.obj_val
        nIterIpopt[i] = iterArr[1]
        retCodeIpopt[i] = ipoptProb.status
        runTimeIpopt[i] = benchRunTime

        # Ipopt without hessian 
        Ipopt.AddIpoptStrOption(ipoptProb, "hessian_approximation", "limited-memory")
        ipoptProb.x = deepcopy(p0)
        benchRunTime = @elapsed  sol_opt = Ipopt.IpoptSolve(ipoptProb)
        # Save relevant values 
        costIpoptBFGS[i] = ipoptProb.obj_val
        nIterIpoptBFGS[i] = iterArr[1]
        retCodeIpoptBFGS[i] = ipoptProb.status
        runTimeIpoptBFGS[i] = benchRunTime

        # Optim with Hessian 
        res = Optim.optimize(dfOpt, dfcOpt, p0, IPNewton(), Optim.Options(iterations = 1000, show_trace = false))
        costOptim[i] = res.minimum
        nIterOptim[i] = res.iterations
        retCodeOptim[i] = res.f_converged
        runTimeOptim[i] = res.time_run

        # Save after each iteration (do not loose data)
        dataSave = ["IpoptHess" costIpopt[i] runTimeIpopt[i] retCodeIpopt[i] nIterIpopt[i] i;
                "IpoptBGFS" costIpoptBFGS[i] runTimeIpoptBFGS[i] retCodeIpoptBFGS[i] nIterIpoptBFGS[i] i;
                "OptimInt" costOptim[i] runTimeOptim[i] retCodeOptim[i] nIterOptim[i] i]
        dataSave = DataFrame(dataSave, ["Alg", "Cost", "Run_time", "Ret_code", "N_iter", "Start_guess"])
        shouldAppend = isfile(fileSave) ? true : false
        CSV.write(fileSave, dataSave, append=shouldAppend)
    end
end


function genIpoptHess(modelName, solver)

    # Set up paths for data files 
    modelPath = joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels")
    modelNameShort = modelName[7:end]
    modelFile = modelName * ".jl" 
    methodPath = joinpath(pwd(), "Pipeline_ModelParameterEstimation", modelName)
    # Reading data for model
    readDataPath = joinpath(pwd(), "Pipeline_ModelParameterEstimation", "Data", modelName)
    dataEnding = modelNameShort * ".tsv"
    experimentalConditions = CSV.read(joinpath(readDataPath, "experimentalCondition_" * dataEnding), DataFrame)
    measurementData = CSV.read(joinpath(readDataPath, "measurementData_" * dataEnding), DataFrame)
    observables = CSV.read(joinpath(readDataPath, "observables_" * dataEnding), DataFrame)
    parameterBounds = CSV.read(joinpath(readDataPath, "parameters_" * dataEnding), DataFrame)
    allModelFiles = getModelFiles(modelPath)
    usedModelFunction = allModelFunctionVector[[allModelFile in [modelFile] for allModelFile in allModelFiles]][1]

    # Generate a cost function, inplace gradient function, in place hessian and lower and upper bounds 
    fPre, fGradPre, fHessPre, lowerBounds, upperBounds = forwardAutomaticDifferentiation_hessian_proto_model_Boehm_JProteomeRes2014(usedModelFunction, solver, experimentalConditions, measurementData, observables, parameterBounds)
    nParam = length(lowerBounds)

    ipoptProb, iterArr = createIpoptProb(fPre, fGradPre, fHessPre, lowerBounds, upperBounds)
    return ipoptProb, iterArr
end

