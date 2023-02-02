function checkGradientResiduals(peTabModel::PeTabModel, solver, tol; verbose::Bool=true)

    # Process PeTab files into type-stable Julia structs 
    experimentalConditionsFile, measurementDataFile, parameterDataFile, observablesDataFile = readPEtabFiles(peTabModel.dirModel, readObservables=true)
    parameterData = processParameters(parameterDataFile)
    measurementData = processMeasurements(measurementDataFile, observablesDataFile) 
    simulationInfo = processSimulationInfo(peTabModel, measurementData)

    # Indices for mapping parameter-estimation vector to dynamic, observable and sd parameters correctly when calculating cost
    paramEstIndices = computeIndicesθ(parameterData, measurementData, peTabModel.odeSystem, experimentalConditionsFile)
        
    # Set model parameter values to those in the PeTab parameter data ensuring correct value of constant parameters 
    setParamToFileValues!(peTabModel.paramMap, peTabModel.stateMap, parameterData)
    priorInfo::PriorInfo = processPriors(paramEstIndices, parameterDataFile)

    # The time-span 5e3 is overwritten when performing actual forward simulations 
    odeProb = ODEProblem(peTabModel.odeSystem, peTabModel.stateMap, (0.0, 5e3), peTabModel.paramMap, jac=true, sparse=false)
    odeProb = remake(odeProb, p = convert.(Float64, odeProb.p), u0 = convert.(Float64, odeProb.u0))
    # Functions to map experimental conditions and parameters correctly to the ODE model 
    changeToExperimentalCondUse! = (pVec, u0Vec, expID, dynParamEst) -> _changeExperimentalCondition!(pVec, u0Vec, expID, dynParamEst, peTabModel, paramEstIndices)
    changeModelParamUse! = (pVec, u0Vec, paramEst) -> changeODEProblemParameters!(pVec, u0Vec, paramEst, paramEstIndices, peTabModel)
    solveOdeModelAllCondUse! = (solArrayArg, odeProbArg, dynParamEst, expIDSolveArg) -> solveODEAllExperimentalConditions!(solArrayArg, odeProbArg, dynParamEst, changeToExperimentalCondUse!, simulationInfo, solver, tol, tol, peTabModel.getTStops, onlySaveAtObservedTimes=true, expIDSolve=expIDSolveArg)
    solveOdeModelAllCondGuassNewtonForwardEq! = (solArrayArg, SMat, odeProbArg, dynParamEst, expIDSolveArg) -> solveODEAllExperimentalConditions!(solArrayArg, SMat, odeProbArg, dynParamEst, changeToExperimentalCondUse!, changeModelParamUse!, simulationInfo, solver, tol, tol, peTabModel.getTStops, onlySaveAtObservedTimes=true, expIDSolve=expIDSolveArg)                                           
    evalResiduals = (paramVecEst) -> computeCost(paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, priorInfo, computeResiduals=true)
    evalJacResiduals = (out, paramVecEst) -> computeGaussNewtonHessianApproximation!(out, paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondGuassNewtonForwardEq!, priorInfo, returnJacobian=true)      

    # Extract parameter vector 
    namesParamEst = paramEstIndices.θ_estNames
    paramVecNominal = [parameterData.nominalValue[findfirst(x -> x == namesParamEst[i], parameterData.parameterId)] for i in eachindex(namesParamEst)]
    paramVec = transformθ(paramVecNominal, namesParamEst, paramEstIndices, reverseTransform=true)

    jacOut = zeros(length(paramVec), length(measurementData.time))
    residualGrad = ForwardDiff.gradient(evalResiduals, paramVec)
    evalJacResiduals(jacOut, paramVec)
    sqDiff = sum((sum(jacOut, dims=2) - residualGrad))

    if !(sqDiff ≤ 1e-5)
        @printf("sqDiff = %.3e\n", sqDiff)
        return false
    end
    if verbose == true
        @printf("sqDiff = %.3e\n", sqDiff)
    end
    return true
end