function checkGradientResiduals(peTabModel::PeTabModel, solver, tol; verbose::Bool=true)

    # Process PeTab files into type-stable Julia structs 
    experimentalConditionsFile, measurementDataFile, parameterDataFile, observablesDataFile = readDataFiles(peTabModel.dirModel, readObs=true)
    parameterData = processParameterData(parameterDataFile)
    measurementData = processMeasurementData(measurementDataFile, observablesDataFile) 
    simulationInfo = getSimulationInfo(peTabModel, measurementDataFile, measurementData)

    # Indices for mapping parameter-estimation vector to dynamic, observable and sd parameters correctly when calculating cost
    paramEstIndices = getIndicesParam(parameterData, measurementData, peTabModel.odeSystem, experimentalConditionsFile)
        
    # Set model parameter values to those in the PeTab parameter data ensuring correct value of constant parameters 
    setParamToFileValues!(peTabModel.paramMap, peTabModel.stateMap, parameterData)
    priorInfo::PriorInfo = getPriorInfo(paramEstIndices, parameterDataFile)

    # The time-span 5e3 is overwritten when performing actual forward simulations 
    odeProb = ODEProblem(peTabModel.odeSystem, peTabModel.stateMap, (0.0, 5e3), peTabModel.paramMap, jac=true, sparse=false)
    odeProb = remake(odeProb, p = convert.(Float64, odeProb.p), u0 = convert.(Float64, odeProb.u0))
    # Functions to map experimental conditions and parameters correctly to the ODE model 
    changeToExperimentalCondUse! = (pVec, u0Vec, expID, dynParamEst) -> changeExperimentalCondEst!(pVec, u0Vec, expID, dynParamEst, peTabModel, paramEstIndices)
    changeToExperimentalCondSenseEqUse! = (pVec, u0Vec, expID, dynParamEst) -> changeExperimentalCondEstSenseEq!(pVec, u0Vec, expID, dynParamEst, peTabModel, paramEstIndices)
    changeModelParamUse! = (pVec, u0Vec, paramEst) -> changeODEProblemParameters!(pVec, u0Vec, paramEst, paramEstIndices, peTabModel)
    solveOdeModelAllCondUse! = (solArrayArg, odeProbArg, dynParamEst, expIDSolveArg) -> solveOdeModelAllExperimentalCond!(solArrayArg, odeProbArg, dynParamEst, changeToExperimentalCondUse!, simulationInfo, solver, tol, tol, peTabModel.getTStops, onlySaveAtTobs=true, expIDSolve=expIDSolveArg)
    solveOdeModelAllCondGuassNewtonForwardEq! = (solArrayArg, SMat, odeProbArg, dynParamEst, expIDSolveArg) -> solveOdeModelAllExperimentalCond!(solArrayArg, SMat, odeProbArg, dynParamEst, changeToExperimentalCondUse!, changeModelParamUse!, simulationInfo, solver, tol, tol, peTabModel.getTStops, onlySaveAtTobs=true, expIDSolve=expIDSolveArg)                                           
    evalResiduals = (paramVecEst) -> computeCost(paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, priorInfo, computeResiduals=true)
    evalJacResiduals = (out, paramVecEst) -> computeGaussNewtonHessianApproximation!(out, paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondGuassNewtonForwardEq!, priorInfo, returnJacobian=true)      

    # Extract parameter vector 
    namesParamEst = paramEstIndices.namesParamEst
    paramVecNominal = [parameterData.paramVal[findfirst(x -> x == namesParamEst[i], parameterData.parameterID)] for i in eachindex(namesParamEst)]
    paramVec = transformθ(paramVecNominal, namesParamEst, parameterData, reverseTransform=true)

    jacOut = zeros(length(paramVec), length(measurementData.tObs))
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