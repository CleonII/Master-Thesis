include(joinpath(pwd(), "src", "PeTab_importer", "Common.jl"))
include(joinpath(pwd(), "src", "PeTab_importer", "Map_parameters.jl"))
include(joinpath(pwd(), "src", "PeTab_importer", "Create_obs_u0_sd_functions.jl"))
include(joinpath(pwd(), "src", "PeTab_importer", "Create_obs_u0_sd_derivatives.jl"))
include(joinpath(pwd(), "src", "PeTab_importer", "Process_PeTab_files.jl"))
include(joinpath(pwd(), "src", "PeTab_importer", "Distributed.jl"))
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
                           solver::SciMLAlgorithm, 
                           tol::Float64; 
                           sensealg=ForwardDiffSensitivity(),
                           sparseJac::Bool=false, 
                           absTolSS::Float64=1e-8, 
                           relTolSS::Float64=1e-6, 
                           nProcs::Signed=1,
                           adjSolver=Rodas5P(), 
                           adjTol::Float64=1e-6, 
                           adjSensealg::SciMLSensitivity.AbstractAdjointSensitivityAlgorithm=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(false)), 
                           adjSensealgSS::SciMLSensitivity.AbstractAdjointSensitivityAlgorithm=SteadyStateAdjoint())::PeTabOpt

    if !(typeof(adjSensealgSS) <: SteadyStateAdjoint)
        println("If you are using adjoint sensitivity analysis for a model with PreEq-criteria the most 
                 the most efficient adjSensealgSS is usually SteadyStateAdjoint. The algorithm you have 
                 provided, ", adjSensealgSS, "might not work (as there are some bugs here). In case it does 
                 not work, and SteadyStateAdjoint fails (because a dependancy on time) a good choice might 
                 be QuadratureAdjoint(autodiff=false, autojacvec=false)")
    end

    # Process PeTab files into type-stable Julia structs 
    experimentalConditionsFile, measurementDataFile, parameterDataFile, observablesDataFile = readDataFiles(peTabModel.dirModel, readObs=true)
    parameterData = processParameterData(parameterDataFile)
    measurementData = processMeasurementData(measurementDataFile, observablesDataFile) 
    simulationInfo = getSimulationInfo(measurementDataFile, measurementData, absTolSS=absTolSS, relTolSS=relTolSS)

    # Indices for mapping parameter-estimation vector to dynamic, observable and sd parameters correctly when calculating cost
    paramEstIndices = getIndicesParam(parameterData, measurementData, peTabModel.odeSystem, experimentalConditionsFile)
    
    # Set up potential prior for the parameters to estimate 
    priorInfo::PriorInfo = getPriorInfo(paramEstIndices, parameterDataFile)

    # Set model parameter values to those in the PeTab parameter data ensuring correct value of constant parameters 
    setParamToFileValues!(peTabModel.paramMap, peTabModel.stateMap, parameterData)

    # The time-span 5e3 is overwritten when performing actual forward simulations 
    odeProb = ODEProblem(peTabModel.odeSystem, peTabModel.stateMap, (0.0, 5e3), peTabModel.paramMap, jac=true, sparse=sparseJac)
    odeProb = remake(odeProb, p = convert.(Float64, odeProb.p), u0 = convert.(Float64, odeProb.u0))

    # Functions to map experimental conditions and parameters correctly to the ODE model 
    changeToExperimentalCondUse! = (pVec, u0Vec, expID, dynParamEst) -> changeExperimentalCondEst!(pVec, u0Vec, expID, dynParamEst, peTabModel, paramEstIndices)
    changeToExperimentalCondUse = (pVec, u0Vec, expID, dynParamEst) -> changeExperimentalCondEst(pVec, u0Vec, expID, dynParamEst, peTabModel, paramEstIndices)
    changeModelParamUse! = (pVec, u0Vec, paramEst) -> changeModelParam!(pVec, u0Vec, paramEst, paramEstIndices, peTabModel)
    changeModelParamUse = (pVec, paramEst) -> changeModelParam(pVec, paramEst, paramEstIndices, peTabModel)

    # Set up function which solves the ODE model for all conditions and stores result 
    solveOdeModelAllCondUse! = (solArrayArg, odeProbArg, dynParamEst, expIDSolveArg) -> solveOdeModelAllExperimentalCond!(solArrayArg, odeProbArg, dynParamEst, changeToExperimentalCondUse!, simulationInfo, solver, tol, tol, onlySaveAtTobs=true, expIDSolve=expIDSolveArg)
    solveOdeModelAllCondAdjUse! = (solArrayArg, odeProbArg, dynParamEst, expIDSolveArg) -> solveOdeModelAllExperimentalCond!(solArrayArg, odeProbArg, dynParamEst, changeToExperimentalCondUse!, simulationInfo, solver, tol, tol, denseSol=true, expIDSolve=expIDSolveArg)
    solveOdeModelAtCondZygoteUse = (odeProbArg, conditionId, dynParamEst, t_max) -> solveOdeModelAtExperimentalCondZygote(odeProbArg, conditionId, dynParamEst, t_max, changeToExperimentalCondUse, measurementData, simulationInfo, solver, tol, tol, sensealg)

    if nProcs > 1 && nprocs() != nProcs
        println("Error : PEtab importer was set to build the cost, grad and hessian with $nProcs processes, 
                 however, Julia is currently running with ", nprocs(), " processes which does not match input 
                 value. Input argument nProcs must match nprocs()")
    elseif nProcs == 1
        evalF = (paramVecEst) -> calcCost(paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, priorInfo)
        evalGradF = (grad, paramVecEst) -> calcGradCost!(grad, paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, priorInfo)    
        evalHessApprox = (hessianMat, paramVecEst) -> calcHessianApprox!(hessianMat, paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, priorInfo)
        # Sometimes even though the cost can be computed the ODE solver with Dual number can fail (as the error is 
        # monitored slightly differently). When this happens we error out and the code crash, hence the need of a catch statement 
        _evalHess = (paramVecEst) -> calcCost(paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, priorInfo, calcHessian=true)
        evalHess = (hessianMat, paramVec) -> begin 
                                                if all([simulationInfo.solArray[i].retcode == :Success for i in eachindex(simulationInfo.solArray)])
                                                    try 
                                                        hessianMat .= Symmetric(ForwardDiff.hessian(_evalHess, paramVec))
                                                    catch
                                                        hessianMat .= 0.0
                                                    end
                                                else
                                                    hessianMat .= 0.0
                                                end
                                            end
    elseif nProcs > 1 && nprocs() == nProcs
        evalF, evalGradF, evalHess, evalHessApprox = setUpPEtabOptDistributed(peTabModel, solver, tol, parameterData, measurementData, 
                                                                              simulationInfo, paramEstIndices, priorInfo, odeProb)
    end

    evalFZygote = (paramVecEst) -> calcCostZygote(paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse, solveOdeModelAtCondZygoteUse, priorInfo)
    evalGradFZygote = (grad, paramVecEst) -> calcGradZygote!(grad, paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse, solveOdeModelAtCondZygoteUse, priorInfo)
    
    # Gradient computed via adjoint sensitivity analysis (not Zygote interface)
    evalGradFAdjoint = (grad, paramVecEst) -> calcGradCostAdj!(grad, paramVecEst, adjSolver, adjSensealg, adjSensealgSS, adjTol, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondAdjUse!, priorInfo) 
    
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
                        evalFZygote,
                        evalGradF, 
                        evalGradFZygote,
                        evalGradFAdjoint,
                        evalHess,
                        evalHessApprox, 
                        Int64(length(namesParamEst)), 
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
function calcCost(paramVecEst::AbstractVector,
                  odeProb::ODEProblem,  
                  peTabModel::PeTabModel,
                  simulationInfo::SimulationInfo,
                  paramIndices::ParameterIndices,
                  measurementData::MeasurementData,
                  parameterData::ParamData,
                  changeModelParamUse!::Function,
                  solveOdeModelAllCondUse!::Function, 
                  priorInfo::PriorInfo;
                  expIDSolve::Array{String, 1} = ["all"],
                  calcHessian::Bool=false)::Real

    # Correctly map paramVecEst to dynmaic, observable and sd param. The new vectors 
    # are all distinct copies.
    dynamicParamEst = paramVecEst[paramIndices.iDynParam]
    obsParEst = paramVecEst[paramIndices.iObsParam]
    sdParamEst = paramVecEst[paramIndices.iSdParam]
    nonDynamicParamEst = paramVecEst[paramIndices.iNonDynParam]

    logLik = calcLogLikSolveODE(dynamicParamEst, sdParamEst, obsParEst, nonDynamicParamEst, odeProb, peTabModel, simulationInfo, paramIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, priorInfo, calcHessDynParam=calcHessian, expIDSolve=expIDSolve)

    if priorInfo.hasPriors == true
        paramVecEstTransformed = transformParamVec(paramVecEst, paramIndices.namesParamEst, parameterData)
        logLik += evalPriors(paramVecEstTransformed, paramVecEst, paramIndices.namesParamEst, paramIndices, priorInfo)
    end

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
function calcGradCost!(grad::Vector{Float64}, 
                       paramVecEst::Vector{Float64}, 
                       odeProb::ODEProblem,  
                       peTabModel::PeTabModel,
                       simulationInfo::SimulationInfo,
                       paramIndices::ParameterIndices,
                       measurementData::MeasurementData,
                       parameterData::ParamData, 
                       changeModelParamUse!::Function, 
                       solveOdeModelAllCondUse!::Function, 
                       priorInfo::PriorInfo;
                       expIDSolve::Array{String, 1} = ["all"])     

    # Split input into observeble and dynamic parameters 
    dynamicParamEst = paramVecEst[paramIndices.iDynParam]
    obsParEst = paramVecEst[paramIndices.iObsParam]
    sdParamEst = paramVecEst[paramIndices.iSdParam]
    noneDynParamEst = paramVecEst[paramIndices.iNonDynParam]
    namesSdParam = paramIndices.namesSdParam
    namesObsParam = paramIndices.namesObsParam
    namesNonDynParam = paramIndices.namesNonDynParam
    namesSdObsNonDynPar = paramIndices.namesSdObsNonDynPar

    # Calculate gradient seperately for dynamic and non dynamic parameter. 

    # I have tried to decrease run time here with chunking without success (deafult value performs best). Might be 
    # worth to look into a parellisation over the chunks (as for larger models each call takes relatively long time). 
    # Also parellisation of the chunks should be faster than paralellisation over experimental condtions.
    calcCostDyn = (x) -> calcLogLikSolveODE(x, sdParamEst, obsParEst, noneDynParamEst, odeProb, peTabModel, simulationInfo, paramIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, priorInfo, calcGradDynParam=true, expIDSolve=expIDSolve)
    try 
        grad[paramIndices.iDynParam] .= ForwardDiff.gradient(calcCostDyn, dynamicParamEst)::Vector{Float64}
    catch
        grad = 1e8
        return
    end

    @inbounds for i in eachindex(simulationInfo.solArrayGrad)
        if expIDSolve[1] == "all" || simulationInfo.conditionIdSol[i] ∈ expIDSolve
            if simulationInfo.solArrayGrad[i].retcode != :Success
                grad .= 1e8
                return 
            end
        end
    end

    # Happens when at least one forward pass fails 
    if all(grad[paramIndices.iDynParam] .== 0.0)
        grad .= 1e8
        return 
    end

    # Here it is crucial to account for that obs- and sd parameter can be overlapping. Thus, a name-map 
    # of both Sd and Obs param is used to account for this. This is not a worry for non-dynamic parameters.
    paramNotOdeSys = paramVecEst[paramIndices.iSdObsNonDynPar]
    iSdUse = [findfirst(x -> x == namesSdParam[i], namesSdObsNonDynPar) for i in eachindex(namesSdParam)]
    iObsUse = [findfirst(x -> x == namesObsParam[i],  namesSdObsNonDynPar) for i in eachindex(namesObsParam)]
    iNonDynUse = [findfirst(x -> x == namesNonDynParam[i],  namesSdObsNonDynPar) for i in eachindex(namesNonDynParam)]

    # TODO : Make choice of gradient availble 
    calcCostNonDyn = (x) -> calcLogLikNotSolveODE(dynamicParamEst, x[iSdUse], x[iObsUse], x[iNonDynUse], peTabModel, simulationInfo, paramIndices, measurementData, parameterData, priorInfo, calcGradObsSdParam=true, expIDSolve=expIDSolve)
    @views ReverseDiff.gradient!(grad[paramIndices.iSdObsNonDynPar], calcCostNonDyn, paramNotOdeSys)
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
function calcHessianApprox!(hessian::Matrix{Float64}, 
                            paramVecEst::Vector{Float64},
                            odeProb::ODEProblem,  
                            peTabModel::PeTabModel,
                            simulationInfo::SimulationInfo,
                            paramIndices::ParameterIndices,
                            measurementData::MeasurementData,
                            parameterData::ParamData, 
                            changeModelParamUse!::Function,
                            solveOdeModelAllCondUse!::Function, 
                            priorInfo::PriorInfo; 
                            expIDSolve::Array{String, 1} = ["all"]) 

    # Avoid incorrect non-zero values 
    hessian .= 0.0

    # Split input into observeble and dynamic parameters 
    dynamicParamEst = paramVecEst[paramIndices.iDynParam]
    obsParEst = paramVecEst[paramIndices.iObsParam]
    sdParamEst = paramVecEst[paramIndices.iSdParam]
    noneDynParamEst = paramVecEst[paramIndices.iNonDynParam]
    namesSdParam = paramIndices.namesSdParam
    namesObsParam = paramIndices.namesObsParam
    namesNonDynParam = paramIndices.namesNonDynParam
    namesSdObsNonDynPar = paramIndices.namesSdObsNonDynPar

    # Calculate gradient seperately for dynamic and non dynamic parameter. 
    calcCostDyn = (x) -> calcLogLikSolveODE(x, sdParamEst, obsParEst, noneDynParamEst, odeProb, peTabModel, simulationInfo, paramIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, priorInfo, calcHessDynParam=true, expIDSolve=expIDSolve)
    try 
        hessian[paramIndices.iDynParam, paramIndices.iDynParam] .= ForwardDiff.hessian(calcCostDyn, dynamicParamEst)::Matrix{Float64}
    catch
        hessian .= 0.0
        return 
    end

    @inbounds for i in eachindex(simulationInfo.solArrayGrad)
        if expIDSolve[1] == "all" || simulationInfo.conditionIdSol[i] ∈ expIDSolve
            if simulationInfo.solArrayGrad[i].retcode != :Success
                hessian .= 0.0
                return 
            end
        end
    end

    # Happens when at least one forward pass fails 
    if all(hessian[paramIndices.iDynParam, paramIndices.iDynParam] .== 0.0)
        return 
    end

    # Here it is crucial to account for that obs- and sd parameter can be overlapping. Thus, a name-map 
    # of both Sd and Obs param is used to account for this 
    paramNotOdeSys = paramVecEst[paramIndices.iSdObsNonDynPar]
    iSdUse = [findfirst(x -> x == namesSdParam[i], namesSdObsNonDynPar) for i in eachindex(namesSdParam)]
    iObsUse = [findfirst(x -> x == namesObsParam[i],  namesSdObsNonDynPar) for i in eachindex(namesObsParam)]
    iNonDynUse = [findfirst(x -> x == namesNonDynParam[i],  namesSdObsNonDynPar) for i in eachindex(namesNonDynParam)]
    
    # Compute hessian for none dynamic parameters 
    calcCostNonDyn = (x) -> calcLogLikNotSolveODE(dynamicParamEst, x[iSdUse], x[iObsUse], x[iNonDynUse], peTabModel, simulationInfo, paramIndices, measurementData, parameterData, priorInfo, expIDSolve=expIDSolve)
    @views ForwardDiff.hessian!(hessian[paramIndices.iSdObsNonDynPar, paramIndices.iSdObsNonDynPar], calcCostNonDyn, paramNotOdeSys)

end


"""
    calcLogLikSolveODE(dynamicParamEst, 
                       sdParamEst, 
                       obsParEst, 
                       nonDynParamEst,
                       odeProb::ODEProblem,
                       peTabModel::PeTabModel,
                       simulationInfo::SimulationInfo,
                       measurementData ::MeasurementData, 
                       parameterData::ParamData, 
                       changeModelParamUse!::Function,
                       solveOdeModelAllCondUse!::Function, 
                       computeGradOrHess::Bool=false, 
                       )

    Helper function computing the likelhood by solving the ODE system for all 
    PeTab-specifed experimental conditions using the dynamic-parameters, 
    sd-parameters and observable parameters. 

    When computing the cost and gradient/hessian for dynamic parameters the ODE 
    system must be solved before getting the likelhood. Besides the different 
    parameter vector the input arguements are the same as for `calcCost`, and 
    everything is setup by `setUpCostGradHess` function.

    See also: [`calcCost`, `setUpCostGradHess`]
"""
function calcLogLikSolveODE(dynamicParamEst::AbstractVector,
                            sdParamEst::AbstractVector,
                            obsParEst::AbstractVector,
                            nonDynParamEst::AbstractVector,
                            odeProb::ODEProblem,
                            peTabModel::PeTabModel,
                            simulationInfo::SimulationInfo,
                            paramIndices::ParameterIndices,
                            measurementData ::MeasurementData, 
                            parameterData::ParamData, 
                            changeModelParamUse!::Function,
                            solveOdeModelAllCondUse!::Function, 
                            priorInfo::PriorInfo;
                            calcHessDynParam::Bool=false, 
                            calcGradDynParam::Bool=false, 
                            expIDSolve::Array{String, 1} = ["all"])::Real

    dynamicParamEstUse = transformParamVec(dynamicParamEst, paramIndices.namesDynParam, parameterData)
    sdParamEstUse = transformParamVec(sdParamEst, paramIndices.namesSdParam, parameterData)
    obsParEstUse = transformParamVec(obsParEst, paramIndices.namesObsParam, parameterData)
    nonDynParamEstUse = transformParamVec(nonDynParamEst, paramIndices.namesNonDynParam, parameterData)

    odeProbUse = remake(odeProb, p = convert.(eltype(dynamicParamEstUse), odeProb.p), u0 = convert.(eltype(dynamicParamEstUse), odeProb.u0))
    changeModelParamUse!(odeProbUse.p, odeProbUse.u0, dynamicParamEstUse)
    
    # If computing hessian or gradient store ODE solution in arrary with dual numbers, else use 
    # solution array with floats
    if calcHessDynParam == true || calcGradDynParam == true
        success = solveOdeModelAllCondUse!(simulationInfo.solArrayGrad, odeProbUse, dynamicParamEstUse, expIDSolve)
    else
        success = solveOdeModelAllCondUse!(simulationInfo.solArray, odeProbUse, dynamicParamEstUse, expIDSolve)
    end
    if success != true
        println("Failed to solve ODE model")
        return Inf
    end

    logLik = calcLogLik(dynamicParamEstUse, sdParamEstUse, obsParEstUse, nonDynParamEstUse, peTabModel, simulationInfo, paramIndices, measurementData, parameterData, expIDSolve, calcHessDynParam=calcHessDynParam, calcGradDynParam=calcGradDynParam)

    if priorInfo.hasPriors == true && (calcHessDynParam || calcGradDynParam)        
        logLik += evalPriors(dynamicParamEstUse, dynamicParamEst, paramIndices.namesDynParam, paramIndices, priorInfo)
    end

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
function calcLogLikNotSolveODE(dynamicParamEst::Vector{Float64}, 
                               sdParamEst::AbstractVector,
                               obsParamEst::AbstractVector,
                               nonDynParamEst::AbstractVector,
                               peTabModel::PeTabModel,
                               simulationInfo::SimulationInfo,
                               paramIndices::ParameterIndices,
                               measurementData::MeasurementData,
                               parameterData::ParamData, 
                               priorInfo::PriorInfo;
                               calcGradObsSdParamAdj::Bool=false,
                               calcGradObsSdParam::Bool=false, 
                               expIDSolve::Array{String, 1} = ["all"])::Real 

    # To be able to use ReverseDiff sdParamEstUse and obsParamEstUse cannot be overwritten. 
    # Hence new vectors have to be created. Minimal overhead.
    dynamicParamEstUse = transformParamVec(dynamicParamEst, paramIndices.namesDynParam, parameterData)
    sdParamEstUse = transformParamVec(sdParamEst, paramIndices.namesSdParam, parameterData)
    obsParamEstUse = transformParamVec(obsParamEst, paramIndices.namesObsParam, parameterData)
    nonDynParamEstUse = transformParamVec(nonDynParamEst, paramIndices.namesNonDynParam, parameterData)

    logLik = calcLogLik(dynamicParamEstUse, sdParamEstUse, obsParamEstUse, nonDynParamEstUse, peTabModel, simulationInfo, paramIndices, 
                        measurementData, parameterData, expIDSolve, calcGradObsSdParam=calcGradObsSdParam, 
                        calcGradObsSdParamAdj=calcGradObsSdParamAdj)
    
    if priorInfo.hasPriors == true 
        logLik += evalPriors(sdParamEstUse, sdParamEst, paramIndices.namesSdParam, paramIndices, priorInfo)
        logLik += evalPriors(obsParamEstUse, obsParamEst, paramIndices.namesObsParam, paramIndices, priorInfo)
        logLik += evalPriors(nonDynParamEstUse, nonDynParamEst, paramIndices.namesNonDynParam, paramIndices, priorInfo)
    end

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
function calcLogLik(dynamicParamEst::AbstractVector,
                    sdParamEst::AbstractVector, 
                    obsPar::AbstractVector, 
                    nonDynParamEst::AbstractVector,
                    peTabModel::PeTabModel,
                    simulationInfo::SimulationInfo,
                    paramIndices::ParameterIndices,
                    measurementData::MeasurementData, 
                    parameterData::ParamData,
                    expIDSolve::Array{String, 1} = ["all"];
                    calcHessDynParam::Bool=false, 
                    calcGradDynParam::Bool=false, 
                    calcGradObsSdParamAdj::Bool=false,
                    calcGradObsSdParam::Bool=false)::Real 

    if calcHessDynParam == true || calcGradDynParam == true || calcGradObsSdParam == true || calcGradObsSdParamAdj == true
        odeSolArray = simulationInfo.solArrayGrad
    else
        odeSolArray = simulationInfo.solArray
    end

    # Compute yMod and sd-val by looping through all experimental conditons. At the end 
    # update the likelihood 
    logLik = 0.0
    for conditionID in keys(measurementData.iPerConditionId)

        if expIDSolve[1] != "all" && conditionID ∉ expIDSolve
            continue
        end

        # Extract the ODE-solution for specific condition ID
        whichForwardSol = findfirst(x -> x == conditionID, simulationInfo.conditionIdSol)
        odeSol = odeSolArray[whichForwardSol]   
        logLik += calcLogLikExpCond(odeSol, dynamicParamEst, sdParamEst, obsPar, 
                                    nonDynParamEst, peTabModel, conditionID, paramIndices,
                                    measurementData, parameterData, calcGradObsSdParam, 
                                    calcGradObsSdParamAdj=calcGradObsSdParamAdj)

        if isinf(logLik)
            return Inf
        end
    end

    return logLik
end


function calcCostZygote(paramVecEst,
                        odeProb::ODEProblem,  
                        peTabModel::PeTabModel,
                        simulationInfo::SimulationInfo,
                        paramIndices::ParameterIndices,
                        measurementData::MeasurementData,
                        parameterData::ParamData,
                        changeModelParamUse::Function,
                        solveOdeModelAllCondZygoteUse::Function, 
                        priorInfo::PriorInfo)
    
    # Correctly map paramVecEst to dynmaic, observable and sd param. The vectors, 
    # e.g dynamicParamEst, are distinct copies so transforming them will not change 
    # paramVecEst.
    dynamicParamEst = paramVecEst[paramIndices.iDynParam]
    obsParEst = paramVecEst[paramIndices.iObsParam]
    sdParamEst = paramVecEst[paramIndices.iSdParam]
    nonDynamicParamEst = paramVecEst[paramIndices.iNonDynParam]

    logLik = calcLogLikZygote(dynamicParamEst,
                              sdParamEst,
                              obsParEst,
                              nonDynamicParamEst,
                              odeProb,
                              peTabModel,
                              simulationInfo,
                              paramIndices,
                              measurementData,
                              parameterData,
                              changeModelParamUse,
                              solveOdeModelAllCondZygoteUse, 
                              priorInfo)

    if priorInfo.hasPriors == true
        paramVecEstTransformed = transformParamVec(paramVecEst, paramIndices.namesParamEst, parameterData)
        logLik += evalPriors(paramVecEstTransformed, paramVecEst, paramIndices.namesParamEst, paramIndices, priorInfo)
    end                                  

    return logLik                          
end


function calcGradZygote!(grad::T1, 
                         paramVecEst::T2, 
                         odeProb::ODEProblem,  
                         peTabModel::PeTabModel,
                         simulationInfo::SimulationInfo,
                         paramIndices::ParameterIndices,
                         measurementData::MeasurementData,
                         parameterData::ParamData, 
                         changeModelParamUse::Function, 
                         solveOdeModelAllCondZygoteUse::Function, 
                         priorInfo::PriorInfo) where {T1<:Array{<:AbstractFloat, 1}, 
                                                      T2<:Vector{<:Real}}
    
    # Split input into observeble and dynamic parameters 
    dynamicParamEst = paramVecEst[paramIndices.iDynParam]
    obsParEst = paramVecEst[paramIndices.iObsParam]
    sdParamEst = paramVecEst[paramIndices.iSdParam]
    noneDynParamEst = paramVecEst[paramIndices.iNonDynParam]
    namesSdParam = paramIndices.namesSdParam
    namesObsParam = paramIndices.namesObsParam
    namesNonDynParam = paramIndices.namesNonDynParam
    namesSdObsNonDynPar = paramIndices.namesSdObsNonDynPar

    # Calculate gradient seperately for dynamic and non dynamic parameter. This seems to considerble help Zygote 
    # if the model is large enough (which I guess has to with how the derivatives are overloaded when calling sol)

    # I have tried to decrease run time here with chunking without success (deafult value performs best). Might be 
    # worth to look into a parellisation over the chunks (as for larger models each call takes relatively long time). 
    # Also parellisation of the chunks should be faster than paralellisation over experimental condtions.
    calcCostDyn = (x) -> calcLogLikZygote(x, sdParamEst, obsParEst, noneDynParamEst, odeProb, peTabModel, simulationInfo, paramIndices, measurementData, parameterData, changeModelParamUse, solveOdeModelAllCondZygoteUse, priorInfo, evalGradDyn=true)
    grad[paramIndices.iDynParam] .= Zygote.gradient(calcCostDyn, dynamicParamEst)[1]

    # Here it is crucial to account for that obs- and sd parameter can be overlapping. Thus, a name-map 
    # of both Sd and Obs param is used to account for this. This is not a worry for non-dynamic parameters.
    paramNotOdeSys = paramVecEst[paramIndices.iSdObsNonDynPar]
    iSdUse = [findfirst(x -> x == namesSdParam[i], namesSdObsNonDynPar) for i in eachindex(namesSdParam)]
    iObsUse = [findfirst(x -> x == namesObsParam[i],  namesSdObsNonDynPar) for i in eachindex(namesObsParam)]
    iNonDynUse = [findfirst(x -> x == namesNonDynParam[i],  namesSdObsNonDynPar) for i in eachindex(namesNonDynParam)]

    # This is subtle. By using Zygote-ignore when solving the ODE the ODE solution is stored in simulationInfo.solArray
    # which can be used to compute the gradient for the non-dynamic parameters without having to resolve the ODE system.
    calcCostNonDyn = (x) -> calcLogLikNotSolveODE(dynamicParamEst, x[iSdUse], x[iObsUse], x[iNonDynUse], peTabModel, simulationInfo, paramIndices, measurementData, parameterData, priorInfo, calcGradObsSdParam=false)
    @views ForwardDiff.gradient!(grad[paramIndices.iSdObsNonDynPar], calcCostNonDyn, paramNotOdeSys)
end


# Computes the likelihood in such a in a Zygote compatible way, which mainly means that no arrays are mutated.
function calcLogLikZygote(dynamicParamEst,
                          sdParamEst,
                          obsParEst,
                          nonDynamicParamEst,
                          odeProb::ODEProblem,  
                          peTabModel::PeTabModel,
                          simulationInfo::SimulationInfo,
                          paramIndices::ParameterIndices,
                          measurementData::MeasurementData,
                          parameterData::ParamData,
                          changeModelParamUse::Function,
                          solveOdeModelAllCondZygoteUse::Function, 
                          priorInfo::PriorInfo; 
                          evalGradDyn::Bool=false)::Real

    # Correctly transform parameter if, for example, they are on the log-scale.
    dynamicParamEstUse = transformParamVec(dynamicParamEst, paramIndices.namesDynParam, parameterData)
    sdParamEstUse = transformParamVec(sdParamEst, paramIndices.namesSdParam, parameterData)
    obsParEstUse = transformParamVec(obsParEst, paramIndices.namesObsParam, parameterData)
    nonDynamicParamEstUse = transformParamVec(nonDynamicParamEst, paramIndices.namesNonDynParam, parameterData)                          
                                                            
    pOdeSysUse, u0Use = changeModelParamUse(odeProb.p, dynamicParamEstUse)
    odeProbUse = remake(odeProb, p = convert.(eltype(dynamicParamEstUse), pOdeSysUse), u0 = convert.(eltype(dynamicParamEstUse), u0Use))
    
    # Compute yMod and sd-val by looping through all experimental conditons. At the end 
    # update the likelihood 
    logLik = convert(eltype(dynamicParamEstUse), 0.0)
    for conditionID in keys(measurementData.iPerConditionId)
        
        # Solve ODE system 
        whichTMax = findfirst(x -> x == conditionID, simulationInfo.conditionIdSol)
        tMax = simulationInfo.tMaxForwardSim[whichTMax]
        
        odeSol, success = solveOdeModelAllCondZygoteUse(odeProbUse, conditionID, dynamicParamEstUse, tMax)
        if success != true
            return Inf
        end

        logLik += calcLogLikExpCond(odeSol, dynamicParamEstUse, sdParamEstUse, obsParEstUse, 
                                    nonDynamicParamEstUse, peTabModel, conditionID, paramIndices,
                                    measurementData, parameterData, false)

        if isinf(logLik)
            return logLik
        end
    end

    if priorInfo.hasPriors == true && evalGradDyn == true
        logLik += evalPriors(dynamicParamEstUse, dynamicParamEst, paramIndices.namesDynParam, paramIndices, priorInfo)
    end

    return logLik
end


function calcGradCostAdj!(grad::Vector{Float64}, 
                          paramVecEst::Vector{Float64}, 
                          adjSolver::SciMLAlgorithm, 
                          sensealg::SciMLSensitivity.AbstractAdjointSensitivityAlgorithm,
                          sensealgSS::SciMLSensitivity.AbstractAdjointSensitivityAlgorithm,
                          tol::Float64,
                          odeProb::ODEProblem,  
                          peTabModel::PeTabModel,
                          simulationInfo::SimulationInfo,
                          paramIndices::ParameterIndices,
                          measurementData::MeasurementData,
                          parameterData::ParamData, 
                          changeModelParamUse!::Function, 
                          solveOdeModelAllCondUse!::Function, 
                          priorInfo::PriorInfo;
                          expIDSolve::Array{String, 1} = ["all"])   
    
    # Split input into observeble and dynamic parameters 
    dynamicParamEst::Vector{Float64} = paramVecEst[paramIndices.iDynParam]
    obsParEst::Vector{Float64} = paramVecEst[paramIndices.iObsParam]
    sdParamEst::Vector{Float64} = paramVecEst[paramIndices.iSdParam]
    noneDynParamEst::Vector{Float64} = paramVecEst[paramIndices.iNonDynParam]
    namesSdParam::Vector{String} = paramIndices.namesSdParam
    namesObsParam::Vector{String} = paramIndices.namesObsParam
    namesNonDynParam::Vector{String} = paramIndices.namesNonDynParam
    namesSdObsNonDynPar::Vector{String} = paramIndices.namesSdObsNonDynPar

    # Calculate gradient seperately for dynamic and non dynamic parameter. 

    gradDynParam::Vector{Float64} = zeros(Float64, length(dynamicParamEst))
    calcGradAdjDynParam!(gradDynParam, dynamicParamEst, sdParamEst, obsParEst,
                         noneDynParamEst, odeProb, adjSolver, tol, sensealg, peTabModel, simulationInfo, paramIndices,
                         measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, priorInfo;
                         expIDSolve=expIDSolve, sensealgSS=sensealgSS)
    grad[paramIndices.iDynParam] .= gradDynParam

    # Happens when at least one forward pass fails and I set the gradient to 1e8 
    if all(grad[paramIndices.iDynParam] .== 1e8)
        return 
    end

    # Here it is crucial to account for that obs- and sd parameter can be overlapping. Thus, a name-map 
    # of both Sd and Obs param is used to account for this. This is not a worry for non-dynamic parameters.
    paramNotOdeSys = paramVecEst[paramIndices.iSdObsNonDynPar]
    iSdUse = [findfirst(x -> x == namesSdParam[i], namesSdObsNonDynPar) for i in eachindex(namesSdParam)]
    iObsUse = [findfirst(x -> x == namesObsParam[i],  namesSdObsNonDynPar) for i in eachindex(namesObsParam)]
    iNonDynUse = [findfirst(x -> x == namesNonDynParam[i],  namesSdObsNonDynPar) for i in eachindex(namesNonDynParam)]

    # TODO : Make choice of gradient availble 
    calcCostNonDyn = (x) -> calcLogLikNotSolveODE(dynamicParamEst, x[iSdUse], x[iObsUse], x[iNonDynUse], peTabModel, simulationInfo, paramIndices, measurementData, parameterData, priorInfo, calcGradObsSdParamAdj=true, expIDSolve=expIDSolve)
    @views ReverseDiff.gradient!(grad[paramIndices.iSdObsNonDynPar], calcCostNonDyn, paramNotOdeSys)

    # Account for priors 
    if priorInfo.hasPriors == true
        evalLogLikPrior = (pVec) ->     begin
                                            paramVecEstTransformed = transformParamVec(pVec, paramIndices.namesParamEst, parameterData)
                                            return evalPriors(paramVecEstTransformed, pVec, paramIndices.namesParamEst, paramIndices, priorInfo)
                                        end
        grad .+= ForwardDiff.gradient(evalLogLikPrior, paramVecEst)
    end
end


function calcYBarSS(simulationInfo::SimulationInfo, 
                    sensealgSS::SteadyStateAdjoint, 
                    solver::SciMLAlgorithm,
                    tol::Float64,
                    expIDSolve::Vector{String})::NamedTuple

    # Extract all unique Pre-equlibrium conditions. If the code is run in parallell 
    # (expIDSolve != [["all]]) the number of preEq cond. might be smaller than the 
    # total number of preEq cond.
    if expIDSolve[1] == "all"
        preEqIds = unique(simulationInfo.firstExpIds)
    else
        whichPreEq = findall(x -> x ∈ simulationInfo.conditionIdSol, expIDSolve)
        preEqIds = unique(simulationInfo.preEqIdSol[whichPreEq])
    end

    yBarF = Array{Function, 1}(undef, length(preEqIds))
    for i in eachindex(preEqIds)
        whichPreEq = findfirst(x -> x == preEqIds[i], simulationInfo.preEqIdSol)
        probUse = simulationInfo.solArrayPreEq[whichPreEq].prob
        ssProb = SteadyStateProblem(probUse)
        ySS, ybarSS = Zygote.pullback((p) ->    (
                                                solve(ssProb, 
                                                      DynamicSS(solver, abstol=simulationInfo.absTolSS, reltol=simulationInfo.relTolSS), 
                                                      abstol=tol, 
                                                      reltol=tol, 
                                                      p=p, 
                                                      sensealg=sensealgSS)[:]), probUse.p)
                                                
        yBarF[i] = ybarSS
    end
    yBarTuple = Tuple(f for f in yBarF)

    return NamedTuple{Tuple(Symbol(preEqId) for preEqId in preEqIds)}(yBarTuple)
end
function calcYBarSS(simulationInfo::SimulationInfo, 
                    sensealgSS::Union{QuadratureAdjoint, InterpolatingAdjoint}, 
                    solver::SciMLAlgorithm,
                    tol::Float64,
                    expIDSolve::Vector{String})::NamedTuple

    # Extract all unique Pre-equlibrium conditions. If the code is run in parallell 
    # (expIDSolve != [["all]]) the number of preEq cond. might be smaller than the 
    # total number of preEq cond.
    if expIDSolve[1] == "all"
        preEqIds = unique(simulationInfo.firstExpIds)
    else
        whichPreEq = findall(x -> x ∈ simulationInfo.conditionIdSol, expIDSolve)
        preEqIds = unique(simulationInfo.preEqIdSol[whichPreEq])
    end

    yBarF::Array{Function, 1} = Array{Function, 1}(undef, length(preEqIds))
    for i in eachindex(preEqIds)
        whichPreEq = findfirst(x -> x == preEqIds[i], simulationInfo.preEqIdSol)
        preEqSol = simulationInfo.solArrayPreEq[whichPreEq]
        probSSPullback = remake(preEqSol.prob, tspan=(0.0, preEqSol.t[end]))
        ssVal, ybarSS = Zygote.pullback((p) -> solve(probSSPullback, solver, p=p, abstol=tol, reltol=tol, sensealg=sensealgSS)[:, end], preEqSol.prob.p)
                                                
        yBarF[i] = ybarSS
    end
    yBarTuple = Tuple(f for f in yBarF)

    return NamedTuple{Tuple(Symbol(preEqId) for preEqId in preEqIds)}(yBarTuple)
end


# Compute the adjoint gradient across all experimental conditions 
function calcGradAdjDynParam!(gradAdj::Vector{Float64},
                              dynamicParamEst::Vector{Float64},
                              sdParamEst::Vector{Float64},
                              obsParEst::Vector{Float64},
                              nonDynParamEst::Vector{Float64},
                              odeProb::ODEProblem,
                              adjSolver::SciMLAlgorithm, 
                              tol::Float64,
                              sensealg::SciMLSensitivity.AbstractAdjointSensitivityAlgorithm,
                              peTabModel::PeTabModel,
                              simulationInfo::SimulationInfo,
                              paramIndices::ParameterIndices,
                              measurementData ::MeasurementData, 
                              parameterData::ParamData, 
                              changeModelParamUse!::Function,
                              solveOdeModelAllCondUse!::Function, 
                              priorInfo::PriorInfo;
                              sensealgSS=SteadyStateAdjoint(),
                              expIDSolve::Array{String, 1} = ["all"])

    
    dynamicParamEstUse = transformParamVec(dynamicParamEst, paramIndices.namesDynParam, parameterData)
    sdParamEstUse = transformParamVec(sdParamEst, paramIndices.namesSdParam, parameterData)
    obsParEstUse = transformParamVec(obsParEst, paramIndices.namesObsParam, parameterData)
    nonDynParamEstUse = transformParamVec(nonDynParamEst, paramIndices.namesNonDynParam, parameterData)

    odeProbUse = remake(odeProb, p = convert.(eltype(dynamicParamEstUse), odeProb.p), u0 = convert.(eltype(dynamicParamEstUse), odeProb.u0))
    changeModelParamUse!(odeProbUse.p, odeProbUse.u0, dynamicParamEstUse)
    success = solveOdeModelAllCondUse!(simulationInfo.solArrayGrad, odeProbUse, dynamicParamEstUse, expIDSolve)
    if success != true
        grad .= 1e8
        return
    end

    # In case of PreEq-critera we need to compute the pullback function at tSS to compute the VJP between 
    # λ_t0 and the sensitivites at t_ss
    if simulationInfo.simulateSS == true
        yBarSS = calcYBarSS(simulationInfo, sensealgSS, adjSolver, tol, expIDSolve)
    else
        yBarSSUse = identity
    end

    gradAdj .= 0.0
    # Compute the gradient by looping through all experimental conditions.
    for conditionID in keys(measurementData.iPerConditionId)

        if expIDSolve[1] != "all" && conditionID ∉ expIDSolve
            continue
        end

        whichForwardSol = findfirst(x -> x == conditionID, simulationInfo.conditionIdSol)
        if simulationInfo.simulateSS == true
            yBarSSUse = yBarSS[Symbol(simulationInfo.preEqIdSol[whichForwardSol])]
        end

        whichForwardSol = findfirst(x -> x == conditionID, simulationInfo.conditionIdSol)
        sol = simulationInfo.solArrayGrad[whichForwardSol]
        calcGradAdjExpCond!(gradAdj, sol, sensealg, tol, adjSolver, dynamicParamEstUse,
                            sdParamEstUse, obsParEstUse, nonDynParamEstUse, conditionID, 
                            peTabModel, paramIndices, measurementData, parameterData, 
                            simulationInfo.simulateSS, yBarSSUse)
    end

end


# For a given experimental condition compute the gradient using adjoint sensitivity analysis 
function calcGradAdjExpCond!(grad::Vector{Float64},
                             sol::ODESolution,
                             sensealg::SciMLSensitivity.AbstractAdjointSensitivityAlgorithm,
                             tol::Float64,
                             adjSolver::SciMLAlgorithm, 
                             dynParam::Vector{Float64},
                             sdParam::Vector{Float64}, 
                             obsParam::Vector{Float64}, 
                             nonDynParam::Vector{Float64},
                             conditionID::String,
                             peTabModel::PeTabModel,
                             paramIndices::ParameterIndices,
                             measurementData::MeasurementData, 
                             parameterData::ParamData, 
                             simulateSS::Bool,
                             yBarSS::Function;
                             expIDSolve::Array{String, 1} = ["all"])
                            
    if expIDSolve[1] != "all" && conditionID ∉ expIDSolve
        return
    end

    iGroupedTObs = measurementData.iGroupedTObs[conditionID]
    # Pre allcoate vectors needed for computations 
    dYmodDu = zeros(Float64, length(sol.prob.u0))
    dSdDu = zeros(Float64, length(sol.prob.u0))
    dYmodDp = zeros(Float64, length(sol.prob.p))
    dSdDp = zeros(Float64, length(sol.prob.p))

    # Functions needed by the lower level interface 
    calcDgDuDiscrete = (out, u, p, t, i) -> begin calcdGd_Discrete(out, u, p, t, i, iGroupedTObs, 
                                                                   measurementData, parameterData, 
                                                                   paramIndices, peTabModel, 
                                                                   dynParam, sdParam, obsParam, nonDynParam, 
                                                                   dYmodDu, dSdDu, calcdGdU=true)
                                            end
    calcDgDpDiscrete = (out, u, p, t, i) -> begin calcdGd_Discrete(out, u, p, t, i, iGroupedTObs, 
                                                                   measurementData, parameterData, 
                                                                   paramIndices, peTabModel, 
                                                                   dynParam, sdParam, obsParam, nonDynParam, 
                                                                   dYmodDp, dSdDp, calcdGdU=false)
                                                                                            end
                                                                                        
    # Time points for which we have observed data. The standard allow cases where we only observe data at t0, that 
    # is we do not solve the ODE. Here adjoint_sensitivities fails (naturally). In this case we compute the gradient 
    # via ∇G_p = dp + du*J(u(t_0)) where du is the cost function differentiated with respect to the states at time zero, 
    # dp is the cost function differentiated with respect to the parameters at time zero and J is sensititvites at time 
    # zero. Overall, the only workflow that changes below is that we compute du outside of the adjoint interface 
    # and use sol[:] as we no longer can interpolate from the forward solution.
    tSaveAt = measurementData.tVecSave[conditionID]
    onlyObsAtZero::Bool = false
    if !(length(tSaveAt) == 1 && tSaveAt[1] == 0.0)
        du, dp = adjoint_sensitivities(sol, 
                                    adjSolver,
                                    dgdp_discrete=nothing,
                                    dgdu_discrete=calcDgDuDiscrete, 
                                    t=tSaveAt, 
                                    sensealg=sensealg, 
                                    abstol=tol, 
                                    reltol=tol)
    else
        du = zeros(Float64, length(sol.prob.u0))
        calcDgDuDiscrete(du, sol[1], sol.prob.p, 0.0, 1)
        dp = zeros(Float64, length(sol.prob.p))'
        onlyObsAtZero = true
    end
    # Technically we can pass calcDgDpDiscrete above to dgdp_discrete. However, odeProb.p often contain 
    # constant parameters which are not a part ode the parameter estimation problem. Sometimes 
    # the gradient for these evaluate to NaN (as they where never thought to be estimated) which 
    # results in the entire gradient evaluating to NaN. Hence, we perform this calculation outside 
    # of the lower level interface. 
    dgDpOut = zeros(Float64, length(sol.prob.p))
    for i in eachindex(tSaveAt)     
        if onlyObsAtZero == false                                                                                   
            calcDgDpDiscrete(dgDpOut, sol(tSaveAt[i]), sol.prob.p, tSaveAt[i], i)
        else
            calcDgDpDiscrete(dgDpOut, sol[1], sol.prob.p, tSaveAt[i], i)
        end
        dp .+= dgDpOut'
    end

    if simulateSS == false
        # In case we do not simulate the ODE for a steady state first we can compute 
        # the initial sensitivites easily via automatic differantitatiom
        sMatAtT0::Matrix{Float64} = Matrix{Float64}(undef, (length(sol.prob.u0), length(sol.prob.p)))
        ForwardDiff.jacobian!(sMatAtT0, peTabModel.evalU0, sol.prob.p)
        gradTot = dp .+ du'*sMatAtT0

    else
        # In case we simulate to a stady state we need to compute a VJP. We use 
        # Zygote pullback to avoid having to having build the Jacobian, rather 
        # we create the yBar function required for the vector Jacobian product.
        gradTot = (dp .+ (yBarSS(du)[1])')[:]
    end

    # Thus far have have computed dY/dθ, but for parameters on the log-scale we 
    # want dY/dθ_log. We can adjust via;
    # dY/dθ_log = log(10) * θ * dY/dθ
    # Account for parameters which are not condition specific 
    grad[paramIndices.mapDynParEst.iDynParamInVecEst] .+= transformParamVecGrad(gradTot[paramIndices.mapDynParEst.iDynParamInSys], 
                                                                                dynParam[paramIndices.mapDynParEst.iDynParamInVecEst], 
                                                                                paramIndices.namesDynParam[paramIndices.mapDynParEst.iDynParamInVecEst], 
                                                                                parameterData)
    # For parameters which are specific to an experimental condition 
    whichExpMap = findfirst(x -> x == conditionID, [paramIndices.mapExpCond[i].condID for i in eachindex(paramIndices.mapExpCond)])
    expMap = paramIndices.mapExpCond[whichExpMap]                                          
    grad[expMap.iDynEstVec] .+= transformParamVecGrad(gradTot[expMap.iOdeProbDynParam], 
                                                      dynParam[expMap.iDynEstVec], 
                                                      paramIndices.namesDynParam[expMap.iDynEstVec], 
                                                      parameterData)                                   
end


function calcdGd_Discrete(out,
                          u::Vector{Float64}, 
                          p::Vector{Float64}, # odeProb.p
                          t::Float64, 
                          i::Integer, 
                          iGroupedTObs::Vector{Vector{Int64}}, 
                          measurementData::MeasurementData, 
                          parameterData::ParamData,
                          paramIndices::ParameterIndices,
                          peTabModel::PeTabModel, 
                          dynParam::Vector{Float64},
                          sdParam::Vector{Float64}, 
                          obsParam::Vector{Float64}, 
                          nonDynParam::Vector{Float64}, 
                          dYmodD_::Vector{Float64}, 
                          dSdD_::Vector{Float64}; 
                          calcdGdU::Bool=true)

    out .= 0.0
    for iMeasurementData in iGroupedTObs[i]
        dYmodD_ .= 0.0
        dSdD_ .= 0.0

        mapObsParam = paramIndices.mapArrayObsParam[paramIndices.indexObsParamMap[iMeasurementData]]
        yMod = peTabModel.evalYmod(u, t, dynParam, obsParam, nonDynParam, parameterData, measurementData.observebleID[iMeasurementData], mapObsParam)
        yModTrans = transformObsOrData(yMod, measurementData.transformData[iMeasurementData])

        # Compute associated SD-value or extract said number if it is known 
        mapSdParam = paramIndices.mapArraySdParam[paramIndices.indexSdParamMap[iMeasurementData]]
        if typeof(measurementData.sdParams[iMeasurementData]) <: AbstractFloat
            sdVal = measurementData.sdParams[iMeasurementData]
        else
            sdVal = peTabModel.evalSd!(u, t, sdParam, dynParam, nonDynParam, parameterData, measurementData.observebleID[iMeasurementData], mapSdParam)
        end

        if calcdGdU == true
            peTabModel.evalDYmodDu(u, t, p, obsParam, nonDynParam, measurementData.observebleID[iMeasurementData], mapObsParam, dYmodD_)
            peTabModel.evalDSdDu!(u, t, sdParam, p, nonDynParam, parameterData, measurementData.observebleID[iMeasurementData], mapSdParam, dSdD_)
        else
            peTabModel.evalDYmodDp(u, t, p, obsParam, nonDynParam, measurementData.observebleID[iMeasurementData], mapObsParam, dYmodD_)
            peTabModel.evalDSdDp!(u, t, sdParam, p, nonDynParam, parameterData, measurementData.observebleID[iMeasurementData], mapSdParam, dSdD_)
        end

        if measurementData.transformData[iMeasurementData] == :log10
            yObs = measurementData.yObsTransformed[iMeasurementData]
            dYmodD_ .*= 1 / (log(10) * exp10(yModTrans))
        elseif measurementData.transformData[iMeasurementData] == :lin
            yObs = measurementData.yObsNotTransformed[iMeasurementData]
        end

        dLogLikDyModTrans::Float64 = ( yModTrans - yObs ) / sdVal^2
        dLogLikDsdVal::Float64 = 1/sdVal - (( yModTrans - yObs )^2 / sdVal^3)
        to_out = dLogLikDyModTrans*dYmodD_ .+ dLogLikDsdVal*dSdD_
        @inbounds for i in eachindex(to_out)
            out[i] += to_out[i]
        end
    end
    return
end


function calcdGDiscrete(u::AbstractVector,
                        p::AbstractVector, # odeProb.p
                        t::AbstractFloat, 
                        i::Integer, 
                        iGroupedTObs::Vector{Vector{Int64}}, 
                        measurementData::MeasurementData, 
                        parameterData::ParamData,
                        paramIndices::ParameterIndices,
                        peTabModel::PeTabModel, 
                        sdParam::Vector{Float64}, 
                        obsParam::Vector{Float64}, 
                        nonDynParam::Vector{Float64})

    dynParam = p[paramIndices.mapDynParEst.iDynParamInSys]
    logLik = 0.0

    for iMeasurementData in iGroupedTObs[i]
        
        mapObsParam = paramIndices.mapArrayObsParam[paramIndices.indexObsParamMap[iMeasurementData]]
        yMod = peTabModel.evalYmod(u, t, dynParam, obsParam, nonDynParam, parameterData, measurementData.observebleID[iMeasurementData], mapObsParam)
        yModTrans = transformObsOrData(yMod, measurementData.transformData[iMeasurementData])

        # Compute associated SD-value or extract said number if it is known 
        mapSdParam = paramIndices.mapArraySdParam[paramIndices.indexSdParamMap[iMeasurementData]]
        if typeof(measurementData.sdParams[iMeasurementData]) <: AbstractFloat
            sdVal = measurementData.sdParams[iMeasurementData]
        else
            sdVal = peTabModel.evalSd!(u, t, sdParam, dynParam, nonDynParam, parameterData, measurementData.observebleID[iMeasurementData], mapSdParam)
        end

        # Update log-likelihood 
        if measurementData.transformData[iMeasurementData] == :lin
            logLik += log(sdVal) + 0.5*log(2*pi) + 0.5*((yModTrans - measurementData.yObsNotTransformed[iMeasurementData]) / sdVal)^2
        elseif measurementData.transformData[iMeasurementData] == :log10
            logLik += log(sdVal) + 0.5*log(2*pi) + log(log(10)) + log(exp10(measurementData.yObsTransformed[iMeasurementData])) + 0.5*( ( log(exp10(yModTrans)) - log(exp10(measurementData.yObsTransformed[iMeasurementData])) ) / (log(10)*sdVal))^2
        else
            println("Transformation ", measurementData.transformData[iMeasurementData], "not yet supported.")
            return Inf
        end   
        
    end
    return logLik
end


function calcLogLikExpCond(odeSol::ODESolution,
                           dynamicParamEst::AbstractVector,
                           sdParamEst::AbstractVector, 
                           obsPar::AbstractVector, 
                           nonDynParamEst::AbstractVector,
                           peTabModel::PeTabModel,
                           conditionID::String,
                           paramIndices::ParameterIndices,
                           measurementData::MeasurementData,
                           parameterData::ParamData, 
                           calcGradObsSdParam::Bool;
                           calcGradObsSdParamAdj=false)::Real

    if !(odeSol.retcode == :Success || odeSol.retcode == :Terminated)
        return Inf
    end

    # Compute yMod and sd for all observations having id conditionID 
    logLik = 0.0
    for i in measurementData.iPerConditionId[conditionID]
        # Compute Y-mod value 
        t = measurementData.tObs[i]
        if calcGradObsSdParam == true
            odeSolAtT = dualVecToFloatVec(odeSol[:, measurementData.iTObs[i]]) #u
        elseif calcGradObsSdParamAdj == true
            # In case we only have sol.t = 0.0 (or similar) interpolation does not work
            if length(odeSol.t) > 1 
                odeSolAtT = odeSol(t)
            else
                odeSolAtT = odeSol[1]
            end
        else
            odeSolAtT = odeSol[:, measurementData.iTObs[i]]
        end
        mapObsParam = paramIndices.mapArrayObsParam[paramIndices.indexObsParamMap[i]]
        yMod = peTabModel.evalYmod(odeSolAtT, t, dynamicParamEst, obsPar, nonDynParamEst, parameterData, measurementData.observebleID[i], mapObsParam) 

        # Compute associated SD-value or extract said number if it is known 
        if typeof(measurementData.sdParams[i]) <: AbstractFloat
            sdVal = measurementData.sdParams[i]
        else
            mapSdParam = paramIndices.mapArraySdParam[paramIndices.indexSdParamMap[i]]
            sdVal = peTabModel.evalSd!(odeSolAtT, t, sdParamEst, dynamicParamEst, nonDynParamEst, parameterData, measurementData.observebleID[i], mapSdParam)
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


# Evaluate contribution of potential prior to the final cost function value. Not mutating so works with both Zygote 
# and forwardiff.
function evalPriors(paramVecTransformed, 
                    paramVecNotTransformed,
                    namesParamVec::Array{String, 1}, 
                    paramEstIndices::ParameterIndices, 
                    priorInfo::PriorInfo)::Real

    if priorInfo.hasPriors == false
        return 0.0
    end

    k = 0
    priorContribution = 0.0
    for i in eachindex(paramVecNotTransformed)
        iParam = findfirst(x -> x == namesParamVec[i], paramEstIndices.namesParamEst)
        if priorInfo.priorOnParamScale[iParam] == true
            pInput = paramVecNotTransformed[i]
        else
            pInput = paramVecTransformed[i]
        end
        priorContribution += priorInfo.logpdf[iParam](pInput)

        if priorInfo.logpdf[iParam](pInput) != 0
            k += 1
        end
    end

    return priorContribution
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
function changeModelParam!(pOdeSys::AbstractVector, 
                           u0::AbstractVector,
                           paramVecEst::AbstractVector,
                           paramIndices::ParameterIndices,
                           peTabModel::PeTabModel)

    mapDynParam = paramIndices.mapDynParEst
    pOdeSys[mapDynParam.iDynParamInSys] .= paramVecEst[mapDynParam.iDynParamInVecEst]
    peTabModel.evalU0!(u0, pOdeSys) 
    
    return nothing
end


function changeModelParam(pOdeSys::AbstractVector,
                          paramVecEst::AbstractVector,
                          paramIndices::ParameterIndices,
                          peTabModel::PeTabModel)

    # Helper function to not-inplace map parameters 
    function mapParamToEst(iUse::Integer, mapDynParam::MapDynParEst)
        whichIndex = findfirst(x -> x == iUse, mapDynParam.iDynParamInSys)
        return mapDynParam.iDynParamInVecEst[whichIndex]
    end

    mapDynParam = paramIndices.mapDynParEst
    pOdeSysRet = [i ∈ mapDynParam.iDynParamInSys ? paramVecEst[mapParamToEst(i, mapDynParam)] : pOdeSys[i] for i in eachindex(pOdeSys)]
    u0Ret = peTabModel.evalU0(pOdeSysRet) 
    
    return pOdeSysRet, u0Ret   
end


"""
    transformParamVec!(paramVec, namesParam::Array{String, 1}, paramData::ParamData; revTransform::Bool=false)

    Helper function which transforms in-place a parameter vector with parameters specied in namesParam according to the 
    transformation for said parameter specifid in paramData.shouldTransform. In case revTransform is true 
    performs the inverse parameter transformation (e.g exp10 instead of log10)
"""
function transformParamVec!(paramVec::AbstractVector, 
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
    performs the inverse parameter transformation (e.g exp10 instead of log10). 

    The function is fully compatible with Zygote.
"""
function transformParamVec(paramVec::AbstractVector, 
                           namesParam::Array{String, 1}, 
                           paramData::ParamData; 
                           revTransform::Bool=false)::AbstractVector
    
    iParam = [findfirst(x -> x == namesParam[i], paramData.parameterID) for i in eachindex(namesParam)]
    shouldTransform = [paramData.logScale[i] == true ? true : false for i in iParam]
    shouldNotTransform = .!shouldTransform

    if revTransform == false
        return exp10.(paramVec) .* shouldTransform .+ paramVec .* shouldNotTransform
    else
        return log10.(paramVec) .* shouldTransform .+ paramVec .* shouldNotTransform
    end
end


function transformParamVecGrad(gradientVec::AbstractVector, 
                               dynamicParam::AbstractVector,
                               namesParam::Array{String, 1}, 
                               paramData::ParamData)::AbstractVector
    
    iParam = [findfirst(x -> x == namesParam[i], paramData.parameterID) for i in eachindex(namesParam)]
    shouldTransform = [paramData.logScale[i] == true ? true : false for i in iParam]
    shouldNotTransform = .!shouldTransform

    return log(10) .* gradientVec .* dynamicParam .* shouldTransform .+ gradientVec .* shouldNotTransform
end


function dualVecToFloatVec(dualVec::T1)::Vector{Float64} where T1<:Vector{<:ForwardDiff.Dual}
    return [dualVec[i].value for i in eachindex(dualVec)]
end
