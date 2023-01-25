# PEtab structs 
include(joinpath(pwd(), "src", "PeTab_structs.jl"))

# Files related to computing the cost (likelihood)
include(joinpath(pwd(), "src", "Compute_cost", "Compute_priors.jl"))
include(joinpath(pwd(), "src", "Compute_cost", "Compute_cost.jl"))
include(joinpath(pwd(), "src", "Compute_cost", "Compute_cost_zygote.jl"))

# Files related to computing derivatives 
include(joinpath(pwd(), "src", "Derivatives", "Hessian.jl"))
include(joinpath(pwd(), "src", "Derivatives", "Gradient.jl"))
include(joinpath(pwd(), "src", "Derivatives", "Adjoint_sensitivity_analysis.jl"))
include(joinpath(pwd(), "src", "Derivatives", "Forward_sensitivity_equations.jl"))
include(joinpath(pwd(), "src", "Derivatives", "Gauss_newton.jl"))
include(joinpath(pwd(), "src", "Derivatives", "Common.jl"))

# Files related to solving the ODE-system 
include(joinpath(pwd(), "src", "Solve_ODE", "Change_experimental_condition.jl"))
include(joinpath(pwd(), "src", "Solve_ODE", "Solve_ode_Zygote.jl"))
include(joinpath(pwd(), "src", "Solve_ODE", "Solve_ode_model.jl"))

# Files related to distributed computing 
include(joinpath(pwd(), "src", "Distributed", "Distributed.jl"))

# Files related to processing PEtab files 
include(joinpath(pwd(), "src", "Process_PEtab_files", "Common.jl"))
include(joinpath(pwd(), "src", "Process_PEtab_files", "Get_simulation_info.jl"))
include(joinpath(pwd(), "src", "Process_PEtab_files", "Get_parameter_indices.jl"))
include(joinpath(pwd(), "src", "Process_PEtab_files", "Process_measurements.jl"))
include(joinpath(pwd(), "src", "Process_PEtab_files", "Process_parameters.jl"))
include(joinpath(pwd(), "src", "Process_PEtab_files", "Create_obs_u0_sd_common.jl"))
include(joinpath(pwd(), "src", "Process_PEtab_files", "Create_obs_u0_sd_functions.jl"))
include(joinpath(pwd(), "src", "Process_PEtab_files", "Create_obs_u0_sd_derivatives.jl"))
include(joinpath(pwd(), "src", "Process_PEtab_files", "Process_PeTab_files.jl"))
include(joinpath(pwd(), "src", "Common.jl"))


"""
    setUpPeTabModel(modelName::String, dirModel::String)::PeTabModel

    Given a model directory (dirModel) containing the PeTab files and a 
    xml-file on format modelName.xml will return a PeTabModel struct holding 
    paths to PeTab files, ode-system in ModellingToolkit format, functions for 
    evaluating yMod, u0 and standard deviations, and a parameter and state maps 
    for how parameters and states are mapped in the ModellingToolkit ODE system
    along with state and parameter names.

    dirModel must contain a SBML file named modelName.xml, and files starting with 
    measurementInfo, experimentalCondition, parameter, and observables (tsv-files).
    The latter files must be unique (e.g only one file starting with measurementInfo)

    TODO : Example  
"""
function setUpPeTabModel(modelName::String, dirModel::String; forceBuildJlFile::Bool=false, verbose::Bool=true, ifElseToEvent=true)::PeTabModel

    # Sanity check user input 
    modelFileXml = dirModel * modelName * ".xml"
    modelFileJl = dirModel * modelName * ".jl"
    if !isdir(dirModel)
        if verbose
            @printf("Model directory %s does not exist\n", dirModel)
        end
    end
    if !isfile(modelFileXml)
        if verbose
            @printf("Model directory does not contain xml-file with name %s\n", modelName * "xml")
        end
    end
    # If Julia model file does exists build it 
    if !isfile(modelFileJl) && forceBuildJlFile == false
        if verbose
            @printf("Julia model file does not exist - will build it\n")
        end
        modelDict = XmlToModellingToolkit(modelFileXml, modelName, dirModel, ifElseToEvent=ifElseToEvent)
    elseif isfile(modelFileJl) && forceBuildJlFile == false
        if verbose
            @printf("Julia model file exists at %s - will not rebuild it\n", modelFileJl)
        end
    elseif forceBuildJlFile == true
        if verbose
            @printf("By user option rebuilds Julia model file\n")
        end
        if isfile(modelFileJl)
            rm(modelFileJl)
        end
        modelDict = XmlToModellingToolkit(modelFileXml, modelName, dirModel, ifElseToEvent=ifElseToEvent)
    end

    # Extract ODE-system and mapping of maps of how to map parameters to states and model parmaeters 
    include(modelFileJl)
    expr = Expr(:call, Symbol("getODEModel_" * modelName))
    odeSys, stateMap, paramMap = eval(expr)
    #odeSysUse = ode_order_lowering(odeSys)
    odeSysUse = structural_simplify(odeSys)
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
    pathDObsSdU0 = dirModel * modelName * "DObsSdU0.jl"
    if !isfile(pathObsSdU0) || forceBuildJlFile == true
        if verbose && forceBuildJlFile == false
            @printf("File for yMod, U0 and Sd does not exist - building it\n")
        end
        if verbose && forceBuildJlFile == true
            @printf("By user option will rebuild Ymod, Sd and u0\n")
        end
        if !@isdefined(modelDict)
            modelDict = XmlToModellingToolkit(modelFileXml, modelName, dirModel, writeToFile=false, ifElseToEvent=ifElseToEvent)
        end
        createFileYmodSdU0(modelName, dirModel, odeSysUse, stateMap, modelDict)
        createFileDYmodSdU0(modelName, dirModel, odeSysUse, stateMap, modelDict)
    else
        if verbose
            @printf("File for yMod, U0 and Sd does exist - will not rebuild it\n")
        end
    end
    include(pathObsSdU0)
    include(pathDObsSdU0)    

    pathCallback = dirModel * "/" * modelName * "Callbacks_time_piecewise.jl"
    if !isfile(pathCallback) || forceBuildJlFile == true
        if verbose && forceBuildJlFile == false
            @printf("File for callback does not exist - building it\n")
        end
        if verbose && forceBuildJlFile == true
            @printf("By user option will rebuild callback file\n")
        end
        if !@isdefined(modelDict)
            modelDict = XmlToModellingToolkit(modelFileXml, modelName, dirModel, writeToFile=false, ifElseToEvent=ifElseToEvent)
        end
        getCallbacksForTimeDepedentPiecewise(odeSysUse, modelDict, modelName, dirModel)
    end
    include(pathCallback)
    exprCallback = Expr(:call, Symbol("getCallbacks_" * modelName))
    cbSet::CallbackSet, checkCbActive::Vector{Function} = eval(exprCallback)    

    peTabModel = PeTabModel(modelName,
                            evalYmod,
                            evalU0!,
                            evalU0,
                            evalSd!,
                            evalDYmodDu,
                            evalDSdDu!,
                            evalDYmodDp,
                            evalDSdDp!,
                            getTstops,
                            odeSysUse,
                            paramMap,
                            stateMap,
                            parameterNames, 
                            stateNames,
                            dirModel,
                            pathMeasurementData,
                            pathExperimentalCond,
                            pathObservables, 
                            pathParameters, 
                            cbSet, 
                            checkCbActive)

    return peTabModel
end


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
                           sensealgForward::Union{Symbol, SciMLSensitivity.AbstractForwardSensitivityAlgorithm}=ForwardSensitivity(),
                           sparseJac::Bool=false, 
                           absTolSS::Float64=1e-8, 
                           relTolSS::Float64=1e-6, 
                           nProcs::Signed=1,
                           adjSolver=Rodas5P(), 
                           solverForward=Rodas5P(autodiff=false),
                           adjTol::Float64=1e-6, 
                           adjSensealg::SciMLSensitivity.AbstractAdjointSensitivityAlgorithm=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(false)), 
                           adjSensealgSS::SciMLSensitivity.AbstractAdjointSensitivityAlgorithm=SteadyStateAdjoint())::PeTabOpt

    if !(typeof(adjSensealgSS) <: SteadyStateAdjoint)
        println("If you are using adjoint sensitivity analysis for a model with PreEq-criteria the most the most efficient adjSensealgSS is usually SteadyStateAdjoint. The algorithm you have provided, ", adjSensealgSS, "might not work (as there are some bugs here). In case it does not work, and SteadyStateAdjoint fails (because a dependancy on time) a good choice might be QuadratureAdjoint(autodiff=false, autojacvec=false)")
    end

    # Process PeTab files into type-stable Julia structs 
    experimentalConditionsFile, measurementDataFile, parameterDataFile, observablesDataFile = readPEtabFiles(peTabModel.dirModel, readObservables=true)
    parameterData = processParameters(parameterDataFile)
    measurementInfo = processMeasurements(measurementDataFile, observablesDataFile) 
    simulationInfo = processSimulationInfo(peTabModel, measurementInfo, sensealg=sensealg, absTolSS=absTolSS, relTolSS=relTolSS)

    # Indices for mapping parameter-estimation vector to dynamic, observable and sd parameters correctly when calculating cost
    paramEstIndices = computeIndicesθ(parameterData, measurementInfo, peTabModel.odeSystem, experimentalConditionsFile)
    
    # Set up potential prior for the parameters to estimate 
    priorInfo::PriorInfo = processPriors(paramEstIndices, parameterDataFile)

    # Set model parameter values to those in the PeTab parameter data ensuring correct value of constant parameters 
    setParamToFileValues!(peTabModel.paramMap, peTabModel.stateMap, parameterData)

    # The time-span 5e3 is overwritten when performing actual forward simulations 
    odeProb = ODEProblem(peTabModel.odeSystem, peTabModel.stateMap, (0.0, 5e3), peTabModel.paramMap, jac=true, sparse=sparseJac)
    odeProb = remake(odeProb, p = convert.(Float64, odeProb.p), u0 = convert.(Float64, odeProb.u0))
    if sensealgForward == :AutoDiffForward
        odeProbSenseEq = deepcopy(odeProb)
    else
        odeProbSenseEq = ODEForwardSensitivityProblem(odeProb.f, odeProb.u0, odeProb.tspan, odeProb.p, 
                                                      sensealg=sensealgForward)
    end

    # Functions to map experimental conditions and parameters correctly to the ODE model 
    changeToExperimentalCondUse! = (pVec, u0Vec, expID, dynParamEst) -> _changeExperimentalCondition!(pVec, u0Vec, expID, dynParamEst, peTabModel, paramEstIndices)
    changeToExperimentalCondSenseEqUse! = (pVec, u0Vec, expID, dynParamEst) -> _changeExperimentalCondition!(pVec, u0Vec, expID, dynParamEst, peTabModel, paramEstIndices, computeForwardSensitivites=true)
    changeToExperimentalCondUse = (pVec, u0Vec, expID, dynParamEst) -> _changeExperimentalCondition(pVec, u0Vec, expID, dynParamEst, peTabModel, paramEstIndices)
    changeModelParamUse! = (pVec, u0Vec, paramEst) -> changeODEProblemParameters!(pVec, u0Vec, paramEst, paramEstIndices, peTabModel)
    changeModelParamUse = (pVec, paramEst) -> changeODEProblemParameters(pVec, paramEst, paramEstIndices, peTabModel)

    # Set up function which solves the ODE model for all conditions and stores result 
    solveOdeModelAllCondUse! = (solArrayArg, odeProbArg, dynParamEst, expIDSolveArg) -> solveODEAllExperimentalConditions!(solArrayArg, odeProbArg, dynParamEst, changeToExperimentalCondUse!, simulationInfo, solver, tol, tol, peTabModel.getTStops, onlySaveAtObservedTimes=true, expIDSolve=expIDSolveArg)
    solveOdeModelAllCondAdjUse! = (solArrayArg, odeProbArg, dynParamEst, expIDSolveArg) -> solveODEAllExperimentalConditions!(solArrayArg, odeProbArg, dynParamEst, changeToExperimentalCondUse!, simulationInfo, solver, tol, tol, peTabModel.getTStops, denseSolution=true, expIDSolve=expIDSolveArg, trackCallback=true)
    solveOdeModelAtCondZygoteUse = (odeProbArg, conditionId, dynParamEst, t_max) -> solveOdeModelAtExperimentalCondZygote(odeProbArg, conditionId, dynParamEst, t_max, changeToExperimentalCondUse, measurementInfo, simulationInfo, solver, tol, tol, sensealg, peTabModel.getTStops)
    if sensealgForward == :AutoDiffForward
        solveOdeModelAllCondForwardEq! = (solArrayArg, SMat, odeProbArg, dynParamEst, expIDSolveArg) -> solveODEAllExperimentalConditions!(solArrayArg, SMat, odeProbArg, dynParamEst, changeToExperimentalCondUse!, changeModelParamUse!, simulationInfo, solverForward, tol, tol, peTabModel.getTStops, onlySaveAtObservedTimes=true, expIDSolve=expIDSolveArg)                                           
    else
        solveOdeModelAllCondForwardEq! = (solArrayArg, odeProbArg, dynParamEst, expIDSolveArg) -> solveODEAllExperimentalConditions!(solArrayArg, odeProbArg, dynParamEst, changeToExperimentalCondSenseEqUse!, simulationInfo, solverForward, tol, tol, peTabModel.getTStops, onlySaveAtObservedTimes=true, expIDSolve=expIDSolveArg)
    end
    solveOdeModelAllCondGuassNewtonForwardEq! = (solArrayArg, SMat, odeProbArg, dynParamEst, expIDSolveArg) -> solveODEAllExperimentalConditions!(solArrayArg, SMat, odeProbArg, dynParamEst, changeToExperimentalCondUse!, changeModelParamUse!, simulationInfo, solverForward, tol, tol, peTabModel.getTStops, onlySaveAtObservedTimes=true, expIDSolve=expIDSolveArg)                                           

    if nProcs > 1 && nprocs() != nProcs
        println("Error : PEtab importer was set to build the cost, grad and hessian with $nProcs processes, 
                 however, Julia is currently running with ", nprocs(), " processes which does not match input 
                 value. Input argument nProcs must match nprocs()")
    elseif nProcs == 1
        evalF = (paramVecEst) -> computeCost(paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementInfo, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, priorInfo)
        evalGradF = (grad, paramVecEst) -> computeGradientAutoDiff!(grad, paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementInfo, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, priorInfo)    
        evalGradFAdjoint = (grad, paramVecEst) -> computeGradientAdjointEquations!(grad, paramVecEst, adjSolver, adjSensealg, adjSensealgSS, adjTol, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementInfo, parameterData, changeModelParamUse!, solveOdeModelAllCondAdjUse!, priorInfo) 
        evalGradFForwardEq = (grad, paramVecEst) -> computeGradientForwardEquations!(grad, paramVecEst, peTabModel, odeProbSenseEq, sensealgForward, simulationInfo, paramEstIndices, measurementInfo, parameterData, changeModelParamUse!, solveOdeModelAllCondForwardEq!, priorInfo) 
        evalHessApprox = (hessianMat, paramVecEst) -> computeHessianBlockApproximation!(hessianMat, paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementInfo, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, priorInfo)
        evalHess = (hessianMat, paramVecEst) -> computeHessian(hessianMat, paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementInfo, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, priorInfo)
    elseif nProcs > 1 && nprocs() == nProcs
        evalF, evalGradF, evalGradFForwardEq, evalGradFAdjoint, evalHess, evalHessApprox = setUpPEtabOptDistributed(peTabModel, solver, tol, 
                                                                                                                    adjSolver, adjSensealg, adjSensealgSS, adjTol,
                                                                                                                    solverForward, sensealgForward, 
                                                                                                                    parameterData, measurementInfo, 
                                                                                                                    simulationInfo, paramEstIndices, priorInfo, odeProb)
    end
    evalFZygote = (paramVecEst) -> computeCostZygote(paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementInfo, parameterData, changeModelParamUse, solveOdeModelAtCondZygoteUse, priorInfo)
    evalGradFZygote = (grad, paramVecEst) -> computeGradientZygote(grad, paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementInfo, parameterData, changeModelParamUse, solveOdeModelAtCondZygoteUse, priorInfo)
    evalHessGaussNewton = (hessian, paramVecEst) -> computeGaussNewtonHessianApproximation!(hessian, paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementInfo, parameterData, changeModelParamUse!, solveOdeModelAllCondGuassNewtonForwardEq!, priorInfo)      

    # Lower and upper bounds for parameters to estimate 
    namesParamEst = paramEstIndices.θ_estNames
    lowerBounds = [parameterData.lowerBound[findfirst(x -> x == namesParamEst[i], parameterData.parameterId)] for i in eachindex(namesParamEst)] 
    upperBounds = [parameterData.upperBound[findfirst(x -> x == namesParamEst[i], parameterData.parameterId)] for i in eachindex(namesParamEst)] 
    # Parameter with nominal values in PeTab file 
    paramVecNominal = [parameterData.nominalValue[findfirst(x -> x == namesParamEst[i], parameterData.parameterId)] for i in eachindex(namesParamEst)]

    # Transform upper and lower bounds if the case 
    transformθ!(lowerBounds, namesParamEst, parameterData, reverseTransform=true)
    transformθ!(upperBounds, namesParamEst, parameterData, reverseTransform=true)
    paramVecNominalTransformed = transformθ(paramVecNominal, namesParamEst, parameterData, reverseTransform=true)

    peTabOpt = PeTabOpt(evalF, 
                        evalFZygote,
                        evalGradF, 
                        evalGradFZygote,
                        evalGradFAdjoint,
                        evalGradFForwardEq,
                        evalHess,
                        evalHessApprox, 
                        evalHessGaussNewton,
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
