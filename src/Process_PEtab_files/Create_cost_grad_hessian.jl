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
include(joinpath(pwd(), "src", "Process_PEtab_files", "Process_callbacks.jl"))
include(joinpath(pwd(), "src", "Process_PEtab_files", "Observables", "Common.jl"))
include(joinpath(pwd(), "src", "Process_PEtab_files", "Observables", "Create_h_sigma_derivatives.jl"))
include(joinpath(pwd(), "src", "Process_PEtab_files", "Observables", "Create_u0_h_sigma.jl"))
include(joinpath(pwd(), "src", "Process_PEtab_files", "Process_PeTab_files.jl"))
include(joinpath(pwd(), "src", "Common.jl"))


"""
    setUpPeTabModel(modelName::String, dirModel::String)::PEtabModel

    Given a model directory (dirModel) containing the PeTab files and a 
    xml-file on format modelName.xml will return a PEtabModel struct holding 
    paths to PeTab files, ode-system in ModellingToolkit format, functions for 
    evaluating yMod, u0 and standard deviations, and a parameter and state maps 
    for how parameters and states are mapped in the ModellingToolkit ODE system
    along with state and parameter names.

    dirModel must contain a SBML file named modelName.xml, and files starting with 
    measurementData, experimentalCondition, parameter, and observables (tsv-files).
    The latter files must be unique (e.g only one file starting with measurementData)

    TODO : Example  
"""
function readPEtabModel(pathYAML::String; 
                        forceBuildJuliaFiles::Bool=false, 
                        verbose::Bool=true, 
                        ifElseToEvent::Bool=true, 
                        jlFile=false)::PEtabModel

    pathSBML, pathParameters, pathConditions, pathObservables, pathMeasurements, dirJulia, dirModel, modelName = readPEtabYamlFile(pathYAML)                    
    
    if jlFile == false
        
        pathModelJlFile = joinpath(dirJulia, modelName * ".jl")
        
        if !isfile(pathModelJlFile) && forceBuildJuliaFiles == false
            verbose == true && @printf("Julia model file does not exist, will build it\n")
            modelDict = XmlToModellingToolkit(pathSBML, pathModelJlFile, modelName, ifElseToEvent=ifElseToEvent)

        elseif isfile(pathModelJlFile) && forceBuildJuliaFiles == false
            verbose == true && @printf("Julia model file exists at %s - will not rebuilt\n", pathModelJlFile)
            
        elseif forceBuildJuliaFiles == true
            verbose == true && @printf("By user option will rebuild Julia model file\n")
            isfile(pathModelJlFile) == true && rm(pathModelJlFile)
            modelDict = XmlToModellingToolkit(pathSBML, pathModelJlFile, modelName, ifElseToEvent=ifElseToEvent)
        end

    else
        modelDict, modelFileJl = JLToModellingToolkit(modelName, dirModel, ifElseToEvent=ifElseToEvent)
    end

    # Load model ODE-system
    include(pathModelJlFile)
    expr = Expr(:call, Symbol("getODEModel_" * modelName))
    _odeSystem, stateMap, parameterMap = eval(expr)
    odeSystem = structural_simplify(_odeSystem)
    # TODO : Make these to strings here to save conversions
    parameterNames = parameters(odeSystem)
    stateNames = states(odeSystem)

    # Build functions for observables, sd and u0 if does not exist and include
    path_u0_h_sigma = joinpath(dirJulia, modelName * "_h_sd_u0.jl")
    path_D_h_sd = joinpath(dirJulia, modelName * "_D_h_sd.jl")
    if !isfile(path_u0_h_sigma) || !isfile(path_D_h_sd) || forceBuildJuliaFiles == true
        verbose && forceBuildJuliaFiles == false && @printf("File for h, u0 and σ does not exist will build it\n")
        verbose && forceBuildJuliaFiles == true && @printf("By user option will rebuild h, σ and u0\n")
            
        if !@isdefined(modelDict)
            modelDict = XmlToModellingToolkit(modelFileXml, pathModelJlFile, modelName, writeToFile=false, ifElseToEvent=ifElseToEvent)
        end
        create_σ_h_u0_File(modelName, pathYAML, dirJulia, odeSystem, stateMap, modelDict, verbose=verbose)
        createDerivative_σ_h_File(modelName, pathYAML, dirJulia, odeSystem, modelDict, verbose=verbose)
    else
        verbose == true && @printf("File for h, u0 and σ exists will not rebuild it\n")
    end
    include(path_u0_h_sigma)
    include(path_D_h_sd)    

    pathCallback = joinpath(dirJulia, modelName * "_callbacks.jl")
    if !isfile(pathCallback) || forceBuildJuliaFiles == true
        verbose && forceBuildJuliaFiles == false && @printf("File for callback does not exist will build it\n")
        verbose && forceBuildJuliaFiles == true && @printf("By user option will rebuild callback file\n")
        
        if !@isdefined(modelDict)
            modelDict = XmlToModellingToolkit(modelFileXml, pathModelJlFile, modelName, writeToFile=false, ifElseToEvent=ifElseToEvent)
        end
        createCallbacksForTimeDepedentPiecewise(odeSystem, modelDict, modelName, pathYAML, dirJulia)
    end
    include(pathCallback)
    exprCallback = Expr(:call, Symbol("getCallbacks_" * modelName))
    cbSet::CallbackSet, checkCbActive::Vector{Function} = eval(exprCallback)    

    petabModel = PEtabModel(modelName,
                            compute_h,
                            compute_u0!,
                            compute_u0,
                            compute_σ,
                            compute_∂h∂u!,
                            compute_∂σ∂σu!,
                            compute_∂h∂p!,
                            compute_∂σ∂σp!,
                            computeTstops,
                            odeSystem,
                            parameterMap,
                            stateMap,
                            parameterNames, 
                            stateNames,
                            dirModel,
                            dirJulia,
                            pathMeasurements,
                            pathConditions,
                            pathObservables, 
                            pathParameters, 
                            cbSet, 
                            checkCbActive)

    return petabModel
end


function readPEtabYamlFile(pathYAML::AbstractString)

    if !isfile(pathYAML)
        throw(PEtabFileError("Model YAML file does not exist in the model directory"))
    end
    fileYAML = YAML.load_file(pathYAML)

    dirModel = dirname(pathYAML)

    pathSBML = joinpath(dirModel, fileYAML["problems"][1]["sbml_files"][1])
    if !isfile(pathSBML)
        throw(PEtabFileError("SBML file does not exist in the model directory"))
    end

    pathMeasurements = joinpath(dirModel, fileYAML["problems"][1]["measurement_files"][1])
    if !isfile(pathMeasurements)
        throw(PEtabFileError("Measurements file does not exist in the model directory"))
    end

    pathObservables = joinpath(dirModel, fileYAML["problems"][1]["observable_files"][1])
    if !isfile(pathObservables)
        throw(PEtabFileError("Observables file does not exist in the models directory"))
    end

    pathConditions = joinpath(dirModel, fileYAML["problems"][1]["condition_files"][1])
    if !isfile(pathConditions)
        throw(PEtabFileError("Conditions file does not exist in the models directory"))
    end

    pathParameters = joinpath(dirModel, fileYAML["parameter_file"])
    if !isfile(pathParameters)
        throw(PEtabFileError("Parameter file does not exist in the models directory"))
    end

    # Extract YAML directory and use directory name as model name and build directory for Julia files
    dirJulia = joinpath(dirModel, "Julia_model_files") 
    modelName = splitdir(dirModel)[end]
    if !isdir(dirJulia)
        mkdir(dirJulia)
    end

    return pathSBML, pathParameters, pathConditions, pathObservables, pathMeasurements, dirJulia, dirModel, modelName
end


"""
    setUpCostGradHess(petabModel::PEtabModel, solver, tol::Float64)

    For a PeTab-model set up functions for computing i) the likelihood, ii) likelhood gradient, 
    and iii) likelhood Hessian block approximation. The functions are stored in PeTabOpt-struct 
    that can be used as input to the optimizers. 

    Currently the gradient for dynamic parameters (part of ODE-system) is computed via ForwardDiff, 
    and ReverseDiff is used for observable and sd parameters. The hessian approximation assumes the 
    interaction betweeen dynamic and (observable, sd) parameters is zero.
"""
function setUpPEtabODEProblem(petabModel::PEtabModel, 
                              odeSolver::SciMLAlgorithm;
                              solverAbsTol::Float64=1e-8,
                              solverRelTol::Float64=1e-8,
                              solverSSRelTol::Float64=1e-6,
                              solverSSAbsTol::Float64=1e-8,
                              sparseJacobian::Bool=false,
                              numberOfprocesses::Signed=1,
                              sensealgZygote=ForwardDiffSensitivity(),
                              odeSolverForwardEquations::SciMLAlgorithm=Rodas5(autodiff=false),
                              sensealgForwardEquations::Union{Symbol, SciMLSensitivity.AbstractForwardSensitivityAlgorithm}=ForwardSensitivity(),                           
                              odeSolverAdjoint::SciMLAlgorithm=KenCarp4(), 
                              solverAdjointAbsTol::Float64=1e-8, 
                              solverAdjointRelTol::Float64=1e-8,
                              sensealgAdjoint::SciMLSensitivity.AbstractAdjointSensitivityAlgorithm=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(false)), 
                              sensealgAdjointSS::SciMLSensitivity.AbstractAdjointSensitivityAlgorithm=SteadyStateAdjoint())::PEtabODEProblem

    if !(typeof(sensealgAdjointSS) <: SteadyStateAdjoint)
        println("If you are using adjoint sensitivity analysis for a model with PreEq-criteria the most the most efficient adjSensealgSS is usually SteadyStateAdjoint. The algorithm you have provided, ", sensealgAdjointSS, "might not work (as there are some bugs here). In case it does not work, and SteadyStateAdjoint fails (because a dependancy on time or a singular Jacobian) a good choice might be QuadratureAdjoint(autodiff=false, autojacvec=false)")
    end

    experimentalConditions, measurementsData, parametersData, observablesData = readPEtabFiles(petabModel)
    parameterInfo = processParameters(parametersData) 
    measurementInfo = processMeasurements(measurementsData, observablesData) 
    simulationInfo = processSimulationInfo(petabModel, measurementInfo, sensealg=sensealgAdjoint, absTolSS=solverSSAbsTol, relTolSS=solverSSRelTol)
    θ_indices = computeIndicesθ(parameterInfo, measurementInfo, petabModel.odeSystem, experimentalConditions)
    
    # Set up potential prior for the parameters to estimate 
    priorInfo = processPriors(θ_indices, parametersData)

    # Set model parameter values to those in the PeTab parameter to ensure correct value for constant parameters  
    setParamToFileValues!(petabModel.parameterMap, petabModel.stateMap, parameterInfo)

    # The time-span 5e3 is overwritten when performing forward simulations. As we solve an expanded system with the forward 
    # equations, we need a seperate problem for it 
    _odeProblem = ODEProblem(petabModel.odeSystem, petabModel.stateMap, (0.0, 5e3), petabModel.parameterMap, jac=true, sparse=sparseJacobian)
    odeProblem = remake(_odeProblem, p = convert.(Float64, _odeProblem.p), u0 = convert.(Float64, _odeProblem.u0))
    odeProblemForwardEquations = getODEProblemForwardEquations(odeProblem, sensealgForwardEquations)

    # The cost (likelihood) can either be computed in the standard way or the Zygote way. The second consumes more 
    # memory as in-place mutations are not compatible with Zygote 
    expIdSolve = [:all]
    computeCost = setUpCost(:Standard, odeProblem, odeSolver, solverAbsTol, solverRelTol, petabModel, 
                            simulationInfo, θ_indices, measurementInfo, parameterInfo, priorInfo, expIdSolve)
    computeCostZygote = setUpCost(:Zygote, odeProblem, odeSolver, solverAbsTol, solverRelTol, 
                                  petabModel, simulationInfo, θ_indices, measurementInfo, parameterInfo, 
                                  priorInfo, expIdSolve, sensealgZygote=sensealgZygote)          
                                         
    # The gradient can either be computed via autodiff, forward sensitivity equations, adjoint sensitivity equations 
    # and Zygote 
    computeGradientAutoDiff = setUpGradient(:AutoDiff, odeProblem, odeSolver, solverAbsTol, solverRelTol, petabModel, 
                                            simulationInfo, θ_indices, measurementInfo, parameterInfo, priorInfo, expIdSolve)
    computeGradientForwardEquations = setUpGradient(:ForwardEquations, odeProblemForwardEquations, odeSolverForwardEquations, solverAbsTol, 
                                                    solverRelTol, petabModel, simulationInfo, θ_indices, measurementInfo, 
                                                    parameterInfo, priorInfo, expIdSolve, sensealg=sensealgForwardEquations)                                                   
    computeGradientAdjoint = setUpGradient(:Adjoint, odeProblem, odeSolverAdjoint, solverAdjointAbsTol, solverAdjointRelTol, 
                                           petabModel, simulationInfo, θ_indices, measurementInfo, parameterInfo, priorInfo, 
                                           expIdSolve, sensealg=sensealgAdjoint, sensealgSS=sensealgAdjointSS)   
    computeGradientZygote = setUpGradient(:Zygote, odeProblem, odeSolver, solverAbsTol, solverRelTol, petabModel, 
                                          simulationInfo, θ_indices, measurementInfo, parameterInfo, priorInfo, 
                                          expIdSolve, sensealg=sensealgZygote)                 
            
    # The Hessian can either be computed via automatic differentation, or approximated via a block approximation or the 
    # Gauss Newton method                                       
    computeHessian = setUpHessian(:AutoDiff, odeProblem, odeSolver, solverAbsTol, solverRelTol, petabModel, simulationInfo,
                                  θ_indices, measurementInfo, parameterInfo, priorInfo, expIdSolve)
    computeHessianBlock = setUpHessian(:BlockAutoDiff, odeProblem, odeSolver, solverAbsTol, solverRelTol, petabModel, simulationInfo,
                                        θ_indices, measurementInfo, parameterInfo, priorInfo, expIdSolve)                                  
    computeHessianGN = setUpHessian(:GaussNewton, odeProblem, odeSolver, solverAbsTol, solverRelTol, petabModel, simulationInfo,
                                    θ_indices, measurementInfo, parameterInfo, priorInfo, expIdSolve)                                                                          

    # Extract nominal parameter vector and parameter bounds. If needed transform parameters
    θ_estNames = θ_indices.θ_estNames               
    lowerBounds = [parameterInfo.lowerBound[findfirst(x -> x == θ_estNames[i], parameterInfo.parameterId)] for i in eachindex(θ_estNames)] 
    upperBounds = [parameterInfo.upperBound[findfirst(x -> x == θ_estNames[i], parameterInfo.parameterId)] for i in eachindex(θ_estNames)] 
    θ_nominal = [parameterInfo.nominalValue[findfirst(x -> x == θ_estNames[i], parameterInfo.parameterId)] for i in eachindex(θ_estNames)]
    transformθ!(lowerBounds, θ_estNames, θ_indices, reverseTransform=true)
    transformθ!(upperBounds, θ_estNames, θ_indices, reverseTransform=true)
    θ_nominalT = transformθ(θ_nominal, θ_estNames, θ_indices, reverseTransform=true)

    petabProblem = PEtabODEProblem(computeCost, 
                                   computeCostZygote,
                                   computeGradientAutoDiff,
                                   computeGradientZygote,
                                   computeGradientAdjoint, 
                                   computeGradientForwardEquations,
                                   computeHessian,
                                   computeHessianBlock, 
                                   computeHessianGN,
                                   Int64(length(θ_estNames)), 
                                   θ_estNames, 
                                   θ_nominal, 
                                   θ_nominalT, 
                                   lowerBounds, 
                                   upperBounds, 
                                   joinpath(petabModel.dirJulia, "Cube" * petabModel.modelName * ".csv"),
                                   petabModel)
    return petabProblem
end


function setUpCost(whichMethod::Symbol, 
                   odeProblem::ODEProblem,  
                   odeSolver::SciMLAlgorithm, 
                   solverAbsTol::Float64, 
                   solverRelTol::Float64,
                   petabModel::PEtabModel,
                   simulationInfo::SimulationInfo,
                   θ_indices::ParameterIndices,
                   measurementInfo::MeasurementsInfo,
                   parameterInfo::ParametersInfo,
                   priorInfo::PriorInfo,
                   expIDSolveArg::Vector{Symbol}; 
                   sensealgZygote=ForwardDiffSensitivity())

    # Functions needed for mapping θ_est to the ODE problem, and then for solving said ODE-system                           
    if whichMethod == :Standard
        _changeODEProblemParameters! = (pODEProblem, u0, θ_est) -> changeODEProblemParameters!(pODEProblem, u0, θ_est, θ_indices, petabModel)                          
        changeExperimentalCondition! = (pODEProblem, u0, conditionId, θ_dynamic) -> _changeExperimentalCondition!(pODEProblem, u0, conditionId, θ_dynamic, petabModel, θ_indices)
        _solveODEAllExperimentalConditions! = (odeSolutions, odeProblem, θ_dynamic, _expIDSolve) -> solveODEAllExperimentalConditions!(odeSolutions, odeProblem, θ_dynamic, changeExperimentalCondition!, simulationInfo, odeSolver, solverAbsTol, solverRelTol, petabModel.computeTStops, onlySaveAtObservedTimes=true, expIDSolve=_expIDSolve)
        __computeCost = (θ_est) -> computeCost(θ_est, 
                                            odeProblem, 
                                            petabModel, 
                                            simulationInfo, 
                                            θ_indices, 
                                            measurementInfo, 
                                            parameterInfo, 
                                            _changeODEProblemParameters!, 
                                            _solveODEAllExperimentalConditions!, 
                                            priorInfo, 
                                            expIDSolve=expIDSolveArg)
                                            
    elseif whichMethod == :Zygote
        changeExperimentalCondition = (pODEProblem, u0, conditionId, θ_dynamic) -> _changeExperimentalCondition(pODEProblem, u0, conditionId, θ_dynamic, petabModel, θ_indices)
        _changeODEProblemParameters = (pODEProblem, θ_est) -> changeODEProblemParameters(pODEProblem, θ_est, θ_indices, petabModel)
        solveODEExperimentalCondition = (odeProblem, conditionId, θ_dynamic, tMax) -> solveOdeModelAtExperimentalCondZygote(odeProblem, conditionId, θ_dynamic, tMax, changeExperimentalCondition, measurementInfo, simulationInfo, odeSolver, solverAbsTol, solverRelTol, sensealgZygote, petabModel.computeTStops)
        __computeCost = (θ_est) -> computeCostZygote(θ_est, 
                                                   odeProblem, 
                                                   petabModel, 
                                                   simulationInfo, 
                                                   θ_indices, 
                                                   measurementInfo, 
                                                   parameterInfo, 
                                                   _changeODEProblemParameters, 
                                                   solveODEExperimentalCondition, 
                                                   priorInfo)

    end

    return __computeCost
end


function setUpGradient(whichMethod::Symbol, 
                       odeProblem::ODEProblem,  
                       odeSolver::SciMLAlgorithm, 
                       solverAbsTol::Float64, 
                       solverRelTol::Float64,
                       petabModel::PEtabModel,
                       simulationInfo::SimulationInfo,
                       θ_indices::ParameterIndices,
                       measurementInfo::MeasurementsInfo,
                       parameterInfo::ParametersInfo,
                       priorInfo::PriorInfo, 
                       expIDSolve::Vector{Symbol}; 
                       sensealg=nothing, 
                       sensealgSS=nothing)

    # Functions needed for mapping θ_est to the ODE problem, and then for solving said ODE-system                           
    if whichMethod == :AutoDiff
        _changeODEProblemParameters! = (pODEProblem, u0, θ_est) -> changeODEProblemParameters!(pODEProblem, u0, θ_est, θ_indices, petabModel)                          
        changeExperimentalCondition! = (pODEProblem, u0, conditionId, θ_dynamic) -> _changeExperimentalCondition!(pODEProblem, u0, conditionId, θ_dynamic, petabModel, θ_indices)
        _solveODEAllExperimentalConditions! = (odeSolutions, odeProblem, θ_dynamic, _expIDSolve) -> solveODEAllExperimentalConditions!(odeSolutions, odeProblem, θ_dynamic, changeExperimentalCondition!, simulationInfo, odeSolver, solverAbsTol, solverRelTol, petabModel.computeTStops, onlySaveAtObservedTimes=true, expIDSolve=_expIDSolve)
        _computeGradient! = (gradient, θ_est) -> computeGradientAutoDiff!(gradient, 
                                                                         θ_est, 
                                                                         odeProblem, 
                                                                         petabModel, 
                                                                         simulationInfo, 
                                                                         θ_indices, 
                                                                         measurementInfo, 
                                                                         parameterInfo, 
                                                                         _changeODEProblemParameters!, 
                                                                         _solveODEAllExperimentalConditions!, 
                                                                         priorInfo, 
                                                                         expIDSolve=expIDSolve)    
                                            
    elseif whichMethod == :ForwardEquations
        _changeODEProblemParameters! = (pODEProblem, u0, θ_est) -> changeODEProblemParameters!(pODEProblem, u0, θ_est, θ_indices, petabModel)
        if sensealg == :AutoDiffForward            
            changeExperimentalCondition! = (pODEProblem, u0, conditionId, θ_dynamic) -> _changeExperimentalCondition!(pODEProblem, u0, conditionId, θ_dynamic, petabModel, θ_indices)
            _solveODEAllExperimentalConditions! = (odeSolutions, S, odeProblem, θ_dynamic, _expIDSolve) -> solveODEAllExperimentalConditions!(odeSolutions, S, odeProblem, θ_dynamic, changeExperimentalCondition!, _changeODEProblemParameters!, simulationInfo, odeSolver, solverAbsTol, solverRelTol, petabModel.computeTStops, onlySaveAtObservedTimes=true, expIDSolve=_expIDSolve)                                           
        
        else
            changeExperimentalCondition! = (pODEProblem, u0, conditionId, θ_dynamic) -> _changeExperimentalCondition!(pODEProblem, u0, conditionId, θ_dynamic, petabModel, θ_indices, computeForwardSensitivites=true)
            _solveODEAllExperimentalConditions! = (odeSolutions, odeProblem, θ_dynamic, _expIDSolve) -> solveODEAllExperimentalConditions!(odeSolutions, odeProblem, θ_dynamic, changeExperimentalCondition!, simulationInfo, odeSolver, solverAbsTol, solverRelTol, petabModel.computeTStops, onlySaveAtObservedTimes=true, expIDSolve=_expIDSolve)
        end
        _computeGradient! = (gradient, θ_est) -> computeGradientForwardEquations!(gradient, 
                                                                                 θ_est, 
                                                                                 petabModel, 
                                                                                 odeProblem, 
                                                                                 sensealg, 
                                                                                 simulationInfo, 
                                                                                 θ_indices, 
                                                                                 measurementInfo, 
                                                                                 parameterInfo, 
                                                                                 _changeODEProblemParameters!, 
                                                                                 _solveODEAllExperimentalConditions!, 
                                                                                 priorInfo, 
                                                                                 expIDSolve=expIDSolve) 

    elseif whichMethod == :Adjoint                                                                                 
        _changeODEProblemParameters! = (pODEProblem, u0, θ_est) -> changeODEProblemParameters!(pODEProblem, u0, θ_est, θ_indices, petabModel)                          
        changeExperimentalCondition! = (pODEProblem, u0, conditionId, θ_dynamic) -> _changeExperimentalCondition!(pODEProblem, u0, conditionId, θ_dynamic, petabModel, θ_indices)
        _solveODEAllExperimentalConditions! = (odeSolutions, odeProblem, θ_dynamic, _expIDSolve) -> solveODEAllExperimentalConditions!(odeSolutions, odeProblem, θ_dynamic, changeExperimentalCondition!, simulationInfo, odeSolver, solverAbsTol, solverRelTol, petabModel.computeTStops, denseSolution=true, expIDSolve=_expIDSolve, trackCallback=true)
        _computeGradient! = (gradient, θ_est) -> computeGradientAdjointEquations!(gradient, 
                                                                                 θ_est, 
                                                                                 odeSolver, 
                                                                                 sensealg, 
                                                                                 sensealgSS, 
                                                                                 solverAbsTol, 
                                                                                 solverRelTol, 
                                                                                 odeProblem, 
                                                                                 petabModel, 
                                                                                 simulationInfo, 
                                                                                 θ_indices, 
                                                                                 measurementInfo, 
                                                                                 parameterInfo, 
                                                                                 _changeODEProblemParameters!, 
                                                                                 _solveODEAllExperimentalConditions!, 
                                                                                 priorInfo, 
                                                                                 expIDSolve=expIDSolve) 

    elseif whichMethod == :Zygote

        changeExperimentalCondition = (pODEProblem, u0, conditionId, θ_dynamic) -> _changeExperimentalCondition(pODEProblem, u0, conditionId, θ_dynamic, petabModel, θ_indices)
        _changeODEProblemParameters = (pODEProblem, θ_est) -> changeODEProblemParameters(pODEProblem, θ_est, θ_indices, petabModel)
        solveODEExperimentalCondition = (odeProblem, conditionId, θ_dynamic, tMax) -> solveOdeModelAtExperimentalCondZygote(odeProblem, conditionId, θ_dynamic, tMax, changeExperimentalCondition, measurementInfo, simulationInfo, odeSolver, solverAbsTol, solverRelTol, sensealg, petabModel.computeTStops)
        _computeGradient! = (gradient, θ_est) -> computeGradientZygote(gradient, 
                                                                      θ_est, 
                                                                      odeProblem, 
                                                                      petabModel, 
                                                                      simulationInfo, 
                                                                      θ_indices, 
                                                                      measurementInfo, 
                                                                      parameterInfo, 
                                                                      _changeODEProblemParameters, 
                                                                      solveODEExperimentalCondition, 
                                                                      priorInfo)
    end

    return _computeGradient!
end


function setUpHessian(whichMethod::Symbol, 
                      odeProblem::ODEProblem,  
                      odeSolver::SciMLAlgorithm, 
                      solverAbsTol::Float64, 
                      solverRelTol::Float64,
                      petabModel::PEtabModel,
                      simulationInfo::SimulationInfo,
                      θ_indices::ParameterIndices,
                      measurementInfo::MeasurementsInfo,
                      parameterInfo::ParametersInfo,
                      priorInfo::PriorInfo, 
                      expIDSolve::Vector{Symbol})

    # Functions needed for mapping θ_est to the ODE problem, and then for solving said ODE-system                           
    if whichMethod == :AutoDiff || whichMethod == :BlockAutoDiff
        _changeODEProblemParameters! = (pODEProblem, u0, θ_est) -> changeODEProblemParameters!(pODEProblem, u0, θ_est, θ_indices, petabModel)                          
        changeExperimentalCondition! = (pODEProblem, u0, conditionId, θ_dynamic) -> _changeExperimentalCondition!(pODEProblem, u0, conditionId, θ_dynamic, petabModel, θ_indices)
        _solveODEAllExperimentalConditions! = (odeSolutions, odeProblem, θ_dynamic, _expIDSolve) -> solveODEAllExperimentalConditions!(odeSolutions, odeProblem, θ_dynamic, changeExperimentalCondition!, simulationInfo, odeSolver, solverAbsTol, solverRelTol, petabModel.computeTStops, onlySaveAtObservedTimes=true, expIDSolve=_expIDSolve)

        if whichMethod == :AutoDiff
            _computeHessian = (hessian, θ_est) -> computeHessian!(hessian, 
                                                                θ_est, 
                                                                odeProblem, 
                                                                petabModel, 
                                                                simulationInfo, 
                                                                θ_indices, 
                                                                measurementInfo, 
                                                                parameterInfo, 
                                                                _changeODEProblemParameters!, 
                                                                _solveODEAllExperimentalConditions!, 
                                                                priorInfo, 
                                                                expIDSolve=expIDSolve)
        else
            _computeHessian = (hessian, θ_est) -> computeHessianBlockApproximation!(hessian, 
                                                                                  θ_est, 
                                                                                  odeProblem, 
                                                                                  petabModel, 
                                                                                  simulationInfo, 
                                                                                  θ_indices, 
                                                                                  measurementInfo, 
                                                                                  parameterInfo, 
                                                                                  _changeODEProblemParameters!, 
                                                                                  _solveODEAllExperimentalConditions!, 
                                                                                  priorInfo, 
                                                                                  expIDSolve=expIDSolve)            
        end
                                            
    elseif whichMethod == :GaussNewton
        _changeODEProblemParameters! = (pODEProblem, u0, θ_est) -> changeODEProblemParameters!(pODEProblem, u0, θ_est, θ_indices, petabModel)
        changeExperimentalCondition! = (pODEProblem, u0, conditionId, θ_dynamic) -> _changeExperimentalCondition!(pODEProblem, u0, conditionId, θ_dynamic, petabModel, θ_indices)
        _solveODEAllExperimentalConditions! = (odeSolutions, S, odeProblem, θ_dynamic, _expIDSolve) -> solveODEAllExperimentalConditions!(odeSolutions, S, odeProblem, θ_dynamic, changeExperimentalCondition!, _changeODEProblemParameters!, simulationInfo, odeSolver, solverAbsTol, solverRelTol, petabModel.computeTStops, onlySaveAtObservedTimes=true, expIDSolve=_expIDSolve)                                           
        _computeHessian = (hessian, θ_est) -> computeGaussNewtonHessianApproximation!(hessian, 
                                                                                    θ_est, 
                                                                                    odeProblem, 
                                                                                    petabModel, 
                                                                                    simulationInfo, 
                                                                                    θ_indices, 
                                                                                    measurementInfo, 
                                                                                    parameterInfo, 
                                                                                    _changeODEProblemParameters!, 
                                                                                    _solveODEAllExperimentalConditions!, 
                                                                                    priorInfo, 
                                                                                    expIDSolve=expIDSolve) 

    end

    return _computeHessian
end


function getODEProblemForwardEquations(odeProblem::ODEProblem, 
                                       sensealgForwardEquations::SciMLSensitivity.AbstractForwardSensitivityAlgorithm)::ODEProblem
    return ODEForwardSensitivityProblem(odeProblem.f, odeProblem.u0, odeProblem.tspan, odeProblem.p, sensealg=sensealgForwardEquations)
end
function getODEProblemForwardEquations(odeProblem::ODEProblem, 
                                       sensealgForwardEquations::Symbol)::ODEProblem
    return deepcopy(odeProblem)
end
    