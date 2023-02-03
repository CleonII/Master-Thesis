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
                              sensealgAdjointSS::SciMLSensitivity.AbstractAdjointSensitivityAlgorithm=SteadyStateAdjoint(), 
                              chunkSize::Union{Nothing, Int64}=nothing, 
                              terminateSSMethod::Symbol=:Norm)::PEtabODEProblem

    if !(typeof(sensealgAdjointSS) <: SteadyStateAdjoint)
        println("If you are using adjoint sensitivity analysis for a model with PreEq-criteria the most the most efficient adjSensealgSS is usually SteadyStateAdjoint. The algorithm you have provided, ", sensealgAdjointSS, "might not work (as there are some bugs here). In case it does not work, and SteadyStateAdjoint fails (because a dependancy on time or a singular Jacobian) a good choice might be QuadratureAdjoint(autodiff=false, autojacvec=false)")
    end

    experimentalConditions, measurementsData, parametersData, observablesData = readPEtabFiles(petabModel)
    parameterInfo = processParameters(parametersData) 
    measurementInfo = processMeasurements(measurementsData, observablesData) 
    simulationInfo = processSimulationInfo(petabModel, measurementInfo, sensealg=sensealgAdjoint, absTolSS=solverSSAbsTol, relTolSS=solverSSRelTol, terminateSSMethod=terminateSSMethod)
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
                                            simulationInfo, θ_indices, measurementInfo, parameterInfo, priorInfo, expIdSolve, 
                                            chunkSize=chunkSize)
    computeGradientForwardEquations = setUpGradient(:ForwardEquations, odeProblemForwardEquations, odeSolverForwardEquations, solverAbsTol, 
                                                    solverRelTol, petabModel, simulationInfo, θ_indices, measurementInfo, 
                                                    parameterInfo, priorInfo, expIdSolve, sensealg=sensealgForwardEquations, chunkSize=chunkSize)                                                   
    computeGradientAdjoint = setUpGradient(:Adjoint, odeProblem, odeSolverAdjoint, solverAdjointAbsTol, solverAdjointRelTol, 
                                           petabModel, simulationInfo, θ_indices, measurementInfo, parameterInfo, priorInfo, 
                                           expIdSolve, sensealg=sensealgAdjoint, sensealgSS=sensealgAdjointSS)   
    computeGradientZygote = setUpGradient(:Zygote, odeProblem, odeSolver, solverAbsTol, solverRelTol, petabModel, 
                                          simulationInfo, θ_indices, measurementInfo, parameterInfo, priorInfo, 
                                          expIdSolve, sensealg=sensealgZygote)                 
            
    # The Hessian can either be computed via automatic differentation, or approximated via a block approximation or the 
    # Gauss Newton method                                       
    computeHessian = setUpHessian(:AutoDiff, odeProblem, odeSolver, solverAbsTol, solverRelTol, petabModel, simulationInfo,
                                  θ_indices, measurementInfo, parameterInfo, priorInfo, chunkSize, expIdSolve)
    computeHessianBlock = setUpHessian(:BlockAutoDiff, odeProblem, odeSolver, solverAbsTol, solverRelTol, petabModel, simulationInfo,
                                        θ_indices, measurementInfo, parameterInfo, priorInfo, chunkSize, expIdSolve)                                  
    computeHessianGN = setUpHessian(:GaussNewton, odeProblem, odeSolver, solverAbsTol, solverRelTol, petabModel, simulationInfo,
                                    θ_indices, measurementInfo, parameterInfo, priorInfo, chunkSize, expIdSolve)                                                                          

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
                       chunkSize::Union{Nothing, Int64}=nothing,
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
                                                                         chunkSize,
                                                                         expIDSolve=expIDSolve)    
                                            
    elseif whichMethod == :ForwardEquations
        _changeODEProblemParameters! = (pODEProblem, u0, θ_est) -> changeODEProblemParameters!(pODEProblem, u0, θ_est, θ_indices, petabModel)
        if sensealg == :AutoDiffForward  
            nTimePointsSaveAt = sum(length(simulationInfo.timeObserved[experimentalConditionId]) for experimentalConditionId in simulationInfo.experimentalConditionId)                                                    
            nModelStates = length(odeProblem.u0)                                                                                                      
            odeSolutionValues = zeros(Float64, nModelStates, nTimePointsSaveAt)          
            changeExperimentalCondition! = (pODEProblem, u0, conditionId, θ_dynamic) -> _changeExperimentalCondition!(pODEProblem, u0, conditionId, θ_dynamic, petabModel, θ_indices)
            _solveODEAllExperimentalConditions! = (odeSolutions, S, odeProblem, θ_dynamic, _expIDSolve) -> solveODEAllExperimentalConditions!(odeSolutions, S, odeProblem, θ_dynamic, changeExperimentalCondition!, _changeODEProblemParameters!, simulationInfo, odeSolver, solverAbsTol, solverRelTol, petabModel.computeTStops, odeSolutionValues, onlySaveAtObservedTimes=true, expIDSolve=_expIDSolve, chunkSize=chunkSize)                                           
        
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
                      chunkSize::Union{Nothing, Int64},
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
                                                                chunkSize,
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
                                                                                  chunkSize,
                                                                                  expIDSolve=expIDSolve)            
        end
                                            
    elseif whichMethod == :GaussNewton
        _changeODEProblemParameters! = (pODEProblem, u0, θ_est) -> changeODEProblemParameters!(pODEProblem, u0, θ_est, θ_indices, petabModel)
        changeExperimentalCondition! = (pODEProblem, u0, conditionId, θ_dynamic) -> _changeExperimentalCondition!(pODEProblem, u0, conditionId, θ_dynamic, petabModel, θ_indices)
        nTimePointsSaveAt = sum(length(simulationInfo.timeObserved[experimentalConditionId]) for experimentalConditionId in simulationInfo.experimentalConditionId)                                                    
        nModelStates = length(odeProblem.u0)                                                                                                      
        odeSolutionValues = zeros(Float64, nModelStates, nTimePointsSaveAt)          
        _solveODEAllExperimentalConditions! = (odeSolutions, S, odeProblem, θ_dynamic, _expIDSolve) -> solveODEAllExperimentalConditions!(odeSolutions, S, odeProblem, θ_dynamic, changeExperimentalCondition!, _changeODEProblemParameters!, simulationInfo, odeSolver, solverAbsTol, solverRelTol, petabModel.computeTStops, odeSolutionValues, onlySaveAtObservedTimes=true, expIDSolve=_expIDSolve, chunkSize=chunkSize)                                           
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