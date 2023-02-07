#=
    The top-level functions for computing the gradient via i) exactly via forward-mode autodiff, ii) forward sensitivty 
    eqations, iii) adjoint sensitivity analysis and iv) Zygote interface. 

    Due to it slow speed Zygote does not have full support for all models, e.g, models with priors and pre-eq criteria.
=#


# Compute the gradient via forward mode automatic differentitation 
function computeGradientAutoDiff!(gradient::Vector{Float64}, 
                                  θ_est::Vector{Float64}, 
                                  odeProblem::ODEProblem,  
                                  petabModel::PEtabModel,
                                  simulationInfo::SimulationInfo,
                                  θ_indices::ParameterIndices,
                                  measurementInfo::MeasurementsInfo,
                                  parameterInfo::ParametersInfo, 
                                  changeODEProblemParameters!::Function,
                                  solveOdeModelAllConditions!::Function, 
                                  priorInfo::PriorInfo, 
                                  chunkSize::Union{Nothing, Int64};
                                  expIDSolve::Vector{Symbol} = [:all])     

    θ_dynamic, θ_observable, θ_sd, θ_nonDynamic = splitParameterVector(θ_est, θ_indices) 
    # Compute hessian for parameters which are a part of the ODE-system (dynamic parameters)
    computeCostDynamicθ = (x) -> computeCostSolveODE(x, θ_sd, θ_observable, θ_nonDynamic, odeProblem, petabModel, 
                                                     simulationInfo, θ_indices, measurementInfo, parameterInfo, 
                                                     changeODEProblemParameters!, solveOdeModelAllConditions!, 
                                                     computeGradientDynamicθ=true, expIDSolve=expIDSolve)

    if !isnothing(chunkSize)                                                     
        cfg = ForwardDiff.GradientConfig(computeCostDynamicθ, θ_dynamic, ForwardDiff.Chunk(chunkSize))
    else
        cfg = ForwardDiff.GradientConfig(computeCostDynamicθ, θ_dynamic, ForwardDiff.Chunk(θ_dynamic))
    end

    try 
        @views ForwardDiff.gradient!(gradient[θ_indices.iθ_dynamic], computeCostDynamicθ, θ_dynamic, cfg)
    catch
        gradient .= 1e8
        return
    end


    # Check if we could solve the ODE (first), and if Inf was returned (second)
    if couldSolveODEModel(simulationInfo, expIDSolve) == false
        gradient .= 0.0
        return
    end
    if all(gradient[θ_indices.iθ_dynamic] .== 0.0)
        gradient .= 1e8
        return 
    end

    # Compute hessian for parameters which are not in ODE-system. Important to keep in mind that Sd- and observable 
    # parameters can overlap in θ_est.
    iθ_sd, iθ_observable, iθ_nonDynamic, iθ_notOdeSystem = getIndicesParametersNotInODESystem(θ_indices)
    computeCostNotODESystemθ = (x) -> computeCostNotSolveODE(θ_dynamic, x[iθ_sd], x[iθ_observable], x[iθ_nonDynamic], 
                                                             petabModel, simulationInfo, θ_indices, measurementInfo, 
                                                             parameterInfo, expIDSolve=expIDSolve, 
                                                             computeGradientNotSolveAutoDiff=true)
    @views ReverseDiff.gradient!(gradient[iθ_notOdeSystem], computeCostNotODESystemθ, θ_est[iθ_notOdeSystem])

    # If we have prior contribution its gradient is computed via autodiff for all parameters 
    if priorInfo.hasPriors == true
        computeGradientPrior!(gradient, θ_est, θ_indices, priorInfo)
    end
end


# Compute the gradient via forward sensitivity equations 
function computeGradientForwardEquations!(gradient::Vector{Float64}, 
                                          θ_est::Vector{Float64}, 
                                          petabModel::PEtabModel,
                                          odeProblem::ODEProblem,
                                          sensealg::Union{Symbol, SciMLSensitivity.AbstractForwardSensitivityAlgorithm},
                                          simulationInfo::SimulationInfo,
                                          θ_indices::ParameterIndices,
                                          measurementInfo::MeasurementsInfo,
                                          parameterInfo::ParametersInfo, 
                                          changeODEProblemParameters!::Function,
                                          solveOdeModelAllConditions!::Function, 
                                          priorInfo::PriorInfo;
                                          expIDSolve::Vector{Symbol} = [:all])   
    
    θ_dynamic, θ_observable, θ_sd, θ_nonDynamic = splitParameterVector(θ_est, θ_indices) 

    # Calculate gradient seperately for dynamic and non dynamic parameter. 
    gradientDyanmicθ::Vector{Float64} = zeros(Float64, length(θ_dynamic))
    computeGradientForwardEqDynamicθ!(gradientDyanmicθ, θ_dynamic, θ_sd, θ_observable, θ_nonDynamic, simulationInfo.S, petabModel, 
                                      sensealg, odeProblem, simulationInfo, θ_indices, measurementInfo, parameterInfo, 
                                      changeODEProblemParameters!, solveOdeModelAllConditions!, expIDSolve=expIDSolve)                            
    gradient[θ_indices.iθ_dynamic] .= gradientDyanmicθ

    # Happens when at least one forward pass fails and I set the gradient to 1e8 
    if all(gradientDyanmicθ .== 1e8)
        gradient .= 1e8
        return 
    end

    # Compute gradient for parameters which are not in ODE-system. Important to keep in mind that Sd- and observable 
    # parameters can overlap in θ_est.
    iθ_sd, iθ_observable, iθ_nonDynamic, iθ_notOdeSystem = getIndicesParametersNotInODESystem(θ_indices)
    computeCostNotODESystemθ = (x) -> computeCostNotSolveODE(θ_dynamic, x[iθ_sd], x[iθ_observable], x[iθ_nonDynamic], 
                                                             petabModel, simulationInfo, θ_indices, measurementInfo, 
                                                             parameterInfo, expIDSolve=expIDSolve, 
                                                             computeGradientNotSolveForward=true)
    @views ReverseDiff.gradient!(gradient[iθ_notOdeSystem], computeCostNotODESystemθ, θ_est[iθ_notOdeSystem])

    if priorInfo.hasPriors == true
        computeGradientPrior!(gradient, θ_est, θ_indices, priorInfo)
    end
end


# Compute gradient via adjoint sensitivity analysis 
function computeGradientAdjointEquations!(gradient::Vector{Float64}, 
                                          θ_est::Vector{Float64}, 
                                          adjointODESolver::SciMLAlgorithm, 
                                          sensealg::SciMLSensitivity.AbstractAdjointSensitivityAlgorithm,
                                          sensealgSS::SciMLSensitivity.AbstractAdjointSensitivityAlgorithm,
                                          solverAbsTol::Float64,
                                          solverRelTol::Float64,
                                          odeProblem::ODEProblem,  
                                          petabModel::PEtabModel,
                                          simulationInfo::SimulationInfo,
                                          θ_indices::ParameterIndices,
                                          measurementInfo::MeasurementsInfo,
                                          parameterInfo::ParametersInfo, 
                                          changeODEProblemParameters!::Function,
                                          solveOdeModelAllConditions!::Function, 
                                          priorInfo::PriorInfo;
                                          expIDSolve::Vector{Symbol} = [:all])   
    
    # Split input into observeble and dynamic parameters 
    θ_dynamic, θ_observable, θ_sd, θ_nonDynamic = splitParameterVector(θ_est, θ_indices) 

    # Calculate gradient seperately for dynamic and non dynamic parameter. 
    gradientDyanmicθ::Vector{Float64} = zeros(Float64, length(θ_dynamic))
    computeGradientAdjointDynamicθ(gradientDyanmicθ, θ_dynamic, θ_sd, θ_observable, θ_nonDynamic, odeProblem, adjointODESolver, 
                                   solverAbsTol, solverRelTol, sensealg, petabModel, simulationInfo, θ_indices, measurementInfo, parameterInfo, 
                                   changeODEProblemParameters!, solveOdeModelAllConditions!; expIDSolve=expIDSolve, 
                                   sensealgSS=sensealgSS)
    gradient[θ_indices.iθ_dynamic] .= gradientDyanmicθ

    # Happens when at least one forward pass fails and I set the gradient to 1e8 
    if all(gradientDyanmicθ .== 1e8)
        gradient .= 1e8
        return 
    end

    # Compute gradient for parameters which are not in ODE-system. Important to keep in mind that Sd- and observable 
    # parameters can overlap in θ_est.
    iθ_sd, iθ_observable, iθ_nonDynamic, iθ_notOdeSystem = getIndicesParametersNotInODESystem(θ_indices)
    computeCostNotODESystemθ = (x) -> computeCostNotSolveODE(θ_dynamic, x[iθ_sd], x[iθ_observable], x[iθ_nonDynamic], 
                                                             petabModel, simulationInfo, θ_indices, measurementInfo, 
                                                             parameterInfo, expIDSolve=expIDSolve, 
                                                             computeGradientNotSolveAdjoint=true)
    @views ReverseDiff.gradient!(gradient[iθ_notOdeSystem], computeCostNotODESystemθ, θ_est[iθ_notOdeSystem])

    if priorInfo.hasPriors == true
        computeGradientPrior!(gradient, θ_est, θ_indices, priorInfo)
    end
end


# Compute gradient using Zygote 
function computeGradientZygote(gradient::Vector{Float64}, 
                               θ_est::Vector{Float64}, 
                               odeProblem::ODEProblem,  
                               petabModel::PEtabModel,
                               simulationInfo::SimulationInfo,
                               θ_indices::ParameterIndices,
                               measurementInfo::MeasurementsInfo,
                               parameterInfo::ParametersInfo, 
                               changeODEProblemParameters::Function,
                               solveOdeModelAllConditions::Function, 
                               priorInfo::PriorInfo)
    
    # Split input into observeble and dynamic parameters 
    θ_dynamic, θ_observable, θ_sd, θ_nonDynamic = splitParameterVector(θ_est, θ_indices) 

    # For Zygote the code must be out-of place. Hence a special likelihood funciton is needed.
    computeGradientZygoteDynamicθ = (x) -> _computeCostZygote(x, θ_sd, θ_observable, θ_nonDynamic, odeProblem, 
                                                              petabModel, simulationInfo, θ_indices, measurementInfo, 
                                                              parameterInfo, changeODEProblemParameters, 
                                                              solveOdeModelAllConditions)
    gradient[θ_indices.iθ_dynamic] .= Zygote.gradient(computeGradientZygoteDynamicθ, θ_dynamic)[1]

    # Compute gradient for parameters which are not in ODE-system. Important to keep in mind that Sd- and observable 
    # parameters can overlap in θ_est.
    iθ_sd, iθ_observable, iθ_nonDynamic, iθ_notOdeSystem = getIndicesParametersNotInODESystem(θ_indices)
    computeCostNotODESystemθ = (x) -> computeCostNotSolveODE(θ_dynamic, x[iθ_sd], x[iθ_observable], x[iθ_nonDynamic], 
                                                             petabModel, simulationInfo, θ_indices, measurementInfo, 
                                                             parameterInfo)
    @views ReverseDiff.gradient!(gradient[iθ_notOdeSystem], computeCostNotODESystemθ, θ_est[iθ_notOdeSystem])
end    


# Compute prior contribution to log-likelihood 
function computeGradientPrior!(gradient::AbstractVector, 
                               θ::AbstractVector, 
                               θ_indices::ParameterIndices, 
                               priorInfo::PriorInfo)
            

    _evalPriors = (θ_est) -> begin
                                θ_estT = transformθ(θ_est, θ_indices.θ_estNames, θ_indices)
                                return computePriors(θ_est, θ_estT, θ_indices.θ_estNames, priorInfo)
                            end
    gradient .+= ForwardDiff.gradient(_evalPriors, θ)                                
end