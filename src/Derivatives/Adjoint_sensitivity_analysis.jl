# Compute the adjoint gradient across all experimental conditions 
function computeGradientAdjointDynamicθ(gradient::Vector{Float64},
                                        θ_dynamic::Vector{Float64},
                                        θ_sd::Vector{Float64},
                                        θ_observable::Vector{Float64},
                                        θ_nonDynamic::Vector{Float64},
                                        odeProblem::ODEProblem,
                                        odeSolver::SciMLAlgorithm, 
                                        tolerance::Float64,
                                        sensealg::SciMLSensitivity.AbstractAdjointSensitivityAlgorithm,
                                        peTabModel::PeTabModel,
                                        simulationInfo::SimulationInfo,
                                        θ_indices::ParameterIndices,
                                        measurementInfo ::MeasurementsInfo, 
                                        parameterInfo::ParametersInfo, 
                                        changeODEProblemParameters!::Function,
                                        solveOdeModelAllConditions!::Function;
                                        sensealgSS=SteadyStateAdjoint(),
                                        expIDSolve::Vector{Symbol} = [:all])

    θ_dynamicT = transformθ(θ_dynamic, θ_indices.θ_dynamicNames, θ_indices)
    θ_sdT = transformθ(θ_sd, θ_indices.θ_sdNames, θ_indices)
    θ_observableT = transformθ(θ_observable, θ_indices.θ_observableNames, θ_indices)
    θ_nonDynamicT = transformθ(θ_nonDynamic, θ_indices.θ_nonDynamicNames, θ_indices)

    _odeProblem = remake(odeProblem, p = convert.(eltype(θ_dynamicT), odeProblem.p), u0 = convert.(eltype(θ_dynamicT), odeProblem.u0))
    changeODEProblemParameters!(_odeProblem.p, _odeProblem.u0, θ_dynamicT)
    success = solveOdeModelAllConditions!(simulationInfo.odeSolutionsDerivatives, _odeProblem, θ_dynamicT, expIDSolve)
    if success != true
        gradient .= 1e8
        return
    end

    # In case of PreEq-critera we need to compute the pullback function at tSS to compute the VJP between 
    # λ_t0 and the sensitivites at steady state time
    if simulationInfo.haspreEquilibrationConditionId == true
        evalVJPSSVec = generateVJPSSFunction(simulationInfo, sensealgSS, odeSolver, tolerance, expIDSolve)
    end

    gradient .= 0.0
    # Compute the gradient by looping through all experimental conditions.
    for i in eachindex(simulationInfo.experimentalConditionId)
        experimentalConditionId = simulationInfo.experimentalConditionId[i]
        simulationConditionId = simulationInfo.simulationConditionId[i]
        
        if expIDSolve[1] != :all && experimentalConditionId ∉ experimentalConditionId
            continue
        end

        if (haspreEquilibrationConditionId = simulationInfo.haspreEquilibrationConditionId) == true
            evalVJPSS = evalVJPSSVec[simulationInfo.preEquilibrationConditionId[i]]
        else
            evalVJPSS = identity
        end

        # In case the model is simulated first to a steady state we need to keep track of the post-equlibrium experimental 
        # condition Id to identify parameters specific to an experimental condition.
        sol = simulationInfo.odeSolutionsDerivatives[experimentalConditionId]
        success = computeGradientAdjointExpCond!(gradient, sol, sensealg, tolerance, odeSolver, θ_dynamicT,
                                                 θ_sdT, θ_observableT, θ_nonDynamicT, experimentalConditionId, 
                                                 simulationConditionId, simulationInfo,
                                                 peTabModel, θ_indices, measurementInfo, parameterInfo, evalVJPSS)

        if success == false
            gradient .= 1e8
            return
        end
    end
    return 
end


function generateVJPSSFunction(simulationInfo::SimulationInfo, 
                               sensealgSS::SteadyStateAdjoint, 
                               odeSolver::SciMLAlgorithm,
                               tolerance::Float64,
                               expIDSolve::Vector{Symbol})::NamedTuple

    # Extract all unique Pre-equlibrium conditions. If the code is run in parallell 
    # (expIDSolve != [["all]]) the number of preEq cond. might be smaller than the 
    # total number of preEq cond.
    if expIDSolve[1] == :all
        preEquilibrationConditionId = unique(simulationInfo.preEquilibrationConditionId)
    else
        whichIds  = findall(x -> x ∈ simulationInfo.experimentalConditionId, expIDSolve)
        preEquilibrationConditionId = unique(simulationInfo.preEquilibrationConditionId[whichIds])
    end

    _evalVJPSS = Vector{Function}(undef, length(preEquilibrationConditionId))
    for i in eachindex(preEquilibrationConditionId)
        
        odeProblem = simulationInfo.odePreEqulibriumSolutions[preEquilibrationConditionId[i]].prob
        ssOdeProblem = SteadyStateProblem(odeProblem)
        ySS, _evalVJPSSi = Zygote.pullback((p) ->    (
                                                      solve(ssOdeProblem, 
                                                            DynamicSS(odeSolver, abstol=simulationInfo.absTolSS, reltol=simulationInfo.relTolSS), 
                                                            abstol=tolerance, 
                                                            reltol=tolerance, 
                                                            p=p, 
                                                            sensealg=sensealgSS)[:]), odeProblem.p)
                                                
        _evalVJPSS[i] = _evalVJPSSi
    end

    evalVJPSS = Tuple(f for f in _evalVJPSS)
    return NamedTuple{Tuple(name for name in preEquilibrationConditionId)}(evalVJPSS)
end
function generateVJPSSFunction(simulationInfo::SimulationInfo, 
                               sensealgSS::Union{QuadratureAdjoint, InterpolatingAdjoint}, 
                               odeSolver::SciMLAlgorithm,
                               tolerance::Float64,
                               expIDSolve::Vector{Symbol})::NamedTuple

    # Extract all unique Pre-equlibrium conditions. If the code is run in parallell 
    # (expIDSolve != [["all]]) the number of preEq cond. might be smaller than the 
    # total number of preEq cond.
    if expIDSolve[1] == :all
        preEquilibrationConditionId = unique(simulationInfo.preEquilibrationConditionId)
    else
        whichIds  = findall(x -> x ∈ simulationInfo.experimentalConditionId, expIDSolve)
        preEquilibrationConditionId = unique(simulationInfo.preEquilibrationConditionId[whichIds])
    end

    _evalVJPSS = Vector{Function}(undef, length(preEquilibrationConditionId))
    for i in eachindex(preEquilibrationConditionId)
        
        # As we already have solved for steady state once and know the steady state time we here 
        # build a problem where we simulate exactly to said steady state.
        preEqulibriumOdeSolution = simulationInfo.odePreEqulibriumSolutions[preEquilibrationConditionId[i]]
        odeProblemPullback = remake(preEqulibriumOdeSolution.prob, tspan=(0.0, preEqulibriumOdeSolution.t[end]))
        ySS, _evalVJPSSi = Zygote.pullback((p) -> solve(odeProblemPullback, odeSolver, p=p, abstol=tolerance, reltol=tolerance, sensealg=sensealgSS)[:, end], preEqulibriumOdeSolution.prob.p)
                                                
        _evalVJPSS[i] = _evalVJPSSi
    end

    evalVJPSS = Tuple(f for f in _evalVJPSS)
    return NamedTuple{Tuple(name for name in preEquilibrationConditionId)}(evalVJPSS)
end


# For a given experimental condition compute the gradient using adjoint sensitivity analysis 
# for a funciton on the form G = (h - yObs)^2 / σ^2
# TODO : Important function - improve documentation.
function computeGradientAdjointExpCond!(gradient::Vector{Float64},
                                        sol::ODESolution,
                                        sensealg::SciMLSensitivity.AbstractAdjointSensitivityAlgorithm,
                                        tolerance::Float64,
                                        odeSolver::SciMLAlgorithm, 
                                        θ_dynamic::Vector{Float64},
                                        θ_sd::Vector{Float64}, 
                                        θ_observable::Vector{Float64}, 
                                        θ_nonDynamic::Vector{Float64},
                                        experimentalConditionId::Symbol,
                                        simulationConditionId::Symbol,
                                        simulationInfo::SimulationInfo,
                                        peTabModel::PeTabModel,
                                        θ_indices::ParameterIndices,
                                        measurementInfo::MeasurementsInfo, 
                                        parameterInfo::ParametersInfo, 
                                        evalVJPSS::Function)::Bool
                                        
    # Extract experimetnalCondition specific parameter required to solve the 
    # adjoitn ODE                                        
    iPerTimePoint = simulationInfo.iPerTimePoint[experimentalConditionId]                                        
    timeObserved = simulationInfo.timeObserved[experimentalConditionId]
    callback = simulationInfo.callbacks[experimentalConditionId]

    # Pre allcoate vectors needed for computations 
    ∂h∂u, ∂σ∂u, ∂h∂p, ∂σ∂p = allocateObservableFunctionDerivatives(sol, peTabModel) 
    
    compute∂G∂u = (out, u, p, t, i) -> begin compute∂G∂_(out, u, p, t, i, iPerTimePoint, 
                                                         measurementInfo, parameterInfo, 
                                                         θ_indices, peTabModel, 
                                                         θ_dynamic, θ_sd, θ_observable, θ_nonDynamic, 
                                                         ∂h∂u, ∂σ∂u, compute∂G∂U=true)
                                            end
    compute∂G∂p = (out, u, p, t, i) -> begin compute∂G∂_(out, u, p, t, i, iPerTimePoint, 
                                                         measurementInfo, parameterInfo, 
                                                         θ_indices, peTabModel, 
                                                         θ_dynamic, θ_sd, θ_observable, θ_nonDynamic, 
                                                         ∂h∂p, ∂σ∂p, compute∂G∂U=false)
                                        end                                        
                                                                   
    # The standard allow cases where we only observe data at t0, that is we do not solve the ODE. Here adjoint_sensitivities fails (naturally). In this case we compute the gradient 
    # via ∇G_p = dp + du*J(u(t_0)) where du is the cost function differentiated with respect to the states at time zero, 
    # dp is the cost function  with respect to the parameters at time zero and J is sensititvites at time 
    # zero. Overall, the only workflow that changes below is that we compute du outside of the adjoint interface 
    # and use sol[:] as we no longer can interpolate from the forward solution.
    onlyObsAtZero::Bool = false
    if !(length(timeObserved) == 1 && timeObserved[1] == 0.0)

        # adjoint_sensitivities does not return a retcode. Hence upon integration error only a warning is thrown. To capture 
        # this stderr is redirected to a read only stream. Thus upon warning an error is triggered, upon which we return 
        # false to signify that we could not compute the gradient.
        stderrOld = stderr
        redirect_stderr(open(touch(tempname()), "r"))
        local du, dp
            try
            du, dp = adjoint_sensitivities(sol, 
                                           odeSolver,
                                           dgdp_discrete=nothing,
                                           dgdu_discrete=compute∂G∂u, 
                                           callback=callback,
                                           t=timeObserved, 
                                           sensealg=sensealg, 
                                           abstol=tolerance, 
                                           reltol=tolerance)
            catch
                redirect_stderr(stderrOld)
                return false
            end
        redirect_stderr(stderrOld)
    else
        du = zeros(Float64, length(sol.prob.u0))
        dp = zeros(Float64, length(sol.prob.p))'
        compute∂G∂u(du, sol[1], sol.prob.p, 0.0, 1)
        onlyObsAtZero = true
    end
    # Technically we can pass compute∂G∂p above to dgdp_discrete. However, odeProb.p often contain 
    # constant parameters which are not a part ode the parameter estimation problem. Sometimes 
    # the gradient for these evaluate to NaN (as they where never thought to be estimated) which 
    # results in the entire gradient evaluating to NaN. Hence, we perform this calculation outside 
    # of the lower level interface. 
    ∂G∂p_ = zeros(Float64, length(sol.prob.p))
    for i in eachindex(timeObserved)     
        if onlyObsAtZero == false                                                                                   
            compute∂G∂p(∂G∂p_, sol(timeObserved[i]), sol.prob.p, timeObserved[i], i)
        else
            compute∂G∂p(∂G∂p_, sol[1], sol.prob.p, timeObserved[i], i)
        end
        dp .+= ∂G∂p_'
    end

    if simulationInfo.haspreEquilibrationConditionId == false
        # In case we do not simulate the ODE for a steady state first we can compute 
        # the initial sensitivites easily via automatic differantitatiom
        St0::Matrix{Float64} = Matrix{Float64}(undef, (length(sol.prob.u0), length(sol.prob.p)))
        ForwardDiff.jacobian!(St0, peTabModel.evalU0, sol.prob.p)
        _gradient = dp .+ du'*St0

    else
        # In case we simulate to a stady state we need to compute a VJP. We use 
        # Zygote pullback to avoid having to having build the Jacobian, rather 
        # we create the yBar function required for the vector Jacobian product.
        _gradient = (dp .+ (evalVJPSS(du)[1])')[:]
    end

    # Thus far have have computed dY/dθ, but for parameters on the log-scale we want dY/dθ_log. We can adjust via;
    # dY/dθ_log = log(10) * θ * dY/dθ
    adjustGradientTransformedParameters!(gradient, _gradient[:], nothing, θ_dynamic, θ_indices, 
                                         simulationConditionId, adjoint=true)                           
    return true                                         
end
