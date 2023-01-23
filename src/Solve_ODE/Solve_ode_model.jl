function solveOdeModelAllExperimentalCond!(solArray::Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1},
                                           prob::ODEProblem, 
                                           changeToExperimentalCondUse!::Function, 
                                           simulationInfo::SimulationInfo,
                                           solver::Union{SciMLAlgorithm, Vector{Symbol}},
                                           absTol::Float64,
                                           relTol::Float64, 
                                           calcTStops::Function;
                                           expIDSolve::Array{String, 1} = ["all"],
                                           nTSave::Int64=0, 
                                           onlySaveAtTobs::Bool=false,
                                           denseSol::Bool=true, 
                                           trackCallback::Bool=false)::Bool

    local sucess::Bool = true
    # In case the model is first simulated to a steady state 
    if simulationInfo.simulateSS == true

        # Extract all unique Pre-equlibrium conditions. If the code is run in parallell 
        # (expIDSolve != [["all]]) the number of preEq cond. might be smaller than the 
        # total number of preEq cond.
        if expIDSolve[1] == "all"
            preEqIds = unique(simulationInfo.firstExpIds)
        else
            whichPreEq = findall(x -> x ∈ simulationInfo.conditionIdSol, expIDSolve)
            preEqIds = unique(simulationInfo.preEqIdSol[whichPreEq])
        end
        
        # Arrays to store steady state (pre-eq) values. When computing the sensitivity equations we are solving an 
        # expanded ODE system where along the original ode we solve on of dimension length(u0)*length(p)
        uAtSS = Matrix{eltype(prob.p)}(undef, (length(prob.u0), length(preEqIds)))
        u0PreSimSS = Matrix{eltype(prob.p)}(undef, (length(prob.u0), length(preEqIds)))

        for i in eachindex(preEqIds)
            # Sometimes due to strongly ill-conditioned Jacobian the linear-solve runs 
            # into a domain error. This is treated as integration error.
            try
                whichPreEq = findfirst(x -> x == preEqIds[i], simulationInfo.preEqIdSol)
                simulationInfo.solArrayPreEq[whichPreEq] = solveODEPreEqulibrium!((@view uAtSS[:, i]), 
                                                                                  (@view u0PreSimSS[:, i]), 
                                                                                  prob, 
                                                                                  changeToExperimentalCondUse!, 
                                                                                  preEqIds[i], 
                                                                                  absTol, 
                                                                                  relTol, 
                                                                                  solver, 
                                                                                  simulationInfo.absTolSS, 
                                                                                  simulationInfo.relTolSS)
                if simulationInfo.solArrayPreEq[whichPreEq].retcode != :Terminated
                    return false
                end
            catch e
                if e isa BoundsError
                    println("Bounds error on PreEq solution")
                elseif e isa DomainError
                    println("Bounds error on PreEq solution")
                else
                    rethrow(e)
                end
            end
        end

        @inbounds for i in eachindex(simulationInfo.conditionIdSol)
            
            if expIDSolve[1] != "all" && simulationInfo.conditionIdSol[i] ∉ expIDSolve                
                continue
            end

            whichPreEq = findfirst(x -> x == simulationInfo.preEqIdSol[i], preEqIds)

            # Whether or not we only want to save solution at observed time-points 
            if onlySaveAtTobs == true
                nTSave = 0
                # Extract t-save point for specific condition ID 
                tSave = simulationInfo.tVecSave[simulationInfo.conditionIdSol[i]]
            else
                tSave=Float64[]
            end

            t_max_ss = simulationInfo.tMaxForwardSim[i]
            # See comment above on domain error
            try 
                solArray[i] = solveODEPostEqulibrium(prob, 
                                                    (@view uAtSS[:, whichPreEq]),
                                                    (@view u0PreSimSS[:, whichPreEq]), 
                                                    i,
                                                    changeToExperimentalCondUse!,
                                                    simulationInfo, 
                                                    simulationInfo.postEqIdSol[i],
                                                    absTol,
                                                    relTol,
                                                    t_max_ss,
                                                    solver,
                                                    calcTStops,
                                                    tSave=tSave, 
                                                    nTSave=nTSave, 
                                                    denseSol=denseSol, 
                                                    trackCallback=trackCallback)

                if solArray[i].retcode != :Success
                    sucess = false
                end
            catch e
                if e isa BoundsError
                    println("Bounds error ODE solve")
                    sucess = false
                elseif e isa DomainError
                    println("Bounds error on PreEq solution")
                else
                    rethrow(e)
                end
            end
            if sucess == false
                return false
            end
        end

    # In case the model is not first simulated to a steady state 
    elseif simulationInfo.simulateSS == false

        @inbounds for i in eachindex(simulationInfo.firstExpIds)

            if expIDSolve[1] != "all" && simulationInfo.conditionIdSol[i] ∉ expIDSolve                
                continue
            end

            # Whether or not we only want to save solution at observed time-points 
            firstExpId = simulationInfo.firstExpIds[i]
            if onlySaveAtTobs == true
                nTSave = 0
                # Extract t-save point for specific condition ID 
                tSave = simulationInfo.tVecSave[firstExpId]
            else
                tSave=Float64[]
            end

            t_max = simulationInfo.tMaxForwardSim[i]
            # See comment on domain error above 
            try
                solArray[i] = solveOdeNoSS(prob, 
                                           changeToExperimentalCondUse!, 
                                           i,
                                           simulationInfo,
                                           firstExpId, 
                                           absTol,
                                           relTol, 
                                           solver, 
                                           t_max, 
                                           calcTStops,
                                           nTSave=nTSave, 
                                           tSave=tSave,
                                           denseSol=denseSol, 
                                           trackCallback=trackCallback)
                if !(solArray[i].retcode == :Success || solArray[i].retcode == :Terminated)
                    sucess = false
                end
            catch e
                if e isa BoundsError
                    println("Bounds error ODE solve")
                    sucess = false
                elseif e isa DomainError
                    println("Bounds error on PreEq solution")
                else
                    rethrow(e)
                end
            end
            if sucess == false
                return false
            end                                       
        end
    end

    return sucess
end
function solveOdeModelAllExperimentalCond!(solArray::Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1},
                                           prob::ODEProblem, 
                                           dynParamEst::AbstractVector,
                                           changeToExperimentalCondUsePre!::Function, 
                                           simulationInfo::SimulationInfo,
                                           solver::Union{SciMLAlgorithm, Vector{Symbol}},
                                           absTol::Float64,
                                           relTol::Float64, 
                                           calcTStops::Function;
                                           nTSave::Int64=0, 
                                           onlySaveAtTobs::Bool=false,
                                           expIDSolve::Array{String, 1} = ["all"],
                                           denseSol::Bool=true, 
                                           trackCallback::Bool=false)::Bool

    changeToExperimentalCondUse! = (pVec, u0Vec, expID) -> changeToExperimentalCondUsePre!(pVec, u0Vec, expID, dynParamEst)
    sucess = solveOdeModelAllExperimentalCond!(solArray, 
                                               prob, 
                                               changeToExperimentalCondUse!, 
                                               simulationInfo, 
                                               solver, 
                                               absTol, 
                                               relTol,
                                               calcTStops,
                                               nTSave=nTSave, 
                                               denseSol=denseSol, 
                                               onlySaveAtTobs=onlySaveAtTobs, 
                                               expIDSolve=expIDSolve, 
                                               trackCallback=trackCallback)

    return sucess
end
function solveOdeModelAllExperimentalCond!(solArray::Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1},
                                           SMat::Matrix{Float64},
                                           prob::ODEProblem, 
                                           dynParamEst::AbstractVector,
                                           changeToExperimentalCondUsePre!::Function, 
                                           changeToModelParamUse!::Function,
                                           simulationInfo::SimulationInfo,
                                           solver::Union{SciMLAlgorithm, Vector{Symbol}},
                                           absTol::Float64,
                                           relTol::Float64, 
                                           calcTStops::Function;
                                           nTSave::Int64=0, 
                                           onlySaveAtTobs::Bool=false,
                                           expIDSolve::Array{String, 1} = ["all"],
                                           denseSol::Bool=true, 
                                           trackCallback::Bool=false)::Bool

    function getSenseMat(paramEst::AbstractVector)

        odeProbUse = remake(prob, p = convert.(eltype(paramEst), prob.p), u0 = convert.(eltype(paramEst), prob.u0))
        
        changeToModelParamUse!(odeProbUse.p, odeProbUse.u0, paramEst)
        changeToExperimentalCondUse! = (pVec, u0Vec, expID) -> changeToExperimentalCondUsePre!(pVec, u0Vec, expID, paramEst)

        sucess = solveOdeModelAllExperimentalCond!(solArray, 
                                                   odeProbUse, 
                                                   changeToExperimentalCondUse!, 
                                                   simulationInfo, 
                                                   solver, 
                                                   absTol, 
                                                   relTol,
                                                   calcTStops,
                                                   nTSave=nTSave, 
                                                   denseSol=denseSol, 
                                                   onlySaveAtTobs=onlySaveAtTobs, 
                                                   expIDSolve=expIDSolve, 
                                                   trackCallback=trackCallback)

        nStates = length(prob.u0)                                                                                                      
        if sucess != true
            nObs = sum(length(simulationInfo.tVecSave[simulationInfo.conditionIdSol[i]]) for i in eachindex(simulationInfo.conditionIdSol))
            return zeros(nStates, nObs)
        end
        
        arrRet = Array{eltype(prob.u0), 2}(undef, (nStates, 0))
        for i in eachindex(simulationInfo.conditionIdSol) 
            if expIDSolve[1] != "all" && simulationInfo.conditionIdSol[i] ∉ expIDSolve        
                tSave = simulationInfo.tVecSave[simulationInfo.conditionIdSol[i]]        
                arrAdd = zeros(eltype(prob.u0), (nStates, length(tSave)))
            else
                arrAdd = Array(solArray[i])
            end
            arrRet = hcat(arrRet, arrAdd)
        end
        return arrRet
    end
    # Compute sensitivity matrix via forward mode automatic differentation 
    ForwardDiff.jacobian!(SMat, getSenseMat, dynParamEst)

    # Check retcode 
    sucess = true
    for i in eachindex(simulationInfo.conditionIdSol) 
        if expIDSolve[1] != "all" && simulationInfo.conditionIdSol[i] ∉ expIDSolve        
            continue
        end
        if !(solArray[i].retcode == :Success || solArray[i].retcode == :Terminated)
            sucess = false
        end
    end
    
    return sucess
end


function solveOdeModelAllExperimentalCond(prob::ODEProblem, 
                                          changeToExperimentalCondUse!::Function, 
                                          simulationInfo::SimulationInfo,
                                          solver::Union{SciMLAlgorithm, Vector{Symbol}},
                                          absTol::Float64,
                                          relTol::Float64,
                                          calcTStops::Function;
                                          nTSave::Int64=0, 
                                          onlySaveAtTobs::Bool=false,
                                          denseSol::Bool=true, 
                                          trackCallback::Bool=false)

    local solArray
    # Compute the number of ODE-solutions to cover all experimental conditions 
    if simulationInfo.simulateSS == true
        nShiftId = Int(sum([length(simulationInfo.shiftExpIds[i]) for i in eachindex(simulationInfo.shiftExpIds)]))
        solArray = Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1}(undef, nShiftId)
    else
        nExperimentalCond = Int64(length(simulationInfo.firstExpIds))
        solArray = Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1}(undef, nExperimentalCond)
    end
        
    success = solveOdeModelAllExperimentalCond!(solArray, 
                                                prob, 
                                                changeToExperimentalCondUse!, 
                                                simulationInfo, 
                                                solver, 
                                                absTol, 
                                                relTol,
                                                calcTStops,
                                                nTSave=nTSave, 
                                                denseSol=denseSol, 
                                                onlySaveAtTobs=onlySaveAtTobs, 
                                                trackCallback=trackCallback)

    return solArray, success
end


# Solve the ODE system for one experimental conditions in a Zygote compatible manner.
function solveOdeModelAtExperimentalCondZygote(prob::ODEProblem, 
                                               conditionId::String,
                                               dynParamEst,
                                               t_max,
                                               changeToExperimentalCondUsePre::Function, 
                                               measurementData::MeasurementData,
                                               simulationInfo::SimulationInfo,
                                               solver::Union{SciMLAlgorithm, Vector{Symbol}}, 
                                               absTol::Float64, 
                                               relTol::Float64, 
                                               sensealg, 
                                               calcTStops::Function)

    changeToExperimentalCondUse = (pVec, u0Vec, expID) -> changeToExperimentalCondUsePre(pVec, u0Vec, expID, dynParamEst)                                               

    # For storing ODE solution (required for split gradient computations)
    whichCondID = findfirst(x -> x == conditionId, simulationInfo.conditionIdSol)

    # In case the model is first simulated to a steady state 
    local success = true
    if simulationInfo.simulateSS == true

        firstExpId = measurementData.preEqCond[measurementData.iPerConditionId[conditionId][1]]
        shiftExpId = measurementData.simCond[measurementData.iPerConditionId[conditionId][1]]
        tSave = simulationInfo.tVecSave[conditionId]            
    
        u0Pre = prob.u0[:]                                        
        pUsePre, u0UsePre = changeToExperimentalCondUse(prob.p, prob.u0, firstExpId)
        probUsePre = remake(prob, tspan=(0.0, 1e8), u0 = convert.(eltype(dynParamEst), u0UsePre), p = convert.(eltype(dynParamEst), pUsePre))
        ssProb = SteadyStateProblem(probUsePre)
        solSS = solve(ssProb, DynamicSS(solver, abstol=simulationInfo.absTolSS, reltol=simulationInfo.relTolSS), abstol=absTol, reltol=relTol)

        # Terminate if a steady state was not reached in preequilibration simulations 
        if solSS.retcode != :Success
            return sol_pre, false
        end

        # Change to parameters for the post steady state parameters 
        pUsePost, u0UsePostTmp = changeToExperimentalCondUse(prob.p, prob.u0, shiftExpId)
        
        # Given the standard the experimentaCondition-file can change the initial values for a state 
        # whose value was changed in the preequilibration-simulation. The experimentalCondition
        # value is prioritized by only changing u0 to the steady state value for those states  
        # that were not affected by change to shiftExpId.
        hasNotChanged = (u0UsePostTmp .== u0Pre)
        u0UsePost = [hasNotChanged[i] == true ? solSS[i] : u0UsePostTmp[i] for i in eachindex(u0UsePostTmp)]
        probUsePost = remake(prob, tspan=(0.0, t_max), u0 = convert.(eltype(dynParamEst), u0UsePost), p = convert.(eltype(dynParamEst), pUsePost))
        
        # Different funcion calls to solve are required if a solver or a Alg-hint are provided. 
        # The preequilibration simulations are terminated upon a steady state using the TerminateSteadyState callback.
        tStops = calcTStops(probUsePost.u0, probUsePost.p)
        solveCallPost = (prob) -> solve(prob, 
                                        solver, 
                                        abstol=absTol, 
                                        reltol=relTol, 
                                        saveat=tSave,
                                        sensealg=sensealg, 
                                        callback=simulationInfo.callbacks[whichCondID], 
                                        tstops=tStops)
        
        
        sol = solveCallPost(probUsePost)

        if sol.retcode != :Success
            sucess = false
        end

    # In case the model is not first simulated to a steady state 
    elseif simulationInfo.simulateSS == false

        firstExpId = measurementData.simCond[measurementData.iPerConditionId[conditionId][1]]
        tSave = simulationInfo.tVecSave[conditionId]
        t_max_use = isinf(t_max) ? 1e8 : t_max

        pUse, u0Use = changeToExperimentalCondUse(prob.p, prob.u0, firstExpId)
        probUse = remake(prob, tspan=(0.0, t_max_use), u0 = convert.(eltype(dynParamEst), u0Use), p = convert.(eltype(dynParamEst), pUse))

        # Different funcion calls to solve are required if a solver or a Alg-hint are provided. 
        # If t_max = inf the model is simulated to steady state using the TerminateSteadyState callback.
        tStops = calcTStops(probUse.u0, probUse.p)
        if !(typeof(solver) <: Vector{Symbol}) && isinf(t_max)
            solveCall = (probArg) -> solve(probArg, 
                                           solver, 
                                           abstol=absTol, 
                                           reltol=relTol, 
                                           save_on=false,
                                           save_end=true, 
                                           dense=dense, 
                                           callback=TerminateSteadyState(absTolSS, relTolSS))

        elseif !(typeof(solver) <: Vector{Symbol}) && !isinf(t_max)
            solveCall = (probArg) -> solve(probArg, 
                                           solver, 
                                           p = pUse,
                                           u0 = u0Use,
                                           abstol=absTol, 
                                           reltol=relTol, 
                                           saveat=tSave, 
                                           sensealg=sensealg, 
                                           callback=simulationInfo.callbacks[whichCondID], 
                                           tstops=tStops)
        else
            println("Error : Solver option does not exist")        
        end
        
        sol = solveCall(probUse)

        Zygote.@ignore simulationInfo.solArray[whichCondID] = sol

        if typeof(sol) <: ODESolution && !(sol.retcode == :Success || sol.retcode == :Terminated)
            sucess = false
        end
    end
    
    return sol, success
end


function solveODEPreEqulibrium!(uAtSSVec::AbstractVector, 
                                u0PreSimSS::AbstractVector,
                                prob::ODEProblem, 
                                changeToExperimentalCondUse!::Function,
                                firstExpId::String, 
                                absTol::Float64, 
                                relTol::Float64, 
                                solver::Union{SciMLAlgorithm, Vector{Symbol}},
                                absTolSS::Float64, 
                                relTolSS::Float64)::Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}

    # Change to parameters for the preequilibration simulations 
    changeToExperimentalCondUse!(prob.p, prob.u0, firstExpId)
    prob = remake(prob, tspan = (0.0, 1e8), p = prob.p[:], u0 = prob.u0[:])
    u0PreSimSS .= prob.u0

    # Terminate if a steady state was not reached in preequilibration simulations 
    sol_pre = getSolPreEq(prob, solver, absTol, relTol, absTolSS, relTolSS)
    if sol_pre.retcode == :Terminated
        uAtSSVec .= sol_pre.u[end]
    end
    return sol_pre
end


function getCallbackSet(prob::ODEProblem, 
                        simulationInfo::SimulationInfo, 
                        whichForwardSol::Int64, 
                        trackCallback::Bool)::SciMLBase.DECallback
    
    if trackCallback == true
        cbSet = SciMLSensitivity.track_callbacks(simulationInfo.callbacks[whichForwardSol], prob.tspan[1], 
                                                 prob.u0, prob.p, simulationInfo.sensealg) 
        simulationInfo.callbacks[whichForwardSol] = cbSet
    end
    return simulationInfo.callbacks[whichForwardSol]
end


function solveODEPostEqulibrium(prob::ODEProblem, 
                                uAtSS::AbstractVector,
                                u0PreSimSS::AbstractVector,
                                whichForwardSol::Int64,
                                changeToExperimentalCondUse!::Function, 
                                simulationInfo::SimulationInfo,
                                shiftExpId::String,
                                absTol::Float64, 
                                relTol::Float64,
                                t_max_ss::Float64,
                                solver::Union{SciMLAlgorithm, Vector{Symbol}},
                                calcTStops::Function;
                                tSave=Float64[], 
                                nTSave=0, 
                                denseSol::Bool=true, 
                                trackCallback::Bool=false)::Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}

    # Sanity check input. Can only provide either nTsave (points to save solution at) or tSave (number of points to save)
    if length(tSave) != 0 && nTSave != 0
        println("Error : Can only provide tSave (vector to save at) or nTSave as saveat argument to solvers")
    elseif nTSave != 0
        saveAtVec = collect(LinRange(0.0, t_max_ss, nTSave))
    else
        saveAtVec = tSave
    end

    # Sanity check input. Both tSave and nTsave must be empty or zero in order to be able to output a dense solution.
    if (isempty(tSave) && nTSave == 0) && denseSol == true
        dense = true
    else
        dense = false
    end                                

    changeToExperimentalCondUse!(prob.p, prob.u0, shiftExpId)
    # Sometimes the experimentaCondition-file changes the initial values for a state 
    # whose value was changed in the preequilibration-simulation. The experimentaCondition
    # value is prioritized by only changing u0 to the steady state value for those states  
    # that were not affected by change to shiftExpId.
    has_not_changed = (prob.u0 .== u0PreSimSS)
    prob.u0[has_not_changed] .= uAtSS[has_not_changed]

    # Here it is IMPORTANT that we copy prob.p[:] else different experimental conditions will 
    # share the same parameter vector p. This will, for example, cause the lower level adjoint 
    # sensitivity interface to fail.
    t_max_ss = isinf(t_max_ss) ? 1e8 : t_max_ss
    probUse = remake(prob, tspan = (0.0, t_max_ss), u0=prob.u0[:], p=prob.p[:])     
    probUse.u0 .= prob.u0[:] # This is needed due to remake for ODE-sensitivity problem does not correctly map u0

    # If case of adjoint sensitivity analysis we need to track the callback to get correct gradients 
    tStops = calcTStops(probUse.u0, probUse.p)
    callBackSet = getCallbackSet(probUse, simulationInfo, whichForwardSol, trackCallback)
    sol = getSolSolveOdeNoSS(probUse, solver, absTol, relTol, simulationInfo.absTolSS, simulationInfo.relTolSS, 
                             t_max_ss, saveAtVec, dense, callBackSet, tStops)

    return sol                                     
end


function solveOdeNoSS(prob::ODEProblem, 
                      changeToExperimentalCondUse!::Function, 
                      whichForwardSol::Int64, 
                      simulationInfo::SimulationInfo,
                      firstExpId::String, 
                      absTol::Float64, 
                      relTol::Float64,
                      solver::Union{SciMLAlgorithm, Vector{Symbol}},
                      t_max::Float64, 
                      calcTStops::Function; 
                      tSave=Float64[], 
                      nTSave::Int64=0, 
                      denseSol::Bool=true, 
                      trackCallback::Bool=false)::Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}

    # Sanity check input. Can only provide either nTsave (points to save solution at) or tSave (number of points to save)
    if length(tSave) != 0 && nTSave != 0
        println("Error : Can only provide tSave (vector to save at) or nTSave as saveat argument to solvers")
    elseif nTSave != 0
        saveAtVec = collect(LinRange(0.0, t_max, nTSave))
    else
        saveAtVec = tSave
    end

    # Sanity check input. Both tSave and nTsave must be empty or zero in order to be able to output a dense solution.
    if (isempty(tSave) && nTSave == 0) && denseSol == true
        dense = true
    else
        dense = false
    end                                

    # Change experimental condition 
    t_max_use = isinf(t_max) ? 1e8 : t_max
    changeToExperimentalCondUse!(prob.p, prob.u0, firstExpId)
    probUse = remake(prob, tspan=(0.0, t_max_use), u0 = prob.u0[:], p = prob.p[:])
    probUse.u0 .= prob.u0[:] # Required as due to a bub in remake it does not handle Senstivity-problems well 
    
    # If case of adjoint sensitivity analysis we need to track the callback 
    callBackSet = getCallbackSet(probUse, simulationInfo, whichForwardSol, trackCallback)
    tStops = calcTStops(probUse.u0, probUse.p)

    sol = getSolSolveOdeNoSS(probUse, solver, absTol, relTol, simulationInfo.absTolSS, simulationInfo.relTolSS, 
                             t_max_use, saveAtVec, dense, callBackSet, tStops)

    return sol
end


function getSolSolveOdeNoSS(prob::ODEProblem, 
                            solver::Vector{Symbol}, 
                            absTol::Float64, 
                            relTol::Float64,
                            absTolSS::Float64, 
                            relTolSS::Float64,
                            t_max::Float64, 
                            saveAtVec::Vector{Float64}, 
                            dense::Bool, 
                            callBackSet::SciMLBase.DECallback, 
                            tStops::Vector{<:AbstractFloat})::Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}

    # Different funcion calls to solve are required if a solver or a Alg-hint are provided. 
    # If t_max = inf the model is simulated to steady state using the TerminateSteadyState callback.
    if isinf(t_max) || t_max == 1e8
        sol = solve(prob, alg_hints=solver, abstol=absTol, reltol=relTol, save_on=false, save_end=true, dense=dense, callback=TerminateSteadyState(absTolSS, relTolSS))
    else 
        sol = solve(prob, alg_hints=solver, abstol=absTol, reltol=relTol, saveat=saveAtVec, dense=dense, tstops=tStops, callback=callBackSet)   
    end
    return sol
end
function getSolSolveOdeNoSS(prob::ODEProblem, 
                            solver::SciMLAlgorithm, 
                            absTol::Float64, 
                            relTol::Float64,
                            absTolSS::Float64, 
                            relTolSS::Float64,
                            t_max::Float64, 
                            saveAtVec::Vector{Float64}, 
                            dense::Bool, 
                            callBackSet::SciMLBase.DECallback, 
                            tStops::Vector{<:AbstractFloat})::Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}

    # Different funcion calls to solve are required if a solver or a Alg-hint are provided. 
    # If t_max = inf the model is simulated to steady state using the TerminateSteadyState callback.
    if isinf(t_max) || t_max == 1e8
        sol = solve(prob, solver, abstol=absTol, reltol=relTol, save_on=false, save_end=true, dense=dense, callback=TerminateSteadyState(absTolSS, relTolSS))
    else
        sol = solve(prob, solver, abstol=absTol, reltol=relTol, saveat=saveAtVec, dense=dense, tstops=tStops, callback=callBackSet)
    end
    return sol
end


function getSolPreEq(prob::ODEProblem,
                     solver::Vector{Symbol},
                     absTol::Float64, 
                     relTol::Float64, 
                     absTolSS::Float64, 
                     relTolSS::Float64)::Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}

    sol = solve(prob, alg_hints=solver, abstol=absTol, reltol=relTol, dense=false, callback=TerminateSteadyState(absTolSS, relTolSS))
    return sol
end
function getSolPreEq(prob::ODEProblem,
                     solver::SciMLAlgorithm,
                     absTol::Float64, 
                     relTol::Float64, 
                     absTolSS::Float64, 
                     relTolSS::Float64)::Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}

    sol = solve(prob, solver, abstol=absTol, reltol=relTol, dense=false, callback=TerminateSteadyState(absTolSS, relTolSS))
    return sol
end


# Change experimental condition when running parameter estimation. A lot of heavy lifting here is done by 
# an index which correctly maps parameters for a conditionId to the ODEProblem.
function changeExperimentalCondEst!(pODEProblem::AbstractVector, 
                                    u0::AbstractVector, 
                                    conditionId::String, 
                                    θ_dynamic::AbstractVector,
                                    peTabModel::PeTabModel, 
                                    θ_indices::ParameterIndices;
                                    computeForwardSensitivites::Bool=false)

    mapConditionId = θ_indices.mapsConiditionId[Symbol(conditionId)] 

    # Constant parameters 
    pODEProblem[mapConditionId.iODEProblemConstantParameters] .= mapConditionId.constantParameters

    # Parameters to estimate 
    pODEProblem[mapConditionId.iODEProblemθDynamic] .= θ_dynamic[mapConditionId.iθDynamic]

    # Given changes in parameters initial values might have to be re-evaluated 
    nModelStates = length(peTabModel.stateNames)
    peTabModel.evalU0!((@view u0[1:nModelStates]), pODEProblem) 

    # Account for any potential events (callbacks) which are active at time zero
    for f! in peTabModel.checkCallbackActive
        f!(u0, pODEProblem)
    end

    # In case an experimental condition maps directly to the initial value of a state. 
    if !isempty(mapConditionId.iODEProblemConstantStates)
        u0[mapConditionId.iODEProblemConstantStates] .= mapConditionId.constantsStates
    end

    # In case we solve the forward sensitivity equations we must adjust the initial sensitives 
    # by computing the jacobian at t0
    if computeForwardSensitivites == true
        St0::Matrix{Float64} = Matrix{Float64}(undef, (nModelStates, length(pODEProblem)))
        ForwardDiff.jacobian!(St0, peTabModel.evalU0, pODEProblem)
        u0[(nModelStates+1):end] .= vec(St0)
    end

    return nothing
end                                 
                             

function changeExperimentalCondEst(pODEProblem::AbstractVector, 
                                   u0::AbstractVector, 
                                   conditionId::String, 
                                   θ_dynamic::AbstractVector,
                                   peTabModel::PeTabModel, 
                                   θ_indices::ParameterIndices)

    mapConditionId = θ_indices.mapsConiditionId[Symbol(conditionId)] 
    
    # For a non-mutating way of mapping constant parameters 
    function iConstantParam(iUse::Integer)
        whichIndex = findfirst(x -> x == iUse, mapConditionId.iODEProblemConstantParameters)
        return whichIndex
    end
    # For a non-mutating mapping of parameters to estimate 
    function iParametersEst(iUse::Integer)
        whichIndex = findfirst(x -> x == iUse, mapConditionId.iODEProblemθDynamic)
        return mapConditionId.iθDynamic[whichIndex]
    end
    
    # Constant parameters 
    _pODEProblem = [i ∈ mapConditionId.iODEProblemConstantParameters ? mapConditionId.constantParameters[iConstantParam(i)] : pODEProblem[i] for i in eachindex(pODEProblem)]
    
    # Parameters to estimate 
    __pODEProblem = [i ∈ mapConditionId.iODEProblemθDynamic ? θ_dynamic[iParametersEst(i)] : _pODEProblem[i] for i in eachindex(_pODEProblem)]    
    
    # When using AD as Zygote we must use the non-mutating version of evalU0
    _u0 = peTabModel.evalU0(__pODEProblem) 

    # In case an experimental condition maps directly to the initial value of a state. 
    # Very rare, will fix if ever needed for slow Zygote 
    if !isempty(mapConditionId.iODEProblemConstantStates)
        u0[mapConditionId.iODEProblemConstantStates] .= mapConditionId.constantsStates
    end

    # Account for any potential events 
    for f! in peTabModel.checkCallbackActive
        f!(_u0, __pODEProblem)
    end

    return __pODEProblem, _u0
end     


function getFileODEvalues(peTabModel::PeTabModel)
  
    # Change model parameters 
    experimentalConditionsFile, measurementDataFile, parameterDataFile, observablesDataFile = readDataFiles(peTabModel.dirModel, readObs=true)
    parameterInfo = processParameterData(parameterDataFile)
    measurementData = processMeasurementData(measurementDataFile, observablesDataFile) 
    θ_indices = computeIndicesθ(parameterInfo, measurementData, peTabModel.odeSystem, experimentalConditionsFile)

    θ_estNames = θ_indices.θ_estNames
    θ_est = parameterInfo.paramVal[findall(x -> x ∈ θ_estNames, parameterInfo.parameterID)]

    return θ_est[θ_indices.iθ_dynamic]
end