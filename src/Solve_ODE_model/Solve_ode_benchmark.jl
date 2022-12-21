# Solve the ODE models for all experimental condiditions and return the run-time for all solve-calls.
function solveOdeModelAllExperimentalCondBench(prob::ODEProblem, 
                                               changeToExperimentalCondUse!::Function, 
                                               simulationInfo::SimulationInfo,
                                               solver, 
                                               absTol::Float64,
                                               relTol::Float64;
                                               nTSave::Int64=0, 
                                               onlySaveAtTobs::Bool=false,
                                               denseSol::Bool=true, 
                                               savePreEqTime::Bool=false)

    
    absTolSS, relTolSS = simulationInfo.absTolSS, simulationInfo.relTolSS

    bPreEq = 0.0
    bSim = 0.0
    local sucess::Bool = true   
    # In case the model is first simulated to a steady state 
    if simulationInfo.simulateSS == true

        preEqIds = unique(simulationInfo.firstExpIds)
        # Arrays to store steady state (pre-eq) values 
        uAtSS = Matrix{eltype(prob.p)}(undef, (length(prob.u0), length(preEqIds)))
        u0PreSimSS = Matrix{eltype(prob.p)}(undef, (length(prob.u0), length(preEqIds)))

        for i in eachindex(preEqIds)
            
            uAtSSVec = @view uAtSS[:, i]
            u0PreSimSS = @view u0PreSimSS[:, i]

            changeToExperimentalCondUse!(prob.p, prob.u0, preEqIds[i])
            prob = remake(prob, tspan = (0.0, 1e8), p = prob.p[:], u0 = prob.u0[:])
            u0PreSimSS .= prob.u0

            # Terminate if a steady state was not reached in preequilibration simulations 
            bPreEq += @elapsed sol_pre = getSolPreEq(prob, solver, absTol, relTol, absTolSS, relTolSS)
            if sol_pre.retcode == :Terminated
                uAtSSVec .= sol_pre.u[end]
            else
                return false, Inf
            end
        end

        @inbounds for i in eachindex(simulationInfo.conditionIdSol)

            whichPreEq = findfirst(x -> x == simulationInfo.preEqIdSol[i], preEqIds)
            uAtSSVec = @view uAtSS[:, whichPreEq]
            u0PreSimSSVec = @view u0PreSimSS[:, whichPreEq]

            t_max_ss = simulationInfo.tMaxForwardSim[i]
            # Sanity check input. Can only provide either nTsave (points to save solution at) or tSave (number of points to save)
            if nTSave != 0
                saveAtVec = collect(LinRange(0.0, t_max, nTSave))
            else
                saveAtVec = Float64[]
            end
            # Sanity check input. Both tSave and nTsave must be empty or zero in order to be able to output a dense solution.
            if (isempty(saveAtVec) && nTSave == 0) && denseSol == true
                dense = true
            else
                dense = false
            end                                

            changeToExperimentalCondUse!(prob.p, prob.u0, simulationInfo.postEqIdSol[i])
            has_not_changed = (prob.u0 .== u0PreSimSSVec)
            prob.u0[has_not_changed] .= uAtSSVec[has_not_changed]
            probUse = remake(prob, tspan = (0.0, t_max_ss), u0=prob.u0[:], p=prob.p[:])     
            bSim += @elapsed sol = getSolSolveOdeNoSS(probUse, solver, absTol, relTol, absTolSS, relTolSS, t_max_ss, saveAtVec, dense)

            if !(sol.retcode == :Success || sol.retcode == :Terminated)
                sucess = false, Inf
            end
        end

    # In case the model is not first simulated to a steady state 
    elseif simulationInfo.simulateSS == false

        @inbounds for i in eachindex(simulationInfo.firstExpIds)

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
            
            bSim += @elapsed sol = getSolSolveOdeNoSS(probUse, solver, absTol, relTol, absTolSS, relTolSS, t_max_use, saveAtVec, dense)

            if !(sol.retcode == :Success || sol.retcode == :Terminated)
                sucess = false, Inf
            end
        end
    end

    if savePreEqTime == true
        return sucess, bSim + bPreEq
    else
        return sucess, bSim 
    end
end
