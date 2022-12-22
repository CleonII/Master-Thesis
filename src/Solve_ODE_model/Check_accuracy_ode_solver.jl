"""
    calcHighAccOdeSolution(prob::ODEProblem, 
                           changeToExperimentalCondUse!::Function, 
                           measurementData::DataFrame,
                           simulationInfo::SimulationInfo,
                           tol::Float64=1e-15) 

    For a PeTab ODE model with the parameter values in prob compute a high accuracy ODE solution 
    with absTol=relTol=tol using BigFloats for the ode-problem. The ODE is solved for all experimental 
    conditions specified in the PeTab files. Returns an array with ODE-solution and a bool which is 
    true if the ODE model could be solved. 

    By default the function first tries to solve the ODE problem using a non-stiff solver AutoVern9(Rodas4P()). 
    In case the non-stiff solver fails the Rodas4P() is used to compute the high accuracy solution. If both solvers
    fail status fail is returned.
"""
function calcHighAccOdeSolution(prob::ODEProblem, 
                                changeToExperimentalCondUse!::Function, 
                                simulationInfo::SimulationInfo;
                                absTol::Float64=1e-15, 
                                relTol::Float64=1e-15,
                                nTSave=100) 
                                
    bigFloatOdeProb = createBigFloatODEProblem(prob)

    solverNonStiff = AutoVern9(Rodas4P())
    solverStiff = Rodas4P()
    solverStiff = KenCarp58()
    
    local solArrayHighAcc
    local sucessSolver
    try 
        solArrayHighAcc, sucessSolver = solveOdeModelAllExperimentalCond(bigFloatOdeProb, changeToExperimentalCondUse!, simulationInfo, solverNonStiff, absTol, relTol; nTSave=nTSave)                                      
    catch 
        solArrayHighAcc, sucessSolver = solveOdeModelAllExperimentalCond(bigFloatOdeProb, changeToExperimentalCondUse!, simulationInfo, solverStiff, absTol, relTol; nTSave=nTSave)                                      
    end
    GC.gc()

    # In cases the non-stiff solver just fails 
    if sucessSolver != true
        println("Failed with composite solver - moving on to stiff solver")
        try 
            solArrayHighAcc, sucessSolver = solveOdeModelAllExperimentalCond(bigFloatOdeProb, changeToExperimentalCondUse!, simulationInfo, solverStiff, absTol, relTol; nTSave=nTSave)                                      
        catch
            sucessSolver = false
            solArrayHighAcc = nothing
        end
    end

    return solArrayHighAcc, sucessSolver
end


"""
    calcAccuracyOdeSolver(prob::ODEProblem, 
                          solArrayHighAccuracy::Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1},
                          changeToExperimentalCondUse!::Function, 
                          measurementData::DataFrame,
                          simulationInfo::SimulationInfo,
                          solver,
                          tol::Float64)::Float64

    Check the accuracy of an ODE solver at specific tol=absTol=relTol for a PeTab ODE model (odeProb) by 
    for each experimental condition computing the squared sum difference against a high accuracy ODE 
    solution (for each experimental condition) stored in solArrayHighAccuracy. 

    Recomended to compute high accuracy solution with small tolerances (1e-15) using a high accuracy solver
    and BigFloat. 
"""
function calcAccuracyOdeSolver(prob::ODEProblem, 
                               solArrayHighAccuracy::Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1},
                               changeToExperimentalCondUse!::Function, 
                               simulationInfo::SimulationInfo,
                               solver,
                               absTol::Float64, 
                               relTol::Float64)::Float64
                   
    # Check if model can be solved (without using forced stops for integrator 
    # as this can make the solver converge).
    solArrayTmp, sucess = solveOdeModelAllExperimentalCond(prob, changeToExperimentalCondUse!, simulationInfo, solver, absTol, relTol, nTSave=100)
    solArrayTmp = 0
    if sucess == false
        return Inf 
    end

    # Compute the square error by comparing the solver against the high accuracy solArray. The comparison 
    # is done at the time-points in the high accuracy solver. 
    sqErr::Float64 = 0.0
    local couldSolve = true

    if simulationInfo.simulateSS == true

        preEqIds = unique(simulationInfo.firstExpIds)    
        # Arrays to store steady state (pre-eq) values 
        uAtSS = Matrix{eltype(prob.p)}(undef, (length(prob.u0), length(preEqIds)))
        u0PreSimSS = Matrix{eltype(prob.p)}(undef, (length(prob.u0), length(preEqIds)))

        for i in eachindex(preEqIds)
            # Sometimes due to strongly ill-conditioned Jacobian the linear-solve runs 
            # into a domain error.
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
                couldSolve == false
                break
            end        
        end

        for i in eachindex(simulationInfo.conditionIdSol)

            if couldSolve == false
                break
            end
                
            whichPreEq = findfirst(x -> x == simulationInfo.preEqIdSol[i], preEqIds)
            t_max_ss = simulationInfo.tMaxForwardSim[i]
            # See comment above on domain error
            solCompare = solveODEPostEqulibrium(prob, 
                                                (@view uAtSS[:, whichPreEq]),
                                                (@view u0PreSimSS[:, whichPreEq]), 
                                                changeToExperimentalCondUse!,
                                                simulationInfo.postEqIdSol[i], 
                                                absTol,
                                                relTol,
                                                t_max_ss,
                                                simulationInfo.absTolSS, 
                                                simulationInfo.relTolSS,
                                                solver,
                                                tSave=solArrayHighAccuracy[i].t)

            if solCompare.retcode != :Success
                couldSolve = false
                break
            end

            sqErr += calcSqErrVal(solArrayHighAccuracy[i], solCompare, t_max_ss)
        end

    elseif simulationInfo.simulateSS == false

        for i in eachindex(simulationInfo.firstExpIds)
            
            firstExpId = simulationInfo.firstExpIds[i]
            t_max = simulationInfo.tMaxForwardSim[i]
            solCompare = solveOdeNoSS(prob, changeToExperimentalCondUse!, firstExpId, absTol, relTol, solver, t_max, tSave=solArrayHighAccuracy[i].t)

            # In case t_max = Inf only calculate sqErr if the model could reach a steady state 
            if isinf(t_max) && solCompare.retcode != :Terminated
                couldSolve = false
                break
            end
            # In case t_max != Inf compute sqErr if Success retcode, and if solCompare made it to the end-point.
            if !isinf(t_max) && !(solCompare.retcode == :Success && solCompare.t[end] == solArrayHighAccuracy[i].t[end])
                couldSolve = false
                break 
            end

            sqErr += calcSqErrVal(solArrayHighAccuracy[i], solCompare, t_max)
        end
    end

    if couldSolve == false
        return Inf
    end

    return sqErr
end
    

"""
    calcSqErrVal(solHighAccuracy, sol, t_max::Float64)::Float64

    Helper function to compute the squared sum difference betweeen two 
    ODESolutions solved to t_max. Here solHighAccuracy is meant to 
    be a high accuracy solution. 
"""
function calcSqErrVal(solHighAccuracy, sol, t_max::Float64)::Float64

    sqErr::Float64 = 0.0
    
    # If t_max = Inf only compare at the end-point (steady state value)
    if isinf(t_max)
        sqErr = sum((solHighAccuracy.u[end] .- sol.u[end]).^2)
    
    # In case solHighAccuracy and sol do not have equal number of time-points 
    # there are events in the solution causing problem (events forces saving). 
    # Hence sol is interpolated to solHighAccuracy values, and to avoid interpolation
    # errors around events the even sqErr is subtracted when solHighAccuracy is interpolated 
    # upon itself.
    elseif length(solHighAccuracy.t) != length(sol.t)
        sqErr = sum((Array(sol(solHighAccuracy.t)) - Array(solHighAccuracy)).^2)
        sqErr -= sum((Array(solHighAccuracy(solHighAccuracy.t)) - Array(solHighAccuracy)).^2) 
    
    # Same strategy as above to avoid interpolation errors with events.
    else
        sqErr = sum((sol[:,:] - solHighAccuracy[:,:]).^2)
        sqErr -= sum((Array(solHighAccuracy(solHighAccuracy.t)) - Array(solHighAccuracy)).^2)
    end

    return sqErr
end


"""
    createBigFloatODEProblem(peTabModel::PeTabModel)::ODEProblem

    From a PeTab model create its corresponding ODE-problem with 
    BigFloat (long double). 
"""
function createBigFloatODEProblem(peTabModel::PeTabModel)::ODEProblem
    
    paramMapUse = convert(Vector{Pair{Num, BigFloat}}, peTabModel.paramMap)

    if typeof(peTabModel.stateMap) == Vector{Pair{Num, Num}}
        # This means that u0 has some initializations that depend on par, 
        # since the values of par is BigFloat, u0 will be BigFloat
        stateMapUse = peTabModel.stateMap
    elseif typeof(peTabModel.stateMap) == Vector{Pair{Num, Float64}}
        stateMapUse = convert(Vector{Pair{Num, BigFloat}}, peTabModel.stateMap)
    else
        println("Error: Could not parse PeTab model into BigFloat problem")
    end

    prob = ODEProblem(peTabModel.odeSystem, stateMapUse, (0.0, 1e8), paramMapUse, jac=true)
    return prob
end
"""
    createBigFloatODEProblem(odeProb::ODEProblem)::ODEProblem

    Convert an ODE problem with arbitrary float to one an ODE-problem with BigFloat.
"""
function createBigFloatODEProblem(odeProb::ODEProblem)::ODEProblem
    return remake(odeProb, p = convert.(BigFloat, odeProb.p), u0 = convert.(BigFloat, odeProb.u0))
end