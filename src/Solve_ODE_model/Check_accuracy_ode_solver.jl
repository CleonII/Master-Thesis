# Functions for checking accuracy by comparing ODE-solution against high accuracy solution. 
# Must change before doing benchmark - this file will go through a lot of changes. 

function calcSqErrVal(solHigh, solCompare, t_max)

    sqErr::Float64 = 0.0
    
    if isinf(t_max)
        # Only match end point when data i matched against steady state value 
        # TODO : Fix if some data points for experiment are not inf 
        sqErr = sum((solHigh.u[end] .- solCompare.u[end]).^2)
    elseif length(solHigh.t) != length(solCompare.t)
        # This likelly happens due to events in the solution 
        sqErr = sum((Array(solCompare(solHigh.t)) - Array(solHigh)).^2)
        sqErr -= sum((Array(solHigh(solHigh.t)) - Array(solHigh)).^2) # Protect against interpolation errors with events 
    else
        sqErr = sum((solCompare[:,:] - solHigh[:,:]).^2)
        sqErr -= sum((Array(solHigh(solHigh.t)) - Array(solHigh)).^2) # Protect against interpolation errors with events 
    end

    return sqErr
end


# Solve and ODE-model (in prob) for all experimental conditions specifed by 
# first and shift experimental IDs using a given solver tolerance and 
# ode solver and compute the squared error against a reference soluation obtained 
# solveOdeModelAllCond and stored in solArray. Based on whether or not the model 
# should be simulated to a steady state it calls on, of two, possible solvers. 
# changeToCondUse! ensures that ODE-parameters are changed correctly. 
# Returns inf if the solver fails. 
# Note : only for benchmarking ODE-solvers. 
function calcSqErr(prob::ODEProblem, 
                   changeToCondUse!::Function, 
                   solArray,
                   simulateSS::Bool, 
                   measurementData::DataFrame,
                   firstExpIds, 
                   shiftExpIds, 
                   tol::Float64, 
                   solver)

    # Check if model can be solved (without using forced stops for integrator)
    solArrayTmp, sucess = solveOdeModelAllCond(prob, changeToCondUse!, simulateSS, measurementData, firstExpIds, shiftExpIds, tol, solver, denseArg=false)
    solArrayTmp = 0
    if sucess == false
        return Inf 
    end

    sqErr::Float64 = 0.0
    local couldSolve = true
    if simulateSS == true
        nShiftId = Int(sum([length(shiftExpIds[i]) for i in eachindex(shiftExpIds)]))
        k = 1

        for i in eachindex(firstExpIds)
            for j in eachindex(shiftExpIds[i])
                firstExpId = firstExpIds[i]
                shiftExpId = shiftExpIds[i][j]
                t_max_ss = getTimeMax(measurementData, shiftExpId)
                solCompare = solveOdeSS(prob, changeToCondUse!, firstExpId, shiftExpId, tol, t_max_ss, solver, tSave=solArray[k].t)
                
                if solCompare.retcode != :Success && solCompare.t[end] != solArray[k].t[end]
                    couldSolve = false
                    break 
                end

                sqErr += calcSqErrVal(solArray[k], solCompare, t_max_ss)
                k += 1
            end
        end

    elseif simulateSS == false
        nExperimentalCond = Int64(length(firstExpIds))

        for i in eachindex(firstExpIds)
            firstExpId = firstExpIds[i]
            t_max = getTimeMax(measurementData, firstExpId)
            solCompare = solveOdeNoSS(prob, changeToCondUse!, firstExpId, tol, solver, t_max, tSave=solArray[i].t)

            if isinf(t_max) && solCompare.retcode != :Terminated
                couldSolve = false
                break
            end
            if !isinf(t_max) && !(solCompare.retcode == :Success && solCompare.t[end] == solArray[i].t[end])
                couldSolve = false
                break 
            end

            sqErr += calcSqErrVal(solArray[i], solCompare, t_max)
        end
    end

    if couldSolve == false
        return Inf
    end

    return sqErr
end
    