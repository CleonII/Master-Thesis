# In a DataFrame data with colSearch find which row of colSearch 
# has the value ExpId. Helper function to navigate experimentalCondition 
# files. 
function getRowExpId(expId, data::DataFrame; colSearch="conditionId") 
    return findfirst(x -> x == expId, data[!, colSearch])
end


# Check if a string is a number. 
function isNumber(a::String)
    return tryparse(Float64, a) !== nothing
end


# Change ODE parameter, paramVec, and initial values, stateVec, to the input 
# parameters for expID. In case the values experimentalConditions cannot be 
# be found in experimentalConditions file will search for the values 
# in the parameterBounds (parameter values) file. To ensure that the parameter 
# values are mapped correctly parameterNames = parameters(odeSys) and 
# stateNames = states(odeSys) is combined paramMap and stateMap (parameters in 
# ModellingToolkit format) to ensuare a correct mapping. 
# As the function currently rounds values, e.g does not allow for different 
# parameters for different conditions, a bit of tweaking is required to make 
# it compatible with parameter estimation. 
function changeToCond!(paramVec, 
                       stateVec, 
                       expID, 
                       parameterBounds::DataFrame,
                       experimentalConditions::DataFrame,
                       parameterNames, 
                       stateNames, 
                       paramMap, 
                       stateMap)

    # Set up simulation conditions (add to build SS-sim funciton)
    colNames = names(experimentalConditions)
    i_start = "conditionName" in colNames ? 3 : 2
    stateParamChange = colNames[i_start:end]
    if isempty(stateParamChange)
        return 
    end
    valsChangeTo = Vector(experimentalConditions[getRowExpId(expID, experimentalConditions), i_start:end])
    valsChangeTo = string.(valsChangeTo)

    parameterNamesStr = string.([paramMap[i].first for i in eachindex(paramMap)])
    stateNamesStr = replace.(string.([stateMap[i].first for i in eachindex(stateMap)]), "(t)" => "")

    # Change parameters (and states)
    for i in eachindex(stateParamChange)
        param = stateParamChange[i]

        # Get value to change to 
        valChangeTo::Float64 = 0.0
        if isNumber(valsChangeTo[i])
            valChangeTo = parse(Float64, valsChangeTo[i])
        elseif getRowExpId(valsChangeTo[i], parameterBounds, colSearch="parameterId") != nothing
            iRow = getRowExpId(valsChangeTo[i], parameterBounds, colSearch="parameterId") 
            valChangeTo = parameterBounds[iRow, "nominalValue"]
        else
            println("Error : Simulation parameter not found for experimental condition $expID")
            println("valsChangeTo[i] = ", valsChangeTo[i])
        end

        # Check for value to change to in parameter file 
        i_param = findfirst(x -> x == param, parameterNamesStr)
        i_state = findfirst(x -> x == param, stateNamesStr)

        if !isnothing(i_param)
            paramMap[i_param] = Pair(paramMap[i_param].first, valChangeTo) 
        elseif !isnothing(i_state)
            stateMap[i_state] = Pair(stateMap[i_state].first, valChangeTo)
        else
            println("Error : Simulation parameter to change not found for experimental condition $expID")
        end
    end

    n_states = length(stateNames)
    newVal = ModelingToolkit.varmap_to_vars(vcat(stateMap, paramMap), vcat(stateNames, parameterNames))
    stateVec .= newVal[1:n_states]
    paramVec .= newVal[(n_states+1):end]

    return nothing
end


# Given a model name, e.g model_Beer_MolBioSystems2014.jl, read the associated PeTab files 
# for the measurements, parameters and experimental conditions. 
function readDataFiles(modelName)

    # Reading data for model
    modelNameShort = modelName[7:end]
    readDataPath = joinpath(pwd(), "Pipeline_ModelParameterEstimation", "Data", modelName)
    dataEnding = modelNameShort * ".tsv"
    experimentalConditions = CSV.read(joinpath(readDataPath, "experimentalCondition_" * dataEnding), DataFrame)
    measurementData = CSV.read(joinpath(readDataPath, "measurementData_" * dataEnding), DataFrame)
    parameterData = CSV.read(joinpath(readDataPath, "parameters_" * dataEnding), DataFrame)

    return experimentalConditions, measurementData, parameterData
end


# For a model in the form of an Ode-system get simulation information. Specifcially 
# using the measurment data file wheter the model should be simulated to a steady-state 
# first (simulateSS), names of the model-parameters, name of the model stataes, the
# experimentalConditionId for the pre-equlibration conditions (if simulating to a steady state),
# and the experimental id:s for the experimental conditions where matching against data occurs.
# If the model is not simulated to a steady state shiftExp is empty, and the experimental 
# id:s are stored firstExpIds. 
function getSimulationInfo(measurementData::DataFrame, 
                           odeSystem::ODESystem)

    # Check if model first should be simulated to a steady state 
    colNames = names(measurementData)
    if !("preequilibrationConditionId" in colNames)
        preEqIDs = []
    else
        preEqIDs = unique(filter(x -> !ismissing(x), measurementData[!, "preequilibrationConditionId"]))
    end
    simulateSS = length(preEqIDs) > 0

    # If simulated to ss choose one of the preequilibrationConditionId as first ID 
    if simulateSS == true
        firstExpIds = preEqIDs
        shiftExpIds = Any[]
        for firstExpId in firstExpIds
            iRows = findall(x -> x == firstExpId, measurementData[!, "preequilibrationConditionId"])
            shiftExpId = unique(measurementData[iRows, "simulationConditionId"])
            push!(shiftExpIds, shiftExpId)
        end
    end

    if simulateSS == false
        firstExpIds = unique(measurementData[!, "simulationConditionId"])
        shiftExpIds = nothing
    end

    parameterNames = parameters(odeSystem)
    stateNames = states(odeSystem)

    return firstExpIds, shiftExpIds, simulateSS, parameterNames, stateNames
end


# For a simulationConditionId (expId) get the end-time to simulate to from the measurementData
function getTimeMax(measurementData::DataFrame, expId)
    return Float64(maximum(measurementData[findall(x -> x == expId, measurementData[!, "simulationConditionId"]), "time"]))
end


# Set the values in paramMap and stateMap to the nominal values in the parametersFile. 
# This allows constant values to be set correctly when mapping paramMap and stateMap 
# to the ODE-problem.  
function setParamToParamFileVal!(paramMap, stateMap, parameterBounds::DataFrame)

    parameterNames = parameterBounds[!, "parameterId"]
    parameterNamesStr = string.([paramMap[i].first for i in eachindex(paramMap)])
    stateNamesStr = replace.(string.([stateMap[i].first for i in eachindex(stateMap)]), "(t)" => "")
    for i in eachindex(parameterNames)
        
        parameterName = parameterNames[i]
        valChangeTo = parameterBounds[i, "nominalValue"]
        
        # Check for value to change to in parameter file 
        i_param = findfirst(x -> x == parameterName, parameterNamesStr)
        i_state = findfirst(x -> x == parameterName, stateNamesStr)

        if !isnothing(i_param)
            paramMap[i_param] = Pair(paramMap[i_param].first, valChangeTo) 
        elseif !isnothing(i_state)
            stateMap[i_state] = Pair(stateMap[i_state].first, valChangeTo)
        end
    end

end


# Solve and ODE-model (in prob) for all experimental conditions specifed by 
# first and shift experimental IDs using a given solver tolerance and 
# ode solver. Based on whether or not the model should be simulated to a 
# steady state it calls on, of two, possible solvers. changeToCondUse! ensures 
# that ODE-parameters are changed correctly. 
# Returns solArray, an array of the ode solutions across all conditions, 
# and succes::Bool. If the solver did not return Succesfull retcode 
# succes is returned as false. 
function solveOdeModelAllCond(prob::ODEProblem, 
                              changeToCondUse!::Function, 
                              simulateSS::Bool, 
                              measurementData::DataFrame,
                              firstExpIds, 
                              shiftExpIds, 
                              tol::Float64, 
                              solver)

    local solArray
    local sucess = true
    if simulateSS == true
        nShiftId = Int(sum([length(shiftExpIds[i]) for i in eachindex(shiftExpIds)]))
        solArray = Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1}(undef, nShiftId)
        k = 1

        for i in eachindex(firstExpIds)
            for j in eachindex(shiftExpIds[i])
                firstExpId = firstExpIds[i]
                shiftExpId = shiftExpIds[i][j]
                t_max_ss = getTimeMax(measurementData, shiftExpId)
                solArray[k] = solveOdeSS(prob, changeToCondUse!, firstExpId, shiftExpId, tol, t_max_ss, solver)
                
                if solArray[k].retcode != :Success
                    sucess = false
                    break 
                end

                k += 1
            end
        end

    elseif simulateSS == false
        nExperimentalCond = Int64(length(firstExpIds))
        solArray = Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1}(undef, nExperimentalCond)

        for i in eachindex(firstExpIds)
            firstExpId = firstExpIds[i]
            t_max = getTimeMax(measurementData, firstExpId)
            solArray[i] = solveOdeNoSS(prob, changeToCondUse!, firstExpId, tol, solver, t_max)

            if !(solArray[i].retcode == :Success || solArray[i].retcode == :Terminated)
                sucess = false
                break 
            end
        end
    end

    return solArray, sucess
end
# For an experimental condition solve the where first a pre-simulation is performed 
# (simulate to a stady state). If the pre-simulation fails, return failed solver. 
# Ultimately, if provided by the user the solution can be saved at the points tSave. 
function solveOdeSS(prob::ODEProblem, 
                      changeToCondUse!::Function, 
                      firstExpId, 
                      shiftExpId,
                      tol::Float64, 
                      t_max_ss,
                      solver;
                      tSave=Float64[])

    changeToCondUse!(prob.p, prob.u0, firstExpId)
    u0_pre = deepcopy(prob.u0)
    prob = remake(prob, tspan = (0.0, 1e6))
    sol_pre = solve(prob, solver, abstol=tol, reltol=tol, callback=TerminateSteadyState())
    # In case a steady state was not reached 
    if sol_pre.retcode != :Terminated
        return sol_pre
    end

    # Post ss 
    changeToCondUse!(prob.p, prob.u0, shiftExpId)
    # Make sure to add steady state vector for u0-val that has not changed with applaying shift-condition 
    has_not_changed = prob.u0 .== u0_pre
    prob.u0[has_not_changed] .= sol_pre.u[end][has_not_changed]
    
    prob = remake(prob, tspan = (0.0, t_max_ss))
    sol = solve(prob, solver, abstol=tol, reltol=tol, saveat=tSave)

    return sol
end
# For an experimental condition solve the where first a pre-simulation is not performed. 
# Ultimately, if provided by the user the solution can be saved at the points tSave. 
function solveOdeNoSS(prob, changeToCondUse!::Function, firstExpId, tol::Float64, solver, t_max::Float64; tSave=Float64[])
    
    changeToCondUse!(prob.p, prob.u0, firstExpId)
    prob = remake(prob, tspan=(0.0, t_max))
    t_max = prob.tspan[2]
    if isinf(t_max) 
        # Is case the model is simulated to and matched against steady state data (not standard)
        sol = solve(prob, solver, abstol=tol, reltol=tol, callback=TerminateSteadyState(), save_on=false, save_end=true)
    else
        sol = solve(prob, solver, abstol=tol, reltol=tol, saveat=tSave)
    end

    return sol
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

    
    local solArraySolver
    local couldSolve = true
    if simulateSS == true
        nShiftId = Int(sum([length(shiftExpIds[i]) for i in eachindex(shiftExpIds)]))
        solArraySolver = Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1}(undef, nShiftId)
        k = 1

        for i in eachindex(firstExpIds)
            for j in eachindex(shiftExpIds[i])
                firstExpId = firstExpIds[i]
                shiftExpId = shiftExpIds[i][j]
                t_max_ss = getTimeMax(measurementData, shiftExpId)
                solArraySolver[k] = solveOdeSS(prob, changeToCondUse!, firstExpId, shiftExpId, tol, t_max_ss, solver, tSave=solArray[k].t)
                
                if solArraySolver[k].retcode != :Success && solArraySolver[i].t[end] != solArray[i].t[end]
                    couldSolve = false
                    break 
                end

                k += 1
            end
        end

    elseif simulateSS == false
        nExperimentalCond = Int64(length(firstExpIds))
        solArraySolver = Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1}(undef, nExperimentalCond)

        for i in eachindex(firstExpIds)
            firstExpId = firstExpIds[i]
            t_max = getTimeMax(measurementData, firstExpId)
            solArraySolver[i] = solveOdeNoSS(prob, changeToCondUse!, firstExpId, tol, solver, t_max, tSave=solArray[i].t)

            if !(solArraySolver[i].retcode == :Success || solArraySolver[i].retcode == :Terminated) || solArraySolver[i].t[end] != solArray[i].t[end]
                couldSolve = false
                break 
            end
        end
    end

    if couldSolve == false
        return Inf
    end

    sqErr::Float64 = 0.0
    for i in eachindex(solArray)
        solI = solArraySolver[i]
        solIHigh = solArray[i]

        if length(solI.t) != length(solIHigh.t)
            # This likelly happens due to events in the solution 
            sqErrI = sum((Array(solI(solIHigh.t)) - Array(solIHigh)).^2)
            sqErrI -= sum((Array(solIHigh(solIHigh.t)) - Array(solIHigh)).^2) # Protect against interpolation errors with events 
        else
            sqErrI = sum((solI[:,:] - solIHigh[:,:]).^2)
            sqErrI -= sum((Array(solIHigh(solIHigh.t)) - Array(solIHigh)).^2) # Protect against interpolation errors with events 
        end

        sqErr += sqErrI
    end

    return sqErr
end
    