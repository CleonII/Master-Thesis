struct ParamData{T1<:Array{<:AbstractFloat}, 
                 T2<:Array{<:String, 1}, 
                 T3<:Array{Bool, 1}, 
                 T4<:Signed}

    paramVal::T1
    lowerBounds::T1
    upperBounds::T1
    parameterID::T2
    logScale::T3
    shouldEst::T3
    nParamEst::T4
end


struct ObservedData{T1<:Array{<:AbstractFloat, 1}, 
                    T2<:Array{<:String, 1}, 
                    T3<:Array{<:Symbol, 1}}

    yObs::T1
    noiseObs::T1
    tObs::T1
    observebleDd::T2
    conditionId::T2  # Sum of pre-eq + simulation-cond id 
    sdParamId::T2
    transformData::T3 # Only done once 
end


function processObsData(measurementData::DataFrame)::ObservedData

    # Process measurment data 
    yObs::Array{Float64, 1} = convert(Array{Float64, 1}, measurementData[!, "measurement"])
    tObs::Array{Float64, 1} = convert(Array{Float64, 1}, measurementData[!, "time"])
    nObs = length(yObs)

    conditionId::Array{String, 1} = Array{String, 1}(undef, nObs)
    preEq = measurementData[!, "preequilibrationConditionId"]
    simCond = measurementData[!, "simulationConditionId"]
    for i in eachindex(conditionId)
        if ismissing(preEq[i])
            conditionId[i] = String(simCond[i])
        else
            conditionId[i] = String(preEq[i]) * String(simCond[i])
        end
    end

    obsID::Array{String, 1} = String.(measurementData[!, "observableId"])

    # Noise parameters are either mapped to have a parameter ID, or set directly to their value if known.
    sdObs = String.(measurementData[!, "noiseParameters"])
    sdId::Array{String, 1} = Array{String, 1}(undef, nObs)
    sdVal::Array{Float64, 1} = Array{Float64, 1}(undef, nObs)
    for i in eachindex(sdObs)
        if isNumber(sdObs[i])
            sdVal[i] = parse(Float64, sdObs[i])
            sdId[i] = ""
        else
            sdVal[i] = 0.0
            sdId[i] = sdObs[i]
        end
    end

    # Whether or not data should be transformed or not 
    if !("observableTransformation" in names(measurementData))
        transformArr = [:lin for i in 1:nObs]
    else
        transformArr = [Symbol(measurementData[i, "observableTransformation"]) for i in 1:nObs]
    end

    # Transform measurment data if given 
    transformObsOrDataArr!(yObs, transformArr)

    return ObservedData(yObs, sdVal, tObs, obsID, conditionId, sdId, transformArr)
end

# Create paramData from parameter data. Putting everything into a struct 
# will make subsequent simulations more efficient. 
function processParameterData(parameterData::DataFrame)

    nParam = length(parameterData[!, "estimate"])
    lb::Array{Float64, 1} = zeros(Float64, nParam) 
    ub::Array{Float64, 1} = zeros(Float64, nParam)
    paramVal::Array{Float64, 1} = zeros(Float64, nParam)
    logScale::Array{Bool, 1} = Array{Bool, 1}(undef, nParam)
    paramId::Array{String, 1} = Array{String, 1}(undef, nParam)
    shouldEst::Array{Bool, 1} = Array{Bool, 1}(undef, nParam)

    for i in eachindex(shouldEst)
        if ismissing(parameterData[i, "lowerBound"])
            lb[i] = -Inf
        else
            lb[i] = parameterData[i, "lowerBound"]
        end
        if ismissing(parameterData[i, "upperBound"])
            ub[i] = Inf
        else
            ub[i] = parameterData[i, "upperBound"]
        end
        paramVal[i] = parameterData[i, "nominalValue"]
        paramId[i] = parameterData[i, "parameterId"]
        logScale[i] = parameterData[i, "parameterScale"] == "log10" ? true : false
        shouldEst[i] = parameterData[i, "estimate"] == 1 ? true : false
    end
    nParamEst::Int = Int(sum(shouldEst))

    return ParamData(paramVal, lb, ub, paramId, logScale, shouldEst, nParamEst)
end


function transformObsOrDataArr!(vals, transformArr::Array{Symbol, 1})
    for i in eachindex(vals)
        vals[i] -= (vals[i] - transformObsOrData(vals[i], transformArr[i]))
    end
end


function transformObsOrData(val, transform::Symbol)
    if transform == :lin
        return val
    elseif transform == :log10
        return log10(val)
    elseif transform == :log
        return log(val)
    else
        println("Error : Not a legitimate transform")
    end
end


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
# in the paramData (parameter values) struct. To ensure that the parameter 
# values are mapped correctly parameterNames = parameters(odeSys) and 
# stateNames = states(odeSys) is combined paramMap and stateMap (parameters in 
# ModellingToolkit format) to ensuare a correct mapping. 
# As the function currently rounds values, e.g does not allow for different 
# parameters for different conditions, a bit of tweaking is required to make 
# it compatible with parameter estimation. 
function changeToCond!(paramVec, 
                       stateVec, 
                       expID, 
                       paramData::ParamData,
                       experimentalConditions::DataFrame,
                       parameterNames, 
                       stateNames, 
                       paramMap, 
                       stateMap)

    # Allow the code to propegate dual numbers for gradients 
    paramMapUse = convert.(Pair{Num, eltype(paramVec)}, paramMap)
    stateMapUse = convert.(Pair{Num, Num}, stateMap)

    # Set up simulation conditions (add to build SS-sim funciton)
    colNames = names(experimentalConditions)
    i_start = "conditionName" in colNames ? 3 : 2
    stateParamChange = colNames[i_start:end]
    if isempty(stateParamChange)
        return 
    end
    valsChangeTo = Vector(experimentalConditions[getRowExpId(expID, experimentalConditions), i_start:end])
    valsChangeTo = string.(valsChangeTo)

    parameterNamesStr = string.([paramMapUse[i].first for i in eachindex(paramMapUse)])
    stateNamesStr = replace.(string.([stateMapUse[i].first for i in eachindex(stateMapUse)]), "(t)" => "")

    # Change parameters (and states)
    for i in eachindex(stateParamChange)
        param = stateParamChange[i]

        # Get value to change to 
        valChangeTo::Float64 = 0.0
        if isNumber(valsChangeTo[i])
            valChangeTo = parse(Float64, valsChangeTo[i])
        elseif findfirst(x -> x == valsChangeTo[i], paramData.parameterID) != nothing
            iVal = findfirst(x -> x == valsChangeTo[i], paramData.parameterID)
            valChangeTo = paramData.paramVal[iVal]
        else
            println("Error : Simulation parameter not found for experimental condition $expID")
            println("valsChangeTo[i] = ", valsChangeTo[i])
        end

        # Check for value to change to in parameter file 
        i_param = findfirst(x -> x == param, parameterNamesStr)
        i_state = findfirst(x -> x == param, stateNamesStr)

        if !isnothing(i_param)
            paramMapUse[i_param] = Pair(paramMapUse[i_param].first, valChangeTo) 
        elseif !isnothing(i_state)
            stateMapUse[i_state] = Pair(stateMapUse[i_state].first, valChangeTo)
        else
            println("Error : Simulation parameter to change not found for experimental condition $expID")
        end
    end

    n_states = length(stateNames)
    newVal = ModelingToolkit.varmap_to_vars(vcat(stateMapUse, paramMapUse), vcat(stateNames, parameterNames))
    stateVec .-= (stateVec - newVal[1:n_states])
    paramVec .-= (paramVec - newVal[(n_states+1):end])

    return nothing
end


# Change the model parameters that should be estimated 
function changeModelParam!(paramVec, 
                           stateVec,
                           paramEstVec,
                           paramIdChange::Array{String, 1},
                           paramData::ParamData,
                           parameterNames, 
                           paramMap)

    # Allow the code to propegate dual numbers for gradients 
    paramMapUse = convert.(Pair{Num, eltype(paramVec)}, paramMap)
    parameterNamesStr = string.([paramMapUse[i].first for i in eachindex(paramMapUse)])

    # Change parameters (and states)
    for i in eachindex(paramIdChange)
        param = paramIdChange[i]

        # Get value to change to 
        valChangeTo= 0.0
        
        if findfirst(x -> x == paramIdChange[i], paramData.parameterID) != nothing
            iVal = findfirst(x -> x == paramIdChange[i], paramData.parameterID)
            valChangeTo = paramEstVec[iVal]
        else
            println("Error : Could not find the correct model parameter to change")
            println("paramIdChange[i] = ", paramIdChange[i])
        end

        # Check for value to change to in parameter file 
        i_param = findfirst(x -> x == param, parameterNamesStr)

        if !isnothing(i_param)
            paramMapUse[i_param] = Pair(paramMapUse[i_param].first, valChangeTo) 
        else
            println("Error : Simulation parameter to change not found for experimental condition $expID")
        end
    end

    newVal = ModelingToolkit.varmap_to_vars(paramMapUse, parameterNames)
    paramVec .= newVal

    # Update stateVec according to mapping in availble file 
    Boehm_JProteomeRes2014_t0!(stateVec, paramVec) 
    
    return nothing
end
function changeToCondEst!(paramVec, 
                          stateVec, 
                          expID, 
                          paramData::ParamData,
                          experimentalConditions::DataFrame,
                          parameterNames, 
                          paramMap)

    # Allow the code to propegate dual numbers for gradients 
    paramMapUse = convert.(Pair{Num, eltype(paramVec)}, paramMap)

    # Set up simulation conditions (add to build SS-sim funciton)
    colNames = names(experimentalConditions)
    i_start = "conditionName" in colNames ? 3 : 2
    stateParamChange = colNames[i_start:end]
    if isempty(stateParamChange)
        return 
    end
    valsChangeTo = Vector(experimentalConditions[getRowExpId(expID, experimentalConditions), i_start:end])
    valsChangeTo = string.(valsChangeTo)

    parameterNamesStr = string.([paramMapUse[i].first for i in eachindex(paramMapUse)])

    # Change parameters (and states)
    for i in eachindex(stateParamChange)
        param = stateParamChange[i]

        # Get value to change to 
        valChangeTo::Float64 = 0.0
        if isNumber(valsChangeTo[i])
            valChangeTo = parse(Float64, valsChangeTo[i])
        elseif findfirst(x -> x == valsChangeTo[i], paramData.parameterID) != nothing
            iVal = findfirst(x -> x == valsChangeTo[i], paramData.parameterID)
            valChangeTo = paramData.paramVal[iVal]
        else
            println("Error : Simulation parameter not found for experimental condition $expID")
            println("valsChangeTo[i] = ", valsChangeTo[i])
        end

        # Check for value to change to in parameter file 
        i_param = findfirst(x -> x == param, parameterNamesStr)

        if !isnothing(i_param)
            paramMapUse[i_param] = Pair(paramMapUse[i_param].first, valChangeTo) 
        else
            println("Error : Simulation parameter to change not found for experimental condition $expID")
        end
    end

    newVal = ModelingToolkit.varmap_to_vars(paramMapUse, parameterNames)
    paramVec .= newVal
    Boehm_JProteomeRes2014_t0!(stateVec, paramVec) 

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
function setParamToParamFileVal!(paramMap, stateMap, paramData::ParamData)

    parameterNames = paramData.parameterID
    parameterNamesStr = string.([paramMap[i].first for i in eachindex(paramMap)])
    stateNamesStr = replace.(string.([stateMap[i].first for i in eachindex(stateMap)]), "(t)" => "")
    for i in eachindex(parameterNames)
        
        parameterName = parameterNames[i]
        valChangeTo = paramData.paramVal[i]
        
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
                              solver;
                              nTSave::Int64=0, 
                              denseArg=true)

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
                solArray[k] = solveOdeSS(prob, changeToCondUse!, firstExpId, shiftExpId, tol, t_max_ss, solver, nTSave=nTSave, denseArg=denseArg)
                
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
            solArray[i] = solveOdeNoSS(prob, changeToCondUse!, firstExpId, tol, solver, t_max, nTSave=nTSave, denseArg=denseArg)

            if !(solArray[i].retcode == :Success || solArray[i].retcode == :Terminated)
                sucess = false
                break 
            end
        end
    end

    return solArray, sucess
end
function solveOdeModelAllCond!(solArray::Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1},
                              conditionId::Array{String, 1},
                              prob::ODEProblem, 
                              changeToCondUse!::Function, 
                              simulateSS::Bool, 
                              measurementData::DataFrame,
                              firstExpIds, 
                              shiftExpIds, 
                              tol::Float64, 
                              solver;
                              nTSave::Int64=0, 
                              denseArg=true)

    local sucess = true
    if simulateSS == true
        k = 1
        for i in eachindex(firstExpIds)
            for j in eachindex(shiftExpIds[i])
                firstExpId = firstExpIds[i]
                shiftExpId = shiftExpIds[i][j]
                t_max_ss = getTimeMax(measurementData, shiftExpId)
                solArray[k] = solveOdeSS(prob, changeToCondUse!, firstExpId, shiftExpId, tol, t_max_ss, solver, nTSave=nTSave, denseArg=denseArg)
                
                if solArray[k].retcode != :Success
                    sucess = false
                    break 
                end

                conditionId[k] = firstExpId * shiftExpId

                k += 1
            end
        end

    elseif simulateSS == false
        for i in eachindex(firstExpIds)
            firstExpId = firstExpIds[i]
            t_max = getTimeMax(measurementData, firstExpId)
            solArray[i] = solveOdeNoSS(prob, changeToCondUse!, firstExpId, tol, solver, t_max, nTSave=nTSave, denseArg=denseArg)

            conditionId[i] = firstExpId

            if !(solArray[i].retcode == :Success || solArray[i].retcode == :Terminated)
                sucess = false
                break 
            end
        end
    end

    return sucess
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
                      tSave=Float64[], 
                      nTSave=0, 
                      denseArg=true)

    # Check that for no clash between number of data-points to save tSave (only one can be activate)    
    if length(tSave) != 0 && nTSave != 0
        println("Error : Can only provide tSave (vector to save at) or nTSave as saveat argument to solvers")
    elseif nTSave != 0
        saveAtVec = collect(LinRange(0.0, t_max_ss, nTSave))
    else
        saveAtVec = tSave
    end

    # Whether or not producing a dense soluation 
    if (isempty(tSave) && nTSave == 0) && denseArg == true
        dense = true
    else
        dense = false
    end
    
    # Check wheter a solver or alg-hint is provided 
    if typeof(solver) <: Vector{Symbol} # In case an Alg-hint is provided 
        solveCallPre = (prob) -> solve(prob, alg_hints=solver, abstol=tol, reltol=tol, callback=TerminateSteadyState(), dense=false)
        solveCallPost = (prob) -> solve(prob, alg_hints=solver, abstol=tol, reltol=tol, tstops=saveAtVec, saveat=saveAtVec, dense=dense)
    else
        solveCallPre = (prob) -> solve(prob, solver, abstol=tol, reltol=tol, callback=TerminateSteadyState(), dense=false)
        solveCallPost = (prob) -> solve(prob, solver, abstol=tol, reltol=tol, tstops=saveAtVec, saveat=saveAtVec, dense=dense)
    end

    changeToCondUse!(prob.p, prob.u0, firstExpId)
    u0_pre = deepcopy(prob.u0)
    prob = remake(prob, tspan = (0.0, 1e6))
    sol_pre = solveCallPre(prob)
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
    sol = solveCallPost(prob)

    return sol
end
# For an experimental condition solve the where first a pre-simulation is not performed. 
# Ultimately, if provided by the user the solution can be saved at the points tSave. 
function solveOdeNoSS(prob, changeToCondUse!::Function, firstExpId, tol::Float64, solver, t_max::Float64; tSave=Float64[], nTSave::Int64=0, denseArg=denseArg)

    changeToCondUse!(prob.p, prob.u0, firstExpId)
    prob = remake(prob, tspan=(0.0, t_max))
    t_max = prob.tspan[2]

    # Check that for no clash between number of data-points to save tSave (only one can be activate)    
    if length(tSave) != 0 && nTSave != 0
        println("Error : Can only provide tSave (vector to save at) or nTSave as saveat argument to solvers")
    elseif nTSave != 0
        saveAtVec = collect(LinRange(0.0, t_max, nTSave))
    else
        saveAtVec = tSave
    end

    # Whether or not producing a dense soluation 
    # Whether or not producing a dense soluation 
    if (isempty(tSave) && nTSave == 0) && denseArg == true
        dense = true
    else
        dense = false
    end

    if typeof(solver) <: Vector{Symbol} && isinf(t_max)
        solveCall = (prob) -> solve(prob, alg_hints=solver, abstol=tol, reltol=tol, callback=TerminateSteadyState(), save_on=false, save_end=true, dense=dense)
    elseif typeof(solver) <: Vector{Symbol} && !isinf(t_max)
        solveCall = (prob) -> solve(prob, alg_hints=solver, abstol=tol, reltol=tol, tstops=saveAtVec, saveat=saveAtVec, dense=dense)
    elseif !(typeof(solver) <: Vector{Symbol}) && isinf(t_max)
        solveCall = (prob) -> solve(prob, solver, abstol=tol, reltol=tol, callback=TerminateSteadyState(), save_on=false, save_end=true, dense=dense)
    elseif !(typeof(solver) <: Vector{Symbol}) && !isinf(t_max)
        solveCall = (prob) -> solve(prob, solver, abstol=tol, reltol=tol, tstops=saveAtVec, saveat=saveAtVec, dense=dense)
    else
        println("Error : Solver option does not exist")        
    end
    
    sol = solveCall(prob)
    return sol
end


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
    GC.gc(); GC.gc(); GC.gc(); GC.gc()
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
                GC.gc(); GC.gc(); GC.gc()
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
            if !isinf(t_max) && !(solCompare[i].retcode == :Success && solCompare[i].t[end] == solArray[i].t[end])
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
    