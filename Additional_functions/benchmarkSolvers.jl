struct PeTabModel{T1<:Vector{<:Pair{Num, <:Union{AbstractFloat, Num}}}, 
                  T2<:Vector{<:Pair{Num, <:Union{AbstractFloat, Num}}},
                  T3<:Vector{Sym{Real, Base.ImmutableDict{DataType, Any}}}, 
                  T4<:Vector{Any}}
    modelName::String
    evalYmod::Function 
    evalU0!::Function
    evalSd!::Function
    odeSystem::ODESystem 
    paramMap::T1
    stateMap::T2
    paramNames::T3
    stateNames::T4
    dirModel::String
    pathMeasurementData::String
    pathExperimentalConditions::String
    pathObservables::String
    pathParameters::String
end


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
    tObs::T1
    observebleDd::T2
    conditionId::T2  # Sum of pre-eq + simulation-cond id 
    sdParams::T2
    transformData::T3 # Only done once 
    obsParam::T2
end

struct SimulationInfo{T1<:Array{<:String, 1}, 
                      T2<:Bool,
                      T3<:Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1}}
    firstExpIds::T1
    shiftExpIds::T1
    conditionIdSol::T1
    simulateSS::T2
    solArray::T3
end


function processObsData(measurementData::DataFrame, observableData::DataFrame)::ObservedData

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
    sdParams::Array{String, 1} = Array{String, 1}(undef, nObs)
    for i in eachindex(sdObs)
        if ismissing(sdObs[i])
            sdParams[i] = ""
        else
            sdParams[i] = sdObs[i]
        end
    end

    # Handle observeble parameters (e.g scale, offset parmeters)
    obsParamFile = measurementData[!, "observableParameters"]
    obsParam = Array{String, 1}(undef, nObs)
    for i in 1:nObs
        if ismissing(obsParamFile[i])
            obsParam[i] = ""
        else
            obsParam[i] = obsParamFile[i]
        end
    end

    # Whether or not data should be transformed or not 
    transformArr = Array{Symbol, 1}(undef, nObs)
    if !("observableTransformation" in names(observableData))
        transformArr .= [:lin for i in 1:nObs]
    else
        for i in 1:nObs
            iRow = findfirst(x -> x == obsID[i], observableData[!, "observableId"])
            transformArr[i] = Symbol(observableData[iRow, "observableTransformation"]) 
        end
    end

    # Transform measurment data if given 
    transformObsOrDataArr!(yObs, transformArr)

    return ObservedData(yObs, tObs, obsID, conditionId, sdParams, transformArr, obsParam)
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
function isNumber(a::SubString{String})
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
function changeExperimentalCond!(paramVec, 
                                 stateVec, 
                                 expID, 
                                 parameterData::ParamData,
                                 experimentalConditions::DataFrame,
                                 peTabModel::PeTabModel)

    # Allow the code to propegate dual numbers for gradients 
    paramMapUse = convert.(Pair{Num, eltype(paramVec)}, peTabModel.paramMap)

    # Set up simulation conditions (add to build SS-sim funciton)
    # TODO : Process experimental conditions into Julia struct 
    colNames = names(experimentalConditions)
    i_start = "conditionName" in colNames ? 3 : 2
    stateParamChange = colNames[i_start:end]
    if isempty(stateParamChange)
        return 
    end
    valsChangeTo = Vector(experimentalConditions[getRowExpId(expID, experimentalConditions), i_start:end])
    valsChangeTo = string.(valsChangeTo)

    # Parameter and state names in string-format for easier comparisons 
    parameterNamesStr = string.([paramMapUse[i].first for i in eachindex(paramMapUse)])
    
    # Change parameters (and states)
    iParamChange = Array{Int64, 1}(undef, length(stateParamChange))
    for i in eachindex(stateParamChange)
        param = stateParamChange[i]
        iParamChange[i] = findfirst(x -> x == param, string.(peTabModel.paramNames))

        # Get value to change to 
        valChangeTo::Float64 = 0.0
        if isNumber(valsChangeTo[i])
            valChangeTo = parse(Float64, valsChangeTo[i])
        elseif findfirst(x -> x == valsChangeTo[i], parameterData.parameterID) != nothing
            iVal = findfirst(x -> x == valsChangeTo[i], parameterData.parameterID)
            valChangeTo = parameterData.paramVal[iVal]
        else
            println("Error : Simulation parameter not found for experimental condition $expID")
            println("valsChangeTo[i] = ", valsChangeTo[i])
        end

        # Propegate dual numbers correctly. In the reference paramMap want to have floats (not duals)
        valChangeToFloat = typeof(valChangeTo) <: Union{ForwardDiff.Dual, ForwardDiff.Dual{<:ForwardDiff.Dual}} ? valChangeTo.value : valChangeTo

        # Check for value to change to in parameter file 
        i_param = findfirst(x -> x == param, parameterNamesStr)

        if !isnothing(i_param)
            paramMapUse[i_param] = Pair(paramMapUse[i_param].first, valChangeTo) 
            peTabModel.paramMap[i_param] = Pair(peTabModel.paramMap[i_param].first, valChangeToFloat) 
        else
            println("Error : Simulation parameter to change not found for experimental condition $expID")
        end
    end

    # To avoid a reset of parameter values only change the parameters that gouvern the behaviour of an experimental state. 
    newVal = ModelingToolkit.varmap_to_vars(paramMapUse, peTabModel.paramNames)
    paramVec[iParamChange] .= newVal[iParamChange]
    peTabModel.evalU0!(stateVec, paramVec) 

    return nothing
end


# Change the model parameters that should be estimated 
function changeModelParam!(paramVec, 
                           stateVec,
                           paramEst,
                           paramEstNames::Array{String, 1},
                           paramData::ParamData,
                           peTabModel::PeTabModel)

    # Allow the code to propegate dual numbers for gradients 
    paramMapUse = convert.(Pair{Num, eltype(paramVec)}, peTabModel.paramMap)
    # TODO: Precompute this step 
    parameterNamesStr = string.([paramMapUse[i].first for i in eachindex(paramMapUse)])

    # Change parameters (and states) to current iterations in parameter estimation vector 
    for i in eachindex(paramEstNames)
        
        paramChangeName = paramEstNames[i]
        valChangeTo = paramEst[i]
        # Propegate dual numbers correctly. In the reference paramMap want to have floats (not duals)
        valChangeToFloat = typeof(valChangeTo) <: Union{ForwardDiff.Dual, ForwardDiff.Dual{<:ForwardDiff.Dual}} ? valChangeTo.value : valChangeTo
        
        i_param = findfirst(x -> x == paramChangeName, parameterNamesStr)
        if !isnothing(i_param)
            paramMapUse[i_param] = Pair(paramMapUse[i_param].first, valChangeTo) 
            # Update reference paramMap 
            if typeof(valChangeToFloat) <: AbstractFloat
                peTabModel.paramMap[i_param] = Pair(peTabModel.paramMap[i_param].first, valChangeToFloat) 
            end
        else
            println("Error : Simulation parameter to change not found for experimental condition $expID")
        end
    end

    # Use ModellingToolkit and u0 function to correctly map parameters to ODE-system 
    newVal = ModelingToolkit.varmap_to_vars(paramMapUse, peTabModel.paramNames)
    paramVec .= newVal
    peTabModel.evalU0!(stateVec, paramVec) 
    
    return nothing
end


# Given a model name, e.g model_Beer_MolBioSystems2014.jl, read the associated PeTab files 
# for the measurements, parameters and experimental conditions. 
function readDataFiles(dirModel; readObs::Bool=false)

    # Extract PeTab-files 
    pathMeasurementData = checkForPeTabFile("measurementData", dirModel)
    pathExperimentalCond = checkForPeTabFile("experimentalCondition", dirModel)
    pathParameters = checkForPeTabFile("parameters", dirModel)
    pathObservables = checkForPeTabFile("observables", dirModel)

    # Reading data for model
    experimentalConditions = CSV.read(pathExperimentalCond, DataFrame)
    measurementData = CSV.read(pathMeasurementData, DataFrame)
    parameterData = CSV.read(pathParameters, DataFrame)

    if readObs == true
        observableData = CSV.read(pathObservables, DataFrame)
        return experimentalConditions, measurementData, parameterData, observableData
    else
        return experimentalConditions, measurementData, parameterData
    end
end


# For a model in the form of an Ode-system get simulation information. Specifcially 
# using the measurment data file wheter the model should be simulated to a steady-state 
# first (simulateSS), names of the model-parameters, name of the model stataes, the
# experimentalConditionId for the pre-equlibration conditions (if simulating to a steady state),
# and the experimental id:s for the experimental conditions where matching against data occurs.
# If the model is not simulated to a steady state shiftExp is empty, and the experimental 
# id:s are stored firstExpIds. 
function getSimulationInfo(measurementData::DataFrame, 
                           odeSystem::ODESystem)::SimulationInfo

    # Check if model first should be simulated to a steady state 
    colNames = names(measurementData)
    if !("preequilibrationConditionId" in colNames)
        preEqIDs = Array{String, 1}(undef, 0)
    else
        preEqIDs = convert(Array{String, 1}, unique(filter(x -> !ismissing(x), measurementData[!, "preequilibrationConditionId"])))
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
        shiftExpIds = convert(Array{String, 1}, shiftExpIds)
    end

    if simulateSS == false
        firstExpIds = convert(Array{String, 1}, unique(measurementData[!, "simulationConditionId"]))
        shiftExpIds = Array{String, 1}(undef, 0)
    end

    # Obtain number of forward simulations to perform 
    if simulateSS == true
        nForwardSol = Int64(sum([length(shiftExpIds[i]) for i in eachindex(shiftExpIds)]))
    else
        nForwardSol = Int64(length(firstExpIds))
    end
    # Array with conition-ID for a specific forward simulation 
    conditionIdSol = Array{String, 1}(undef, nForwardSol)
    # Foward ODE-solutions array 
    solArray = Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1}(undef, nForwardSol)

    simulationInfo = SimulationInfo(firstExpIds, 
                                    shiftExpIds, 
                                    conditionIdSol, 
                                    simulateSS,
                                    solArray)
    return simulationInfo
end


# For a simulationConditionId (expId) get the end-time to simulate to from the measurementData
function getTimeMax(measurementData::DataFrame, expId)::Float64
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
                              changeToExperimentalCondUse!::Function, 
                              measurementData::DataFrame,
                              simulationInfo::SimulationInfo,
                              solver, 
                              tol::Float64;
                              nTSave::Int64=0, 
                              denseArg=true)

    local solArray
    if simulationInfo.simulateSS == true
        nShiftId = Int(sum([length(simulationInfo.shiftExpIds[i]) for i in eachindex(simulationInfo.shiftExpIds)]))
        solArray = Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1}(undef, nShiftId)
    else
        nExperimentalCond = Int64(length(simulationInfo.firstExpIds))
        solArray = Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1}(undef, nExperimentalCond)
    end
        
    success = solveOdeModelAllCond!(solArray, 
                                    prob, 
                                    changeToExperimentalCondUse!, 
                                    measurementData, 
                                    simulationInfo, 
                                    solver, 
                                    tol, 
                                    nTSave=nTSave, 
                                    denseArg=denseArg)

    return solArray, success
end
function solveOdeModelAllCond!(solArray::Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1},
                               prob::ODEProblem, 
                               changeToExperimentalCondUse!::Function, 
                               measurementData::DataFrame,
                               simulationInfo::SimulationInfo,
                               solver,
                               tol::Float64;
                               nTSave::Int64=0, 
                               denseArg=true)

    local sucess = true
    if simulationInfo.simulateSS == true
        k = 1
        for i in eachindex(simulationInfo.firstExpIds)
            for j in eachindex(simulationInfo.shiftExpIds[i])
                firstExpId = simulationInfo.firstExpIds[i]
                shiftExpId = simulationInfo.shiftExpIds[i][j]
                t_max_ss = getTimeMax(measurementData, shiftExpId)
                solArray[k] = solveOdeSS(prob, changeToExperimentalCondUse!, firstExpId, shiftExpId, tol, t_max_ss, solver, nTSave=nTSave, denseArg=denseArg)
                
                if solArray[k].retcode != :Success
                    sucess = false
                    break 
                end

                simulationInfo.conditionIdSol[k] = firstExpId * shiftExpId

                k += 1
            end
        end

    elseif simulationInfo.simulateSS == false
        for i in eachindex(simulationInfo.firstExpIds)
            firstExpId = simulationInfo.firstExpIds[i]
            t_max = getTimeMax(measurementData, firstExpId)
            solArray[i] = solveOdeNoSS(prob, changeToExperimentalCondUse!, firstExpId, tol, solver, t_max, nTSave=nTSave, denseArg=denseArg)

            simulationInfo.conditionIdSol[i] = firstExpId

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
        solveCallPost = (prob) -> solve(prob, alg_hints=solver, abstol=tol, reltol=tol, saveat=saveAtVec, dense=dense)
    else
        solveCallPre = (prob) -> solve(prob, solver, abstol=tol, reltol=tol, callback=TerminateSteadyState(), dense=false)
        solveCallPost = (prob) -> solve(prob, solver, abstol=tol, reltol=tol, saveat=saveAtVec, dense=dense)
    end

    changeToCondUse!(prob.p, prob.u0, firstExpId)
    u0_pre = deepcopy(prob.u0)
    prob = remake(prob, tspan = (0.0, 1e6), p = prob.p[:], u0 = prob.u0[:])
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
function solveOdeNoSS(prob, changeToCondUse!::Function, firstExpId, tol::Float64, solver, t_max::Float64; tSave=Float64[], nTSave::Int64=0, denseArg=true)

    changeToCondUse!(prob.p, prob.u0, firstExpId)
    probUse = remake(prob, tspan=(0.0, t_max), u0 = prob.u0[:], p = prob.p[:])
    t_max = probUse.tspan[2]

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
        solveCall = (probArg) -> solve(probArg, alg_hints=solver, abstol=tol, reltol=tol, callback=TerminateSteadyState(), save_on=false, save_end=true, dense=dense)
    elseif typeof(solver) <: Vector{Symbol} && !isinf(t_max)
        solveCall = (probArg) -> solve(probArg, alg_hints=solver, abstol=tol, reltol=tol, saveat=saveAtVec, dense=dense)
    elseif !(typeof(solver) <: Vector{Symbol}) && isinf(t_max)
        solveCall = (probArg) -> solve(probArg, solver, abstol=tol, reltol=tol, callback=TerminateSteadyState(), save_on=false, save_end=true, dense=dense)
    elseif !(typeof(solver) <: Vector{Symbol}) && !isinf(t_max)
        solveCall = (probArg) -> solve(probArg, solver, abstol=tol, reltol=tol, saveat=saveAtVec, dense=dense)
    else
        println("Error : Solver option does not exist")        
    end

    sol = solveCall(probUse)

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
    

# For a model extract the relevant observeble parameters or sd-parameter for a given 
# observation-id, simulation-condition and time-point.
function getObsOrSdParam(vecParam,
                         paramData::ParamData, 
                         obsData::ObservedData, 
                         observableId::String,
                         simulationConditionId::String,
                         t::Float64;
                         getObsPar::Bool=true)

    if getObsPar == true
        idVec = getIdEst(obsData.obsParam, paramData)
    else
        idVec = getIdEst(obsData.sdParams, paramData)
    end

    nObs = length(obsData.tObs)
    # Acquire which observation is being handled (for really large parameter files smarter data structs might be needed)
    iUse = findfirst(x -> obsData.observebleDd[x] == observableId && obsData.conditionId[x] == simulationConditionId && obsData.tObs[x] == t, 1:nObs)
    if isnothing(iUse)
        println("Warning : Cannot identify an observation in getObsParam ")
    end
    if getObsPar == true && isempty(obsData.obsParam[iUse])
        return nothing
    elseif getObsPar == false && isempty(obsData.sdParams[iUse])
        return nothing
    end

    if getObsPar == true
        paramsRet = split(obsData.obsParam[iUse], ';')
    else
        paramsRet = split(obsData.sdParams[iUse], ';')
    end

    output = Array{eltype(vecParam), 1}(undef, length(paramsRet))
    for i in eachindex(paramsRet)
        # Hard coded constant number 
        if isNumber(paramsRet[i])
            output[i] = parse(Float64, paramsRet[i])

        # Parameter to be estimated 
        elseif paramsRet[i] in idVec
            output[i] = vecParam[findfirst(x -> x == paramsRet[i], idVec)]
        
        # Constant parameter 
        elseif paramsRet[i] in paramData.parameterID
            output[i] = paramData.paramVal[findfirst(x -> x == paramsRet[i], paramData.parameterID)]
        else
            println("Warning : cannot find matching for ", paramsRet[i])
        end
    end

    if length(output) == 1
        return output[1]
    else
        return output
    end
end


# Get the id of the obs-param which should be estimated (not treated as constants)
function getIdEst(idsInStr::Array{String, 1}, paramData::ParamData)
    idVec = String[]

    for i in eachindex(idsInStr)
        if isempty(idsInStr[i])
            continue
        else
            idsInStrSplit = split(idsInStr[i], ';')
            for idInStr in idsInStrSplit

                # Disregard Id if parameters should not be estimated, or 
                iParam = findfirst(x -> x == idInStr, paramData.parameterID)
                if isNumber(idInStr)
                    continue
                elseif isnothing(iParam)                    
                    println("Warning : param $idInStr could not be found in parameter file")
                elseif idInStr in idVec
                    continue
                elseif paramData.shouldEst[iParam] == false
                    continue
                else
                    idVec = vcat(idVec, idInStr)
                end
            end
        end
    end

    return idVec
end


struct ParameterIndices{T1<:Array{<:Integer, 1}, 
                        T2<:Array{<:String, 1}}

    iDynParam::T1
    iObsParam::T1
    iSdParam::T1
    namesDynParam::T2
    namesObsParam::T2
    namesSdParam::T2
    namesParamEst::T2
end


# Create indices for dynamic, observed and sd parameters 
function getIndicesParam(paramData::ParamData, obsData::ObservedData)::ParameterIndices

    idObsParam = getIdEst(obsData.obsParam, paramData)
    isObsParam = [paramData.parameterID[i] in idObsParam for i in eachindex(paramData.parameterID)]

    # Set up to bookeeping of indices between observed and noise parameters in the large input vector 
    idSdParam = getIdEst(obsData.sdParams, paramData)
    
    isSd = [paramData.parameterID[i] in idSdParam for i in eachindex(paramData.parameterID)]

    isDynamic = (paramData.shouldEst .&& .!isSd .&& .!isObsParam)
    idParamDyn = paramData.parameterID[isDynamic]
    
    # Index vector for the dynamic and sd parameters as UInt32 vectors 
    iDynPar::Array{UInt32, 1} = convert(Array{UInt32, 1}, collect(1:length(idParamDyn)))
    iSdPar::Array{UInt32, 1} = convert(Array{UInt32, 1}, collect((length(idParamDyn)+1):(length(idParamDyn) + length(idSdParam))))
    iObsPar::Array{UInt32, 1} = convert(Array{UInt32, 1}, collect((length(idParamDyn) + length(idSdParam) + 1):(length(idParamDyn) + length(idSdParam) + length(idObsParam))))
    namesParamEst::Array{String, 1} = string.(vcat(string.(idParamDyn), string.(idSdParam), string.(idObsParam)))

    paramIndicies = ParameterIndices(iDynPar, 
                                     iObsPar, 
                                     iSdPar, 
                                     string.(idParamDyn), 
                                     string.(idObsParam), 
                                     string.(idSdParam),
                                     namesParamEst)

    return paramIndicies
end