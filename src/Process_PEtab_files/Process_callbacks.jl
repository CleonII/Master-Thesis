# Function generating callbacksets for time-depedent SBML piecewise expressions, as callbacks are more efficient than 
# using ifelse (e.g better integration stabillity)
function createCallbacksForTimeDepedentPiecewise(odeSystem::ODESystem, 
                                                 SBMLDict::Dict, 
                                                 modelName::String, 
                                                 dirModel::String)

    pODEProblemNames = string.(parameters(odeSystem))
    modelStateNames = replace.(string.(states(odeSystem)), "(t)" => "")
    
    # Compute indices tracking parameters (needed as down the line we need to know if a parameter should be estimated 
    # or not, as if such a parameter triggers a callback we must let it be a continious callback)
    experimentalConditions, measurementsData, parametersData, observablesData = readPEtabFiles(dirModel, readObservables=true)
    parameterInfo = processParameters(parametersData) 
    measurementInfo = processMeasurements(measurementsData, observablesData) 
    θ_indices = computeIndicesθ(parameterInfo, measurementInfo, odeSystem, experimentalConditions)                                                  

    # In case of no-callbacks the function for getting callbacks will be empty, likewise for the function 
    # which compute tstops (callback-times)     
    stringWriteCallbacks = "function getCallbacks_" * modelName * "()\n"
    stringWriteFunctions = ""
    stringWriteTstops = "function computeTstops(u::AbstractVector, p::AbstractVector)\n"

    # In case we do not have any events 
    if isempty(SBMLDict["boolVariables"])
        callbackNames = ""
        checkIfActivatedT0Names = ""
        stringWriteTstops *= "\t return Float64[]\nend"
    else
        for key in keys(SBMLDict["boolVariables"])
            functionsStr, callbackStr =  createCallback(key, SBMLDict, pODEProblemNames, modelStateNames, θ_indices) 
            stringWriteCallbacks *= callbackStr * "\n"
            stringWriteFunctions *= functionsStr * "\n"
        end
        callbackNames = prod(["cb_" * key * ", " for key in keys(SBMLDict["boolVariables"])])[1:end-2]
        checkIfActivatedT0Names = prod(["isActiveAtTime0_" * key * "!, " for key in keys(SBMLDict["boolVariables"])])[1:end-2]

        stringWriteTstops *= "\t return " * createFuncionForTstops(SBMLDict, modelStateNames, pODEProblemNames, θ_indices) * "\nend"
    end

    stringWriteCallbacks *= "\treturn CallbackSet(" * callbackNames * "), [" * checkIfActivatedT0Names * "]\nend"
    fileWrite = dirModel * "/" * modelName * "_callbacks.jl"
    if isfile(fileWrite)
        rm(fileWrite)
    end
    io = open(fileWrite, "w")    

    write(io, stringWriteCallbacks * "\n\n")
    write(io, stringWriteFunctions)
    write(io, stringWriteTstops)
    
    close(io)
end


function createCallback(callbackName::String, 
                        SBMLDict::Dict, 
                        pODEProblemNames::Vector{String}, 
                        modelStateNames::Vector{String}, 
                        θ_indices::ParameterIndices)

    # Check if the event trigger depend on parameters which are to be i) estimated, or ii) if it depend on models state. 
    # For i) it must be a cont. event in order for us to be able to compute the gradient. For ii) we cannot compute 
    # tstops (the event times) prior to starting to solve the ODE so it most be cont. callback
    _conditionFormula = SBMLDict["boolVariables"][callbackName][1]
    hasModelStates = conditionHasStates(_conditionFormula, modelStateNames)
    hasParametersToEstimate = conditionHasParametersToEstimate(_conditionFormula, pODEProblemNames, θ_indices)
    discreteEvent = hasParametersToEstimate == true || hasModelStates == true ? false : true
        
    # Replace any state or parameter with their corresponding index in the ODE system to be comaptible with event 
    # syntax 
    for i in eachindex(modelStateNames)
        _conditionFormula = replaceWholeWord(_conditionFormula, modelStateNames[i], "u["*string(i)*"]")
    end
    for i in eachindex(pODEProblemNames)
        _conditionFormula = replaceWholeWord(_conditionFormula, pODEProblemNames[i], "integrator.p["*string(i)*"]")
    end

    # Replace inequality with - (root finding cont. event) or with == in case of 
    # discrete event 
    replaceWith = discreteEvent == true ? "==" : "-"
    conditionFormula = replace(_conditionFormula, "<=" => replaceWith)
    conditionFormula = replace(conditionFormula, ">=" => replaceWith)
    conditionFormula = replace(conditionFormula, ">" => replaceWith)
    conditionFormula = replace(conditionFormula, "<" => replaceWith)

    # Build the condition statement used in the jl function 
    conditionStr = "\nfunction condition_" * callbackName * "(u, t, integrator)\n"
    conditionStr *= "\t" * conditionFormula * "\nend\n"

    # Build the affect function 
    whichParameter = findfirst(x -> x == callbackName, pODEProblemNames)
    affectStr = "function affect_" * callbackName * "!(integrator)\n"
    affectStr *= "\tintegrator.p[" * string(whichParameter) * "] = 1.0\nend\n"

    # Build the callback 
    if discreteEvent == false
        callbackStr = "\tcb_" * callbackName * " = ContinuousCallback(" * "condition_" * callbackName * ", " * "affect_" * callbackName * "!, "
    else
        callbackStr = "\tcb_" * callbackName * " = DiscreteCallback(" * "condition_" * callbackName * ", " * "affect_" * callbackName * "!, "
    end
    callbackStr *= "save_positions=(false, false))\n" # So we do not get problems with saveat in the ODE solver 

    # Building a function which check if a callback is activated at time zero (as this is not something Julia will 
    # check for us)
    sideInequality = SBMLDict["boolVariables"][callbackName][2] == "right" ? "!" : "" # Check if true or false evaluates expression to true
    activeAtT0Str = "function isActiveAtTime0_" * callbackName * "!(u, p)\n"
    activeAtT0Str *= "\tt = 0.0 # Used to check conditions activated at t0=0\n" * "\tp[" * string(whichParameter) * "] = 0.0 # Default to being off\n"
    conditionFormula = replace(_conditionFormula, "integrator." => "")
    conditionFormula = replace(conditionFormula, "<=" => "≤")
    conditionFormula = replace(conditionFormula, ">=" => "≥")
    activeAtT0Str *= "\tif " * sideInequality *"(" * conditionFormula * ")\n" * "\t\tp[" * string(whichParameter) * "] = 1.0\n\tend\nend\n"

    functionsStr = conditionStr * '\n' * affectStr * '\n' * activeAtT0Str * '\n'

    return functionsStr, callbackStr
end


# Function computing t-stops (time for events) for piecewise expressions using the symbolics package 
# to symboically solve for where the condition is zero
function createFuncionForTstops(SBMLDict::Dict, modelStateNames::Vector{String}, pODEProblemNames::Vector{String}, θ_indices::ParameterIndices)

    tstopsStr = Vector{String}(undef, length(keys(SBMLDict["boolVariables"])))
    i = 1
    for key in keys(SBMLDict["boolVariables"])

        conditionFormula = SBMLDict["boolVariables"][key][1]
        # In case the activation formula contains a state we cannot precompute the t-stop time as it depends on 
        # the actual ODE solution.
        if conditionHasStates(conditionFormula, modelStateNames) || conditionHasParametersToEstimate(conditionFormula, pODEProblemNames, θ_indices)
            tstopsStr[i] = ""
            i += 1
            continue
        end
        
        # We need to make the parameters and states symbolic in order to solve the condition expression 
        # using the Symbolics package.
        variablesStr = "@variables t, "
        variablesStr *= prod(string.(collect(keys(SBMLDict["parameters"]))) .* ", " )[1:end-2] * " "
        variablesStr *= prod(string.(collect(keys(SBMLDict["states"]))) .* ", " )[1:end-2] 
        variablesSymbolic = eval(Meta.parse(variablesStr))

        # Note - below order counts (e.g having < first results in ~= incase what actually stands is <=)
        conditionFormula = replace(conditionFormula, "<=" => "~")
        conditionFormula = replace(conditionFormula, ">=" => "~")
        conditionFormula = replace(conditionFormula, "<" => "~")
        conditionFormula = replace(conditionFormula, ">" => "~")
        conditionSymbolic = eval(Meta.parse(conditionFormula))

        # Expression for the time at which the condition is triggered
        expressionForTime = string.(Symbolics.solve_for(conditionSymbolic, variablesSymbolic[1], simplify=true))

        # Make compatible with the PEtab importer syntax 
        for i in eachindex(modelStateNames)
            expressionForTime = replaceWholeWord(expressionForTime, modelStateNames[i], "u["*string(i)*"]")
        end
        for i in eachindex(pODEProblemNames)
            expressionForTime = replaceWholeWord(expressionForTime, pODEProblemNames[i], "p["*string(i)*"]")
        end
        # dualToFloat is needed as tstops for the integrator cannot be of type Dual 
        tstopsStr[i] = "dualToFloat(" * expressionForTime * ")"
        i += 1
    end

    return "Float64[" * prod([isempty(tstopsStr[i]) ? "" : tstopsStr[i] * ", " for i in eachindex(tstopsStr)])[1:end-2] * "]"
end


function conditionHasStates(conditionFormula::AbstractString, modelStateNames::Vector{String})::Bool
    for i in eachindex(modelStateNames)
        _conditionFormula = replaceWholeWord(conditionFormula, modelStateNames[i], "")
        if _conditionFormula != conditionFormula
            return true
        end
    end
    return false
end


function conditionHasParametersToEstimate(conditionFormula::AbstractString, 
                                          pODEProblemNames::Vector{String}, 
                                          θ_indices::ParameterIndices)::Bool

    # Parameters which are present for each experimental condition, and condition specific parameters 
    iODEProblemθConstantDynamic = θ_indices.mapODEProblem.iODEProblemθDynamic
    iODEProblemθDynamicCondition = reduce(vcat, [θ_indices.mapsConiditionId[i].iODEProblemθDynamic for i in keys(θ_indices.mapsConiditionId)])

    for i in eachindex(pODEProblemNames)
        _conditionFormula = replaceWholeWord(conditionFormula, pODEProblemNames[i], "integrator.p["*string(i)*"]")
        if _conditionFormula != conditionFormula
            if i ∈ iODEProblemθConstantDynamic || i ∈ iODEProblemθDynamicCondition
                return true
            end
        end
    end
    return false
end