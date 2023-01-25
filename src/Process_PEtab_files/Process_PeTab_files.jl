"""
    checkForPeTabFile(fileSearchFor::String, dirModel::String)::String

    Helper function to check in dirModel if a file starting with fileSearchFor exists. 
    If true return file path.
"""
function checkForPeTabFile(fileSearchFor::String, dirModel::String)::String

    filesDirModel = readdir(dirModel)
    iUse = findall(x -> occursin(fileSearchFor, x), filesDirModel)
    if length(iUse) > 1 
        @printf("Error : More than 1 file starting with %s in %s\n", fileSearchFor, filesDirModel)
    end
    if length(iUse) == 0
        @printf("Error : No file starting with %s in %s\n", fileSearchFor, filesDirModel)
    end

    return dirModel * filesDirModel[iUse[1]]
end


function getPriorInfo(θ_indices::ParameterIndices, parameterDataFile::DataFrame)::PriorInfo

    if "objectivePriorType" ∉ names(parameterDataFile)
        return PriorInfo(Array{Function, 1}(undef, 0), Bool[], false)
    end

    namesParamEst = string.(θ_indices.θ_estNames)
    priorLogPdf = Array{Function, 1}(undef, length(namesParamEst))
    priorOnParamScale = Array{Bool, 1}(undef, length(namesParamEst))
    paramID = string.(parameterDataFile[!, "parameterId"])

    contPrior = 0.0
    for i in eachindex(namesParamEst)

        iUse = findfirst(x -> x == namesParamEst[i], paramID)

        priorF = parameterDataFile[iUse, "objectivePriorType"]
        if ismissing(priorF)
            priorLogPdf[i] = noPrior
            priorOnParamScale[i] = false
            continue
        end

        priorVal = parse.(Float64, split(parameterDataFile[iUse, "objectivePriorParameters"], ";"))

        if priorF == "parameterScaleNormal"
            contPrior += logpdf(Normal(priorVal[1], priorVal[2]), log10(parameterDataFile[iUse, "nominalValue"]))
            priorLogPdf[i] = (p) -> logpdf(Normal(priorVal[1], priorVal[2]), p) 
            priorOnParamScale[i] = true
        elseif priorF == "parameterScaleLaplace"
            priorLogPdf[i] = (p) -> logpdf(Laplace(priorVal[1], priorVal[2]), p) 
            priorOnParamScale[i] = true
        elseif priorF == "normal"
            priorLogPdf[i] = (p) -> logpdf(Normal(priorVal[1], priorVal[2]), p) 
            priorOnParamScale[i] = false
        elseif priorF == "laplace"
            priorLogPdf[i] = (p) -> logpdf(Laplace(priorVal[1], priorVal[2]), p) 
            priorOnParamScale[i] = false
        elseif priorF == "logNormal"
            priorLogPdf[i] = (p) -> logpdf(LogNormal(priorVal[1], priorVal[2]), p) 
            priorOnParamScale[i] = false
        elseif priorF == "logLaplace"
            println("Error : Julia does not yet have support for log-laplace")
        else
            println("Error : PeTab standard does not support a prior of type ", priorF)
        end

    end

    return PriorInfo(priorLogPdf, priorOnParamScale, true)
end
# Helper function in case there is not any parameter priors 
function noPrior(p::Real)::Real
    return 0.0
end


# Function generating callbacksets for time-depedent SBML piecewise expressions
function getCallbacksForTimeDepedentPiecewise(odeSys::ODESystem, modelDict::Dict, modelName::String, dirModel::String)

    # ParamEstIndices is needed to see if event-triggers contain parameters we want to estimate 
    experimentalConditionsFile, measurementDataFile, parameterDataFile, observablesDataFile = readPEtabFiles(dirModel, readObservables=true)
    parameterData = processParameters(parameterDataFile)
    measurementInfo = processMeasurements(measurementDataFile, observablesDataFile) 
    θ_indices = computeIndicesθ(parameterData, measurementInfo, odeSys, experimentalConditionsFile)

    parameterNames = parameters(odeSys)
    stateNames = string.(states(odeSys))
    parameterNames = string.(parameters(odeSys))
    
    stringWrite = "function getCallbacks_" * modelName * "()\n"
    stringWriteFunctions = ""
    stringWriteTstops = "function getTstops(u, p)\n"

    # In case we do not have any events 
    if isempty(modelDict["boolVariables"])
        callbackNames = ""
        checkActivatedNames = ""
        stringWriteTstops *= "\t return Float64[]\nend"
    else
        # In case we have events loop over each variable 
        for key in keys(modelDict["boolVariables"])
            stringFunctions, stringSet =  createCallbackCont(key, modelDict, parameterNames, stateNames, θ_indices) 
            stringWrite *= stringSet * "\n"
            stringWriteFunctions *= stringFunctions * "\n"
        end
        callbackNames = prod(["cb_" * key * ", " for key in keys(modelDict["boolVariables"])])[1:end-2]
        checkActivatedNames = prod(["activeAtTime0_" * key * "!, " for key in keys(modelDict["boolVariables"])])[1:end-2]

        stringWriteTstops *= "\t return " * createFuncForTstops(modelDict, stateNames, parameterNames) * "\nend"
    end


    # Put everything together
    stringWrite *= "\treturn CallbackSet(" * callbackNames * "), [" * checkActivatedNames * "]\nend"
    fileWrite = dirModel * "/" * modelName * "Callbacks_time_piecewise.jl"
    if isfile(fileWrite)
        rm(fileWrite)
    end
    io = open(fileWrite, "w")    

    write(io, stringWrite * "\n\n")
    write(io, stringWriteFunctions)
    write(io, stringWriteTstops)
    
    close(io)
end


function createCallbackCont(paramName::String, 
                            modelDict::Dict, 
                            parameterNames::Vector{String}, 
                            stateNames::Vector{String}, 
                            θ_indices::ParameterIndices)

    stringCallbackFunctions = ""
    stringCallbackSet = ""

    # Build the condition statement 
    stringCallbackFunctions *= "\nfunction condition_" * paramName * "(u, t, integrator)\n"
    activationFormula = modelDict["boolVariables"][paramName][1]
    sideActivated = modelDict["boolVariables"][paramName][2]

    # Check if the event trigger depend on parameters which are to be i) estimated, or ii) if it depend on models state. 
    # For i) it must be a cont. event in order for us to be able to compute the gradient. For ii) we cannot compute 
    # tstops (the event times) prior to starting to solve the ODE.
    hasStates = triggerHasStates(activationFormula, stateNames)
    hasParametersEst = triggerHasParametersToEst(activationFormula, parameterNames, θ_indices)
    discreteEvent = true
    if hasParametersEst == true || hasStates == true
        discreteEvent = false
    end
    
    # Replace any state or parameter with their corresponding index in the ODE system 
    for i in eachindex(stateNames)
        activationFormula = replaceWholeWord(activationFormula, stateNames[i], "u["*string(i)*"]")
    end
    for i in eachindex(parameterNames)
        activationFormula = replaceWholeWord(activationFormula, parameterNames[i], "integrator.p["*string(i)*"]")
    end

    # Replace inequality 
    if discreteEvent == false
        activationFormulaCond = replace(activationFormula, "<=" => "-")
        activationFormulaCond = replace(activationFormulaCond, ">=" => "-")
        activationFormulaCond = replace(activationFormulaCond, ">" => "-")
        activationFormulaCond = replace(activationFormulaCond, "<" => "-")
    else
        activationFormulaCond = replace(activationFormula, "<=" => "==")
        activationFormulaCond = replace(activationFormulaCond, ">=" => "==")
        activationFormulaCond = replace(activationFormulaCond, ">" => "==")
        activationFormulaCond = replace(activationFormulaCond, "<" => "==")
    end
    stringCallbackFunctions *= "\t" * activationFormulaCond * "\nend\n"

    # Build the affect function 
    whichParam = findfirst(x -> x == paramName, parameterNames)
    stringCallbackFunctions *= "function affect_" * paramName * "!(integrator)\n"
    stringCallbackFunctions *= "\tintegrator.p[" * string(whichParam) * "] = 1.0\nend\n"

    # Build the callback 
    if discreteEvent == false
        stringCallbackSet *= "\tcb_" * paramName * " = ContinuousCallback(" * "condition_" * paramName * ", " * "affect_" * paramName * "!, "
    else
        stringCallbackSet *= "\tcb_" * paramName * " = DiscreteCallback(" * "condition_" * paramName * ", " * "affect_" * paramName * "!, "
    end
    stringCallbackSet *= "save_positions=(false, false))\n"

    # Build a function to check if a condition should be true at time zero (as in Julia events activated at time zero 
    # are usually not triggered by default)
    sideInequality = sideActivated == "right" ? "!" : ""
    stringCallbackFunctions *= "function activeAtTime0_" * paramName * "!(u, p)\n"
    stringCallbackFunctions *= "\tt = 0.0 # Used to check conditions activated at t0=0\n"
    stringCallbackFunctions *= "\tp[" * string(whichParam) * "] = 0.0 # Default to being off\n"
    activationFormulaBool = replace(activationFormula, "integrator." => "")
    activationFormulaBool = replace(activationFormulaBool, "<=" => "≤")
    activationFormulaBool = replace(activationFormulaBool, ">=" => "≥")
    stringCallbackFunctions *= "\tif " * sideInequality *"(" * activationFormulaBool * ")\n"
    stringCallbackFunctions *= "\t\tp[" * string(whichParam) * "] = 1.0\n\tend\nend\n"

    return stringCallbackFunctions, stringCallbackSet
end


function triggerHasStates(activationFormula::AbstractString, stateNames::Vector{String})::Bool
    for i in eachindex(stateNames)
        activationFormulaNew = replaceWholeWord(activationFormula, stateNames[i], "u["*string(i)*"]")
        if activationFormulaNew != activationFormula
            return true
        end
    end
    return false
end


function triggerHasParametersToEst(activationFormula::AbstractString, parameterNames::Vector{String}, θ_indices::ParameterIndices)::Bool

    # Parameters which are present for each experimental condition, and condition specific parameters 
    indexParamEstODESys = θ_indices.mapODEProblem.iODEProblemθDynamic
    indexParamEstODESysCond = reduce(vcat, [θ_indices.mapsConiditionId[i].iODEProblemθDynamic for i in keys(θ_indices.mapsConiditionId)])

    for i in eachindex(parameterNames)
        activationFormulaNew = replaceWholeWord(activationFormula, parameterNames[i], "integrator.p["*string(i)*"]")
        if activationFormulaNew != activationFormula
            if i ∈ indexParamEstODESys || i ∈ indexParamEstODESysCond
                return true
            end
        end
    end

    return false
end


# Function computing t-stops (time for events) for piecewise expressions using the symbolics package 
function createFuncForTstops(modelDict::Dict, stateNames::Vector{String}, parameterNames::Vector{String})

    tStopExpressions = Array{String, 1}(undef, length(keys(modelDict["boolVariables"])))

    i = 1
    for key in keys(modelDict["boolVariables"])

        condFormula = modelDict["boolVariables"][key][1]
        # In case the activation formula contains a state we cannot precompute the t-stop time as it depends on 
        # the actual ODE solution.
        if triggerHasStates(condFormula, stateNames) == true
            i += 1
            continue
        end
        
        stringVariables = "@variables t, "
        stringVariables *= prod(string.(collect(keys(modelDict["parameters"]))) .* ", " )[1:end-2] * " "
        stringVariables *= prod(string.(collect(keys(modelDict["states"]))) .* ", " )[1:end-2] 
        symbolicVariables = eval(Meta.parse(stringVariables))

        # Note - below order counts (e.g having < first results in ~= incase what actually stands is <=)
        condFormula = replace(condFormula, "<=" => "~")
        condFormula = replace(condFormula, ">=" => "~")
        condFormula = replace(condFormula, "<" => "~")
        condFormula = replace(condFormula, ">" => "~")
        expSymbolic = eval(Meta.parse(condFormula))

        # Expression for the time at which the condition is triggered
        expForTime = string.(Symbolics.solve_for(expSymbolic, symbolicVariables[1], simplify=true))

        for i in eachindex(stateNames)
            expForTime = replaceWholeWord(expForTime, stateNames[i], "u["*string(i)*"]")
        end
        for i in eachindex(parameterNames)
            expForTime = replaceWholeWord(expForTime, parameterNames[i], "p["*string(i)*"]")
        end
        tStopExpressions[i] = "dualToFloat(" * expForTime * ")"
        i += 1
    end

    stringTstop = "[" * prod([str * ", " for str in tStopExpressions])[1:end-2] * "]"
    return stringTstop
end

