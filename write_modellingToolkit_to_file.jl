using PyCall

#ENV["PYTHON"] = "H:\\Program\\Python37\\python.exe"
#Pkg.build("PyCall")

libsbml = pyimport("libsbml")
reader = libsbml.SBMLReader()
modelName = "model_Crauste_CellSystems2017"
document = reader[:readSBML]("./SBML/Benchmark/" * modelName * ".xml")
model = document[:getModel]() # Get the model

# Splits strings by comma (,), but only if the comma is not inside a function / parenthesis.
function splitByComma(stringToSplit)::Vector{SubString{String}}
    parts = Vector{SubString{String}}(undef, trunc(Int, length(stringToSplit)))
    numParts = 0
    inParenthesis = 0
    startPart = 1
    endPart = 1
    for i = 1:length(stringToSplit)
        if stringToSplit[i] == '('
            inParenthesis += 1
        elseif stringToSplit[i] == ')'
            inParenthesis -= 1
        end
        if stringToSplit[i] == ',' && inParenthesis == 0
            endPart = i-1
            numParts += 1
            parts[numParts] = stringToSplit[startPart:endPart]
            startPart = i+2
        end
    end
    numParts += 1
    parts[numParts] = stringToSplit[startPart:end]
    parts = parts[1:numParts]
end

# replaces the word "time" in functions with "t". 
# Makes sure not to change for example "time1" or "shift_time"
function replaceTime(timeString)
    if timeString == "time"
        return "t"
    end
    newTimeString = timeString
    i = 1
    while i < length(newTimeString)
        indices = findnext("time", newTimeString, i)
        if indices === nothing
            break
        end
        if indices[1] == 1
            if newTimeString[indices[end]+1] in [' ', ')']
                newTimeString = newTimeString[1:indices[1]-1] * "t" * newTimeString[indices[end]+1:end]
            end
        elseif indices[end] == length(newTimeString)
            if newTimeString[indices[1]-1] in [' ', '(']
                newTimeString = newTimeString[1:indices[1]-1] * "t" * newTimeString[indices[end]+1:end]
            end
        else
            if newTimeString[indices[1]-1] in [' ', '('] && newTimeString[indices[end]+1] in [' ', ')']
                newTimeString = newTimeString[1:indices[1]-1] * "t" * newTimeString[indices[end]+1:end]
            end
        end
        i = indices[1] + 1
    end
    return newTimeString
end

# Goes through a condition for an event written for SBML and creates the appropriate events in ModelingToolkit.
# May call itself to handle more complicated events.
# Makes sure there are no doublett events. 
# Also adds constants and parameters used to a dummy variable so that ModelingToolkit does not think they are unused.
function goToBottomEvent(condition, variable, valActive, valInactive, eventDict, dicts)   
    if "and" == condition[1:3]
        strippedCondition = condition[5:end-1]
        parts = splitByComma(strippedCondition)
        for part in parts
            goToBottomEvent(part, variable, valActive, valInactive, eventDict, dicts)
        end
    elseif "or" == condition[1:2]
        strippedCondition = condition[4:end-1]
        parts = splitByComma(strippedCondition)
        for part in parts
            goToBottomEvent(part, variable, valActive, valInactive, eventDict, dicts)
        end
    elseif "geq" == condition[1:3] || "gt" == condition[1:2]
        if "geq" == condition[1:3]
            strippedCondition = condition[5:end-1]
        else
            strippedCondition = condition[4:end-1]
        end
        # if the trigger or event contains parameters or constant, add them to the dummyVariable to keep defined
        args = split(getArguments(strippedCondition), ", ")
        for arg in args
            if arg in keys(dicts["parameters"])
                dicts["dummyVariable"][arg] = dicts["parameters"][arg]
            elseif arg in keys(dicts["constants"])
                dicts["dummyVariable"][arg] = dicts["constants"][arg]
            end
        end

        parts = splitByComma(strippedCondition)
        if occursin("time", parts[1])
            parts[1] = replaceTime(parts[1])
        end
        trigger = "[" * parts[1] * " ~ " * parts[2] * "]"
        event = "[" * variable * " ~ " * valActive * "]"

        eventDict[trigger] = event
    elseif "leq" == condition[1:3] || "lt" == condition[1:2]
        if "leq" == condition[1:3]
            strippedCondition = condition[5:end-1]
        else
            strippedCondition = condition[4:end-1]
        end
        # if the trigger or event contains parameters or constant, add them to the dummyVariable to keep defined
        args = split(getArguments(strippedCondition), ", ")
        for arg in args
            if arg in keys(dicts["parameters"])
                dicts["dummyVariable"][arg] = dicts["parameters"][arg]
            elseif arg in keys(dicts["constants"])
                dicts["dummyVariable"][arg] = dicts["constants"][arg]
            end
        end

        parts = splitByComma(strippedCondition)
        if occursin("time", parts[1])
            parts[1] = replaceTime(parts[1])
        end
        trigger = "[" * parts[1] * " ~ " * parts[2] * "]"
        event = "[" * variable * " ~ " * valInactive * "]"

        if ~(trigger in keys(eventDict))
            eventDict[trigger] = event
        end
    end
    nothing
end

# Handles piecewise functions that are to be redefined as events. 
# will add the variable affected by the event as a variable in ModelingToolkit if it is not already.
# Calls goToBottomEvent to create the actual events and then combines them.
function piecewiseToEvent(piecewiseString, variable, dicts)
    piecewiseString = piecewiseString[11:end-1]

    # if the variable in the event is set as a parameter or constant, set it as a variable and delete the previous form
    if variable in keys(dicts["parameters"])
        dicts["variableParameters"][variable] = dicts["parameters"][variable]
        delete!(dicts["parameters"], variable)
    elseif variable in keys(dicts["constants"])
        dicts["variableParameters"][variable] = dicts["constants"][variable]
        delete!(dicts["constants"], variable)
    end

    eventDict = Dict()
    args = splitByComma(piecewiseString)
    vals = args[1:2:end]
    conds = args[2:2:end]

    # if any values are a/a set of parameter(s) or constant(s), add it/them to the dummyVariable to keep defined, also insert function definitions as needed
    for (vIndex, val) in enumerate(vals)
        valArgs = getArguments(val, dicts["functions"])
        if valArgs[2]
            vals[vIndex] = rewriteFunctionOfFunction(val, dicts["functions"])
        end
        for valArg in split(valArgs[1], ", ")
            if valArg in keys(dicts["parameters"])
                dicts["dummyVariable"][valArg] = dicts["parameters"][valArg]
            elseif valArg in keys(dicts["constants"])
                dicts["dummyVariable"][valArg] = dicts["constants"][valArg]
            end
        end
    end    

    for (cIndex, condition) in enumerate(conds)
        valActive = vals[cIndex]
        valInactive = vals[end]
        goToBottomEvent(condition, variable, valActive, valInactive, eventDict, dicts)
    end

    eventString = ""
    for (trigger, event) in eventDict
        if length(eventString) == 0
            eventString = trigger * " => " * event
        else
            eventString = eventString * ", " * trigger * " => " * event
        end
    end
    return eventString
end

# Rewrites triggers in events to propper form for ModelingToolkit
function asTrigger(triggerFormula)
    if "geq" == triggerFormula[1:3]
        strippedFormula = triggerFormula[5:end-1]
    elseif "gt" == triggerFormula[1:2]
        strippedFormula = triggerFormula[4:end-1]
    elseif "leq" == triggerFormula[1:3]
        strippedFormula = triggerFormula[5:end-1]
    elseif "lt" == triggerFormula[1:2]
        strippedFormula = triggerFormula[4:end-1]
    end
    parts = splitByComma(strippedFormula)
    if occursin("time", parts[1])
        parts[1] = replaceTime(parts[1])
    end
    expression = "[" * parts[1] * " ~ " * parts[2] * "]"
    return expression
end

# Extracts the argument from a function.
# If a dictionary (with functions) is supplied, will also check if there are nested functions and will 
# include the arguments of these nested functions as arguments of the first function.
# The returned string will only contain unique arguments.
function getArguments(functionAsString)
    parts = split(functionAsString, ['(', ')', '/', '+', '-', '*', ' ', '~', '>', '<', '=', ','], keepempty = false)
    arguments = Dict()
    for part in parts
        if isdigit(part[1])
            nothing
        else
            if (part in values(arguments)) == false
                arguments[length(arguments)+1] = part
            end
        end
    end
    if length(arguments) > 0
        argumentString = arguments[1]
        for i = 2:length(arguments)
            argumentString = argumentString * ", " * arguments[i]
        end
    else
        argumentString = ""
    end
    return argumentString
end
function getArguments(functionAsString, dictionary::Dict)
    parts = split(functionAsString, ['(', ')', '/', '+', '-', '*', ' ', '~', '>', '<', '=', ','], keepempty = false)
    existingFunctions = keys(dictionary)
    includesFunction = false
    arguments = Dict()
    for part in parts
        if isdigit(part[1])
            nothing
        else
            if part in existingFunctions
                includesFunction = true
                funcArgs = dictionary[part][1]
                funcArgs = split(funcArgs, [',', ' '], keepempty = false)
                for arg in funcArgs
                    if (arg in values(arguments)) == false
                        arguments[length(arguments)+1] = arg
                    end
                end
            else
                if (part in values(arguments)) == false
                    arguments[length(arguments)+1] = part
                end
            end
        end
    end
    if length(arguments) > 0
        argumentString = arguments[1]
        for i = 2:length(arguments)
            argumentString = argumentString * ", " * arguments[i]
        end
    else
        argumentString = ""
    end
    return [argumentString, includesFunction]
end

# Rewrites functions including variables that have been defined as functions in funcNameArgFormula.
# Replaces the names of the variables with the new function definitions.
# Makes sure not to replace variables with simmilar names.
function rewriteFunctionOfFunction(functionAsString, funcNameArgFormula)
    parts = split(functionAsString, ['(', ')', '/', '+', '-', '*', ' ', ','], keepempty = false)
    existingFunctions = keys(funcNameArgFormula)
    newFunctionString = functionAsString
    for part in parts
        if part in existingFunctions
            funcArgs = "(" * funcNameArgFormula[part][1] * ")"
            i = 1
            while i < length(newFunctionString)
                indices = findnext(part, newFunctionString, i)
                if indices[1] == 1 && indices[end] == length(newFunctionString)
                    break
                end
                if indices[1] == 1
                    if newFunctionString[indices[end]+1] in [' ', ')']
                        newFunctionString = newFunctionString[1:indices[end]] * funcArgs * newFunctionString[indices[end]+1:end]
                        break
                    end
                elseif indices[end] == length(newFunctionString)
                    if newFunctionString[indices[1]-1] in [' ', '(']
                        newFunctionString = newFunctionString[1:indices[end]] * funcArgs * newFunctionString[indices[end]+1:end]
                        break
                    end
                else
                    if newFunctionString[indices[1]-1] in [' ', '('] && newFunctionString[indices[end]+1] in [' ', ')']
                        newFunctionString = newFunctionString[1:indices[end]] * funcArgs * newFunctionString[indices[end]+1:end]
                        break
                    end
                end
                i = indices[end] + 1
            end
        end
    end
    return newFunctionString
end

# Differs from rewriteFunctionOfFunction in that it checks if the function is defined by the model.
# Will substitute the function definition for the formula given by the model with the arguments 
# given to the function used instead of the general arguments
function insertModelDefineFunctions(functionsAsString, modelFuncNameArgFormula)
    newFunctionsAsString = functionsAsString
    possibleModelFunctions = split(getArguments(functionsAsString), ", ")
    for possibleModelFunction in possibleModelFunctions 
        if possibleModelFunction in keys(modelFuncNameArgFormula)
            modelFuncArguments, modelFuncFormula = modelFuncNameArgFormula[possibleModelFunction]
            modelFuncArguments = split(modelFuncArguments, ", ")
            inParenthesis = 0
            started = 0
            startIndex = 1
            endIndex = 1
            i = 1
            while i <= length(functionsAsString)
                if started == 0
                    next = findnext(possibleModelFunction*"(", functionsAsString, i)
                    if next === nothing
                        break
                    end
                    startIndex = next[1]
                    i = startIndex + length(possibleModelFunction) + 1 
                    started = 1
                    inParenthesis = 1
                else 
                    if functionsAsString[i] == '('
                        inParenthesis += 1
                    elseif functionsAsString[i] == ')'
                        inParenthesis -= 1
                    end
                    if inParenthesis == 0
                        endIndex = i
                        modelFunction = functionsAsString[startIndex:endIndex]
                        argumentString = modelFunction[length(possibleModelFunction)+2:end-1]
                        arguments = split(argumentString, ", ", keepempty = false)
                        newFunction = modelFuncFormula
                        for (aIndex, arg) in enumerate(arguments)
                            j = 1
                            while j < length(newFunction)
                                indices = findnext(modelFuncArguments[aIndex], newFunction, j)
                                if indices === nothing
                                    break
                                end
                                if indices[1] == 1
                                    if newFunction[indices[end]+1] in [' ', '(']
                                        if occursin(" ", arg)
                                            newFunction = newFunction[1:indices[1]-1] * "(" * arg * ")" * newFunction[indices[end]+1:end]
                                        else
                                            newFunction = newFunction[1:indices[1]-1] * arg * newFunction[indices[end]+1:end]
                                        end
                                    end
                                elseif indices[end] == length(newFunction)
                                    if newFunction[indices[1]-1] in [' ', ')']
                                        if occursin(" ", arg)
                                            newFunction = newFunction[1:indices[1]-1] * "(" * arg * ")" * newFunction[indices[end]+1:end]
                                        else
                                            newFunction = newFunction[1:indices[1]-1] * arg * newFunction[indices[end]+1:end]
                                        end
                                    end
                                else
                                    if newFunction[indices[1]-1] in [' ', ')'] && newFunction[indices[end]+1] in [' ', '(']
                                        if occursin(" ", arg)
                                            newFunction = newFunction[1:indices[1]-1] * "(" * arg * ")" * newFunction[indices[end]+1:end]
                                        else
                                            newFunction = newFunction[1:indices[1]-1] * arg * newFunction[indices[end]+1:end]
                                        end
                                    end
                                end
                                j = indices[end] + 1
                            end
                        end
                        newFunctionsAsString = newFunctionsAsString[1:startIndex-1] * "(" * newFunction * ")" * newFunctionsAsString[startIndex + length(modelFunction):end]
                        startIndex = i
                        started = 0
                    end
                    i += 1
                end
            end
        end
    end
    return newFunctionsAsString
end

# Rewrites pow(a,b) into a^b, which Julia can handle
function removePowFunctions(functionAsString)
    newFunctionAsString = functionAsString
    inParenthesis = 0
    started = 0
    startIndex = 1
    endIndex = 1
    i = 1
    while i <= length(functionAsString)
        if started == 0
            next = findnext("pow(", functionAsString, i)
            if next === nothing
                break
            end
            startIndex = next[1]
            i = startIndex + 4
            started = 1
            inParenthesis = 1
        else 
            if functionAsString[i] == '('
                inParenthesis += 1
            elseif functionAsString[i] == ')'
                inParenthesis -= 1
            end
            if inParenthesis == 0
                endIndex = i
                powFunction = functionAsString[startIndex:endIndex]
                powArguments = powFunction[5:end-1]
                parts = split(powArguments, ", ", keepempty = false)
                newPowFunction = "(" * parts[1] * ")^(" * parts[2] * ")"
                newFunctionAsString = replace(newFunctionAsString, powFunction => newPowFunction, count = 1)
                startIndex = i
                started = 0
            end
            i += 1
        end
    end
    return newFunctionAsString
end

# Handles conditions of piecewise functions in derivatives and rewrites them into a propper format.
# may call itself to handle more complicated functions.
# returns both a trigger and an event to be combined in rewritePiecewiseInDerivative
function goToBottomDerivative(condition, valActive, valInactive)
    if valInactive == ""
        event = valActive
    else
        event = "(" * valActive * " - " * valInactive * ")"
    end
    if "and" == condition[1:3]
        strippedCondition = condition[5:end-1]
        parts = splitByComma(strippedCondition)
        combinedTrigger = ""
        for part in parts
            trigger,  = goToBottomDerivative(part, valActive, valInactive)
            combinedTrigger = combinedTrigger * " * " * trigger
        end
    elseif "or" == condition[1:2]
        strippedCondition = condition[4:end-1]
        parts = splitByComma(strippedCondition)
        combinedTrigger = ""
        for (pIndex, part) in enumerate(parts)
            trigger,  = goToBottomDerivative(part, valActive, valInactive)
            if pIndex == 1
                combinedTrigger = "(" * trigger
            else
                combinedTrigger = combinedTrigger * " + " * trigger
            end
            combinedTrigger = combinedTrigger * ")"
        end
    elseif "geq" == condition[1:3] || "gt" == condition[1:2]
        if "geq" == condition[1:3]
            strippedCondition = condition[5:end-1]
            operator = " >= "
        else
            strippedCondition = condition[4:end-1]
            operator = " > "
        end
        parts = splitByComma(strippedCondition)
        if occursin("time", parts[1])
            parts[1] = replaceTime(parts[1])
        end

        trigger = "(" * parts[1] * operator * parts[2] * ")"

    elseif "leq" == condition[1:3] || "lt" == condition[1:2]
        if "leq" == condition[1:3]
            strippedCondition = condition[5:end-1]
            operator = " <= "
        else
            strippedCondition = condition[4:end-1]
            operator = " < " 
        end
        parts = splitByComma(strippedCondition)
        if occursin("time", parts[1])
            parts[1] = replaceTime(parts[1])
        end

        trigger = "(" * parts[1] * operator * parts[2] * ")"

    end
    return trigger, event
end

# Handel piecewise functions in derivatives and rewrites them into a propper format.
# Calls goToBottomDerivative to handle each condition. 
function rewritePiecewiseInDerivative(functionAsString)
    newFunctionAsString = functionAsString
    inParenthesis = 0
    started = 0
    startIndex = 1
    endIndex = 1
    i = 1
    while i <= length(functionAsString)
        if started == 0
            next = findnext("piecewise(", functionAsString, i)
            if next === nothing
                break
            end
            startIndex = next[1]
            i = startIndex + 10
            started = 1
            inParenthesis = 1
        else 
            if functionAsString[i] == '('
                inParenthesis += 1
            elseif functionAsString[i] == ')'
                inParenthesis -= 1
            end
            if inParenthesis == 0
                endIndex = i
                piecewiseFunction = functionAsString[startIndex:endIndex]
                piecewiseString = piecewiseFunction[11:end-1]
                args = splitByComma(piecewiseString)
                vals = args[1:2:end]
                conds = args[2:2:end]

                if length(vals) > length(conds)
                    expression = vals[end]
                else
                    expression = ""
                end
                for (cIndex, condition) in enumerate(conds)
                    valActive = vals[cIndex]
                    if length(vals) > length(conds)
                        valInactive = vals[end]
                    else
                        valInactive = ""
                    end
                    trigger, event = goToBottomDerivative(condition, valActive, valInactive)
                    if expression == ""
                        expression = trigger * " * " * event
                    else
                        expression = "(" * expression * " + " * trigger * " * " * event * ")"
                    end
                end

                newFunctionAsString = replace(newFunctionAsString, piecewiseFunction => expression, count = 1)
                startIndex = i
                started = 0
            end
            i += 1
        end
    end
    return newFunctionAsString
end

# Insert the function definitions for newly defined functions. 
function insertFunctionDefinitions(functionAsString, funcNameArgFormula)
    parts = split(functionAsString, ['(', ')', '/', '+', '-', '*', ' ', ','], keepempty = false)
    existingFunctions = keys(funcNameArgFormula)
    newFunctionString = functionAsString
    for part in parts
        if part in existingFunctions
            funcArgs = "(" * funcNameArgFormula[part][1] * ")"
            i = 1
            while i < length(newFunctionString)
                indices = findnext(part, newFunctionString, i)
                if indices === nothing
                    break
                end
                if indices[1] == 1
                    if newFunctionString[indices[end]+1] in [' ', ')']
                        newFunctionString = newFunctionString[1:indices[end]] * funcArgs * newFunctionString[indices[end]+1:end]
                    end
                elseif indices[end] == length(newFunctionString)
                    if newFunctionString[indices[1]-1] in [' ', '(']
                        newFunctionString = newFunctionString[1:indices[end]] * funcArgs * newFunctionString[indices[end]+1:end]
                    end
                else
                    if newFunctionString[indices[1]-1] in [' ', '('] && newFunctionString[indices[end]+1] in [' ', ')']
                        newFunctionString = newFunctionString[1:indices[end]] * funcArgs * newFunctionString[indices[end]+1:end]
                    end
                end
                i = indices[1] + 1
            end
        end
    end
    return newFunctionString
end

# Rewrites derivatives using insertModelDefineFunctions, removePowFunctions, rewritePiecewiseInDerivative and insertFunctionDefinitions
function rewriteDerivatives(derivativeAsString, dicts)
    newDerivativeAsString = derivativeAsString
    newDerivativeAsString = insertModelDefineFunctions(newDerivativeAsString, dicts["modelFunctions"])
    if occursin("pow(", newDerivativeAsString)
        newDerivativeAsString = removePowFunctions(newDerivativeAsString)
    end
    if occursin("piecewise(", newDerivativeAsString)
        newDerivativeAsString = rewritePiecewiseInDerivative(newDerivativeAsString)
    end
    newDerivativeAsString = insertFunctionDefinitions(newDerivativeAsString, dicts["functions"])
    return newDerivativeAsString
end

function writeODEModelToFile()
    dicts = Dict()
    variableDict = Dict()
    dicts["variables"] = variableDict
    parameterDict = Dict()
    dicts["parameters"] = parameterDict
    constantsDict = Dict()
    dicts["constants"] = constantsDict
    variableParameterDict = Dict()
    dicts["variableParameters"] = variableParameterDict
    dummyVariableDict = Dict()
    dicts["dummyVariable"] = dummyVariableDict

    # for functions defined by the model, these functions will be rewritten
    modelFuncNameArgFormula = Dict()
    dicts["modelFunctions"] = modelFuncNameArgFormula
    
    # for functions rewritten from rules etc.
    funcNameArgFormula = Dict()
    dicts["functions"] = funcNameArgFormula
    dus = Dict()
    dicts["derivatives"] = dus
    

    ### Define extra functions

    ### Define variables and read initial value
    for spec in model[:getListOfSpecies]()
        variableDict[spec[:getId]()] = spec[:getInitialAmount]() == 0 ? string(spec[:getInitialConcentration]()) : string(spec[:getInitialAmount]())
    end
    
    ### Define parameters and read true parameter values
    for par in model[:getListOfParameters]()
        parameterDict[par[:getId]()] = string(par[:getValue]())
    end

    ### Define constants and read constants values
    for comp in model[:getListOfCompartments]()
        constantsDict[comp[:getId]()] = string(comp[:getSize]())
    end

    ### Define functions
    for (fIndex, functionDefinition) in enumerate(model[:getListOfFunctionDefinitions]())
        math = functionDefinition[:getMath]()
        functionName = functionDefinition[:getId]()
        args = "("
        if math[:getNumChildren]() > 1
            args = args * math[:getLeftChild]()[:getName]()
            for n in range(1, math[:getNumChildren]()-1, step = 1)
                arg = math[:getChild](n)[:getName]()
                if arg !== nothing
                    args = args * ", " * arg
                end
            end
            args = args * ")"
        end
        mathAsString = libsbml[:formulaToString](math)
        expressionStart = findfirst('(', mathAsString)+1
        expressionStop = findlast(')', mathAsString)-1
        StrippedMathAsString = mathAsString[expressionStart:expressionStop]
        functionFormula = lstrip(split(StrippedMathAsString, ',')[end])
        functionFormula = removePowFunctions(functionFormula)
        modelFuncNameArgFormula[functionName] = [args[2:end-1], functionFormula]
    end

    ### Define events
    stringOfEvents = ""
    for (eIndex, event) in enumerate(model[:getListOfEvents]())
        eventName = event[:getName]()
        trigger = event[:getTrigger]()
        triggerMath = trigger[:getMath]()
        triggerFormula = asTrigger(libsbml[:formulaToString](triggerMath))
        # if the trigger contains parameters or constants, add them to the dummy variables
        triggerArgs = split(getArguments(triggerFormula[2:end-1]), ", ")
        for triggerArg in triggerArgs
            if triggerArg in keys(dicts["parameters"])
                dicts["dummyVariable"][triggerArg] = dicts["parameters"][triggerArg]
            elseif arg in keys(dicts["constants"])
                dicts["dummyVariable"][triggerArg] = dicts["constants"][triggerArg]
            end
        end
        eventAsString = ""
        for (eaIndex, eventAssignment) in enumerate(event[:getListOfEventAssignments]())
            variableName = eventAssignment[:getVariable]()
            eventMath = eventAssignment[:getMath]()
            eventMathAsString = libsbml[:formulaToString](eventMath)
            if eaIndex == 1
                eventAsString = "[" * variableName * " ~ " * eventMathAsString
            else
                eventAsString = eventAsString * ", " * variableName * " ~ " * eventMathAsString
            end
            # if the variable in the event is not set as a variable, make it so and remove it as a parameter or constant
            if variableName in keys(parameterDict)
                variableParameterDict[lhsArg] = parameterDict[lhsArg]
                delete!(parameterDict, lhsArg)
            elseif variableName in keys(constantsDict)
                variableParameterDict[lhsArg] = constantsDict[lhsArg]
                delete!(constantsDict, lhsArg)
            end
            # if the event math contains parameters or constants, add them to the dummy variable so that they stay defined 
            eventMathArgs = split(getArguments(eventMathAsString, funcNameArgFormula)[1], ", ")
            for eventMathArg in eventMathArgs
                if eventMathArg in keys(parameterDict)
                    dummyVariableDict[eventMathArg] = parameterDict[eventMathArg]
                elseif eventMathArg in keys(constantsDict)
                    dummyVariableDict[eventMathArg] = constantsDict[eventMathArg]
                end
            end
        end
        eventAsString = eventAsString * "]"
        fullEvent = triggerFormula * " => " * eventAsString
        if eIndex == 1
            stringOfEvents = fullEvent
        else
            stringOfEvents = stringOfEvents * ", " * fullEvent
        end
    end

    ### Define rules
    for rule in model[:getListOfRules]()
        ruleType = rule[:getElementName]() # type
        if ruleType == "assignmentRule"
            ruleVariable = rule[:getVariable]() # variable
            ruleFormula = rule[:getFormula]() # need to add piecewise function (and others)
            if occursin("piecewise", ruleFormula)
                if ruleFormula[1:9] == "piecewise" && findnext("piecewise", ruleFormula, 2) === nothing
                    # adds rule as event
                    localEventString = piecewiseToEvent(ruleFormula, ruleVariable, dicts)
                    if length(stringOfEvents) == 0
                        stringOfEvents = localEventString
                    else
                        stringOfEvents = stringOfEvents * ", " * localEventString
                    end
                else
                    # adds rule as function
                    ruleFormula = rewritePiecewiseInDerivative(ruleFormula)
                    ruleFormula = replaceTime(ruleFormula)
                    arguments, includesFunction = getArguments(ruleFormula, funcNameArgFormula)
                    ruleFormula = removePowFunctions(ruleFormula)
                    if includesFunction
                        ruleFormula = rewriteFunctionOfFunction(ruleFormula, funcNameArgFormula)
                    end
                    funcNameArgFormula[ruleVariable] = [arguments, ruleFormula]
                    if ruleVariable in keys(variableDict)
                        variableDict = delete!(variableDict, ruleVariable)
                    end
                    if ruleVariable in keys(parameterDict)
                        parameterDict = delete!(parameterDict, ruleVariable)
                    end
                end
                
            else
                # adds rule as function
                arguments, includesFunction = getArguments(ruleFormula, funcNameArgFormula)
                if arguments == ""
                    parameterDict[ruleVariable] = ruleFormula
                else
                    ruleFormula = removePowFunctions(ruleFormula)
                    if includesFunction
                        ruleFormula = rewriteFunctionOfFunction(ruleFormula, funcNameArgFormula)
                    end
                    funcNameArgFormula[ruleVariable] = [arguments, ruleFormula]
                    if ruleVariable in keys(variableDict)
                        variableDict = delete!(variableDict, ruleVariable)
                    end
                    if ruleVariable in keys(parameterDict)
                        parameterDict = delete!(parameterDict, ruleVariable)
                    end
                end
            end
        elseif ruleType == "algebraicRule"
            # TODO
        elseif ruleType == "rateRule"
            # TODO
        end
    end

    ### Implement Initial assignments 
    # Positioned after rules since some assignments may include functions
    for initAssign in model[:getListOfInitialAssignments]()
        assignName = initAssign[:getId]()
        assignMath = initAssign[:getMath]()
        assignFormula = libsbml[:formulaToString](assignMath)
        _, inFunction = getArguments(assignFormula, funcNameArgFormula)
        if inFunction
            assignFormula = rewriteFunctionOfFunction(assignFormula, funcNameArgFormula)
        end
        if assignName in keys(variableDict)
            variableDict[assignName] = assignFormula
        elseif assignName in keys(parameterDict)
            parameterDict[assignName] = assignFormula
        elseif assignName in keys(constantsDict)
            constantsDict[assignName] = assignFormula
        else
            println("Error: could not find assigned variable/parameter")
        end
    end

    ### Defining derivatives
    for spec in model[:getListOfSpecies]()
        dus[spec[:getId]()] = "D(" * spec[:getId]() * ") ~ "
    end

    reactions = [(r, r[:getKineticLaw]()[:getFormula]()) for r in model[:getListOfReactions]()]
    for (reac, formula) in reactions
        products = [(p[:species], p[:getStoichiometry]()) for p in reac[:getListOfProducts]()]
        reactants = [(r[:species], r[:getStoichiometry]()) for r in reac[:getListOfReactants]()]
        formula = rewriteDerivatives(formula, dicts)
        for (rName, rStoich) in reactants
            dus[rName] = dus[rName] * "-" * string(rStoich) * " * (" * formula * ")"
        end
        for (pName, pStoich) in products
            dus[pName] = dus[pName] * "+" * string(pStoich) * " * (" * formula * ")"
        end
    end


    ### Writing to file 
    modelFile = open("./RewrittenModels/" * modelName * ".jl", "w")

    println(modelFile, "# Model name: " * modelName)

    println(modelFile, "# Number of parameters: " * string(length(model[:getListOfParameters]())))
    println(modelFile, "# Number of species: " * string(length(model[:getListOfSpecies]())))

    println(modelFile, "")
    println(modelFile, "### Define independent and dependent variables")
    defineVariables = "@variables t"
    for key in keys(variableDict)
        defineVariables = defineVariables * " " * key * "(t)"
    end
    println(modelFile, defineVariables)

    println(modelFile, "")
    println(modelFile, "### Define variable parameters")
    if length(variableParameterDict) > 0
        defineVariableParameters = "@variables"
        for key in keys(variableParameterDict)
            defineVariableParameters = defineVariableParameters * " " * key * "(t)"
        end
        println(modelFile, defineVariableParameters)
    end
    
    println(modelFile, "")
    println(modelFile, "### Define dummy variable")
    if length(dummyVariableDict) > 0
        defineDummyVariables = "@variables dummyVariable(t)"
        println(modelFile, defineDummyVariables)
    end
    
    println(modelFile, "")
    println(modelFile, "### Define parameters")
    defineParameters = "@parameters"
    for key in keys(parameterDict)
        defineParameters = defineParameters * " " * key
    end
    println(modelFile, defineParameters)

    println(modelFile, "")
    println(modelFile, "### Define constants")
    defineConstants = "@parameters"
    for key in keys(constantsDict)
        defineConstants = defineConstants * " " * key
    end
    println(modelFile, defineConstants)

    println(modelFile, "")
    println(modelFile, "### Define an operator for the differentiation w.r.t. time")
    println(modelFile, "D = Differential(t)")

    println(modelFile, "")
    println(modelFile, "### Function definitions ###")
    if length(funcNameArgFormula) > 0
        for (funcName, (funcArg, funcFormula)) in funcNameArgFormula
            functionAsString = funcName * "(" * funcArg * ") = " * funcFormula 
            functionDefinition = funcName * "(" * funcArg * ")"
            println(modelFile, functionAsString)
            println(modelFile, "@register " * functionDefinition)
        end
    end

    println(modelFile, "")
    println(modelFile, "### Events ###")
    if length(stringOfEvents) > 0
        println(modelFile, "continuous_events = [")
        println(modelFile, stringOfEvents)
        println(modelFile, "]")
    end
    
    println(modelFile, "")
    println(modelFile, "### Derivatives ###")
    println(modelFile, "eqs = [")
    for (sIndex, key) in enumerate(keys(dus))
        if sIndex == 1
            print(modelFile, dus[key])
        else
            print(modelFile, ",\n" * dus[key])
        end
    end
    for key in keys(variableParameterDict)
        print(modelFile, ",\nD(" * key * ") ~ 0")
    end
    if length(dummyVariableDict) > 0
        dummyVariableDerivative = ",\nD(dummyVariable) ~ "
        for key in keys(dummyVariableDict)
            dummyVariableDerivative = dummyVariableDerivative * "+" * key
        end
        println(modelFile, dummyVariableDerivative)
    end
    
    println(modelFile, "]")

    println(modelFile, "")
    if length(stringOfEvents) > 0
        println(modelFile, "@named sys = ODESystem(eqs, t, continuous_events = continuous_events)")
    else
        println(modelFile, "@named sys = ODESystem(eqs)")
    end

    println(modelFile, "")
    println(modelFile, "### Initial species concentrations ###")
    println(modelFile, "initialSpeciesValues = [")
    for (index, (key, value)) in enumerate(variableDict)
        if index == 1
            assignString = key * " => " * value
        else
            assignString = ",\n" * key * " => " * value
        end
        print(modelFile, assignString)
    end
    for (key, value) in variableParameterDict
        assignString = ",\n" * key * " => " * value
        print(modelFile, assignString)
    end
    if length(dummyVariableDict) > 0
        assignString = ",\ndummyVariable => 0.0"
        print(modelFile, assignString)
    end
    println(modelFile, "]")

    println(modelFile, "")
    println(modelFile, "### True parameter values ###")
    println(modelFile, "trueParameterValues = [")
    for (index, (key, value)) in enumerate(parameterDict)
        if index == 1
            assignString = key * " => " * value
        else
            assignString = ",\n" * key * " => " * value
        end
        print(modelFile, assignString)
    end
    println(modelFile, "]")
    println(modelFile, "")
    println(modelFile, "trueConstantsValues = [")
    for (index, (key, value)) in enumerate(constantsDict)
        if index == 1
            assignString = key * " => " * value
        else
            assignString = ",\n" * key * " => " * value
        end
        print(modelFile, assignString)
    end
    println(modelFile, "]")
    
    close(modelFile)

end




writeODEModelToFile()
