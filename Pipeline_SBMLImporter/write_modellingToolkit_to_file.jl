


# Splits strings by a given delimiter, but only if the delimiter is not inside a function / parenthesis.
function splitBetween(stringToSplit, delimiter)
    parts = Vector{SubString{String}}(undef, length(stringToSplit))
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
        if stringToSplit[i] == delimiter && inParenthesis == 0
            endPart = i-1
            numParts += 1
            parts[numParts] = stringToSplit[startPart:endPart]
            parts[numParts] = strip(parts[numParts])
            startPart = i+1
        end
    end
    numParts += 1
    parts[numParts] = stringToSplit[startPart:end]
    parts[numParts] = strip(parts[numParts])
    parts = parts[1:numParts]
end

# Finds the ending parenthesis of an argument and returns its position. 
function findEndOffunction(functionAsString, iStart)
    inParenthesis = 1
    i = iStart
    endIndex = i
    while i <= length(functionAsString)
        if functionAsString[i] == '('
            inParenthesis += 1
        elseif functionAsString[i] == ')'
            inParenthesis -= 1
        end
        if inParenthesis == 0
            endIndex = i
            break
        end
        i += 1
    end
    return endIndex
end

# Extracts the argument from a function.
# If a dictionary (with functions) is supplied, will also check if there are nested functions and will 
# include the arguments of these nested functions as arguments of the first function.
# The returned string will only contain unique arguments.
function getArguments(functionAsString, baseFunctions::Dict)
    parts = split(functionAsString, ['(', ')', '/', '+', '-', '*', ' ', '~', '>', '<', '=', ','], keepempty = false)
    arguments = Dict()
    for part in parts
        if isdigit(part[1])
            nothing
        else
            if (part in values(arguments)) == false && ~(part in keys(baseFunctions))
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
function getArguments(functionAsString, dictionary::Dict, baseFunctions::Dict)
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
                if (part in values(arguments)) == false && ~(part in keys(baseFunctions))
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

# replaces a word, "toReplace" in functions with another word, "replacer. 
# Often used to change "time" to "t"
# Makes sure not to change for example "time1" or "shift_time"
function replaceWith(oldString, toReplace, replacer)
    if oldString == toReplace
        return replacer
    end
    newString = oldString
    i = 1
    while i < length(newString)
        indices = findnext(toReplace, newString, i)
        if indices === nothing
            break
        end
        if indices[1] == 1
            if newString[indices[end]+1] in [' ', ',', ')']
                newString = newString[1:indices[1]-1] * replacer * newString[indices[end]+1:end]
            end
        elseif indices[end] == length(newString)
            if newString[indices[1]-1] in [' ', '(', '-', '+']
                newString = newString[1:indices[1]-1] * replacer * newString[indices[end]+1:end]
            end
        else
            if newString[indices[1]-1] in [' ', ',', '(', '-', '+'] && newString[indices[end]+1] in [' ', ')', ',']
                newString = newString[1:indices[1]-1] * replacer * newString[indices[end]+1:end]
            end
        end
        i = indices[1] + 1
    end
    return newString
end

# Goes through a condition for an event written for SBML and creates the appropriate events in ModelingToolkit.
# May call itself to handle more complicated events.
# Makes sure there are no doublett events. 
# Also adds parameters used to a dummy variable so that ModelingToolkit does not think they are unused.
function goToBottomPiecewiseToEvent(condition, variable, valActive, valInactive, eventDict, dicts)   
    if "and" == condition[1:3]
        strippedCondition = condition[5:end-1]
        parts = splitBetween(strippedCondition, ',')
        for part in parts
            goToBottomPiecewiseToEvent(part, variable, valActive, valInactive, eventDict, dicts)
        end
    elseif "or" == condition[1:2]
        strippedCondition = condition[4:end-1]
        parts = splitBetween(strippedCondition, ',')
        for part in parts
            goToBottomPiecewiseToEvent(part, variable, valActive, valInactive, eventDict, dicts)
        end
    elseif "geq" == condition[1:3] || "gt" == condition[1:2]
        if "geq" == condition[1:3]
            strippedCondition = condition[5:end-1]
        else
            strippedCondition = condition[4:end-1]
        end
        # if the trigger or event contains parameters or constant, add them to the dummyVariable to keep defined
        args = split(getArguments(strippedCondition, dicts["baseFunctions"]), ", ")
        for arg in args
            if arg in keys(dicts["parameters"])
                dicts["dummyVariable"][arg] = dicts["parameters"][arg]
            end
        end

        parts = splitBetween(strippedCondition, ',')
        if occursin("time", parts[1])
            parts[1] = replaceWith(parts[1], "time", "t")
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
        args = split(getArguments(strippedCondition, dicts["baseFunctions"]), ", ")
        for arg in args
            if arg in keys(dicts["parameters"])
                dicts["dummyVariable"][arg] = dicts["parameters"][arg]
            end
        end

        parts = splitBetween(strippedCondition, ',')
        if occursin("time", parts[1])
            parts[1] = replaceWith(parts[1], "time", "t")
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
# Calls goToBottomPiecewiseToEvent to create the actual events and then combines them.
function rewritePiecewiseToEvent(piecewiseString, variable, dicts)
    piecewiseString = piecewiseString[11:end-1]

    # if the variable in the event is set as a parameter or constant, set it as a variable and delete the previous form
    if variable in keys(dicts["parameters"])
        dicts["variableParameters"][variable] = dicts["parameters"][variable]
        delete!(dicts["parameters"], variable)
    end

    eventDict = Dict()
    args = splitBetween(piecewiseString, ',')
    vals = args[1:2:end]
    conds = args[2:2:end]

    # if any values are a/a set of parameter(s) or constant(s), add it/them to the dummyVariable to keep defined, also insert function definitions as needed
    for (vIndex, val) in enumerate(vals)
        valArgs = getArguments(val, dicts["functions"], dicts["baseFunctions"])
        if valArgs[2]
            vals[vIndex] = rewriteFunctionOfFunction(val, dicts["functions"])
        end
        for valArg in split(valArgs[1], ", ")
            if valArg in keys(dicts["parameters"])
                dicts["dummyVariable"][valArg] = dicts["parameters"][valArg]
            end
        end
    end    

    for (cIndex, condition) in enumerate(conds)
        valActive = vals[cIndex]
        valInactive = vals[end]
        goToBottomPiecewiseToEvent(condition, variable, valActive, valInactive, eventDict, dicts)
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
    parts = splitBetween(strippedFormula, ',')
    if occursin("time", parts[1])
        parts[1] = replaceWith(parts[1], "time", "t")
    end
    expression = "[" * parts[1] * " ~ " * parts[2] * "]"
    return expression
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
            funcFormula = "(" * funcNameArgFormula[part][2] * ")"
            i = 1
            while i < length(newFunctionString)
                indices = findnext(part, newFunctionString, i)
                if indices[1] == 1 && indices[end] == length(newFunctionString)
                    newFunctionString = newFunctionString[1:indices[1]-1] * funcFormula * newFunctionString[indices[end]+1:end]
                    break
                end
                if indices[1] == 1
                    if newFunctionString[indices[end]+1] in [' ', ')']
                        newFunctionString = newFunctionString[1:indices[1]-1] * funcFormula * newFunctionString[indices[end]+1:end]
                        break
                    end
                elseif indices[end] == length(newFunctionString)
                    if newFunctionString[indices[1]-1] in [' ', '(']
                        newFunctionString = newFunctionString[1:indices[1]-1] * funcFormula * newFunctionString[indices[end]+1:end]
                        break
                    end
                else
                    if newFunctionString[indices[1]-1] in [' ', '('] && newFunctionString[indices[end]+1] in [' ', ')']
                        newFunctionString = newFunctionString[1:indices[1]-1] * funcFormula * newFunctionString[indices[end]+1:end]
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
function insertModelDefineFunctions(functionsAsString, modelFuncNameArgFormula, baseFunctions)
    newFunctionsAsString = functionsAsString
    possibleModelFunctions = split(getArguments(functionsAsString, baseFunctions), ", ")
    for possibleModelFunction in possibleModelFunctions 
        if possibleModelFunction in keys(modelFuncNameArgFormula)
            modelFuncArguments, modelFuncFormula = modelFuncNameArgFormula[possibleModelFunction]
            modelFuncArguments = split(modelFuncArguments, ", ")
            i = 1
            while i <= length(newFunctionsAsString)
                next = findnext(possibleModelFunction*"(", newFunctionsAsString, i)
                if next === nothing
                    break
                end
                startIndex = next[1]
                iStart = next[end] + 1
                endIndex = findEndOffunction(newFunctionsAsString, iStart)
                modelFunction = newFunctionsAsString[startIndex:endIndex]
                argumentString = modelFunction[length(possibleModelFunction)+2:end-1]
                arguments = split(argumentString, ", ", keepempty = false)
                newFunction = modelFuncFormula
                for (aIndex, arg) in enumerate(arguments)
                    j = 1
                    while j <= length(newFunction)
                        indices = findnext(modelFuncArguments[aIndex], newFunction, j)
                        if indices === nothing
                            break
                        end
                        if indices[1] == 1 && indices[end] == length(newFunction)
                            if occursin(" ", arg)
                                newFunction = newFunction[1:indices[1]-1] * "(" * arg * ")" * newFunction[indices[end]+1:end]
                            else
                                newFunction = newFunction[1:indices[1]-1] * arg * newFunction[indices[end]+1:end]
                            end
                        elseif indices[1] == 1
                            if newFunction[indices[end]+1] in [' ', ')']
                                if occursin(" ", arg)
                                    newFunction = newFunction[1:indices[1]-1] * "(" * arg * ")" * newFunction[indices[end]+1:end]
                                else
                                    newFunction = newFunction[1:indices[1]-1] * arg * newFunction[indices[end]+1:end]
                                end
                            end
                        elseif indices[end] == length(newFunction)
                            if newFunction[indices[1]-1] in [' ', '(']
                                if occursin(" ", arg)
                                    newFunction = newFunction[1:indices[1]-1] * "(" * arg * ")" * newFunction[indices[end]+1:end]
                                else
                                    newFunction = newFunction[1:indices[1]-1] * arg * newFunction[indices[end]+1:end]
                                end
                            end
                        else
                            if newFunction[indices[1]-1] in [' ', '('] && newFunction[indices[end]+1] in [' ', ')']
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
                i = startIndex + 1
            end
        end
    end
    return newFunctionsAsString
end

# Rewrites pow(a,b) into a^b, which Julia can handle
function removePowFunctions(functionAsString)
    newFunctionAsString = functionAsString
    i = 1
    while i <= length(newFunctionAsString)
        next = findnext("pow(", newFunctionAsString, i)
        if next === nothing
            break
        end
        startIndex = next[1]
        iStart = next[end] + 1
        endIndex = findEndOffunction(newFunctionAsString, iStart)
        powFunction = newFunctionAsString[startIndex:endIndex]
        powArguments = powFunction[5:end-1]
        parts = splitBetween(powArguments, ',')
        newPowFunction = "(" * parts[1] * ")^(" * parts[2] * ")"
        newFunctionAsString = replace(newFunctionAsString, powFunction => newPowFunction, count = 1)
        i = startIndex + 1
    end
    return newFunctionAsString
end

# Handles conditions of piecewise functions in derivatives and rewrites them into a propper format.
# may call itself to handle more complicated functions.
# returns both a trigger and an event to be combined in rewritePiecewiseToFunction
function goToBottomPiecewiseToFunction(condition, valActive, valInactive)
    if valInactive == ""
        event = valActive
    else
        event = "((" * valActive * ") - (" * valInactive * "))"
    end
    if "and" == condition[1:3]
        strippedCondition = condition[5:end-1]
        parts = splitBetween(strippedCondition, ',')
        combinedTrigger = ""
        for part in parts
            trigger,  = goToBottomPiecewiseToFunction(part, valActive, valInactive)
            if combinedTrigger == ""
                combinedTrigger = trigger
            else
                combinedTrigger = combinedTrigger * " * " * trigger
            end
            trigger = combinedTrigger
        end
    elseif "or" == condition[1:2]
        strippedCondition = condition[4:end-1]
        parts = splitBetween(strippedCondition, ',')
        combinedTrigger = ""
        for (pIndex, part) in enumerate(parts)
            trigger,  = goToBottomPiecewiseToFunction(part, valActive, valInactive)
            if pIndex == 1
                combinedTrigger = "(" * trigger
            else
                combinedTrigger = combinedTrigger * " + " * trigger
            end
            combinedTrigger = combinedTrigger * ")"
        end
        trigger = combinedTrigger
    elseif "geq" == condition[1:3] || "gt" == condition[1:2]
        if "geq" == condition[1:3]
            strippedCondition = condition[5:end-1]
            operator = " >= "
        else
            strippedCondition = condition[4:end-1]
            operator = " > "
        end
        parts = splitBetween(strippedCondition, ',')
        if occursin("time", parts[1])
            parts[1] = replaceWith(parts[1], "time", "t")
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
        parts = splitBetween(strippedCondition, ',')
        if occursin("time", parts[1])
            parts[1] = replaceWith(parts[1], "time", "t")
        end

        trigger = "(" * parts[1] * operator * parts[2] * ")"

    end
    return trigger, event
end

# Handel piecewise functions in derivatives and rewrites them into a propper format.
# Calls goToBottomPiecewiseToFunction to handle each condition. 
function rewritePiecewiseToFunction(functionAsString)
    newFunctionAsString = functionAsString
    i = 1
    while i <= length(newFunctionAsString)
        next = findnext("piecewise(", newFunctionAsString, i)
        if next === nothing
            break
        end
        startIndex = next[1]
        iStart = next[end] + 1
        endIndex = findEndOffunction(newFunctionAsString, iStart)
        piecewiseFunction = newFunctionAsString[startIndex:endIndex]
        piecewiseString = piecewiseFunction[11:end-1]
        args = splitBetween(piecewiseString, ',')
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
            trigger, event = goToBottomPiecewiseToFunction(condition, valActive, valInactive)
            if expression == ""
                expression = trigger * " * " * event
            else
                expression = "(" * expression * " + " * trigger * " * " * event * ")"
            end
        end
        newFunctionAsString = replace(newFunctionAsString, piecewiseFunction => expression, count = 1)
        i = startIndex + 1
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
            funcFormula = "(" * funcNameArgFormula[part][2] * ")"
            i = 1
            while i < length(newFunctionString)
                indices = findnext(part, newFunctionString, i)
                if indices === nothing
                    break
                end
                if indices[1] == 1 && indices[end] == length(newFunctionString)
                    newFunctionString = newFunctionString[1:indices[1]-1] * funcFormula * newFunctionString[indices[end]+1:end]
                elseif indices[1] == 1
                    if newFunctionString[indices[end]+1] in [' ', ')']
                        newFunctionString = newFunctionString[1:indices[1]-1] * funcFormula * newFunctionString[indices[end]+1:end]
                    end
                elseif indices[end] == length(newFunctionString)
                    if newFunctionString[indices[1]-1] in [' ', '(', '-', '+']
                        newFunctionString = newFunctionString[1:indices[1]-1] * funcFormula * newFunctionString[indices[end]+1:end]
                    end
                else
                    if newFunctionString[indices[1]-1] in [' ', '(', '-', '+'] && newFunctionString[indices[end]+1] in [' ', ')']
                        newFunctionString = newFunctionString[1:indices[1]-1] * funcFormula * newFunctionString[indices[end]+1:end]
                    end
                end
                i = indices[1] + 1
            end
        end
    end
    return newFunctionString
end

# Rewrites derivatives using insertModelDefineFunctions, removePowFunctions, rewritePiecewiseToFunction and insertFunctionDefinitions
function rewriteDerivatives(derivativeAsString, dicts)
    newDerivativeAsString = derivativeAsString
    newDerivativeAsString = insertModelDefineFunctions(newDerivativeAsString, dicts["modelFunctions"], dicts["baseFunctions"])
    if occursin("pow(", newDerivativeAsString)
        newDerivativeAsString = removePowFunctions(newDerivativeAsString)
    end
    if occursin("piecewise(", newDerivativeAsString)
        newDerivativeAsString = rewritePiecewiseToFunction(newDerivativeAsString)
    end
    newDerivativeAsString = insertFunctionDefinitions(newDerivativeAsString, dicts["functions"])
    return newDerivativeAsString
end

function writeODEModelToFile(libsbml, model, modelName, path, useData, wrapped; experimentalConditions = [], parameterBounds = [])
    dicts = Dict()
    variableDict = Dict()
    dicts["variables"] = variableDict
    parameterDict = Dict()
    dicts["parameters"] = parameterDict
    variableParameterDict = Dict()
    dicts["variableParameters"] = variableParameterDict
    dummyVariableDict = Dict()
    dicts["dummyVariable"] = dummyVariableDict

    if useData
        constantParameterDict = Dict()
        dicts["constantParameter"] = constantParameterDict
    end

    # for functions defined by the model, these functions will be rewritten
    modelFuncNameArgFormula = Dict()
    dicts["modelFunctions"] = modelFuncNameArgFormula
    
    # for functions rewritten from rules etc.
    funcNameArgFormula = Dict()
    dicts["functions"] = funcNameArgFormula
    dus = Dict()
    dicts["derivatives"] = dus

    baseFunctions = Dict()
    dicts["baseFunctions"] = baseFunctions
    baseFunctions["exp"] = ""
    baseFunctions["log"] = ""
    baseFunctions["log2"] = ""
    baseFunctions["log10"] = ""
    baseFunctions["sin"] = ""
    baseFunctions["cos"] = ""
    baseFunctions["tan"] = ""
    baseFunctions["pi"] = ""

    ### Define extra functions

    ### Define variables and read initial value
    for spec in model[:getListOfSpecies]()
        variableDict[spec[:getId]()] = spec[:getInitialAmount]() == 0 ? spec[:getInitialConcentration]() === NaN ? 
            string(spec[:getInitialAmount]()) : string(spec[:getInitialConcentration]()) : string(spec[:getInitialAmount]())
    end
    ### Defining derivatives
    for spec in model[:getListOfSpecies]()
        dus[spec[:getId]()] = "D(" * spec[:getId]() * ") ~ "
    end
    
    ### Define parameters and read true parameter values
    for par in model[:getListOfParameters]()
        parameterDict[par[:getId]()] = string(par[:getValue]())
    end

    ### Define compartments and read compartments values
    for comp in model[:getListOfCompartments]()
        parameterDict[comp[:getId]()] = string(comp[:getSize]())
    end

    ### Read functions defined by the model
    for functionDefinition in model[:getListOfFunctionDefinitions]()
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
        # if the trigger contains parameters, add them to the dummy variables
        triggerArgs = split(getArguments(triggerFormula[2:end-1], baseFunctions), ", ")
        for triggerArg in triggerArgs
            if triggerArg in keys(dicts["parameters"])
                dicts["dummyVariable"][triggerArg] = dicts["parameters"][triggerArg]
            end
        end
        eventAsString = ""
        for (eaIndex, eventAssignment) in enumerate(event[:getListOfEventAssignments]())
            variableName = eventAssignment[:getVariable]()
            # if the variable in the event is not set as a variable, make it so and remove it as a parameter or constant
            if variableName in keys(parameterDict)
                variableParameterDict[variableName] = parameterDict[variableName]
                delete!(parameterDict, variableName)
            end

            eventMath = eventAssignment[:getMath]()
            eventMathAsString = libsbml[:formulaToString](eventMath)
            # if the event math contains parameters, add them to the dummy variable so that they stay defined 
            eventMathArgs = split(getArguments(eventMathAsString, baseFunctions), ", ")
            for eventMathArg in eventMathArgs
                if eventMathArg in keys(parameterDict)
                    dummyVariableDict[eventMathArg] = parameterDict[eventMathArg]
                end
            end

            # Add the event 
            if eaIndex == 1
                eventAsString = "[" * variableName * " ~ " * eventMathAsString
            else
                eventAsString = eventAsString * ", " * variableName * " ~ " * eventMathAsString
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
        if ruleType == "assignmentRule" # TODO: fix for functions with too many arguments
            ruleVariable = rule[:getVariable]() # variable
            ruleFormula = rule[:getFormula]() 
            ruleFormula = replaceWith(ruleFormula, "time", "t") 
            ruleFormula = removePowFunctions(ruleFormula)
            if occursin("piecewise(", ruleFormula)
                next = findfirst("piecewise(", ruleFormula)
                if next[1] == 1 && findEndOffunction(ruleFormula, next[end]+1) == length(ruleFormula) && findnext("piecewise(", ruleFormula, 2) === nothing
                    # adds rule as event
                    localEventString = rewritePiecewiseToEvent(ruleFormula, ruleVariable, dicts)
                    if length(stringOfEvents) == 0
                        stringOfEvents = localEventString
                    else
                        stringOfEvents = stringOfEvents * ", " * localEventString
                    end
                else
                    # adds rule as function
                    ruleFormula = rewritePiecewiseToFunction(ruleFormula)
                    arguments, includesFunction = getArguments(ruleFormula, funcNameArgFormula, baseFunctions)
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
                arguments, includesFunction = getArguments(ruleFormula, funcNameArgFormula, baseFunctions)
                if arguments == ""
                    parameterDict[ruleVariable] = ruleFormula
                else
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
            println("Has algebraicRule!")
        elseif ruleType == "rateRule"
            ruleVariable = rule[:getVariable]() # variable
            ruleFormula = rule[:getFormula]() 
            ruleFormula = replaceWith(ruleFormula, "time", "t") 
            ruleFormula = removePowFunctions(ruleFormula)
            if occursin("piecewise(", ruleFormula)
                # set as derivative
                ruleFormula = rewritePiecewiseToFunction(ruleFormula)
                arguments, includesFunction = getArguments(ruleFormula, funcNameArgFormula, baseFunctions)
                if includesFunction
                    ruleFormula = rewriteFunctionOfFunction(ruleFormula, funcNameArgFormula)
                end
            else
                # set as derivative
                arguments, includesFunction = getArguments(ruleFormula, funcNameArgFormula, baseFunctions)
                if arguments == ""

                else
                    if includesFunction
                        ruleFormula = rewriteFunctionOfFunction(ruleFormula, funcNameArgFormula)
                    end
                end
            end
            if ruleVariable in keys(variableDict)
                dus[ruleVariable] = "D(" * ruleVariable * ") ~ " * ruleFormula
            elseif ruleVariable in keys(variableParameterDict)
                variableDict[ruleVariable] = variableParameterDict[ruleVariable]
                delete!(variableParameterDict, ruleVariable)
                dus[ruleVariable] = "D(" * ruleVariable * ") ~ " * ruleFormula
            else
                if ruleVariable in keys(parameterDict)
                    variableDict[ruleVariable] = parameterDict[ruleVariable]
                    delete!(parameterDict, ruleVariable)
                    dus[ruleVariable] = "D(" * ruleVariable * ") ~ " * ruleFormula
                else
                    variableDict[ruleVariable] = "0"
                    dus[ruleVariable] = "D(" * ruleVariable * ") ~ " * ruleFormula
                end
            end
        end
    end

    ### Implement Initial assignments 
    # Positioned after rules since some assignments may include functions
    initallyAssignedVariable = Dict{String, String}()
    initallyAssignedParameter = Dict{String, String}()
    for initAssign in model[:getListOfInitialAssignments]()
        assignName = initAssign[:getId]()
        assignMath = initAssign[:getMath]()
        assignFormula = libsbml[:formulaToString](assignMath)
        assignFormula = rewriteDerivatives(assignFormula, dicts)
        if assignName in keys(variableDict)
            variableDict[assignName] = assignFormula
            initallyAssignedVariable[assignName] = "variableDict"
        elseif assignName in keys(variableParameterDict)
            variableParameterDict[assignName] = assignFormula
            initallyAssignedVariable[assignName] = "variableParameterDict"
        elseif assignName in keys(parameterDict)
            parameterDict[assignName] = assignFormula
            initallyAssignedParameter[assignName] = "parameterDict"
        else
            println("Error: could not find assigned variable/parameter")
        end
        for part in split(getArguments(assignFormula, baseFunctions), ", ")
            if part in keys(parameterDict)
                dummyVariableDict[part] = parameterDict[part]
            end
        end
    end

    while true
        nestedVariables = false
        for (variable, dictName) in initallyAssignedVariable
            if dictName == "variableDict"
                variableValue = variableDict[variable]
                args = split(getArguments(variableValue, baseFunctions))
                for arg in args
                    if arg in keys(variableDict)
                        nestedVariables = true
                        variableValue = replaceWith(variableValue, arg, variableDict[arg])
                    elseif arg in keys(variableParameterDict)
                        nestedVariables = true
                        variableValue = replaceWith(variableValue, arg, variableParameterDict[arg])
                    end
                end
                variableDict[variable] = variableValue
            else
                variableValue = variableParameterDict[variable]
                args = split(getArguments(variableValue, baseFunctions))
                for arg in args
                    if arg in keys(variableDict)
                        nestedVariables = true
                        variableValue = replaceWith(variableValue, arg, variableDict[arg])
                    elseif arg in keys(variableParameterDict)
                        nestedVariables = true
                        variableValue = replaceWith(variableValue, arg, variableParameterDict[arg])
                    end
                end
                variableParameterDict[variable] = variableValue
            end
        end
        nestedVariables || break
    end

    while true
        nestedParameter = false
        for (parameter, dictName) in initallyAssignedParameter
            parameterValue = parameterDict[parameter]
            args = split(getArguments(parameterValue, baseFunctions))
            for arg in args
                if arg in keys(parameterDict)
                    nestedParameter = true
                    parameterValue = replaceWith(parameterValue, arg, parameterDict[arg])
                end
            end
            parameterDict[parameter] = parameterValue
        end
        nestedParameter || break
    end

    if useData
        # sort out the parameters that are assumed to be known 

        optParameters = parameterBounds[:,:parameterId]
        inputParameters = names(experimentalConditions)

        for key in keys(parameterDict)
            if key ∉ inputParameters && key ∉ optParameters 
                constantParameterDict[key] = parameterDict[key]
                delete!(parameterDict, key)
            end
        end
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

    isInODESys = falses(length(parameterDict))
    for du in values(dus)
        for (i, pars) in enumerate(keys(parameterDict))
            if replaceWith(du, pars, "") !== du
                isInODESys[i] = true
            end
        end
    end
    for (i, pars) in enumerate(keys(parameterDict))
        if !isInODESys[i]
            dummyVariableDict[pars] = parameterDict[pars]
        end
    end

    ### Writing to file 
    modelFile = open(path * "/" * modelName * ".jl", "w")

    println(modelFile, "# Model name: " * modelName)

    println(modelFile, "# Number of parameters: " * string(length(model[:getListOfParameters]())))
    println(modelFile, "# Number of species: " * string(length(model[:getListOfSpecies]())))

    if wrapped 
        println(modelFile, "function getODEModel_" * modelName * "()")
    end

    if useData

        println(modelFile, "")
        println(modelFile, "    ### Define constant parameters")
        for key in keys(constantParameterDict)
            println(modelFile, "    " * key * " = " * constantParameterDict[key])
        end

    end

    println(modelFile, "")
    println(modelFile, "    ### Define independent and dependent variables")
    defineVariables = "    ModelingToolkit.@variables t"
    for key in keys(variableDict)
        defineVariables = defineVariables * " " * key * "(t)"
    end
    println(modelFile, defineVariables)

    println(modelFile, "")
    println(modelFile, "    ### Define variable parameters")
    if length(variableParameterDict) > 0
        defineVariableParameters = "    ModelingToolkit.@variables"
        for key in keys(variableParameterDict)
            defineVariableParameters = defineVariableParameters * " " * key * "(t)"
        end
        println(modelFile, defineVariableParameters)
    end
    
    println(modelFile, "")
    println(modelFile, "    ### Define dummy variable")
    if length(dummyVariableDict) > 0
        defineDummyVariables = "    ModelingToolkit.@variables dummyVariable(t)"
        println(modelFile, defineDummyVariables)
    end
    
    println(modelFile, "")
    println(modelFile, "    ### Define parameters")
    defineParameters = "    ModelingToolkit.@parameters"
    for key in keys(parameterDict)
        defineParameters = defineParameters * " " * key
    end
    println(modelFile, defineParameters)

    println(modelFile, "")
    println(modelFile, "    ### Define an operator for the differentiation w.r.t. time")
    println(modelFile, "    D = Differential(t)")

    println(modelFile, "")
    println(modelFile, "    ### Events ###")
    if length(stringOfEvents) > 0
        println(modelFile, "    continuous_events = [")
        println(modelFile, "    " * stringOfEvents)
        println(modelFile, "    ]")
    end
    
    println(modelFile, "")
    println(modelFile, "    ### Derivatives ###")
    println(modelFile, "    eqs = [")
    for (sIndex, key) in enumerate(keys(variableDict))
        if sIndex == 1
            print(modelFile, "    " * dus[key])
        else
            print(modelFile, ",\n    " * dus[key])
        end
    end
    for key in keys(variableParameterDict)
        print(modelFile, ",\n    D(" * key * ") ~ 0")
    end
    if length(dummyVariableDict) > 0
        dummyVariableDerivative = ",\n    D(dummyVariable) ~ "
        for key in keys(dummyVariableDict)
            dummyVariableDerivative = dummyVariableDerivative * "+" * key
        end
        println(modelFile, dummyVariableDerivative)
    end
    
    println(modelFile, "    ]")

    println(modelFile, "")
    if length(stringOfEvents) > 0
        println(modelFile, "    @named sys = ODESystem(eqs, t, continuous_events = continuous_events)")
    else
        println(modelFile, "    @named sys = ODESystem(eqs)")
    end

    println(modelFile, "")
    println(modelFile, "    ### Initial species concentrations ###")
    println(modelFile, "    initialSpeciesValues = [")
    for (index, (key, value)) in enumerate(variableDict)
        if index == 1
            assignString = "    " * key * " => " * value
        else
            assignString = ",\n    " * key * " => " * value
        end
        print(modelFile, assignString)
    end
    for (key, value) in variableParameterDict
        assignString = ",\n    " * key * " => " * value
        print(modelFile, assignString)
    end
    if length(dummyVariableDict) > 0
        assignString = ",\n    dummyVariable => 0.0"
        print(modelFile, assignString)
    end
    println(modelFile, "]")

    println(modelFile, "")
    println(modelFile, "    ### True parameter values ###")
    println(modelFile, "    trueParameterValues = [")
    for (index, (key, value)) in enumerate(parameterDict)
        if index == 1
            assignString = "    " * key * " => " * value
        else
            assignString = ",\n    " * key * " => " * value
        end
        print(modelFile, assignString)
    end
    println(modelFile, "]")
    println(modelFile, "")

    if wrapped
        println(modelFile, "    return sys, initialSpeciesValues, trueParameterValues")

        println(modelFile, "")

        println(modelFile, "end")
    end
    
    close(modelFile)
end
