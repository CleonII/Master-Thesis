using PyCall

#ENV["PYTHON"] = "H:\\Program\\Python37\\python.exe"
#Pkg.build("PyCall")

libsbml = pyimport("libsbml")
reader = libsbml.SBMLReader()
modelName = "model_Beer_MolBioSystems2014"
document = reader[:readSBML]("./SBML/Benchmark/" * modelName * ".xml")
model = document[:getModel]() # Get the model

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

function replaceTime(timeString)
    if timeString == "time"
        return "t"
    end
    newTimeString = timeString
    i = 1
    while i < length(timeString)
        indices = findnext("time", timeString, i)
        if indices === nothing
            break
        end
        if indices[1] == 1
            if timeString[indices[end]+1] == ' '
                newTimeString = replace(newTimeString, "time " => "t ", count = 1)
            end
        elseif indices[end] == length(timeString)
            if timeString[indices[1]-1] == ' '
                newTimeString = replace(newTimeString, " time" => " t", count = 1)
            end
        else
            if timeString[indices[1]-1] == ' ' && timeString[indices[end]+1] == ' '
                newTimeString = replace(newTimeString, " time " => " t ", count = 1)
            end
        end
        i = indices[end]
    end
    return newTimeString
end

function goToBottomEvent(condition, variable, valActive, valInactive, eventDict)
    if "and" == condition[1:3]
        strippedCondition = condition[5:end-1]
        parts = splitByComma(strippedCondition)
        for part in parts
            goToBottomEvent(part, variable, valActive, valInactive, eventDict)
        end
    elseif "or" == condition[1:2]
        strippedCondition = condition[4:end-1]
        parts = splitByComma(strippedCondition)
        for part in parts
            goToBottomEvent(part, variable, valActive, valInactive, eventDict)
        end
    elseif "geq" == condition[1:3] || "gt" == condition[1:2]
        if "geq" == condition[1:3]
            strippedCondition = condition[5:end-1]
        else
            strippedCondition = condition[4:end-1]
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

function piecewiseToEvent(piecewiseString, variable)
    piecewiseString = piecewiseString[11:end-1]
    eventDict = Dict()
    args = splitByComma(piecewiseString)
    vals = args[1:2:end]
    conds = args[2:2:end]

    for (cIndex, condition) in enumerate(conds)
        valActive = vals[cIndex]
        valInactive = vals[end]
        goToBottomEvent(condition, variable, valActive, valInactive, eventDict)
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

function asTrigger(triggerFormula)
    activates = false
    if "geq" == triggerFormula[1:3]
        activates = true
        strippedFormula = triggerFormula[5:end-1]
    elseif "gt" == triggerFormula[1:2]
        activates = true
        strippedFormula = triggerFormula[4:end-1]
    elseif "leq" == triggerFormula[1:3]
        activates = false
        strippedFormula = triggerFormula[5:end-1]
    elseif "lt" == triggerFormula[1:2]
        activates = false
        strippedFormula = triggerFormula[4:end-1]
    end
    parts = splitByComma(strippedFormula)
    if occursin("time", parts[1])
        parts[1] = replaceTime(parts[1])
    end
    expression = "[" * parts[1] * " ~ " * parts[2] * "]"
    return [expression, activates]
end

function getArguments(functionAsString, funcNameArgFormula)
    parts = split(functionAsString, ['(', ')', '/', '+', '-', '*', ' ', '~', ','], keepempty = false)
    existingFunctions = keys(funcNameArgFormula)
    includesFunction = false
    arguments = Dict()
    for part in parts
        if isdigit(part[1])
            nothing
        else
            if part in existingFunctions
                includesFunction = true
                funcArgs = funcNameArgFormula[part][1]
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

function rewriteFunctionOfFunction(functionAsString, funcNameArgFormula)
    parts = split(functionAsString, ['(', ')', '/', '+', '-', '*', ' ', ','], keepempty = false)
    existingFunctions = keys(funcNameArgFormula)
    newFunctionString = functionAsString
    for part in parts
        if isdigit(part[1])
            nothing
        else
            if part in existingFunctions
                funcArgs = funcNameArgFormula[part][1]
                funcDef = part * "(" * funcArgs * ")"
                newFunctionString = replace(newFunctionString, part => funcDef)
            end
        end
    end
    return newFunctionString
end

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

function rewriteMaxFunctions(functionAsString)
    newFunctionAsString = functionAsString
    inParenthesis = 0
    started = 0
    startIndex = 1
    endIndex = 1
    i = 1
    while i <= length(functionAsString)
        if started == 0
            next = findnext("max(", functionAsString, i)
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
                newPowFunction = "0.5*(" * parts[1] * " + " * parts[2] * " + abs(" * parts[1] * " - " * parts[2] * "))"
                newFunctionAsString = replace(newFunctionAsString, powFunction => newPowFunction, count = 1)
                startIndex = i
                started = 0
            end
            i += 1
        end
    end
    return newFunctionAsString
end

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
                        expression = expression * " + " * trigger * " * " * event
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

function insertFunctionDefinitions(functionAsString, funcNameArgFormula)
    parts = split(functionAsString, ['(', ')', '/', '+', '-', '*', ' ', ','], keepempty = false)
    existingFunctions = keys(funcNameArgFormula)
    newFunctionString = functionAsString
    for part in parts
        if isdigit(part[1])
            nothing
        else
            if part in existingFunctions
                funcArgs = funcNameArgFormula[part][1]
                funcDef = part * "(" * funcArgs * ")"
                newFunctionString = replace(newFunctionString, part => funcDef)
            end
        end
    end
    return newFunctionString
end

function rewriteDerivatives(derivativeAsString, funcNameArgFormula)
    newDerivativeAsString = derivativeAsString
    if occursin("pow(", newDerivativeAsString)
        newDerivativeAsString = removePowFunctions(newDerivativeAsString)
    end
    if occursin("max(", newDerivativeAsString)
        newDerivativeAsString = rewriteMaxFunctions(newDerivativeAsString)
    end
    if occursin("piecewise(", newDerivativeAsString)
        newDerivativeAsString = rewritePiecewiseInDerivative(newDerivativeAsString)
    end
    newDerivativeAsString = insertFunctionDefinitions(newDerivativeAsString, funcNameArgFormula)
    return newDerivativeAsString
end

function writeODEModelToFile()
    funcNameArgFormula = Dict()

    ### Define extra functions

    ### Define variables and read initial value
    variableDict = Dict()
    for spec in model[:getListOfSpecies]()
        variableDict[spec[:getId]()] = spec[:getInitialAmount]() == 0 ? string(spec[:getInitialConcentration]()) : string(spec[:getInitialAmount]())
    end
    
    ### Define parameters and read true parameter values
    parameterDict = Dict()
    for par in model[:getListOfParameters]()
        parameterDict[par[:getId]()] = string(par[:getValue]())
    end

    ### Define constants and read constants values
    constantsDict = Dict()
    for comp in model[:getListOfCompartments]()
        constantsDict[comp[:getId]()] = string(comp[:getSize]())
    end

    ### Implement Initial assignments 
    for initAssign in model[:getListOfInitialAssignments]()
        assignName = initAssign[:getId]()
        if assignName in keys(variableDict)
            assignMath = initAssign[:getMath]()
            assignFormula = libsbml[:formulaToString](assignMath)
            variableDict[assignName] = assignFormula
        elseif assignName in keys(parameterDict)
            assignMath = initAssign[:getMath]()
            assignFormula = libsbml[:formulaToString](assignMath)
            parameterDict[assignName] = assignFormula
        elseif assignName in keys(constantsDict)
            assignMath = initAssign[:getMath]()
            assignFormula = libsbml[:formulaToString](assignMath)
            constantsDict[assignName] = assignFormula
        else
            println("Error: could not find assigned variable/parameter")
        end
        
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
        funcNameArgFormula[functionName] = [args[2:end-1], functionFormula]
    end

    ### Define events
    stringOfEvents = ""
    for (eIndex, event) in enumerate(model[:getListOfEvents]())
        eventName = event[:getName]()
        trigger = event[:getTrigger]()
        triggerMath = trigger[:getMath]()
        triggerFormula, activates = asTrigger(libsbml[:formulaToString](triggerMath))
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
                # adds rule as event
                localEventString = piecewiseToEvent(ruleFormula, ruleVariable)
                if length(stringOfEvents) == 0
                    stringOfEvents = localEventString
                else
                    stringOfEvents = stringOfEvents * ", " * localEventString
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

    ### Defining derivatives
    dus = Dict()
    for spec in model[:getListOfSpecies]()
        dus[spec[:getId]()] = "D(" * spec[:getId]() * ") ~ "
    end

    reactions = [(r, r[:getKineticLaw]()[:getFormula]()) for r in model[:getListOfReactions]()]
    for (reac, formula) in reactions
        products = [(p[:species], p[:getStoichiometry]()) for p in reac[:getListOfProducts]()]
        reactants = [(r[:species], r[:getStoichiometry]()) for r in reac[:getListOfReactants]()]
        formula = rewriteDerivatives(formula, funcNameArgFormula)
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
    println(modelFile, "### Events ###")
    if length(stringOfEvents) > 0
        println(modelFile, "continuous_events = [")
        println(modelFile, stringOfEvents)
        println(modelFile, "]")
    end

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
    println(modelFile, "### Derivatives ###")
    println(modelFile, "eqs = [")
    for (sIndex, key) in enumerate(keys(dus))
        if sIndex == 1
            print(modelFile, dus[key])
        else
            print(modelFile, ",\n" * dus[key])
        end
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
