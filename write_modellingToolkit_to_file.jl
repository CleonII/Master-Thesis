using PyCall

#ENV["PYTHON"] = "H:\\Program\\Python37\\python.exe"
#Pkg.build("PyCall")

libsbml = pyimport("libsbml")
reader = libsbml.SBMLReader()
document = reader[:readSBML]("./b1.xml")
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
            startPart = i
        end
    end
    numParts += 1
    parts[numParts] = stringToSplit[startPart:end]
    parts = parts[1:numParts]
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
        if parts[1] == "time"
            trigger = "[t ~ " * parts[2] * "]"
            event = "[" * variable * " ~ " * valActive * "]"
        else
            trigger = "[" * parts[1] * " ~ " * parts[2] * "]"
            event = "[" * variable * " ~ " * valActive * "]"
        end
        eventDict[trigger] = event
    elseif "leq" == condition[1:3] || "lt" == condition[1:2]
        if "leq" == condition[1:3]
            strippedCondition = condition[5:end-1]
        else
            strippedCondition = condition[4:end-1]
        end
        parts = splitByComma(strippedCondition)
        if parts[1] == "time"
            trigger = "[t ~ " * parts[2] * "]"
            event = "[" * variable * " ~ " * valInactive * "]"
        else
            trigger = "[" * parts[1] * " ~ " * parts[2] * "]"
            event = "[" * variable * " ~ " * valInactive * "]"
        end
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

function rewritePiecewise(piecewiseString, variable)
    piecewiseString = piecewiseString[11:end-1]
    args = splitByComma(piecewiseString)
    vals = args[1:2:end]
    conds = args[2:2:end]

    pieces = ""
    for (index, condition) in enumerate(conds)
        expression = goToBottom(condition)
        if index == 1
            piece = "if " * expression * "\n" * variable * " = " * vals[index] 
            pieces = piece
        else
            piece = "elseif " * expression * "\n" * variable * " = " * vals[index] 
            pieces = pieces * "\n" * piece
        end
    end
    piece = "else\n" * variable * " = " * vals[end] * "\nend"
    pieces = pieces * "\n" * piece

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
    if parts[1] == "time"
        expression = "[t ~ " * parts[2] * "]"
    else
        expression = "[" * parts[1] * " ~ " * parts[2] * "]"
    end
    return [expression, activates]
end

function getArguments(functionAsString, funcNameArgFormula)
    parts = split(functionAsString, ['(', ')', '/', '+', '-', '*', ' '], keepempty = false)
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

function writeODEModelToFile()
    funcNameArgFormula = Dict()

    ### Define extra functions

    ### Define variables and read initial value
    variableDict = Dict()
    defineVariables = "@variables t"
    for spec in model[:getListOfSpecies]()
        variableDict[spec[:getId]()] = spec[:getInitialAmount]() == 0 ? string(spec[:getInitialConcentration]()) : string(spec[:getInitialAmount]())
        defineVariables = defineVariables * " " * spec[:getId]()
    end
    
    ### Define parameters and read true parameter values
    parameterDict = Dict()
    defineParameters = "@parameters"
    for par in model[:getListOfParameters]()
        parameterDict[par[:getId]()] = string(par[:getValue]())
        defineParameters = defineParameters * " " * par[:getId]()
    end

    ### Define constants and read constants values
    constantsDict = Dict()
    defineConstants = "@parameters"
    for comp in model[:getListOfCompartments]()
        constantsDict[comp[:getId]()] = string(comp[:getSize]())
        defineConstants = defineConstants * " " * comp[:getId]()
    end

    ### Read true parameter values
    trueParameterValues = ""
    for (pIndex, par) in enumerate(model[:getListOfParameters]())
        if pIndex == 1
            trueParameterValues = "[" * string(par[:getValue]())
        else
            trueParameterValues = trueParameterValues * ", " * string(par[:getValue]())
        end
    end
    trueParameterValueString = "trueParameterValues = " * trueParameterValues * "]"

    ### Read constants values
    trueConstantsValues = ""
    for (cIndex, comp) in enumerate(model[:getListOfCompartments]())
        if cIndex == 1
            trueConstantsValues = "[" * string(comp[:getSize]())
        else
            trueConstantsValues = trueConstantsValues * ", " * string(comp[:getSize]())
        end
    end
    trueConstantsValueString = "trueConstantsValues = " * trueConstantsValues * "]"

    ### Read initial values for species
    initialSpeciesValue = ""
    for (sIndex, spec) in enumerate(model[:getListOfSpecies]())
        if sIndex == 1
            initialSpeciesValue = spec[:getInitialAmount]() == 0 ? "[" * string(spec[:getInitialConcentration]()) : "[" * string(spec[:getInitialAmount]())
        else
            initialSpeciesValue = spec[:getInitialAmount]() == 0 ? initialSpeciesValue * ", " * string(spec[:getInitialConcentration]()) : initialSpeciesValue * ", " * string(spec[:getInitialAmount]())
        end
    end
    initialSpeciesValuesString = "initialSpeciesValue = " * initialSpeciesValue * "]"

    ### Define functions
    stringOfFunctions = ""
    stringOfFunctionDefinitions = ""
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
        functionAsString = functionName * args * " = " * functionFormula
        funcNameArgFormula[functionName] = [args[2:end-1], functionFormula]
        if fIndex == 1
            stringOfFunctions = functionAsString
            stringOfFunctionDefinitions = functionName * args
        else
            stringOfFunctions = stringOfFunctions * " # " * functionAsString
            stringOfFunctionDefinitions = stringOfFunctionDefinitions * " # " * functionName * args
        end
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
        if occursin("max(", formula)
            println(formula)
            println("")
        end
        formula = removePowFunctions(formula)
        #formula = rewriteFunctionOfFunction(formula, funcNameArgFormula)
        for (rName, rStoich) in reactants
            dus[rName] = dus[rName] * "-" * string(rStoich) * " * (" * formula * ")"
        end
        for (pName, pStoich) in products
            dus[pName] = dus[pName] * "+" * string(pStoich) * " * (" * formula * ")"
        end
    end


    ### Writing to file 
    modelFile = open("./testMTK.jl", "w")

    println(modelFile, "# Number of parameters: " * string(length(model[:getListOfParameters]())))
    println(modelFile, "# Number of species: " * string(length(model[:getListOfSpecies]())))

    println(modelFile, "")
    println(modelFile, "### Define independent and dependent variables")
    defineVariables = "@variables t"
    for key in keys(variableDict)
        defineVariables = defineVariables * " " * key
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
