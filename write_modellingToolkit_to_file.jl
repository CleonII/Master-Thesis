using PyCall

#ENV["PYTHON"] = "H:\\Program\\Python37\\python.exe"
#Pkg.build("PyCall")

libsbml = pyimport("libsbml")
reader = libsbml.SBMLReader()
document = reader[:readSBML]("./b3.xml")
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
            i += 2
            startPart = i
        end
    end
    numParts += 1
    parts[numParts] = stringToSplit[startPart:end]
    parts = parts[1:numParts]
end

function goToBottom(condition)
    if cmp("and", condition[1:3]) == 0
        strippedCondition = condition[5:end-1]
        parts = splitByComma(strippedCondition)
        expression = ""
        for (index, part) in enumerate(parts)
            expr = goToBottom(part)
            if index == 1
                expression = expr
            else
                expression = expression * " && " * expr
            end
        end
    elseif cmp("or", condition[1:2]) == 0
        strippedCondition = condition[4:end-1]
        parts = splitByComma(strippedCondition)
        expression = ""
        for (index, part) in enumerate(parts)
            expr = goToBottom(part)
            if index == 1
                expression = expr
            else
                expression = expression * " || " * expr
            end
        end
    elseif cmp("geq", condition[1:3]) == 0
        strippedCondition = condition[5:end-1]
        parts = splitByComma(strippedCondition)
        if cmp(parts[1], "time") == 0
            expression = "t >= " * parts[2]
        else
            expression = parts[1] * " >= " * parts[2]
        end
    elseif cmp("gt", condition[1:2]) == 0
        strippedCondition = condition[4:end-1]
        parts = splitByComma(strippedCondition)
        if cmp(parts[1], "time") == 0
            expression = "t > " * parts[2]
        else
            expression = parts[1] * " > " * parts[2]
        end
    elseif cmp("leq", condition[1:3]) == 0
        strippedCondition = condition[5:end-1]
        parts = splitByComma(strippedCondition)
        if cmp(parts[1], "time") == 0
            expression = "t <= " * parts[2]
        else
            expression = parts[1] * " <= " * parts[2]
        end
    elseif cmp("lt", condition[1:2]) == 0
        strippedCondition = condition[4:end-1]
        parts = splitByComma(strippedCondition)
        if cmp(parts[1], "time") == 0
            expression = "t < " * parts[2]
        else
            expression = parts[1] * " < " * parts[2]
        end
    end
    return expression
end

function piecewiseToEvent(piecewiseString, variable)
    piecewiseString = piecewiseString[11:end-1]
    arguments = splitByComma(piecewiseString)
    values = arguments[1:2:end]
    conditions = arguments[2:2:end]

    pieces = ""
    for (index, condition) in enumerate(conditions)
        expression = goToBottom(condition)
        if index == 1
            piece = "if " * expression * "\n" * variable * " = " * values[index] 
            pieces = piece
        else
            piece = "elseif " * expression * "\n" * variable * " = " * values[index] 
            pieces = pieces * "\n" * piece
        end
    end
    piece = "else\n" * variable * " = " * values[end] * "\nend"
    pieces = pieces * "\n" * piece

end

function rewritePiecewise(piecewiseString, variable)
    piecewiseString = piecewiseString[11:end-1]
    arguments = splitByComma(piecewiseString)
    values = arguments[1:2:end]
    conditions = arguments[2:2:end]

    pieces = ""
    for (index, condition) in enumerate(conditions)
        expression = goToBottom(condition)
        if index == 1
            piece = "if " * expression * "\n" * variable * " = " * values[index] 
            pieces = piece
        else
            piece = "elseif " * expression * "\n" * variable * " = " * values[index] 
            pieces = pieces * "\n" * piece
        end
    end
    piece = "else\n" * variable * " = " * values[end] * "\nend"
    pieces = pieces * "\n" * piece

end

function asTrigger(triggerFormula)
    if cmp("geq", triggerFormula[1:3]) == 0
        strippedFormula = triggerFormula[5:end-1]
    elseif cmp("gt", triggerFormula[1:2]) == 0
        strippedFormula = triggerFormula[4:end-1]
    elseif cmp("leq", triggerFormula[1:3]) == 0
        strippedFormula = triggerFormula[5:end-1]
    elseif cmp("lt", triggerFormula[1:2]) == 0
        strippedFormula = triggerFormula[4:end-1]
    end
    parts = splitByComma(strippedFormula)
    if cmp(parts[1], "time") == 0
        expression = "[t ~ " * parts[2] * "]"
    else
        expression = "[" * parts[1] * " ~ " * parts[2] * "]"
    end
    return expression
end

function writeODEModelToFile()
    numFunctions = 0

    ### Define variables
    defineVariables = "@variables t"
    for spec in model[:getListOfSpecies]()
        defineVariables = defineVariables * " " * spec[:getId]()
    end
    
    ### Define parameters
    defineParameters = "@parameters"
    for par in model[:getListOfParameters]()
        defineParameters = defineParameters * " " * par[:getId]()
    end

    ### Define constants
    defineConstants = "@parameters"
    for comp in model[:getListOfCompartments]()
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
        functionAsString = functionName * args * " = " * functionFormula
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
        triggerFormula = asTrigger(libsbml[:formulaToString](triggerMath))
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
    println(modelFile, "")
    println(modelFile, "### Rules ###")
    for rule in model[:getListOfRules]()
        ruleType = rule[:getElementName]() # type
        if ruleType == "assignmentRule"
            ruleVariable = rule[:getVariable]() # variable
            ruleFormula = rule[:getFormula]() # need to add piecewise function (and others)
            if occursin("piecewise", ruleFormula)
                ruleFormulaAsString = rewritePiecewise(ruleFormula, ruleVariable)
            else
                ruleFormulaAsString = ruleVariable * " = " * ruleFormula
            end
            println(modelFile, ruleFormulaAsString)
        elseif cmp(ruleType, "algebraicRule") == 0
            # TODO
        elseif cmp(ruleType, "rateRule") == 0
            # TODO
        end
    end

    println(modelFile, "")
    println(modelFile, "### Derivatives ###")
    dus = Dict()
    for (index, spec) in enumerate(model[:getListOfSpecies]())
        dus[spec[:getId]()] = "du[" * string(index) * "] = "
    end

    reactions = [(r, r[:getKineticLaw]()[:getFormula]()) for r in model[:getListOfReactions]()]
    for (reac, formula) in reactions
        products = [(p[:species], p[:getStoichiometry]()) for p in reac[:getListOfProducts]()]
        reactants = [(r[:species], r[:getStoichiometry]()) for r in reac[:getListOfReactants]()]
        for (rName, rStoich) in reactants
            dus[rName] = dus[rName] * "-" * string(rStoich) * " * (" * formula * ")"
        end
        for (pName, pStoich) in products
            dus[pName] = dus[pName] * "+" * string(pStoich) * " * (" * formula * ")"
        end
    end

    for spec in model[:getListOfSpecies]()
        println(modelFile, dus[spec[:getId]()])
    end

    ### Writing to file 
    modelFile = open("./testMTK.jl", "w")

    println(modelFile, "# Number of parameters: " * string(length(model[:getListOfParameters]())))
    println(modelFile, "# Number of species: " * string(length(model[:getListOfSpecies]())))

    println(modelFile, "")
    println(modelFile, "### Extra functions")
    println(modelFile, "pow(a,b) = a^b")
    println(modelFile, "@register pow(a,b)")

    println(modelFile, "")
    println(modelFile, "### Define independent and dependent variables")
    println(modelFile, defineVariables)
    
    println(modelFile, "")
    println(modelFile, "### Define parameters")
    println(modelFile, defineParameters)

    println(modelFile, "")
    println(modelFile, "### Define constants")
    println(modelFile, defineConstants)

    println(modelFile, "")
    println(modelFile, "### Define an operator for the differentiation w.r.t. time")
    println(modelFile, "D = Differential(t)")

    println(modelFile, "")
    println(modelFile, "### True parameter values ###")
    println(modelFile, trueParameterValueString)
    println(modelFile, trueConstantsValueString)

    println(modelFile, "")
    println(modelFile, "### Initial species concentrations ###")
    println(modelFile, initialSpeciesValuesString)

    println(modelFile, "")
    println(modelFile, "### Function definitions ###")
    fParts = split(stringOfFunctions, " # ")
    fdParts = split(stringOfFunctionDefinitions, " # ")
    for (functionAsString, functionDefinition) in zip(fParts, fdParts)
        println(modelFile, functionAsString)
        println(modelFile, "@register " * functionDefinition)
    end
    
    println(modelFile, "")
    println(modelFile, "### Events ###")
    println(modelFile, "continuous_events = [")
    println(modelFile, stringOfEvents)
    println(modelFile, "]")





    println(modelFile, "")
    println(modelFile, "@named sys = ODESystem(eqs)")

    println(modelFile, "")
    println(modelFile, "nothing")
    println(modelFile, "end")
    
    close(modelFile)

end




writeODEModelToFile()
