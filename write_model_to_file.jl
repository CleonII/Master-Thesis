using PyCall

#ENV["PYTHON"] = "H:\\Program\\Python37\\python.exe"
#Pkg.build("PyCall")

libsbml = pyimport("libsbml")
reader = libsbml.SBMLReader()
document = reader[:readSBML]("./model_Alkan_SciSignal2018.xml")
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
    if length(values) > length(conditions)
        piece = "else\n" * variable * " = " * values[end] * "\nend"
        pieces = pieces * "\n" * piece
    end
end

function fixPiecewiseFormula(formula)
    parts = split(formula, "piecewise")
    expression = ""
    for (pIndex, part) in enumerate(parts[2:end])
        inPar = 0
        piecewiseStop = 0
        for (i, c) in enumerate(part)
            if c == '('
                inPar += 1
            elseif c == ')'
                inPar -= 1
            end
            if inPar == 0
                piecewiseStop = i
                break
            end
        end
        piecewisePart = part[2:piecewiseStop-1]
        restPart = part[piecewiseStop+1:end]
        arguments = splitByComma(piecewisePart)
        values = arguments[1:2:end]
        conditions = arguments[2:2:end]
        pieces = ""
        for (index, condition) in enumerate(conditions)
            expression = goToBottom(condition)
            if index == 1
                piece = "if " * expression * "; " * values[index] 
                pieces = piece
            else
                piece = "elseif " * expression * "; " * values[index] 
                pieces = pieces * "; " * piece
            end
        end
        if length(values) > length(conditions)
            piece = "else; " * values[end] 
            pieces = pieces * "; " * piece
        end
        if pIndex == 1
            expression = "begin; " * pieces * "; end; end" * restPart
        else
            expression = expression * "begin; " * pieces * "; end" * restPart
        end
        expression = parts[1] * expression
    end
    return expression
end

function writeODEModelToFile()

    modelFile = open("./test.jl", "w")

    println(modelFile, "# Number of parameters: " * string(length(model[:getListOfParameters]())))
    println(modelFile, "# Number of species: " * string(length(model[:getListOfSpecies]())))

    println(modelFile, "")
    println(modelFile, "### Extra functions")
    println(modelFile, "pow(a,b) = a^b")

    println(modelFile, "")
    println(modelFile, "### True parameter values ###")
    trueParameterValues = ""
    for (pIndex, par) in enumerate(model[:getListOfParameters]())
        if pIndex == 1
            trueParameterValues = "[" * string(par[:getValue]())
        else
            trueParameterValues = trueParameterValues * ", " * string(par[:getValue]())
        end
    end
    trueParameterValueString = "trueParameterValues = " * trueParameterValues * "]"
    println(modelFile, trueParameterValueString)

    println(modelFile, "")
    println(modelFile, "### Initial species concentrations ###")
    initialSpeciesValue = ""
    for (sIndex, spec) in enumerate(model[:getListOfSpecies]())
        if sIndex == 1
            initialSpeciesValue = spec[:getInitialAmount]() == 0 ? "[" * string(spec[:getInitialConcentration]()) : "[" * string(spec[:getInitialAmount]())
        else
            initialSpeciesValue = spec[:getInitialAmount]() == 0 ? initialSpeciesValue * ", " * string(spec[:getInitialConcentration]()) : initialSpeciesValue * ", " * string(spec[:getInitialAmount]())
        end
    end
    initialSpeciesValuesString = "initialSpeciesValue = " * initialSpeciesValue * "]"
    println(modelFile, initialSpeciesValuesString)

    println(modelFile, "")
    println(modelFile, "function ODE_model!(du, u, p, t)")

    println(modelFile, "")
    println(modelFile, "### Parameters ###")
    parameters = ""
    for (pIndex, par) in enumerate(model[:getListOfParameters]())
        if pIndex == 1
            parameters = par[:getId]()
        else
            parameters = parameters * ", " * par[:getId]()
        end
    end
    parametersValue = parameters * " = p"
    println(modelFile, parametersValue)

    println(modelFile, "")
    println(modelFile, "### Compartments ###")
    compartmentsSize = [c[:getId]() * " = " * string(c[:getSize]()) for c in model[:getListOfCompartments]()]
    for compSize in compartmentsSize
        println(modelFile, compSize)
    end

    println(modelFile, "")
    println(modelFile, "### Species ###")
    species = ""
    for (sIndex, spec) in enumerate(model[:getListOfSpecies]())
        if sIndex == 1
            species = spec[:getId]()
        else
            species = species * ", " * spec[:getId]()
        end
    end
    speciesValue = species * " = u"
    println(modelFile, speciesValue)

    println(modelFile, "")
    println(modelFile, "### Function definitions ###")
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
            args = args * ") = "
        end
        mathAsString = libsbml[:formulaToString](math)
        expressionStart = findfirst('(', mathAsString)+1
        expressionStop = findlast(')', mathAsString)-1
        StrippedMathAsString = mathAsString[expressionStart:expressionStop]
        functionFormula = lstrip(split(StrippedMathAsString, ',')[end])
        functionAsString = functionName * args * functionFormula
        println(modelFile, functionAsString)
    end

    println(modelFile, "")
    println(modelFile, "### Events ###")
    for event in model[:getListOfEvents]()
        eventName = event[:getName]()
        println(modelFile, "# " * eventName)
        trigger = event[:getTrigger]()
        triggerMath = trigger[:getMath]()
        triggerFormula = goToBottom(libsbml[:formulaToString](triggerMath))
        println(modelFile, "if " * triggerFormula)
        for eventAssignment in event[:getListOfEventAssignments]()
            variableName = eventAssignment[:getVariable]()
            eventMath = eventAssignment[:getMath]()
            eventMathAsString = libsbml[:formulaToString](eventMath)
            functionAsString = variableName * " = " * eventMathAsString
            println(modelFile, functionAsString)
        end
        println(modelFile, "end")
    end

    println(modelFile, "")
    println(modelFile, "### Rules ###")
    for rule in model[:getListOfRules]()
        ruleType = rule[:getElementName]() # type
        ruleVariable = rule[:getVariable]() # variable
        ruleFormula = rule[:getFormula]() # need to add piecewise function (and others)
        if occursin("piecewise", ruleFormula)
            ruleFormulaAsString = rewritePiecewise(ruleFormula, ruleVariable)
        else
            ruleFormulaAsString = ruleVariable * " = " * ruleFormula
        end
        println(modelFile, ruleFormulaAsString)
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
        if occursin("piecewise", formula)
            formula = fixPiecewiseFormula(formula)
        end
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

    println(modelFile, "")
    println(modelFile, "nothing")
    println(modelFile, "end")
    
    close(modelFile)

end




writeODEModelToFile()


##### Functions to add
# piecewise
# geq
# gt
# leq
# lt
# and
# or

