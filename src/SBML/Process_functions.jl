# Extracts the argument from a function.
# If a dictionary (with functions) is supplied, will also check if there are nested functions and will 
# include the arguments of these nested functions as arguments of the first function.
# The returned string will only contain unique arguments.
function getArguments(functionAsString, baseFunctions::Array{String, 1})
    parts = split(functionAsString, ['(', ')', '/', '+', '-', '*', ' ', '~', '>', '<', '=', ','], keepempty = false)
    arguments = Dict()
    for part in parts
        if isdigit(part[1])
            nothing
        else
            if (part in values(arguments)) == false && !(part in baseFunctions)
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
function getArguments(functionAsString, dictionary::Dict, baseFunctions::Vector{String})
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
                if (part in values(arguments)) == false && !(part in baseFunctions)
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


# When processing rules into function this function rewrites any function defined rules inside into 
# a function. 
function rewriteFunctionOfFunction(functionAsString, funcNameArgFormula)

    newFunctionString = functionAsString
    for (key,value) in funcNameArgFormula
        varFrom = Regex("(\\b" *key* "\\b)")
        funcFormula = "(" * value[2] * ")"
        newFunctionString = replace(newFunctionString, varFrom => funcFormula)
    end
    return newFunctionString

end


# Will substitute the function definition for the formula given by the model with the arguments 
# when the funciton is a rule that has been rewritten into a rule.
# given to the function used instead of the general arguments in contrast to rewriteFunctionOfFunction. 
# Main goal, inser models formulas when producing the model equations.
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


# For a SBML function extract the function arguments 
function getSBMLFuncArg(mathSBML)::String
    args = "("
    if mathSBML[:getNumChildren]() > 1
        args = args * mathSBML[:getLeftChild]()[:getName]()
        for n in range(1, mathSBML[:getNumChildren]()-1, step = 1)
            arg = mathSBML[:getChild](n)[:getName]()
            if arg !== nothing
                args = args * ", " * arg
            end
        end
        args = args * ")"
    end

    return args
end


function getSBMLFuncFormula(mathSBML, libsbml)

    mathAsString = libsbml[:formulaToString](mathSBML)
        
    # Remove any lambda and only keep the arguments and actual function formulation 
    expressionStart = findfirst('(', mathAsString)+1
    expressionStop = findlast(')', mathAsString)-1
    StrippedMathAsString = mathAsString[expressionStart:expressionStop]
        
    # The actual math formula is given between the last , and end paranthesis 
    functionFormula = lstrip(split(StrippedMathAsString, ',')[end])

    functionFormula = removePowFunctions(functionFormula)

    return functionFormula
end


# Insert the function definitions for newly defined functions. 
function insertFunctionDefinitions(functionAsString, funcNameArgFormula)

    newFunctionString = functionAsString
    for (key,value) in funcNameArgFormula
        varFrom = Regex("(\\b" *key* "\\b)")
        funcFormula = "(" * value[2] * ")"
        newFunctionString = replace(newFunctionString, varFrom => funcFormula)
    end

    return newFunctionString
end