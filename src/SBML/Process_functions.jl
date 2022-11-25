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


# Replaces a word, "replaceFrom" in functions with another word, "replaceTo". 
# Often used to change "time" to "t"
# Makes sure not to change for example "time1" or "shift_time"
function replaceWholeWord(oldString, replaceFrom, replaceTo)
    
    replaceFromRegex = Regex("(\\b" * replaceFrom * "\\b)")
    newString = replace(oldString, replaceFromRegex => replaceTo)
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


# Replaces words in oldString given a dictionary replaceDict.
# In the Dict, the key  is the word to replace and the second
# value is the value to replace with.
# Makes sure to only change whole words.
function replaceWholeWordDict(oldString, replaceDict)

    newString = oldString
    for (key,value) in replaceDict
        newString = replaceWholeWord(newString, key, "(" * value[2] * ")")
    end
    return newString

end


# Substitutes the function with the formula given by the model, but replaces
# the names of the variables in the formula with the input variable names.
# e.g. If fun(a) = a^2 then "constant * fun(b)" will be rewritten as 
# "constant * b^2"
# Main goal, insert model formulas when producing the model equations.
function replaceFunctionWithFormula(functionAsString, funcNameArgFormula)

    newFunctionsAsString = functionAsString
    
    for (key,value) in funcNameArgFormula
        # The string we wish to insert when the correct replacement has been made.
        replaceStr = value[2]
    
        # Finds the old input arguments, removes spaces and puts them in a list
        replaceFrom = split(replace(value[1]," "=>""),",")
        
        # Finds the first time the function is present in the string
        # and captures the name and input arguments used there.
        # Assumes that the function is not called with different 
        # arguments in different parts of the same line.
        varFrom = Regex("\\b(" * key * ")+\\((.*?)\\)")
        mtc = match(varFrom, functionAsString)
        
        if !isnothing(mtc)
            # capture contains two elements 
            # 1) the captured function name 
            # 2) the captured arguments of the function
            captures = mtc.captures
            # Takes the new input arguments removes spaces and puts them in a list
            replaceTo = split(replace(captures[2]," "=>""),",")
            
            # Replace each variable used in the formula with the 
            # variable name used as input for the function.
            for ind in eachindex(replaceTo)
                replaceStr = replaceWholeWord(replaceStr, replaceFrom[ind], replaceTo[ind])
            end

            # Replace function(input) with formula where each variable in formula has the correct name.
            newFunctionsAsString = replace(newFunctionsAsString, captures[1] * "(" * captures[2] * ")" => "(" * replaceStr * ")")

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