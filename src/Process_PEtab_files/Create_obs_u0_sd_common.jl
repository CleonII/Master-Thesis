
"""
    peTabFormulaToJulia(formula::String, stateNames, paramData::ParametersInfo, namesParamDyn::Array{String, 1}, namesNonDynParam::Array{String, 1})::String
    Translate a peTab formula (e.g for observable or for sd-parameter) into Julia syntax and output the result 
    as a string.

    State-names, namesParamDyn and paramData are all required to correctly identify states and parameters in the formula.
"""
function peTabFormulaToJulia(formula::String, stateNames, paramData::ParamData, namesParamDyn::Array{String, 1}, namesNonDynParam::Array{String, 1})::String

    # Characters directly translate to Julia and characters that also are assumed to terminate a word (e.g state and parameter)
    charDirectTranslate = ['(', ')', '+', '-', '/', '*', '^'] 
    lenFormula = length(formula)
    formulaJulia = ""
    i = 1
    while i <= lenFormula
        # In case character i of the string can be translated directly 
        if formula[i] in charDirectTranslate
            formulaJulia *= formula[i] * " "
            i += 1

        # In case character i cannot be translated directly (is part of a word)
        else
            # Get word (e.g param, state, math-operation or number)
            word, iNew = getWord(formula, i, charDirectTranslate)
            # Translate word to Julia syntax 
            formulaJulia *= wordToJuliaSyntax(word, stateNames, paramData, namesParamDyn, namesNonDynParam)
            i = iNew

            # Special case where we have multiplication
            if isNumber(word) && i <= lenFormula && isletter(formula[i])
                formulaJulia *= "* "
            end
        end
    end

    return formulaJulia
end


"""
    getWord(str::String, iStart, charListTerm)

    In a string starting from position iStart extract the next "word", which is the longest 
    concurent occurance of characters that are not in the character list with word termination 
    characters. Returns the word and iEnd (the position where the word ends).
    
    For example, if charListTerm = ['(', ')', '+', '-', '/', '*', '^'] abc123 is 
    considered a word but not abc123*.  
"""
function getWord(str::String, iStart::Int, charListTerm::Array{Char, 1})
    
    wordStr = ""
    iEnd = iStart

    # If the first character is a numberic the termination occurs when 
    # the first non-numeric character (not digit or dot .) is reached. 
    isNumericStart = isnumeric(str[iEnd])

    while iEnd <= length(str)
        if !(str[iEnd] in charListTerm) 

            # Parase sciencetific notation for number 
            if isNumericStart == true && str[iEnd] == 'e'
                if length(str) > iEnd && (str[iEnd+1] == '-' || isnumeric(str[iEnd+1]))
                    if str[iEnd+1] == '-'
                        iEnd += 2
                        wordStr *= "e-"
                    else
                        iEnd += 1
                        wordStr *= "e"
                    end

                else
                    break 
                end
            end

            if isNumericStart == true && !(isnumeric(str[iEnd]) || str[iEnd] == '.')
                break
            end
            wordStr *= str[iEnd]
        else
            break 
        end
        iEnd += 1
    end
    # Remove all spaces from the word
    wordStr = replace(wordStr, " " => "")
    return wordStr, iEnd
end


""""
    wordToJuliaSyntax(wordTranslate::String, 
                           stateNames,
                           paramData::ParamData, 
                           namesParamDyn::Array{String, 1})::String

    Translate a word (state, parameter, math-expression or number) into Julia syntax 
    when building Ymod, U0 and Sd functions.

"""
function wordToJuliaSyntax(wordTranslate::String, 
                           stateNames,
                           paramData::ParametersInfo, 
                           namesParamDyn::Array{String, 1}, 
                           namesNonDynParam::Array{String, 1})::String

    # List of mathemathical operations that are accpeted and will be translated 
    # into Julia syntax (t is assumed to be time)
    listOperations = ["exp", "sin", "cos", "t"]

    stateNamesStr = replace.(string.(stateNames), "(t)" => "")
    wordJuliaSyntax = ""
    
    # If wordTranslate is a constant parameter (is not paramDyn - list of parameter to estimate)
    if wordTranslate in string.(paramData.parameterId) && !(wordTranslate in namesParamDyn) && !(wordTranslate in namesNonDynParam)
        # Constant parameters get a _C appended to tell them apart 
        wordJuliaSyntax *= wordTranslate * "_C"
    end

    # If wordTranslate is a dynamic parameters 
    if wordTranslate in namesParamDyn
        wordJuliaSyntax *= wordTranslate
    end

    if wordTranslate in namesNonDynParam
        wordJuliaSyntax *= wordTranslate
    end

    if wordTranslate in stateNamesStr
        wordJuliaSyntax *= wordTranslate
    end

    if isNumber(wordTranslate)
        wordJuliaSyntax *= wordTranslate
    end

    if wordTranslate in listOperations
        wordJuliaSyntax *= listOperations[wordTranslate .== listOperations]
    end

    # If word is one of the observable parameters 
    if length(wordTranslate) >= 19 && wordTranslate[1:19] == "observableParameter"
        wordJuliaSyntax *= wordTranslate
    end

    # If word is one of the noise parameters 
    if length(wordTranslate) >= 14 && wordTranslate[1:14] == "noiseParameter"
        wordJuliaSyntax *= wordTranslate
    end

    if isempty(wordTranslate)
        println("Warning : When creating observation function $wordTranslate could not be processed")
    end

    wordJuliaSyntax *= " "

    return wordJuliaSyntax
end


"""
    getObsParamStr(measurmentFormula::String)::String

    Helper function to extract all observableParameter in the observableFormula in the PeTab-file. 
"""
function getObsParamStr(measurmentFormula::String)::String
    
    # Find all words on the form observableParameter
    obsWords = sort(unique([ match.match for match in eachmatch(r"observableParameter[0-9]_\w+", measurmentFormula) ]))
    obsWordStr = ""
    for i in eachindex(obsWords)
        if i != length(obsWords) 
            obsWordStr *= obsWords[i] * ", "
        else
            obsWordStr *= obsWords[i] 
        end
    end

    return obsWordStr
end


"""
    getNoiseParamStr(sdFormula::String)::String

    Helper function to extract all the noiseParameter in noiseParameter formula in the PeTab file.
"""
function getNoiseParamStr(sdFormula::String)::String
    
    # Find all words on the form observableParameter
    sdWords = [ match.match for match in eachmatch(r"noiseParameter[0-9]_\w+", sdFormula) ]
    sdWordStr = ""
    for i in eachindex(sdWords)
        if i != length(sdWords) 
            sdWordStr *= sdWords[i] * ", "
        else
            sdWordStr *= sdWords[i] 
        end
    end

    return sdWordStr
end


"""
    replaceVariablesWithArrayIndex(formula,stateNames,parameterNames,namesNonDynParam,paramData)::String

    Replaces any state or parameter from formula with their corresponding index in the ODE system 
    Symbolics can return strings without multiplication sign, e.g. 100.0STAT5 instead of 100.0*STAT5 
    so replaceWholeWord cannot be used here
"""
function replaceVariablesWithArrayIndex(formula,stateNames,parameterNames,namesNonDynParam,paramData)::String
    stateNamesStr = replace.(string.(stateNames), "(t)" => "")
    
    for i in eachindex(stateNamesStr)
        formula=replaceWholeWordWithNumberPrefix(formula, stateNamesStr[i], "u["*string(i)*"]")
    end
    for i in eachindex(parameterNames)
        formula=replaceWholeWordWithNumberPrefix(formula, parameterNames[i], "dynPar["*string(i)*"]")
    end
    for i in eachindex(namesNonDynParam)
        formula=replaceWholeWordWithNumberPrefix(formula, namesNonDynParam[i], "nonDynParam["*string(i)*"]")
    end
    for i in eachindex(paramData.parameterID)
        if paramData.shouldEst[i] == false
            formula=replaceWholeWordWithNumberPrefix(formula, paramData.parameterID[i]*"_C", "paramData.paramVal[" * string(i) *"]")
        end
    end
    return formula
end



"""
    replaceExplicitVariableWithRule(formula,stateNames,parameterNames,namesNonDynParam,paramData)::String

    Replace the explicit rule variable with the explicit rule
"""
function replaceExplicitVariableWithRule(formula, modelDict)::String
    for (key,value) in modelDict["modelRuleFunctions"]            
        formula = replaceWholeWord(formula, key, "(" * value[2] * ")")
    end
    return formula
end


"""
replaceWholeWordWithNumberPrefix(formula,from,to)::String
    Replaces variables that can be prefixed with numbers.
    e.g., replaceWholeWordWithNumberPrefix("4STAT5 + 100.0STAT5 + RE*STAT5 + STAT5","STAT5","u[1]") gives
    4u[1] + 100.0u[1] + RE*u[1] + u[1]   
"""
function replaceWholeWordWithNumberPrefix(oldString, replaceFrom, replaceTo)
    replaceFromRegex = Regex("\\b(\\d+\\.?\\d*)*(" * replaceFrom * ")\\b")
    replaceToRegex = SubstitutionString("\\1" * replaceTo )
    newString = replace(oldString, replaceFromRegex => replaceToRegex)
    return newString

end



