function getRuleFormula(rule)
    
    ruleFormula = rule[:getFormula]()
    ruleFormula = replaceWholeWord(ruleFormula, "time", "t") 
    ruleFormula = removePowFunctions(ruleFormula)
    
    return ruleFormula
end


function processAssignmentRule!(modelDict::Dict, ruleFormula::String, ruleVariable::String, baseFunctions)

    # If piecewise occurs in the rule we are looking at a time-based event which is encoded as an 
    # event into the model to ensure proper evaluation of the gradient.
    if occursin("piecewise(", ruleFormula)
        rewritePiecewiseToIfElse(ruleFormula, ruleVariable, modelDict, baseFunctions)
           
    # If the rule does not involve a piecewise expression simply encode it as a function which downsteram 
    # is integrated into the equations. 
    else
        # Extract the parameters and states which make out the rule, further check if rule consists of an 
        # existing function (nesting).
        arguments, includesFunction = getArguments(ruleFormula, modelDict["modelRuleFunctions"], baseFunctions)
        if isempty(arguments)
            modelDict["parameters"][ruleVariable] = ruleFormula
        else
            if includesFunction == true
                ruleFormula = replaceWholeWordDict(ruleFormula, modelDict["modelRuleFunctions"])
            end
            modelDict["modelRuleFunctions"][ruleVariable] = [arguments, ruleFormula]

            # As we hard-code the rule variable into the equation remove it as state or model parameter.
            # TODO : Add this as expression into the model eq. and allow structurally simplify to act on it.
            if ruleVariable in keys(modelDict["states"])
                modelDict["states"] = delete!(modelDict["states"], ruleVariable)
            end
            if ruleVariable in keys(modelDict["parameters"])
                modelDict["parameters"] = delete!(modelDict["parameters"], ruleVariable)
            end
        end        
    end

end


function processRateRule!(modelDict::Dict, ruleFormula::String, ruleVariable::String, baseFunctions)

    # Rewrite rule to function if there are not any piecewise, eles rewrite to formula with ifelse
    if occursin("piecewise(", ruleFormula)
        ruleFormula = rewritePiecewiseToIfElse(ruleFormula, ruleVariable, modelDict, baseFunctions, retFormula=true)
    else                
        arguments, includesFunction = getArguments(ruleFormula, modelDict["modelRuleFunctions"], baseFunctions)
        if arguments != "" && includesFunction == true
            ruleFormula = replaceWholeWordDict(ruleFormula, modelDict["modelRuleFunctions"])
        end
    end

    # Add rate rule as part of model derivatives and remove from parameters dict if rate rule variable 
    # is a parameter  
    if ruleVariable in keys(modelDict["states"])
        modelDict["derivatives"][ruleVariable] = "D(" * ruleVariable * ") ~ " * ruleFormula
    
    elseif ruleVariable in keys(modelDict["nonConstantParameters"])
        modelDict["states"][ruleVariable] = modelDict["nonConstantParameters"][ruleVariable]
        delete!(modelDict["nonConstantParameters"], ruleVariable)
        modelDict["derivatives"][ruleVariable] = "D(" * ruleVariable * ") ~ " * ruleFormula
    
    elseif ruleVariable in keys(parameterDict)
        modelDict["states"][ruleVariable] = parameterDict[ruleVariable]
        delete!(parameterDict, ruleVariable)
        modelDict["derivatives"][ruleVariable] = "D(" * ruleVariable * ") ~ " * ruleFormula
    else
        println("Warning : Cannot find rate rule variable in either model states or parameters")
    end
end


# Rewrites time-dependent ifElse-statements to depend on a boolean variable. This makes it possible to treat piecewise 
# as events, allowing us to properly handle discontinious. Does not rewrite ifElse if the activation criteria depends 
# on a state.                                   
function timeDependentIfElseToBool!(modelDict::Dict)                                  
                                  
    # Extrad picewise 
    for key in keys(modelDict["inputFunctions"])

        formulaWithIfelse = modelDict["inputFunctions"][key]
        indexIfElse = getIndexPiecewise(formulaWithIfelse)
        
        for i in eachindex(indexIfElse)

            ifelseFormula = formulaWithIfelse[indexIfElse[i]][8:end-1]
            activationRule, leftSide, rightSide = split(ifelseFormula, ",")

            # Find inequality 
            iLt = findfirst(x -> x == '<', activationRule) 
            iGt = findfirst(x -> x == '>', activationRule) 
            if isnothing(iGt) && !isnothing(iLt)
                signUsed = "lt"
                splitBy = "<"
            elseif !isnothing(iGt) && isnothing(iLt)
                signUsed = "lt"
                splitBy = "<"
            end
            lhsRule, rhsRule = split(activationRule, string(splitBy))

            # Identify which side of ifelse expression is activated with time 
            timeRight = checkForTime(string(rhsRule))
            timeLeft = checkForTime(string(lhsRule))
            rewriteIfElse = true
            if timeLeft == false && timeLeft == false
                println("Have ifelse statements which does not contain time. Hence we do not 
                        rewrite as event, but rather keep it as an ifelse.")
                rewriteIfElse = false
                continue
            elseif timeLeft == true
                signTime = checkSignTime(string(lhsRule))
                if (signTime == 1 && signUsed == "lt") || (signTime == -1 && signUsed == "gt")
                    sideActivatedWithTime = "right"
                elseif (signTime == 1 && signUsed == "gt") || (signTime == -1 && signUsed == "lt")
                    sideActivatedWithTime = "left"
                end
            elseif timeRight == true
                signTime = checkSignTime(string(rhsRule))
                if (signTime == 1 && signUsed == "lt") || (signTime == -1 && signUsed == "gt")
                    sideActivatedWithTime = "left"
                elseif (signTime == 1 && signUsed == "gt") || (signTime == -1 && signUsed == "lt")
                    sideActivatedWithTime = "right"
                end
            end

            if rewriteIfElse == true
                varName = string(key) * "_bool" * string(i)
                activatedWithTime = sideActivatedWithTime == "left" ? leftSide : rightSide
                deActivatedWithTime = sideActivatedWithTime == "left" ? rightSide : leftSide
                formulaInModel = "((1 - " * varName * ")*" * "(" * deActivatedWithTime *") + " * varName * "*(" * activatedWithTime * "))"
                modelDict["parameters"][varName] = "0.0"
                modelDict["inputFunctions"][key] = replace(modelDict["inputFunctions"][key], formulaWithIfelse[indexIfElse[i]] => formulaInModel)
                modelDict["boolVariables"][varName] = [activationRule, sideActivatedWithTime]
            end
        end
    end
end


function getIndexPiecewise(str::String)
    
    ret = Array{Any, 1}(undef, 0)

    iStart, iEnd = 0, 0
    i = 1
    while i < length(str)

        println("i = $i")

        if !(length(str) > i+6)
            break
        end

        if str[i:(i+5)] == "ifelse"
            iStart = i
            paranthesisLevel = 1
            for j in (i+7):length(str)
                if str[j] == '('
                    paranthesisLevel += 1
                elseif str[j] == ')'
                    paranthesisLevel -= 1
                end
                if paranthesisLevel == 0
                    iEnd = j
                    break
                end
            end
            ret = push!(ret, collect(iStart:iEnd))
            i = iEnd + 1
            continue
        end
        i += 1
    end

    return ret
end

