using SBMLToolkit, ModelingToolkit, OrdinaryDiffEq

function xmlToOdeSys(xmlPath)

    mdl = readSBML(xmlPath)
    odesys = ODESystem(mdl)

    # Replace the entries of odesys where the parameter 
    # is calculated from a formula
    for key in keys(odesys.defaults)
        index = replace(string(key) , "(t)"=>"")
        if haskey(mdl.initial_assignments,index)
            SBMLToJulia = SBMLToolkit.interpret_as_num(mdl.initial_assignments[index])
            odesys.defaults[key] = SBMLToJulia
        end
    end

    return odesys

end

function OdeSysToModelDict(odesys, modelName, modelPath)

    ### Extract all information from the ODE System
    indVar = string(odesys.iv)
    statesArray = states(odesys)
    parameterArray = parameters(odesys)
    eqsArray = odesys.eqs
    defaultValues = odesys.defaults

    ### Create all strings related to the variables/states
    variableStr = "ModelingToolkit.@variables " * indVar
    statesArrayStr = "stateArray = ["
    initialStatesStr = "initialSpeciesValues = [" * "\n"
    for i in eachindex(statesArray)
        tmpState = string(statesArray[i])
        tmpStateNoIndVar = replace(tmpState,"(" * indVar * ")"=>"")
        tmpInitVal = string(defaultValues[statesArray[i]])
        variableStr *= " " * tmpState 
        statesArrayStr *= tmpStateNoIndVar * ", "
        initialStatesStr *= "    " * tmpStateNoIndVar * " => " * tmpInitVal * ",\n"
    end
    ### Remove last comma in arrays and end them with ]
    statesArrayStr = statesArrayStr[1:end-2] * "]"
    initialStatesStr = initialStatesStr[1:end-2] * "\n    ]"

    ### Create all strings related to the parameters
    parameterStr = "ModelingToolkit.@parameters"
    parameterArrayStr = "parameterArray = ["
    initialparameterStr = "trueParameterValues = [" * "\n"
    for i in eachindex(parameterArray)
        tmpParameter = string(parameterArray[i])
        tmpInitVal = string(defaultValues[parameterArray[i]])
        parameterStr *= " " * tmpParameter
        parameterArrayStr *= tmpParameter * ", "
        initialparameterStr *= "    " * tmpParameter * " => " * tmpInitVal * ",\n"
    end
    ### Remove last comma in arrays and end them with ]
    parameterArrayStr = parameterArrayStr[1:end-2] * "]"
    initialparameterStr = initialparameterStr[1:end-2] * "\n    ]"

    ### Create all strings related to the differentials and equations
    differentialLHS = "D"
    differentialRHS = "Differential(" * indVar * ")"
    differentialStr = differentialLHS * " = " * differentialRHS
    eqsArrayStr = "eqs = [\n"
    for i in eachindex(eqsArray)
        tmpEq = string(eqsArray[i])
        tmpEq = replace(tmpEq, differentialRHS => differentialLHS)
        tmpEq = replace(tmpEq, "(" * indVar * ")"=>"")
        eqsArrayStr *= "    " * tmpEq * ",\n"
    end
    ### Remove last comma in array and end it with ]
    eqsArrayStr = eqsArrayStr[1:end-2] * "\n    ]"

    ### Store all strings in a model dictionary
    modelDict = Dict()
    modelDict["modelName"] = modelName
    modelDict["modelPath"] = modelPath
    modelDict["numberOfParameters"] = string(length(parameterArray))
    modelDict["numberOfStates"] = string(length(statesArray))
    modelDict["variableStr"] = variableStr
    modelDict["statesArrayStr"] = statesArrayStr
    modelDict["parameterStr"] = parameterStr
    modelDict["parameterArrayStr"] = parameterArrayStr
    modelDict["differentialStr"] = differentialStr
    modelDict["eqsArrayStr"] = eqsArrayStr
    modelDict["initialStatesStr"] = initialStatesStr
    modelDict["initialparameterStr"] = initialparameterStr

    return modelDict

end

function writeToFile(modelDict)

    modelFilePath = joinpath(modelDict["modelPath"], modelDict["modelName"] * ".jl")
    modelFile = open(modelFilePath, "w")
    
    println(modelFile, "# Model name: " * modelDict["modelName"])
    println(modelFile, "# Number of parameters: " * modelDict["numberOfParameters"])
    println(modelFile, "# Number of species: " * modelDict["numberOfStates"])
    println(modelFile, "function getODEModel_" * modelDict["modelName"] * "()")
    println(modelFile, "")
    println(modelFile, "    " * "### Define independent and dependent variables")
    println(modelFile, "    " * modelDict["variableStr"])
    println(modelFile, "")
    println(modelFile, "    " * "### Store dependent variables in array for ODESystem command")
    println(modelFile, "    " * modelDict["statesArrayStr"])
    println(modelFile, "")
    println(modelFile, "    " * "### Define variable parameters")
    println(modelFile, "")
    println(modelFile, "    " * "### Define potential algebraic variables")
    println(modelFile, "")
    println(modelFile, "    " * "### Define parameters")
    println(modelFile, "    " * modelDict["parameterStr"])
    println(modelFile, "")
    println(modelFile, "    " * "### Store parameters in array for ODESystem command")
    println(modelFile, "    " * modelDict["parameterArrayStr"])
    println(modelFile, "")
    println(modelFile, "    " * "### Define an operator for the differentiation w.r.t. time")
    println(modelFile, "    " * modelDict["differentialStr"])
    println(modelFile, "")
    println(modelFile, "    " * "### Continious events ###")    
    println(modelFile, "")
    println(modelFile, "    " * "### Discrete events ###")    
    println(modelFile, "")
    println(modelFile, "    " * "### Derivatives ###")
    println(modelFile, "    " * modelDict["eqsArrayStr"])
    println(modelFile, "")
    println(modelFile, "    " * "@named sys = ODESystem(eqs, t, stateArray, parameterArray)")
    println(modelFile, "")
    println(modelFile, "    " * "### Initial species concentrations ###")
    println(modelFile, "    " * modelDict["initialStatesStr"])
    println(modelFile, "")
    println(modelFile, "    " * "### SBML file parameter values ###")
    println(modelFile, "    " * modelDict["initialparameterStr"])
    println(modelFile, "")
    println(modelFile, "    " * "return sys, initialSpeciesValues, trueParameterValues")
    println(modelFile, "")
    println(modelFile, "end")

    close(modelFile)

end

modelName = "model_Bachmann_MSB2011"
modelPath = joinpath(pwd(), "Intermediate", "PeTab_models", modelName)
xmlPath = joinpath(modelPath, modelName * ".xml")

odesys = xmlToOdeSys(xmlPath)

modelDict = OdeSysToModelDict(odesys, modelName, modelPath)
writeToFile(modelDict)