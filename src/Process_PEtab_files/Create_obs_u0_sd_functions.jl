# TODO: Refactor yMod and SD functions to avoid redundant code 


"""
    createFileYmodSdU0(modelName::String, 
                       dirModel::String, 
                       odeSys::ODESystem, 
                       stateMap,
                       modelDict::Dict)

    For a PeTab model with name modelName with all PeTab-files in dirModel and associated 
    ModellingToolkit ODESystem (with its stateMap) build a file containing a functions for 
    i) computing the observable model value (yMod) ii) compute the initial value u0 (by using the 
    stateMap) and iii) computing the standard error (sd) for each observableFormula in the 
    observables PeTab file.

    Note - The produced Julia file will go via the JIT-compiler.
"""
function createFileYmodSdU0(modelName::String, 
                            dirModel::String, 
                            odeSys::ODESystem, 
                            stateMap,
                            modelDict::Dict)
                            
    parameterNames = parameters(odeSys)
    stateNames = states(odeSys)

    # Read data on experimental conditions and parameter values 
    experimentalConditionsFile, measurementDataFile, parameterBoundsFile, observablesDataFile = readDataFiles(dirModel, readObs=true)
    paramData = processParameterData(parameterBoundsFile) # Model data in convient structure 
    measurementData = processMeasurementData(measurementDataFile, observablesDataFile) 
    
    # Indices for mapping parameter-estimation vector to dynamic, observable and sd parameters correctly when calculating cost
    paramIndices = computeIndicesθ(paramData, measurementData, odeSys, experimentalConditionsFile)
    
    createYmodFunction(modelName, dirModel, stateNames, paramData, paramIndices.θ_dynamicNames, paramIndices.θ_nonDynamicNames, observablesDataFile, modelDict)
    println("Done with Ymod function")
    println("")
    
    createU0Function(modelName, dirModel, paramData, string.(parameterNames), stateMap, inPlace=true)
    println("Done with u0 function in-place")
    println("")

    createU0Function(modelName, dirModel, paramData, string.(parameterNames), stateMap, inPlace=false)
    println("Done with u0 function not in-place")
    println("")

    createSdFunction(modelName, dirModel, paramData, stateNames, paramIndices.θ_dynamicNames, paramIndices.θ_nonDynamicNames, observablesDataFile, modelDict)
    println("Done with sd function")
    println("")
end


"""
    createTopOfFun(stateNames, 
        paramData::ParameterInfo, 
        namesParamODEProb::Array{String, 1}, 
        namesNonDynParam::Array{String, 1},
        observablesData::DataFrame,
        modelDict::Dict,
        obsFun::Bool)
    Extracts all variables needed for the functions. 
    Also adds them as variables for Symbolics.jl
"""
function createTopOfFun(stateNames, 
                        paramData::ParameterInfo, 
                        namesParamODEProb::Array{String, 1}, 
                        namesNonDynParam::Array{String, 1},
                        observablesData::DataFrame,
                        modelDict::Dict,
                        obsFun::Bool)
    # Extract names of model states 
    stateNamesShort = replace.(string.(stateNames), "(t)" => "")
    stateStr = "\n\t"
    for i in eachindex(stateNamesShort)
        stateStr *= stateNamesShort[i] * ", "
    end
    stateStr = stateStr[1:end-2]
    stateStr *= "= u \n"
    
    # Extract name of dynamic parameter 
    paramDynStr = "\t"
    for i in eachindex(namesParamODEProb)
        paramDynStr *= namesParamODEProb[i] * ", "
    end
    paramDynStr = paramDynStr[1:end-2]
    paramDynStr *= " = dynPar \n"

    paramNonDynStr = ""
    # Extract name of non-dynamic parameter
    if !isempty(namesNonDynParam)
        paramNonDynStr = "\t"
        for i in eachindex(namesNonDynParam)
            paramNonDynStr *= namesNonDynParam[i] * ", "
        end
        paramNonDynStr = paramNonDynStr[1:end-2]
        paramNonDynStr *= " = nonDynParam \n"
    end
    
    # Extract constant parameters. To avoid cluttering the function only constant parameters that are used to 
    # to compute yMod are written to file.
    paramConstStr = ""
    for i in eachindex(paramData.parameterID)
        if paramData.shouldEst[i] == false
            paramConstStr *= "\t" * paramData.parameterID[i] * "_C = paramData.paramVal[" * string(i) *"] \n" 
        end
    end
    paramConstStr *= "\n"

    # Create namesExplicitRules, if length(modelDict["modelRuleFunctions"]) = 0 this becomes a String[].
    namesExplicitRules = Array{String, 1}(undef,length(modelDict["modelRuleFunctions"]))
    # Extract the keys from the modelRuleFunctions
    for (index, key) in enumerate(keys(modelDict["modelRuleFunctions"]))
        namesExplicitRules[index] = key
    end
       
    # Extracts the explicit rules that have keys among the observables.
    explicitRules = ""
    tmpFormulaConCat = ""
    
    observableIDs = String.(observablesData[!, "observableId"])
    for i in eachindex(observableIDs)
        # Concatenates all formulae and afterwards checks if any explicit rule variable is used there.
        if obsFun == true
            tmpFormula = filter(x -> !isspace(x), String(observablesData[i, "observableFormula"]))
            tmpFormulaConCat *= tmpFormula * "\n"
        else
            tmpFormula = filter(x -> !isspace(x), String(observablesData[i, "noiseFormula"]))
            tmpFormulaConCat *= tmpFormula * "\n"
        end
    end

    for (key,value) in modelDict["modelRuleFunctions"]
        if occursin(Regex("\\b" * key * "\\b"), tmpFormulaConCat)
            explicitRules *= "\t" * key * " = " * value[2] * "\n"
        end
    end

    if length(explicitRules)>0
        explicitRules *= "\n"
    end

    return stateStr, paramDynStr, paramNonDynStr, paramConstStr, explicitRules, namesExplicitRules
end

"""
    createYmodFunction(modelName::String, 
                       dirModel::String, 
                       stateNames, 
                       paramData::ParameterInfo, 
                       namesParamDyn::Array{String, 1}, 
                       observablesData::DataFrame,
                       modelDict::Dict)

    For modelName create a function for computing yMod by translating the observablesData
    PeTab-file into Julia syntax. 

    To correctly create the function the state-names, names of dynamic parameters to estiamte 
    (namesDynParam) and PeTab parameter-file (to get constant parameters) data are needed. 

    The modelDict is used to define explicit rules to the Ymod file.
"""
function createYmodFunction(modelName::String, 
                            dirModel::String, 
                            stateNames, 
                            paramData::ParameterInfo, 
                            namesParamDyn::Array{String, 1}, 
                            namesNonDynParam::Array{String, 1},
                            observablesData::DataFrame,
                            modelDict::Dict)

    io = open(dirModel * "/" * modelName * "ObsSdU0.jl", "w")
    stateStr, paramDynStr, paramNonDynStr, paramConstStr, explicitRules, namesExplicitRules = createTopOfFun(stateNames, paramData, namesParamDyn, namesNonDynParam, observablesData, modelDict, true)

    # Write the formula of each observable to file
    observableIDs = String.(observablesData[!, "observableId"])
    strObserveble = ""
    for i in eachindex(observableIDs)
        # Each observebleID falls below its own if-statement 
        strObserveble *= "\tif observableId == " * "\"" * observableIDs[i] * "\"" * " \n"
        tmpFormula = filter(x -> !isspace(x), String(observablesData[i, "observableFormula"]))

        # Extract observable parameters 
        obsParam = getObsParamStr(tmpFormula)
        if !isempty(obsParam)
            strObserveble *= "\t\t" * obsParam * " = getObsOrSdParam(obsPar, mapObsParam)\n" 
        end 

        # Translate the formula for the observable to Julia syntax 
        juliaFormula = peTabFormulaToJulia(tmpFormula, stateNames, paramData, namesParamDyn, namesNonDynParam, namesExplicitRules)
        strObserveble *= "\t\t" * "return " * juliaFormula * "\n"
        strObserveble *= "\tend\n\n"
    end

    write(io, "function evalYmod(u, t, dynPar, obsPar, nonDynParam, paramData, observableId, mapObsParam) \n")
    write(io, stateStr)
    write(io, paramDynStr)
    write(io, paramNonDynStr)
    write(io, paramConstStr)
    write(io, explicitRules)
    write(io, strObserveble)
    strClose = "end\n\n"
    write(io, strClose)
    close(io)
end


"""
    createU0Function(modelName::String, 
                         dirModel::String, 
                         paramData::ParameterInfo, 
                         namesParameter::Array{String, 1}, 
                         stateMap)

    For modelName create a function for computing initial value by translating the stateMap 
    into Julia syntax.

    To correctly create the function the name of all parameters, paramData (to get constant parameters)
    are required.
"""
function createU0Function(modelName::String, 
                          dirModel::String, 
                          paramData::ParameterInfo, 
                          namesParameter::Array{String, 1}, 
                          stateMap;
                          inPlace::Bool=true)

    io = open(dirModel * "/" * modelName * "ObsSdU0.jl", "a")

    if inPlace == true
        write(io, "function evalU0!(u0Vec, paramVec) \n\n")
    else
        write(io, "function evalU0(paramVec) \n\n")
    end

    # Extract all model parameters (constant and none-constant) and write to file 
    paramDynStr = "\t"
    for i in eachindex(namesParameter)
        paramDynStr *= namesParameter[i] * ", "
    end
    paramDynStr = paramDynStr[1:end-2]
    paramDynStr *= " = paramVec \n\n"
    write(io, paramDynStr)    

    # Write the formula for each initial condition to file 
    stateNames = [replace.(string.(stateMap[i].first), "(t)" => "") for i in eachindex(stateMap)]
    stateExpWrite = ""
    for i in eachindex(stateMap)
        stateName = stateNames[i]
        stateExp = replace(string(stateMap[i].second), " " => "")
        stateFormula = peTabFormulaToJulia(stateExp, stateNames, paramData, namesParameter, String[], String[])
        stateExpWrite *= "\t" * stateName * " = " * stateFormula * "\n"
    end
    write(io, stateExpWrite * "\n")

    # Ensure the states in correct order are written to u0 
    # In place version where we mutate stateVec 
    if inPlace == true
        stateStr = "\tu0Vec .= "
        for i in eachindex(stateNames)
            stateStr *= stateNames[i] * ", "
        end
        stateStr = stateStr[1:end-2]
        write(io, stateStr)

    # Where we return the entire initial value vector 
    elseif inPlace == false
        stateStr = "\t return ["
        for i in eachindex(stateNames)
            stateStr *= stateNames[i] * ", "
        end
        stateStr = stateStr[1:end-2]
        stateStr *= "]"
        write(io, stateStr)
    end

    strClose = "\nend\n\n"
    write(io, strClose)
    close(io)
end


"""
    createSdFunction(modelName::String, 
                          dirModel::String, 
                          paramData::ParameterInfo, 
                          stateNames, 
                          namesParamDyn::Array{String, 1}, 
                          observablesData::DataFrame,
                          modelDict::Dict)

    For modelName create a function for computing the standard deviation by translating the observablesData
    PeTab-file into Julia syntax. 

    To correctly create the function the state-names, names of dynamic parameters to estiamte 
    (namesDynParam) and PeTab parameter-file (to get constant parameters) data are needed. 
"""
function createSdFunction(modelName::String, 
                          dirModel::String, 
                          paramData::ParameterInfo, 
                          stateNames, 
                          namesParamDyn::Array{String, 1}, 
                          namesNonDynParam::Array{String, 1},
                          observablesData::DataFrame,
                          modelDict::Dict)


    io = open(dirModel * "/" * modelName * "ObsSdU0.jl", "a")
    stateStr, paramDynStr, paramNonDynStr, paramConstStr, explicitRules, namesExplicitRules = createTopOfFun(stateNames, paramData, namesParamDyn, namesNonDynParam, observablesData, modelDict, false)

    # Write the formula for standard deviations to file
    observableIDs = String.(observablesData[!, "observableId"])
    strObserveble = ""
    for i in eachindex(observableIDs)
        # Each observebleID falls below its own if-statement 
        strObserveble *= "\tif observableId == " * "\"" * observableIDs[i] * "\"" * " \n"
        tmpFormula = filter(x -> !isspace(x), String(observablesData[i, "noiseFormula"]))

        # Extract noise parameters 
        noiseParam = getNoiseParamStr(tmpFormula)
        if !isempty(noiseParam)
            strObserveble *= "\t\t" * noiseParam * " = getObsOrSdParam(sdPar, mapSdParam)\n" 
        end 

        juliaFormula = peTabFormulaToJulia(tmpFormula, stateNames, paramData, namesParamDyn, namesNonDynParam, namesExplicitRules)
        strObserveble *= "\t\t" * "return " * juliaFormula * "\n"
        strObserveble *= "\tend\n\n"
    end

    write(io, "function evalSd!(u, t, sdPar, dynPar, nonDynParam, paramData, observableId, mapSdParam) \n")
    write(io, stateStr)
    write(io, paramDynStr)
    write(io, paramNonDynStr)
    write(io, paramConstStr)
    write(io, explicitRules)
    write(io, strObserveble)
    strClose = "end"
    write(io, strClose)
    close(io)
end