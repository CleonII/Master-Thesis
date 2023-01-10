
"""
    createFileDYmodSdU0(modelName::String, 
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
function createFileDYmodSdU0(modelName::String, 
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
    paramIndices = getIndicesParam(paramData, measurementData, odeSys, experimentalConditionsFile)
    
    createDYmodFunction(modelName, dirModel, stateNames, paramData, string.(parameterNames), paramIndices.namesNonDynParam, observablesDataFile, modelDict)
    println("Done with DYmod function")
    println("")
    
    createDSdFunction(modelName, dirModel, paramData, stateNames, string.(parameterNames), paramIndices.namesNonDynParam, observablesDataFile, modelDict)
    println("Done with Dsd function")
    println("")

end

"""
    createTopOfDFun(modelName::String, 
                       dirModel::String, 
                       stateNames, 
                       paramData::ParamData, 
                       namesParamODEProb::Array{String, 1}, 
                       observablesData::DataFrame,
                       modelDict::Dict)
    Extracts all variables needed for the functions. 
    Also adds them as variables for Symbolics.jl
"""
function createTopOfDFun(stateNames, 
                        paramData::ParamData, 
                        namesParamODEProb::Array{String, 1}, 
                        namesNonDynParam::Array{String, 1},
                        observablesData::DataFrame,
                        modelDict::Dict)

    variablesStr = "@variables "
    parameterArray = namesParamODEProb
    statesArray = Array{String, 1}(undef,length(stateNames))

    # Extract names of model states 
    stateNamesShort = replace.(string.(stateNames), "(t)" => "")
    stateStr = "\n\t"
    for i in eachindex(stateNamesShort)
        stateStr *= stateNamesShort[i] * ", "
        variablesStr *= stateNamesShort[i] * ", "
        statesArray[i] = stateNamesShort[i]
    end
    stateStr = stateStr[1:end-2]
    stateStr *= "= u \n"
    
    # Extract name of dynamic parameter 
    paramDynStr = "\t"
    for i in eachindex(namesParamODEProb)
        paramDynStr *= namesParamODEProb[i] * ", "
        variablesStr *= namesParamODEProb[i] * ", "
    end
    paramDynStr = paramDynStr[1:end-2]
    paramDynStr *= " = p \n"

    paramNonDynStr = ""
    # Extract name of non-dynamic parameter
    if !isempty(namesNonDynParam)
        paramNonDynStr = "\t"
        for i in eachindex(namesNonDynParam)
            paramNonDynStr *= namesNonDynParam[i] * ", "
            variablesStr *= namesNonDynParam[i] * ", "
        end
        paramNonDynStr = paramNonDynStr[1:end-2]
        paramNonDynStr *= " = nonDynParam \n"
    end
    
    # Create namesExplicitRules, if length(modelDict["modelRuleFunctions"]) = 0 this becomes a String[].
    namesExplicitRules = Array{String, 1}(undef,length(modelDict["modelRuleFunctions"]))
    # Extract the keys from the modelRuleFunctions
    for (index, key) in enumerate(keys(modelDict["modelRuleFunctions"]))
        namesExplicitRules[index] = key
        variablesStr *= key * ", "
    end
       
    # Extracts the explicit rules that have keys among the observables.
    # also extracts all observable- and noise-parameters
    explicitRules = ""
    strObserveble = ""
    
    observableIDs = String.(observablesData[!, "observableId"])
    for i in eachindex(observableIDs)
        tmpFormula = filter(x -> !isspace(x), String(observablesData[i, "observableFormula"]))
        strObserveble *= tmpFormula * "\n"
        # Extract observable parameters 
        obsParam = getObsParamStr(tmpFormula)
        if !isempty(obsParam)
            variablesStr *= obsParam * ", "   
        end 

        tmpFormula = filter(x -> !isspace(x), String(observablesData[i, "noiseFormula"]))
        strObserveble *= tmpFormula * "\n"
        noiseParam = getNoiseParamStr(tmpFormula)
        if !isempty(noiseParam)
            variablesStr *= noiseParam * ", "   
        end 

    end

    # Remove last "," and Run @variables ... string
    variablesStr = variablesStr[1:end-2]
    eval(Meta.parse(variablesStr))
    for (key,value) in modelDict["modelRuleFunctions"]
        if occursin(Regex("\\b" * key * "\\b"), strObserveble)
            explicitRules *= "\t" * key * " = " * value[2] * "\n"
        end
    end

    if length(explicitRules)>0
        explicitRules *= "\n"
    end

    return stateStr, paramDynStr, paramNonDynStr, explicitRules, namesExplicitRules, statesArray, parameterArray
end


"""
createDYmodFunction(modelName::String, 
                       dirModel::String, 
                       stateNames, 
                       paramData::ParamData, 
                       namesParamDyn::Array{String, 1}, 
                       observablesData::DataFrame,
                       modelDict::Dict)
For modelName create a function for computing DyMod/Du and DyMod/Dp
"""
function createDYmodFunction(modelName::String, 
                            dirModel::String, 
                            stateNames, 
                            paramData::ParamData, 
                            namesParamDyn::Array{String, 1}, 
                            namesNonDynParam::Array{String, 1},
                            observablesData::DataFrame,
                            modelDict::Dict)

    io = open(dirModel * "/" * modelName * "DObsSdU0.jl", "w")
    stateStr, paramDynStr, paramNonDynStr, explicitRules, namesExplicitRules, statesArray, parameterArray = createTopOfDFun(stateNames, paramData, namesParamDyn, namesNonDynParam, observablesData, modelDict)
  
    # Store the formula of each observable in string
    observableIDs = String.(observablesData[!, "observableId"])
    strObservebleU = ""
    strObservebleX = ""
    for i in eachindex(observableIDs)
        # Each observebleID falls below its own if-statement 
        strObservebleU *= "\tif observableId == " * "\"" * observableIDs[i] * "\"" * " \n"
        strObservebleX *= "\tif observableId == " * "\"" * observableIDs[i] * "\"" * " \n"
        tmpFormula = filter(x -> !isspace(x), String(observablesData[i, "observableFormula"]))


        # Replace the explicit rule variable with the explicit rule
        for (key,value) in modelDict["modelRuleFunctions"]            
            tmpFormula = replace(tmpFormula, key => "(" * value[2] * ")")
        end
        # Translate the formula for the observable to Julia syntax 
        juliaFormula = peTabFormulaToJulia(tmpFormula, stateNames, paramData, namesParamDyn, namesNonDynParam, namesExplicitRules)
        
        printInd = 0

        for stateInd in eachindex(statesArray)
            if occursin(Regex("\\b" * statesArray[stateInd] * "\\b"), juliaFormula)
                # Extract observable parameters 
                obsParam = getObsParamStr(tmpFormula)
                if !isempty(obsParam) && printInd == 0
                    strObservebleU *= "\t\t" * obsParam * " = getObsOrSdParam(obsPar, mapObsParam)\n" 
                    printInd = 1
                end 
                
                juliaFormulaSym =  eval(Meta.parse(juliaFormula))
                tmpVar = eval(Meta.parse(statesArray[stateInd]))
                dYdu = Symbolics.derivative(juliaFormulaSym, tmpVar; simplify=true)
                strObservebleU *= "\t\tout[" * string(stateInd) * "] = " * string(dYdu) * "\n"
            end
        end

        printInd = 0
        for paramInd in eachindex(parameterArray)
            if occursin(Regex("\\b" * parameterArray[paramInd] * "\\b"), juliaFormula)
                # Extract observable parameters 
                obsParam = getObsParamStr(tmpFormula)
                if !isempty(obsParam) && printInd == 0
                    strObservebleX *= "\t\t" * obsParam * " = getObsOrSdParam(obsPar, mapObsParam)\n" 
                    printInd = 1
                end 
        
                juliaFormulaSym =  eval(Meta.parse(juliaFormula))
                tmpVar = eval(Meta.parse(parameterArray[paramInd]))
                dYdp = Symbolics.derivative(juliaFormulaSym, tmpVar; simplify=true)
                strObservebleX *= "\t\tout[" * string(paramInd) * "] = " * string(dYdp) * "\n"
            end
        end
        
        strObservebleU *= "\t\t" * "return nothing\n"
        strObservebleU *= "\tend\n\n"
        strObservebleX *= "\t\t" * "return nothing\n"
        strObservebleX *= "\tend\n\n"
    end
    

    write(io, "function evalDYmodDu(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) \n")
    write(io, stateStr)
    write(io, paramDynStr)
    write(io, paramNonDynStr)
    write(io, explicitRules)
    write(io, strObservebleU)
    strClose = "end\n\n"
    write(io, strClose)

    write(io, "function evalDYmodDp(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) \n")
    write(io, stateStr)
    write(io, paramDynStr)
    write(io, paramNonDynStr)
    write(io, explicitRules)
    write(io, strObservebleX)
    strClose = "end\n\n"
    write(io, strClose)

    close(io)
end


"""
    createDSdFunction(modelName::String, 
                          dirModel::String, 
                          paramData::ParamData, 
                          stateNames, 
                          namesParamDyn::Array{String, 1}, 
                          observablesData::DataFrame,
                          modelDict::Dict)
    For modelName create a function for computing the standard deviation by translating the observablesData
    PeTab-file into Julia syntax. 
    To correctly create the function the state-names, names of dynamic parameters to estiamte 
    (namesDynParam) and PeTab parameter-file (to get constant parameters) data are needed. 
"""
function createDSdFunction(modelName::String, 
                          dirModel::String, 
                          paramData::ParamData, 
                          stateNames, 
                          namesParamDyn::Array{String, 1}, 
                          namesNonDynParam::Array{String, 1},
                          observablesData::DataFrame,
                          modelDict::Dict)

    io = open(dirModel * "/" * modelName * "DObsSdU0.jl", "a")
    
    stateStr, paramDynStr, paramNonDynStr, explicitRules, namesExplicitRules, statesArray, parameterArray = createTopOfDFun(stateNames, paramData, namesParamDyn, namesNonDynParam, observablesData, Dict("modelRuleFunctions" => Dict()))
    
    # Store the formula for standard deviations in string
    observableIDs = String.(observablesData[!, "observableId"])
    strObservebleU = ""
    strObservebleX = ""
    for i in eachindex(observableIDs)
        # Each observebleID falls below its own if-statement 
        strObservebleU *= "\tif observableId == " * "\"" * observableIDs[i] * "\"" * " \n"
        strObservebleX *= "\tif observableId == " * "\"" * observableIDs[i] * "\"" * " \n"
        tmpFormula = filter(x -> !isspace(x), String(observablesData[i, "noiseFormula"]))

        # Replace the explicit rule variable with the explicit rule
        for (key,value) in modelDict["modelRuleFunctions"]            
            tmpFormula = replace(tmpFormula, key => "(" * value[2] * ")")
        end
        
        juliaFormula = peTabFormulaToJulia(tmpFormula, stateNames, paramData, namesParamDyn, namesNonDynParam, namesExplicitRules)

        printInd = 0

        for stateInd in eachindex(statesArray)
            if occursin(Regex("\\b" * statesArray[stateInd] * "\\b"), juliaFormula)
                # Extract noise parameters 
                noiseParam = getNoiseParamStr(tmpFormula)
                if !isempty(noiseParam) && printInd == 0
                    strObservebleU *= "\t\t" * noiseParam * " = getObsOrSdParam(sdPar, mapSdParam)\n" 
                    printInd = 1
                end 
        
                juliaFormulaSym =  eval(Meta.parse(juliaFormula))
                tmpVar = eval(Meta.parse(statesArray[stateInd]))
                dYdu = Symbolics.derivative(juliaFormulaSym, tmpVar; simplify=true)
                strObservebleU *= "\t\tout[" * string(stateInd) * "] = " * string(dYdu) * "\n"
            end
        end

        printInd = 0
        for paramInd in eachindex(parameterArray)
            if occursin(Regex("\\b" * parameterArray[paramInd] * "\\b"), juliaFormula)
                juliaFormulaSym =  eval(Meta.parse(juliaFormula))
                # Extract noise parameters 
                noiseParam = getNoiseParamStr(tmpFormula)
                if !isempty(noiseParam) && printInd == 0
                    strObservebleX *= "\t\t" * noiseParam * " = getObsOrSdParam(sdPar, mapSdParam)\n" 
                    printInd = 1
                end 
                
                tmpVar = eval(Meta.parse(parameterArray[paramInd]))
                dYdp = Symbolics.derivative(juliaFormulaSym, tmpVar; simplify=true)
                strObservebleX *= "\t\tout[" * string(paramInd) * "] = " * string(dYdp) * "\n"
            end
        end

        strObservebleU *= "\t\t" * "return nothing\n"
        strObservebleU *= "\tend\n\n"
        strObservebleX *= "\t\t" * "return nothing\n"
        strObservebleX *= "\tend\n\n"
    end

    write(io, "function evalDSdDu!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) \n")
    write(io, stateStr)
    write(io, paramDynStr)
    write(io, paramNonDynStr)
    write(io, explicitRules)
    write(io, strObservebleU)
    strClose = "end\n\n"
    write(io, strClose)

    write(io, "function evalDSdDp!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) \n")
    write(io, stateStr)
    write(io, paramDynStr)
    write(io, paramNonDynStr)
    write(io, explicitRules)
    write(io, strObservebleX)
    strClose = "end\n\n"
    write(io, strClose)
    
    close(io)
end