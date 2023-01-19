
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
    createTopOfDFun(stateNames, 
        namesParamODEProb::Array{String, 1}, 
        namesNonDynParam::Array{String, 1},
        observablesData::DataFrame,
        modelDict::Dict)
    Extracts all variables needed for the functions. 
    Also adds them as variables for Symbolics.jl
"""
function createTopOfDFun(stateNames, 
                        namesParamODEProb::Array{String, 1}, 
                        namesNonDynParam::Array{String, 1},
                        observablesData::DataFrame)

    variablesStr = "@variables "

    # Extract names of model states 
    stateNamesShort = replace.(string.(stateNames), "(t)" => "")
    stateStr = "#"
    for i in eachindex(stateNamesShort)
        stateStr *= "u[" * string(i) * "] = " * stateNamesShort[i] * ", "
        variablesStr *= stateNamesShort[i] * ", "
    end
    stateStr = stateStr[1:end-2]
    stateStr *= "\n"
    
    # Extract name of dynamic parameter 
    paramDynStr = "#"
    for i in eachindex(namesParamODEProb)
        paramDynStr *= "dynPar[" * string(i) * "] = " * namesParamODEProb[i] * ", "
        variablesStr *= namesParamODEProb[i] * ", "
    end
    paramDynStr = paramDynStr[1:end-2]
    paramDynStr *= "\n"

    paramNonDynStr = ""
    # Extract name of non-dynamic parameter
    if !isempty(namesNonDynParam)
        paramNonDynStr = "#"
        for i in eachindex(namesNonDynParam)
            paramNonDynStr *= "nonDynParam[" * string(i)* "] = " * namesNonDynParam[i] * ", "
            variablesStr *= namesNonDynParam[i] * ", "
        end
        paramNonDynStr = paramNonDynStr[1:end-2]
        paramNonDynStr *= "\n"
    end

    # Extracts all observable- and noise-parameters    
    observableIDs = String.(observablesData[!, "observableId"])
    for i in eachindex(observableIDs)

        # Extract observable parameters 
        tmpFormula = filter(x -> !isspace(x), String(observablesData[i, "observableFormula"]))
        obsParam = getObsParamStr(tmpFormula)
        if !isempty(obsParam)
            variablesStr *= obsParam * ", "   
        end 

        # Extract noise parameters 
        tmpFormula = filter(x -> !isspace(x), String(observablesData[i, "noiseFormula"]))
        noiseParam = getNoiseParamStr(tmpFormula)
        if !isempty(noiseParam)
            variablesStr *= noiseParam * ", "   
        end 

    end

    # Remove last "," and Run @variables ... string
    variablesStr = variablesStr[1:end-2]
    eval(Meta.parse(variablesStr))

    return stateStr, paramDynStr, paramNonDynStr
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
    stateStr, paramDynStr, paramNonDynStr = createTopOfDFun(stateNames, namesParamDyn, namesNonDynParam, observablesData)
    statesArray = replace.(string.(stateNames), "(t)" => "")

    # Store the formula of each observable in string
    observableIDs = String.(observablesData[!, "observableId"])
    strObservebleU = ""
    strObservebleX = ""
    for i in eachindex(observableIDs)
        # Each observebleID falls below its own if-statement 
        strObservebleU *= "\tif observableId == " * "\"" * observableIDs[i] * "\"" * " \n"
        strObservebleX *= "\tif observableId == " * "\"" * observableIDs[i] * "\"" * " \n"
        tmpFormula = filter(x -> !isspace(x), String(observablesData[i, "observableFormula"]))
        tmpFormula = replaceExplicitVariableWithRule(tmpFormula, modelDict)

        # Translate the formula for the observable to Julia syntax 
        juliaFormula = peTabFormulaToJulia(tmpFormula, stateNames, paramData, namesParamDyn, namesNonDynParam)
        
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
                dYdu = string(dYdu)
                dYdu = replaceVariablesWithArrayIndex(dYdu, stateNames, namesParamDyn, namesNonDynParam, paramData)

                strObservebleU *= "\t\tout[" * string(stateInd) * "] = " * dYdu * "\n"
            end
        end

        printInd = 0
        for paramInd in eachindex(namesParamDyn)
            if occursin(Regex("\\b" * namesParamDyn[paramInd] * "\\b"), juliaFormula)
                # Extract observable parameters 
                obsParam = getObsParamStr(tmpFormula)
                if !isempty(obsParam) && printInd == 0
                    strObservebleX *= "\t\t" * obsParam * " = getObsOrSdParam(obsPar, mapObsParam)\n" 
                    printInd = 1
                end 
        
                juliaFormulaSym =  eval(Meta.parse(juliaFormula))
                tmpVar = eval(Meta.parse(namesParamDyn[paramInd]))
                dYdp = Symbolics.derivative(juliaFormulaSym, tmpVar; simplify=true)
                dYdp = string(dYdp)
                dYdp = replaceVariablesWithArrayIndex(dYdp, stateNames, namesParamDyn, namesNonDynParam, paramData)
                strObservebleX *= "\t\tout[" * string(paramInd) * "] = " * dYdp * "\n"
            end
        end
        
        strObservebleU *= "\t\t" * "return nothing\n"
        strObservebleU *= "\tend\n\n"
        strObservebleX *= "\t\t" * "return nothing\n"
        strObservebleX *= "\tend\n\n"
    end
    

    write(io, stateStr)
    write(io, paramDynStr)
    write(io, paramNonDynStr)
    write(io, "\n")
    write(io, "function evalDYmodDu(u, t, dynPar, obsPar, nonDynParam, observableId, mapObsParam, out) \n")
    write(io, strObservebleU)
    strClose = "end\n\n"
    write(io, strClose)

    write(io, "function evalDYmodDp(u, t, dynPar, obsPar, nonDynParam, observableId, mapObsParam, out) \n")
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
    
    statesArray = replace.(string.(stateNames), "(t)" => "")


    # Store the formula for standard deviations in string
    observableIDs = String.(observablesData[!, "observableId"])
    strObservebleU = ""
    strObservebleX = ""
    for i in eachindex(observableIDs)
        # Each observebleID falls below its own if-statement 
        strObservebleU *= "\tif observableId == " * "\"" * observableIDs[i] * "\"" * " \n"
        strObservebleX *= "\tif observableId == " * "\"" * observableIDs[i] * "\"" * " \n"
        tmpFormula = filter(x -> !isspace(x), String(observablesData[i, "noiseFormula"]))
        tmpFormula = replaceExplicitVariableWithRule(tmpFormula, modelDict)

        # Translate the formula for the noise to Julia syntax 
        juliaFormula = peTabFormulaToJulia(tmpFormula, stateNames, paramData, namesParamDyn, namesNonDynParam)

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
                dYdu = string(dYdu)
                dYdu = replaceVariablesWithArrayIndex(dYdu, stateNames, namesParamDyn, namesNonDynParam, paramData)
                strObservebleU *= "\t\tout[" * string(stateInd) * "] = " * dYdu * "\n"
            end
        end

        printInd = 0
        for paramInd in eachindex(namesParamDyn)
            if occursin(Regex("\\b" * namesParamDyn[paramInd] * "\\b"), juliaFormula)
                juliaFormulaSym =  eval(Meta.parse(juliaFormula))
                # Extract noise parameters 
                noiseParam = getNoiseParamStr(tmpFormula)
                if !isempty(noiseParam) && printInd == 0
                    strObservebleX *= "\t\t" * noiseParam * " = getObsOrSdParam(sdPar, mapSdParam)\n" 
                    printInd = 1
                end 
                
                tmpVar = eval(Meta.parse(namesParamDyn[paramInd]))
                dYdp = Symbolics.derivative(juliaFormulaSym, tmpVar; simplify=true)
                dYdp = string(dYdp)
                dYdp = replaceVariablesWithArrayIndex(dYdp, stateNames, namesParamDyn, namesNonDynParam, paramData)
                strObservebleX *= "\t\tout[" * string(paramInd) * "] = " * dYdp * "\n"
            end
        end

        strObservebleU *= "\t\t" * "return nothing\n"
        strObservebleU *= "\tend\n\n"
        strObservebleX *= "\t\t" * "return nothing\n"
        strObservebleX *= "\tend\n\n"
    end

    write(io, "function evalDSdDu!(u, t, sdPar, dynPar, nonDynParam, paramData, observableId, mapSdParam, out) \n")
    write(io, strObservebleU)
    strClose = "end\n\n"
    write(io, strClose)

    write(io, "function evalDSdDp!(u, t, sdPar, dynPar, nonDynParam, paramData, observableId, mapSdParam, out) \n")
    write(io, strObservebleX)
    strClose = "end\n\n"
    write(io, strClose)
    
    close(io)
end