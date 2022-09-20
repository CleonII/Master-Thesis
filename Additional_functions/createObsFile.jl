using ModelingToolkit, DifferentialEquations, BenchmarkTools, DataFrames, CSV, LSODA, Sundials, ODEInterface, ODEInterfaceDiffEq

include(joinpath(pwd(), "Additional_functions", "additional_tools.jl"))
include(joinpath(pwd(), "Additional_functions", "Solver_info.jl"))
include(joinpath(pwd(), "Additional_functions", "benchmarkSolvers.jl"))

# Ipopt wrapper 
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "CommonParameterEstimationMethods", "CreateIpoptProb.jl"))

# PeTab importer 
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "PeTabImporter.jl"))

# Create a file with a function for the initial values and the observation function for a model. 
# Must extend to observeble parameters appearing only in the observation function, and that aren't 
# scale and offset parameters. 
function createFilesModel(modelName::String, 
                          dirModel::String, 
                          odeSys::ODESystem, 
                          paramMap, 
                          stateMap)
    
    # For correct mapping of parameter and states when building related functions 
    parameterNames = parameters(odeSys)
    stateNames = states(odeSys)

    # Read data on experimental conditions and parameter values 
    experimentalConditionsFile, measurementDataFile, parameterBoundsFile, observablesDataFile = readDataFiles(dirModel, readObs=true)
    paramData = processParameterData(parameterBoundsFile) # Model data in convient structure 
    simulationInfo = getSimulationInfo(measurementDataFile, peTabModel.odeSystem)
    measurementData = processObsData(measurementDataFile, observablesDataFile) 
    
    # Set up a stateMap and paramMap to map parameter correctly to the ODE model via ModellingToolkit
    paramNames = [string(paramMap[i].first) for i in eachindex(paramMap)]

    # Indices for mapping parameter-estimation vector to dynamic, observable and sd parameters correctly when calculating cost
    paramIndices = getIndicesParam(paramData, measurementData)
    
    createObservebleFile(modelName, dirModel, stateNames, paramData, paramIndices.namesDynParam, observablesDataFile)
    println("Done with observation file")
    println("")
    
    createU0File(modelName, dirModel, paramData, string.(parameterNames), stateMap)
    println("Done with u0 observation file")
    println("")

    createSdFile(modelName, dirModel, paramData, stateNames, paramIndices.namesDynParam, observablesDataFile)
    println("Done with sd observation file")
    println("")
end


# Helper function to identify all words with observable in a formula string 
function getObsParamStr(measurmentFormula)
    
    # Find all words on the form observableParameter
    obsWords = [ match.match for match in eachmatch(r"observableParameter[0-9]_\w+", measurmentFormula) ]
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


function getNoiseParamStr(measurmentFormula)
    
    # Find all words on the form observableParameter
    sdWords = [ match.match for match in eachmatch(r"noiseParameter[0-9]_\w+", measurmentFormula) ]
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


# Create a Julia file with the observeble format provided in the PeTab files. 
function createObservebleFile(modelName::String, 
                              dirModel::String, 
                              stateNames, 
                              paramData::ParamData, 
                              idParamDyn::Array{String, 1}, 
                              observablesData::DataFrame)

    io = open(dirModel * "/" * modelName * "ObsSdU0.jl", "w")
    
    # Write funciton header 
    write(io, "function evalYmod(u, t, dynPar, obsPar, paramData, obsData, observableId, simulationId) \n")

    # Construct string to extract states
    stateNamesShort = replace.(string.(stateNames), "(t)" => "")
    stateStr = "\n\t"
    for i in eachindex(stateNamesShort)
        stateStr *= stateNamesShort[i] * ", "
    end
    stateStr = stateStr[1:end-2]
    stateStr *= "= u \n"
    write(io, stateStr)

    # Construct string extracting dynPar (dynamic parameter to be estimated)
    paramDynStr = "\t"
    for i in eachindex(idParamDyn)
        paramDynStr *= idParamDyn[i] * ", "
    end
    paramDynStr = paramDynStr[1:end-2]
    paramDynStr *= " = dynPar \n"
    write(io, paramDynStr)

    # String to extract constant parameters. To prevent cluttering only the constant parameters 
    # are extracted in the file. 
    paramConstStr = ""
    for i in eachindex(paramData.parameterID)
        if paramData.shouldEst[i] == false
            paramConstStr *= "\t" * paramData.parameterID[i] * "_C = paramData.paramVal[" * string(i) *"] \n" 
        end
    end
    paramConstStr *= "\n"
    write(io, paramConstStr)

    # Fix the observeble formulas 
    observableIDs = String.(observablesData[!, "observableId"])
    strObserveble = ""
    for i in eachindex(observableIDs)
        # Each observebleID falls below its own if-statement 
        strObserveble *= "\tif observableId == " * "\"" * observableIDs[i] * "\"" * " \n"
        tmpFormula = filter(x -> !isspace(x), String(observablesData[i, "observableFormula"]))

        # Extract observable parameters 
        obsParam = getObsParamStr(tmpFormula)
        if !isempty(obsParam)
            strObserveble *= "\t\t" * obsParam * " = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)\n" 
        end 

        juliaFormula = obsFormulaToJulia(tmpFormula, stateNames, paramData, idParamDyn)
        strObserveble *= "\t\t" * "return " * juliaFormula * "\n"
        strObserveble *= "\tend\n\n"
    end
    write(io, strObserveble)

    # Close file 
    strClose = "end"
    write(io, strClose)
    close(io)
end


function createU0File(modelName::String, 
                      dirModel::String, 
                      paramData::ParamData, 
                      idParam, 
                      stateMap)

    io = open(dirModel * "/" * modelName * "ObsSdU0.jl", "a")

    # Write funciton header 
    write(io, "\n\nfunction evalU0!(u0Vec, paramVec) \n\n")

    # Extract model parameters 
    paramDynStr = "\t"
    for i in eachindex(idParam)
        paramDynStr *= idParam[i] * ", "
    end
    paramDynStr = paramDynStr[1:end-2]
    paramDynStr *= " = paramVec \n\n"
    write(io, paramDynStr)    

    # Write initial values 
    stateNames = [replace.(string.(stateMap[i].first), "(t)" => "") for i in eachindex(stateMap)]
    stateExpWrite = ""
    for i in eachindex(stateMap)
        stateName = stateNames[i]
        stateExp = replace(string(stateMap[i].second), " " => "")
        stateFormula = obsFormulaToJulia(stateExp, stateNames, paramData, idParam)
        stateExpWrite *= "\t" * stateName * " = " * stateFormula * "\n"
    end
    write(io, stateExpWrite * "\n")

    # Map initial values to u0 vector 
    stateStr = "\tu0Vec .= "
    for i in eachindex(stateNames)
        stateStr *= stateNames[i] * ", "
    end
    stateStr = stateStr[1:end-2]
    write(io, stateStr)

    # Close file 
    strClose = "\nend"
    write(io, strClose)
    close(io)
end


# Create a function for any noise parameters 
function createSdFile(modelName::String, 
                      dirModel::String, 
                      paramData::ParamData, 
                      stateNames, 
                      idParamDyn, 
                      observablesData::DataFrame)


    io = open(dirModel * "/" * modelName * "ObsSdU0.jl", "a")

    # Write funciton header 
    write(io, "\n\nfunction evalSd!(u, t, sdPar, dynPar, paramData, obsData, observableId, simulationId) \n")

    # Construct string to extract states
    stateNamesShort = replace.(string.(stateNames), "(t)" => "")
    stateStr = "\n\t"
    for i in eachindex(stateNamesShort)
        stateStr *= stateNamesShort[i] * ", "
    end
    stateStr = stateStr[1:end-2]
    stateStr *= "= u \n"
    write(io, stateStr)

    # Construct string extracting dynPar (dynamic parameter to be estimated)
    paramDynStr = "\t"
    for i in eachindex(idParamDyn)
        paramDynStr *= idParamDyn[i] * ", "
    end
    paramDynStr = paramDynStr[1:end-2]
    paramDynStr *= " = dynPar \n"
    write(io, paramDynStr)

    # String to extract constant parameters. To prevent cluttering only the constant parameters 
    # are extracted in the file. 
    paramConstStr = ""
    for i in eachindex(paramData.parameterID)
        if paramData.shouldEst[i] == false
            paramConstStr *= "\t" * paramData.parameterID[i] * "_C = paramData.paramVal[" * string(i) *"] \n" 
        end
    end
    paramConstStr *= "\n"
    write(io, paramConstStr)

    # Fix the observeble formulas 
    observableIDs = String.(observablesData[!, "observableId"])
    strObserveble = ""
    for i in eachindex(observableIDs)
        # Each observebleID falls below its own if-statement 
        strObserveble *= "\tif observableId == " * "\"" * observableIDs[i] * "\"" * " \n"
        tmpFormula = filter(x -> !isspace(x), String(observablesData[i, "noiseFormula"]))

        # Extract noise parameters 
        noiseParam = getNoiseParamStr(tmpFormula)
        if !isempty(noiseParam)
            strObserveble *= "\t\t" * noiseParam * " = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)\n" 
        end 

        juliaFormula = obsFormulaToJulia(tmpFormula, stateNames, paramData, idParamDyn)
        strObserveble *= "\t\t" * "return " * juliaFormula * "\n"
        strObserveble *= "\tend\n\n"
    end
    write(io, strObserveble)

    # Close file 
    strClose = "end"
    write(io, strClose)
    close(io)
end


function translateTermObsFunc(termTranslate::String, 
                              stateNames,
                              paramData::ParamData, 
                              idParamDyn::Array{String, 1})

    listOperations = ["exp", "sin", "cos", "t"]
    stateNamesShort = replace.(string.(stateNames), "(t)" => "")
    
    strAdd = ""
    if termTranslate in paramData.parameterID && !(termTranslate in idParamDyn)
        # Constant parameters get a _C appended to tell them apart 
        strAdd *= termTranslate * "_C"
    end
    if termTranslate in idParamDyn
        strAdd *= termTranslate
    end
    if termTranslate in stateNamesShort
        strAdd *= termTranslate
    end
    if isNumber(termTranslate)
        strAdd *= termTranslate
    end
    if termTranslate in listOperations
        strAdd *= listOperations[termTranslate .== listOperations]
    end
    if length(termTranslate) >= 19 && termTranslate[1:19] == "observableParameter"
        strAdd *= termTranslate
    end
    if length(termTranslate) >= 14 && termTranslate[1:14] == "noiseParameter"
        strAdd *= termTranslate
    end

    if isempty(strAdd)
        println(idParamDyn)
        println("Warning : When creating observation function $termTranslate could not be processed")
    end

    strAdd *= " "

    return strAdd
end



# In a string extract the next "word", which is the longest concurent occurance 
# of characters that re not in the character list with word termination characters. 
# Also, in the case 102var 102 is considered a word. 
function getWord(str, iStart, charListTerm)
    
    wordStr = ""
    i = iStart

    # If the first character is a numberic the termination occurs when 
    # the first non-numeric character (or not dot) is reached. 
    isNumericStart = isnumeric(str[i])

    while i <= length(str)
        if !(str[i] in charListTerm) 
            if isNumericStart == true && !(isnumeric(str[i]) || str[i] == '.')
                break
            end
            wordStr *= str[i]
        else
            break 
        end
        i += 1
    end

    return wordStr, i 
end


# Translate a string on the PeTAB formula format into Julia code. 
function obsFormulaToJulia(formulaObs::String, stateNames, paramData::ParamData, idParamDyn::Array{String, 1})

    charListTerm = ['(', ')', '+', '-', '/', '*', '^']
    lenFormula = length(formulaObs)
    strFormula = ""
    i = 1
    while i <= lenFormula
        if formulaObs[i] in charListTerm
            strFormula *= formulaObs[i] * " "
            i += 1
        else
            termTmp, iNew = getWord(formulaObs, i, charListTerm)
            strFormula *= translateTermObsFunc(termTmp, stateNames, paramData, idParamDyn)
            i = iNew

            # Special case where we have multiplication
            if isNumber(termTmp) && i <= lenFormula && isletter(formulaObs[i])
                strFormula *= "* "
            end
        end
    end

    return strFormula
end





