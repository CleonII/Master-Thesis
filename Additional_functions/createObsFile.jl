using ModelingToolkit, DifferentialEquations, BenchmarkTools, DataFrames, CSV, LSODA, Sundials, ODEInterface, ODEInterfaceDiffEq

include(joinpath(pwd(), "Pipeline_ModelSolver", "BigFloatODEProblem.jl"))
include(joinpath(pwd(), "Additional_functions", "additional_tools.jl"))
include(joinpath(pwd(), "Additional_functions", "Solver_info.jl"))
include(joinpath(pwd(), "Additional_functions", "benchmarkSolvers.jl"))

# Ipopt wrapper 
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "CommonParameterEstimationMethods", "CreateIpoptProb.jl"))

# PeTab importer 
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "PeTabImporter.jl"))

allModelFunctionVector = includeAllModels(getModelFiles(joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels")), 
        joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels"))


# Create a file with a function for the initial values and the observation function for a model. 
# Must extend to observeble parameters appearing only in the observation function, and that aren't 
# scale and offset parameters. 
function crateFilesModel(modelName::String)
    
    modelFile = modelName * ".jl"
    modelNameShort = modelName[7:end]
    readPath = joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels")
    readDataPath = joinpath(pwd(), "Pipeline_ModelParameterEstimation", "Data", modelName)

    # Extract relevant model file (ode-system)
    allModelFiles = getModelFiles(readPath)    
    usedModelFunctionVector = allModelFunctionVector[[allModelFile .== modelFile for allModelFile in allModelFiles]][1]
    
    # Create an ODE-system and run via ModellingToolkit
    sys, initialSpeciesValues, trueParameterValues = usedModelFunctionVector()

    # Read data on experimental conditions and parameter values 
    experimentalConditions, measurementData, parameterBounds, observableData = readDataFiles(modelName, readObs=true)
    paramData = processParameterData(parameterBounds) # Model data in convient structure 
    firstExpIds, shiftExpIds, simulateSS, parameterNames, stateNames = getSimulationInfo(measurementData, sys)

    # Read file with observed values and put into struct 
    obsData = processObsData(measurementData, observableData)
    
    # Set up a stateMap and paramMap to map parameter correctly to the ODE model via ModellingToolkit
    paramMap = trueParameterValues    
    stateMap = initialSpeciesValues 
    paramNames = [string(paramMap[i].first) for i in eachindex(paramMap)]

    # Set up to bookeeping of indices between observed and noise parameters in the large input vector 
    iDynPar, iSdPar, iObsPar, idParamDyn, idSdParam, idObsParam = getIndicesParam(paramData, obsData)
    
    createObservebleFile(modelNameShort, readDataPath, stateNames, paramData, idParamDyn)

    println("Done with observation file")
    println("")
    
    createU0File(modelNameShort, readDataPath, paramData, string.(parameterNames), stateMap)

    println("Done with u0 observation file")
    println("")

    createSdFile(modelNameShort, readDataPath, paramData, stateNames, idParamDyn)

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
function createObservebleFile(modelNameShort, readDataPath, stateNames, paramData::ParamData, idParamDyn::Array{String, 1})

    observableDatas = CSV.read(joinpath(readDataPath, "observables_" * modelNameShort * ".tsv"), DataFrame)

    io = open(readDataPath * "/" * modelNameShort * "Obs.jl", "w")
    
    # Write funciton header 
    write(io, "function " * modelNameShort * "(u, t, dynPar, obsPar, paramData, obsData, observableId, simulationId) \n")

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
    observableIDs = String.(observableDatas[!, "observableId"])
    strObserveble = ""
    for i in eachindex(observableIDs)
        # Each observebleID falls below its own if-statement 
        strObserveble *= "\tif observableId == " * "\"" * observableIDs[i] * "\"" * " \n"
        tmpFormula = filter(x -> !isspace(x), String(observableDatas[i, "observableFormula"]))

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


function createU0File(modelNameShort, readDataPath, paramData::ParamData, idParam, stateMap)

    io = open(readDataPath * "/" * modelNameShort * "Obs.jl", "a")

    # Write funciton header 
    write(io, "\n\nfunction " * modelNameShort * "_t0!(u0Vec, paramVec) \n\n")

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
function createSdFile(modelNameShort, readDataPath, paramData::ParamData, stateNames, idParamDyn)

    observableDatas = CSV.read(joinpath(readDataPath, "observables_" * modelNameShort * ".tsv"), DataFrame)

    io = open(readDataPath * "/" * modelNameShort * "Obs.jl", "a")

    # Write funciton header 
    write(io, "\n\nfunction " * modelNameShort * "_sd!(u, t, sdPar, dynPar, paramData, obsData, observableId, simulationId) \n")

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
    observableIDs = String.(observableDatas[!, "observableId"])
    strObserveble = ""
    for i in eachindex(observableIDs)
        # Each observebleID falls below its own if-statement 
        strObserveble *= "\tif observableId == " * "\"" * observableIDs[i] * "\"" * " \n"
        tmpFormula = filter(x -> !isspace(x), String(observableDatas[i, "noiseFormula"]))

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


# "model_Bachmann_MSB2011"
# "model_Boehm_JProteomeRes2014"
crateFilesModel("model_Bachmann_MSB2011")



