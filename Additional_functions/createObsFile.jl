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
    experimentalConditions, measurementData, parameterBounds = readDataFiles(modelName)
    paramData = processParameterData(parameterBounds) # Model data in convient structure 
    firstExpIds, shiftExpIds, simulateSS, parameterNames, stateNames = getSimulationInfo(measurementData, sys)

    # Read file with observed values and put into struct 
    obsData = processObsData(measurementData)
    
    # Set up a stateMap and paramMap to map parameter correctly to the ODE model via ModellingToolkit
    paramMap = trueParameterValues    
    stateMap = initialSpeciesValues 
    paramNames = [string(paramMap[i].first) for i in eachindex(paramMap)]

    # Set up to bookeeping of indices between observed and noise parameters in the large input vector 
    isSd = [paramData.parameterID[i] in obsData.sdParamId for i in eachindex(paramData.parameterID)]
    isDynamic = (paramData.shouldEst .&& .!isSd)
    idParamDyn = paramData.parameterID[isDynamic]
    idParamSd = paramData.parameterID[isSd]
    
    createObservebleFile(modelNameShort, readDataPath, stateNames, paramData, idParamDyn)
    
    createU0File(modelNameShort, readDataPath, paramData, string.(parameterNames), stateMap)
end


# Create a Julia file with the observeble format provided in the PeTab files. 
function createObservebleFile(modelNameShort, readDataPath, stateNames, paramData::ParamData, idParamDyn::Array{String, 1})

    observableDatas = CSV.read(joinpath(readDataPath, "observables_" * modelNameShort * ".tsv"), DataFrame)

    io = open(readDataPath * "/" * modelNameShort * "Obs.jl", "w")
    
    # Write funciton header 
    write(io, "function " * modelNameShort * "(u, t, dynPar, obsPar, paramData, observableId) \n")

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
        juliaFormula = obsFormulaToJulia(tmpFormula, stateNames, paramData, idParamDyn)
        strObserveble *= "\t\t" * "yMod = " * juliaFormula * "\n"
        strObserveble *= "\tend\n\n"
    end
    write(io, strObserveble)

    # Close file 
    strClose = "\treturn yMod\nend"
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
        stateFormula = obsFormulaToJulia(stateExp, stateNames, paramData, paramData.parameterID)
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

    if isempty(strAdd)
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

    println("str = $str")

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

crateFilesModel("model_Boehm_JProteomeRes2014")




#=
function Boehm_JProteomeRes2014(u, t, dynPar, obsPar, paramData, observableId) 

	STAT5A, pApA, nucpApB, nucpBpB, STAT5B, pApB, nucpApA, pBpB, dummyVariable= u 
	Epo_degradation_BaF3, k_exp_hetero, k_exp_homo, k_imp_hetero, k_imp_homo, k_phos= dynPar 
	ratio_C = paramData.paramVal[7] 
	specC17_C = paramData.paramVal[11] 

	if observableId == "pSTAT5A_rel" 
		yMod = ( 100 * pApB + 200 * pApA * specC17_C ) / ( pApB + STAT5A * specC17_C + 2 * pApA * specC17_C ) 
	end

	if observableId == "pSTAT5B_rel" 
		yMod = - ( 100 * pApB - 200 * pBpB * ( specC17_C - 1 ) ) / ( ( STAT5B * ( specC17_C - 1 ) - pApB ) + 2 * pBpB * ( specC17_C - 1 ) ) 
	end

	if observableId == "rSTAT5A_rel" 
		yMod = ( 100 * pApB + 100 * STAT5A * specC17_C + 200 * pApA * specC17_C ) / ( 2 * pApB + STAT5A * specC17_C + 2 * pApA * specC17_C - STAT5B * ( specC17_C - 1 ) - 2 * pBpB * ( specC17_C - 1 ) ) 
	end

	return yMod
end


function Boehm_JProteomeRes2014_t0!(u0Vec, paramVec) 

	k_exp_hetero, Epo_degradation_BaF3, k_exp_homo, k_phos, k_imp_homo, k_imp_hetero, specC17, ratio = paramVec

	STAT5A = 207.6 * ratio
    pApA = 0.0
    nucpApB = 0.0
    nucpBpB = 0.0
    STAT5B = 207.6 - 207.6 * ratio
    pApB = 0.0
    nucpApA = 0.0
    pBpB = 0.0
    dummyVariable = 0.0
		
	u0Vec .= STAT5A, pApA, nucpApB, nucpBpB, STAT5B, pApB, nucpApA, pBpB, dummyVariable
end

=#