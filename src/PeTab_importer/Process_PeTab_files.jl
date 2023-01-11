"""
    setUpPeTabModel(modelName::String, dirModel::String)::PeTabModel

    Given a model directory (dirModel) containing the PeTab files and a 
    xml-file on format modelName.xml will return a PeTabModel struct holding 
    paths to PeTab files, ode-system in ModellingToolkit format, functions for 
    evaluating yMod, u0 and standard deviations, and a parameter and state maps 
    for how parameters and states are mapped in the ModellingToolkit ODE system
    along with state and parameter names.

    dirModel must contain a SBML file named modelName.xml, and files starting with 
    measurementData, experimentalCondition, parameter, and observables (tsv-files).
    The latter files must be unique (e.g only one file starting with measurementData)

    TODO : Example  
"""
function setUpPeTabModel(modelName::String, dirModel::String; forceBuildJlFile::Bool=false, verbose::Bool=true, ifElseToEvent=true)::PeTabModel

    # Sanity check user input 
    modelFileXml = dirModel * modelName * ".xml"
    modelFileJl = dirModel * modelName * ".jl"
    if !isdir(dirModel)
        if verbose
            @printf("Model directory %s does not exist\n", dirModel)
        end
    end
    if !isfile(modelFileXml)
        if verbose
            @printf("Model directory does not contain xml-file with name %s\n", modelName * "xml")
        end
    end
    # If Julia model file does exists build it 
    if !isfile(modelFileJl) && forceBuildJlFile == false
        if verbose
            @printf("Julia model file does not exist - will build it\n")
        end
        modelDict = XmlToModellingToolkit(modelFileXml, modelName, dirModel, ifElseToEvent=ifElseToEvent)
    elseif isfile(modelFileJl) && forceBuildJlFile == false
        if verbose
            @printf("Julia model file exists at %s - will not rebuild it\n", modelFileJl)
        end
    elseif forceBuildJlFile == true
        if verbose
            @printf("By user option rebuilds Julia model file\n")
        end
        if isfile(modelFileJl)
            rm(modelFileJl)
        end
        modelDict = XmlToModellingToolkit(modelFileXml, modelName, dirModel, ifElseToEvent=ifElseToEvent)
    end

    # Extract ODE-system and mapping of maps of how to map parameters to states and model parmaeters 
    include(modelFileJl)
    expr = Expr(:call, Symbol("getODEModel_" * modelName))
    odeSys, stateMap, paramMap = eval(expr)
    #odeSysUse = ode_order_lowering(odeSys)
    odeSysUse = structural_simplify(odeSys)
    # Parameter and state names for ODE-system 
    parameterNames = parameters(odeSysUse)
    stateNames = states(odeSysUse)

    # Sanity check for presence of all PeTab-files 
    pathMeasurementData = checkForPeTabFile("measurementData", dirModel)
    pathExperimentalCond = checkForPeTabFile("experimentalCondition", dirModel)
    pathParameters = checkForPeTabFile("parameters", dirModel)
    pathObservables = checkForPeTabFile("observables", dirModel)

    # Build functions for observables, sd and u0 if does not exist and include
    pathObsSdU0 = dirModel * modelName * "ObsSdU0.jl"
    pathDObsSdU0 = dirModel * modelName * "DObsSdU0.jl"
    if !isfile(pathObsSdU0) || forceBuildJlFile == true
        if verbose && forceBuildJlFile == false
            @printf("File for yMod, U0 and Sd does not exist - building it\n")
        end
        if verbose && forceBuildJlFile == true
            @printf("By user option will rebuild Ymod, Sd and u0\n")
        end
        if !@isdefined(modelDict)
            modelDict = XmlToModellingToolkit(modelFileXml, modelName, dirModel, writeToFile=false, ifElseToEvent=ifElseToEvent)
        end
        createFileYmodSdU0(modelName, dirModel, odeSysUse, stateMap, modelDict)
        createFileDYmodSdU0(modelName, dirModel, odeSysUse, stateMap, modelDict)
    else
        if verbose
            @printf("File for yMod, U0 and Sd does exist - will not rebuild it\n")
        end
    end
    include(pathObsSdU0)
    include(pathDObsSdU0)    

    pathCallback = dirModel * "/" * modelName * "Callbacks_time_piecewise.jl"
    if !isfile(pathCallback) || forceBuildJlFile == true
        if verbose && forceBuildJlFile == false
            @printf("File for callback does not exist - building it\n")
        end
        if verbose && forceBuildJlFile == true
            @printf("By user option will rebuild callback file\n")
        end
        if !@isdefined(modelDict)
            modelDict = XmlToModellingToolkit(modelFileXml, modelName, dirModel, writeToFile=false, ifElseToEvent=ifElseToEvent)
        end
        getCallbacksForTimeDepedentPiecewise(odeSysUse, modelDict, modelName, dirModel)
    end
    include(pathCallback)
    exprCallback = Expr(:call, Symbol("getCallbacks_" * modelName))
    cbSet::CallbackSet, checkCbActive::Vector{Function} = eval(exprCallback)    

    peTabModel = PeTabModel(modelName,
                            evalYmod,
                            evalU0!,
                            evalU0,
                            evalSd!,
                            evalDYmodDu,
                            evalDSdDu!,
                            evalDYmodDp,
                            evalDSdDp!,
                            getTstops,
                            odeSysUse,
                            paramMap,
                            stateMap,
                            parameterNames, 
                            stateNames,
                            dirModel,
                            pathMeasurementData,
                            pathExperimentalCond,
                            pathObservables, 
                            pathParameters, 
                            cbSet, 
                            checkCbActive)

    return peTabModel
end


"""
    readDataFiles(dirModel::String; readObs::Bool=false)

    Given a directory for a model, e.g ./Beer_MolBioSystems2014, read the associated PeTab files 
    for the measurements, parameters, experimental conditions and (if true) the observables.
"""
function readDataFiles(dirModel::String; readObs::Bool=false)

    # Check if PeTab files exist and get their path 
    pathMeasurementData = checkForPeTabFile("measurementData", dirModel)
    pathExperimentalCond = checkForPeTabFile("experimentalCondition", dirModel)
    pathParameters = checkForPeTabFile("parameters", dirModel)
    pathObservables = checkForPeTabFile("observables", dirModel)

    experimentalConditions = CSV.read(pathExperimentalCond, DataFrame)
    measurementData = CSV.read(pathMeasurementData, DataFrame)
    parameterData = CSV.read(pathParameters, DataFrame)
    if readObs == true
        observableData = CSV.read(pathObservables, DataFrame)
        return experimentalConditions, measurementData, parameterData, observableData
    else
        return experimentalConditions, measurementData, parameterData
    end
end


"""
    processParameterData(parameterData::DataFrame)::ParamData

    Process the PeTab parameterData file into a type-stable Julia struct.
"""
function processParameterData(parameterData::DataFrame)::ParamData

    nParam = length(parameterData[!, "estimate"])

    # Pre-allocate arrays to hold data 
    lb::Array{Float64, 1} = zeros(Float64, nParam) 
    ub::Array{Float64, 1} = zeros(Float64, nParam)
    paramVal::Array{Float64, 1} = zeros(Float64, nParam) # Vector with Nominal value in PeTab-file
    logScale::Array{Bool, 1} = Array{Bool, 1}(undef, nParam)
    paramId::Array{String, 1} = Array{String, 1}(undef, nParam)
    shouldEst::Array{Bool, 1} = Array{Bool, 1}(undef, nParam)

    for i in eachindex(shouldEst)

        # If upper or lower bounds are missing assume +Inf and -Inf respectively.
        if ismissing(parameterData[i, "lowerBound"])
            lb[i] = -Inf
        else
            lb[i] = parameterData[i, "lowerBound"]
        end
        if ismissing(parameterData[i, "upperBound"])
            ub[i] = Inf
        else
            ub[i] = parameterData[i, "upperBound"]
        end

        paramVal[i] = parameterData[i, "nominalValue"]
        paramId[i] = parameterData[i, "parameterId"]
        # Currently only supports parameters on log10-scale -> TODO: Change this 
        logScale[i] = parameterData[i, "parameterScale"] == "log10" ? true : false
        shouldEst[i] = parameterData[i, "estimate"] == 1 ? true : false
    end
    nParamEst::Int = Int(sum(shouldEst))

    return ParamData(paramVal, lb, ub, paramId, logScale, shouldEst, nParamEst)
end


"""
    processMeasurementData(measurementData::DataFrame, observableData::DataFrame)::MeasurementData

    Process the PeTab measurementData file into a type-stable Julia struct.
"""
function processMeasurementData(measurementData::DataFrame, observableData::DataFrame)::MeasurementData

    # Arrays for storing observed time-points and measurment values (yObs)
    yObs::Array{Float64, 1} = convert(Array{Float64, 1}, measurementData[!, "measurement"])
    tObs::Array{Float64, 1} = convert(Array{Float64, 1}, measurementData[!, "time"])
    nObs = length(yObs)

    # Get the experimental condition ID describing the experimental conditions for each observed time-point. 
    # In case of preequilibration simulation the condition ID is stored in a single-string as the 
    # concatenation of the pre and post equlibration ID:s.
    conditionId::Array{String, 1} = Array{String, 1}(undef, nObs)
    if !("preequilibrationConditionId" in names(measurementData))
        preEq = [missing for i in 1:nObs]
    else
        preEq = measurementData[!, "preequilibrationConditionId"]
    end
    simCond = String.(measurementData[!, "simulationConditionId"])
    for i in eachindex(conditionId)
        if ismissing(preEq[i])
            conditionId[i] = String(simCond[i])
        else
            conditionId[i] = String(preEq[i]) * String(simCond[i])
        end
    end
    if any(x -> ismissing(x), preEq)
        preEq = Array{String, 1}(undef, 0)
    else
        preEq = String.(preEq)
    end

    # PeTab observable ID for each measurment 
    obsID::Array{String, 1} = string.(measurementData[!, "observableId"])

    # Noise parameters in the PeTab file either have a parameter ID, or they have 
    # a value (fixed). Here regardless the values are mapped to the sdParams vector 
    # as string. If sdObs[i] is numeric is the parsed before computing the cost. 
    if !("noiseParameters" in names(measurementData))
        sdObs = [missing for i in 1:nObs]
    else
        sdObs = measurementData[!, "noiseParameters"]
    end
    sdParams::Array{Union{String, Float64}, 1} = Array{Union{String, Float64}, 1}(undef, nObs)
    for i in eachindex(sdObs)
        if ismissing(sdObs[i]) 
            sdParams[i] = ""
        elseif typeof(sdObs[i]) <:AbstractString && isNumber(sdObs[i])
            sdParams[i] = parse(Float64, sdObs[i])
        elseif typeof(sdObs[i]) <:Real
            sdParams[i] = convert(Float64, sdObs[i])
        else
            sdParams[i] = string(sdObs[i])
        end
    end

    # obsParamFile[i] can store more than one parameter. This is parsed
    # when computing the likelihood. 
    if !("observableParameters" in names(measurementData))
        obsParamFile = [missing for i in 1:nObs]
    else
        obsParamFile = measurementData[!, "observableParameters"]
    end
    obsParam = Array{String, 1}(undef, nObs)
    for i in 1:nObs
        if ismissing(obsParamFile[i]) 
            obsParam[i] = ""
        else
            obsParam[i] = string(obsParamFile[i])
        end
    end

    # Currently supports log10 transformation of measurementData
    transformArr = Array{Symbol, 1}(undef, nObs)
    if !("observableTransformation" in names(observableData))
        transformArr .= [:lin for i in 1:nObs]
    else
        for i in 1:nObs
            iRow = findfirst(x -> x == obsID[i], observableData[!, "observableId"])
            transformArr[i] = Symbol(observableData[iRow, "observableTransformation"]) 
        end
    end

    # To avoid repeating calculations yObs is stored in a transformed and untransformed format 
    yObsTransformed::Array{Float64, 1} = deepcopy(yObs)
    transformYobsOrYmodArr!(yObsTransformed, transformArr)

    # To make everything easier for adjoint sensitivity analysis we can sort the observations for each experimental 
    # condition based on the observed time, that is we do not follow the format in the measurementData file. 
    uniqueConditionID = unique(conditionId)
    indexVec::Array{Int64, 1} = Array{Int64, 1}(undef, length(tObs))
    k = 1
    for i in eachindex(uniqueConditionID)
        obsWithId = findall(x -> x == uniqueConditionID[i], conditionId)
        nObsWithId = length(obsWithId)
        iSort = sortperm(tObs[obsWithId])
        indexVec[k:(nObsWithId+k-1)] .= obsWithId[iSort]
        k += nObsWithId
    end
    yObs = yObs[indexVec]
    yObsTransformed = yObsTransformed[indexVec]
    tObs = tObs[indexVec]
    conditionId = conditionId[indexVec]
    sdParams = sdParams[indexVec]
    transformArr = transformArr[indexVec]
    obsParam = obsParam[indexVec]
    simCond = simCond[indexVec]
    obsID = obsID[indexVec]
    if !isempty(preEq)
        preEq = preEq[indexVec]
    end

    # For each experimental condition we want to know the vector of time-points to save the ODE solution at 
    # for each experimental condition. For each t-obs we also want to know which index in t-save vector 
    # it corresponds to. We do not need to sort here as it is done above
    iTimePoint = Array{Int64, 1}(undef, nObs)
    iPerConditionId = Dict{String, Vector{Int64}}() # Index in measurment data corresponding to specific condition id 
    tVecSave = Dict{String, Vector{Float64}}()
    for i in eachindex(uniqueConditionID)
        iConditionId::Vector{Int64} = findall(x -> x == uniqueConditionID[i], conditionId)
        # Sorting is needed so that when extracting time-points when computing the cost 
        # we extract the correct index.
        tVecSave[uniqueConditionID[i]] = sort(unique(tObs[iConditionId]))
        iPerConditionId[uniqueConditionID[i]] = iConditionId
        for j in iConditionId
            iTimePoint[j] = findfirst(x -> x == tObs[j], tVecSave[uniqueConditionID[i]])
        end
    end

    # When doing adjoint sensitivity analysis via the lower level interface we need at each time 
    # points (for a specific experimental condition) need to know how many observations we have, 
    # and which index they have in the measurementData-struct. This computes a dictionary such 
    # that iGroupedTObs["expId"] returns a vector of vector where iGroupedTObs["expId"][i] returns 
    # the indices (in measurmentData) for all observations we have for unique(tVecExpId)[i], which 
    # is then used by the lower level interface for adjoint sensitivity analysis.
    iGroupedTObs = Dict{String, Vector{Vector{Int64}}}()
    for key in keys(iPerConditionId)
        indexObs = iPerConditionId[key]
        uniqueT = unique(tObs[indexObs])
        vecSave = Array{Array{Int64, 1}, 1}(undef, length(uniqueT))
        for j in eachindex(vecSave)
            saveTmp = findall(x -> x == uniqueT[j], tObs)
            vecSave[j] = saveTmp[findall(x -> conditionId[x] == key, saveTmp)]
        end
        iGroupedTObs[key] = vecSave
    end

    return MeasurementData(yObs, yObsTransformed, tObs, obsID, conditionId, sdParams, 
                           transformArr, obsParam, tVecSave, iTimePoint, iPerConditionId, 
                           preEq, simCond, iGroupedTObs)
end


"""
    getTimeMax(measurementData::DataFrame, expId::String)::Float64

    Small helper function to get the time-max value for a specific simulationConditionId when simulating 
    the PeTab ODE-model 
"""
function getTimeMax(measurementData::DataFrame, expId::String)::Float64
    return Float64(maximum(measurementData[findall(x -> x == expId, measurementData[!, "simulationConditionId"]), "time"]))
end


"""
    getSimulationInfo(measurementData::DataFrame)::SimulationInfo

    Using the PeTab measurementData-file extract information on the foward ODE simulations.

    Specifcially extract the experimental ID:s from the experimentalCondition - PeTab file;
    firstExpIds (preequilibration ID:s), the shiftExpIds (postequilibration), and
    simulateSS (whether or not to simulate ODE-model to steady state). Further 
    stores a solArray with the ODE solution where conditionIdSol of the ID for 
    each forward solution

    TODO: Compute t-vec save at from measurementDataFile (instead of providing another struct)
"""
function getSimulationInfo(peTabModel::PeTabModel,
                           measurementDataFile::DataFrame,
                           measurementData::MeasurementData;
                           absTolSS::Float64=1e-8,
                           relTolSS::Float64=1e-6,
                           sensealg::Union{SciMLSensitivity.AbstractForwardSensitivityAlgorithm, SciMLSensitivity.AbstractAdjointSensitivityAlgorithm}=InterpolatingAdjoint())::SimulationInfo 

    # If preequilibrationConditionId column is not empty the model should 
    # first be simulated to a stady state 
    preEqIDs = unique(measurementData.preEqCond)
    simulateSS = length(preEqIDs) > 0

    # In case the the model is simulated to steday state get pre and post equlibration experimental conditions 
    if simulateSS == true
        firstExpIds = preEqIDs
        shiftExpIds = Any[]
        for firstExpId in firstExpIds
            iRows = findall(x -> x == firstExpId, measurementData.preEqCond)
            shiftExpId = unique(measurementData.simCond[iRows])
            push!(shiftExpIds, shiftExpId)
        end
        shiftExpIds = convert(Vector{Vector{String}}, shiftExpIds)
    end

    # In case the the model is not simulated to steday state store experimental condition in firstExpIds
    if simulateSS == false
        firstExpIds = unique(measurementData.simCond)
        shiftExpIds = Array{Array{String, 1}, 1}(undef, 0)
    end

    # Compute number of foward simulations to cover all experimental conditions and allocate array for them
    if simulateSS == true
        nForwardSol = Int64(sum([length(shiftExpIds[i]) for i in eachindex(shiftExpIds)]))
    else
        nForwardSol = Int64(length(firstExpIds))
    end
    # When computing the gradient and hessian the ODE-system needs to be resolved to compute the gradient 
    # of the dynamic parameters, while for the observable/sd parameters the system should not be resolved. 
    # Rather, an ODE solution without dual numbers is required and this solution can be the same which is 
    # used when computing the cost.
    solArray::Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1} = Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1}(undef, nForwardSol)
    solArrayGrad::Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1} = Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1}(undef, nForwardSol)
    if simulateSS
        lenPreEq = length(unique(firstExpIds))
    else
        lenPreEq = 0
    end
    solArrayPreEq::Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1} = Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1}(undef, lenPreEq)

    # Array with conition-ID for each foward simulations. As we always solve the ODE in the same order this can 
    # be pre-computed.
    conditionIdSol::Array{String, 1} = Array{String, 1}(undef, nForwardSol)
    preEqIdSol::Array{String, 1} = Array{String, 1}(undef, nForwardSol)
    postEqIdSol::Array{String, 1} = Array{String, 1}(undef, nForwardSol)
    tMaxForwardSim::Array{Float64, 1} = Array{Float64, 1}(undef, nForwardSol)
    if simulateSS == true
        k = 1
        for i in eachindex(firstExpIds)
            for j in eachindex(shiftExpIds[i])
                firstExpId = firstExpIds[i]
                shiftExpId = shiftExpIds[i][j]
                conditionIdSol[k] = firstExpId * shiftExpId
                preEqIdSol[k] = firstExpId
                postEqIdSol[k] = shiftExpId
                tMaxForwardSim[k] = getTimeMax(measurementDataFile, shiftExpId)
                k +=1
            end
        end
    else
        for i in eachindex(firstExpIds)    
            firstExpId = firstExpIds[i]
            conditionIdSol[i] = firstExpId
            tMaxForwardSim[i] = getTimeMax(measurementDataFile, firstExpId)
        end
    end

    # Some models, e.g those with time dependent piecewise statements, have callbacks encoded. To accurately 
    # solve these we have callbacks. When doing adjoint sensitivity analysis we need to track these callbacks 
    # for the forward and reverse solution, hence we store within an arrary. 
    modelCallbacks::Vector{SciMLBase.DECallback} = [deepcopy(peTabModel.callbackSet) for i in 1:nForwardSol]

    tVecSave::Dict{String, Vector{Float64}} = deepcopy(measurementData.tVecSave)
    # Set up map to identify for each experimental condition in which order (when solving all experimental condition) a time 
    # point occurs at. 
    posInSolArray = getPosInSolArray(conditionIdSol, tVecSave)
    simulationInfo = SimulationInfo(firstExpIds, 
                                    shiftExpIds, 
                                    preEqIdSol,
                                    postEqIdSol,
                                    conditionIdSol, 
                                    tMaxForwardSim,
                                    simulateSS,
                                    solArray, 
                                    solArrayGrad, 
                                    solArrayPreEq,
                                    absTolSS, 
                                    relTolSS, 
                                    tVecSave, 
                                    posInSolArray,
                                    modelCallbacks, 
                                    sensealg)
    return simulationInfo
end


"""
    checkForPeTabFile(fileSearchFor::String, dirModel::String)::String

    Helper function to check in dirModel if a file starting with fileSearchFor exists. 
    If true return file path.
"""
function checkForPeTabFile(fileSearchFor::String, dirModel::String)::String

    filesDirModel = readdir(dirModel)
    iUse = findall(x -> occursin(fileSearchFor, x), filesDirModel)
    if length(iUse) > 1 
        @printf("Error : More than 1 file starting with %s in %s\n", fileSearchFor, filesDirModel)
    end
    if length(iUse) == 0
        @printf("Error : No file starting with %s in %s\n", fileSearchFor, filesDirModel)
    end

    return dirModel * filesDirModel[iUse[1]]
end


function getPriorInfo(paramEstIndices::ParameterIndices, parameterDataFile::DataFrame)::PriorInfo

    if "objectivePriorType" ∉ names(parameterDataFile)
        return PriorInfo(Array{Function, 1}(undef, 0), Bool[], false)
    end

    namesParamEst = paramEstIndices.namesParamEst
    priorLogPdf = Array{Function, 1}(undef, length(namesParamEst))
    priorOnParamScale = Array{Bool, 1}(undef, length(namesParamEst))
    paramID = string.(parameterDataFile[!, "parameterId"])

    contPrior = 0.0
    for i in eachindex(namesParamEst)

        iUse = findfirst(x -> x == namesParamEst[i], paramID)

        priorF = parameterDataFile[iUse, "objectivePriorType"]
        if ismissing(priorF)
            priorLogPdf[i] = noPrior
            priorOnParamScale[i] = false
            continue
        end

        priorVal = parse.(Float64, split(parameterDataFile[iUse, "objectivePriorParameters"], ";"))

        if priorF == "parameterScaleNormal"
            contPrior += logpdf(Normal(priorVal[1], priorVal[2]), log10(parameterDataFile[iUse, "nominalValue"]))
            priorLogPdf[i] = (p) -> logpdf(Normal(priorVal[1], priorVal[2]), p) 
            priorOnParamScale[i] = true
        elseif priorF == "parameterScaleLaplace"
            priorLogPdf[i] = (p) -> logpdf(Laplace(priorVal[1], priorVal[2]), p) 
            priorOnParamScale[i] = true
        elseif priorF == "normal"
            priorLogPdf[i] = (p) -> logpdf(Normal(priorVal[1], priorVal[2]), p) 
            priorOnParamScale[i] = false
        elseif priorF == "laplace"
            priorLogPdf[i] = (p) -> logpdf(Laplace(priorVal[1], priorVal[2]), p) 
            priorOnParamScale[i] = false
        elseif priorF == "logNormal"
            priorLogPdf[i] = (p) -> logpdf(LogNormal(priorVal[1], priorVal[2]), p) 
            priorOnParamScale[i] = false
        elseif priorF == "logLaplace"
            println("Error : Julia does not yet have support for log-laplace")
        else
            println("Error : PeTab standard does not support a prior of type ", priorF)
        end

    end

    return PriorInfo(priorLogPdf, priorOnParamScale, true)
end
# Helper function in case there is not any parameter priors 
function noPrior(p::Real)::Real
    return 0.0
end


# Function generating callbacksets for time-depedent SBML piecewise expressions
function getCallbacksForTimeDepedentPiecewise(odeSys::ODESystem, modelDict::Dict, modelName::String, dirModel::String)

    # ParamEstIndices is needed to see if event-triggers contain parameters we want to estimate 
    experimentalConditionsFile, measurementDataFile, parameterDataFile, observablesDataFile = readDataFiles(dirModel, readObs=true)
    parameterData = processParameterData(parameterDataFile)
    measurementData = processMeasurementData(measurementDataFile, observablesDataFile) 
    paramEstIndices = getIndicesParam(parameterData, measurementData, odeSys, experimentalConditionsFile)

    parameterNames = parameters(odeSys)
    stateNames = string.(states(odeSys))
    parameterNames = string.(parameters(odeSys))
    
    stringWrite = "function getCallbacks_" * modelName * "()\n"
    stringWriteFunctions = ""
    stringWriteTstops = "function getTstops(u, p)\n"

    # In case we do not have any events 
    if isempty(modelDict["boolVariables"])
        callbackNames = ""
        checkActivatedNames = ""
        stringWriteTstops *= "\t return Float64[]\nend"
    else
        # In case we have events loop over each variable 
        for key in keys(modelDict["boolVariables"])
            stringFunctions, stringSet =  createCallbackCont(key, modelDict, parameterNames, stateNames, paramEstIndices) 
            stringWrite *= stringSet * "\n"
            stringWriteFunctions *= stringFunctions * "\n"
        end
        callbackNames = prod(["cb_" * key * ", " for key in keys(modelDict["boolVariables"])])[1:end-2]
        checkActivatedNames = prod(["activeAtTime0_" * key * "!, " for key in keys(modelDict["boolVariables"])])[1:end-2]

        stringWriteTstops *= "\t return " * createFuncForTstops(modelDict, stateNames, parameterNames) * "\nend"
    end


    # Put everything together
    stringWrite *= "\treturn CallbackSet(" * callbackNames * "), [" * checkActivatedNames * "]\nend"
    fileWrite = dirModel * "/" * modelName * "Callbacks_time_piecewise.jl"
    if isfile(fileWrite)
        rm(fileWrite)
    end
    io = open(fileWrite, "w")    

    write(io, stringWrite * "\n\n")
    write(io, stringWriteFunctions)
    write(io, stringWriteTstops)
    
    close(io)
end


function createCallbackCont(paramName::String, 
                            modelDict::Dict, 
                            parameterNames::Vector{String}, 
                            stateNames::Vector{String}, 
                            paramEstIndices::ParameterIndices)

    stringCallbackFunctions = ""
    stringCallbackSet = ""

    # Build the condition statement 
    stringCallbackFunctions *= "\nfunction condition_" * paramName * "(u, t, integrator)\n"
    activationFormula = modelDict["boolVariables"][paramName][1]
    sideActivated = modelDict["boolVariables"][paramName][2]

    # Check if the event trigger depend on parameters which are to be i) estimated, or ii) if it depend on models state. 
    # For i) it must be a cont. event in order for us to be able to compute the gradient. For ii) we cannot compute 
    # tstops (the event times) prior to starting to solve the ODE.
    hasStates = triggerHasStates(activationFormula, stateNames)
    hasParametersEst = triggerHasParametersToEst(activationFormula, parameterNames, paramEstIndices)
    discreteEvent = true
    if hasParametersEst == true || hasStates == true
        discreteEvent = false
    end
    
    # Replace any state or parameter with their corresponding index in the ODE system 
    for i in eachindex(stateNames)
        activationFormula = replaceWholeWord(activationFormula, stateNames[i], "u["*string(i)*"]")
    end
    for i in eachindex(parameterNames)
        activationFormula = replaceWholeWord(activationFormula, parameterNames[i], "integrator.p["*string(i)*"]")
    end

    # Replace inequality 
    if discreteEvent == false
        activationFormulaCond = replace(activationFormula, "<=" => "-")
        activationFormulaCond = replace(activationFormulaCond, ">=" => "-")
        activationFormulaCond = replace(activationFormulaCond, ">" => "-")
        activationFormulaCond = replace(activationFormulaCond, "<" => "-")
    else
        activationFormulaCond = replace(activationFormula, "<=" => "==")
        activationFormulaCond = replace(activationFormulaCond, ">=" => "==")
        activationFormulaCond = replace(activationFormulaCond, ">" => "==")
        activationFormulaCond = replace(activationFormulaCond, "<" => "==")
    end
    stringCallbackFunctions *= "\t" * activationFormulaCond * "\nend\n"

    # Build the affect function 
    whichParam = findfirst(x -> x == paramName, parameterNames)
    stringCallbackFunctions *= "function affect_" * paramName * "!(integrator)\n"
    stringCallbackFunctions *= "\tintegrator.p[" * string(whichParam) * "] = 1.0\nend\n"

    # Build the callback 
    if discreteEvent == false
        stringCallbackSet *= "\tcb_" * paramName * " = ContinuousCallback(" * "condition_" * paramName * ", " * "affect_" * paramName * "!, "
    else
        stringCallbackSet *= "\tcb_" * paramName * " = DiscreteCallback(" * "condition_" * paramName * ", " * "affect_" * paramName * "!, "
    end
    stringCallbackSet *= "save_positions=(false, false))\n"

    # Build a function to check if a condition should be true at time zero (as in Julia events activated at time zero 
    # are usually not triggered by default)
    sideInequality = sideActivated == "right" ? "!" : ""
    stringCallbackFunctions *= "function activeAtTime0_" * paramName * "!(u, p)\n"
    stringCallbackFunctions *= "\tt = 0.0 # Used to check conditions activated at t0=0\n"
    stringCallbackFunctions *= "\tp[" * string(whichParam) * "] = 0.0 # Default to being off\n"
    activationFormulaBool = replace(activationFormula, "integrator." => "")
    activationFormulaBool = replace(activationFormulaBool, "<=" => "≤")
    activationFormulaBool = replace(activationFormulaBool, ">=" => "≥")
    stringCallbackFunctions *= "\tif " * sideInequality *"(" * activationFormulaBool * ")\n"
    stringCallbackFunctions *= "\t\tp[" * string(whichParam) * "] = 1.0\n\tend\nend\n"

    return stringCallbackFunctions, stringCallbackSet
end


function triggerHasStates(activationFormula::AbstractString, stateNames::Vector{String})::Bool
    for i in eachindex(stateNames)
        activationFormulaNew = replaceWholeWord(activationFormula, stateNames[i], "u["*string(i)*"]")
        if activationFormulaNew != activationFormula
            return true
        end
    end
    return false
end


function triggerHasParametersToEst(activationFormula::AbstractString, parameterNames::Vector{String}, paramEstIndices::ParameterIndices)::Bool

    # Parameters which are present for each experimental condition, and condition specific parameters 
    indexParamEstODESys = paramEstIndices.mapDynParEst.iDynParamInSys
    indexParamEstODESysCond = reduce(vcat, [paramEstIndices.mapExpCond[i].iOdeProbDynParam for i in eachindex(paramEstIndices.mapExpCond)])

    for i in eachindex(parameterNames)
        activationFormulaNew = replaceWholeWord(activationFormula, parameterNames[i], "integrator.p["*string(i)*"]")
        if activationFormulaNew != activationFormula
            if i ∈ indexParamEstODESys || i ∈ indexParamEstODESysCond
                return true
            end
        end
    end

    return false
end


# Function computing t-stops (time for events) for piecewise expressions using the symbolics package 
function createFuncForTstops(modelDict::Dict, stateNames::Vector{String}, parameterNames::Vector{String})

    tStopExpressions = Array{String, 1}(undef, length(keys(modelDict["boolVariables"])))

    i = 1
    for key in keys(modelDict["boolVariables"])

        condFormula = modelDict["boolVariables"][key][1]
        # In case the activation formula contains a state we cannot precompute the t-stop time as it depends on 
        # the actual ODE solution.
        if triggerHasStates(condFormula, stateNames) == true
            i += 1
            continue
        end
        
        stringVariables = "@variables t, "
        stringVariables *= prod(string.(collect(keys(modelDict["parameters"]))) .* ", " )[1:end-2] * " "
        stringVariables *= prod(string.(collect(keys(modelDict["states"]))) .* ", " )[1:end-2] 
        symbolicVariables = eval(Meta.parse(stringVariables))

        # Note - below order counts (e.g having < first results in ~= incase what actually stands is <=)
        condFormula = replace(condFormula, "<=" => "~")
        condFormula = replace(condFormula, ">=" => "~")
        condFormula = replace(condFormula, "<" => "~")
        condFormula = replace(condFormula, ">" => "~")
        expSymbolic = eval(Meta.parse(condFormula))

        # Expression for the time at which the condition is triggered
        expForTime = string.(Symbolics.solve_for(expSymbolic, symbolicVariables[1], simplify=true))

        for i in eachindex(stateNames)
            expForTime = replaceWholeWord(expForTime, stateNames[i], "u["*string(i)*"]")
        end
        for i in eachindex(parameterNames)
            expForTime = replaceWholeWord(expForTime, parameterNames[i], "p["*string(i)*"]")
        end
        tStopExpressions[i] = "dualToFloat(" * expForTime * ")"
        i += 1
    end

    stringTstop = "[" * prod([str * ", " for str in tStopExpressions])[1:end-2] * "]"
    return stringTstop
end


# Function taking a Dual and making it into the underlaying float via recursion. This is needed 
# if having events where t-stops depend on some of the parameter in odeProb.p, as here even though 
# the parameter might be constant it can still be a Dual.
function dualToFloat(x::AbstractFloat)
    return x
end
function dualToFloat(x::ForwardDiff.Dual)
    return dualToFloat(x.value)
end