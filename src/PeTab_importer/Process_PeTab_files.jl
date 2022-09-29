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
function setUpPeTabModel(modelName::String, dirModel::String; forceBuildJlFile::Bool=false)::PeTabModel

    # Sanity check user input 
    modelFileXml = dirModel * modelName * ".xml"
    modelFileJl = dirModel * modelName * ".jl"
    if !isdir(dirModel)
        @printf("Model directory %s does not exist\n", dirModel)
    end
    if !isfile(modelFileXml)
        @printf("Model directory does not contain xml-file with name %s\n", modelName * "xml")
    end
    # If Julia model file does exists build it 
    if !isfile(modelFileJl) && forceBuildJlFile == false
        @printf("Julia model file does not exist - will build it\n")
        XmlToModellingToolkit(modelFileXml, modelName, dirModel)
    elseif isfile(modelFileJl) && forceBuildJlFile == false
        @printf("Julia model file exists at %s - will not rebuild it\n", modelFileJl)
    elseif forceBuildJlFile == true
        @printf("By user option rebuilds Julia model file\n")
        if isfile(modelFileJl)
            rm(modelFileJl)
        end
        XmlToModellingToolkit(modelFileXml, modelName, dirModel)
    end

    # Extract ODE-system and mapping of maps of how to map parameters to states and model parmaeters 
    include(modelFileJl)
    expr = Expr(:call, Symbol("getODEModel_" * modelName))
    odeSys, stateMap, paramMap = eval(expr)
    odeSysUse = ode_order_lowering(odeSys)
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
    if !isfile(pathObsSdU0)
        @printf("File for yMod, U0 and Sd does not exist - building it\n")
        createFileYmodSdU0(modelName, dirModel, odeSysUse, stateMap)
    else
        @printf("File for yMod, U0 and Sd does exist - will not rebuild it\n")
    end
    include(pathObsSdU0)

    peTabModel = PeTabModel(modelName,
                            evalYmod,
                            evalU0!,
                            evalSd!,
                            odeSysUse,
                            paramMap,
                            stateMap,
                            parameterNames, 
                            stateNames,
                            dirModel,
                            pathMeasurementData,
                            pathExperimentalCond,
                            pathObservables, 
                            pathParameters)

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

    # Get the experimental condition ID describing the experimental conditiosn for each observed time-point. 
    # In case of preequilibration simulation the condition ID is stored in a single-string as the 
    # concatenation of the pre and post equlibration ID:s.
    conditionId::Array{String, 1} = Array{String, 1}(undef, nObs)
    if !("preequilibrationConditionId" in names(measurementData))
        preEq = [missing for i in 1:nObs]
    else
        preEq = measurementData[!, "preequilibrationConditionId"]
    end
    simCond = measurementData[!, "simulationConditionId"]
    for i in eachindex(conditionId)
        if ismissing(preEq[i])
            conditionId[i] = String(simCond[i])
        else
            conditionId[i] = String(preEq[i]) * String(simCond[i])
        end
    end

    # PeTab observable ID for each measurment 
    obsID::Array{String, 1} = string.(measurementData[!, "observableId"])

    # Noise parameters in the PeTab file either have a parameter ID, or they have 
    # a value (fixed). Here regardless the values are mapped to the sdParams vector 
    # as string. If sdObs[i] is numeric is the parsed before computing the cost. 
    if !("noiseParameters" in names(measurementData))
        sdObs = [missing for i in 1:nObs]
    else
        sdObs = string.(measurementData[!, "noiseParameters"])
    end
    sdParams::Array{String, 1} = Array{String, 1}(undef, nObs)
    for i in eachindex(sdObs)
        if ismissing(sdObs[i])
            sdParams[i] = ""
        else
            sdParams[i] = sdObs[i]
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
            obsParam[i] = obsParamFile[i]
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

    return MeasurementData(yObs, yObsTransformed, tObs, obsID, conditionId, sdParams, transformArr, obsParam)
end



"""
    getSimulationInfo(measurementData::DataFrame)::SimulationInfo

    Using the PeTab measurementData-file extract information on the foward ODE simulations.

    Specifcially extract the experimental ID:s from the experimentalCondition - PeTab file;
    firstExpIds (preequilibration ID:s), the shiftExpIds (postequilibration), and
    simulateSS (whether or not to simulate ODE-model to steady state). Further 
    stores a solArray with the ODE solution where conditionIdSol of the ID for 
    each forward solution
"""
function getSimulationInfo(measurementData::DataFrame;
                           absTolSS::Float64=1e-8,
                           relTolSS::Float64=1e-6)::SimulationInfo

    # If preequilibrationConditionId column is not empty the model should 
    # first be simulated to a stady state 
    colNames = names(measurementData)
    if !("preequilibrationConditionId" in colNames)
        preEqIDs = Array{String, 1}(undef, 0)
    else
        preEqIDs = convert(Array{String, 1}, unique(filter(x -> !ismissing(x), measurementData[!, "preequilibrationConditionId"])))
    end
    simulateSS = length(preEqIDs) > 0

    # In case the the model is simulated to steday state get pre and post equlibration experimental conditions 
    if simulateSS == true
        firstExpIds = preEqIDs
        shiftExpIds = Any[]
        for firstExpId in firstExpIds
            iRows = findall(x -> x == firstExpId, measurementData[!, "preequilibrationConditionId"])
            shiftExpId = unique(measurementData[iRows, "simulationConditionId"])
            push!(shiftExpIds, shiftExpId)
        end
        shiftExpIds = convert(Vector{Vector{String}}, shiftExpIds)
    end

    # In case the the model is mpt simulated to steday state store experimental condition in firstExpIds
    if simulateSS == false
        firstExpIds = convert(Array{String, 1}, unique(measurementData[!, "simulationConditionId"]))
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
    solArray = Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1}(undef, nForwardSol)
    solArrayGrad = Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1}(undef, nForwardSol)

    # Array with conition-ID for each foward simulations 
    conditionIdSol = Array{String, 1}(undef, nForwardSol)

    simulationInfo = SimulationInfo(firstExpIds, 
                                    shiftExpIds, 
                                    conditionIdSol, 
                                    simulateSS,
                                    solArray, 
                                    solArrayGrad, 
                                    absTolSS, 
                                    relTolSS)
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