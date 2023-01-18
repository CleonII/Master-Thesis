function getNDynParam(peTabModel::PeTabModel)
    experimentalConditionsFile, measurementDataFile, parameterDataFile, observablesDataFile = readDataFiles(peTabModel.dirModel, readObs=true)
    parameterData = processParameterData(parameterDataFile)
    measurementData = processMeasurementData(measurementDataFile, observablesDataFile) 
    # Indices for mapping parameter-estimation vector to dynamic, observable and sd parameters correctly when calculating cost
    paramEstIndices = getIndicesParam(parameterData, measurementData, peTabModel.odeSystem, experimentalConditionsFile)
    return length(paramEstIndices.namesDynParam)
end


# For a PEtab model create a tempory directory where the PEtab files is changed such that nParamFixate are 
# fixated.
function getPEtabModelNparamFixed(peTabModel::PeTabModel, nParamFixate::Integer)::PeTabModel
    
    modelName = peTabModel.modelName
    dirNew = peTabModel.dirModel * "/Fewer_param/"
    if !isdir(dirNew)
        mkdir(dirNew)
    end
    # Copy over PEtab-files 
    cp(peTabModel.pathParameters, dirNew * "parameters_" * modelName * ".tsv", force=true) 
    cp(peTabModel.pathObservables, dirNew * "observables_" * modelName * ".tsv", force=true)
    cp(peTabModel.pathMeasurementData, dirNew * "measurementData_" * modelName * ".tsv", force=true)
    cp(peTabModel.pathExperimentalConditions, dirNew * "experimentalCondition_" * modelName * ".tsv", force=true)
    cp(dirModel * modelName * ".xml", dirNew * modelName * ".xml", force=true)

    # Fixate model parameters in PEtab file 
    experimentalConditionsFile, measurementDataFile, parameterDataFile, observablesDataFile = readDataFiles(dirNew, readObs=true)
    parameterData = processParameterData(parameterDataFile)
    measurementData = processMeasurementData(measurementDataFile, observablesDataFile) 
    # Indices for mapping parameter-estimation vector to dynamic, observable and sd parameters correctly when calculating cost
    paramEstIndices = getIndicesParam(parameterData, measurementData, peTabModel.odeSystem, experimentalConditionsFile)

    paramFixate = sample(paramEstIndices.namesDynParam, nParamFixate)
    for i in 1:nrow(parameterDataFile)
        if parameterDataFile[i, :parameterId] ∈ paramFixate
            parameterDataFile[i, :estimate] = 0
        end
    end
    rm(dirNew * "parameters_" * modelName * ".tsv")
    CSV.write(dirNew * "parameters_" * modelName * ".tsv", parameterDataFile, delim = '\t')

    return setUpPeTabModel(modelName, dirNew)
end
