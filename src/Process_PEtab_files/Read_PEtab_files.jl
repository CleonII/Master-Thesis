function readPEtabYamlFile(pathYAML::AbstractString)

    if !isfile(pathYAML)
        throw(PEtabFileError("Model YAML file does not exist in the model directory"))
    end
    fileYAML = YAML.load_file(pathYAML)

    dirModel = dirname(pathYAML)

    pathSBML = joinpath(dirModel, fileYAML["problems"][1]["sbml_files"][1])
    if !isfile(pathSBML)
        throw(PEtabFileError("SBML file does not exist in the model directory"))
    end

    pathMeasurements = joinpath(dirModel, fileYAML["problems"][1]["measurement_files"][1])
    if !isfile(pathMeasurements)
        throw(PEtabFileError("Measurements file does not exist in the model directory"))
    end

    pathObservables = joinpath(dirModel, fileYAML["problems"][1]["observable_files"][1])
    if !isfile(pathObservables)
        throw(PEtabFileError("Observables file does not exist in the models directory"))
    end

    pathConditions = joinpath(dirModel, fileYAML["problems"][1]["condition_files"][1])
    if !isfile(pathConditions)
        throw(PEtabFileError("Conditions file does not exist in the models directory"))
    end

    pathParameters = joinpath(dirModel, fileYAML["parameter_file"])
    if !isfile(pathParameters)
        throw(PEtabFileError("Parameter file does not exist in the models directory"))
    end

    # Extract YAML directory and use directory name as model name and build directory for Julia files
    dirJulia = joinpath(dirModel, "Julia_model_files") 
    modelName = splitdir(dirModel)[end]
    if !isdir(dirJulia)
        mkdir(dirJulia)
    end

    return pathSBML, pathParameters, pathConditions, pathObservables, pathMeasurements, dirJulia, dirModel, modelName
end


function readPEtabFiles(pathYAML::String)

    pathSBML, pathParameters, pathConditions, pathObservables, pathMeasurements, dirJulia, dirModel, modelName = readPEtabYamlFile(pathYAML)

    experimentalConditions = CSV.read(pathConditions, DataFrame, stringtype=String)
    measurementsData = CSV.read(pathMeasurements, DataFrame, stringtype=String)
    parametersData = CSV.read(pathParameters, DataFrame, stringtype=String)
    observablesData = CSV.read(pathObservables, DataFrame, stringtype=String)
    
    return experimentalConditions, measurementsData, parametersData, observablesData
end
function readPEtabFiles(petabModel::PEtabModel)

    experimentalConditions = CSV.read(petabModel.pathConditions, DataFrame, stringtype=String)
    measurementsData = CSV.read(petabModel.pathMeasurements, DataFrame, stringtype=String)
    parametersData = CSV.read(petabModel.pathParameters, DataFrame, stringtype=String)
    observablesData = CSV.read(petabModel.pathObservables, DataFrame, stringtype=String)
    
    return experimentalConditions, measurementsData, parametersData, observablesData
end