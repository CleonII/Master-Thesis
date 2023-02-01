function readPEtabFiles(pathYAML::String)

    pathSBML, pathParameters, pathConditions, pathObservables, pathMeasurements, dirJulia, dirModel, modelName = readPEtabYamlFile(pathYAML)

    experimentalConditions = CSV.read(pathConditions, DataFrame)
    measurementsData = CSV.read(pathMeasurements, DataFrame)
    parametersData = CSV.read(pathParameters, DataFrame)
    observablesData = CSV.read(pathObservables, DataFrame)
    
    return experimentalConditions, measurementsData, parametersData, observablesData
end
function readPEtabFiles(petabModel::PEtabModel)

    experimentalConditions = CSV.read(petabModel.pathConditions, DataFrame)
    measurementsData = CSV.read(petabModel.pathMeasurements, DataFrame)
    parametersData = CSV.read(petabModel.pathParameters, DataFrame)
    observablesData = CSV.read(petabModel.pathObservables, DataFrame)
    
    return experimentalConditions, measurementsData, parametersData, observablesData
end


"""
    isNumber(x::String)::Bool

    Check if a string x is a number (Float) taking sciencetific notation into account.
"""
function isNumber(x::AbstractString)::Bool
    re1 = r"^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)$" # Picks up scientific notation
    re2 = r"^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)$"
    return (occursin(re1, x) || occursin(re2, x))
end
"""
    isNumber(x::SubString{String})::Bool
"""
function isNumber(x::SubString{String})::Bool
    re1 = r"^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)$" # Picks up scientific notation
    re2 = r"^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)$"
    return (occursin(re1, x) || occursin(re2, x))
end
