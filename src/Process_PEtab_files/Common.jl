"""
    readDataFiles(dirModel::String; readObs::Bool=false)

    Given a directory for a model, e.g ./Beer_MolBioSystems2014, read the associated PeTab files 
    for the measurements, parameters, experimental conditions and (if true) the observables.
"""
function readPEtabFiles(dirModel::String; readObservables::Bool=false)

    # Check if PeTab files exist and get their path 
    pathMeasurements = checkForPeTabFile("measurementData", dirModel)
    pathExperimentalCondidtions = checkForPeTabFile("experimentalCondition", dirModel)
    pathParameters = checkForPeTabFile("parameters", dirModel)

    experimentalConditions = CSV.read(pathExperimentalCondidtions, DataFrame)
    measurementsData = CSV.read(pathMeasurements, DataFrame)
    parametersData = CSV.read(pathParameters, DataFrame)
    if readObservables == true
        pathObservables = checkForPeTabFile("observables", dirModel)
        observablesData = CSV.read(pathObservables, DataFrame)
        return experimentalConditions, measurementsData, parametersData, observablesData
    else
        return experimentalConditions, measurementsData, parametersData
    end
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
