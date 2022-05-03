using PyCall, DataFrames, CSV

include(joinpath(pwd(), "Pipeline_SBMLImporter", "write_modellingToolkit_to_file.jl"))

function getSBMLFiles(path)
    files = readdir(path)
    return files
end

function fixDirectories(path)
    if ~isdir(path)
        mkdir(path)
    end
end


function main()
    libsbml = pyimport("libsbml")
    reader = libsbml.SBMLReader()
    readPath = joinpath(pwd(), "Pipeline_SBMLImporter", "SBML")
    writePath = joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels")
    fixDirectories(writePath)

    files = getSBMLFiles(readPath)
    for file in files
        modelName = file[1:end-4]
        modelNameShort = modelName[7:end]
        dataEnding = modelNameShort * ".tsv"
        readDataPath = joinpath(pwd(), "Pipeline_ModelParameterEstimation", "Data", modelName)
        experimentalConditions = CSV.read(joinpath(readDataPath, "experimentalCondition_" * dataEnding), DataFrame)
        parameterBounds = CSV.read(joinpath(readDataPath, "parameters_" * dataEnding), DataFrame)

        document = reader[:readSBML](joinpath(pwd(), "Pipeline_SBMLImporter", "SBML", file))
        model = document[:getModel]() # Get the model
        writeODEModelToFile(libsbml, model, modelName, writePath, experimentalConditions, parameterBounds)
    end
end

main()