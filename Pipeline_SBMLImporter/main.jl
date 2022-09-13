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


function main(; usedFiles = ["all"]::Vector{String}, useData = false, wrapped = true )
    libsbml = pyimport("libsbml")
    reader = libsbml.SBMLReader()
    readPath = joinpath(pwd(), "Pipeline_SBMLImporter", "SBML")
    writePath = joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels")
    fixDirectories(writePath)

    if usedFiles == ["all"]
        files = getSBMLFiles(readPath)
    else
        files = usedFiles .* ".xml"
    end
    files = files[files .!= "model_Chen_MSB2009.xml"]

    if useData
        for file in files
            modelName = file[1:end-4]
            modelNameShort = modelName[7:end]
            dataEnding = modelNameShort * ".tsv"
            readDataPath = joinpath(pwd(), "Pipeline_ModelParameterEstimation", "Data", modelName)

            # Check for experimental conditions file and parameter bounds (why?)
            pathExpCond = joinpath(readDataPath, "experimentalCondition_" * dataEnding)
            if isfile(pathExpCond)
                experimentalConditions = CSV.read(pathExpCond, DataFrame)
            else
                experimentalConditions = []
            end
            pathBounds = joinpath(readDataPath, "parameters_" * dataEnding)
            if isfile(pathBounds) 
                parameterBounds = CSV.read(pathBounds, DataFrame)
            else
                parameterBounds = []
                useData = false
            end
            document = reader[:readSBML](joinpath(pwd(), "Pipeline_SBMLImporter", "SBML", file))
            model = document[:getModel]() # Get the model
            writeODEModelToFile(libsbml, model, modelName, writePath, useData, wrapped, experimentalConditions = experimentalConditions, parameterBounds = parameterBounds)
        end
    else
        for file in files
            modelName = file[1:end-4]

            document = reader[:readSBML](joinpath(pwd(), "Pipeline_SBMLImporter", "SBML", file))
            model = document[:getModel]() # Get the model
            writeODEModelToFile(libsbml, model, modelName, writePath, useData, wrapped)
        end
    end
end

main(usedFiles = ["all"], useData = true, wrapped = true)