using PyCall

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
        document = reader[:readSBML](joinpath(pwd(), "Pipeline_SBMLImporter", "SBML", file))
        model = document[:getModel]() # Get the model
        writeODEModelToFile(libsbml, model, file[1:end-4], writePath)
    end
end

main()