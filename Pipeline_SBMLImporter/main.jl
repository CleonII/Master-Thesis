using PyCall

include(pwd() * "\\Pipeline_SBMLImporter\\write_modellingToolkit_to_file.jl")

function getSBMLFiles(path)
    files = readdir(path)
end

function fixDirectories(path)
    if ~isdir(path)
        mkdir(path)
    end
end


function main()
    libsbml = pyimport("libsbml")
    reader = libsbml.SBMLReader()
    readPath = pwd() * "\\Pipeline_SBMLImporter\\SBML"
    writePath = pwd() * "\\Pipeline_SBMLImporter\\JuliaModels"
    fixDirectories(writePath)
    

    files = getSBMLFiles(readPath)
    for file in files
        document = reader[:readSBML](pwd() * "\\Pipeline_SBMLImporter\\SBML\\" * file)
        model = document[:getModel]() # Get the model
        writeODEModelToFile(model, file[1:end-4], writePath)
    end
end

main()