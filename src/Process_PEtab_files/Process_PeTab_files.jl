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
