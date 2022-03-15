function getModelFiles(path)
    if isdir(path)
        modelFiles = readdir(path)
        return modelFiles
    else
        println("No such directory")
        return nothing
    end
end

function getTolerances(; onlyMaxTol = false)
    tolList = [1e-6, 1e-9, 1e-12] # 1e-6 is standard
    if onlyMaxTol
        return [tolList[1]]
    else
        return tolList
    end
end

function fixDirectories(path)
    if ~isdir(path)
        println("Create directory: " * path)
        mkdir(path)
    end
end

function getNumberOfFiles(path)
    if isdir(path)
        files = readdir(path)
        return length(files)
    else
        mkdir(path)
        return 0
    end
end

function getAutojacvecs()
    autojacvecs = [ReverseDiffVJP(true),
                   ReverseDiffVJP(),
                   ZygoteVJP(),
                   TrackerVJP(),
                   EnzymeVJP(),
                   nothing,
                   true,
                   false]
    return autojacvecs
end

function getSenseAlgs()
    autojacvecs = getAutojacvecs()
    senseAlgs = [[BacksolveAdjoint(autojacvec = autojacvecs[i]) for i in 1:length(autojacvecs)],
                 [InterpolatingAdjoint(checkpointing=true, autojacvec = autojacvecs[i]) for i in 1:length(autojacvecs)], 
                 [InterpolatingAdjoint(autojacvec = autojacvecs[i]) for i in 1:length(autojacvecs)], 
                 [QuadratureAdjoint(autojacvec = autojacvecs[i]) for i in 1:length(autojacvecs)],
                 [ReverseDiffAdjoint()],
                 [TrackerAdjoint()],
                 [ZygoteAdjoint()],
                 [SensitivityADPassThrough()]
                 ]
    return senseAlgs
end