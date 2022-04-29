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

function fixDirectories(paths::T1) where T1 <: Union{String, AbstractVector{String}}
    if typeof(paths) <: AbstractVector
        for path in paths
            if ~isdir(path)
                println("Create directory: " * path)
                mkdir(path)
            end
        end
    else
        path = paths
        if ~isdir(path)
            println("Create directory: " * path)
            mkdir(path)
        end
    end

    nothing
end

function getNumberOfFiles(path::String)::Int64
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
                   EnzymeVJP()]
    return autojacvecs
end

function getSenseAlgs()
    autojacvecs = getAutojacvecs()
    senseAlgs = [[BacksolveAdjoint(autojacvec = autojacvecs[i]) for i in 1:length(autojacvecs)],
                 [InterpolatingAdjoint(checkpointing=true, autojacvec = autojacvecs[i]) for i in 1:length(autojacvecs)], 
                 [InterpolatingAdjoint(autojacvec = autojacvecs[i]) for i in 1:length(autojacvecs)], 
                 [InterpolatingAdjoint(autodiff = true, autojacvec = false)], 
                 [QuadratureAdjoint(autojacvec = autojacvecs[i]) for i in 1:length(autojacvecs)]]
    senseAlgs = reduce(vcat, senseAlgs)
    return senseAlgs
end

function get_type(lists)
    type = Float64
    if eltype(lists) <: AbstractVector
        for list in lists
            if eltype(list) !== Float64
                type = eltype(list)
                break
            end
        end
    else 
        type = eltype(lists)
    end
    return type
end

function updateAllParameterVectors_proto(modelParameters, modelData)
    allParameters = modelParameters.allParameters
    parameterIndices = modelParameters.parameterIndices
    pars = view(allParameters, parameterIndices)
    dynPar = modelParameters.dynamicParametersVector
    dynPar[modelData.optParameterIndices] .= pars

    scaleIndices = modelParameters.scaleIndices
    scale = view(allParameters, scaleIndices)
    scaleVector = modelParameters.scaleVector
    scaleVector[1:modelParameters.numScale] .= scale

    offsetIndices = modelParameters.offsetIndices
    offset = view(allParameters, offsetIndices)
    offsetVector = modelParameters.offsetVector
    offsetVector[1:modelParameters.numOffset] .= offset

    varianceIndices = modelParameters.varianceIndices
    variance = view(allParameters, varianceIndices)
    varianceVector = modelParameters.varianceVector
    varianceVector[1:modelParameters.numVariance] .= variance
end