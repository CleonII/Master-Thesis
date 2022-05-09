# using LatinHypercubeSampling, CSV, DataFrames, Random

function createSamples(nDims)
    Random.seed!(123)
    dirName = joinpath(pwd(), "Pipeline_ModelParameterEstimation", "LatinHyperCubeParameters", )
    fileName = string(nDims) * ".csv"
    plan = LHCoptim(50, nDims, 1000)[1]
    CSV.write(joinpath(dirName, fileName), DataFrame(plan, :auto))
end

function getSamples(nDims, lowerBounds, upperBounds, index)
    Random.seed!(123)
    dirName = joinpath(pwd(), "Pipeline_ModelParameterEstimation", "LatinHyperCubeParameters")
    fileName = string(nDims) * ".csv"
    if fileName in readdir(dirName)
        plan = Matrix(CSV.read(joinpath(dirName, fileName), DataFrame))
    else
        plan = LHCoptim(50, nDims, 1000)[1]
        CSV.write(joinpath(dirName, fileName), DataFrame(plan, :auto))
    end
    bounds = Array{Tuple{Float64, Float64}}(undef, nDims)
    for i in eachindex(bounds)
        bounds[i] = (lowerBounds[i], upperBounds[i])
    end

    scaled_plan = vec(Matrix(scaleLHC(plan, bounds))[index, :])
    return scaled_plan
end
