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


function createCube(nSamples, lowerBounds, upperBounds, fileSave, fCost::Function; seed=123)

    Random.seed!(seed)
    println("Building cube")

    if isfile(fileSave)
        println("Cube already exists")
        return 
    end

    plan = randomLHC(nSamples*10, nDims)
    bounds = [(lowerBounds[i], upperBounds[i]) for i in eachindex(lowerBounds)] 
    scaledPlan = Matrix(scaleLHC(plan, bounds))
    paramSave = zeros(Float64, (nSamples, nDims))
    k = 1
    for i in 1:nSamples*10
        local cost
        paramI = scaledPlan[i, :]
        try 
            cost = fCost(paramI)
        catch
            cost = Inf
        end
        
        if !(isinf(cost) || isnan(cost))
            paramSave[k, :] .= paramI
            k += 1
        end

        if k + 1 == nSamples
            break
        end
    end

    if k + 1 != nSamples
        println("Error : Did not find $nSamples start guesses for the estimation")
    end

    CSV.write(fileSave, DataFrame(paramSave, :auto))
end
