using LatinHypercubeSampling, CSV, DataFrames, Random
using Distributions

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

    nDims = length(lowerBounds)
    paramSave = zeros(Float64, (nSamples, nDims))
    # The goal, as widely distributed distances as possible 
    useCube = true
    k, maxIter = 1, 1
    while k - 1 != nSamples && maxIter < nSamples*10

        # If model could not be solved for some start-guesses generate a new Cube
        # A new cube is best to maximise spread of start guesses 
        iMax = nSamples - (k - 1)
        println("iMax = ", iMax)
        println("nDims = ", nDims)
        if iMax > 15
            plan = LHCoptim(iMax, nDims, 10)[1]
            bounds = [(lowerBounds[i], upperBounds[i]) for i in eachindex(lowerBounds)] 
            scaledPlan = Matrix(scaleLHC(plan, bounds))
        else
            iMax = 1
            useCube = false
            paramI = [rand(Uniform(lowerBounds[i], upperBounds[i])) for i in eachindex(lowerBounds)]
        end
        

        for i in 1:iMax
            if useCube
                paramI = scaledPlan[i, :]
            end

            local cost
            try 
                cost = fCost(paramI)
            catch
                cost = Inf
            end
            
            if !(isinf(cost) || isnan(cost))
                paramSave[k, :] .= paramI
                k += 1
            end

            if k % 50 == 0
                println("Have found $k start guesses")
            end

            if k - 1 == nSamples
                break
            end

        maxIter += 1
        end
    end

    if k - 1 != nSamples
        println("Error : Did not find $nSamples start guesses for the estimation")
    end

    CSV.write(fileSave, DataFrame(paramSave, :auto))
end
