using ModelingToolkit, DifferentialEquations, BenchmarkTools, DataFrames, CSV
using JuMP, NLopt, LinearAlgebra, DiffEqSensitivity, ForwardDiff
using ModelingToolkit: varmap_to_vars
using ForwardDiff: GradientConfig, Chunk

println("Done loading modules")

include(joinpath(pwd(), "Additional_functions", "additional_tools.jl"))
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "LatinHyperCubeParameters", "LatinHyperCubeSampledParameters.jl"))
include(joinpath(pwd(), "Additional_functions", "importModelInfo.jl"))

function optModelSaveResults(model, p, doLogSearch, option1, option2, option3, iStartPar, iterations, methodFile, write)

    benchRunTime = Vector{Float64}(undef, iterations)
    benchMemory = Vector{Float64}(undef, iterations)
    benchAllocs = Vector{Float64}(undef, iterations)
    success = true
    p_opt = Vector{Float64}(undef, length(p))
    cost_opt = 0.0 

    p_save = copy(p)

    try
        JuMP.optimize!(model)
        p_opt[:] = [value(p[i]) for i=1:length(p)]
        view(p_opt, doLogSearch) .= exp10.(view(p_opt, doLogSearch))
        cost_opt = objective_value(model)
        success = true 
    catch err
        println(err)
        success = false
        p_opt .= NaN
        cost_opt = Inf
    end

    terminationStatus = termination_status(model)
    primalStatus = primal_status(model)

    if success
        for iter in 1:iterations
            b = @benchmark begin
                $p = copy($p_save)
                JuMP.optimize!($model) 
            end samples=1 evals=1 
            
            bMin = minimum(b)
            benchRunTime[iter] = bMin.time # microsecond
            benchMemory[iter] = bMin.memory # bytes
            benchAllocs[iter] = bMin.allocs # number of allocations
        end
    else
        benchRunTime .= NaN
        benchMemory .= NaN
        benchAllocs .= NaN
    end

    data = DataFrame(method = methodFile[1:end-3], option1 = option1, option2 = option2, option3 = option3, startParameterIndex = iStartPar, 
            success = success, terminationStatus = terminationStatus, primalStatus = primalStatus, cost = cost_opt, 
            runTime = benchRunTime, memory = benchMemory, allocs = benchAllocs, iteration = 1:iterations)
    if isfile(write)
        CSV.write(write, data, append = true)
    else
        CSV.write(write, data)
    end

    data_pars = DataFrame(p_opt', :auto)
    insertcols!(data_pars, 1, :method => methodFile[1:end-3])
    insertcols!(data_pars, 2, :option1 => option1)
    insertcols!(data_pars, 3, :option2 => option2)
    insertcols!(data_pars, 4, :option3 => option3)
    insertcols!(data_pars, 5, :startParameterIndex => iStartPar)
    insertcols!(data_pars, 6, :iteration => 1)
    write_pars = write[1:end-4] * "_pars.csv"
    if isfile(write_pars)
        CSV.write(write_pars, data_pars, append = true)
    else
        CSV.write(write_pars, data_pars)
    end
            
end


function optAdamSaveResults(step, adam_opt, doLogSearch, n_it, b2, stepRange, iStartPar, iterations, methodFile, write)
    loss_val = zeros(Float64, n_it)
    theta_save = copy(adam_opt.theta)

    benchRunTime = Vector{Float64}(undef, iterations)
    benchMemory = Vector{Float64}(undef, iterations)
    benchAllocs = Vector{Float64}(undef, iterations)
    success = Vector{Bool}(undef, iterations)
    success .= true
    p_opt = Vector{Float64}(undef, length(theta_save))
    cost_opt = Vector{Float64}(undef, iterations)
    
    for iter in 1:iterations 
        
        b = @benchmark begin
            $adam_opt.theta = copy($theta_save)
            $adam_opt.loss = 0.0
            $adam_opt.m = zeros(length($theta_save))
            $adam_opt.v = zeros(length($theta_save))
            $adam_opt.t = 0
            $adam_opt.β = 1.0
            $adam_opt.fail = 0
            for i in 1:$n_it
                $step()
                $loss_val[i] = $adam_opt.loss
    
                if $adam_opt.fail >= 10
                    $success[$iter] = false
                    break
                end
            end
        end samples=1 evals=1

        cost_opt[iter] = adam_opt.loss
        p_opt[:] = adam_opt.theta
        view(p_opt, doLogSearch) .= exp10.(view(p_opt, doLogSearch))

        bMin = minimum(b)
        benchRunTime[iter] = bMin.time # microsecond
        benchMemory[iter] = bMin.memory # bytes
        benchAllocs[iter] = bMin.allocs # number of allocations

        data_pars = DataFrame(p_opt', :auto)
        insertcols!(data_pars, 1, :method => methodFile[1:end-3])
        insertcols!(data_pars, 2, :option1 => n_it)
        insertcols!(data_pars, 3, :option2 => b2)
        insertcols!(data_pars, 4, :option3 => string(stepRange))
        insertcols!(data_pars, 5, :startParameterIndex => iStartPar)
        insertcols!(data_pars, 6, :iteration => iter)
        write_pars = write[1:end-4] * "_pars.csv"
        if isfile(write_pars)
            CSV.write(write_pars, data_pars, append = true)
        else
            CSV.write(write_pars, data_pars)
        end

        data_cost = DataFrame(loss_val', :auto)
        insertcols!(data_cost, 1, :method => methodFile[1:end-3])
        insertcols!(data_cost, 2, :option1 => n_it)
        insertcols!(data_cost, 3, :option2 => b2)
        insertcols!(data_cost, 4, :option3 => string(stepRange))
        insertcols!(data_cost, 5, :startParameterIndex => iStartPar)
        insertcols!(data_cost, 6, :iteration => iter)
        write_cost = write[1:end-4] * "_forwGrad_cost.csv"
        if isfile(write_cost)
            CSV.write(write_cost, data_cost, append = true)
        else
            CSV.write(write_cost, data_cost)
        end

    end

    data = DataFrame(method = methodFile[1:end-3], option1 = n_it, option2 = b2, option3 = string(stepRange), startParameterIndex = iStartPar, 
            success = success, terminationStatus = "-", primalStatus = "-", cost = cost_opt,
            runTime = benchRunTime, memory = benchMemory, allocs = benchAllocs, iteration = 1:iterations)
    if isfile(write)
        CSV.write(write, data, append = true)
    else
        CSV.write(write, data)
    end

end


function benchMethod(filesAndPaths, methodFile, iStartPar, iterations, adjointSensitivities, forwardAutomaticDifferentiation, forwardGradient)
    write = joinpath(filesAndPaths.writePath, filesAndPaths.writeFile)
    
    # runtime, objective_value (at each step), end parameters, iterations 
    senseAlgs = getSenseAlgs()
    optAlgs = [:LD_MMA, :LD_SLSQP, :LD_LBFGS]
    chunkSizes = [23, 60, 115]
    stepRanges = [[0.0, -3.0], [-1.0, -4.0], [-2.0, -5.0], [-3.0, -3.0]]
    b2s = [0.999, 0.9]
    n_its = [100, 500, 1000, 5000, 10000]
    # r?
    # c?

    if methodFile == "adjointSensitivities.jl"
        println("Running: adjointSensitivities.jl \n")

        for senseAlg in senseAlgs
            model, p, doLogSearch = adjointSensitivities(iStartPar, senseAlg, optAlgs[1])
            optModelSaveResults(model, p, doLogSearch, senseAlg, optAlgs[1], "-", iStartPar, iterations, methodFile, write)
        end
        for optAlg in optAlgs
            model, p, doLogSearch = adjointSensitivities(iStartPar, senseAlgs[1], optAlg)
            optModelSaveResults(model, p, doLogSearch, senseAlgs[1], optAlg, "-", iStartPar, iterations, methodFile, write)
        end
        
    elseif methodFile == "forwardAutomaticDifferentiation.jl"
        println("Running: forwardAutomaticDifferentiation.jl \n")

        for chunkSize in chunkSizes
            model, p, doLogSearch = forwardAutomaticDifferentiation(iStartPar, chunkSize, optAlgs[1])
            optModelSaveResults(model, p, doLogSearch, chunkSize, optAlgs[1], "-", iStartPar, iterations, methodFile, write)
        end
        for optAlg in optAlgs
            model, p, doLogSearch = forwardAutomaticDifferentiation(iStartPar, chunkSizes[1], optAlg)
            optModelSaveResults(model, p, doLogSearch, chunkSizes[1], optAlg, "-", iStartPar, iterations, methodFile, write)
        end
        
    elseif methodFile == "forwardGradient.jl"
        println("Running: forwardGradient.jl \n")

        for n_it in n_its
            step, adam_opt, doLogSearch = forwardGradient(iStartPar, n_it, b2s[1], stepRanges[2])
            optAdamSaveResults(step, adam_opt, doLogSearch, n_it, b2s[1], stepRanges[2], iStartPar, iterations, methodFile, write)
        end

        for b2 in b2s
            step, adam_opt, doLogSearch = forwardGradient(iStartPar, n_its[3], b2, stepRanges[2])
            optAdamSaveResults(step, adam_opt, doLogSearch, n_its[3], b2, stepRanges[2], iStartPar, iterations, methodFile, write)
        end

        for stepRange in stepRanges
            step, adam_opt, doLogSearch = forwardGradient(iStartPar, n_its[3], b2s[1], stepRange)
            optAdamSaveResults(step, adam_opt, doLogSearch, n_its[3], b2s[1], stepRange, iStartPar, iterations, methodFile, write)
        end

    end

end





function iterateMethods(filesAndPaths, methodFiles, iterations, adjointSensitivities, forwardAutomaticDifferentiation, forwardGradient)
    
    for iStartPar=1:50
        for methodFile in methodFiles
            benchMethod(filesAndPaths, methodFile, iStartPar, iterations, adjointSensitivities, forwardAutomaticDifferentiation, forwardGradient)
        end
    end
end



function main(; modelName = "model_Bachmann_MSB2011")
    
    modelPath = joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels")
    modelNameShort = modelName[7:end]
    modelFile = modelName * ".jl" 
    methodPath = joinpath(pwd(), "Pipeline_ModelParameterEstimation", modelName)
    writePath = joinpath(pwd(), "Pipeline_ModelParameterEstimation", "IntermediaryResults")
    writeFile = joinpath("benchmark_" * modelName * "_" * string(getNumberOfFiles(writePath) + 1) * ".csv")

    filesAndPaths = FilesAndPaths(modelName, modelPath, modelFile, methodPath, writePath, writeFile)

    fixDirectories(writePath)

    # Reading data for model
    readDataPath = joinpath(pwd(), "Pipeline_ModelParameterEstimation", "Data", modelName)
    dataEnding = modelNameShort * ".tsv"
    experimentalConditions = CSV.read(joinpath(readDataPath, "experimentalCondition_" * dataEnding), DataFrame)
    measurementData = CSV.read(joinpath(readDataPath, "measurementData_" * dataEnding), DataFrame)
    timeEnd = maximum(measurementData[!, 5])
    observables = CSV.read(joinpath(readDataPath, "observables_" * dataEnding), DataFrame)
    parameterBounds = CSV.read(joinpath(readDataPath, "parameters_" * dataEnding), DataFrame)
    #

    methodFiles = readdir(methodPath)
    for methodFile in methodFiles
        if methodFile !== "sensitivityEquations.jl"
            include(joinpath(methodPath, methodFile))
        end
    end

    adjointSensitivities = (iStartPar, senseAlg, optAlg) -> adjointSensitivities_proto(iStartPar, senseAlg, optAlg, 
            filesAndPaths, timeEnd, experimentalConditions, measurementData, observables, parameterBounds)

    forwardAutomaticDifferentiation = (iStartPar, chunkSize, optAlg) -> forwardAutomaticDifferentiation_proto(iStartPar, chunkSize, optAlg, 
            filesAndPaths, timeEnd, experimentalConditions, measurementData, observables, parameterBounds)

    forwardGradient = (iStartPar, n_it, b2, stepRange) -> forwardGradient_proto(iStartPar, n_it, b2, stepRange, 
            filesAndPaths, timeEnd, experimentalConditions, measurementData, observables, parameterBounds)

    println("Starting benchmark")
    
    iterations = 5
    iterateMethods(filesAndPaths, methodFiles, iterations, adjointSensitivities, forwardAutomaticDifferentiation, forwardGradient)

end

main()