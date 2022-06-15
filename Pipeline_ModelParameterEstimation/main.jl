using ModelingToolkit, DifferentialEquations, DataFrames, CSV, Random
using JuMP, NLopt, Ipopt, LinearAlgebra, DiffEqSensitivity, ForwardDiff
using ModelingToolkit: varmap_to_vars
using ForwardDiff: GradientConfig, Chunk
using LatinHypercubeSampling
using DelimitedFiles
using Optim 


println("Done loading modules")

include(joinpath(pwd(), "Additional_functions", "additional_tools.jl"))
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "LatinHyperCubeParameters", "LatinHyperCubeSampledParameters.jl"))
include(joinpath(pwd(), "Additional_functions", "importModelInfo.jl"))
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "optAndSave.jl"))

include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "CommonParameterEstimationMethods", "adjointSensitivities.jl"))
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "CommonParameterEstimationMethods", "forwardAutomaticDifferentiation.jl"))
include(joinpath(pwd(), "Pipeline_ModelParameterEstimation", "CommonParameterEstimationMethods", "forwardGradient.jl"))


allModelFunctionVector = includeAllModels(getModelFiles(joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels")), 
        joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels"))

modelsWithWorkingMethods = ["model_Alkan_SciSignal2018", "model_Bachmann_MSB2011", "model_Boehm_JProteomeRes2014"]
allMethodFunctionArray = includeAllMethods(joinpath.(pwd(), "Pipeline_ModelParameterEstimation", modelsWithWorkingMethods))

# Include with hessian 
include(pwd() * "/Pipeline_ModelParameterEstimation/With_Full_Hessian/model_Boehm_JProteomeRes2014/forwardAutomaticDifferentiation_with_hessian.jl")

function benchMethod_proto(usedModelFunction, usedMethodFunctions, iStartPar, filesAndPaths, timeEnd, experimentalConditions, measurementData, observables, parameterBounds,
        optAlg, method)
    write = joinpath(filesAndPaths.writePath, filesAndPaths.writeFile)
    
    solver = QNDF()
    
    if method == "adjointSensitivities" || method == "all"
        println("Running: adjointSensitivities \n")
        senseAlgs = getSenseAlgs()
        solver_adjoint = KenCarp4()

        for senseAlg in senseAlgs
            model, p, doLogSearch = usedMethodFunctions[1](usedModelFunction, iStartPar, senseAlg, optAlg, solver, solver_adjoint, 
                    timeEnd, experimentalConditions, measurementData, observables, parameterBounds)
            optModelSaveResults(model, p, doLogSearch, solver, optAlg, senseAlg, solver_adjoint, iStartPar, "adjointSensitivities", write)
        end
    end
    if method == "forwardAutomaticDifferentiation" || method == "all"
        println("Running: forwardAutomaticDifferentiation \n")

        model, p, doLogSearch = usedMethodFunctions[2](usedModelFunction, iStartPar, optAlg, solver, 
                timeEnd, experimentalConditions, measurementData, observables, parameterBounds)
        optModelSaveResults(model, p, doLogSearch, solver, optAlg, "-", "-", iStartPar, "forwardAutomaticDifferentiation", write)

    end
    if method == "forwardGradient" || method == "all"
        println("Running: forwardGradient \n")
        stepRanges = [[0.0, -3.0], [-1.0, -4.0], [-2.0, -5.0], [-3.0, -3.0]]
        b2s = [0.999, 0.9]
        n_its = [100, 500, 1000, 5000, 10000]
        iterations = 4

        for n_it in n_its
            for b2 in b2s
                for stepRange in stepRanges

                    step, adam_opt, doLogSearch = usedMethodFunctions[3](usedModelFunction, iStartPar, n_it, b2, stepRange, solver, 
                            timeEnd, experimentalConditions, measurementData, observables, parameterBounds)
                    optAdamSaveResults(step, adam_opt, doLogSearch, solver, n_it, b2, stepRange, iStartPar, iterations, "forwardGradient", write)

                end
            end
        end
    end

    nothing

end





function iterateMethods(benchMethod, numStartGuesses)
    
    for iStartPar=1:numStartGuesses
        benchMethod(iStartPar)
    end

    nothing

end



function main(; modelName = "model_Bachmann_MSB2011", optAlg = :Ipopt, method = "all")
    
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
    
    allModelFiles = getModelFiles(modelPath)
    usedModelFunction = allModelFunctionVector[[allModelFile in [modelFile] for allModelFile in allModelFiles]][1]

    usedMethodFunctions = allMethodFunctionArray[findfirst(modelsWithWorkingMethods .== modelName), :]

    benchMethod = (iStartPar) -> benchMethod_proto(usedModelFunction, usedMethodFunctions, iStartPar, filesAndPaths, timeEnd, experimentalConditions, 
            measurementData, observables, parameterBounds, optAlg, method)

    println("Starting benchmark")
    
    numStartGuesses = 10
    iterateMethods(benchMethod, numStartGuesses)

    nothing

end

#main(modelName = "model_Boehm_JProteomeRes2014", method = "forwardAutomaticDifferentiation")
#optAlgs = [:LD_MMA, :LD_LBFGS, :Ipopt]
#methods = ["adjointSensitivities", "forwardAutomaticDifferentiation", "forwardGradient", "all"]


function test_hessian(;get_functions=false, i=1)

    modelName = "model_Boehm_JProteomeRes2014"
    solver = QNDF()
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
        
    allModelFiles = getModelFiles(modelPath)
    usedModelFunction = allModelFunctionVector[[allModelFile in [modelFile] for allModelFile in allModelFiles]][1]
    usedMethodFunctions = allMethodFunctionArray[findfirst(modelsWithWorkingMethods .== modelName), :]

    if get_functions == true
        return forwardAutomaticDifferentiation_hessian_proto_model_Boehm_JProteomeRes2014(usedModelFunction, i, solver, timeEnd, experimentalConditions, measurementData, observables, parameterBounds)
    end

    cost_array = zeros(50)
    for i in 1:50
        println("i = $i")
        f_pre, f_grad_pre, f_hess_pre, p, lb, ub = forwardAutomaticDifferentiation_hessian_proto_model_Boehm_JProteomeRes2014(usedModelFunction, i, solver, timeEnd, experimentalConditions, measurementData, observables, parameterBounds)

        p0 = deepcopy(p)
        f_use = (x) -> f_pre(x...)
        f_grad_use = (grad, x) -> f_grad_pre(grad, x...)
        f_hess_use = (hess, x) -> f_hess_pre(hess, x...)
        df = TwiceDifferentiable(f_use, f_grad_use, f_hess_use, p0)
        dfc = TwiceDifferentiableConstraints(lb .- 0.01, ub .+ 0.01)
        res = Optim.optimize(df, dfc, p0, IPNewton(), Optim.Options(iterations = 1000, show_trace = true))
        cost_array[i] = res.minimum
    end

    return cost_array
end


function eval_h(x_arg::Vector{Float64}, 
                rows::Vector{Int32}, 
                cols::Vector{Int32}, 
                obj_factor::Float64, 
                lambda::Vector{Float64}, 
                values::Union{Nothing,Vector{Float64}}, 
                n_param, 
                calc_hess::Function)

    if values === nothing
        # Symmetric matrix, fill the lower left triangle only
        idx = 1
        for row in 1:n_param
            for col in 1:row
                rows[idx] = row
                cols[idx] = col
                idx += 1
            end
        end
    else
        # Again, only lower left triangle
        # Objective
        hessian_mat = zeros(n_param, n_param)
        calc_hess(hessian_mat, x_arg)
        idx = 1
        for row in 1:n_param
            for col in 1:row
                values[idx] = hessian_mat[row, col] * obj_factor
                idx += 1
            end
        end
    end
    return
end
function eval_jac_g(x::Vector{Float64}, rows::Vector{Int32}, cols::Vector{Int32}, values::Union{Nothing,Vector{Float64}})
    return 
end
function eval_g(x::Vector{Float64}, g::Vector{Float64})
    return 
end


function create_ipopt_prob(i)

    modelName = "model_Boehm_JProteomeRes2014"
    solver = QNDF()
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
        
    allModelFiles = getModelFiles(modelPath)
    usedModelFunction = allModelFunctionVector[[allModelFile in [modelFile] for allModelFile in allModelFiles]][1]
    usedMethodFunctions = allMethodFunctionArray[findfirst(modelsWithWorkingMethods .== modelName), :]

    f_pre, f_grad_pre, f_hess_pre, p, lb, ub = forwardAutomaticDifferentiation_hessian_proto_model_Boehm_JProteomeRes2014(usedModelFunction, i, solver, timeEnd, experimentalConditions, measurementData, observables, parameterBounds)
    p0 = deepcopy(p)
    n_param = length(p0)
    f_use = (x) -> begin cost = f_pre(x...); return cost end
    f_grad_use = (x, grad) -> begin f_grad_pre(grad, x...); return end
    f_hess_use = (hess, x) -> f_hess_pre(hess, x...)
    eval_h_use = (x_arg, rows, cols, obj_factor, lambda, values) -> eval_h(x_arg, rows, cols, obj_factor, lambda, values, n_param, f_hess_use)
    
    m = 0
    n_param_hess = Int(n_param*(n_param + 1) / 2)
    g_L = Float64[]
    g_U = Float64[]
    prob = Ipopt.CreateIpoptProblem(n_param, lb, ub, m, g_L, g_U, 0, n_param_hess, f_use, eval_g, f_grad_use, eval_jac_g, eval_h_use)

    return prob, p0
end


function test_otpimizers()

    modelName = "model_Boehm_JProteomeRes2014"
    solver = QNDF()
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
        
    allModelFiles = getModelFiles(modelPath)
    usedModelFunction = allModelFunctionVector[[allModelFile in [modelFile] for allModelFile in allModelFiles]][1]
    usedMethodFunctions = allMethodFunctionArray[findfirst(modelsWithWorkingMethods .== modelName), :]

    cost_array_newton = zeros(50)
    cost_array_ipopt = zeros(50)
    cost_array_ipopt_bfgs = zeros(50)
    for i in 1:50
        println("i = $i")

        if i == 40
            continue
        end

        f_pre, f_grad_pre, f_hess_pre, p, lb, ub = forwardAutomaticDifferentiation_hessian_proto_model_Boehm_JProteomeRes2014(usedModelFunction, i, solver, timeEnd, experimentalConditions, measurementData, observables, parameterBounds)

        p0 = deepcopy(p)
        f_use = (x) -> f_pre(x...)
        f_grad_use = (grad, x) -> f_grad_pre(grad, x...)
        f_hess_use = (hess, x) -> f_hess_pre(hess, x...)
        df = TwiceDifferentiable(f_use, f_grad_use, f_hess_use, p0)
        dfc = TwiceDifferentiableConstraints(lb .- 0.01, ub .+ 0.01)
        res = Optim.optimize(df, dfc, p0, IPNewton(), Optim.Options(iterations = 1000, show_trace = true))
        cost_array_newton[i] = res.minimum

        ipopt_prob, p0 = create_ipopt_prob(i)
        Ipopt.AddIpoptIntOption(ipopt_prob, "print_level", 5)
        Ipopt.AddIpoptIntOption(ipopt_prob, "max_iter", 1000)
        Ipopt.AddIpoptNumOption(ipopt_prob, "tol", 1e-8)
        Ipopt.AddIpoptNumOption(ipopt_prob, "acceptable_tol", 1e-6)
        ipopt_prob.x = deepcopy(p0)
        sol_opt = Ipopt.IpoptSolve(ipopt_prob)
        cost_array_ipopt[i] = ipopt_prob.obj_val

        #Ipopt.AddIpoptStrOption(ipopt_prob, "hessian_approximation", "limited-memory")
        ipopt_prob.x = deepcopy(p0)
        Ipopt.AddIpoptStrOption(ipopt_prob, "hessian_approximation", "limited-memory")
        sol_opt = Ipopt.IpoptSolve(ipopt_prob)
        cost_array_ipopt_bfgs[i] = ipopt_prob.obj_val
    end

    return cost_array_newton, cost_array_ipopt, cost_array_ipopt_bfgs
end

c_newton, c_ipopt, c_ipopt_bfgs = test_otpimizers()

open(pwd() * "/Array_newton.csv", "w") do io
    writedlm(io, c_newton)
end
open(pwd() * "/Array_ipopt.csv", "w") do io
    writedlm(io, c_ipopt)
end
open(pwd() * "/Array_ipopt_bfgs.csv", "w") do io
    writedlm(io, c_ipopt_bfgs)
end