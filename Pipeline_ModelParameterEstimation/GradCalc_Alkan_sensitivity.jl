using ModelingToolkit, DifferentialEquations, BenchmarkTools, DataFrames, CSV, LSODA, Sundials, ODEInterface, ODEInterfaceDiffEq
using JuMP, NLopt, LinearAlgebra, DiffEqSensitivity, DiffEqFlux
using ModelingToolkit: varmap_to_vars

println("Done loading modules")

include(pwd() * "/Additional_functions/additional_tools.jl")
include(pwd() * "/Additional_functions/Solver_info.jl")


mutable struct StaticParameters
    scale::Vector{Float64}
    shift::Vector{Float64}
    variance::Vector{Float64}
end

struct ConstantValues
    observedObservable::Vector{Int64}
    measurementFor::Vector{Vector{Float64}}
    minVariance::Float64
    observedAtIndex::Vector{Vector{Int64}}
    optParameterIndices::Vector{Int64}
    ts::Vector{Float64}
    variablesObservedBy::Vector{Vector{Int64}}
    #
    numObservables::Int64
    numUsedParameters::Int64
    numTimePointsFor::Vector{Int64}
    numTimeSteps::Int64
    numVariables::Int64
    numVariablesObservedBy::Vector{Int64}
end

struct MutatingArrays
    dfdu::Vector{Float64}
    dhdp::Matrix{Float64}
    dp::Vector{Matrix{Float64}}
    h_bar::Vector{Vector{Float64}}
    h_hat::Vector{Vector{Float64}}
    costFor::Vector{Float64}
    p_vector::Vector{Float64}
    u::Matrix{Float64}
end

struct FilesAndPaths
    modelFile::String
    modelPath::String
    writefile_solver::String
    writefile_sensealg::String
end


function f_proto(prob, staticParameters, constantValues, mutatingArrays, p_tuple...)
    println("in f()")

    p_vector = mutatingArrays.p_vector
    p_vector[constantValues.optParameterIndices] .= p_tuple
    _prob = remake(prob, p = p_vector) 
    
    sol = OrdinaryDiffEq.solve(_prob, AutoVern7(Rodas5()), reltol=1e-8, abstol=1e-8, saveat = constantValues.ts)

    h_bar = mutatingArrays.h_bar
    observedAtIndex = constantValues.observedAtIndex

    h_bar[1] = (sol[Symbol("Cells_Apo4(t)")] + sol[Symbol("Cells_Apo_ReOx(t)")])[observedAtIndex[1]]
    h_bar[2] = (sol[Symbol("Cells(t)")] + sol[Symbol("Cells_Cycle_S(t)")] + sol[Symbol("Cells_Cycle_G2(t)")] + 
                sol[Symbol("Cells_SSBDamage_S(t)")] + sol[Symbol("Cells_DSBDamage_S(t)")] + sol[Symbol("Cells_DSBDamage_G2(t)")] + 
                sol[Symbol("Cells_Apo(t)")] + sol[Symbol("Cells_Apo1(t)")] + sol[Symbol("Cells_Apo2(t)")] + 
                sol[Symbol("Cells_Apo3(t)")] + sol[Symbol("Cells_Apo4(t)")] +  + sol[Symbol("Cells_Apo_ReOx(t)")])[observedAtIndex[2]]
    h_bar[3] = (sol[Symbol("Wip1_mRNA_G2(t)")] + sol[Symbol("Wip1_mRNA_S(t)")])[observedAtIndex[3]]
    h_bar[4] = (sol[Symbol("p21_mRNA_S(t)")] + sol[Symbol("p21_mRNA_G2(t)")])[observedAtIndex[4]]
    h_bar[5] = (sol[Symbol("pATM_G2(t)")] + sol[Symbol("pATM_S(t)")])[observedAtIndex[5]]
    h_bar[6] = (sol[Symbol("pChk1_G2(t)")] + sol[Symbol("pChk1_S(t)")])[observedAtIndex[6]]
    h_bar[7] = (sol[Symbol("pChk2_G2(t)")] + sol[Symbol("pChk2_S(t)")])[observedAtIndex[7]]
    h_bar[8] = (sol[Symbol("pDNAPK_G2(t)")] + sol[Symbol("pDNAPK_S(t)")])[observedAtIndex[8]]
    h_bar[9] = (sol[Symbol("pp53_S(t)")] + sol[Symbol("pp53_G2(t)")])[observedAtIndex[9]]
    h_bar[10] = (sol[Symbol("p21_S(t)")] + sol[Symbol("p21_G2(t)")])[observedAtIndex[10]]
    h_bar[11] = (sol[Symbol("pp53_S(t)")] + sol[Symbol("pp53_G2(t)")])[observedAtIndex[11]]
    h_bar[12] = (sol[Symbol("yH2AX_S(t)")] + sol[Symbol("yH2AX_G2(t)")])[observedAtIndex[12]]

    h_hat = mutatingArrays.h_hat
    scale = staticParameters.scale
    shift = staticParameters.shift
    variance = staticParameters.variance
    costFor = mutatingArrays.costFor
    measurementFor = constantValues.measurementFor
    num_t = constantValues.numTimePointsFor
    numObservables = constantValues.numObservables
    minVariance = constantValues.minVariance

    for i = 3:numObservables
        scale[i] = (dot(h_bar[i], measurementFor[i]) - shift[i] * sum(h_bar[i])) / dot(h_bar[i], h_bar[i])
    end
    replace!(scale, NaN=>0.0) # necessary when h_bar is zeros 
    # Mayby replace 0.0 (for non zero h_bar) with a lower bound

    for i in constantValues.observedObservable
        h_hat[i] = scale[i] * h_bar[i] .+ shift[i]
        variance[i] = (dot(measurementFor[i], measurementFor[i]) - 2 * dot(measurementFor[i], h_hat[i]) + dot(h_hat[i], h_hat[i])) / num_t[i]
        variance[i] = max(variance[i], minVariance)
        costFor[i] = (num_t[i]/2) * log(2*pi*variance[i]) + (dot(measurementFor[i], measurementFor[i]) - 2*dot(measurementFor[i], h_hat[i]) + dot(h_hat[i], h_hat[i])) / (2 * variance[i])
    end

    println(sum(costFor[constantValues.observedObservable]))
    return sum(costFor[constantValues.observedObservable])
end

function f_prime_proto(grad, senProb, staticParameters, constantValues, mutatingArrays, p_tuple...)
    println("In f_prime()")

    p_vector = mutatingArrays.p_vector
    p_vector[constantValues.optParameterIndices] .= p_tuple

    _senProb = remake(senProb, p = p_vector) 

    senSol = OrdinaryDiffEq.solve(_senProb, DP8(), reltol=1e-8, abstol=1e-8, saveat = constantValues.ts)
    u = mutatingArrays.u
    dp = mutatingArrays.dp
    u[:], dp[:] = extract_local_sensitivities(senSol)
    # dp[i] containst the sensitivities for all states with regards to parameter i

    dhdp = mutatingArrays.dhdp
    variablesObservedBy = constantValues.variablesObservedBy
    numVariablesObservedBy = constantValues.numVariablesObservedBy
    observedAtIndex = constantValues.observedAtIndex
    numUsedParameters = constantValues.numUsedParameters
    numTimePointsFor = constantValues.numTimePointsFor
    for i in 1:numUsedParameters
        cumulativeLen = 0
        for j in constantValues.observedObservable
            lenj = numTimePointsFor[j]
            dhdp[i, cumulativeLen+1:cumulativeLen+lenj] = sum([@view dp[i][variablesObservedBy[j][k], observedAtIndex[j]] for k in 1:numVariablesObservedBy[j]])
            cumulativeLen += lenj
        end
    end
        
    dfdu = mutatingArrays.dfdu
    scale = staticParameters.scale
    variance = staticParameters.variance
    h_hat = mutatingArrays.h_hat
    measurementFor = constantValues.measurementFor
    dfdu[:] = reduce(vcat, (scale ./ variance) .* (h_hat .- measurementFor))
    grad[:] = (dhdp * dfdu)[constantValues.optParameterIndices]
    println(grad)   

    nothing
end



function GradCalc_sensitivity_eq(filesAndPaths, timeEnd, relevantMeasurementData, observables, inputParameterValues)
    minVariance = 1e-2

    observableIDs = relevantMeasurementData[:,1]
    ts = sort(unique(relevantMeasurementData[:,5]))
    numTimeSteps = length(ts)
    numObservables = length(observables)

    isObserved = [observables[i] in unique(observableIDs) for i=1:numObservables]
    observedObservable = collect(1:numObservables)[isObserved]
    
    observedAt = Array{Array{Float64, 1}, 1}(undef, numObservables)
    for (i, obsId) = enumerate(observables)
        observedAt[i] = relevantMeasurementData[observableIDs .== obsId, 5]
    end
    observedAtIndex = Array{Array{Int64, 1}, 1}(undef, numObservables)
    for (i, obsId) = enumerate(observables)
        observedAtIndex[i] = indexin(observedAt[i], ts)
    end
    measurementFor = Array{Array{Float64, 1}, 1}(undef, numObservables)
    for (i, obsId) = enumerate(observables)
        measurementFor[i] = relevantMeasurementData[observableIDs .== obsId, 4]
    end

    include(filesAndPaths.modelPath)
    new_sys = ode_order_lowering(sys)
    u0 = initialSpeciesValues
    p = trueParameterValues 

    c = trueConstantsValues
    pars = [p;c]
    tspan = (0.0, timeEnd)
    prob = ODEProblem(new_sys,u0,tspan,pars,jac=true)

    numUsedParameters = length(prob.p)
    problemParameterSymbols = parameters(new_sys)
    inputParameterSymbols = [:Dox_level, :Gem_level, :SN38_level]
    inputParameterIndices = collect(1:numUsedParameters)[[pPS.name ∈ inputParameterSymbols for pPS in problemParameterSymbols]]
    numInputParameters = length(inputParameterIndices)
    optParameterIndices = collect(1:numUsedParameters)[[pPS.name ∉ inputParameterSymbols for pPS in problemParameterSymbols]]
    numOptParameters = length(optParameterIndices)

    orderedParameters = varmap_to_vars(pars, parameters(new_sys))[optParameterIndices]
    #orderedParameters = varmap_to_vars(pars, parameters(new_sys))
    orderedParameters .+= (rand() - 0.5)*2
    test_orderedParameters = [0.5963242072135031, 10.493146364550325, 0.5272709883024861, 0.6341879767569332, 0.5419473609958273, 0.8495974979765691, 0.4931463645503251, 0.5833646215561278, 0.4985240061626153, 0.49315636455032524, 0.49315636489693954, 0.49450081003042373, 3.8941626652463848, 1000.4931463644693, 0.8495975440479882, 1.493146364550249, 1.493146364550246, 0.49315636455032513, 0.5931463645503251, 0.6166479339949781, 0.4931463645503251, 50.363193375832125, 174.92788346257532, 1000.4928206967563, 1000.4931463644702, 0.5236459913737985, 0.49315636455032513, 0.5443854987023068, 0.6658139068659021, 0.49315636455032513, 0.4931563645503252, 1.4931463645503251, 0.5000133675306864, 0.9931463645503251, 0.8931463645503251, 0.8173439927147161, 215.2534279479523, 0.49315636455032513, 0.6805384575877701, 0.6025007692235561, 3.904789288616315, 0.49315636455032513, 1000.4931463645503, 0.49315636455032513, 7.386228342688855, 10.493146364549526, 0.5723297840375486, 100.49314636455033, 10.493146353401634, 8.925645615519846, 7.081215591635155, 10.493144013792245]
    #orderedParameters[inputParameterIndices] .= 0.0
    orderedParameters = test_orderedParameters[optParameterIndices]

    senProb = ODELocalSensitivityProblem(new_sys, prob.u0, tspan, orderedParameters)

    numVariables = length(u0)
    numTimePointsFor = Vector{Int64}(undef, numObservables)
    numTimePointsFor .= length.(observedAt)

    variablesObservedBy = [[13, 18], [2, 8, 11, 12, 13, 18, 19, 20, 25, 29, 30, 36], [0], [0], [6, 10], [31, 33], [21, 28], [32, 34], [9, 35], [7, 24], [9, 35], [4, 14]]
    numVariablesObservedBy = [2, 12, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2]
    
    constantValues = ConstantValues(observedObservable, measurementFor, minVariance, observedAtIndex, optParameterIndices,
                                    ts, variablesObservedBy, numObservables, numUsedParameters, numTimePointsFor, numTimeSteps, numVariables, numVariablesObservedBy)

    scale = Array{Float64, 1}(undef, numObservables)
    scale[[1,2]] .= 1.0
    shift = zeros(Float64, numObservables)
    variance = Array{Float64, 1}(undef, numObservables)
    shift[[3,4]] .= 1.0
    staticParameters = StaticParameters(scale, shift, variance)


    dfdu = Vector{Float64}(undef, sum(numTimePointsFor))
    dhdp = Array{Float64, 2}(undef, numUsedParameters, sum(numTimePointsFor))
    dp = Vector{Array{Float64, 2}}(undef, numUsedParameters)
    for i = 1:numUsedParameters
        dp[i] = Array{Float64, 2}(undef, numVariables, numTimeSteps)
    end
    h_bar = Array{Array{Float64, 1}, 1}(undef, numObservables)
    h_hat = Array{Array{Float64, 1}, 1}(undef, numObservables)
    h_hat[3] = Float64[]
    h_hat[4] = Float64[]
    costFor = Array{Float64, 1}(undef, numObservables)
    p_vector = Vector{Float64}(undef, numUsedParameters)
    p_vector[inputParameterIndices] = inputParameterValues
    u = Array{Float64, 2}(undef, numVariables, numTimeSteps)
    mutatingArrays = MutatingArrays(dfdu, dhdp, dp, h_bar, h_hat, costFor, p_vector, u)

    f = (p_tuple...) -> f_proto(prob, staticParameters, constantValues, mutatingArrays, p_tuple...)   

    f_prime = (grad, p_tuple...) -> f_prime_proto(grad, senProb, staticParameters, constantValues, mutatingArrays, p_tuple...)


    # test first cost and grad at test_orderedParameters
    cost = f(orderedParameters...)
    println(cost)
    # 29.216144533943794
    grad = zeros(Float64, numOptParameters)
    f_prime(grad, orderedParameters...)
    println(grad)
    # [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.9510894754496647e-16, 0.0, 10.276324511313856, 12.216211010949758, 0.0, 0.0, 0.0, -1.7268644853198856, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.2758121304483095e-5, 0.0, 0.0, 0.0, 0.0]
    
    
    #=
    model = Model(NLopt.Optimizer)
    #JuMP.register(model::Model, s::Symbol, dimension::Integer, f::Function, ∇f::Function, ∇²f::Function)
    register(model, :f, numOptParameters, f, f_prime)
    set_optimizer_attribute(model, "algorithm", :LD_MMA)
    @variable(model, p[1:numOptParameters] >= 0) # fix lower and upper bounds
    for i in 1:numOptParameters
        set_start_value(p[i], orderedParameters[i])
    end
    @NLobjective(model, Min, f(p...))
    JuMP.optimize!(model)
    #@show termination_status(model)
    #@show primal_status(model)
    println("Done!")
    p_opt = [value(p[i]) for i=1:numOptParameters]
    println(p_opt)
    cost_opt = objective_value(model)
    println(cost_opt)
    =#
    
end



function main()
    
    readModelPath = pwd() * "/Pipeline_SBMLImporter/JuliaModels"
    modelFile = "model_Alkan_SciSignal2018.jl"
    modelPath = readModelPath * "/" * modelFile
    writePath = pwd() * "/Pipeline_ModelParameterEstimation/IntermediaryResults"
    fixDirectories(writePath)
    writefile_solver = writePath * "/benchmark_solver_" * string(getNumberOfFiles(writePath) + 1) * ".csv"
    writefile_sensealg = writePath * "/benchmark_sensealg_" * string(getNumberOfFiles(writePath) + 1) * ".csv"

    filesAndPaths = FilesAndPaths(modelFile, modelPath, writefile_solver, writefile_sensealg)
    
    readTimePath = pwd() * "/Pipeline_ModelSolver/IntermediaryResults"
    timeEnds = CSV.read(readTimePath * "/timeScales.csv", DataFrame)
    timeEnd = timeEnds[timeEnds[:,1] .== modelFile, 2][1]
    readDataPath = pwd() * "/Pipeline_ModelParameterEstimation/Data" * "/" * modelFile[1:end-3]
    experimentalConditions = CSV.read(readDataPath * "/experimentalCondition_Alkan_SciSignal2018.tsv", DataFrame)
    measurementData = CSV.read(readDataPath * "/measurementData_Alkan_SciSignal2018.tsv", DataFrame)
    relevantMeasurementData = measurementData[measurementData[:,3] .== experimentalConditions[1,1], :]
    observables = CSV.read(readDataPath * "/observables_Alkan_SciSignal2018.tsv", DataFrame)[:,1]

    inputParameterValues = [0.0, 0.0, 0.0]

    println("Starting optimisation")
    GradCalc_sensitivity_eq(filesAndPaths, timeEnd, relevantMeasurementData, observables, inputParameterValues)

end

main()