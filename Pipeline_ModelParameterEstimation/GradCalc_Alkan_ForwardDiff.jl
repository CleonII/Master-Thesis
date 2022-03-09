using ModelingToolkit, DifferentialEquations, BenchmarkTools, DataFrames, CSV, LSODA, Sundials, ODEInterface, ODEInterfaceDiffEq
using JuMP, NLopt, LinearAlgebra, DiffEqSensitivity, DiffEqFlux, ForwardDiff
using ModelingToolkit: varmap_to_vars

println("Done loading modules")

include(pwd() * "/Additional_functions/additional_tools.jl")
include(pwd() * "/Additional_functions/Solver_info.jl")

struct StaticParameters
    scale::Vector{ForwardDiff.Dual}
    shift::Vector{Float64}
    variance::Vector{ForwardDiff.Dual}
end

struct ConstantValues
    observedObservable::Vector{Int64}
    measurementFor::Vector{Vector{Float64}}
    minVariance::Float64
    observablesTimeIndexIndices::Array{Vector{Int64}, 2}
    observedAtIndex::Vector{Vector{Int64}}
    ts::Vector{Float64}
    #
    numObservables::Int64
    numParameters::Int64
    numUsedParameters::Int64
    numTimePointsFor::Vector{Int64}
    numTimeSteps::Int64
    numVariables::Int64
end

struct MutatingArrays
    h_bar::Vector{Vector{ForwardDiff.Dual}}
    h_hat::Vector{Vector{ForwardDiff.Dual}}
    costFor::Vector{ForwardDiff.Dual}
    p_vector::Vector{Float64}
end

struct FilesAndPaths
    modelFile::String
    modelPath::String
    writefile_solver::String
    writefile_sensealg::String
end

function f_proto(prob, staticParameters, constantValues, mutatingArrays, p)
    println("in f()")
    #println(p, "\n")

    #p_vector = mutatingArrays.p_vector
    #p_vector .= p_tuple
    #_prob = remake(prob, p = p_vector) 
    _prob = remake(prob, u0 = convert.(eltype(p), prob.u0), p = p)

    sol = OrdinaryDiffEq.solve(_prob, AutoVern7(Rodas5()), reltol=1e-8, abstol=1e-8, saveat = constantValues.ts)
    
    #println(sol.t)
    #println(sol[Symbol("Cells_Apo4(t)")])
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

    for i = 1:numObservables
        h_hat[i] = scale[i] * h_bar[i] .+ shift[i]
        variance[i] = (dot(measurementFor[i], measurementFor[i]) - 2 * dot(measurementFor[i], h_hat[i]) + dot(h_hat[i], h_hat[i])) / num_t[i]
        variance[i] = max(variance[i], minVariance)
        costFor[i] = (num_t[i]/2) * log(2*pi*variance[i]) + (dot(measurementFor[i], measurementFor[i]) - 2*dot(measurementFor[i], h_hat[i]) + dot(h_hat[i], h_hat[i])) / (2 * variance[i])
    end

    return sum(costFor[constantValues.observedObservable])
end



function GradCalc_forwardDiff(filesAndPaths, timeEnd, relevantMeasurementData, observables)
    minVariance = 1e-2

    observableIDs = relevantMeasurementData[:,1]
    ts = sort(unique(relevantMeasurementData[:,5]))
    numTimeSteps = length(ts)
    numObservables = length(observables)

    isObserved = [observables[i] in unique(observableIDs) for i=1:numObservables]
    observedObservable = collect(1:numObservables)[isObserved]

    observedAt = Vector{Vector{Float64}}(undef, numObservables)
    for (i, obsId) = enumerate(observables)
        observedAt[i] = relevantMeasurementData[observableIDs .== obsId, 5]
    end
    observedAtIndex = Vector{Vector{Int64}}(undef, numObservables)
    for i = 1:numObservables
        observedAtIndex[i] = indexin(observedAt[i], ts)
    end
    measurementFor = Vector{Vector{Float64}}(undef, numObservables)
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
    numParameters = length(pars)
    orderedParameters = varmap_to_vars(pars, parameters(new_sys))
    orderedParameters .+= rand()
    test_orderedParameters = [0.9035043754031366, 10.800326532739959, 0.8344511564921195, 0.9413681449465665, 0.8491275291854607, 1.1567776661662026, 0.8003265327399586, 0.8905447897457612, 0.8057041743522487, 0.8003365327399586, 0.800336533086573, 0.8016809782200571, 4.201342833436018, 1000.800326532659, 1.1567777122376215, 1.8003265327398825, 1.8003265327398794, 0.8003365327399585, 0.9003265327399586, 0.9238281021846115, 0.8003265327399586, 50.67037354402176, 175.23506363076496, 1000.800000864946, 1000.8003265326599, 0.830826159563432, 0.8003365327399585, 0.8515656668919402, 0.9729940750555356, 0.8003365327399585, 0.8003365327399586, 1.8003265327399585, 0.8071935357203198, 1.3003265327399585, 1.2003265327399586, 1.1245241609043495, 215.56060811614194, 0.8003365327399585, 0.9877186257774035, 0.9096809374131896, 4.2119694568059485, 0.8003365327399585, 1000.80032653274, 0.8003365327399585, 7.693408510878489, 10.80032653273916, 0.8795099522271821, 100.80032653273996, 10.800326521591268, 9.23282578370948, 7.388395759824789, 10.800324181981878]

    numVariables = length(u0)
    numTimePointsFor = Vector{Int64}(undef, numObservables)
    numTimePointsFor .= length.(observedAt)

    observablesTimeIndexIndices = Array{Vector{Int64}, 2}(undef, numObservables, numTimeSteps)
    for i = 1:numObservables
        for j = 1:numTimeSteps
            observablesTimeIndexIndices[i,j] = findall(x->x==j, observedAtIndex[i])
        end
    end

    constantValues = ConstantValues(observedObservable, measurementFor, minVariance, observablesTimeIndexIndices, observedAtIndex, 
                                    ts, numObservables, numParameters, numUsedParameters, numTimePointsFor, numTimeSteps, numVariables)

    scale = Array{ForwardDiff.Dual, 1}(undef, numObservables)
    #scale[[1,2]] .= 1.0
    shift = zeros(Float64, numObservables)# Array{Float64, 1}(undef, numObservables)
    shift[[3,4]] .+= 1.0
    variance = Array{ForwardDiff.Dual, 1}(undef, numObservables)
    #shift[[3,4]] .= 1.0
    staticParameters = StaticParameters(scale, shift, variance)

    #=
    h_bar = Array{Array{Float64, 1}, 1}(undef, numObservables)
    h_hat = Array{Array{Float64, 1}, 1}(undef, numObservables)
    costFor = Array{Float64, 1}(undef, numObservables)
    p_vector = Vector{Float64}(undef, numUsedParameters)
    mutatingArrays = MutatingArrays(h_bar, h_hat, costFor, p_vector)
    =#
    h_bar = Array{Array{ForwardDiff.Dual, 1}, 1}(undef, numObservables)
    h_hat = Array{Array{ForwardDiff.Dual, 1}, 1}(undef, numObservables)
    costFor = Array{ForwardDiff.Dual, 1}(undef, numObservables)
    p_vector = Vector{Float64}(undef, numUsedParameters)
    mutatingArrays = MutatingArrays(h_bar, h_hat, costFor, p_vector)

    f = (p) -> f_proto(prob, staticParameters, constantValues, mutatingArrays, p)

    chunkSize = 11

    for i in 1:numObservables
        scale[i] = zero(ForwardDiff.Dual{ForwardDiff.Tag{typeof(f), Float64},Float64, chunkSize})
        variance[i] = zero(ForwardDiff.Dual{ForwardDiff.Tag{typeof(f), Float64},Float64, chunkSize})
    end
    scale[[1,2]] .+= 1.0


    for i in 1:numObservables
        h_bar[i] = zeros(ForwardDiff.Dual{ForwardDiff.Tag{typeof(f), Float64},Float64, chunkSize}, numTimePointsFor[i])
        h_hat[i] = zeros(ForwardDiff.Dual{ForwardDiff.Tag{typeof(f), Float64},Float64, chunkSize}, numTimePointsFor[i])
        costFor[i] = zero(ForwardDiff.Dual{ForwardDiff.Tag{typeof(f), Float64},Float64, chunkSize})
    end


    #grad = ForwardDiff.gradient(f, test_orderedParameters)
    #println(grad)
    out = Vector{Float64}(undef, numUsedParameters)
    ForwardDiff.gradient!(out, f, test_orderedParameters)
    println(out)
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


    println("Starting optimisation")
    GradCalc_forwardDiff(filesAndPaths, timeEnd, relevantMeasurementData, observables)

end

main()