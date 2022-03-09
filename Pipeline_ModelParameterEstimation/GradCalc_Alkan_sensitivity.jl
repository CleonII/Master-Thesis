using ModelingToolkit, DifferentialEquations, BenchmarkTools, DataFrames, CSV, LSODA, Sundials, ODEInterface, ODEInterfaceDiffEq
using JuMP, NLopt, LinearAlgebra, DiffEqSensitivity, DiffEqFlux
using ModelingToolkit: varmap_to_vars

println("Done loading modules")

include(pwd() * "/Additional_functions/additional_tools.jl")
include(pwd() * "/Additional_functions/Solver_info.jl")


struct StaticParameters
    scale::Vector{Float64}
    shift::Vector{Float64}
    variance::Vector{Float64}
end

struct ConstantValues
    observedObservable::Vector{Int64}
    measurementFor::Vector{Vector{Float64}}
    minVariance::Float64
    observablesTimeIndexIndices::Array{Vector{Int64}, 2}
    observedAtIndex::Vector{Vector{Int64}}
    ts::Vector{Float64}
    variablesObservedBy::Vector{Vector{Int64}}
    #
    numObservables::Int64
    numParameters::Int64
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
    println(p_tuple)

    p_vector = mutatingArrays.p_vector
    p_vector .= p_tuple
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
    p_vector .= p_tuple

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
    grad[:] = dhdp * dfdu
    println(grad)   

    nothing
end



function GradCalc_sensitivity_eq(filesAndPaths, timeEnd, relevantMeasurementData, observables)
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
    numParameters = length(pars)
    orderedParameters = varmap_to_vars(pars, parameters(new_sys))
    orderedParameters .+= rand()
    test_orderedParameters = [0.9035043754031366, 10.800326532739959, 0.8344511564921195, 0.9413681449465665, 0.8491275291854607, 1.1567776661662026, 0.8003265327399586, 0.8905447897457612, 0.8057041743522487, 0.8003365327399586, 0.800336533086573, 0.8016809782200571, 4.201342833436018, 1000.800326532659, 1.1567777122376215, 1.8003265327398825, 1.8003265327398794, 0.8003365327399585, 0.9003265327399586, 0.9238281021846115, 0.8003265327399586, 50.67037354402176, 175.23506363076496, 1000.800000864946, 1000.8003265326599, 0.830826159563432, 0.8003365327399585, 0.8515656668919402, 0.9729940750555356, 0.8003365327399585, 0.8003365327399586, 1.8003265327399585, 0.8071935357203198, 1.3003265327399585, 1.2003265327399586, 1.1245241609043495, 215.56060811614194, 0.8003365327399585, 0.9877186257774035, 0.9096809374131896, 4.2119694568059485, 0.8003365327399585, 1000.80032653274, 0.8003365327399585, 7.693408510878489, 10.80032653273916, 0.8795099522271821, 100.80032653273996, 10.800326521591268, 9.23282578370948, 7.388395759824789, 10.800324181981878]
    #[1.8666279863282882, 10.800325983245347, 1.4849914433851856, 0.0, 0.07882578708199628, 2.2765182792654164, 0.0, 0.03537292965290244, 0.8057041743522487, 0.8003365327399586, 0.800336533086573, 0.8016809782200571, 4.208261180909193, 1000.8002040862448, 1.5241810171634373, 2.2661656037562476, 1.9231351564113515, 0.7941765327140669, 0.8963112375082627, 2.084685334577463, 0.7866795993355498, 50.67030144394519, 175.23499891290757, 1000.8017826414157, 1000.7963109698443, 0.8326079360332286, 0.8003381610913134, 1.1741817873312708, 0.575669389212187, 0.7734650612695332, 0.6871828727236781, 1.800326532739958, 0.9975373454097706, 1.2963401508219894, 1.3046421389407894, 0.9186651322476836, 215.56151063598094, 0.7603790855837445, 0.851366860474014, 1.8696322766123301, 4.218789878042009, 0.4546347191828532, 1000.7998405537494, 0.0, 7.519032585820315, 11.196735858380208, 0.17357304856409217, 100.803954578622, 10.779598392784923, 8.908264050383941, 7.489795067966958, 10.79468160240232]


    senProb = ODELocalSensitivityProblem(new_sys, prob.u0, tspan, orderedParameters)

    numVariables = length(u0)
    numTimePointsFor = Vector{Int64}(undef, numObservables)
    numTimePointsFor .= length.(observedAt)

    observablesTimeIndexIndices = Array{Array{Int64, 1}, 2}(undef, numObservables, numTimeSteps)
    for i = 1:numObservables
        for j = 1:numTimeSteps
            observablesTimeIndexIndices[i,j] = findall(x->x==j, observedAtIndex[i])
        end
    end

    variablesObservedBy = [[13, 18], [2, 8, 11, 12, 13, 18, 19, 20, 25, 29, 30, 36], [0], [0], [6, 10], [31, 33], [21, 28], [32, 34], [9, 35], [7, 24], [9, 35], [4, 14]]
    numVariablesObservedBy = [2, 12, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2]
    
    constantValues = ConstantValues(observedObservable, measurementFor, minVariance, observablesTimeIndexIndices, observedAtIndex, 
                                    ts, variablesObservedBy, numObservables, numParameters, numUsedParameters, numTimePointsFor, numTimeSteps, numVariables, numVariablesObservedBy)

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
    u = Array{Float64, 2}(undef, numVariables, numTimeSteps)
    mutatingArrays = MutatingArrays(dfdu, dhdp, dp, h_bar, h_hat, costFor, p_vector, u)

    f = (p_tuple...) -> f_proto(prob, staticParameters, constantValues, mutatingArrays, p_tuple...)   

    f_prime = (grad, p_tuple...) -> f_prime_proto(grad, senProb, staticParameters, constantValues, mutatingArrays, p_tuple...)


    # test first cost and grad at test_orderedParameters
    # cost: 47.333303744050696
    # grad: [-5.3182143165552676, 5.515530130306332e-7, -5.82086852469523, 3.7521742178916146, 0.44418628571104485, -1.762966208682973, 2.7063868066930308, -0.08352076007252987, 4.336808689942018e-19, 2.168404344971009e-19, -1.3010426069826053e-18, 3.0567087759227537e-19, -5.8784745398482716e-5, 2.2111025691617808e-5, -0.10343374323538441, 0.10080976815337757, 0.012841677960984607, 0.0011469332691993667, 0.004470247451825353, 0.6855534324271183, 0.01545960769334584, 7.942870966394309e-5, 7.060649830222212e-5, 5.159319581707668e-5, -0.00023485864224371874, 5.1593195817076664e-5, -2.3789424486946293e-6, 0.2000771484241956, 0.1306018374948179, 0.04250949832814352, 0.08605708929116086, -6.092607821075757e-19, 0.008445227773165528, -0.00032214328237251897, -0.010848558916109821, -0.00810496166007523, 3.23963827930824e-5, -0.00023211705648472909, 0.2408354729850729, -0.004237850096504894, -5.8438497677064455e-5, 0.34360585205990357, -4.969438210718985e-5, 1.0095047751672284, 0.0008281019919793915, -0.00030820500687614595, -0.06336619008821796, -0.00238476288533326, -3.0566021566127903e-5, -0.017195498749074372, -0.00253013601676205, -0.0015803689278361642]
    #=
    cost = f(test_orderedParameters...)
    println(cost)
    grad = zeros(Float64, length(test_orderedParameters))
    f_prime(grad, test_orderedParameters...)
    println(grad)
    =#
    
    model = Model(NLopt.Optimizer)
    #JuMP.register(model::Model, s::Symbol, dimension::Integer, f::Function, ∇f::Function, ∇²f::Function)
    register(model, :f, numUsedParameters, f, f_prime)
    set_optimizer_attribute(model, "algorithm", :LD_MMA)
    @variable(model, p[1:numUsedParameters] >= 0) # fix lower and upper bounds
    for i in 1:numUsedParameters
        set_start_value(p[i], test_orderedParameters[i])
    end
    @NLobjective(model, Min, f(p...))
    JuMP.optimize!(model)
    #@show termination_status(model)
    #@show primal_status(model)
    p_opt = [value(p[i]) for i=1:numUsedParameters]
    println("Done!")
    println(p_opt)
    cost_opt = objective_value(model)
    println(cost_opt)
    
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
    GradCalc_sensitivity_eq(filesAndPaths, timeEnd, relevantMeasurementData, observables)

end

main()