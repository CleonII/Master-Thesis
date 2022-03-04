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
    isObserved::Vector{Bool}
    iterations::Int64
    measurementFor::Vector{Vector{Float64}}
    minVariance::Float64
    observablesTimeIndexIndices::Array{Vector{Int64}, 2}
    observedAtIndex::Vector{Vector{Int64}}
    ts::Vector{Float64}
    #
    numObservables::Int64
    numParameters::Int64
    numUsedParameters::Int64
    numTimePointsFor::Vector{Float64}
    numTimeSteps::Int64
    numVariables::Int64
end

struct MutatingArrays
    h_bar::Vector{Vector{Float64}}
    h_hat::Vector{Vector{Float64}}
    costFor::Vector{Float64}
    p_vector::Vector{Float64}
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

    for i = 1:numObservables
        h_hat[i] = scale[i] * h_bar[i] .+ shift[i]
        variance[i] = (dot(measurementFor[i], measurementFor[i]) - 2 * dot(measurementFor[i], h_hat[i]) + dot(h_hat[i], h_hat[i])) / num_t[i]
        variance[i] = max(variance[i], minVariance)
        costFor[i] = (num_t[i]/2) * log(2*pi*variance[i]) + (dot(measurementFor[i], measurementFor[i]) - 2*dot(measurementFor[i], h_hat[i]) + dot(h_hat[i], h_hat[i])) / (2 * variance[i])
    end

    println(sum(costFor[constantValues.isObserved]))
    return sum(costFor[constantValues.isObserved])
end

# the total cost, G, is the integral/sum of the costs at each time step, g
# here we have the derivative, dg, of the cost at each time step
# Note: adjoint_sensitivities seems to use the negative gradient due to a sign error
# Each out[i] corresponds with the gradient of the cost with regards to state i at time t
function dg_proto!(out,u,p,t,i, staticParameters, constantValues, mutatingArrays) 
    scale = staticParameters.scale
    variance = staticParameters.variance

    measurementFor = constantValues.measurementFor
    observablesTimeIndexIndices = constantValues.observablesTimeIndexIndices

    h_hat = mutatingArrays.h_hat
    #out .= 0.0
    out[1] =  0.0                                                                                                                                # Cells_Dead
    out[2] =  (scale[2]/variance[2])   * (sum(@view measurementFor[2][observablesTimeIndexIndices[2,i]])   - sum(@view h_hat[2][observablesTimeIndexIndices[2,i]]))            # Cells_Apo1: 2
    out[3] =  0.0                                                                                                                                # pATR_G2
    out[4] =  (scale[12]/variance[12]) * (sum(@view measurementFor[12][observablesTimeIndexIndices[12,i]]) - sum(@view h_hat[12][observablesTimeIndexIndices[12,i]]))      # yH2AX_G2: 12
    out[5] =  0.0#(scale[3]/variance[3]) * (sum(@view measurementFor[3][observablesTimeIndexIndices[3,i]])   - sum(@view h_hat[3][observablesTimeIndexIndices[3,i]]))            # Wip1_mRNA_S: 3
    out[6] =  (scale[5]/variance[5])   * (sum(@view measurementFor[5][observablesTimeIndexIndices[5,i]])   - sum(@view h_hat[5][observablesTimeIndexIndices[5,i]]))            # pATM_S: 5
    out[7] =  (scale[10]/variance[10]) * (sum(@view measurementFor[10][observablesTimeIndexIndices[10,i]]) - sum(@view h_hat[10][observablesTimeIndexIndices[10,i]]))      # p21_G2: 10
    out[8] =  (scale[2]/variance[2])   * (sum(@view measurementFor[2][observablesTimeIndexIndices[2,i]])   - sum(@view h_hat[2][observablesTimeIndexIndices[2,i]]))            # Cells_SSBDamage_S: 2
    out[9] =  (scale[9]/variance[9])   * (sum(@view measurementFor[9][observablesTimeIndexIndices[9,i]])   - sum(@view h_hat[9][observablesTimeIndexIndices[9,i]])) + 
              (scale[11]/variance[11]) * (sum(@view measurementFor[11][observablesTimeIndexIndices[11,i]]) - sum(@view h_hat[11][observablesTimeIndexIndices[11,i]]))      # pp53_S: 9, 11
    out[10] = (scale[5]/variance[5])   * (sum(@view measurementFor[5][observablesTimeIndexIndices[5,i]])   - sum(@view h_hat[5][observablesTimeIndexIndices[5,i]]))           # pATM_G2: 5
    out[11] = (scale[2]/variance[2])   * (sum(@view measurementFor[2][observablesTimeIndexIndices[2,i]])   - sum(@view h_hat[2][observablesTimeIndexIndices[2,i]]))           # Cells_Apo2: 2
    out[12] = (scale[2]/variance[2])   * (sum(@view measurementFor[2][observablesTimeIndexIndices[2,i]])   - sum(@view h_hat[2][observablesTimeIndexIndices[2,i]]))           # Cells_Cycle_G2: 2
    out[13] = (scale[1]/variance[1])   * (sum(@view measurementFor[1][observablesTimeIndexIndices[1,i]])   - sum(@view h_hat[1][observablesTimeIndexIndices[1,i]])) + 
              (scale[2]/variance[2])   * (sum(@view measurementFor[2][observablesTimeIndexIndices[2,i]])   - sum(@view h_hat[2][observablesTimeIndexIndices[2,i]]))           # Cells_Apo_ReOx: 1,2
    out[14] = (scale[12]/variance[12]) * (sum(@view measurementFor[12][observablesTimeIndexIndices[12,i]]) - sum(@view h_hat[12][observablesTimeIndexIndices[12,i]]))     # yH2AX_S: 12
    out[15] = 0.0#(scale[4]/variance[4]) * (sum(@view measurementFor[4][observablesTimeIndexIndices[4,i]])   - sum(@view h_hat[4][observablesTimeIndexIndices[4,i]]))           # p21_mRNA_S: 4
    out[16] = 0.0                                                                                                                               # Wip1_S
    out[17] = 0.0#(scale[3]/variance[3]) * (sum(@view measurementFor[3][observablesTimeIndexIndices[3,i]])   - sum(@view h_hat[3][observablesTimeIndexIndices[3,i]]))           # Wip1_mRNA_G2: 3
    out[18] = (scale[1]/variance[1])   * (sum(@view measurementFor[1][observablesTimeIndexIndices[1,i]])   - sum(@view h_hat[1][observablesTimeIndexIndices[1,i]])) + 
              (scale[2]/variance[2])   * (sum(@view measurementFor[2][observablesTimeIndexIndices[2,i]])   - sum(@view h_hat[2][observablesTimeIndexIndices[2,i]]))           # Cells_Apo4: 1,2
    out[19] = (scale[2]/variance[2])   * (sum(@view measurementFor[2][observablesTimeIndexIndices[2,i]])   - sum(@view h_hat[2][observablesTimeIndexIndices[2,i]]))           # Cells_DSBDamage_G2: 2
    out[20] = (scale[2]/variance[2])   * (sum(@view measurementFor[2][observablesTimeIndexIndices[2,i]])   - sum(@view h_hat[2][observablesTimeIndexIndices[2,i]]))           # Cells_Apo: 2
    out[21] = (scale[7]/variance[7])   * (sum(@view measurementFor[7][observablesTimeIndexIndices[7,i]])   - sum(@view h_hat[7][observablesTimeIndexIndices[7,i]]))           # pChk2_G2: 7
    out[22] = 0.0                                                                                                                               # Wip1_G2 
    out[23] = 0.0                                                                                                                               # pATR_S
    out[24] = (scale[10]/variance[10]) * (sum(@view measurementFor[10][observablesTimeIndexIndices[10,i]]) - sum(@view h_hat[10][observablesTimeIndexIndices[10,i]]))     # p21_S: 10
    out[25] = (scale[2]/variance[2])   * (sum(@view measurementFor[2][observablesTimeIndexIndices[2,i]])   - sum(@view h_hat[2][observablesTimeIndexIndices[2,i]]))           # Cells_DSBDamage_S: 2
    out[26] = 0.0#(scale[4]/variance[4]) * (sum(@view measurementFor[4][observablesTimeIndexIndices[4,i]])   - sum(@view h_hat[4][observablesTimeIndexIndices[4,i]]))           # p21_mRNA_G2: 4
    out[27] = 0.0                                                                                                                               # Space
    out[28] = (scale[7]/variance[7])   * (sum(@view measurementFor[7][observablesTimeIndexIndices[7,i]])   - sum(@view h_hat[7][observablesTimeIndexIndices[7,i]]))           # pChk2_S: 7
    out[29] = (scale[2]/variance[2])   * (sum(@view measurementFor[2][observablesTimeIndexIndices[2,i]])   - sum(@view h_hat[2][observablesTimeIndexIndices[2,i]]))           # Cells_Cycle_S: 2
    out[30] = (scale[2]/variance[2])   * (sum(@view measurementFor[2][observablesTimeIndexIndices[2,i]])   - sum(@view h_hat[2][observablesTimeIndexIndices[2,i]]))           # Cells_Apo3: 2
    out[31] = (scale[6]/variance[6])   * (sum(@view measurementFor[6][observablesTimeIndexIndices[6,i]])   - sum(@view h_hat[6][observablesTimeIndexIndices[6,i]]))           # pChk1_G2: 6
    out[32] = (scale[8]/variance[8])   * (sum(@view measurementFor[8][observablesTimeIndexIndices[8,i]])   - sum(@view h_hat[8][observablesTimeIndexIndices[8,i]]))           # pDNAPK_S: 8
    out[33] = (scale[6]/variance[6])   * (sum(@view measurementFor[6][observablesTimeIndexIndices[6,i]])   - sum(@view h_hat[6][observablesTimeIndexIndices[6,i]]))           # pChk1_S: 6
    out[34] = (scale[8]/variance[8])   * (sum(@view measurementFor[8][observablesTimeIndexIndices[8,i]])   - sum(@view h_hat[8][observablesTimeIndexIndices[8,i]]))           # pDNAPK_G2: 8
    out[35] = (scale[9]/variance[9])   * (sum(@view measurementFor[9][observablesTimeIndexIndices[9,i]])   - sum(@view h_hat[9][observablesTimeIndexIndices[9,i]])) + 
              (scale[11]/variance[11]) * (sum(@view measurementFor[11][observablesTimeIndexIndices[11,i]]) - sum(@view h_hat[11][observablesTimeIndexIndices[11,i]]))     # pp53_G2: 9, 11
    out[36] = (scale[2]/variance[2])   * (sum(@view measurementFor[2][observablesTimeIndexIndices[2,i]])   - sum(@view h_hat[2][observablesTimeIndexIndices[2,i]]))           # Cells: 2
end

function f_prime_proto(grad, prob, constantValues, mutatingArrays, dg!, p_tuple...)
    println("In f_prime()")
    p_vector = mutatingArrays.p_vector
    p_vector .= p_tuple
    _prob = remake(prob, p = p_vector) 

    sol = OrdinaryDiffEq.solve(_prob, AutoVern7(Rodas5()), reltol=1e-8, abstol=1e-8, saveat = constantValues.ts, 
            sensealg=InterpolatingAdjoint(checkpointing=true))

    ~, grad[:] = adjoint_sensitivities(sol, Trapezoid(), dg!, constantValues.ts, sensalg = BacksolveAdjoint(autojacvec = ReverseDiffVJP(true)))

    nothing
end


function GradCalc_adjoint(filesAndPaths, timeEnd, relevantMeasurementData, observables, iterations)
    minVariance = 1e-2

    observableIDs = relevantMeasurementData[:,1]
    ts = sort(unique(relevantMeasurementData[:,5]))
    numTimeSteps = length(ts)
    numObservables = length(observables)

    isObserved = [observables[i] in unique(observableIDs) for i=1:numObservables]

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

    numVariables = length(u0)
    numTimePointsFor = Vector{Float64}(undef, numObservables)
    numTimePointsFor .= length.(observedAt)

    observablesTimeIndexIndices = Array{Vector{Int64}, 2}(undef, numObservables, numTimeSteps)
    for i = 1:numObservables
        for j = 1:numTimeSteps
            observablesTimeIndexIndices[i,j] = findall(x->x==j, observedAtIndex[i])
        end
    end

    constantValues = ConstantValues(isObserved, iterations, measurementFor, minVariance, observablesTimeIndexIndices, observedAtIndex, 
                                    ts, numObservables, numParameters, numUsedParameters, numTimePointsFor, numTimeSteps, numVariables)

    scale = Array{Float64, 1}(undef, numObservables)
    scale[[1,2]] .= 1.0
    shift = zeros(Float64, numObservables)
    variance = Array{Float64, 1}(undef, numObservables)
    shift[[3,4]] .= 1.0
    staticParameters = StaticParameters(scale, shift, variance)

    h_bar = Array{Array{Float64, 1}, 1}(undef, numObservables)
    h_hat = Array{Array{Float64, 1}, 1}(undef, numObservables)
    costFor = Array{Float64, 1}(undef, numObservables)
    p_vector = Vector{Float64}(undef, numUsedParameters)
    mutatingArrays = MutatingArrays(h_bar, h_hat, costFor, p_vector)

    #_pars = [Pair(pars[i].first, pars[i].second+abs(rand())) for i in 1:numParameters]

    #p_tuple = Tuple(orderedParameters)

    f = (p_tuple...) -> f_proto(prob, staticParameters, constantValues, mutatingArrays, p_tuple...)    
    #totcost = costFunc(p_tuple...)
    #println(totcost)

    dg! = (out, u, p, t, i) -> dg_proto!(out, u, p, t, i, staticParameters, constantValues, mutatingArrays)

    #grad = Vector{Float64}(undef, numUsedParameters)

    f_prime = (grad, p_tuple...) -> f_prime_proto(grad, prob, constantValues, mutatingArrays, dg!, p_tuple...)
    #f_prime(grad, p_tuple...)

    #println(grad)


    model = Model(NLopt.Optimizer)
    #JuMP.register(model::Model, s::Symbol, dimension::Integer, f::Function, ∇f::Function, ∇²f::Function)
    register(model, :f, numUsedParameters, f, f_prime)
    set_optimizer_attribute(model, "algorithm", :LD_MMA)
    @variable(model, p[1:numUsedParameters] >= 0)
    for i in 1:numUsedParameters
        set_start_value(p[i], orderedParameters[i])
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
    
    #return p_opt, cost_opt
    
    
    nothing
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

    iterations = 5

    println("Starting optimisation")
    GradCalc_adjoint(filesAndPaths, timeEnd, relevantMeasurementData, observables, iterations)

end

main()