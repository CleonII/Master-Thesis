using ModelingToolkit, DifferentialEquations, BenchmarkTools, DataFrames, CSV, LSODA, Sundials, ODEInterface, ODEInterfaceDiffEq
using JuMP, NLopt, LinearAlgebra, DiffEqSensitivity, DiffEqFlux, ForwardDiff
using ModelingToolkit: varmap_to_vars
using ForwardDiff: GradientConfig, Chunk
using ProfileView

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
    observedAtIndex::Vector{Vector{Int64}}
    optParameterIndices::Vector{Int64}
    ts::Vector{Float64}
    #
    numObservables::Int64
    numUsedParameters::Int64
    numTimePointsFor::Vector{Int64}
    numTimeSteps::Int64
end

struct MutatingArrays
    h_bar::Vector{Vector{ForwardDiff.Dual}}
    h_hat::Vector{Vector{ForwardDiff.Dual}}
    costFor::Vector{ForwardDiff.Dual}
    p_vector::Vector{Float64}
    dualP_vector::Vector{ForwardDiff.Dual}
    grad::Vector{Float64}
end

struct FilesAndPaths
    modelFile::String
    modelPath::String
    writefile_solver::String
    writefile_sensealg::String
end

function costFunc(prob::ODEProblem, staticParameters::StaticParameters, constantValues::ConstantValues, mutatingArrays::MutatingArrays, p::Vector{T1})::T1 where T1 <: ForwardDiff.Dual
    
    dualP_vector = mutatingArrays.dualP_vector
    dualP_vector[constantValues.optParameterIndices] .= p

    _prob = remake(prob, u0 = convert.(eltype(p), prob.u0), p = dualP_vector)

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
        scale[i] = ((dot(h_bar[i], measurementFor[i]) - shift[i] * sum(h_bar[i])) / dot(h_bar[i], h_bar[i])).value
    end
    replace!(scale, NaN=>0.0) # necessary when h_bar is zeros 
    # Mayby replace 0.0 (for non zero h_bar) with a lower bound

    for i = 1:numObservables
        h_hat[i] = scale[i] * h_bar[i] .+ shift[i]
        variance[i] = ((dot(measurementFor[i], measurementFor[i]) - 2 * dot(measurementFor[i], h_hat[i]) + dot(h_hat[i], h_hat[i])) / num_t[i]).value
        variance[i] = max(variance[i], minVariance)
        costFor[i] = (num_t[i]/2) * log(2*pi*variance[i]) + (dot(measurementFor[i], measurementFor[i]) - 2*dot(measurementFor[i], h_hat[i]) + dot(h_hat[i], h_hat[i])) / (2 * variance[i])
    end

    return sum(costFor[constantValues.observedObservable])
end


function f_cost_proto(result::DiffResults.MutableDiffResult, f::Function, cfg::GradientConfig, mutatingArrays::MutatingArrays, p_tuple...)::Float64
    p = mutatingArrays.p_vector
    p .= p_tuple

    ForwardDiff.gradient!(result, f, p, cfg)
    mutatingArrays.grad[:] = DiffResults.gradient(result)
    
    println(DiffResults.value(result), "\n")
    return DiffResults.value(result)
end

function f_grad_proto(grad, mutatingArrays, p_tuple...)
    grad[:] = mutatingArrays.grad
    nothing
end


function GradCalc_forwardDiff(filesAndPaths, timeEnd, relevantMeasurementData, observables, inputParameterValues)
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

    numTimePointsFor = Vector{Int64}(undef, numObservables)
    numTimePointsFor .= length.(observedAt)

    constantValues = ConstantValues(observedObservable, measurementFor, minVariance, observedAtIndex, optParameterIndices,
                                    ts, numObservables, numUsedParameters, numTimePointsFor, numTimeSteps)

    
    scale = Vector{Float64}(undef, numObservables)
    scale[[1,2]] .= 1.0
    shift = zeros(Float64, numObservables)
    shift[[3,4]] .= 1.0
    variance = Vector{Float64}(undef, numObservables)
    staticParameters = StaticParameters(scale, shift, variance)

    usedType = ForwardDiff.Dual
    h_bar = Vector{Vector{usedType}}(undef, numObservables)
    h_hat = Vector{Vector{usedType}}(undef, numObservables)
    costFor = Vector{usedType}(undef, numObservables)
    p_vector = Vector{Float64}(undef, numOptParameters)
    dualP_vector = Vector{usedType}(undef, numUsedParameters)
    grad = Vector{Float64}(undef, numOptParameters)
    mutatingArrays = MutatingArrays(h_bar, h_hat, costFor, p_vector, dualP_vector, grad)

    f_cost = (p) -> costFunc(prob, staticParameters, constantValues, mutatingArrays, p)

    chunkSize = 13
    cfg = GradientConfig(f_cost, orderedParameters, Chunk{chunkSize}())

    dualP_vector[inputParameterIndices] .= inputParameterValues .* one(ForwardDiff.Dual{ForwardDiff.Tag{typeof(f_cost), Float64},Float64, chunkSize})

    #
    result = DiffResults.GradientResult(orderedParameters::Vector{Float64})

    f = (p_tuple...) -> f_cost_proto(result, f_cost, cfg, mutatingArrays, p_tuple...)
    f_grad = (grad, p_tuple...) -> f_grad_proto(grad, mutatingArrays, p_tuple...)

    #cost = f(orderedParameters...)
    # 29.21614392324067
    #grad2 = Vector{Float64}(undef, numOptParameters)
    #f_grad(grad2, orderedParameters...)
    #println(grad2)
    # [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.6406411214537998e-16, 0.0, 10.276324153033542, 12.216214038184734, 0.0, 0.0, 0.0, -1.7268649092217891, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.275811836644368e-5, 0.0, 0.0, 0.0, 0.0]

    # test first cost and grad at test_orderedParameters
    # cost: 47.33330405467803
    # grad: [-5.318215299039426, 5.515534973056204e-7, -5.820869002075568, 3.752174833073845, 0.44418583163493813, -1.7629644064642325, 2.7063844789080984, -0.08352082360898885, 0.0, -3.009265538105056e-34, 1.2037062152420224e-33, -1.504632769052528e-34, -5.8784718997286444e-5, 2.2111064355330217e-5, -0.10343392410137235, 0.10080959381932418, 0.01284167627301536, 0.001146966953022865, 0.004470244153038387, 0.6855533519457021, 0.015459596285039507, 7.94286510500958e-5, 7.060644619866742e-5, 5.163092560927721e-5, -0.0002348586275941731, 5.1630925609287945e-5, -2.3425009389426196e-6, 0.20007711244312204, 0.13060181312925243, 0.04250949099587211, 0.0860570749829805, -5.877350297371893e-19, 0.008445234457328957, -0.00032214327389440785, -0.010848559062022471, -0.008104965953355267, 3.239638604151488e-5, -0.00023211706158480746, 0.2408354592904409, -0.00423789046207008, -5.843846087025392e-5, 0.3436057744244619, -4.9694372697670505e-5, 1.0095043641393266, 0.0008280843800659434, -0.0003082048401328952, -0.06336618742551742, -0.002384762259748042, -3.056610602007926e-5, -0.01719550591281074, -0.002530135631328737, -0.0015803688265661744]
    #=
    cost = f(test_orderedParameters...)
    grad2 = Vector{Float64}(undef, numUsedParameters)
    f_grad(grad2, test_orderedParameters...)
    println(cost)
    println(grad2)
    =#
    
    model = Model(NLopt.Optimizer)
    #JuMP.register(model::Model, s::Symbol, dimension::Integer, f::Function, ∇f::Function, ∇²f::Function)
    register(model, :f, numOptParameters, f, f_grad)
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
    GradCalc_forwardDiff(filesAndPaths, timeEnd, relevantMeasurementData, observables, inputParameterValues)

end

main()