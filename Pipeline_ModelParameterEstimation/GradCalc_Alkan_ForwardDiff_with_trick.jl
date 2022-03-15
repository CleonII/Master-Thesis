using ModelingToolkit, DifferentialEquations, BenchmarkTools, DataFrames, CSV, LSODA, Sundials, ODEInterface, ODEInterfaceDiffEq
using LinearAlgebra, DiffEqSensitivity, DiffEqFlux, ForwardDiff, Random
using ModelingToolkit: varmap_to_vars
using ForwardDiff: GradientConfig, Chunk
using ProfileView, ProgressMeter

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
    dualP_vector::Vector{ForwardDiff.Dual}
    grad::Vector{Float64}
end

struct FilesAndPaths
    modelFile::String
    modelPath::String
    writefile_solver::String
    writefile_sensealg::String
end

mutable struct Adam
    theta::AbstractArray{Float64} # Parameter array
    grad::Function                # Gradient function
    loss::Float64
    m::AbstractArray{Float64}     # First moment
    v::AbstractArray{Float64}     # Second moment
    b1::Float64                   # Exp. decay first moment
    b2::Float64                   # Exp. decay second moment
    a::Float64                    # Step size
    eps::Float64                  # Epsilon for stability
    t::Int                        # Time step (iteration)
end
  
# Outer constructor
function Adam(theta::AbstractArray{Float64}, grad::Function; a::Float64=0.005)
    loss = 0
    m   = zeros(size(theta))
    v   = zeros(size(theta))
    b1  = 0.9
    b2  = 0.999
    a   = a
    eps = 1e-8
    t   = 0
    Adam(theta, grad, loss, m, v, b1, b2, a, eps, t)
end

# Step function with optional keyword arguments for the data passed to grad()
function step!(opt::Adam)
    opt.t += 1
    gt, opt.loss    = opt.grad(opt.theta)
    opt.m = opt.b1 .* opt.m + (1 - opt.b1) .* gt
    opt.v = opt.b2 .* opt.v + (1 - opt.b2) .* gt .^ 2
    mhat = opt.m ./ (1 - opt.b1^opt.t)
    vhat = opt.v ./ (1 - opt.b2^opt.t)
    opt.theta -= opt.a .* (mhat ./ (sqrt.(vhat) .+ opt.eps))
end

# F - function on p 
function unbiased_grad(p::T1, f::Function, result)::Tuple{Array{Float64, 1}, Float64} where T1<:Array{<:AbstractFloat, 1}

    v::Array{Float64, 1} = randn(length(p)) # v ~ MV(0, I)

    g = r -> f(p + r*v)
    result = ForwardDiff.derivative!(result, g, 0.0)
    uGrad = DiffResults.gradient(result) * v
    loss = DiffResults.value(result)
    return uGrad, loss
end

function costFunc(prob::ODEProblem, staticParameters::StaticParameters, constantValues::ConstantValues, mutatingArrays::MutatingArrays, p::Vector{T1})::T1 where T1 <: ForwardDiff.Dual

    dualP_vector = mutatingArrays.dualP_vector
    dualP_vector[:] = convert.(eltype(p), prob.p)
    dualP_vector[constantValues.optParameterIndices] .= p

    _prob = remake(prob, u0 = convert.(eltype(p), prob.u0), p = dualP_vector)

    sol = OrdinaryDiffEq.solve(_prob, Rodas5(), reltol=1e-8, abstol=1e-8, saveat = constantValues.ts)
    
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

    pWithInputs = prob.p 
    pWithInputs[inputParameterIndices] .= inputParameterValues
    prob = remake(prob, p=pWithInputs)

    orderedParameters = varmap_to_vars(pars, parameters(new_sys))[optParameterIndices]
    #orderedParameters = varmap_to_vars(pars, parameters(new_sys))
    orderedParameters .+= (rand() - 0.5)*5
    orderedParameters[orderedParameters .< 0.0] .= 0.0
    orderedParameters .= rand()*5
    #test_orderedParameters = [0.5963242072135031, 10.493146364550325, 0.5272709883024861, 0.6341879767569332, 0.5419473609958273, 0.8495974979765691, 0.4931463645503251, 0.5833646215561278, 0.4985240061626153, 0.49315636455032524, 0.49315636489693954, 0.49450081003042373, 3.8941626652463848, 1000.4931463644693, 0.8495975440479882, 1.493146364550249, 1.493146364550246, 0.49315636455032513, 0.5931463645503251, 0.6166479339949781, 0.4931463645503251, 50.363193375832125, 174.92788346257532, 1000.4928206967563, 1000.4931463644702, 0.5236459913737985, 0.49315636455032513, 0.5443854987023068, 0.6658139068659021, 0.49315636455032513, 0.4931563645503252, 1.4931463645503251, 0.5000133675306864, 0.9931463645503251, 0.8931463645503251, 0.8173439927147161, 215.2534279479523, 0.49315636455032513, 0.6805384575877701, 0.6025007692235561, 3.904789288616315, 0.49315636455032513, 1000.4931463645503, 0.49315636455032513, 7.386228342688855, 10.493146364549526, 0.5723297840375486, 100.49314636455033, 10.493146353401634, 8.925645615519846, 7.081215591635155, 10.493144013792245]
    #orderedParameters[inputParameterIndices] .= 0.0
    #orderedParameters = test_orderedParameters[optParameterIndices]

    numTimePointsFor = Vector{Int64}(undef, numObservables)
    numTimePointsFor .= length.(observedAt)

    constantValues = ConstantValues(observedObservable, measurementFor, minVariance, observedAtIndex, optParameterIndices, 
                                    ts, numObservables, numUsedParameters, numTimePointsFor, numTimeSteps)

    usedType = ForwardDiff.Dual
    scale = Vector{Float64}(undef, numObservables)
    scale[[1,2]] .= 1.0
    shift = zeros(Float64, numObservables)
    shift[[3,4]] .= 1.0
    variance = Vector{Float64}(undef, numObservables)
    staticParameters = StaticParameters(scale, shift, variance)

    h_bar = Vector{Vector{usedType}}(undef, numObservables)
    h_hat = Vector{Vector{usedType}}(undef, numObservables)
    costFor = Vector{usedType}(undef, numObservables)
    dualP_vector = Vector{usedType}(undef, numUsedParameters)
    grad = Vector{Float64}(undef, numOptParameters)
    mutatingArrays = MutatingArrays(h_bar, h_hat, costFor, dualP_vector, grad)

    f_cost = (p) -> costFunc(prob, staticParameters, constantValues, mutatingArrays, p)

    result = DiffResults.DiffResult(0.0, 0.0)

    grad_g = (p) -> unbiased_grad(p, f_cost, result)





    n_it = 7500
    theta_save = Array{Float64, 2}(undef, (numOptParameters, n_it))
    loss_val = Array{Float64, 1}(undef, n_it)

    adam_opt = Adam(orderedParameters, grad_g, a=0.001)
    
    @showprogress 1 "Running ADAM... " for i in 1:n_it
        step!(adam_opt)
        theta_save[:, i] .= adam_opt.theta
        loss_val[i] = adam_opt.loss
    end

    return theta_save, loss_val
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
    theta_save, loss_val = GradCalc_forwardDiff(filesAndPaths, timeEnd, relevantMeasurementData, observables, inputParameterValues)

    plotly()
    plot(loss_val)
end

main()