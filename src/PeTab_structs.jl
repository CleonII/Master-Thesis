"""
    PeTabModel

    Struct storing information about a PeTab-model. Create by the `setUpCostFunc` function.

    # Args
    `modelName`: PeTab model name (must match the xml-file name)
    `evalYmod`: Function to evaluate yMod for the log-likelhood. 
    `evalU0!`: Function that computes the initial u0 value for the ODE-system.
    `evalSd!`: Function that computes the standard deviation value for the log-likelhood.
    `odeSystem`: ModellingToolkit ODE-system for the PeTab model.
    `paramMap`: A map to correctly map model parameters to the ODE-system.
    `stateMap`: A map to correctly mapping the parameters to the u0 values.
    `paramNames`: Names of the model parameters (both fixed and those to be estimated).
    `stateNames`: Names of the model states.
    `dirModel`: Directory where the model.xml and PeTab files are stored.
    `pathMeasurementData`: Path to the measurementData PeTab file.
    `pathMeasurementData`: Path to the experimentaCondition PeTab file
    `pathMeasurementData`: Path to the observables PeTab file
    `pathMeasurementData`: Path to the parameters PeTab file

    See also: [`setUpCostFunc`]
"""
struct PeTabModel{F1<:Function, 
                  F2<:Function, 
                  F3<:Function,
                  F4<:Function,
                  F5<:Function,
                  F6<:Function,
                  F7<:Function,
                  F8<:Function,
                  F9<:Function,
                  S<:ODESystem,
                  T1<:Vector{<:Pair{Num, <:Union{AbstractFloat, Num}}}, 
                  T2<:Vector{<:Pair{Num, <:Union{AbstractFloat, Num}}}, 
                  T3<:Vector{<:Any}, 
                  T4<:Vector{<:Any}, 
                  C <: SciMLBase.DECallback,
                  FA<:Vector{<:Function}}
    modelName::String
    evalYmod::F1
    evalU0!::F2
    evalU0::F3
    evalSd!::F4
    evalDYmodDu::F5
    evalDSdDu!::F6
    evalDYmodDp::F7
    evalDSdDp!::F8
    getTStops::F9
    odeSystem::S
    paramMap::T1
    stateMap::T2
    paramNames::T3
    stateNames::T4
    dirModel::String
    pathMeasurementData::String
    pathExperimentalConditions::String
    pathObservables::String
    pathParameters::String
    callbackSet::C
    checkCallbackActive::FA
end


struct PeTabOpt{F1<:Function, 
                F2<:Function, 
                F3<:Function,
                F4<:Function, 
                F5<:Function, 
                F6<:Function, 
                F7<:Function, 
                F8<:Function}

    evalF::F1
    evalFZygote::F2
    evalGradF::F3
    evalGradFZygote::F4
    evalGradFAdjoint::F5
    evalGradFForwardEq::F6
    evalHess::F7
    evalHessApprox::F8
    nParamEst::Int64
    namesParam::Array{String, 1}
    paramVecNotTransformed::Vector{Float64}
    paramVecTransformed::Vector{Float64}
    lowerBounds::Vector{Float64}
    upperBounds::Vector{Float64}
    pathCube::String
    peTabModel::PeTabModel
end


"""
    ParamData

    Struct storing the data in the PeTab parameter-file in type-stable manner.

    Currently logScale notices whether or not parameters are estimated on the 
    log10 scale or not.

    See also: [`processParameterData`]
"""
struct ParamData{T1<:Array{<:AbstractFloat}, 
                 T2<:Array{<:String, 1}, 
                 T3<:Array{Bool, 1}, 
                 T4<:Signed}

    # TODO: logScale make symbol to support more transformations 
    paramVal::T1
    lowerBounds::T1
    upperBounds::T1
    parameterID::T2
    logScale::T3
    shouldEst::T3
    nParamEst::T4
end


"""
    MeasurementData

    Struct storing the data in the PeTab measurementData-file in type-stable manner.

    Transform data supports log and log10 transformations of the data.  

    See also: [`processMeasurementData`]
"""
struct MeasurementData{T6<:Array{<:Union{<:String, <:AbstractFloat}, 1}}
                    
    yObsNotTransformed::Vector{Float64}
    yObsTransformed::Vector{Float64}
    tObs::Vector{Float64}
    observebleID::Vector{String}
    conditionId::Vector{String}  # Sum of pre-eq + simulation-cond id 
    sdParams::T6
    transformData::Vector{Symbol} # Only done once 
    obsParam::Vector{String}
    tVecSave::Dict{String, Vector{Float64}}
    iTObs::Vector{Int64}
    iPerConditionId::Dict{String, Vector{Int64}}
    preEqCond::Vector{String}
    simCond::Vector{String}
    iGroupedTObs::Dict{String, Vector{Vector{Int64}}}
end


"""
    SimulationInfo

    Struct storing simulation (forward ODE-solution) information. Specifcially 
    stores the experimental ID:s from the experimentalCondition - PeTab file;
    firstExpIds (preequilibration ID:s), the shiftExpIds (postequilibration), and
    simulateSS (whether or not to simulate ODE-model to steady state). Further 
    stores a solArray with the ODE solution where conditionIdSol of the ID for 
    each forward solution. It also stores for each experimental condition which 
    time-points we have observed data at

    See also: [`getSimulationInfo`]
"""
struct SimulationInfo{T4<:Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1}, 
                      T5<:Dict{<:String, <:Vector{<:Float64}}, 
                      T6<:Vector{<:SciMLBase.DECallback}, 
                      T7<:Union{<:SciMLSensitivity.AbstractForwardSensitivityAlgorithm, <:SciMLSensitivity.AbstractAdjointSensitivityAlgorithm}}
    firstExpIds::Vector{String}
    shiftExpIds::Vector{Vector{String}}
    preEqIdSol::Vector{String}
    postEqIdSol::Vector{String}
    conditionIdSol::Vector{String}
    tMaxForwardSim::Vector{Float64}
    simulateSS::Bool
    solArray::T4
    solArrayGrad::T4
    solArrayPreEq::T4
    absTolSS::Float64
    relTolSS::Float64
    tVecSave::T5
    callbacks::T6
    sensealg::T7
end


"""
    ParamMap

    Struct which makes out a map to correctly for an observation extract the correct observable 
    or sd-param via the getObsOrSdParam function when computing the likelihood. Correctly built 
    by `buildMapParameters`, and is part of the ParameterIndices-struct.

    For noise or observable parameters belong to an observation, e.g (obsParam1, obsParam2), 
    this struct stores which parameters should be estimtated, and for those parameters which 
    index they correspond to in the parameter estimation vector. For constant parameters 
    the struct stores the values. 

    See also: [`getIndicesParam`, `buildMapParameters`]
"""
struct ParamMap
    shouldEst::Array{Bool, 1}
    indexUse::Array{Int64, 1}
    valuesConst::Vector{Float64}
    nParam::Int64
end


struct MapExpCond
    condID::String
    expCondParamConstVal::Vector{Float64}
    iOdeProbParamConstVal::Vector{Int64}
    expCondStateConstVal::Vector{Float64}
    iOdeProbStateConstVal::Vector{Int64}
    iDynEstVec::Vector{Int64}
    iOdeProbDynParam::Vector{Int64}
end


struct MapDynParEst
    iDynParamInSys::Vector{Int64}
    iDynParamInVecEst::Vector{Int64}
end


"""
    ParameterIndices

    Struct storing names and mapping indices for mapping the parameter provided 
    to the optimizers correctly. 
    
    Optimizers require a single vector input of parameters (pVecEst). However, the PeTab 
    model has three kind of parameters, Dynmaic (part of the ODE-system), 
    Observable (only part of the observation model) and Standard-deviation 
    (only part of the standard deviation expression in the log-likelhood). This 
    struct stores mapping indices (starting with i) to map pVecEst 
    correctly when computing the likelihood (e.g map the SD-parameters in pVecEst
    correctly to a vector of SD-vals). Also stores the name of each parameter.

    Furthermore, when computing yMod or SD the correct observable and sd parameters 
    has to be used for each observation. The mapArrays effectively contains precomputed  
    maps allowing said parameter to be effectively be extracted by the getObsOrSdParam 
    function. 

    See also: [`getIndicesParam`, `ParamMap`]
"""
struct ParameterIndices{T4<:Array{<:ParamMap, 1}, 
                        T5<:MapDynParEst, 
                        T6<:Array{<:MapExpCond, 1}}

    iDynParam::Vector{Int64}
    iObsParam::Vector{Int64}
    iSdParam::Vector{Int64}
    iSdObsNonDynPar::Vector{Int64}
    iNonDynParam::Vector{Int64}
    namesDynParam::Vector{String}
    namesObsParam::Vector{String}
    namesSdParam::Vector{String}
    namesSdObsNonDynPar::Vector{String}
    namesNonDynParam::Vector{String}
    namesParamEst::Vector{String}
    indexObsParamMap::Vector{Int64}
    indexSdParamMap::Vector{Int64}
    mapArrayObsParam::T4
    mapArraySdParam::T4
    mapDynParEst::T5
    mapExpCond::T6
    constParamPerCond::Vector{Vector{Float64}}
end


struct PriorInfo{T1 <: Vector{<:Function}}
    logpdf::T1
    priorOnParamScale::Vector{Bool}
    hasPriors::Bool
end