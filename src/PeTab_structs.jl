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
struct PeTabModel{T1<:Vector{<:Pair{Num, <:Union{AbstractFloat, Num}}}, 
                  T2<:Vector{<:Pair{Num, <:Union{AbstractFloat, Num}}},
                  T3<:Vector{Sym{Real, Base.ImmutableDict{DataType, Any}}}, 
                  T4<:Vector{Any}}
    modelName::String
    evalYmod::Function 
    evalU0!::Function
    evalSd!::Function
    odeSystem::ODESystem 
    paramMap::T1
    stateMap::T2
    paramNames::T3
    stateNames::T4
    dirModel::String
    pathMeasurementData::String
    pathExperimentalConditions::String
    pathObservables::String
    pathParameters::String
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
struct MeasurementData{T1<:Array{<:AbstractFloat, 1}, 
                       T2<:Array{<:String, 1}, 
                       T3<:Array{<:Symbol, 1}}
                    
    yObsNotTransformed::T1
    yObsTransformed::T1
    tObs::T1
    observebleID::T2
    conditionId::T2  # Sum of pre-eq + simulation-cond id 
    sdParams::T2
    transformData::T3 # Only done once 
    obsParam::T2
end


"""
    SimulationInfo

    Struct storing simulation (forward ODE-solution) information. Specifcially 
    stores the experimental ID:s from the experimentalCondition - PeTab file;
    firstExpIds (preequilibration ID:s), the shiftExpIds (postequilibration), and
    simulateSS (whether or not to simulate ODE-model to steady state). Further 
    stores a solArray with the ODE solution where conditionIdSol of the ID for 
    each forward solution

    See also: [`getSimulationInfo`]
"""
struct SimulationInfo{T1<:Array{<:String, 1}, 
                      T2<:Bool,
                      T3<:Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1}}
    firstExpIds::T1
    shiftExpIds::T1
    conditionIdSol::T1
    simulateSS::T2
    solArray::T3
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
    correctly to a vector of SD-vals). Further stores the name of each parameter.

    See also: [`getIndicesParam`]
"""
struct ParameterIndices{T1<:Array{<:Integer, 1}, 
                        T2<:Array{<:String, 1}}

    iDynParam::T1
    iObsParam::T1
    iSdParam::T1
    namesDynParam::T2
    namesObsParam::T2
    namesSdParam::T2
    namesParamEst::T2
end