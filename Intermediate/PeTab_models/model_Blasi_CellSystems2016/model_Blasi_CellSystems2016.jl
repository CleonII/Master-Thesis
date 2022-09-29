# Model name: model_Blasi_CellSystems2016
# Number of parameters: 34
# Number of species: 16
function getODEModel_model_Blasi_CellSystems2016()

    ### Define independent and dependent variables
    ModelingToolkit.@variables t x_k5k12k16(t) x_k8(t) x_k16(t) x_0ac(t) x_k12(t) x_k5k8(t) x_k5k12(t) x_k12k16(t) x_k8k12k16(t) x_k5(t) x_k5k16(t) x_k5k8k12(t) x_k8k12(t) x_4ac(t) x_k8k16(t) x_k5k8k16(t)

    ### Define variable parameters

    ### Define dummy variable
    ModelingToolkit.@variables dummyVariable(t)

    ### Define parameters
    ModelingToolkit.@parameters a_k5_k5k12 a_k8 d default a_k12k16_k8k12k16 a_basal a_k5k12_k5k8k12 a_k8k12k16_4ac a_k12_k5k12 a_k16_k12k16

    ### Define an operator for the differentiation w.r.t. time
    D = Differential(t)

    ### Events ###

    ### Derivatives ###
    eqs = [
    D(x_k5k12k16) ~ +1.0 * ( 1 /default ) * ((a_basal) * x_k5k12 - d * x_k5k12k16)+1.0 * ( 1 /default ) * ((a_basal) * x_k5k16 - d * x_k5k12k16)+1.0 * ( 1 /default ) * ((a_basal) * x_k12k16 - d * x_k5k12k16)-1.0 * ( 1 /default ) * ((a_basal) * x_k5k12k16 - d * x_4ac),
    D(x_k8) ~ +1.0 * ( 1 /default ) * (a_k8 * x_0ac - d * x_k8)-1.0 * ( 1 /default ) * ((a_basal) * x_k8 - d * x_k5k8)-1.0 * ( 1 /default ) * ((a_basal) * x_k8 - d * x_k8k12)-1.0 * ( 1 /default ) * ((a_basal) * x_k8 - d * x_k8k16),
    D(x_k16) ~ +1.0 * ( 1 /default ) * ((a_basal) * x_0ac - d * x_k16)-1.0 * ( 1 /default ) * ((a_basal) * x_k16 - d * x_k5k16)-1.0 * ( 1 /default ) * ((a_basal) * x_k16 - d * x_k8k16)-1.0 * ( 1 /default ) * (a_k16_k12k16 * x_k16 - d * x_k12k16),
    D(x_0ac) ~ -1.0 * ( 1 /default ) * ((a_basal) * x_0ac - d * x_k5)-1.0 * ( 1 /default ) * (a_k8 * x_0ac - d * x_k8)-1.0 * ( 1 /default ) * ((a_basal) * x_0ac - d * x_k12)-1.0 * ( 1 /default ) * ((a_basal) * x_0ac - d * x_k16),
    D(x_k12) ~ +1.0 * ( 1 /default ) * ((a_basal) * x_0ac - d * x_k12)-1.0 * ( 1 /default ) * (a_k12_k5k12 * x_k12 - d * x_k5k12)-1.0 * ( 1 /default ) * ((a_basal) * x_k12 - d * x_k8k12)-1.0 * ( 1 /default ) * ((a_basal) * x_k12 - d * x_k12k16),
    D(x_k5k8) ~ +1.0 * ( 1 /default ) * ((a_basal) * x_k5 - d * x_k5k8)+1.0 * ( 1 /default ) * ((a_basal) * x_k8 - d * x_k5k8)-1.0 * ( 1 /default ) * ((a_basal) * x_k5k8 - d * x_k5k8k12)-1.0 * ( 1 /default ) * ((a_basal) * x_k5k8 - d * x_k5k8k16),
    D(x_k5k12) ~ +1.0 * ( 1 /default ) * (a_k5_k5k12 * x_k5 - d * x_k5k12)+1.0 * ( 1 /default ) * (a_k12_k5k12 * x_k12 - d * x_k5k12)-1.0 * ( 1 /default ) * (a_k5k12_k5k8k12 * x_k5k12 - d * x_k5k8k12)-1.0 * ( 1 /default ) * ((a_basal) * x_k5k12 - d * x_k5k12k16),
    D(x_k12k16) ~ +1.0 * ( 1 /default ) * ((a_basal) * x_k12 - d * x_k12k16)+1.0 * ( 1 /default ) * (a_k16_k12k16 * x_k16 - d * x_k12k16)-1.0 * ( 1 /default ) * ((a_basal) * x_k12k16 - d * x_k5k12k16)-1.0 * ( 1 /default ) * (a_k12k16_k8k12k16 * x_k12k16 - d * x_k8k12k16),
    D(x_k8k12k16) ~ +1.0 * ( 1 /default ) * ((a_basal) * x_k8k12 - d * x_k8k12k16)+1.0 * ( 1 /default ) * ((a_basal) * x_k8k16 - d * x_k8k12k16)+1.0 * ( 1 /default ) * (a_k12k16_k8k12k16 * x_k12k16 - d * x_k8k12k16)-1.0 * ( 1 /default ) * (a_k8k12k16_4ac * x_k8k12k16 - d * x_4ac),
    D(x_k5) ~ +1.0 * ( 1 /default ) * ((a_basal) * x_0ac - d * x_k5)-1.0 * ( 1 /default ) * ((a_basal) * x_k5 - d * x_k5k8)-1.0 * ( 1 /default ) * (a_k5_k5k12 * x_k5 - d * x_k5k12)-1.0 * ( 1 /default ) * ((a_basal) * x_k5 - d * x_k5k16),
    D(x_k5k16) ~ +1.0 * ( 1 /default ) * ((a_basal) * x_k5 - d * x_k5k16)+1.0 * ( 1 /default ) * ((a_basal) * x_k16 - d * x_k5k16)-1.0 * ( 1 /default ) * ((a_basal) * x_k5k16 - d * x_k5k8k16)-1.0 * ( 1 /default ) * ((a_basal) * x_k5k16 - d * x_k5k12k16),
    D(x_k5k8k12) ~ +1.0 * ( 1 /default ) * ((a_basal) * x_k5k8 - d * x_k5k8k12)+1.0 * ( 1 /default ) * (a_k5k12_k5k8k12 * x_k5k12 - d * x_k5k8k12)+1.0 * ( 1 /default ) * ((a_basal) * x_k8k12 - d * x_k5k8k12)-1.0 * ( 1 /default ) * ((a_basal) * x_k5k8k12 - d * x_4ac),
    D(x_k8k12) ~ +1.0 * ( 1 /default ) * ((a_basal) * x_k8 - d * x_k8k12)+1.0 * ( 1 /default ) * ((a_basal) * x_k12 - d * x_k8k12)-1.0 * ( 1 /default ) * ((a_basal) * x_k8k12 - d * x_k5k8k12)-1.0 * ( 1 /default ) * ((a_basal) * x_k8k12 - d * x_k8k12k16),
    D(x_4ac) ~ +1.0 * ( 1 /default ) * ((a_basal) * x_k5k8k12 - d * x_4ac)+1.0 * ( 1 /default ) * ((a_basal) * x_k5k8k16 - d * x_4ac)+1.0 * ( 1 /default ) * ((a_basal) * x_k5k12k16 - d * x_4ac)+1.0 * ( 1 /default ) * (a_k8k12k16_4ac * x_k8k12k16 - d * x_4ac),
    D(x_k8k16) ~ +1.0 * ( 1 /default ) * ((a_basal) * x_k8 - d * x_k8k16)+1.0 * ( 1 /default ) * ((a_basal) * x_k16 - d * x_k8k16)-1.0 * ( 1 /default ) * ((a_basal) * x_k8k16 - d * x_k5k8k16)-1.0 * ( 1 /default ) * ((a_basal) * x_k8k16 - d * x_k8k12k16),
    D(x_k5k8k16) ~ +1.0 * ( 1 /default ) * ((a_basal) * x_k5k8 - d * x_k5k8k16)+1.0 * ( 1 /default ) * ((a_basal) * x_k5k16 - d * x_k5k8k16)+1.0 * ( 1 /default ) * ((a_basal) * x_k8k16 - d * x_k5k8k16)-1.0 * ( 1 /default ) * ((a_basal) * x_k5k8k16 - d * x_4ac),
    D(dummyVariable) ~ 1e-60*( +default)
    ]

    @named sys = ODESystem(eqs)

    ### Initial species concentrations ###
    initialSpeciesValues = [
    x_k5k12k16 => 0.0,
    x_k8 => 0.0,
    x_k16 => 0.0,
    x_0ac => 1.0,
    x_k12 => 0.0,
    x_k5k8 => 0.0,
    x_k5k12 => 0.0,
    x_k12k16 => 0.0,
    x_k8k12k16 => 0.0,
    x_k5 => 0.0,
    x_k5k16 => 0.0,
    x_k5k8k12 => 0.0,
    x_k8k12 => 0.0,
    x_4ac => 0.0,
    x_k8k16 => 0.0,
    x_k5k8k16 => 0.0,
    dummyVariable => 0.0]

    ### True parameter values ###
    trueParameterValues = [
    a_k5_k5k12 => 2.062,
    a_k8 => 0.0273,
    d => 1.0,
    default => 1.0,
    a_k12k16_k8k12k16 => 2.206,
    a_basal => 0.067,
    a_k5k12_k5k8k12 => 0.325,
    a_k8k12k16_4ac => 3.5917,
    a_k12_k5k12 => 0.552,
    a_k16_k12k16 => 0.696]

    return sys, initialSpeciesValues, trueParameterValues

end
