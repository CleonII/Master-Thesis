# Model name: model_Test_model2
# Number of parameters: 2
# Number of species: 2
function getODEModel_model_Test_model2()

    ### Define constant parameters
    default = 1.0

    ### Define independent and dependent variables
    ModelingToolkit.@variables t sebastian(t) damiano(t)

    ### Define variable parameters

    ### Define dummy variable

    ### Define parameters
    ModelingToolkit.@parameters alpha beta

    ### Define an operator for the differentiation w.r.t. time
    D = Differential(t)

    ### Events ###

    ### Derivatives ###
    eqs = [
    D(sebastian) ~ +1.0 * ( 1 /default ) * (alpha * sebastian),
    D(damiano) ~ +1.0 * ( 1 /default ) * (beta * damiano)    ]

    @named sys = ODESystem(eqs)

    ### Initial species concentrations ###
    initialSpeciesValues = [
    sebastian => 8,
    damiano => 4]

    ### True parameter values ###
    trueParameterValues = [
    alpha => 5.0,
    beta => 3.0]

    return sys, initialSpeciesValues, trueParameterValues

end
