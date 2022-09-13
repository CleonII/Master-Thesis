# Model name: Test_model1
# Number of parameters: 4
# Number of species: 2
function getODEModel_Test_model1()

    ### Define independent and dependent variables
    ModelingToolkit.@variables t sebastian(t) damiano(t)

    ### Define variable parameters

    ### Define dummy variable
    ModelingToolkit.@variables dummyVariable(t)

    ### Define parameters
    ModelingToolkit.@parameters default alpha gamma delta beta

    ### Define an operator for the differentiation w.r.t. time
    D = Differential(t)

    ### Events ###

    ### Derivatives ###
    eqs = [
    D(sebastian) ~ +1.0 * ( 1 /default ) * (alpha * sebastian + beta * damiano),
    D(damiano) ~ +1.0 * ( 1 /default ) * (gamma * sebastian + delta * damiano),
    D(dummyVariable) ~ +default
    ]

    @named sys = ODESystem(eqs)

    ### Initial species concentrations ###
    initialSpeciesValues = [
    sebastian => 8,
    damiano => 4,
    dummyVariable => 0.0]

    ### True parameter values ###
    trueParameterValues = [
    default => 1.0,
    alpha => 5.0,
    gamma => 3.0,
    delta => 5.0,
    beta => 3.0]

    return sys, initialSpeciesValues, trueParameterValues

end
