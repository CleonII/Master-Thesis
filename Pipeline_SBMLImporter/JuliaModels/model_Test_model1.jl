# Model name: model_Test_model1
# Number of parameters: 5
# Number of species: 2
function getODEModel_model_Test_model1()

    ### Define constant parameters
    default = 1.0

    ### Define independent and dependent variables
    ModelingToolkit.@variables t sebastian(t) damiano(t)

    ### Define variable parameters

    ### Define dummy variable
    ModelingToolkit.@variables dummyVariable(t)

    ### Define parameters
    ModelingToolkit.@parameters alpha gamma delta beta alpha_scale

    ### Define an operator for the differentiation w.r.t. time
    D = Differential(t)

    ### Events ###

    ### Derivatives ###
    eqs = [
    D(sebastian) ~ +1.0 * ( 1 /default ) * (alpha*alpha_scale * sebastian + beta * damiano),
    D(damiano) ~ +1.0 * ( 1 /default ) * (gamma * sebastian + delta * damiano),
    D(dummyVariable) ~ +alpha
    ]

    @named sys = ODESystem(eqs)

    ### Initial species concentrations ###
    initialSpeciesValues = [
    sebastian => 8,
    damiano => 4,
    dummyVariable => 0.0]

    ### True parameter values ###
    trueParameterValues = [
    alpha => 5.0,
    gamma => 3.0,
    delta => 5.0,
    beta => 3.0, 
    alpha_scale => 1.0]

    return sys, initialSpeciesValues, trueParameterValues

end
