# Model name: model_Bertozzi_PNAS2020
# Number of parameters: 21
# Number of species: 3
function getODEModel_model_Bertozzi_PNAS2020()

    ### Define independent and dependent variables
    ModelingToolkit.@variables t Infected(t) Recovered(t) Susceptible(t)

    ### Store dependent variables in array for ODESystem command
    stateArray = [Infected, Recovered, Susceptible]

    ### Define variable parameters

    ### Define potential algebraic variables
    ModelingToolkit.@variables Ro(t)

    ### Define parameters
    ModelingToolkit.@parameters Lockdown_NY_end Pop_CA Ro_bool1 Ro_bool2 Io_CA Ro_bool4 Pop_NY Io_NY Trigger_NY USA___CA__NY Lockdown_CA_start gamma_NY Trigger_Lockdown Ro_NY Lockdown_CA_end Ro_bool5 gamma_CA Ro_CA Lockdown_NY_start

    ### Store parameters in array for ODESystem command
    parameterArray = [Lockdown_NY_end, Pop_CA, Ro_bool1, Ro_bool2, Io_CA, Ro_bool4, Pop_NY, Io_NY, Trigger_NY, USA___CA__NY, Lockdown_CA_start, gamma_NY, Trigger_Lockdown, Ro_NY, Lockdown_CA_end, Ro_bool5, gamma_CA, Ro_CA, Lockdown_NY_start]

    ### Define an operator for the differentiation w.r.t. time
    D = Differential(t)

    ### Continious events ###

    ### Discrete events ###

    ### Derivatives ###
    eqs = [
    D(Infected) ~ +1.0 * ( 1 /USA___CA__NY ) * (USA___CA__NY * ((Trigger_NY * gamma_NY + (1 - Trigger_NY) * gamma_CA) * Ro * Infected * Susceptible))-1.0 * ( 1 /USA___CA__NY ) * (USA___CA__NY * (Trigger_NY * gamma_NY + (1 - Trigger_NY) * gamma_CA) * Infected),
    D(Recovered) ~ +1.0 * ( 1 /USA___CA__NY ) * (USA___CA__NY * (Trigger_NY * gamma_NY + (1 - Trigger_NY) * gamma_CA) * Infected),
    D(Susceptible) ~ -1.0 * ( 1 /USA___CA__NY ) * (USA___CA__NY * ((Trigger_NY * gamma_NY + (1 - Trigger_NY) * gamma_CA) * Ro * Infected * Susceptible)),
    Ro ~ ((((1 - Ro_bool1)*( 0.0) + Ro_bool1*( 1.0))) * (((1 - Ro_bool2)*( 1.0) + Ro_bool2*( 0.0)))) * (ifelse(Trigger_Lockdown >= 0.5, 1.0, 0.0)) * (Trigger_NY * (0.5 * Ro_NY) + (1 - Trigger_NY) * (0.5 * Ro_CA)) + (1 - ((((1 - Ro_bool1)*( 0.0) + Ro_bool1*( 1.0))) * (((1 - Ro_bool2)*( 1.0) + Ro_bool2*( 0.0)))) * (ifelse(Trigger_Lockdown >= 0.5, 1.0, 0.0))) * (Trigger_NY * Ro_NY + (1 - Trigger_NY) * Ro_CA)
    ]

    @named sys = ODESystem(eqs, t, stateArray, parameterArray)

    ### Initial species concentrations ###
    initialSpeciesValues = [
    Infected => (Trigger_NY * Io_NY + (1 - Trigger_NY) * Io_CA) / (Trigger_NY * Pop_NY + (1 - Trigger_NY) * Pop_CA),
    Recovered => 0.0,
    Susceptible => ((Trigger_NY * Pop_NY + (1 - Trigger_NY) * Pop_CA) - (Trigger_NY * Io_NY + (1 - Trigger_NY) * Io_CA) / (Trigger_NY * Pop_NY + (1 - Trigger_NY) * Pop_CA)) / (Trigger_NY * Pop_NY + (1 - Trigger_NY) * Pop_CA)
    ]

    ### SBML file parameter values ###
    trueParameterValues = [
    Lockdown_NY_end => 105.0,
    Pop_CA => 3.956e7,
    Ro_bool1 => 0.0,
    Ro_bool2 => 0.0,
    Io_CA => 0.1,
    Ro_bool4 => 0.0,
    Pop_NY => 1.954e7,
    Io_NY => 0.005,
    Trigger_NY => 0.0,
    USA___CA__NY => 1.0,
    Lockdown_CA_start => 66.0,
    gamma_NY => 0.1,
    Trigger_Lockdown => 0.0,
    Ro_NY => 4.1,
    Lockdown_CA_end => 105.0,
    Ro_bool5 => 0.0,
    gamma_CA => 0.12,
    Ro_CA => 2.7,
    Lockdown_NY_start => 69.0
    ]

    return sys, initialSpeciesValues, trueParameterValues

end
