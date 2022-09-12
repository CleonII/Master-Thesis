# Model name: model_Rahman_MBS2016
# Number of parameters: 24
# Number of species: 7
function getODEModel_model_Rahman_MBS2016()

    ### Define constant parameters
    c1 = 1.0
    treated_normal_death_rate = 0.0408
    infected_normal_treatment_rate = 0.0
    infected_moderate_treatment_rate = 0.0
    susceptible_death_rate = 0.0288
    infected_normal_death_rate = 0.0888
    treated_weak_death_rate = 0.1752
    treated_moderate_death_rate = 0.0528
    recruitment_rate = 1.032672e6
    infected_moderate_death_rate = 0.1368
    infected_weak_death_rate = 0.3108

    ### Define independent and dependent variables
    ModelingToolkit.@variables t treated_weak(t) infected_moderate(t) treated_normal(t) infected_normal(t) infected_weak(t) treated_moderate(t) susceptible(t)

    ### Define variable parameters

    ### Define dummy variable

    ### Define parameters
    ModelingToolkit.@parameters infected_normal_transmission_rate_relative infected_weak_transmission_rate_relative infected_moderate_transmission_rate infected_normal_worsen_rate treated_moderate_improve_rate infected_weak_treatment_rate behavioural_change_rate infected_moderate_worsen_rate treated_weak_improve_rate

    ### Define an operator for the differentiation w.r.t. time
    D = Differential(t)

    ### Events ###

    ### Derivatives ###
    eqs = [
    D(treated_weak) ~ -1.0 * ( 1 /c1 ) * (treated_weak_improve_rate * treated_weak)+1.0 * ( 1 /c1 ) * (infected_weak_treatment_rate * infected_weak)-1.0 * ( 1 /c1 ) * (treated_weak_death_rate * treated_weak),
    D(infected_moderate) ~ +1.0 * ( 1 /c1 ) * (infected_normal_worsen_rate * infected_normal)-1.0 * ( 1 /c1 ) * (infected_moderate_worsen_rate * infected_moderate)-1.0 * ( 1 /c1 ) * (infected_moderate_treatment_rate * infected_moderate)-1.0 * ( 1 /c1 ) * (infected_moderate_death_rate * infected_moderate),
    D(treated_normal) ~ +1.0 * ( 1 /c1 ) * (treated_moderate_improve_rate * treated_moderate)+1.0 * ( 1 /c1 ) * (infected_normal_treatment_rate * infected_normal)-1.0 * ( 1 /c1 ) * (treated_normal_death_rate * treated_normal),
    D(infected_normal) ~ +1.0 * ( 1 /c1 ) * ((((infected_normal_transmission_rate_relative * infected_moderate_transmission_rate) * infected_normal + infected_moderate_transmission_rate * infected_moderate + (infected_weak_transmission_rate_relative * infected_moderate_transmission_rate) * infected_weak + (0.04 * infected_moderate_transmission_rate) * (treated_normal + treated_moderate + treated_weak)) / (susceptible + infected_normal + infected_moderate + infected_weak + treated_normal + treated_moderate + treated_weak) * exp(-1 * behavioural_change_rate * (infected_normal + infected_moderate + infected_weak + treated_normal + treated_moderate + treated_weak))) * susceptible)-1.0 * ( 1 /c1 ) * (infected_normal_worsen_rate * infected_normal)-1.0 * ( 1 /c1 ) * (infected_normal_treatment_rate * infected_normal)-1.0 * ( 1 /c1 ) * (infected_normal_death_rate * infected_normal),
    D(infected_weak) ~ +1.0 * ( 1 /c1 ) * (infected_moderate_worsen_rate * infected_moderate)-1.0 * ( 1 /c1 ) * (infected_weak_treatment_rate * infected_weak)-1.0 * ( 1 /c1 ) * (infected_weak_death_rate * infected_weak),
    D(treated_moderate) ~ +1.0 * ( 1 /c1 ) * (treated_weak_improve_rate * treated_weak)-1.0 * ( 1 /c1 ) * (treated_moderate_improve_rate * treated_moderate)+1.0 * ( 1 /c1 ) * (infected_moderate_treatment_rate * infected_moderate)-1.0 * ( 1 /c1 ) * (treated_moderate_death_rate * treated_moderate),
    D(susceptible) ~ +1.0 * ( 1 /c1 ) * (recruitment_rate)-1.0 * ( 1 /c1 ) * ((((infected_normal_transmission_rate_relative * infected_moderate_transmission_rate) * infected_normal + infected_moderate_transmission_rate * infected_moderate + (infected_weak_transmission_rate_relative * infected_moderate_transmission_rate) * infected_weak + (0.04 * infected_moderate_transmission_rate) * (treated_normal + treated_moderate + treated_weak)) / (susceptible + infected_normal + infected_moderate + infected_weak + treated_normal + treated_moderate + treated_weak) * exp(-1 * behavioural_change_rate * (infected_normal + infected_moderate + infected_weak + treated_normal + treated_moderate + treated_weak))) * susceptible)-1.0 * ( 1 /c1 ) * (susceptible_death_rate * susceptible)    ]

    @named sys = ODESystem(eqs)

    ### Initial species concentrations ###
    initialSpeciesValues = [
    treated_weak => 0.0,
    infected_moderate => 9000.0,
    treated_normal => 0.0,
    infected_normal => 16300.0,
    infected_weak => 11000.0,
    treated_moderate => 0.0,
    susceptible => 1.794e7]

    ### True parameter values ###
    trueParameterValues = [
    infected_normal_transmission_rate_relative => 12.57,
    infected_weak_transmission_rate_relative => 4.54,
    infected_moderate_transmission_rate => 0.082,
    infected_normal_worsen_rate => 0.33,
    treated_moderate_improve_rate => 0.57,
    infected_weak_treatment_rate => 0.11,
    behavioural_change_rate => 2.4744e-7,
    infected_moderate_worsen_rate => 0.34,
    treated_weak_improve_rate => 0.82]

    return sys, initialSpeciesValues, trueParameterValues

end
