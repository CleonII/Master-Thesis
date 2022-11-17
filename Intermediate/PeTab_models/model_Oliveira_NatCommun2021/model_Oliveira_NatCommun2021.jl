# Model name: model_Oliveira_NatCommun2021
# Number of parameters: 33
# Number of species: 9
function getODEModel_model_Oliveira_NatCommun2021()

    ### Define independent and dependent variables
    ModelingToolkit.@variables t Hospital(t) Symptomatic(t) Cumulative_cases(t) Asymptomatic(t) Exposed(t) ICU(t) Recovered(t) Deaths(t) Susceptible(t)

    ### Store dependent variables in array for ODESystem command
    stateArray = [Hospital, Symptomatic, Cumulative_cases, Asymptomatic, Exposed, ICU, Recovered, Deaths, Susceptible]

    ### Define variable parameters

    ### Define potential algebraic variables
    ModelingToolkit.@variables beta(t)

    ### Define parameters
    ModelingToolkit.@parameters asymptomatic_init_concentration beta_2_multiplier t_2 gamma_u exposed_init_concentration omega_u kappa h_hosp_rate xi delta_ t_1 beta_0 symptomatic_init_concentration Interior mu_u omega_h mu_h beta_2 beta_1 population p_symp_rate gamma_s gamma_h gamma_a

    ### Store parameters in array for ODESystem command
    parameterArray = [asymptomatic_init_concentration, beta_2_multiplier, t_2, gamma_u, exposed_init_concentration, omega_u, kappa, h_hosp_rate, xi, delta_, t_1, beta_0, symptomatic_init_concentration, Interior, mu_u, omega_h, mu_h, beta_2, beta_1, population, p_symp_rate, gamma_s, gamma_h, gamma_a]

    ### Define an operator for the differentiation w.r.t. time
    D = Differential(t)

    ### Continious events ###

    ### Discrete events ###

    ### Derivatives ###
    eqs = [
    D(Hospital) ~ +1.0 * ( 1 /Interior ) * (Interior * (h_hosp_rate * xi * gamma_s) * Symptomatic)-1.0 * ( 1 /Interior ) * (Interior * ((1 - omega_h) * (1 - mu_h) * gamma_h) * Hospital)-1.0 * ( 1 /Interior ) * (Interior * (omega_h * gamma_h) * Hospital)+1.0 * ( 1 /Interior ) * (Interior * (gamma_u * ICU * (1 - mu_u + omega_u * mu_u)))-1.0 * ( 1 /Interior ) * (Interior * ((1 - omega_h) * mu_h * gamma_h) * Hospital),
    D(Symptomatic) ~ +1.0 * ( 1 /Interior ) * (Interior * (kappa * p_symp_rate) * Exposed)-1.0 * ( 1 /Interior ) * (Interior * ((1 - h_hosp_rate) * gamma_s) * Symptomatic)-1.0 * ( 1 /Interior ) * (Interior * (h_hosp_rate * xi * gamma_s) * Symptomatic)-1.0 * ( 1 /Interior ) * (Interior * (h_hosp_rate * (1 - xi) * gamma_s) * Symptomatic),
    D(Cumulative_cases) ~ (kappa * p_symp_rate) * Exposed * Interior,
    D(Asymptomatic) ~ +1.0 * ( 1 /Interior ) * (Interior * (kappa * (1 - p_symp_rate)) * Exposed)-1.0 * ( 1 /Interior ) * (Interior * gamma_a * Asymptomatic),
    D(Exposed) ~ +1.0 * ( 1 /Interior ) * (Interior * (beta * Susceptible * (Symptomatic + delta_ * Asymptomatic) / population))-1.0 * ( 1 /Interior ) * (Interior * (kappa * (1 - p_symp_rate)) * Exposed)-1.0 * ( 1 /Interior ) * (Interior * (kappa * p_symp_rate) * Exposed),
    D(ICU) ~ +1.0 * ( 1 /Interior ) * (Interior * (h_hosp_rate * (1 - xi) * gamma_s) * Symptomatic)+1.0 * ( 1 /Interior ) * (Interior * (omega_h * gamma_h) * Hospital)-1.0 * ( 1 /Interior ) * (Interior * (gamma_u * ICU * (1 - mu_u + omega_u * mu_u)))-1.0 * ( 1 /Interior ) * (Interior * ((1 - omega_u) * mu_u * gamma_u) * ICU),
    D(Recovered) ~ +1.0 * ( 1 /Interior ) * (Interior * gamma_a * Asymptomatic)+1.0 * ( 1 /Interior ) * (Interior * ((1 - h_hosp_rate) * gamma_s) * Symptomatic)+1.0 * ( 1 /Interior ) * (Interior * ((1 - omega_h) * (1 - mu_h) * gamma_h) * Hospital),
    D(Deaths) ~ +1.0 * ( 1 /Interior ) * (Interior * ((1 - omega_h) * mu_h * gamma_h) * Hospital)+1.0 * ( 1 /Interior ) * (Interior * ((1 - omega_u) * mu_u * gamma_u) * ICU),
    D(Susceptible) ~ -1.0 * ( 1 /Interior ) * (Interior * (beta * Susceptible * (Symptomatic + delta_ * Asymptomatic) / population)),
    beta ~ ifelse(t < t_1, beta_0, ifelse(t < t_2, beta_1, beta_2 * beta_2_multiplier + beta_1 * (1 - beta_2_multiplier)))
    ]

    @named sys = ODESystem(eqs, t, stateArray, parameterArray)

    ### Initial species concentrations ###
    initialSpeciesValues = [
    Hospital => 0.0,
    Symptomatic => symptomatic_init_concentration * population,
    Cumulative_cases => symptomatic_init_concentration * population,
    Asymptomatic => asymptomatic_init_concentration * population,
    Exposed => exposed_init_concentration * population,
    ICU => 0.0,
    Recovered => 0.0,
    Deaths => 0.0,
    Susceptible => population * (1 - asymptomatic_init_concentration - exposed_init_concentration - symptomatic_init_concentration)
    ]

    ### SBML file parameter values ###
    trueParameterValues = [
    asymptomatic_init_concentration => 6.4828815498309e-7,
    beta_2_multiplier => 1.0,
    t_2 => 97.0,
    gamma_u => 0.1870201755885,
    exposed_init_concentration => 2.24673170025082e-6,
    omega_u => 0.29,
    kappa => 0.25,
    h_hosp_rate => 0.0,
    xi => 0.53,
    delta_ => 0.6951361054,
    t_1 => 23.4530530945,
    beta_0 => 0.807194090959563,
    symptomatic_init_concentration => 1.99259135560373e-6,
    Interior => 1.0,
    mu_u => 0.4,
    omega_h => 0.14,
    mu_h => 0.15,
    beta_2 => 0.66,
    beta_1 => 0.560641488770886,
    population => 1.2098867e7,
    p_symp_rate => 0.2,
    gamma_s => 0.25,
    gamma_h => 0.131351606436,
    gamma_a => 0.285714
    ]

    return sys, initialSpeciesValues, trueParameterValues

end
