# Model name: model_Oliveira_NatCommun2021
# Number of parameters: 33
# Number of species: 9

### Define independent and dependent variables
@variables t Hospital(t) Symptomatic(t) Cumulative_cases(t) Asymptomatic(t) Exposed(t) ICU(t) Recovered(t) Deaths(t) Susceptible(t)

### Define variable parameters

### Define dummy variable

### Define parameters
@parameters asymptomatic_init_concentration beta_2_multiplier t_2 gamma_u exposed_init_concentration omega_u kappa h_hosp_rate xi delta_ t_1 beta_0 symptomatic_init_concentration mu_u omega_h mu_h beta_2 beta_1 population p_symp_rate gamma_s gamma_h gamma_a

### Define constants
@parameters Interior

### Define an operator for the differentiation w.r.t. time
D = Differential(t)

### Function definitions ###
exposed_to_asymptomatic(kappa, p_symp_rate) = kappa * (1 - p_symp_rate)
@register exposed_to_asymptomatic(kappa, p_symp_rate)
symptomatic_to_icu(h_hosp_rate, xi, gamma_s) = h_hosp_rate * (1 - xi) * gamma_s
@register symptomatic_to_icu(h_hosp_rate, xi, gamma_s)
symptomatic_to_hospitalized(h_hosp_rate, xi, gamma_s) = h_hosp_rate * xi * gamma_s
@register symptomatic_to_hospitalized(h_hosp_rate, xi, gamma_s)
hospitalized_to_deaths(omega_h, mu_h, gamma_h) = (1 - omega_h) * mu_h * gamma_h
@register hospitalized_to_deaths(omega_h, mu_h, gamma_h)
hospitalized_to_recovered(omega_h, mu_h, gamma_h) = (1 - omega_h) * (1 - mu_h) * gamma_h
@register hospitalized_to_recovered(omega_h, mu_h, gamma_h)
icu_to_deaths(omega_u, mu_u, gamma_u) = (1 - omega_u) * mu_u * gamma_u
@register icu_to_deaths(omega_u, mu_u, gamma_u)
hospitalized_to_icu(omega_h, gamma_h) = omega_h * gamma_h
@register hospitalized_to_icu(omega_h, gamma_h)
symptomatic_to_recovered(h_hosp_rate, gamma_s) = (1 - h_hosp_rate) * gamma_s
@register symptomatic_to_recovered(h_hosp_rate, gamma_s)
exposed_to_symptomatic(kappa, p_symp_rate) = kappa * p_symp_rate
@register exposed_to_symptomatic(kappa, p_symp_rate)
beta(beta_2, beta_2_multiplier, beta_1, t, t_2, t_1, beta_0) = ((beta_2 * beta_2_multiplier + beta_1 * (1 - beta_2_multiplier) + (t < t_2) * ((beta_1) - (beta_2 * beta_2_multiplier + beta_1 * (1 - beta_2_multiplier)))) + (t < t_1) * ((beta_0) - ((beta_2 * beta_2_multiplier + beta_1 * (1 - beta_2_multiplier) + (t < t_2) * ((beta_1) - (beta_2 * beta_2_multiplier + beta_1 * (1 - beta_2_multiplier)))))))
@register beta(beta_2, beta_2_multiplier, beta_1, t, t_2, t_1, beta_0)

### Events ###

### Derivatives ###
eqs = [
D(Hospital) ~ +1.0 * (Interior * symptomatic_to_hospitalized(h_hosp_rate, xi, gamma_s) * Symptomatic)-1.0 * (Interior * hospitalized_to_recovered(omega_h, mu_h, gamma_h) * Hospital)-1.0 * (Interior * hospitalized_to_icu(omega_h, gamma_h) * Hospital)+1.0 * (Interior * (gamma_u * ICU * (1 - mu_u + omega_u * mu_u)))-1.0 * (Interior * hospitalized_to_deaths(omega_h, mu_h, gamma_h) * Hospital),
D(Symptomatic) ~ +1.0 * (Interior * exposed_to_symptomatic(kappa, p_symp_rate) * Exposed)-1.0 * (Interior * symptomatic_to_recovered(h_hosp_rate, gamma_s) * Symptomatic)-1.0 * (Interior * symptomatic_to_hospitalized(h_hosp_rate, xi, gamma_s) * Symptomatic)-1.0 * (Interior * symptomatic_to_icu(h_hosp_rate, xi, gamma_s) * Symptomatic),
D(Cumulative_cases) ~ exposed_to_symptomatic(kappa, p_symp_rate) * Exposed * Interior,
D(Asymptomatic) ~ +1.0 * (Interior * exposed_to_asymptomatic(kappa, p_symp_rate) * Exposed)-1.0 * (Interior * gamma_a * Asymptomatic),
D(Exposed) ~ +1.0 * (Interior * (beta(beta_2, beta_2_multiplier, beta_1, t, t_2, t_1, beta_0) * Susceptible * (Symptomatic + delta_ * Asymptomatic) / population))-1.0 * (Interior * exposed_to_asymptomatic(kappa, p_symp_rate) * Exposed)-1.0 * (Interior * exposed_to_symptomatic(kappa, p_symp_rate) * Exposed),
D(ICU) ~ +1.0 * (Interior * symptomatic_to_icu(h_hosp_rate, xi, gamma_s) * Symptomatic)+1.0 * (Interior * hospitalized_to_icu(omega_h, gamma_h) * Hospital)-1.0 * (Interior * (gamma_u * ICU * (1 - mu_u + omega_u * mu_u)))-1.0 * (Interior * icu_to_deaths(omega_u, mu_u, gamma_u) * ICU),
D(Recovered) ~ +1.0 * (Interior * gamma_a * Asymptomatic)+1.0 * (Interior * symptomatic_to_recovered(h_hosp_rate, gamma_s) * Symptomatic)+1.0 * (Interior * hospitalized_to_recovered(omega_h, mu_h, gamma_h) * Hospital),
D(Deaths) ~ +1.0 * (Interior * hospitalized_to_deaths(omega_h, mu_h, gamma_h) * Hospital)+1.0 * (Interior * icu_to_deaths(omega_u, mu_u, gamma_u) * ICU),
D(Susceptible) ~ -1.0 * (Interior * (beta(beta_2, beta_2_multiplier, beta_1, t, t_2, t_1, beta_0) * Susceptible * (Symptomatic + delta_ * Asymptomatic) / population))]

@named sys = ODESystem(eqs)

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
Susceptible => population * (1 - asymptomatic_init_concentration - exposed_init_concentration - symptomatic_init_concentration)]

### True parameter values ###
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
mu_u => 0.4,
omega_h => 0.14,
mu_h => 0.15,
beta_2 => 0.66,
beta_1 => 0.560641488770886,
population => 1.2098867e7,
p_symp_rate => 0.2,
gamma_s => 0.25,
gamma_h => 0.131351606436,
gamma_a => 0.285714]

trueConstantsValues = [
Interior => 1.0]
