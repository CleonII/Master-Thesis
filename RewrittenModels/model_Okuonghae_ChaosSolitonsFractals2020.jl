# Model name: model_Okuonghae_ChaosSolitonsFractals2020
# Number of parameters: 17
# Number of species: 9

### Define independent and dependent variables
@variables t detected_cumulative(t) symptomatic(t) asymptomatic(t) detected(t) exposed(t) deceased(t) recovered(t) susceptible(t)

### Define variable parameters

### Define dummy variable

### Define parameters
@parameters psi nu gamma_i asymptomatic_start delta symptomatic_start sigma theta alpha d_0 transmission_rate_effective eps exposed_start gamma_0 d_D gamma_a

### Define constants
@parameters Lagos

### Define an operator for the differentiation w.r.t. time
D = Differential(t)

### Events ###

### Derivatives ###
eqs = [
D(detected_cumulative) ~ +1.0 * (Lagos * psi * symptomatic)+1.0 * (Lagos * theta * asymptomatic),
D(symptomatic) ~ +1.0 * (Lagos * (sigma * (1 - nu) * exposed))-1.0 * (Lagos * psi * symptomatic)-1.0 * (Lagos * gamma_0 * symptomatic)-1.0 * (Lagos * d_0 * symptomatic),
D(asymptomatic) ~ +1.0 * (Lagos * (sigma * nu * exposed))-1.0 * (Lagos * theta * asymptomatic)-1.0 * (Lagos * gamma_a * asymptomatic),
D(detected) ~ +1.0 * (Lagos * psi * symptomatic)+1.0 * (Lagos * theta * asymptomatic)-1.0 * (Lagos * gamma_i * detected)-1.0 * (Lagos * d_D * detected),
D(exposed) ~ +1.0 * (Lagos * (transmission_rate_effective * (1 - eps) * (1 - delta) * (alpha * asymptomatic + symptomatic) / ((asymptomatic + detected + exposed + recovered + susceptible + symptomatic) - detected)) * susceptible)-1.0 * (Lagos * (sigma * (1 - nu) * exposed))-1.0 * (Lagos * (sigma * nu * exposed)),
D(deceased) ~ +1.0 * (Lagos * d_0 * symptomatic)+1.0 * (Lagos * d_D * detected),
D(recovered) ~ +1.0 * (Lagos * gamma_0 * symptomatic)+1.0 * (Lagos * gamma_a * asymptomatic)+1.0 * (Lagos * gamma_i * detected),
D(susceptible) ~ -1.0 * (Lagos * (transmission_rate_effective * (1 - eps) * (1 - delta) * (alpha * asymptomatic + symptomatic) / ((asymptomatic + detected + exposed + recovered + susceptible + symptomatic) - detected)) * susceptible)]

@named sys = ODESystem(eqs)

### Initial species concentrations ###
initialSpeciesValues = [
detected_cumulative => 1.0,
symptomatic => symptomatic_start,
asymptomatic => asymptomatic_start,
detected => 1.0,
exposed => exposed_start,
deceased => 0.0,
recovered => 0.0,
susceptible => 1.4367982e7]

### True parameter values ###
trueParameterValues = [
psi => 0.0135,
nu => 0.5,
gamma_i => 1 / 15,
asymptomatic_start => 188.0,
delta => 0.0,
symptomatic_start => 212.0,
sigma => 1 / 5.2,
theta => 1.8999 * (10)^(-12),
alpha => 0.5,
d_0 => 0.015,
transmission_rate_effective => 0.4236,
eps => 0.8,
exposed_start => 441.0,
gamma_0 => 0.13978,
d_D => 0.015,
gamma_a => 0.13978]

trueConstantsValues = [
Lagos => 1.0]