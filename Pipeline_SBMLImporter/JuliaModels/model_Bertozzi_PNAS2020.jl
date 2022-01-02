# Model name: model_Bertozzi_PNAS2020
# Number of parameters: 21
# Number of species: 3

### Define independent and dependent variables
@variables t Infected(t) Recovered(t) Susceptible(t)

### Define variable parameters
@variables Ro(t)

### Define dummy variable
@variables dummyVariable(t)

### Define parameters
@parameters Lockdown_NY_end Pop_CA Io_CA Pop_NY Io_NY Trigger_NY Lockdown_CA_start gamma_NY Trigger_Lockdown Ro_NY Lockdown_CA_end gamma_CA Ro_CA Lockdown_NY_start

### Define constants
@parameters USA___CA__NY

### Define an operator for the differentiation w.r.t. time
D = Differential(t)

### Events ###
continuous_events = [
[t ~ Trigger_NY * Lockdown_NY_start + (1 - Trigger_NY) * Lockdown_CA_start] => [Ro ~ Trigger_NY * (0.5 * Ro_NY) + (1 - Trigger_NY) * (0.5 * Ro_CA)], [t ~ Trigger_NY * Lockdown_NY_end + (1 - Trigger_NY) * Lockdown_CA_end] => [Ro ~ Trigger_NY * Ro_NY + (1 - Trigger_NY) * Ro_CA], [Trigger_Lockdown ~ 0.5] => [Ro ~ Trigger_NY * (0.5 * Ro_NY) + (1 - Trigger_NY) * (0.5 * Ro_CA)]
]

### Derivatives ###
eqs = [
D(Infected) ~ +1.0 * (USA___CA__NY * ((Trigger_NY * gamma_NY + (1 - Trigger_NY) * gamma_CA) * Ro * Infected * Susceptible))-1.0 * (USA___CA__NY * (Trigger_NY * gamma_NY + (1 - Trigger_NY) * gamma_CA) * Infected),
D(Recovered) ~ +1.0 * (USA___CA__NY * (Trigger_NY * gamma_NY + (1 - Trigger_NY) * gamma_CA) * Infected),
D(Susceptible) ~ -1.0 * (USA___CA__NY * ((Trigger_NY * gamma_NY + (1 - Trigger_NY) * gamma_CA) * Ro * Infected * Susceptible)),
D(Ro) ~ 0,
D(dummyVariable) ~ +Ro_NY+Lockdown_NY_end+Lockdown_CA_end+Trigger_NY+Ro_CA+Lockdown_CA_start+Lockdown_NY_start+Trigger_Lockdown
]

@named sys = ODESystem(eqs, t, continuous_events = continuous_events)

### Initial species concentrations ###
initialSpeciesValues = [
Infected => (Trigger_NY * Io_NY + (1 - Trigger_NY) * Io_CA) / (Trigger_NY * Pop_NY + (1 - Trigger_NY) * Pop_CA),
Recovered => 0.0,
Susceptible => ((Trigger_NY * Pop_NY + (1 - Trigger_NY) * Pop_CA) - Infected) / (Trigger_NY * Pop_NY + (1 - Trigger_NY) * Pop_CA),
Ro => 2.7,
dummyVariable => 0.0]

### True parameter values ###
trueParameterValues = [
Lockdown_NY_end => 105.0,
Pop_CA => 3.956e7,
Io_CA => 0.1,
Pop_NY => 1.954e7,
Io_NY => 0.005,
Trigger_NY => 0.0,
Lockdown_CA_start => 66.0,
gamma_NY => 0.1,
Trigger_Lockdown => 0.0,
Ro_NY => 4.1,
Lockdown_CA_end => 105.0,
gamma_CA => 0.12,
Ro_CA => 2.7,
Lockdown_NY_start => 69.0]

trueConstantsValues = [
USA___CA__NY => 1.0]
