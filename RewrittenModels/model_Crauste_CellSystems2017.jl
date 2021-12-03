# Model name: model_Crauste_CellSystems2017
# Number of parameters: 12
# Number of species: 5

### Define independent and dependent variables
@variables t Naive(t) Pathogen(t) LateEffector(t) EarlyEffector(t) Memory(t)

### Define variable parameters

### Define dummy variable

### Define parameters
@parameters mu_LL delta_NE mu_PE mu_P mu_PL delta_EL mu_EE mu_N rho_E delta_LM rho_P mu_LE

### Define constants
@parameters default

### Define an operator for the differentiation w.r.t. time
D = Differential(t)

### Function definitions ###

### Events ###

### Derivatives ###
eqs = [
D(Naive) ~ -1.0 * (Naive * mu_N)-1.0 * (Naive * Pathogen * delta_NE),
D(Pathogen) ~ +1.0 * ((Pathogen)^(2) * rho_P)-1.0 * (EarlyEffector * Pathogen * mu_PE)-1.0 * (LateEffector * Pathogen * mu_PL)-1.0 * (Pathogen * mu_P),
D(LateEffector) ~ +1.0 * (EarlyEffector * delta_EL)-1.0 * ((LateEffector)^(2) * mu_LL)-1.0 * (EarlyEffector * LateEffector * mu_LE)-1.0 * (LateEffector * delta_LM),
D(EarlyEffector) ~ +1.0 * (Naive * Pathogen * delta_NE)+1.0 * (EarlyEffector * Pathogen * rho_E)-1.0 * ((EarlyEffector)^(2) * mu_EE)-1.0 * (EarlyEffector * delta_EL),
D(Memory) ~ +1.0 * (LateEffector * delta_LM)]

@named sys = ODESystem(eqs)

### Initial species concentrations ###
initialSpeciesValues = [
Naive => 8090.0,
Pathogen => 1.0,
LateEffector => 0.0,
EarlyEffector => 0.0,
Memory => 0.0]

### True parameter values ###
trueParameterValues = [
mu_LL => 8.11520135326853e-6,
delta_NE => 0.0119307857579241,
mu_PE => 1.36571832778378e-10,
mu_P => 1.00000002976846e-5,
mu_PL => 3.6340308186265e-5,
delta_EL => 0.51794597529254,
mu_EE => 3.91359322673521e-5,
mu_N => 0.739907308603256,
rho_E => 0.507415703707752,
delta_LM => 0.0225806365892933,
rho_P => 0.126382288121756,
mu_LE => 1.00000000000005e-10]

trueConstantsValues = [
default => 1.0]
