# Model name: model_Beer_MolBioSystems2014
# Number of parameters: 8
# Number of species: 4

### Define independent and dependent variables
@variables t Glu(t) cGlu(t) Ind(t) Bac(t) tau(t) lag(t)

### Define parameters
@parameters kdegi Bacmax ksyn kdim init_Bac beta

### Define constants
@parameters medium

### Define an operator for the differentiation w.r.t. time
D = Differential(t)

### Events ###
continuous_events = [
[t ~ tau] => [lag ~ 1]
]

### Function definitions ###

### Derivatives ###
eqs = [
D(Glu) ~ +1.0 * (medium * -Bac * Glu * ksyn),
D(cGlu) ~ +1.0 * (medium * (Bac * Glu * ksyn - (cGlu)^(2) * kdim)),
D(Ind) ~ +1.0 * (medium * ((cGlu)^(2) * kdim - Ind * kdegi)),
D(Bac) ~ +1.0 * (medium * (Bac * beta * lag * (Bacmax + -Bac) / Bacmax)),
D(tau) ~ 0,
D(lag) ~ 0]

@named sys = ODESystem(eqs, t, continuous_events = continuous_events)

### Initial species concentrations ###
initialSpeciesValues = [
Glu => 10.0,
cGlu => 0.0,
Ind => 0.0,
Bac => init_Bac,
tau => 1.0,
lag => 0.0]

### True parameter values ###
trueParameterValues = [
kdegi => 1.0,
Bacmax => 1.0,
ksyn => 1.0,
kdim => 1.0,
init_Bac => 0.0147007946993721,
beta => 1.0]

trueConstantsValues = [
medium => 1.0]
