# Model name: model_Perelson_Science1996
# Number of parameters: 5
# Number of species: 4

### Define independent and dependent variables
ModelingToolkit.@variables t Vni(t) V(t) Vin(t) Tstar(t)

### Define variable parameters

### Define dummy variable

### Define parameters
ModelingToolkit.@parameters c T0 K0 NN delta

### Define constants
ModelingToolkit.@parameters default

### Define an operator for the differentiation w.r.t. time
D = Differential(t)

### Events ###

### Derivatives ###
eqs = [
D(Vni) ~ -1.0 * (c * Vni)+1.0 * (delta * NN * Tstar),
D(V) ~ -1.0 * (c * Vin)-1.0 * (c * Vni)+1.0 * (delta * NN * Tstar),
D(Vin) ~ -1.0 * (c * Vin),
D(Tstar) ~ +1.0 * (K0 * T0 * Vin)-1.0 * (delta * Tstar)]

@named sys = ODESystem(eqs)

### Initial species concentrations ###
initialSpeciesValues = [
Vni => 0.0,
V => 1.86e6,
Vin => 1.86e6,
Tstar => 15061.32075]

### True parameter values ###
trueParameterValues = [
c => 2.06,
T0 => 11000.0,
K0 => 3.9e-7,
NN => 480.0,
delta => 0.53]

trueConstantsValues = [
default => 1.0]
