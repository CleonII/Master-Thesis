# Model name: model_Bruno_JExpBot2016
# Number of parameters: 17
# Number of species: 7

### Define independent and dependent variables
@variables t b10(t) bio(t) ohbio(t) zea(t) bcry(t) ohb10(t) bcar(t)

### Define variable parameters

### Define dummy variable

### Define parameters
@parameters kc2_multiplier init_zea kc4_multiplier k5_multiplier kc1_multiplier init_b10 init_bcry kb1_multiplier kb2_multiplier kc1 kc4 init_ohb10 init_bcar kc2 kb2 k5 kb1

### Define constants
@parameters cyt

### Define an operator for the differentiation w.r.t. time
D = Differential(t)

### Events ###

### Derivatives ###
eqs = [
D(b10) ~ +1.0 * (cyt * bcar * kb1 * kb1_multiplier)-1.0 * (cyt * b10 * kb2 * kb2_multiplier)+1.0 * (cyt * bcry * kc1 * kc1_multiplier),
D(bio) ~ +1.0 * (cyt * bcar * kb1 * kb1_multiplier)+1.0 * (cyt * b10 * kb2 * kb2_multiplier)+1.0 * (cyt * bcry * kc2 * kc2_multiplier),
D(ohbio) ~ +1.0 * (cyt * bcry * kc1 * kc1_multiplier)+1.0 * (cyt * kc4 * kc4_multiplier * ohb10)+1.0 * (cyt * k5 * k5_multiplier * zea),
D(zea) ~ -1.0 * (cyt * k5 * k5_multiplier * zea),
D(bcry) ~ -1.0 * (cyt * bcry * kc1 * kc1_multiplier)-1.0 * (cyt * bcry * kc2 * kc2_multiplier),
D(ohb10) ~ +1.0 * (cyt * bcry * kc2 * kc2_multiplier)-1.0 * (cyt * kc4 * kc4_multiplier * ohb10)+1.0 * (cyt * k5 * k5_multiplier * zea),
D(bcar) ~ -1.0 * (cyt * bcar * kb1 * kb1_multiplier)]

@named sys = ODESystem(eqs)

### Initial species concentrations ###
initialSpeciesValues = [
b10 => init_b10,
bio => 0.0,
ohbio => 0.0,
zea => init_zea,
bcry => init_bcry,
ohb10 => init_ohb10,
bcar => init_bcar]

### True parameter values ###
trueParameterValues = [
kc2_multiplier => 0.521817105884857,
init_zea => 3.25673310603827,
kc4_multiplier => 0.521817105884857,
k5_multiplier => 0.521817105884857,
kc1_multiplier => 0.521817105884857,
init_b10 => 3.25673310603827,
init_bcry => 3.25673310603827,
kb1_multiplier => 0.521817105884857,
kb2_multiplier => 0.521817105884857,
kc1 => 0.00162720870139293,
kc4 => 0.00617966905040397,
init_ohb10 => 3.25673310603827,
init_bcar => 3.25673310603827,
kc2 => 0.00673060855967958,
kb2 => 0.00537146249584483,
k5 => 0.00305007086506138,
kb1 => 0.0164169857330715]

trueConstantsValues = [
cyt => 1.0]
