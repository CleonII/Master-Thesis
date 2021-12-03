# Model name: model_Borghans_BiophysChem1997
# Number of parameters: 20
# Number of species: 3

### Define independent and dependent variables
@variables t A_state(t) Y_state(t) Z_state(t)

### Define variable parameters

### Define dummy variable

### Define parameters
@parameters v0 Ky Vm3 K2 Kz v1 Vm2 beta_par init_Y_state n_par K_par Kd epsilon_par Kp Kf init_A_state Vd init_Z_state Vp Ka

### Define constants
@parameters extracellular cytosol intravesicular

### Define an operator for the differentiation w.r.t. time
D = Differential(t)

### Function definitions ###

### Events ###

### Derivatives ###
eqs = [
D(A_state) ~ +1.0 * (cytosol * Vp * beta_par)-1.0 * (cytosol * ((A_state)^(2) * Vd * (Z_state)^(n_par) / (((Kd)^(n_par) + (Z_state)^(n_par)) * ((A_state)^(2) + (Kp)^(2)))))-1.0 * (cytosol * A_state * epsilon_par),
D(Y_state) ~ +1.0 * (cytosol * (Vm2 * (Z_state)^(2) / ((K2)^(2) + (Z_state)^(2))))-1.0 * (intravesicular * ((A_state)^(4) * Vm3 * (Y_state)^(2) * (Z_state)^(4) / (((A_state)^(4) + (Ka)^(4)) * ((Ky)^(2) + (Y_state)^(2)) * ((Kz)^(4) + (Z_state)^(4)))))-1.0 * (intravesicular * Kf * Y_state),
D(Z_state) ~ +1.0 * (cytosol * (v0 + beta_par * v1))-1.0 * (cytosol * (Vm2 * (Z_state)^(2) / ((K2)^(2) + (Z_state)^(2))))+1.0 * (intravesicular * ((A_state)^(4) * Vm3 * (Y_state)^(2) * (Z_state)^(4) / (((A_state)^(4) + (Ka)^(4)) * ((Ky)^(2) + (Y_state)^(2)) * ((Kz)^(4) + (Z_state)^(4)))))+1.0 * (intravesicular * Kf * Y_state)-1.0 * (cytosol * K_par * Z_state)]

@named sys = ODESystem(eqs)

### Initial species concentrations ###
initialSpeciesValues = [
A_state => init_A_state,
Y_state => init_Y_state,
Z_state => init_Z_state]

### True parameter values ###
trueParameterValues = [
v0 => 2.31778715779187,
Ky => 0.200364133028272,
Vm3 => 22.7025248825351,
K2 => 0.0999485580156493,
Kz => 0.303353547057472,
v1 => 1.00488755696677,
Vm2 => 7.45712445492225,
beta_par => 1.12395230256787,
init_Y_state => 0.999348084438687,
n_par => 4.1025144127497,
K_par => 11.4120804594265,
Kd => 0.392497512107474,
epsilon_par => 0.163232195306189,
Kp => 0.996567464196473,
Kf => 1.13794968361236,
init_A_state => 0.99999999999996,
Vd => 92.8739300577965,
init_Z_state => 0.0879205244255038,
Vp => 2.75685784345759,
Ka => 0.197593940310187]

trueConstantsValues = [
extracellular => 1.0,
cytosol => 1.0,
intravesicular => 1.0]
