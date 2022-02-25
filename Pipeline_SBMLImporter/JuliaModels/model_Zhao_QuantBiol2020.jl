# Model name: model_Zhao_QuantBiol2020
# Number of parameters: 37
# Number of species: 5

### Define independent and dependent variables
ModelingToolkit.@variables t Quarantined_Infected(t) Confirmed_Infected(t) Susceptible(t) Unquarantined_Infected(t)

### Define variable parameters

### Define dummy variable

### Define parameters
ModelingToolkit.@parameters Total_Pop_Hubei R_Stage_I_Wuhan gamma_2_Stage_III_Wuhan Total_Pop_China Trigger_Wuhan gamma_1_Stage_I_China gamma_1_Stage_I_Hubei sigma R_Stage_II_Wuhan R_Stage_II_China gamma_1_Stage_I_Wuhan gamma_2_Stage_II_Wuhan Trigger_China gamma_2_Stage_II_China R_Stage_III_Wuhan Trigger_Stage_I gamma_2_Stage_I_Wuhan R_Stage_I_China Trigger_Stage_III gamma_1_Stage_II_China gamma_1_Stage_III_Wuhan gamma_2_Stage_I_China Trigger_Stage_II R_Stage_II_Hubei Total_Pop_Wuhan gamma_1_Stage_II_Wuhan gamma_2_Stage_I_Hubei R_Stage_I_Hubei gamma_1_Stage_II_Hubei Trigger_Hubei gamma_2_Stage_II_Hubei

### Define constants
ModelingToolkit.@parameters China

### Define an operator for the differentiation w.r.t. time
D = Differential(t)

### Events ###

### Derivatives ###
eqs = [
D(Quarantined_Infected) ~ +1.0 * (China * (Trigger_Wuhan * (Trigger_Stage_I * gamma_1_Stage_I_Wuhan + Trigger_Stage_II * gamma_1_Stage_II_Wuhan + Trigger_Stage_III * gamma_1_Stage_III_Wuhan) + Trigger_Hubei * (Trigger_Stage_I * gamma_1_Stage_I_Hubei + Trigger_Stage_II * gamma_1_Stage_II_Hubei) + Trigger_China * (Trigger_Stage_I * gamma_1_Stage_I_China + Trigger_Stage_II * gamma_1_Stage_II_China)) * Unquarantined_Infected)-1.0 * (China * ((Trigger_Wuhan * (Trigger_Stage_I * gamma_2_Stage_I_Wuhan + Trigger_Stage_II * gamma_2_Stage_II_Wuhan + Trigger_Stage_III * gamma_2_Stage_III_Wuhan) + Trigger_Hubei * (Trigger_Stage_I * gamma_2_Stage_I_Hubei + Trigger_Stage_II * gamma_2_Stage_II_Hubei) + Trigger_China * (Trigger_Stage_I * gamma_2_Stage_I_China + Trigger_Stage_II * gamma_2_Stage_II_China)) + (1 - (Trigger_Wuhan * (Trigger_Stage_I * gamma_2_Stage_I_Wuhan + Trigger_Stage_II * gamma_2_Stage_II_Wuhan + Trigger_Stage_III * gamma_2_Stage_III_Wuhan) + Trigger_Hubei * (Trigger_Stage_I * gamma_2_Stage_I_Hubei + Trigger_Stage_II * gamma_2_Stage_II_Hubei) + Trigger_China * (Trigger_Stage_I * gamma_2_Stage_I_China + Trigger_Stage_II * gamma_2_Stage_II_China))) * sigma) * Quarantined_Infected),
D(Confirmed_Infected) ~ +1.0 * (China * ((Trigger_Wuhan * (Trigger_Stage_I * gamma_2_Stage_I_Wuhan + Trigger_Stage_II * gamma_2_Stage_II_Wuhan + Trigger_Stage_III * gamma_2_Stage_III_Wuhan) + Trigger_Hubei * (Trigger_Stage_I * gamma_2_Stage_I_Hubei + Trigger_Stage_II * gamma_2_Stage_II_Hubei) + Trigger_China * (Trigger_Stage_I * gamma_2_Stage_I_China + Trigger_Stage_II * gamma_2_Stage_II_China)) + (1 - (Trigger_Wuhan * (Trigger_Stage_I * gamma_2_Stage_I_Wuhan + Trigger_Stage_II * gamma_2_Stage_II_Wuhan + Trigger_Stage_III * gamma_2_Stage_III_Wuhan) + Trigger_Hubei * (Trigger_Stage_I * gamma_2_Stage_I_Hubei + Trigger_Stage_II * gamma_2_Stage_II_Hubei) + Trigger_China * (Trigger_Stage_I * gamma_2_Stage_I_China + Trigger_Stage_II * gamma_2_Stage_II_China))) * sigma) * Quarantined_Infected),
D(Susceptible) ~ -1.0 * (China * (((Trigger_Wuhan * (Trigger_Stage_I * R_Stage_I_Wuhan + Trigger_Stage_II * R_Stage_II_Wuhan + Trigger_Stage_III * R_Stage_III_Wuhan) + Trigger_Hubei * (Trigger_Stage_I * R_Stage_I_Hubei + Trigger_Stage_II * R_Stage_II_Hubei) + Trigger_China * (Trigger_Stage_I * R_Stage_I_China + Trigger_Stage_II * R_Stage_II_China)) * (Trigger_Wuhan * (Trigger_Stage_I * gamma_1_Stage_I_Wuhan + Trigger_Stage_II * gamma_1_Stage_II_Wuhan + Trigger_Stage_III * gamma_1_Stage_III_Wuhan) + Trigger_Hubei * (Trigger_Stage_I * gamma_1_Stage_I_Hubei + Trigger_Stage_II * gamma_1_Stage_II_Hubei) + Trigger_China * (Trigger_Stage_I * gamma_1_Stage_I_China + Trigger_Stage_II * gamma_1_Stage_II_China))) * Susceptible * Unquarantined_Infected / (Trigger_Wuhan * Total_Pop_Wuhan + Trigger_Hubei * Total_Pop_Hubei + Trigger_China * Total_Pop_China))),
D(Unquarantined_Infected) ~ +1.0 * (China * (((Trigger_Wuhan * (Trigger_Stage_I * R_Stage_I_Wuhan + Trigger_Stage_II * R_Stage_II_Wuhan + Trigger_Stage_III * R_Stage_III_Wuhan) + Trigger_Hubei * (Trigger_Stage_I * R_Stage_I_Hubei + Trigger_Stage_II * R_Stage_II_Hubei) + Trigger_China * (Trigger_Stage_I * R_Stage_I_China + Trigger_Stage_II * R_Stage_II_China)) * (Trigger_Wuhan * (Trigger_Stage_I * gamma_1_Stage_I_Wuhan + Trigger_Stage_II * gamma_1_Stage_II_Wuhan + Trigger_Stage_III * gamma_1_Stage_III_Wuhan) + Trigger_Hubei * (Trigger_Stage_I * gamma_1_Stage_I_Hubei + Trigger_Stage_II * gamma_1_Stage_II_Hubei) + Trigger_China * (Trigger_Stage_I * gamma_1_Stage_I_China + Trigger_Stage_II * gamma_1_Stage_II_China))) * Susceptible * Unquarantined_Infected / (Trigger_Wuhan * Total_Pop_Wuhan + Trigger_Hubei * Total_Pop_Hubei + Trigger_China * Total_Pop_China)))-1.0 * (China * (Trigger_Wuhan * (Trigger_Stage_I * gamma_1_Stage_I_Wuhan + Trigger_Stage_II * gamma_1_Stage_II_Wuhan + Trigger_Stage_III * gamma_1_Stage_III_Wuhan) + Trigger_Hubei * (Trigger_Stage_I * gamma_1_Stage_I_Hubei + Trigger_Stage_II * gamma_1_Stage_II_Hubei) + Trigger_China * (Trigger_Stage_I * gamma_1_Stage_I_China + Trigger_Stage_II * gamma_1_Stage_II_China)) * Unquarantined_Infected)]

@named sys = ODESystem(eqs)

### Initial species concentrations ###
initialSpeciesValues = [
Quarantined_Infected => Trigger_Wuhan * (Trigger_Stage_I * 0 + Trigger_Stage_II * 0 + Trigger_Stage_III * 5000) + Trigger_Hubei * (Trigger_Stage_I * 0 + Trigger_Stage_II * 1500) + Trigger_China * (Trigger_Stage_I * 0 + Trigger_Stage_II * 2000),
Confirmed_Infected => Trigger_Wuhan * (Trigger_Stage_I * 258 + Trigger_Stage_II * 2000 + Trigger_Stage_III * 36000) + Trigger_Hubei * (Trigger_Stage_I * 0 + Trigger_Stage_II * 1600) + Trigger_China * (Trigger_Stage_I * 0 + Trigger_Stage_II * 4000),
Susceptible => (Trigger_Wuhan * Total_Pop_Wuhan + Trigger_Hubei * Total_Pop_Hubei + Trigger_China * Total_Pop_China),
Unquarantined_Infected => Trigger_Wuhan * (Trigger_Stage_I * 258 + Trigger_Stage_II * 15270 + Trigger_Stage_III * 4000) + Trigger_Hubei * (Trigger_Stage_I * 270 + Trigger_Stage_II * 5700) + Trigger_China * (Trigger_Stage_I * 291 + Trigger_Stage_II * 2800)]

### True parameter values ###
trueParameterValues = [
Total_Pop_Hubei => 4.8e7,
R_Stage_I_Wuhan => 4.7092,
gamma_2_Stage_III_Wuhan => 0.322,
Total_Pop_China => 1.335e9,
Trigger_Wuhan => 1.0,
gamma_1_Stage_I_China => 0.1941,
gamma_1_Stage_I_Hubei => 0.05,
sigma => 0.0,
R_Stage_II_Wuhan => 0.7575,
R_Stage_II_China => 0.5753,
gamma_1_Stage_I_Wuhan => 0.063,
gamma_2_Stage_II_Wuhan => 0.0643,
Trigger_China => 0.0,
gamma_2_Stage_II_China => 0.2189,
R_Stage_III_Wuhan => 0.4797,
Trigger_Stage_I => 1.0,
gamma_2_Stage_I_Wuhan => 0.05,
R_Stage_I_China => 1.5283,
Trigger_Stage_III => 0.0,
gamma_1_Stage_II_China => 0.5157,
gamma_1_Stage_III_Wuhan => 0.6185,
gamma_2_Stage_I_China => 0.05,
Trigger_Stage_II => 0.0,
R_Stage_II_Hubei => 0.6079,
Total_Pop_Wuhan => 9.01e6,
gamma_1_Stage_II_Wuhan => 0.3917,
gamma_2_Stage_I_Hubei => 0.05,
R_Stage_I_Hubei => 5.934,
gamma_1_Stage_II_Hubei => 0.488,
Trigger_Hubei => 0.0,
gamma_2_Stage_II_Hubei => 0.1914]

trueConstantsValues = [
China => 1.0]
