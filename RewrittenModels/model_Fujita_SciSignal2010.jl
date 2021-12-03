# Model name: model_Fujita_SciSignal2010
# Number of parameters: 20
# Number of species: 9

### Define independent and dependent variables
@variables t pAkt_S6(t) pAkt(t) pS6(t) EGFR(t) pEGFR_Akt(t) pEGFR(t) Akt(t) S6(t) EGF_EGFR(t)

### Define variable parameters
@variables EGF(t)

### Define dummy variable
@variables dummyVariable(t)

### Define parameters
@parameters EGF_end reaction_5_k1 reaction_2_k2 init_AKT init_EGFR EGF_rate EGFR_turnover reaction_1_k1 reaction_1_k2 reaction_8_k1 reaction_4_k1 reaction_6_k1 reaction_2_k1 init_S6 reaction_7_k1 reaction_9_k1 reaction_3_k1 reaction_5_k2 EGF_0

### Define constants
@parameters Cell

### Define an operator for the differentiation w.r.t. time
D = Differential(t)

### Function definitions ###

### Events ###
continuous_events = [
[t ~ EGF_end] => [EGF ~ 0]
]

### Derivatives ###
eqs = [
D(pAkt_S6) ~ +1.0 * (Cell * (S6 * pAkt * reaction_5_k1 - pAkt_S6 * reaction_5_k2))-1.0 * (Cell * pAkt_S6 * reaction_6_k1),
D(pAkt) ~ +1.0 * (Cell * pEGFR_Akt * reaction_3_k1)-1.0 * (Cell * (S6 * pAkt * reaction_5_k1 - pAkt_S6 * reaction_5_k2))+1.0 * (Cell * pAkt_S6 * reaction_6_k1)-1.0 * (Cell * pAkt * reaction_7_k1),
D(pS6) ~ +1.0 * (Cell * pAkt_S6 * reaction_6_k1)-1.0 * (Cell * pS6 * reaction_8_k1),
D(EGFR) ~ -1.0 * (Cell * (EGF * EGFR * reaction_1_k1 - EGF_EGFR * reaction_1_k2))-1.0 * (Cell * EGFR * EGFR_turnover)+1.0 * (Cell * 68190 * EGFR_turnover),
D(pEGFR_Akt) ~ +1.0 * (Cell * (Akt * pEGFR * reaction_2_k1 - pEGFR_Akt * reaction_2_k2))-1.0 * (Cell * pEGFR_Akt * reaction_3_k1),
D(pEGFR) ~ -1.0 * (Cell * (Akt * pEGFR * reaction_2_k1 - pEGFR_Akt * reaction_2_k2))+1.0 * (Cell * pEGFR_Akt * reaction_3_k1)-1.0 * (Cell * pEGFR * reaction_4_k1)+1.0 * (Cell * EGF_EGFR * reaction_9_k1),
D(Akt) ~ -1.0 * (Cell * (Akt * pEGFR * reaction_2_k1 - pEGFR_Akt * reaction_2_k2))+1.0 * (Cell * pAkt * reaction_7_k1),
D(S6) ~ -1.0 * (Cell * (S6 * pAkt * reaction_5_k1 - pAkt_S6 * reaction_5_k2))+1.0 * (Cell * pS6 * reaction_8_k1),
D(EGF_EGFR) ~ +1.0 * (Cell * (EGF * EGFR * reaction_1_k1 - EGF_EGFR * reaction_1_k2))-1.0 * (Cell * EGF_EGFR * reaction_9_k1),
D(EGF) ~ 0,
D(dummyVariable) ~ +EGF_end+EGF_rate+EGF_0
]

@named sys = ODESystem(eqs, t, continuous_events = continuous_events)

### Initial species concentrations ###
initialSpeciesValues = [
pAkt_S6 => 0.0,
pAkt => 0.0,
pS6 => 0.0,
EGFR => init_EGFR,
pEGFR_Akt => 0.0,
pEGFR => 0.0,
Akt => init_AKT,
S6 => init_S6,
EGF_EGFR => 0.0,
EGF => 0.0,
dummyVariable => 0.0]

### True parameter values ###
trueParameterValues = [
EGF_end => 0.0,
reaction_5_k1 => 2.9643709900602e-6,
reaction_2_k2 => 41469.6914053245,
init_AKT => 0.00332683237159935,
init_EGFR => 2.26508055977911e7,
EGF_rate => 0.0,
EGFR_turnover => 0.001449799125736,
reaction_1_k1 => 0.00372345533395159,
reaction_1_k2 => 0.00262709856442467,
reaction_8_k1 => 0.000941161525754959,
reaction_4_k1 => 0.0308146966905863,
reaction_6_k1 => 9.20585474645043e-6,
reaction_2_k1 => 0.00103236148008131,
init_S6 => 205.86301335244,
reaction_7_k1 => 0.0119329694583145,
reaction_9_k1 => 0.0273281571867514,
reaction_3_k1 => 0.454840577578597,
reaction_5_k2 => 0.000404055756190126,
EGF_0 => 0.0]

trueConstantsValues = [
Cell => 1.0]
