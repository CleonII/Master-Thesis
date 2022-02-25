# Model name: model_Laske_PLOSComputBiol2019
# Number of parameters: 90
# Number of species: 41

### Define independent and dependent variables
ModelingToolkit.@variables t P_PA(t) Vp_cyt(t) Vp_nuc(t) P_NA(t) R_C_RdRp(t) R_M6(t) P_NEP(t) V_end(t) B_att_Hi(t) R_M5(t) Vp_cyt_M1(t) V_rel(t) R_C(t) Vp_nuc_M1(t) Cp(t) V_att_Lo(t) V_att_Hi(t) B_att_Lo(t) P_HA(t) P_M1(t) R_M1(t) P_M2(t) R_M7(t) P_NP(t) V_ex(t) R_M4(t) R_M8(t) R_V_RdRp(t) R_M2(t) R_V(t) P_B2(t) P_RdRp(t) R_M3(t) P_B1(t)

### Define variable parameters

### Define dummy variable

### Define parameters
ModelingToolkit.@parameters k_deg_R_RdRp ModelValue_104 ModelValue_63 ModelValue_114 ModelValue_108 ModelValue_64 L6 k_end N_P_RdRp D_rib ModelValue_111 ModelValue_107 k_rel ModelValue_105 L1 K_eq_Hi N_P_NP L2 ModelValue_113 k_bind_M1 ModelValue_116 F_Spl7 k_deg_Rnp N_P_NA K_eq_Lo ModelValue_101 k_exp_Vp_nuc_M1 ModelValue_79 ModelValue_69 ModelValue_90 ModelValue_91 ModelValue_103 N_P_M1 k_syn_R_M L8 ModelValue_84 ModelValue_85 ModelValue_115 k_att_Lo k_fus N_P_HA k_deg_R N_P_NEP k_bind_RdRp ModelValue_80 ModelValue_89 k_RdRp L5 N_P_M2 k_syn_R_V L4 ModelValue_87 k_bind_NP L3 k_imp ModelValue_106 F_Spl8 k_deg_R_M k_syn_P ModelValue_82 k_syn_R_C L7 K_V_rel ModelValue_86 ModelValue_102 F_fus k_att_Hi ModelValue_88

### Define constants
ModelingToolkit.@parameters compartment

### Define an operator for the differentiation w.r.t. time
D = Differential(t)

### Events ###

### Derivatives ###
eqs = [
D(P_PA) ~ +1.0 * (compartment * (ModelValue_80 / ModelValue_82) * R_M3)-1.0 * (compartment * k_RdRp * P_B1 * P_B2 * P_PA),
D(Vp_cyt) ~ +8.0 * (compartment * k_fus * V_end)-1.0 * (compartment * k_imp * Vp_cyt),
D(Vp_nuc) ~ +1.0 * (compartment * k_imp * Vp_cyt)-1.0 * (compartment * k_syn_R_C * Vp_nuc)+1.0 * (compartment * k_syn_R_C * Vp_nuc)+1.0 * (compartment * (k_bind_NP * R_V_RdRp * P_NP))-1.0 * (compartment * k_deg_Rnp * Vp_nuc)-1.0 * (compartment * (ModelValue_79 / ModelValue_84 / 8) * Vp_nuc)+1.0 * (compartment * (ModelValue_79 / ModelValue_84 / 8) * Vp_nuc)-1.0 * (compartment * (ModelValue_79 / ModelValue_85 / 8) * Vp_nuc)+1.0 * (compartment * (ModelValue_79 / ModelValue_85 / 8) * Vp_nuc)-1.0 * (compartment * (ModelValue_79 / ModelValue_86 / 8) * Vp_nuc)+1.0 * (compartment * (ModelValue_79 / ModelValue_86 / 8) * Vp_nuc)-1.0 * (compartment * (ModelValue_79 / ModelValue_87 / 8) * Vp_nuc)+1.0 * (compartment * (ModelValue_79 / ModelValue_87 / 8) * Vp_nuc)-1.0 * (compartment * (ModelValue_79 / ModelValue_88 / 8) * Vp_nuc)+1.0 * (compartment * (ModelValue_79 / ModelValue_88 / 8) * Vp_nuc)-1.0 * (compartment * (ModelValue_79 / ModelValue_89 / 8) * Vp_nuc)+1.0 * (compartment * (ModelValue_79 / ModelValue_89 / 8) * Vp_nuc)-1.0 * (compartment * (ModelValue_79 / ModelValue_90 / 8) * Vp_nuc)+1.0 * (compartment * (ModelValue_79 / ModelValue_90 / 8) * Vp_nuc)-1.0 * (compartment * (ModelValue_79 / ModelValue_91 / 8) * Vp_nuc)+1.0 * (compartment * (ModelValue_79 / ModelValue_91 / 8) * Vp_nuc)-1.0 * (compartment * (k_bind_M1 * Vp_nuc * P_M1)),
D(P_NA) ~ +1.0 * (compartment * (ModelValue_80 / ModelValue_82) * R_M6)-100.0 * (compartment * (k_rel * Vp_cyt_M1 * P_RdRp * P_HA * P_NP * P_NA * P_M1 * P_M2 * P_NEP / ((P_RdRp + (ModelValue_107 * ModelValue_101)) * (P_HA + (ModelValue_107 * ModelValue_102)) * (P_NP + (ModelValue_107 * ModelValue_106)) * (P_NA + (ModelValue_107 * ModelValue_104)) * (P_M1 + (ModelValue_107 * ModelValue_103)) * (P_M2 + (ModelValue_107 * ModelValue_113)) * (P_NEP + (ModelValue_107 * ModelValue_105))))),
D(R_C_RdRp) ~ +1.0 * (compartment * k_bind_RdRp * R_C * P_RdRp)-1.0 * (compartment * k_deg_R_RdRp * R_C_RdRp)-1.0 * (compartment * (k_bind_NP * R_C_RdRp * P_NP)),
D(R_M6) ~ +1.0 * (compartment * (ModelValue_79 / ModelValue_89 / 8) * Vp_nuc)-1.0 * (compartment * (ModelValue_80 / ModelValue_82) * R_M6)+1.0 * (compartment * (ModelValue_80 / ModelValue_82) * R_M6)-1.0 * (compartment * k_deg_R_M * R_M6),
D(P_NEP) ~ +1.0 * (compartment * (ModelValue_80 / ModelValue_82 * ModelValue_111) * R_M8)-1.0 * (compartment * k_exp_Vp_nuc_M1 * Vp_nuc_M1 * P_NEP)-157.0 * (compartment * (k_rel * Vp_cyt_M1 * P_RdRp * P_HA * P_NP * P_NA * P_M1 * P_M2 * P_NEP / ((P_RdRp + (ModelValue_107 * ModelValue_101)) * (P_HA + (ModelValue_107 * ModelValue_102)) * (P_NP + (ModelValue_107 * ModelValue_106)) * (P_NA + (ModelValue_107 * ModelValue_104)) * (P_M1 + (ModelValue_107 * ModelValue_103)) * (P_M2 + (ModelValue_107 * ModelValue_113)) * (P_NEP + (ModelValue_107 * ModelValue_105))))),
D(V_end) ~ +1.0 * (compartment * k_end * V_att_Hi)-1.0 * (compartment * ((1 - ModelValue_116) / ModelValue_116 * ModelValue_69) * V_end)-1.0 * (compartment * k_fus * V_end)+1.0 * (compartment * k_end * V_att_Lo),
D(B_att_Hi) ~ -1.0 * (compartment * k_att_Hi * V_ex * B_att_Hi)+1.0 * (compartment * (ModelValue_63 / ModelValue_114) * V_att_Hi)+1.0 * (compartment * k_end * V_att_Hi),
D(R_M5) ~ +1.0 * (compartment * (ModelValue_79 / ModelValue_88 / 8) * Vp_nuc)-1.0 * (compartment * (ModelValue_80 / ModelValue_82) * R_M5)+1.0 * (compartment * (ModelValue_80 / ModelValue_82) * R_M5)-1.0 * (compartment * k_deg_R_M * R_M5),
D(Vp_cyt_M1) ~ +1.0 * (compartment * k_exp_Vp_nuc_M1 * Vp_nuc_M1 * P_NEP)-8.0 * (compartment * (k_rel * Vp_cyt_M1 * P_RdRp * P_HA * P_NP * P_NA * P_M1 * P_M2 * P_NEP / ((P_RdRp + (ModelValue_107 * ModelValue_101)) * (P_HA + (ModelValue_107 * ModelValue_102)) * (P_NP + (ModelValue_107 * ModelValue_106)) * (P_NA + (ModelValue_107 * ModelValue_104)) * (P_M1 + (ModelValue_107 * ModelValue_103)) * (P_M2 + (ModelValue_107 * ModelValue_113)) * (P_NEP + (ModelValue_107 * ModelValue_105)))))-1.0 * (compartment * k_deg_Rnp * Vp_cyt_M1),
D(V_rel) ~ +1.0 * (compartment * (k_rel * Vp_cyt_M1 * P_RdRp * P_HA * P_NP * P_NA * P_M1 * P_M2 * P_NEP / ((P_RdRp + (ModelValue_107 * ModelValue_101)) * (P_HA + (ModelValue_107 * ModelValue_102)) * (P_NP + (ModelValue_107 * ModelValue_106)) * (P_NA + (ModelValue_107 * ModelValue_104)) * (P_M1 + (ModelValue_107 * ModelValue_103)) * (P_M2 + (ModelValue_107 * ModelValue_113)) * (P_NEP + (ModelValue_107 * ModelValue_105))))),
D(R_C) ~ +1.0 * (compartment * k_syn_R_C * Vp_nuc)-1.0 * (compartment * k_deg_R * R_C)-1.0 * (compartment * k_bind_RdRp * R_C * P_RdRp),
D(Vp_nuc_M1) ~ +1.0 * (compartment * (k_bind_M1 * Vp_nuc * P_M1))-1.0 * (compartment * k_exp_Vp_nuc_M1 * Vp_nuc_M1 * P_NEP)-1.0 * (compartment * k_deg_Rnp * Vp_nuc_M1),
D(Cp) ~ +1.0 * (compartment * (k_bind_NP * R_C_RdRp * P_NP))-1.0 * (compartment * k_syn_R_V * Cp)+1.0 * (compartment * k_syn_R_V * Cp)-1.0 * (compartment * k_deg_Rnp * Cp),
D(V_att_Lo) ~ +1.0 * (compartment * k_att_Lo * V_ex * B_att_Lo)-1.0 * (compartment * (ModelValue_64 / ModelValue_115) * V_att_Lo)-1.0 * (compartment * k_end * V_att_Lo),
D(V_att_Hi) ~ +1.0 * (compartment * k_att_Hi * V_ex * B_att_Hi)-1.0 * (compartment * (ModelValue_63 / ModelValue_114) * V_att_Hi)-1.0 * (compartment * k_end * V_att_Hi),
D(B_att_Lo) ~ -1.0 * (compartment * k_att_Lo * V_ex * B_att_Lo)+1.0 * (compartment * (ModelValue_64 / ModelValue_115) * V_att_Lo)+1.0 * (compartment * k_end * V_att_Lo),
D(P_HA) ~ +1.0 * (compartment * (ModelValue_80 / ModelValue_82) * R_M4)-500.0 * (compartment * (k_rel * Vp_cyt_M1 * P_RdRp * P_HA * P_NP * P_NA * P_M1 * P_M2 * P_NEP / ((P_RdRp + (ModelValue_107 * ModelValue_101)) * (P_HA + (ModelValue_107 * ModelValue_102)) * (P_NP + (ModelValue_107 * ModelValue_106)) * (P_NA + (ModelValue_107 * ModelValue_104)) * (P_M1 + (ModelValue_107 * ModelValue_103)) * (P_M2 + (ModelValue_107 * ModelValue_113)) * (P_NEP + (ModelValue_107 * ModelValue_105))))),
D(P_M1) ~ +1.0 * (compartment * (ModelValue_80 / ModelValue_82 * (1 - ModelValue_108)) * R_M7)-8.5 * (compartment * (k_bind_M1 * Vp_nuc * P_M1))-2932.0 * (compartment * (k_rel * Vp_cyt_M1 * P_RdRp * P_HA * P_NP * P_NA * P_M1 * P_M2 * P_NEP / ((P_RdRp + (ModelValue_107 * ModelValue_101)) * (P_HA + (ModelValue_107 * ModelValue_102)) * (P_NP + (ModelValue_107 * ModelValue_106)) * (P_NA + (ModelValue_107 * ModelValue_104)) * (P_M1 + (ModelValue_107 * ModelValue_103)) * (P_M2 + (ModelValue_107 * ModelValue_113)) * (P_NEP + (ModelValue_107 * ModelValue_105))))),
D(R_M1) ~ +1.0 * (compartment * (ModelValue_79 / ModelValue_84 / 8) * Vp_nuc)-1.0 * (compartment * (ModelValue_80 / ModelValue_82) * R_M1)+1.0 * (compartment * (ModelValue_80 / ModelValue_82) * R_M1)-1.0 * (compartment * k_deg_R_M * R_M1),
D(P_M2) ~ -40.0 * (compartment * (k_rel * Vp_cyt_M1 * P_RdRp * P_HA * P_NP * P_NA * P_M1 * P_M2 * P_NEP / ((P_RdRp + (ModelValue_107 * ModelValue_101)) * (P_HA + (ModelValue_107 * ModelValue_102)) * (P_NP + (ModelValue_107 * ModelValue_106)) * (P_NA + (ModelValue_107 * ModelValue_104)) * (P_M1 + (ModelValue_107 * ModelValue_103)) * (P_M2 + (ModelValue_107 * ModelValue_113)) * (P_NEP + (ModelValue_107 * ModelValue_105)))))+1.0 * (compartment * (ModelValue_80 / ModelValue_82 * ModelValue_108) * R_M7),
D(R_M7) ~ +1.0 * (compartment * (ModelValue_79 / ModelValue_90 / 8) * Vp_nuc)-1.0 * (compartment * (ModelValue_80 / ModelValue_82 * (1 - ModelValue_108)) * R_M7)+1.0 * (compartment * (ModelValue_80 / ModelValue_82 * (1 - ModelValue_108)) * R_M7)-1.0 * (compartment * k_deg_R_M * R_M7)-1.0 * (compartment * (ModelValue_80 / ModelValue_82 * ModelValue_108) * R_M7)+1.0 * (compartment * (ModelValue_80 / ModelValue_82 * ModelValue_108) * R_M7),
D(P_NP) ~ -71.0 * (compartment * (k_bind_NP * R_C_RdRp * P_NP))-71.0 * (compartment * (k_bind_NP * R_V_RdRp * P_NP))+1.0 * (compartment * (ModelValue_80 / ModelValue_82) * R_M5)-433.0 * (compartment * (k_rel * Vp_cyt_M1 * P_RdRp * P_HA * P_NP * P_NA * P_M1 * P_M2 * P_NEP / ((P_RdRp + (ModelValue_107 * ModelValue_101)) * (P_HA + (ModelValue_107 * ModelValue_102)) * (P_NP + (ModelValue_107 * ModelValue_106)) * (P_NA + (ModelValue_107 * ModelValue_104)) * (P_M1 + (ModelValue_107 * ModelValue_103)) * (P_M2 + (ModelValue_107 * ModelValue_113)) * (P_NEP + (ModelValue_107 * ModelValue_105))))),
D(V_ex) ~ -1.0 * (compartment * k_att_Hi * V_ex * B_att_Hi)+1.0 * (compartment * (ModelValue_63 / ModelValue_114) * V_att_Hi)-1.0 * (compartment * k_att_Lo * V_ex * B_att_Lo)+1.0 * (compartment * (ModelValue_64 / ModelValue_115) * V_att_Lo),
D(R_M4) ~ +1.0 * (compartment * (ModelValue_79 / ModelValue_87 / 8) * Vp_nuc)-1.0 * (compartment * (ModelValue_80 / ModelValue_82) * R_M4)+1.0 * (compartment * (ModelValue_80 / ModelValue_82) * R_M4)-1.0 * (compartment * k_deg_R_M * R_M4),
D(R_M8) ~ +1.0 * (compartment * (ModelValue_79 / ModelValue_91 / 8) * Vp_nuc)-1.0 * (compartment * (ModelValue_80 / ModelValue_82 * ModelValue_111) * R_M8)+1.0 * (compartment * (ModelValue_80 / ModelValue_82 * ModelValue_111) * R_M8)-1.0 * (compartment * k_deg_R_M * R_M8),
D(R_V_RdRp) ~ +1.0 * (compartment * k_bind_RdRp * R_V * P_RdRp)-1.0 * (compartment * (k_bind_NP * R_V_RdRp * P_NP))-1.0 * (compartment * k_deg_R_RdRp * R_V_RdRp),
D(R_M2) ~ +1.0 * (compartment * (ModelValue_79 / ModelValue_85 / 8) * Vp_nuc)-1.0 * (compartment * (ModelValue_80 / ModelValue_82) * R_M2)+1.0 * (compartment * (ModelValue_80 / ModelValue_82) * R_M2)-1.0 * (compartment * k_deg_R_M * R_M2),
D(R_V) ~ +1.0 * (compartment * k_syn_R_V * Cp)-1.0 * (compartment * k_bind_RdRp * R_V * P_RdRp)-1.0 * (compartment * k_deg_R * R_V),
D(P_B2) ~ +1.0 * (compartment * (ModelValue_80 / ModelValue_82) * R_M1)-1.0 * (compartment * k_RdRp * P_B1 * P_B2 * P_PA),
D(P_RdRp) ~ -1.0 * (compartment * k_bind_RdRp * R_C * P_RdRp)-1.0 * (compartment * k_bind_RdRp * R_V * P_RdRp)+1.0 * (compartment * k_RdRp * P_B1 * P_B2 * P_PA)-37.0 * (compartment * (k_rel * Vp_cyt_M1 * P_RdRp * P_HA * P_NP * P_NA * P_M1 * P_M2 * P_NEP / ((P_RdRp + (ModelValue_107 * ModelValue_101)) * (P_HA + (ModelValue_107 * ModelValue_102)) * (P_NP + (ModelValue_107 * ModelValue_106)) * (P_NA + (ModelValue_107 * ModelValue_104)) * (P_M1 + (ModelValue_107 * ModelValue_103)) * (P_M2 + (ModelValue_107 * ModelValue_113)) * (P_NEP + (ModelValue_107 * ModelValue_105))))),
D(R_M3) ~ +1.0 * (compartment * (ModelValue_79 / ModelValue_86 / 8) * Vp_nuc)-1.0 * (compartment * (ModelValue_80 / ModelValue_82) * R_M3)+1.0 * (compartment * (ModelValue_80 / ModelValue_82) * R_M3)-1.0 * (compartment * k_deg_R_M * R_M3),
D(P_B1) ~ +1.0 * (compartment * (ModelValue_80 / ModelValue_82) * R_M2)-1.0 * (compartment * k_RdRp * P_B1 * P_B2 * P_PA)]

@named sys = ODESystem(eqs)

### Initial species concentrations ###
initialSpeciesValues = [
P_PA => 0.0,
Vp_cyt => 0.0,
Vp_nuc => 0.0,
P_NA => 0.0,
R_C_RdRp => 0.0,
R_M6 => 0.0,
P_NEP => 0.0,
V_end => 0.0,
B_att_Hi => 150.0,
R_M5 => 0.0,
Vp_cyt_M1 => 0.0,
V_rel => 0.0,
R_C => 0.0,
Vp_nuc_M1 => 0.0,
Cp => 0.0,
V_att_Lo => 0.0,
V_att_Hi => 0.0,
B_att_Lo => 1000.0,
P_HA => 0.0,
P_M1 => 0.0,
R_M1 => 0.0,
P_M2 => 0.0,
R_M7 => 0.0,
P_NP => 0.0,
V_ex => 50.0,
R_M4 => 0.0,
R_M8 => 0.0,
R_V_RdRp => 0.0,
R_M2 => 0.0,
R_V => 0.0,
P_B2 => 0.0,
P_RdRp => 0.0,
R_M3 => 0.0,
P_B1 => 0.0]

### True parameter values ###
trueParameterValues = [
k_deg_R_RdRp => 4.25,
ModelValue_104 => 100.0,
ModelValue_63 => 0.0809,
ModelValue_114 => 0.0113,
ModelValue_108 => 0.02,
ModelValue_64 => 0.000455,
L6 => 1392.0,
k_end => 4.8,
N_P_RdRp => 45.0,
D_rib => 160.0,
ModelValue_111 => 0.125,
ModelValue_107 => 10.0,
k_rel => 0.0011,
ModelValue_105 => 165.0,
L1 => 2320.0,
K_eq_Hi => 0.0113,
N_P_NP => 1000.0,
L2 => 2320.0,
ModelValue_113 => 40.0,
k_bind_M1 => 1.82e-6,
ModelValue_116 => 0.51,
F_Spl7 => 0.02,
k_deg_Rnp => 0.09,
N_P_NA => 100.0,
K_eq_Lo => 8.33e-5,
ModelValue_101 => 45.0,
k_exp_Vp_nuc_M1 => 1.0e-6,
ModelValue_79 => 30600.0,
ModelValue_69 => 3.21,
ModelValue_90 => 1005.0,
ModelValue_91 => 868.0,
ModelValue_103 => 3000.0,
N_P_M1 => 3000.0,
k_syn_R_M => 30600.0,
L8 => 868.0,
ModelValue_84 => 2320.0,
ModelValue_85 => 2320.0,
ModelValue_115 => 8.33e-5,
k_att_Lo => 0.000455,
k_fus => 3.21,
N_P_HA => 500.0,
k_deg_R => 36.36,
N_P_NEP => 165.0,
k_bind_RdRp => 1.0,
ModelValue_80 => 64800.0,
ModelValue_89 => 1392.0,
k_RdRp => 1.0,
L5 => 1540.0,
N_P_M2 => 40.0,
k_syn_R_V => 100.93,
L4 => 1757.0,
ModelValue_87 => 1757.0,
k_bind_NP => 0.000301,
L3 => 2211.0,
k_imp => 0.296,
ModelValue_106 => 1000.0,
F_Spl8 => 0.125,
k_deg_R_M => 0.33,
k_syn_P => 64800.0,
ModelValue_82 => 160.0,
k_syn_R_C => 1.53,
L7 => 1005.0,
K_V_rel => 10.0,
ModelValue_86 => 2211.0,
ModelValue_102 => 500.0,
F_fus => 0.51,
k_att_Hi => 0.0809,
ModelValue_88 => 1540.0]

trueConstantsValues = [
compartment => 1.0]
