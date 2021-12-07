# Model name: model_Laske_PLOSComputBiol2019
# Number of parameters: 90
# Number of species: 41

### Define independent and dependent variables
@variables t P_PA(t) Vp_cyt(t) Vp_nuc(t) P_NA(t) R_C_RdRp(t) R_M6(t) P_NEP(t) V_end(t) B_att_Hi(t) R_M5(t) Vp_cyt_M1(t) V_rel(t) R_C(t) Vp_nuc_M1(t) Cp(t) V_att_Lo(t) V_att_Hi(t) B_att_Lo(t) P_HA(t) P_M1(t) R_M1(t) P_M2(t) R_M7(t) P_NP(t) V_ex(t) R_M4(t) R_M8(t) R_V_RdRp(t) R_M2(t) R_V(t) P_B2(t) P_RdRp(t) R_M3(t) P_B1(t)

### Define variable parameters

### Define dummy variable

### Define parameters
@parameters k_deg_R_RdRp ModelValue_104 ModelValue_63 ModelValue_114 ModelValue_108 ModelValue_64 L6 k_end N_P_RdRp D_rib ModelValue_111 ModelValue_107 k_rel ModelValue_105 L1 K_eq_Hi N_P_NP L2 ModelValue_113 k_bind_M1 ModelValue_116 F_Spl7 k_deg_Rnp N_P_NA K_eq_Lo ModelValue_101 k_exp_Vp_nuc_M1 ModelValue_79 ModelValue_69 ModelValue_90 ModelValue_91 ModelValue_103 N_P_M1 k_syn_R_M L8 ModelValue_84 ModelValue_85 ModelValue_115 k_att_Lo k_fus N_P_HA k_deg_R N_P_NEP k_bind_RdRp ModelValue_80 ModelValue_89 k_RdRp L5 N_P_M2 k_syn_R_V L4 ModelValue_87 k_bind_NP L3 k_imp ModelValue_106 F_Spl8 k_deg_R_M k_syn_P ModelValue_82 k_syn_R_C L7 K_V_rel ModelValue_86 ModelValue_102 F_fus k_att_Hi ModelValue_88

### Define constants
@parameters compartment

### Define an operator for the differentiation w.r.t. time
D = Differential(t)

### Function definitions ###
k_dis_Hi(ModelValue_63, ModelValue_114) = ModelValue_63 / ModelValue_114
@register k_dis_Hi(ModelValue_63, ModelValue_114)
k_syn_P_M1_rel(ModelValue_80, ModelValue_82, ModelValue_108) = ModelValue_80 / ModelValue_82 * (1 - ModelValue_108)
@register k_syn_P_M1_rel(ModelValue_80, ModelValue_82, ModelValue_108)
KmB(ModelValue_107, ModelValue_101) = ModelValue_107 * ModelValue_101
@register KmB(ModelValue_107, ModelValue_101)
k_syn_R_M5_rel(ModelValue_79, ModelValue_88) = ModelValue_79 / ModelValue_88 / 8
@register k_syn_R_M5_rel(ModelValue_79, ModelValue_88)
F_rnp_nuc(Vp_nuc_M1, Vp_nuc, V_end, Vp_cyt, Vp_cyt_M1) = RNP_nuc(Vp_nuc_M1, Vp_nuc) / (RNP_nuc(Vp_nuc_M1, Vp_nuc) + RNP_cyt(V_end, Vp_cyt, Vp_cyt_M1)) * 100
@register F_rnp_nuc(Vp_nuc_M1, Vp_nuc, V_end, Vp_cyt, Vp_cyt_M1)
KmD(ModelValue_107, ModelValue_106) = ModelValue_107 * ModelValue_106
@register KmD(ModelValue_107, ModelValue_106)
k_syn_P_NEP_rel(ModelValue_80, ModelValue_82, ModelValue_111) = ModelValue_80 / ModelValue_82 * ModelValue_111
@register k_syn_P_NEP_rel(ModelValue_80, ModelValue_82, ModelValue_111)
KmE(ModelValue_107, ModelValue_104) = ModelValue_107 * ModelValue_104
@register KmE(ModelValue_107, ModelValue_104)
KmG(ModelValue_107, ModelValue_113) = ModelValue_107 * ModelValue_113
@register KmG(ModelValue_107, ModelValue_113)
RNP_nuc(Vp_nuc_M1, Vp_nuc) = Vp_nuc_M1 + Vp_nuc
@register RNP_nuc(Vp_nuc_M1, Vp_nuc)
R_V_seg_tot(V_att_Hi, V_att_Lo, V_end, Vp_cyt_M1, Vp_nuc_M1, Vp_cyt, Vp_nuc, R_V) = R_V_total_0(V_att_Hi, V_att_Lo, V_end, Vp_cyt_M1, Vp_nuc_M1, Vp_cyt, Vp_nuc, R_V) / 8
@register R_V_seg_tot(V_att_Hi, V_att_Lo, V_end, Vp_cyt_M1, Vp_nuc_M1, Vp_cyt, Vp_nuc, R_V)
R_C_seg_tot(Cp, R_C_RdRp, R_C) = R_C_total_0(Cp, R_C_RdRp, R_C) / 8
@register R_C_seg_tot(Cp, R_C_RdRp, R_C)
R_V_total_0(V_att_Hi, V_att_Lo, V_end, Vp_cyt_M1, Vp_nuc_M1, Vp_cyt, Vp_nuc, R_V) = 8 * (V_att_Hi + V_att_Lo + V_end) + Vp_cyt_M1 + Vp_nuc_M1 + Vp_cyt + Vp_nuc + R_V
@register R_V_total_0(V_att_Hi, V_att_Lo, V_end, Vp_cyt_M1, Vp_nuc_M1, Vp_cyt, Vp_nuc, R_V)
RNP_cyt(V_end, Vp_cyt, Vp_cyt_M1) = 8 * V_end + Vp_cyt + Vp_cyt_M1
@register RNP_cyt(V_end, Vp_cyt, Vp_cyt_M1)
k_syn_R_M3_rel(ModelValue_79, ModelValue_86) = ModelValue_79 / ModelValue_86 / 8
@register k_syn_R_M3_rel(ModelValue_79, ModelValue_86)
k_deg_end(ModelValue_116, ModelValue_69) = (1 - ModelValue_116) / ModelValue_116 * ModelValue_69
@register k_deg_end(ModelValue_116, ModelValue_69)
KmC(ModelValue_107, ModelValue_102) = ModelValue_107 * ModelValue_102
@register KmC(ModelValue_107, ModelValue_102)
KmF(ModelValue_107, ModelValue_103) = ModelValue_107 * ModelValue_103
@register KmF(ModelValue_107, ModelValue_103)
k_syn_R_M2_rel(ModelValue_79, ModelValue_85) = ModelValue_79 / ModelValue_85 / 8
@register k_syn_R_M2_rel(ModelValue_79, ModelValue_85)
KmH(ModelValue_107, ModelValue_105) = ModelValue_107 * ModelValue_105
@register KmH(ModelValue_107, ModelValue_105)
k_syn_R_M1_rel(ModelValue_79, ModelValue_84) = ModelValue_79 / ModelValue_84 / 8
@register k_syn_R_M1_rel(ModelValue_79, ModelValue_84)
R_C_total_0(Cp, R_C_RdRp, R_C) = Cp + R_C_RdRp + R_C
@register R_C_total_0(Cp, R_C_RdRp, R_C)
k_syn_P_rel(ModelValue_80, ModelValue_82) = ModelValue_80 / ModelValue_82
@register k_syn_P_rel(ModelValue_80, ModelValue_82)
k_syn_P_M2_rel(ModelValue_80, ModelValue_82, ModelValue_108) = ModelValue_80 / ModelValue_82 * ModelValue_108
@register k_syn_P_M2_rel(ModelValue_80, ModelValue_82, ModelValue_108)
k_syn_R_M6_rel(ModelValue_79, ModelValue_89) = ModelValue_79 / ModelValue_89 / 8
@register k_syn_R_M6_rel(ModelValue_79, ModelValue_89)
k_syn_R_M4_rel(ModelValue_79, ModelValue_87) = ModelValue_79 / ModelValue_87 / 8
@register k_syn_R_M4_rel(ModelValue_79, ModelValue_87)
k_syn_R_M8_rel(ModelValue_79, ModelValue_91) = ModelValue_79 / ModelValue_91 / 8
@register k_syn_R_M8_rel(ModelValue_79, ModelValue_91)
k_dis_Lo(ModelValue_64, ModelValue_115) = ModelValue_64 / ModelValue_115
@register k_dis_Lo(ModelValue_64, ModelValue_115)
k_syn_R_M7_rel(ModelValue_79, ModelValue_90) = ModelValue_79 / ModelValue_90 / 8
@register k_syn_R_M7_rel(ModelValue_79, ModelValue_90)

### Events ###

### Derivatives ###
eqs = [
D(P_PA) ~ +1.0 * (compartment * k_syn_P_rel(ModelValue_80, ModelValue_82) * R_M3)-1.0 * (compartment * k_RdRp * P_B1 * P_B2 * P_PA),
D(Vp_cyt) ~ +8.0 * (compartment * k_fus * V_end)-1.0 * (compartment * k_imp * Vp_cyt),
D(Vp_nuc) ~ +1.0 * (compartment * k_imp * Vp_cyt)-1.0 * (compartment * k_syn_R_C * Vp_nuc)+1.0 * (compartment * k_syn_R_C * Vp_nuc)+1.0 * (compartment * (k_bind_NP * R_V_RdRp * P_NP))-1.0 * (compartment * k_deg_Rnp * Vp_nuc)-1.0 * (compartment * k_syn_R_M1_rel(ModelValue_79, ModelValue_84) * Vp_nuc)+1.0 * (compartment * k_syn_R_M1_rel(ModelValue_79, ModelValue_84) * Vp_nuc)-1.0 * (compartment * k_syn_R_M2_rel(ModelValue_79, ModelValue_85) * Vp_nuc)+1.0 * (compartment * k_syn_R_M2_rel(ModelValue_79, ModelValue_85) * Vp_nuc)-1.0 * (compartment * k_syn_R_M3_rel(ModelValue_79, ModelValue_86) * Vp_nuc)+1.0 * (compartment * k_syn_R_M3_rel(ModelValue_79, ModelValue_86) * Vp_nuc)-1.0 * (compartment * k_syn_R_M4_rel(ModelValue_79, ModelValue_87) * Vp_nuc)+1.0 * (compartment * k_syn_R_M4_rel(ModelValue_79, ModelValue_87) * Vp_nuc)-1.0 * (compartment * k_syn_R_M5_rel(ModelValue_79, ModelValue_88) * Vp_nuc)+1.0 * (compartment * k_syn_R_M5_rel(ModelValue_79, ModelValue_88) * Vp_nuc)-1.0 * (compartment * k_syn_R_M6_rel(ModelValue_79, ModelValue_89) * Vp_nuc)+1.0 * (compartment * k_syn_R_M6_rel(ModelValue_79, ModelValue_89) * Vp_nuc)-1.0 * (compartment * k_syn_R_M7_rel(ModelValue_79, ModelValue_90) * Vp_nuc)+1.0 * (compartment * k_syn_R_M7_rel(ModelValue_79, ModelValue_90) * Vp_nuc)-1.0 * (compartment * k_syn_R_M8_rel(ModelValue_79, ModelValue_91) * Vp_nuc)+1.0 * (compartment * k_syn_R_M8_rel(ModelValue_79, ModelValue_91) * Vp_nuc)-1.0 * (compartment * (k_bind_M1 * Vp_nuc * P_M1)),
D(P_NA) ~ +1.0 * (compartment * k_syn_P_rel(ModelValue_80, ModelValue_82) * R_M6)-100.0 * (compartment * (k_rel * Vp_cyt_M1 * P_RdRp * P_HA * P_NP * P_NA * P_M1 * P_M2 * P_NEP / ((P_RdRp + KmB(ModelValue_107, ModelValue_101)) * (P_HA + KmC(ModelValue_107, ModelValue_102)) * (P_NP + KmD(ModelValue_107, ModelValue_106)) * (P_NA + KmE(ModelValue_107, ModelValue_104)) * (P_M1 + KmF(ModelValue_107, ModelValue_103)) * (P_M2 + KmG(ModelValue_107, ModelValue_113)) * (P_NEP + KmH(ModelValue_107, ModelValue_105))))),
D(R_C_RdRp) ~ +1.0 * (compartment * k_bind_RdRp * R_C * P_RdRp)-1.0 * (compartment * k_deg_R_RdRp * R_C_RdRp)-1.0 * (compartment * (k_bind_NP * R_C_RdRp * P_NP)),
D(R_M6) ~ +1.0 * (compartment * k_syn_R_M6_rel(ModelValue_79, ModelValue_89) * Vp_nuc)-1.0 * (compartment * k_syn_P_rel(ModelValue_80, ModelValue_82) * R_M6)+1.0 * (compartment * k_syn_P_rel(ModelValue_80, ModelValue_82) * R_M6)-1.0 * (compartment * k_deg_R_M * R_M6),
D(P_NEP) ~ +1.0 * (compartment * k_syn_P_NEP_rel(ModelValue_80, ModelValue_82, ModelValue_111) * R_M8)-1.0 * (compartment * k_exp_Vp_nuc_M1 * Vp_nuc_M1 * P_NEP)-157.0 * (compartment * (k_rel * Vp_cyt_M1 * P_RdRp * P_HA * P_NP * P_NA * P_M1 * P_M2 * P_NEP / ((P_RdRp + KmB(ModelValue_107, ModelValue_101)) * (P_HA + KmC(ModelValue_107, ModelValue_102)) * (P_NP + KmD(ModelValue_107, ModelValue_106)) * (P_NA + KmE(ModelValue_107, ModelValue_104)) * (P_M1 + KmF(ModelValue_107, ModelValue_103)) * (P_M2 + KmG(ModelValue_107, ModelValue_113)) * (P_NEP + KmH(ModelValue_107, ModelValue_105))))),
D(V_end) ~ +1.0 * (compartment * k_end * V_att_Hi)-1.0 * (compartment * k_deg_end(ModelValue_116, ModelValue_69) * V_end)-1.0 * (compartment * k_fus * V_end)+1.0 * (compartment * k_end * V_att_Lo),
D(B_att_Hi) ~ -1.0 * (compartment * k_att_Hi * V_ex * B_att_Hi)+1.0 * (compartment * k_dis_Hi(ModelValue_63, ModelValue_114) * V_att_Hi)+1.0 * (compartment * k_end * V_att_Hi),
D(R_M5) ~ +1.0 * (compartment * k_syn_R_M5_rel(ModelValue_79, ModelValue_88) * Vp_nuc)-1.0 * (compartment * k_syn_P_rel(ModelValue_80, ModelValue_82) * R_M5)+1.0 * (compartment * k_syn_P_rel(ModelValue_80, ModelValue_82) * R_M5)-1.0 * (compartment * k_deg_R_M * R_M5),
D(Vp_cyt_M1) ~ +1.0 * (compartment * k_exp_Vp_nuc_M1 * Vp_nuc_M1 * P_NEP)-8.0 * (compartment * (k_rel * Vp_cyt_M1 * P_RdRp * P_HA * P_NP * P_NA * P_M1 * P_M2 * P_NEP / ((P_RdRp + KmB(ModelValue_107, ModelValue_101)) * (P_HA + KmC(ModelValue_107, ModelValue_102)) * (P_NP + KmD(ModelValue_107, ModelValue_106)) * (P_NA + KmE(ModelValue_107, ModelValue_104)) * (P_M1 + KmF(ModelValue_107, ModelValue_103)) * (P_M2 + KmG(ModelValue_107, ModelValue_113)) * (P_NEP + KmH(ModelValue_107, ModelValue_105)))))-1.0 * (compartment * k_deg_Rnp * Vp_cyt_M1),
D(V_rel) ~ +1.0 * (compartment * (k_rel * Vp_cyt_M1 * P_RdRp * P_HA * P_NP * P_NA * P_M1 * P_M2 * P_NEP / ((P_RdRp + KmB(ModelValue_107, ModelValue_101)) * (P_HA + KmC(ModelValue_107, ModelValue_102)) * (P_NP + KmD(ModelValue_107, ModelValue_106)) * (P_NA + KmE(ModelValue_107, ModelValue_104)) * (P_M1 + KmF(ModelValue_107, ModelValue_103)) * (P_M2 + KmG(ModelValue_107, ModelValue_113)) * (P_NEP + KmH(ModelValue_107, ModelValue_105))))),
D(R_C) ~ +1.0 * (compartment * k_syn_R_C * Vp_nuc)-1.0 * (compartment * k_deg_R * R_C)-1.0 * (compartment * k_bind_RdRp * R_C * P_RdRp),
D(Vp_nuc_M1) ~ +1.0 * (compartment * (k_bind_M1 * Vp_nuc * P_M1))-1.0 * (compartment * k_exp_Vp_nuc_M1 * Vp_nuc_M1 * P_NEP)-1.0 * (compartment * k_deg_Rnp * Vp_nuc_M1),
D(Cp) ~ +1.0 * (compartment * (k_bind_NP * R_C_RdRp * P_NP))-1.0 * (compartment * k_syn_R_V * Cp)+1.0 * (compartment * k_syn_R_V * Cp)-1.0 * (compartment * k_deg_Rnp * Cp),
D(V_att_Lo) ~ +1.0 * (compartment * k_att_Lo * V_ex * B_att_Lo)-1.0 * (compartment * k_dis_Lo(ModelValue_64, ModelValue_115) * V_att_Lo)-1.0 * (compartment * k_end * V_att_Lo),
D(V_att_Hi) ~ +1.0 * (compartment * k_att_Hi * V_ex * B_att_Hi)-1.0 * (compartment * k_dis_Hi(ModelValue_63, ModelValue_114) * V_att_Hi)-1.0 * (compartment * k_end * V_att_Hi),
D(B_att_Lo) ~ -1.0 * (compartment * k_att_Lo * V_ex * B_att_Lo)+1.0 * (compartment * k_dis_Lo(ModelValue_64, ModelValue_115) * V_att_Lo)+1.0 * (compartment * k_end * V_att_Lo),
D(P_HA) ~ +1.0 * (compartment * k_syn_P_rel(ModelValue_80, ModelValue_82) * R_M4)-500.0 * (compartment * (k_rel * Vp_cyt_M1 * P_RdRp * P_HA * P_NP * P_NA * P_M1 * P_M2 * P_NEP / ((P_RdRp + KmB(ModelValue_107, ModelValue_101)) * (P_HA + KmC(ModelValue_107, ModelValue_102)) * (P_NP + KmD(ModelValue_107, ModelValue_106)) * (P_NA + KmE(ModelValue_107, ModelValue_104)) * (P_M1 + KmF(ModelValue_107, ModelValue_103)) * (P_M2 + KmG(ModelValue_107, ModelValue_113)) * (P_NEP + KmH(ModelValue_107, ModelValue_105))))),
D(P_M1) ~ +1.0 * (compartment * k_syn_P_M1_rel(ModelValue_80, ModelValue_82, ModelValue_108) * R_M7)-8.5 * (compartment * (k_bind_M1 * Vp_nuc * P_M1))-2932.0 * (compartment * (k_rel * Vp_cyt_M1 * P_RdRp * P_HA * P_NP * P_NA * P_M1 * P_M2 * P_NEP / ((P_RdRp + KmB(ModelValue_107, ModelValue_101)) * (P_HA + KmC(ModelValue_107, ModelValue_102)) * (P_NP + KmD(ModelValue_107, ModelValue_106)) * (P_NA + KmE(ModelValue_107, ModelValue_104)) * (P_M1 + KmF(ModelValue_107, ModelValue_103)) * (P_M2 + KmG(ModelValue_107, ModelValue_113)) * (P_NEP + KmH(ModelValue_107, ModelValue_105))))),
D(R_M1) ~ +1.0 * (compartment * k_syn_R_M1_rel(ModelValue_79, ModelValue_84) * Vp_nuc)-1.0 * (compartment * k_syn_P_rel(ModelValue_80, ModelValue_82) * R_M1)+1.0 * (compartment * k_syn_P_rel(ModelValue_80, ModelValue_82) * R_M1)-1.0 * (compartment * k_deg_R_M * R_M1),
D(P_M2) ~ -40.0 * (compartment * (k_rel * Vp_cyt_M1 * P_RdRp * P_HA * P_NP * P_NA * P_M1 * P_M2 * P_NEP / ((P_RdRp + KmB(ModelValue_107, ModelValue_101)) * (P_HA + KmC(ModelValue_107, ModelValue_102)) * (P_NP + KmD(ModelValue_107, ModelValue_106)) * (P_NA + KmE(ModelValue_107, ModelValue_104)) * (P_M1 + KmF(ModelValue_107, ModelValue_103)) * (P_M2 + KmG(ModelValue_107, ModelValue_113)) * (P_NEP + KmH(ModelValue_107, ModelValue_105)))))+1.0 * (compartment * k_syn_P_M2_rel(ModelValue_80, ModelValue_82, ModelValue_108) * R_M7),
D(R_M7) ~ +1.0 * (compartment * k_syn_R_M7_rel(ModelValue_79, ModelValue_90) * Vp_nuc)-1.0 * (compartment * k_syn_P_M1_rel(ModelValue_80, ModelValue_82, ModelValue_108) * R_M7)+1.0 * (compartment * k_syn_P_M1_rel(ModelValue_80, ModelValue_82, ModelValue_108) * R_M7)-1.0 * (compartment * k_deg_R_M * R_M7)-1.0 * (compartment * k_syn_P_M2_rel(ModelValue_80, ModelValue_82, ModelValue_108) * R_M7)+1.0 * (compartment * k_syn_P_M2_rel(ModelValue_80, ModelValue_82, ModelValue_108) * R_M7),
D(P_NP) ~ -71.0 * (compartment * (k_bind_NP * R_C_RdRp * P_NP))-71.0 * (compartment * (k_bind_NP * R_V_RdRp * P_NP))+1.0 * (compartment * k_syn_P_rel(ModelValue_80, ModelValue_82) * R_M5)-433.0 * (compartment * (k_rel * Vp_cyt_M1 * P_RdRp * P_HA * P_NP * P_NA * P_M1 * P_M2 * P_NEP / ((P_RdRp + KmB(ModelValue_107, ModelValue_101)) * (P_HA + KmC(ModelValue_107, ModelValue_102)) * (P_NP + KmD(ModelValue_107, ModelValue_106)) * (P_NA + KmE(ModelValue_107, ModelValue_104)) * (P_M1 + KmF(ModelValue_107, ModelValue_103)) * (P_M2 + KmG(ModelValue_107, ModelValue_113)) * (P_NEP + KmH(ModelValue_107, ModelValue_105))))),
D(V_ex) ~ -1.0 * (compartment * k_att_Hi * V_ex * B_att_Hi)+1.0 * (compartment * k_dis_Hi(ModelValue_63, ModelValue_114) * V_att_Hi)-1.0 * (compartment * k_att_Lo * V_ex * B_att_Lo)+1.0 * (compartment * k_dis_Lo(ModelValue_64, ModelValue_115) * V_att_Lo),
D(R_M4) ~ +1.0 * (compartment * k_syn_R_M4_rel(ModelValue_79, ModelValue_87) * Vp_nuc)-1.0 * (compartment * k_syn_P_rel(ModelValue_80, ModelValue_82) * R_M4)+1.0 * (compartment * k_syn_P_rel(ModelValue_80, ModelValue_82) * R_M4)-1.0 * (compartment * k_deg_R_M * R_M4),
D(R_M8) ~ +1.0 * (compartment * k_syn_R_M8_rel(ModelValue_79, ModelValue_91) * Vp_nuc)-1.0 * (compartment * k_syn_P_NEP_rel(ModelValue_80, ModelValue_82, ModelValue_111) * R_M8)+1.0 * (compartment * k_syn_P_NEP_rel(ModelValue_80, ModelValue_82, ModelValue_111) * R_M8)-1.0 * (compartment * k_deg_R_M * R_M8),
D(R_V_RdRp) ~ +1.0 * (compartment * k_bind_RdRp * R_V * P_RdRp)-1.0 * (compartment * (k_bind_NP * R_V_RdRp * P_NP))-1.0 * (compartment * k_deg_R_RdRp * R_V_RdRp),
D(R_M2) ~ +1.0 * (compartment * k_syn_R_M2_rel(ModelValue_79, ModelValue_85) * Vp_nuc)-1.0 * (compartment * k_syn_P_rel(ModelValue_80, ModelValue_82) * R_M2)+1.0 * (compartment * k_syn_P_rel(ModelValue_80, ModelValue_82) * R_M2)-1.0 * (compartment * k_deg_R_M * R_M2),
D(R_V) ~ +1.0 * (compartment * k_syn_R_V * Cp)-1.0 * (compartment * k_bind_RdRp * R_V * P_RdRp)-1.0 * (compartment * k_deg_R * R_V),
D(P_B2) ~ +1.0 * (compartment * k_syn_P_rel(ModelValue_80, ModelValue_82) * R_M1)-1.0 * (compartment * k_RdRp * P_B1 * P_B2 * P_PA),
D(P_RdRp) ~ -1.0 * (compartment * k_bind_RdRp * R_C * P_RdRp)-1.0 * (compartment * k_bind_RdRp * R_V * P_RdRp)+1.0 * (compartment * k_RdRp * P_B1 * P_B2 * P_PA)-37.0 * (compartment * (k_rel * Vp_cyt_M1 * P_RdRp * P_HA * P_NP * P_NA * P_M1 * P_M2 * P_NEP / ((P_RdRp + KmB(ModelValue_107, ModelValue_101)) * (P_HA + KmC(ModelValue_107, ModelValue_102)) * (P_NP + KmD(ModelValue_107, ModelValue_106)) * (P_NA + KmE(ModelValue_107, ModelValue_104)) * (P_M1 + KmF(ModelValue_107, ModelValue_103)) * (P_M2 + KmG(ModelValue_107, ModelValue_113)) * (P_NEP + KmH(ModelValue_107, ModelValue_105))))),
D(R_M3) ~ +1.0 * (compartment * k_syn_R_M3_rel(ModelValue_79, ModelValue_86) * Vp_nuc)-1.0 * (compartment * k_syn_P_rel(ModelValue_80, ModelValue_82) * R_M3)+1.0 * (compartment * k_syn_P_rel(ModelValue_80, ModelValue_82) * R_M3)-1.0 * (compartment * k_deg_R_M * R_M3),
D(P_B1) ~ +1.0 * (compartment * k_syn_P_rel(ModelValue_80, ModelValue_82) * R_M2)-1.0 * (compartment * k_RdRp * P_B1 * P_B2 * P_PA)]

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
ModelValue_104 => N_P_NA,
ModelValue_63 => k_att_Hi,
ModelValue_114 => K_eq_Hi,
ModelValue_108 => F_Spl7,
ModelValue_64 => k_att_Lo,
L6 => 1392.0,
k_end => 4.8,
N_P_RdRp => 45.0,
D_rib => 160.0,
ModelValue_111 => F_Spl8,
ModelValue_107 => K_V_rel,
k_rel => 0.0011,
ModelValue_105 => N_P_NEP,
L1 => 2320.0,
K_eq_Hi => 0.0113,
N_P_NP => 1000.0,
L2 => 2320.0,
ModelValue_113 => N_P_M2,
k_bind_M1 => 1.82e-6,
ModelValue_116 => F_fus,
F_Spl7 => 0.02,
k_deg_Rnp => 0.09,
N_P_NA => 100.0,
K_eq_Lo => 8.33e-5,
ModelValue_101 => N_P_RdRp,
k_exp_Vp_nuc_M1 => 1.0e-6,
ModelValue_79 => k_syn_R_M,
ModelValue_69 => k_fus,
ModelValue_90 => L7,
ModelValue_91 => L8,
ModelValue_103 => N_P_M1,
N_P_M1 => 3000.0,
k_syn_R_M => 30600.0,
L8 => 868.0,
ModelValue_84 => L1,
ModelValue_85 => L2,
ModelValue_115 => K_eq_Lo,
k_att_Lo => 0.000455,
k_fus => 3.21,
N_P_HA => 500.0,
k_deg_R => 36.36,
N_P_NEP => 165.0,
k_bind_RdRp => 1.0,
ModelValue_80 => k_syn_P,
ModelValue_89 => L6,
k_RdRp => 1.0,
L5 => 1540.0,
N_P_M2 => 40.0,
k_syn_R_V => 100.93,
L4 => 1757.0,
ModelValue_87 => L4,
k_bind_NP => 0.000301,
L3 => 2211.0,
k_imp => 0.296,
ModelValue_106 => N_P_NP,
F_Spl8 => 0.125,
k_deg_R_M => 0.33,
k_syn_P => 64800.0,
ModelValue_82 => D_rib,
k_syn_R_C => 1.53,
L7 => 1005.0,
K_V_rel => 10.0,
ModelValue_86 => L3,
ModelValue_102 => N_P_HA,
F_fus => 0.51,
k_att_Hi => 0.0809,
ModelValue_88 => L5]

trueConstantsValues = [
compartment => 1.0]
