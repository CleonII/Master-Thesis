function evalDYmodDu(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	P_PA, Vp_cyt, Vp_nuc, P_NA, R_C_RdRp, R_M6, P_NEP, V_end, B_att_Hi, R_M5, Vp_cyt_M1, V_rel, R_C, Vp_nuc_M1, Cp, V_att_Lo, V_att_Hi, B_att_Lo, P_HA, P_M1, R_M1, P_M2, R_M7, P_NP, V_ex, R_M4, R_M8, R_V_RdRp, R_M2, R_V, P_B2, P_RdRp, R_M3, P_B1= u 
	k_deg_R_RdRp, ModelValue_104, ModelValue_63, ModelValue_114, ModelValue_108, ModelValue_64, L6, k_end, N_P_RdRp, D_rib, ModelValue_111, ModelValue_107, k_rel, ModelValue_105, compartment, L1, K_eq_Hi, N_P_NP, L2, ModelValue_113, k_bind_M1, ModelValue_116, F_Spl7, k_deg_Rnp, N_P_NA, K_eq_Lo, ModelValue_101, k_exp_Vp_nuc_M1, ModelValue_79, ModelValue_69, ModelValue_90, ModelValue_91, ModelValue_103, N_P_M1, k_syn_R_M, L8, ModelValue_84, ModelValue_85, ModelValue_115, k_att_Lo, k_fus, N_P_HA, k_deg_R, N_P_NEP, k_bind_RdRp, ModelValue_80, ModelValue_89, k_RdRp, L5, N_P_M2, k_syn_R_V, L4, ModelValue_87, k_bind_NP, L3, k_imp, ModelValue_106, F_Spl8, k_deg_R_M, k_syn_P, ModelValue_82, k_syn_R_C, L7, K_V_rel, ModelValue_86, ModelValue_102, F_fus, k_att_Hi, ModelValue_88 = p 
	F_rnp_nuc = (Vp_nuc_M1 + Vp_nuc) / ((Vp_nuc_M1 + Vp_nuc) + (8 * V_end + Vp_cyt + Vp_cyt_M1)) * 100
	R_V_seg_tot = (8 * (V_att_Hi + V_att_Lo + V_end) + Vp_cyt_M1 + Vp_nuc_M1 + Vp_cyt + Vp_nuc + R_V) / 8
	R_C_seg_tot = (Cp + R_C_RdRp + R_C) / 8

	if observableId == "RM5" 
		out[10] = 1
		return nothing
	end

	if observableId == "RVSegTot" 
		return nothing
	end

	if observableId == "RCSegTot" 
		return nothing
	end

	if observableId == "Vrel" 
		out[12] = 1
		return nothing
	end

	if observableId == "IntNucOffset" 
		return nothing
	end

	if observableId == "FracNucInt_1" 
		return nothing
	end

	if observableId == "FracNucInt_2" 
		return nothing
	end

	if observableId == "FracNucInt_3" 
		return nothing
	end

	if observableId == "FracNucInt_4" 
		return nothing
	end

	if observableId == "FracNucInt_5" 
		return nothing
	end

	if observableId == "FracNucInt_6" 
		return nothing
	end

	if observableId == "FracNucInt_7" 
		return nothing
	end

	if observableId == "FracNucInt_8" 
		return nothing
	end

end

function evalDYmodDp(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	P_PA, Vp_cyt, Vp_nuc, P_NA, R_C_RdRp, R_M6, P_NEP, V_end, B_att_Hi, R_M5, Vp_cyt_M1, V_rel, R_C, Vp_nuc_M1, Cp, V_att_Lo, V_att_Hi, B_att_Lo, P_HA, P_M1, R_M1, P_M2, R_M7, P_NP, V_ex, R_M4, R_M8, R_V_RdRp, R_M2, R_V, P_B2, P_RdRp, R_M3, P_B1= u 
	k_deg_R_RdRp, ModelValue_104, ModelValue_63, ModelValue_114, ModelValue_108, ModelValue_64, L6, k_end, N_P_RdRp, D_rib, ModelValue_111, ModelValue_107, k_rel, ModelValue_105, compartment, L1, K_eq_Hi, N_P_NP, L2, ModelValue_113, k_bind_M1, ModelValue_116, F_Spl7, k_deg_Rnp, N_P_NA, K_eq_Lo, ModelValue_101, k_exp_Vp_nuc_M1, ModelValue_79, ModelValue_69, ModelValue_90, ModelValue_91, ModelValue_103, N_P_M1, k_syn_R_M, L8, ModelValue_84, ModelValue_85, ModelValue_115, k_att_Lo, k_fus, N_P_HA, k_deg_R, N_P_NEP, k_bind_RdRp, ModelValue_80, ModelValue_89, k_RdRp, L5, N_P_M2, k_syn_R_V, L4, ModelValue_87, k_bind_NP, L3, k_imp, ModelValue_106, F_Spl8, k_deg_R_M, k_syn_P, ModelValue_82, k_syn_R_C, L7, K_V_rel, ModelValue_86, ModelValue_102, F_fus, k_att_Hi, ModelValue_88 = p 
	F_rnp_nuc = (Vp_nuc_M1 + Vp_nuc) / ((Vp_nuc_M1 + Vp_nuc) + (8 * V_end + Vp_cyt + Vp_cyt_M1)) * 100
	R_V_seg_tot = (8 * (V_att_Hi + V_att_Lo + V_end) + Vp_cyt_M1 + Vp_nuc_M1 + Vp_cyt + Vp_nuc + R_V) / 8
	R_C_seg_tot = (Cp + R_C_RdRp + R_C) / 8

	if observableId == "RM5" 
		return nothing
	end

	if observableId == "RVSegTot" 
		return nothing
	end

	if observableId == "RCSegTot" 
		return nothing
	end

	if observableId == "Vrel" 
		return nothing
	end

	if observableId == "IntNucOffset" 
		return nothing
	end

	if observableId == "FracNucInt_1" 
		return nothing
	end

	if observableId == "FracNucInt_2" 
		return nothing
	end

	if observableId == "FracNucInt_3" 
		return nothing
	end

	if observableId == "FracNucInt_4" 
		return nothing
	end

	if observableId == "FracNucInt_5" 
		return nothing
	end

	if observableId == "FracNucInt_6" 
		return nothing
	end

	if observableId == "FracNucInt_7" 
		return nothing
	end

	if observableId == "FracNucInt_8" 
		return nothing
	end

end

function evalDSdDu!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	P_PA, Vp_cyt, Vp_nuc, P_NA, R_C_RdRp, R_M6, P_NEP, V_end, B_att_Hi, R_M5, Vp_cyt_M1, V_rel, R_C, Vp_nuc_M1, Cp, V_att_Lo, V_att_Hi, B_att_Lo, P_HA, P_M1, R_M1, P_M2, R_M7, P_NP, V_ex, R_M4, R_M8, R_V_RdRp, R_M2, R_V, P_B2, P_RdRp, R_M3, P_B1= u 
	k_deg_R_RdRp, ModelValue_104, ModelValue_63, ModelValue_114, ModelValue_108, ModelValue_64, L6, k_end, N_P_RdRp, D_rib, ModelValue_111, ModelValue_107, k_rel, ModelValue_105, compartment, L1, K_eq_Hi, N_P_NP, L2, ModelValue_113, k_bind_M1, ModelValue_116, F_Spl7, k_deg_Rnp, N_P_NA, K_eq_Lo, ModelValue_101, k_exp_Vp_nuc_M1, ModelValue_79, ModelValue_69, ModelValue_90, ModelValue_91, ModelValue_103, N_P_M1, k_syn_R_M, L8, ModelValue_84, ModelValue_85, ModelValue_115, k_att_Lo, k_fus, N_P_HA, k_deg_R, N_P_NEP, k_bind_RdRp, ModelValue_80, ModelValue_89, k_RdRp, L5, N_P_M2, k_syn_R_V, L4, ModelValue_87, k_bind_NP, L3, k_imp, ModelValue_106, F_Spl8, k_deg_R_M, k_syn_P, ModelValue_82, k_syn_R_C, L7, K_V_rel, ModelValue_86, ModelValue_102, F_fus, k_att_Hi, ModelValue_88 = p 
	if observableId == "RM5" 
		return nothing
	end

	if observableId == "RVSegTot" 
		return nothing
	end

	if observableId == "RCSegTot" 
		return nothing
	end

	if observableId == "Vrel" 
		return nothing
	end

	if observableId == "IntNucOffset" 
		return nothing
	end

	if observableId == "FracNucInt_1" 
		return nothing
	end

	if observableId == "FracNucInt_2" 
		return nothing
	end

	if observableId == "FracNucInt_3" 
		return nothing
	end

	if observableId == "FracNucInt_4" 
		return nothing
	end

	if observableId == "FracNucInt_5" 
		return nothing
	end

	if observableId == "FracNucInt_6" 
		return nothing
	end

	if observableId == "FracNucInt_7" 
		return nothing
	end

	if observableId == "FracNucInt_8" 
		return nothing
	end

end

function evalDSdDp!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	P_PA, Vp_cyt, Vp_nuc, P_NA, R_C_RdRp, R_M6, P_NEP, V_end, B_att_Hi, R_M5, Vp_cyt_M1, V_rel, R_C, Vp_nuc_M1, Cp, V_att_Lo, V_att_Hi, B_att_Lo, P_HA, P_M1, R_M1, P_M2, R_M7, P_NP, V_ex, R_M4, R_M8, R_V_RdRp, R_M2, R_V, P_B2, P_RdRp, R_M3, P_B1= u 
	k_deg_R_RdRp, ModelValue_104, ModelValue_63, ModelValue_114, ModelValue_108, ModelValue_64, L6, k_end, N_P_RdRp, D_rib, ModelValue_111, ModelValue_107, k_rel, ModelValue_105, compartment, L1, K_eq_Hi, N_P_NP, L2, ModelValue_113, k_bind_M1, ModelValue_116, F_Spl7, k_deg_Rnp, N_P_NA, K_eq_Lo, ModelValue_101, k_exp_Vp_nuc_M1, ModelValue_79, ModelValue_69, ModelValue_90, ModelValue_91, ModelValue_103, N_P_M1, k_syn_R_M, L8, ModelValue_84, ModelValue_85, ModelValue_115, k_att_Lo, k_fus, N_P_HA, k_deg_R, N_P_NEP, k_bind_RdRp, ModelValue_80, ModelValue_89, k_RdRp, L5, N_P_M2, k_syn_R_V, L4, ModelValue_87, k_bind_NP, L3, k_imp, ModelValue_106, F_Spl8, k_deg_R_M, k_syn_P, ModelValue_82, k_syn_R_C, L7, K_V_rel, ModelValue_86, ModelValue_102, F_fus, k_att_Hi, ModelValue_88 = p 
	if observableId == "RM5" 
		return nothing
	end

	if observableId == "RVSegTot" 
		return nothing
	end

	if observableId == "RCSegTot" 
		return nothing
	end

	if observableId == "Vrel" 
		return nothing
	end

	if observableId == "IntNucOffset" 
		return nothing
	end

	if observableId == "FracNucInt_1" 
		return nothing
	end

	if observableId == "FracNucInt_2" 
		return nothing
	end

	if observableId == "FracNucInt_3" 
		return nothing
	end

	if observableId == "FracNucInt_4" 
		return nothing
	end

	if observableId == "FracNucInt_5" 
		return nothing
	end

	if observableId == "FracNucInt_6" 
		return nothing
	end

	if observableId == "FracNucInt_7" 
		return nothing
	end

	if observableId == "FracNucInt_8" 
		return nothing
	end

end

