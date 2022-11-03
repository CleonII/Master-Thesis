function evalYmod(u, t, dynPar, obsPar, paramData, observableId, mapObsParam) 

	P_PA, Vp_cyt, Vp_nuc, P_NA, R_C_RdRp, R_M6, P_NEP, V_end, B_att_Hi, R_M5, Vp_cyt_M1, V_rel, R_C, Vp_nuc_M1, Cp, V_att_Lo, V_att_Hi, B_att_Lo, P_HA, P_M1, R_M1, P_M2, R_M7, P_NP, V_ex, R_M4, R_M8, R_V_RdRp, R_M2, R_V, P_B2, P_RdRp, R_M3, P_B1, dummyVariable= u 
	k_imp, k_syn_R_M, k_syn_R_C, k_syn_R_V, k_bind_M1, k_rel = dynPar 

	F_rnp_nuc = RNP_nuc / (RNP_nuc + RNP_cyt) * 100
	R_C_seg_tot = (Cp + R_C_RdRp + R_C) / 8
	R_V_seg_tot = (8*(V_att_Hi+V_att_Lo+V_end)+Vp_cyt_M1+Vp_nuc_M1+Vp_cyt+Vp_nuc+R_V) / 8.0

	if observableId == "RM5" 
		return R_M5 
	end

	if observableId == "RVSegTot" 
		observableParameter1_RVSegTot = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return observableParameter1_RVSegTot + R_V_seg_tot
	end

	if observableId == "RCSegTot" 
		return  R_C_seg_tot
	end

	if observableId == "Vrel" 
		return V_rel 
	end

	if observableId == "IntNucOffset" 
		observableParameter1_IntNucOffset = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return observableParameter1_IntNucOffset 
	end

	if observableId == "FracNucInt_1" 
		observableParameter1_FracNucInt_1 = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return F_rnp_nuc + observableParameter1_FracNucInt_1 
	end

	if observableId == "FracNucInt_2" 
		observableParameter1_FracNucInt_2 = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return F_rnp_nuc + observableParameter1_FracNucInt_2 
	end

	if observableId == "FracNucInt_3" 
		observableParameter1_FracNucInt_3 = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return F_rnp_nuc + observableParameter1_FracNucInt_3 
	end

	if observableId == "FracNucInt_4" 
		observableParameter1_FracNucInt_4 = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return F_rnp_nuc + observableParameter1_FracNucInt_4 
	end

	if observableId == "FracNucInt_5" 
		observableParameter1_FracNucInt_5 = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return F_rnp_nuc + observableParameter1_FracNucInt_5 
	end

	if observableId == "FracNucInt_6" 
		observableParameter1_FracNucInt_6 = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return F_rnp_nuc + observableParameter1_FracNucInt_6 
	end

	if observableId == "FracNucInt_7" 
		observableParameter1_FracNucInt_7 = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return F_rnp_nuc + observableParameter1_FracNucInt_7 
	end

	if observableId == "FracNucInt_8" 
		observableParameter1_FracNucInt_8 = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return F_rnp_nuc + observableParameter1_FracNucInt_8 
	end

end

function evalU0!(u0Vec, paramVec) 

	ModelValue_82, k_RdRp, ModelValue_80, k_imp, k_fus, k_deg_Rnp, k_bind_NP, k_bind_M1, k_rel, ModelValue_104, ModelValue_105, ModelValue_101, ModelValue_106, ModelValue_103, ModelValue_113, ModelValue_102, ModelValue_107, k_bind_RdRp, k_deg_R_RdRp, ModelValue_89, k_deg_R_M, ModelValue_79, k_exp_Vp_nuc_M1, ModelValue_111, ModelValue_69, k_end, ModelValue_116, ModelValue_114, ModelValue_63, k_att_Hi, ModelValue_88, k_syn_R_C, k_deg_R, ModelValue_115, ModelValue_64, k_att_Lo, ModelValue_108, ModelValue_84, ModelValue_90, ModelValue_87, ModelValue_91, ModelValue_85, k_syn_R_V, ModelValue_86, N_P_NEP, L5, D_rib, N_P_M1, N_P_NA, L4, N_P_RdRp, F_Spl7, N_P_HA, K_eq_Lo, L3, L6, K_V_rel, F_Spl8, L1, L7, N_P_M2, L8, L2, k_syn_P, k_syn_R_M, K_eq_Hi, N_P_NP, F_fus = paramVec 

	P_PA = 0.0 
	Vp_cyt = 0.0 
	Vp_nuc = 0.0 
	P_NA = 0.0 
	R_C_RdRp = 0.0 
	R_M6 = 0.0 
	P_NEP = 0.0 
	V_end = 0.0 
	B_att_Hi = 150.0 
	R_M5 = 0.0 
	Vp_cyt_M1 = 0.0 
	V_rel = 0.0 
	R_C = 0.0 
	Vp_nuc_M1 = 0.0 
	Cp = 0.0 
	V_att_Lo = 0.0 
	V_att_Hi = 0.0 
	B_att_Lo = 1000.0 
	P_HA = 0.0 
	P_M1 = 0.0 
	R_M1 = 0.0 
	P_M2 = 0.0 
	R_M7 = 0.0 
	P_NP = 0.0 
	V_ex = 50.0 
	R_M4 = 0.0 
	R_M8 = 0.0 
	R_V_RdRp = 0.0 
	R_M2 = 0.0 
	R_V = 0.0 
	P_B2 = 0.0 
	P_RdRp = 0.0 
	R_M3 = 0.0 
	P_B1 = 0.0 
	dummyVariable = 0.0 

	u0Vec .= P_PA, Vp_cyt, Vp_nuc, P_NA, R_C_RdRp, R_M6, P_NEP, V_end, B_att_Hi, R_M5, Vp_cyt_M1, V_rel, R_C, Vp_nuc_M1, Cp, V_att_Lo, V_att_Hi, B_att_Lo, P_HA, P_M1, R_M1, P_M2, R_M7, P_NP, V_ex, R_M4, R_M8, R_V_RdRp, R_M2, R_V, P_B2, P_RdRp, R_M3, P_B1, dummyVariable
end

function evalSd!(u, t, sdPar, dynPar, paramData, observableId, mapSdParam) 

	P_PA, Vp_cyt, Vp_nuc, P_NA, R_C_RdRp, R_M6, P_NEP, V_end, B_att_Hi, R_M5, Vp_cyt_M1, V_rel, R_C, Vp_nuc_M1, Cp, V_att_Lo, V_att_Hi, B_att_Lo, P_HA, P_M1, R_M1, P_M2, R_M7, P_NP, V_ex, R_M4, R_M8, R_V_RdRp, R_M2, R_V, P_B2, P_RdRp, R_M3, P_B1, dummyVariable= u 
	k_imp, k_syn_R_M, k_syn_R_C, k_syn_R_V, k_bind_M1, k_rel = dynPar 

	if observableId == "RM5" 
		noiseParameter1_RM5 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_RM5 
	end

	if observableId == "RVSegTot" 
		noiseParameter1_RVSegTot = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_RVSegTot 
	end

	if observableId == "RCSegTot" 
		noiseParameter1_RCSegTot = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_RCSegTot 
	end

	if observableId == "Vrel" 
		noiseParameter1_Vrel = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_Vrel 
	end

	if observableId == "IntNucOffset" 
		noiseParameter1_IntNucOffset = getObsOrSdParam(sdPar, mapSdParam)
		return 13.0334748500416 * noiseParameter1_IntNucOffset 
	end

	if observableId == "FracNucInt_1" 
		noiseParameter1_FracNucInt_1 = getObsOrSdParam(sdPar, mapSdParam)
		return 3.10956052629092 * noiseParameter1_FracNucInt_1 
	end

	if observableId == "FracNucInt_2" 
		noiseParameter1_FracNucInt_2 = getObsOrSdParam(sdPar, mapSdParam)
		return 3.44650910342625 * noiseParameter1_FracNucInt_2 
	end

	if observableId == "FracNucInt_3" 
		noiseParameter1_FracNucInt_3 = getObsOrSdParam(sdPar, mapSdParam)
		return 1.96540708251497 * noiseParameter1_FracNucInt_3 
	end

	if observableId == "FracNucInt_4" 
		noiseParameter1_FracNucInt_4 = getObsOrSdParam(sdPar, mapSdParam)
		return 3.03676774438437 * noiseParameter1_FracNucInt_4 
	end

	if observableId == "FracNucInt_5" 
		noiseParameter1_FracNucInt_5 = getObsOrSdParam(sdPar, mapSdParam)
		return 4.00311961683219 * noiseParameter1_FracNucInt_5 
	end

	if observableId == "FracNucInt_6" 
		noiseParameter1_FracNucInt_6 = getObsOrSdParam(sdPar, mapSdParam)
		return 4.8783911282307 * noiseParameter1_FracNucInt_6 
	end

	if observableId == "FracNucInt_7" 
		noiseParameter1_FracNucInt_7 = getObsOrSdParam(sdPar, mapSdParam)
		return 3.96353798182046 * noiseParameter1_FracNucInt_7 
	end

	if observableId == "FracNucInt_8" 
		noiseParameter1_FracNucInt_8 = getObsOrSdParam(sdPar, mapSdParam)
		return 4.8820897165046 * noiseParameter1_FracNucInt_8 
	end

end