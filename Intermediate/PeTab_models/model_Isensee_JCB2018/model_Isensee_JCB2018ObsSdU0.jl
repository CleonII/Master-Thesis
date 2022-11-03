function evalYmod(u, t, dynPar, obsPar, paramData, observableId, mapObsParam) 

	pAC, Rp8_Br_cAMPS, Rp8_pCPT_cAMPS, PDE, Rp_cAMPS, RII_2, RIIp_Rp8_Br_cAMPS_C_2, cAMP, RIIp_Sp8_Br_cAMPS_C_2, IBMX, AC_Fsk, RIIp_C_2, Sp8_Br_cAMPS, RII_C_2, RIIp_Sp8_Br_cAMPS_2, Csub, Csub_H89, AC, RIIp_cAMP_2, pAC_Fsk, RIIp_Rp8_pCPT_cAMPS_C_2, pPDE, RIIp_Rp_cAMPS_C_2, RIIp_2, RIIp_cAMP_C_2, dummyVariable= u 
	KD_Fsk, KD_H89, KD_IBMX, KD_cAMP, b_Calpha_global, b_pRII_global, kdeg_cAMP, kdeg_cAMP_free, kf_Fsk, kf_H89, kf_RII_2__RII_C_2, kf_RII_C_2__RII_2, kf_RII_C_2__RIIp_C_2, kf_RIIp_2__RII_2, kf_RIIp_C_2__RII_C_2, kf_RIIp_cAMP_C_2__RIIp_2, kf_cAMP, ki_IBMX, ki_Rp8_Br_cAMPS_pAB, ki_Rp8_pCPT_cAMPS_pAB, ki_Rp_cAMPS_pAB, ki_Sp8_Br_cAMPS_AM, ks_AC_cAMP, rel_open, s_Calpha_global, s_pRII_Western, s_pRII_global, rho_Calpha_Microscopy, rho_pRII_Microscopy, rho_pRII_Western, xi_AC_cAMP_Fsk, xi_KD_Rp8_Br_cAMPS, xi_KD_Rp8_pCPT_cAMPS, xi_KD_Rp_cAMPS, xi_KD_Sp8_Br_cAMPS, xi_b_Rp8_Br_cAMPS, xi_b_Rp8_pCPT_cAMPS, xi_b_Rp_cAMPS, xi_b_Sp8_Br_cAMPS, xi_kf_RII_2__RII_C_2, xi_kf_RII_C_2__RII_2, xi_rel_open = dynPar 
	AC_total_C = paramData.paramVal[1] 
	KD_PDE_Csub_C = paramData.paramVal[5] 
	PDE_total_C = paramData.paramVal[7] 
	RII2_total_C = paramData.paramVal[8] 
	kdp_AC_C = paramData.paramVal[13] 
	kf_PDE_Csub_C = paramData.paramVal[16] 
	kp_AC_C = paramData.paramVal[29] 
	s_Calpha_JI09_160201_Drg453_452_CycNuc_C = paramData.paramVal[32] 
	xi_i_Rp8_Br_cAMPS_pAB_C = paramData.paramVal[52] 
	xi_i_Rp8_pCPT_cAMPS_pAB_C = paramData.paramVal[53] 
	xi_i_Rp_cAMPS_pAB_C = paramData.paramVal[54] 
	xi_i_Sp8_Br_cAMPS_AM_C = paramData.paramVal[55] 
	xi_pAC_C = paramData.paramVal[58] 
	xi_pPDE_C = paramData.paramVal[59] 

	if observableId == "pRII_Microscopy" 
		observableParameter1_pRII_Microscopy = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pRII_Microscopy * ( b_pRII_global + s_pRII_global * ( 2 * RIIp_2 + 2 * RIIp_cAMP_2 + 2 * RIIp_Sp8_Br_cAMPS_2 + 2 * rel_open * ( RIIp_C_2 + RIIp_Rp_cAMPS_C_2 + RIIp_Rp8_Br_cAMPS_C_2 + RIIp_Rp8_pCPT_cAMPS_C_2 ) + 2 * ( RIIp_cAMP_C_2 + RIIp_Sp8_Br_cAMPS_C_2 ) * ( rel_open - xi_rel_open * ( rel_open - 1 ) ) ) ) 
	end

	if observableId == "pRII_Western" 
		return s_pRII_Western * ( 2 * RIIp_2 + 2 * RIIp_C_2 + 2 * RIIp_cAMP_2 + 2 * RIIp_cAMP_C_2 + 2 * RIIp_Rp_cAMPS_C_2 + 2 * RIIp_Sp8_Br_cAMPS_2 + 2 * RIIp_Rp8_Br_cAMPS_C_2 + 2 * RIIp_Sp8_Br_cAMPS_C_2 + 2 * RIIp_Rp8_pCPT_cAMPS_C_2 ) 
	end

	if observableId == "Calpha_Microscopy" 
		observableParameter1_Calpha_Microscopy = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_Calpha_Microscopy * ( b_Calpha_global + s_Calpha_global * ( 2 * Csub + 2 * Csub_H89 + 2 * rel_open * ( RIIp_C_2 + RIIp_Rp_cAMPS_C_2 + RIIp_Rp8_Br_cAMPS_C_2 + RIIp_Rp8_pCPT_cAMPS_C_2 ) + 2 * ( RIIp_cAMP_C_2 + RIIp_Sp8_Br_cAMPS_C_2 ) * ( rel_open - xi_rel_open * ( rel_open - 1 ) ) ) ) 
	end

end

function evalU0!(u0Vec, paramVec) 

	KD_Fsk, kf_Fsk, kp_AC, kdp_AC, xi_i_Rp8_Br_cAMPS_pAB, ki_Rp8_Br_cAMPS_pAB, ki_Rp8_pCPT_cAMPS_pAB, xi_i_Rp8_pCPT_cAMPS_pAB, KD_PDE_Csub, kf_PDE_Csub, xi_i_Rp_cAMPS_pAB, ki_Rp_cAMPS_pAB, kf_RIIp_2__RII_2, kf_RII_2__RII_C_2, kf_RII_C_2__RII_2, kf_cAMP, xi_b_Rp8_Br_cAMPS, KD_cAMP, xi_KD_Rp8_Br_cAMPS, xi_AC_cAMP_Fsk, kdeg_cAMP_free, KD_IBMX, xi_pAC, xi_pPDE, ks_AC_cAMP, xi_b_Sp8_Br_cAMPS, kf_RIIp_cAMP_C_2__RIIp_2, xi_KD_Sp8_Br_cAMPS, ki_IBMX, kf_RII_C_2__RIIp_C_2, xi_b_Rp8_pCPT_cAMPS, xi_KD_Rp8_pCPT_cAMPS, kf_RIIp_C_2__RII_C_2, xi_kf_RII_C_2__RII_2, xi_b_Rp_cAMPS, kdeg_cAMP, xi_kf_RII_2__RII_C_2, xi_KD_Rp_cAMPS, ki_Sp8_Br_cAMPS_AM, xi_i_Sp8_Br_cAMPS_AM, KD_H89, kf_H89, Rp_cAMPS_pAB_level, Rp_cAMPS_pAB_incubation_time, H89_time, H89_level, Fsk_time, Fsk_level, IBMX_time, IBMX_level, Rp8_Br_cAMPS_pAB_level, Rp8_Br_cAMPS_pAB_incubation_time, Rp8_pCPT_cAMPS_pAB_level, Rp8_pCPT_cAMPS_pAB_incubation_time, fourABnOH_level, fourABnOH_incubation_time, Sp8_Br_cAMPS_AM_time, Sp8_Br_cAMPS_AM_level, nuc, PDE_total, default, RII2_total, AC_total = paramVec 

	pAC = 0.0 
	Rp8_Br_cAMPS = 0.0 
	Rp8_pCPT_cAMPS = 0.0 
	PDE = 1.0 
	Rp_cAMPS = 0.0 
	RII_2 = 0.057671482854616 
	RIIp_Rp8_Br_cAMPS_C_2 = 0.0 
	cAMP = 0.0575218758977949 
	RIIp_Sp8_Br_cAMPS_C_2 = 0.0 
	IBMX = 0.0 
	AC_Fsk = 0.0 
	RIIp_C_2 = 0.287974352643203 
	Sp8_Br_cAMPS = 0.0 
	RII_C_2 = 0.514562213223039 
	RIIp_Sp8_Br_cAMPS_2 = 0.0 
	Csub = 0.192051405433538 
	Csub_H89 = 0.0 
	AC = 1.0 
	RIIp_cAMP_2 = 0.000400416634006059 
	pAC_Fsk = 0.0 
	RIIp_Rp8_pCPT_cAMPS_C_2 = 0.0 
	pPDE = 0.0 
	RIIp_Rp_cAMPS_C_2 = 0.0 
	RIIp_2 = 0.133979505944916 
	RIIp_cAMP_C_2 = 0.00541202870022029 
	dummyVariable = 0.0 

	u0Vec .= pAC, Rp8_Br_cAMPS, Rp8_pCPT_cAMPS, PDE, Rp_cAMPS, RII_2, RIIp_Rp8_Br_cAMPS_C_2, cAMP, RIIp_Sp8_Br_cAMPS_C_2, IBMX, AC_Fsk, RIIp_C_2, Sp8_Br_cAMPS, RII_C_2, RIIp_Sp8_Br_cAMPS_2, Csub, Csub_H89, AC, RIIp_cAMP_2, pAC_Fsk, RIIp_Rp8_pCPT_cAMPS_C_2, pPDE, RIIp_Rp_cAMPS_C_2, RIIp_2, RIIp_cAMP_C_2, dummyVariable
end

function evalSd!(u, t, sdPar, dynPar, paramData, observableId, mapSdParam) 

	pAC, Rp8_Br_cAMPS, Rp8_pCPT_cAMPS, PDE, Rp_cAMPS, RII_2, RIIp_Rp8_Br_cAMPS_C_2, cAMP, RIIp_Sp8_Br_cAMPS_C_2, IBMX, AC_Fsk, RIIp_C_2, Sp8_Br_cAMPS, RII_C_2, RIIp_Sp8_Br_cAMPS_2, Csub, Csub_H89, AC, RIIp_cAMP_2, pAC_Fsk, RIIp_Rp8_pCPT_cAMPS_C_2, pPDE, RIIp_Rp_cAMPS_C_2, RIIp_2, RIIp_cAMP_C_2, dummyVariable= u 
	KD_Fsk, KD_H89, KD_IBMX, KD_cAMP, b_Calpha_global, b_pRII_global, kdeg_cAMP, kdeg_cAMP_free, kf_Fsk, kf_H89, kf_RII_2__RII_C_2, kf_RII_C_2__RII_2, kf_RII_C_2__RIIp_C_2, kf_RIIp_2__RII_2, kf_RIIp_C_2__RII_C_2, kf_RIIp_cAMP_C_2__RIIp_2, kf_cAMP, ki_IBMX, ki_Rp8_Br_cAMPS_pAB, ki_Rp8_pCPT_cAMPS_pAB, ki_Rp_cAMPS_pAB, ki_Sp8_Br_cAMPS_AM, ks_AC_cAMP, rel_open, s_Calpha_global, s_pRII_Western, s_pRII_global, rho_Calpha_Microscopy, rho_pRII_Microscopy, rho_pRII_Western, xi_AC_cAMP_Fsk, xi_KD_Rp8_Br_cAMPS, xi_KD_Rp8_pCPT_cAMPS, xi_KD_Rp_cAMPS, xi_KD_Sp8_Br_cAMPS, xi_b_Rp8_Br_cAMPS, xi_b_Rp8_pCPT_cAMPS, xi_b_Rp_cAMPS, xi_b_Sp8_Br_cAMPS, xi_kf_RII_2__RII_C_2, xi_kf_RII_C_2__RII_2, xi_rel_open = dynPar 
	AC_total_C = paramData.paramVal[1] 
	KD_PDE_Csub_C = paramData.paramVal[5] 
	PDE_total_C = paramData.paramVal[7] 
	RII2_total_C = paramData.paramVal[8] 
	kdp_AC_C = paramData.paramVal[13] 
	kf_PDE_Csub_C = paramData.paramVal[16] 
	kp_AC_C = paramData.paramVal[29] 
	s_Calpha_JI09_160201_Drg453_452_CycNuc_C = paramData.paramVal[32] 
	xi_i_Rp8_Br_cAMPS_pAB_C = paramData.paramVal[52] 
	xi_i_Rp8_pCPT_cAMPS_pAB_C = paramData.paramVal[53] 
	xi_i_Rp_cAMPS_pAB_C = paramData.paramVal[54] 
	xi_i_Sp8_Br_cAMPS_AM_C = paramData.paramVal[55] 
	xi_pAC_C = paramData.paramVal[58] 
	xi_pPDE_C = paramData.paramVal[59] 

	if observableId == "pRII_Microscopy" 
		noiseParameter1_pRII_Microscopy = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pRII_Microscopy * s_pRII_global * rho_pRII_Microscopy 
	end

	if observableId == "pRII_Western" 
		return s_pRII_Western * rho_pRII_Western 
	end

	if observableId == "Calpha_Microscopy" 
		return s_Calpha_global * rho_Calpha_Microscopy 
	end

end