function evalYmod(u, t, dynPar, obsPar, nonDynParam, paramData, observableId, mapObsParam) 

	pAC, Rp8_Br_cAMPS, Rp8_pCPT_cAMPS, PDE, Rp_cAMPS, RII_2, RIIp_Rp8_Br_cAMPS_C_2, cAMP, RIIp_Sp8_Br_cAMPS_C_2, IBMX, AC_Fsk, RIIp_C_2, Sp8_Br_cAMPS, RII_C_2, RIIp_Sp8_Br_cAMPS_2, Csub, Csub_H89, AC, RIIp_cAMP_2, pAC_Fsk, RIIp_Rp8_pCPT_cAMPS_C_2, pPDE, RIIp_Rp_cAMPS_C_2, RIIp_2, RIIp_cAMP_C_2= u 
	KD_Fsk, KD_H89, KD_IBMX, KD_cAMP, kdeg_cAMP, kdeg_cAMP_free, kf_Fsk, kf_H89, kf_RII_2__RII_C_2, kf_RII_C_2__RII_2, kf_RII_C_2__RIIp_C_2, kf_RIIp_2__RII_2, kf_RIIp_C_2__RII_C_2, kf_RIIp_cAMP_C_2__RIIp_2, kf_cAMP, ki_IBMX, ki_Rp8_Br_cAMPS_pAB, ki_Rp8_pCPT_cAMPS_pAB, ki_Rp_cAMPS_pAB, ki_Sp8_Br_cAMPS_AM, ks_AC_cAMP, xi_AC_cAMP_Fsk, xi_KD_Rp8_Br_cAMPS, xi_KD_Rp8_pCPT_cAMPS, xi_KD_Rp_cAMPS, xi_KD_Sp8_Br_cAMPS, xi_b_Rp8_Br_cAMPS, xi_b_Rp8_pCPT_cAMPS, xi_b_Rp_cAMPS, xi_b_Sp8_Br_cAMPS, xi_kf_RII_2__RII_C_2, xi_kf_RII_C_2__RII_2 = dynPar 
	b_Calpha_global, b_pRII_global, rel_open, s_Calpha_global, s_pRII_Western, s_pRII_global, rho_Calpha_Microscopy, rho_pRII_Microscopy, rho_pRII_Western, xi_rel_open = nonDynParam 
	AC_total_C = paramData.nominalValue[1] 
	KD_PDE_Csub_C = paramData.nominalValue[5] 
	PDE_total_C = paramData.nominalValue[7] 
	RII2_total_C = paramData.nominalValue[8] 
	kdp_AC_C = paramData.nominalValue[13] 
	kf_PDE_Csub_C = paramData.nominalValue[16] 
	kp_AC_C = paramData.nominalValue[29] 
	s_Calpha_JI09_160201_Drg453_452_CycNuc_C = paramData.nominalValue[32] 
	xi_i_Rp8_Br_cAMPS_pAB_C = paramData.nominalValue[52] 
	xi_i_Rp8_pCPT_cAMPS_pAB_C = paramData.nominalValue[53] 
	xi_i_Rp_cAMPS_pAB_C = paramData.nominalValue[54] 
	xi_i_Sp8_Br_cAMPS_AM_C = paramData.nominalValue[55] 
	xi_pAC_C = paramData.nominalValue[58] 
	xi_pPDE_C = paramData.nominalValue[59] 

	if observableId == :pRII_Microscopy 
		observableParameter1_pRII_Microscopy = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pRII_Microscopy * ( b_pRII_global + s_pRII_global * ( 2 * RIIp_2 + 2 * RIIp_cAMP_2 + 2 * RIIp_Sp8_Br_cAMPS_2 + 2 * rel_open * ( RIIp_C_2 + RIIp_Rp_cAMPS_C_2 + RIIp_Rp8_Br_cAMPS_C_2 + RIIp_Rp8_pCPT_cAMPS_C_2 ) + 2 * ( RIIp_cAMP_C_2 + RIIp_Sp8_Br_cAMPS_C_2 ) * ( rel_open - xi_rel_open * ( rel_open - 1 ) ) ) ) 
	end

	if observableId == :pRII_Western 
		return s_pRII_Western * ( 2 * RIIp_2 + 2 * RIIp_C_2 + 2 * RIIp_cAMP_2 + 2 * RIIp_cAMP_C_2 + 2 * RIIp_Rp_cAMPS_C_2 + 2 * RIIp_Sp8_Br_cAMPS_2 + 2 * RIIp_Rp8_Br_cAMPS_C_2 + 2 * RIIp_Sp8_Br_cAMPS_C_2 + 2 * RIIp_Rp8_pCPT_cAMPS_C_2 ) 
	end

	if observableId == :Calpha_Microscopy 
		observableParameter1_Calpha_Microscopy = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_Calpha_Microscopy * ( b_Calpha_global + s_Calpha_global * ( 2 * Csub + 2 * Csub_H89 + 2 * rel_open * ( RIIp_C_2 + RIIp_Rp_cAMPS_C_2 + RIIp_Rp8_Br_cAMPS_C_2 + RIIp_Rp8_pCPT_cAMPS_C_2 ) + 2 * ( RIIp_cAMP_C_2 + RIIp_Sp8_Br_cAMPS_C_2 ) * ( rel_open - xi_rel_open * ( rel_open - 1 ) ) ) ) 
	end

end

function evalU0!(u0Vec, paramVec) 

	ki_Rp8_pCPT_cAMPS_pAB, xi_b_Rp_cAMPS, RII2_total, H89_level, fourABnOH_incubation_time, KD_Fsk, kdeg_cAMP_free, Rp8_Br_cAMPS_pAB_level, xi_KD_Rp8_Br_cAMPS, kf_PDE_Csub, Sp8_Br_cAMPS_AM_level, xi_kf_RII_C_2__RII_2, kf_RIIp_2__RII_2, kf_cAMP, IBMX_time, H89_bool1, Rp_cAMPS_pAB_incubation_time, kf_H89, kf_RII_C_2__RII_2, kf_RII_C_2__RIIp_C_2, xi_i_Rp8_pCPT_cAMPS_pAB, xi_pAC, fourABnOH_level, ki_Sp8_Br_cAMPS_AM, Rp8_pCPT_cAMPS_pAB_bool1, xi_b_Sp8_Br_cAMPS, xi_b_Rp8_Br_cAMPS, H89_time, Rp_cAMPS_pAB_bool1, kdeg_cAMP, xi_AC_cAMP_Fsk, xi_b_Rp8_pCPT_cAMPS, ki_Rp_cAMPS_pAB, xi_i_Rp_cAMPS_pAB, KD_PDE_Csub, ki_IBMX, IBMXex_bool1, Fsk_time, Rp8_pCPT_cAMPS_pAB_incubation_time, PDE_total, Rp8_Br_cAMPS_pAB_bool1, ki_Rp8_Br_cAMPS_pAB, xi_kf_RII_2__RII_C_2, default, xi_i_Rp8_Br_cAMPS_pAB, Sp8_Br_cAMPS_AM_time, xi_KD_Rp8_pCPT_cAMPS, kp_AC, xi_pPDE, xi_KD_Rp_cAMPS, Rp8_Br_cAMPS_pAB_incubation_time, AC_total, kf_RIIp_C_2__RII_C_2, xi_KD_Sp8_Br_cAMPS, nuc, Sp8_Br_cAMPS_AM_bool1, kf_RIIp_cAMP_C_2__RIIp_2, KD_cAMP, KD_IBMX, kf_Fsk, xi_i_Sp8_Br_cAMPS_AM, cyt, KD_H89, ks_AC_cAMP, fourABnOH_bool1, Rp8_pCPT_cAMPS_pAB_level, Rp_cAMPS_pAB_level, kf_RII_2__RII_C_2, Fsk_level, Fsk_bool1, kdp_AC, IBMX_level = paramVec 

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

	u0Vec .= pAC, Rp8_Br_cAMPS, Rp8_pCPT_cAMPS, PDE, Rp_cAMPS, RII_2, RIIp_Rp8_Br_cAMPS_C_2, cAMP, RIIp_Sp8_Br_cAMPS_C_2, IBMX, AC_Fsk, RIIp_C_2, Sp8_Br_cAMPS, RII_C_2, RIIp_Sp8_Br_cAMPS_2, Csub, Csub_H89, AC, RIIp_cAMP_2, pAC_Fsk, RIIp_Rp8_pCPT_cAMPS_C_2, pPDE, RIIp_Rp_cAMPS_C_2, RIIp_2, RIIp_cAMP_C_2
end

function evalU0(paramVec) 

	ki_Rp8_pCPT_cAMPS_pAB, xi_b_Rp_cAMPS, RII2_total, H89_level, fourABnOH_incubation_time, KD_Fsk, kdeg_cAMP_free, Rp8_Br_cAMPS_pAB_level, xi_KD_Rp8_Br_cAMPS, kf_PDE_Csub, Sp8_Br_cAMPS_AM_level, xi_kf_RII_C_2__RII_2, kf_RIIp_2__RII_2, kf_cAMP, IBMX_time, H89_bool1, Rp_cAMPS_pAB_incubation_time, kf_H89, kf_RII_C_2__RII_2, kf_RII_C_2__RIIp_C_2, xi_i_Rp8_pCPT_cAMPS_pAB, xi_pAC, fourABnOH_level, ki_Sp8_Br_cAMPS_AM, Rp8_pCPT_cAMPS_pAB_bool1, xi_b_Sp8_Br_cAMPS, xi_b_Rp8_Br_cAMPS, H89_time, Rp_cAMPS_pAB_bool1, kdeg_cAMP, xi_AC_cAMP_Fsk, xi_b_Rp8_pCPT_cAMPS, ki_Rp_cAMPS_pAB, xi_i_Rp_cAMPS_pAB, KD_PDE_Csub, ki_IBMX, IBMXex_bool1, Fsk_time, Rp8_pCPT_cAMPS_pAB_incubation_time, PDE_total, Rp8_Br_cAMPS_pAB_bool1, ki_Rp8_Br_cAMPS_pAB, xi_kf_RII_2__RII_C_2, default, xi_i_Rp8_Br_cAMPS_pAB, Sp8_Br_cAMPS_AM_time, xi_KD_Rp8_pCPT_cAMPS, kp_AC, xi_pPDE, xi_KD_Rp_cAMPS, Rp8_Br_cAMPS_pAB_incubation_time, AC_total, kf_RIIp_C_2__RII_C_2, xi_KD_Sp8_Br_cAMPS, nuc, Sp8_Br_cAMPS_AM_bool1, kf_RIIp_cAMP_C_2__RIIp_2, KD_cAMP, KD_IBMX, kf_Fsk, xi_i_Sp8_Br_cAMPS_AM, cyt, KD_H89, ks_AC_cAMP, fourABnOH_bool1, Rp8_pCPT_cAMPS_pAB_level, Rp_cAMPS_pAB_level, kf_RII_2__RII_C_2, Fsk_level, Fsk_bool1, kdp_AC, IBMX_level = paramVec 

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

	 return [pAC, Rp8_Br_cAMPS, Rp8_pCPT_cAMPS, PDE, Rp_cAMPS, RII_2, RIIp_Rp8_Br_cAMPS_C_2, cAMP, RIIp_Sp8_Br_cAMPS_C_2, IBMX, AC_Fsk, RIIp_C_2, Sp8_Br_cAMPS, RII_C_2, RIIp_Sp8_Br_cAMPS_2, Csub, Csub_H89, AC, RIIp_cAMP_2, pAC_Fsk, RIIp_Rp8_pCPT_cAMPS_C_2, pPDE, RIIp_Rp_cAMPS_C_2, RIIp_2, RIIp_cAMP_C_2]
end

function evalSd!(u, t, sdPar, dynPar, nonDynParam, paramData, observableId, mapSdParam) 

	pAC, Rp8_Br_cAMPS, Rp8_pCPT_cAMPS, PDE, Rp_cAMPS, RII_2, RIIp_Rp8_Br_cAMPS_C_2, cAMP, RIIp_Sp8_Br_cAMPS_C_2, IBMX, AC_Fsk, RIIp_C_2, Sp8_Br_cAMPS, RII_C_2, RIIp_Sp8_Br_cAMPS_2, Csub, Csub_H89, AC, RIIp_cAMP_2, pAC_Fsk, RIIp_Rp8_pCPT_cAMPS_C_2, pPDE, RIIp_Rp_cAMPS_C_2, RIIp_2, RIIp_cAMP_C_2= u 
	KD_Fsk, KD_H89, KD_IBMX, KD_cAMP, kdeg_cAMP, kdeg_cAMP_free, kf_Fsk, kf_H89, kf_RII_2__RII_C_2, kf_RII_C_2__RII_2, kf_RII_C_2__RIIp_C_2, kf_RIIp_2__RII_2, kf_RIIp_C_2__RII_C_2, kf_RIIp_cAMP_C_2__RIIp_2, kf_cAMP, ki_IBMX, ki_Rp8_Br_cAMPS_pAB, ki_Rp8_pCPT_cAMPS_pAB, ki_Rp_cAMPS_pAB, ki_Sp8_Br_cAMPS_AM, ks_AC_cAMP, xi_AC_cAMP_Fsk, xi_KD_Rp8_Br_cAMPS, xi_KD_Rp8_pCPT_cAMPS, xi_KD_Rp_cAMPS, xi_KD_Sp8_Br_cAMPS, xi_b_Rp8_Br_cAMPS, xi_b_Rp8_pCPT_cAMPS, xi_b_Rp_cAMPS, xi_b_Sp8_Br_cAMPS, xi_kf_RII_2__RII_C_2, xi_kf_RII_C_2__RII_2 = dynPar 
	b_Calpha_global, b_pRII_global, rel_open, s_Calpha_global, s_pRII_Western, s_pRII_global, rho_Calpha_Microscopy, rho_pRII_Microscopy, rho_pRII_Western, xi_rel_open = nonDynParam 
	AC_total_C = paramData.nominalValue[1] 
	KD_PDE_Csub_C = paramData.nominalValue[5] 
	PDE_total_C = paramData.nominalValue[7] 
	RII2_total_C = paramData.nominalValue[8] 
	kdp_AC_C = paramData.nominalValue[13] 
	kf_PDE_Csub_C = paramData.nominalValue[16] 
	kp_AC_C = paramData.nominalValue[29] 
	s_Calpha_JI09_160201_Drg453_452_CycNuc_C = paramData.nominalValue[32] 
	xi_i_Rp8_Br_cAMPS_pAB_C = paramData.nominalValue[52] 
	xi_i_Rp8_pCPT_cAMPS_pAB_C = paramData.nominalValue[53] 
	xi_i_Rp_cAMPS_pAB_C = paramData.nominalValue[54] 
	xi_i_Sp8_Br_cAMPS_AM_C = paramData.nominalValue[55] 
	xi_pAC_C = paramData.nominalValue[58] 
	xi_pPDE_C = paramData.nominalValue[59] 

	if observableId == :pRII_Microscopy 
		noiseParameter1_pRII_Microscopy = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pRII_Microscopy * s_pRII_global * rho_pRII_Microscopy 
	end

	if observableId == :pRII_Western 
		return s_pRII_Western * rho_pRII_Western 
	end

	if observableId == :Calpha_Microscopy 
		return s_Calpha_global * rho_Calpha_Microscopy 
	end

end