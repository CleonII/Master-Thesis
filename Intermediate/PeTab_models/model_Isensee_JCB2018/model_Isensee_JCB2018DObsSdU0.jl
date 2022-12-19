function evalDYmodDu(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	pAC, Rp8_Br_cAMPS, Rp8_pCPT_cAMPS, PDE, Rp_cAMPS, RII_2, RIIp_Rp8_Br_cAMPS_C_2, cAMP, RIIp_Sp8_Br_cAMPS_C_2, IBMX, AC_Fsk, RIIp_C_2, Sp8_Br_cAMPS, RII_C_2, RIIp_Sp8_Br_cAMPS_2, Csub, Csub_H89, AC, RIIp_cAMP_2, pAC_Fsk, RIIp_Rp8_pCPT_cAMPS_C_2, pPDE, RIIp_Rp_cAMPS_C_2, RIIp_2, RIIp_cAMP_C_2= u 
	ki_Rp8_pCPT_cAMPS_pAB, xi_b_Rp_cAMPS, RII2_total, H89_level, fourABnOH_incubation_time, KD_Fsk, kdeg_cAMP_free, Rp8_Br_cAMPS_pAB_level, xi_KD_Rp8_Br_cAMPS, kf_PDE_Csub, Sp8_Br_cAMPS_AM_level, xi_kf_RII_C_2__RII_2, kf_RIIp_2__RII_2, kf_cAMP, IBMX_time, Rp_cAMPS_pAB_incubation_time, kf_H89, kf_RII_C_2__RII_2, kf_RII_C_2__RIIp_C_2, xi_i_Rp8_pCPT_cAMPS_pAB, xi_pAC, fourABnOH_level, ki_Sp8_Br_cAMPS_AM, xi_b_Sp8_Br_cAMPS, xi_b_Rp8_Br_cAMPS, H89_time, kdeg_cAMP, xi_AC_cAMP_Fsk, xi_b_Rp8_pCPT_cAMPS, ki_Rp_cAMPS_pAB, xi_i_Rp_cAMPS_pAB, KD_PDE_Csub, ki_IBMX, Fsk_time, Rp8_pCPT_cAMPS_pAB_incubation_time, PDE_total, ki_Rp8_Br_cAMPS_pAB, xi_kf_RII_2__RII_C_2, default, xi_i_Rp8_Br_cAMPS_pAB, Sp8_Br_cAMPS_AM_time, xi_KD_Rp8_pCPT_cAMPS, kp_AC, xi_pPDE, xi_KD_Rp_cAMPS, Rp8_Br_cAMPS_pAB_incubation_time, AC_total, kf_RIIp_C_2__RII_C_2, xi_KD_Sp8_Br_cAMPS, nuc, kf_RIIp_cAMP_C_2__RIIp_2, KD_cAMP, KD_IBMX, kf_Fsk, xi_i_Sp8_Br_cAMPS_AM, cyt, KD_H89, ks_AC_cAMP, Rp8_pCPT_cAMPS_pAB_level, Rp_cAMPS_pAB_level, kf_RII_2__RII_C_2, Fsk_level, kdp_AC, IBMX_level = p 
	b_Calpha_global, b_pRII_global, rel_open, s_Calpha_global, s_pRII_Western, s_pRII_global, rho_Calpha_Microscopy, rho_pRII_Microscopy, rho_pRII_Western, xi_rel_open = nonDynParam 
	if observableId == "pRII_Microscopy" 
		observableParameter1_pRII_Microscopy = getObsOrSdParam(obsPar, mapObsParam)
		out[7] = 2observableParameter1_pRII_Microscopy*rel_open*s_pRII_global
		out[9] = observableParameter1_pRII_Microscopy*s_pRII_global*(2rel_open - 2xi_rel_open*(rel_open - 1))
		out[12] = 2observableParameter1_pRII_Microscopy*rel_open*s_pRII_global
		out[15] = 2observableParameter1_pRII_Microscopy*s_pRII_global
		out[19] = 2observableParameter1_pRII_Microscopy*s_pRII_global
		out[21] = 2observableParameter1_pRII_Microscopy*rel_open*s_pRII_global
		out[23] = 2observableParameter1_pRII_Microscopy*rel_open*s_pRII_global
		out[24] = 2observableParameter1_pRII_Microscopy*s_pRII_global
		out[25] = observableParameter1_pRII_Microscopy*s_pRII_global*(2rel_open - 2xi_rel_open*(rel_open - 1))
		return nothing
	end

	if observableId == "pRII_Western" 
		out[7] = 2s_pRII_Western
		out[9] = 2s_pRII_Western
		out[12] = 2s_pRII_Western
		out[15] = 2s_pRII_Western
		out[19] = 2s_pRII_Western
		out[21] = 2s_pRII_Western
		out[23] = 2s_pRII_Western
		out[24] = 2s_pRII_Western
		out[25] = 2s_pRII_Western
		return nothing
	end

	if observableId == "Calpha_Microscopy" 
		observableParameter1_Calpha_Microscopy = getObsOrSdParam(obsPar, mapObsParam)
		out[7] = 2observableParameter1_Calpha_Microscopy*rel_open*s_Calpha_global
		out[9] = observableParameter1_Calpha_Microscopy*s_Calpha_global*(2rel_open - 2xi_rel_open*(rel_open - 1))
		out[12] = 2observableParameter1_Calpha_Microscopy*rel_open*s_Calpha_global
		out[16] = 2observableParameter1_Calpha_Microscopy*s_Calpha_global
		out[17] = 2observableParameter1_Calpha_Microscopy*s_Calpha_global
		out[21] = 2observableParameter1_Calpha_Microscopy*rel_open*s_Calpha_global
		out[23] = 2observableParameter1_Calpha_Microscopy*rel_open*s_Calpha_global
		out[25] = observableParameter1_Calpha_Microscopy*s_Calpha_global*(2rel_open - 2xi_rel_open*(rel_open - 1))
		return nothing
	end

end

function evalDYmodDp(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	pAC, Rp8_Br_cAMPS, Rp8_pCPT_cAMPS, PDE, Rp_cAMPS, RII_2, RIIp_Rp8_Br_cAMPS_C_2, cAMP, RIIp_Sp8_Br_cAMPS_C_2, IBMX, AC_Fsk, RIIp_C_2, Sp8_Br_cAMPS, RII_C_2, RIIp_Sp8_Br_cAMPS_2, Csub, Csub_H89, AC, RIIp_cAMP_2, pAC_Fsk, RIIp_Rp8_pCPT_cAMPS_C_2, pPDE, RIIp_Rp_cAMPS_C_2, RIIp_2, RIIp_cAMP_C_2= u 
	ki_Rp8_pCPT_cAMPS_pAB, xi_b_Rp_cAMPS, RII2_total, H89_level, fourABnOH_incubation_time, KD_Fsk, kdeg_cAMP_free, Rp8_Br_cAMPS_pAB_level, xi_KD_Rp8_Br_cAMPS, kf_PDE_Csub, Sp8_Br_cAMPS_AM_level, xi_kf_RII_C_2__RII_2, kf_RIIp_2__RII_2, kf_cAMP, IBMX_time, Rp_cAMPS_pAB_incubation_time, kf_H89, kf_RII_C_2__RII_2, kf_RII_C_2__RIIp_C_2, xi_i_Rp8_pCPT_cAMPS_pAB, xi_pAC, fourABnOH_level, ki_Sp8_Br_cAMPS_AM, xi_b_Sp8_Br_cAMPS, xi_b_Rp8_Br_cAMPS, H89_time, kdeg_cAMP, xi_AC_cAMP_Fsk, xi_b_Rp8_pCPT_cAMPS, ki_Rp_cAMPS_pAB, xi_i_Rp_cAMPS_pAB, KD_PDE_Csub, ki_IBMX, Fsk_time, Rp8_pCPT_cAMPS_pAB_incubation_time, PDE_total, ki_Rp8_Br_cAMPS_pAB, xi_kf_RII_2__RII_C_2, default, xi_i_Rp8_Br_cAMPS_pAB, Sp8_Br_cAMPS_AM_time, xi_KD_Rp8_pCPT_cAMPS, kp_AC, xi_pPDE, xi_KD_Rp_cAMPS, Rp8_Br_cAMPS_pAB_incubation_time, AC_total, kf_RIIp_C_2__RII_C_2, xi_KD_Sp8_Br_cAMPS, nuc, kf_RIIp_cAMP_C_2__RIIp_2, KD_cAMP, KD_IBMX, kf_Fsk, xi_i_Sp8_Br_cAMPS_AM, cyt, KD_H89, ks_AC_cAMP, Rp8_pCPT_cAMPS_pAB_level, Rp_cAMPS_pAB_level, kf_RII_2__RII_C_2, Fsk_level, kdp_AC, IBMX_level = p 
	b_Calpha_global, b_pRII_global, rel_open, s_Calpha_global, s_pRII_Western, s_pRII_global, rho_Calpha_Microscopy, rho_pRII_Microscopy, rho_pRII_Western, xi_rel_open = nonDynParam 
	if observableId == "pRII_Microscopy" 
		return nothing
	end

	if observableId == "pRII_Western" 
		return nothing
	end

	if observableId == "Calpha_Microscopy" 
		return nothing
	end

end

function evalDSdDu!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	pAC, Rp8_Br_cAMPS, Rp8_pCPT_cAMPS, PDE, Rp_cAMPS, RII_2, RIIp_Rp8_Br_cAMPS_C_2, cAMP, RIIp_Sp8_Br_cAMPS_C_2, IBMX, AC_Fsk, RIIp_C_2, Sp8_Br_cAMPS, RII_C_2, RIIp_Sp8_Br_cAMPS_2, Csub, Csub_H89, AC, RIIp_cAMP_2, pAC_Fsk, RIIp_Rp8_pCPT_cAMPS_C_2, pPDE, RIIp_Rp_cAMPS_C_2, RIIp_2, RIIp_cAMP_C_2= u 
	ki_Rp8_pCPT_cAMPS_pAB, xi_b_Rp_cAMPS, RII2_total, H89_level, fourABnOH_incubation_time, KD_Fsk, kdeg_cAMP_free, Rp8_Br_cAMPS_pAB_level, xi_KD_Rp8_Br_cAMPS, kf_PDE_Csub, Sp8_Br_cAMPS_AM_level, xi_kf_RII_C_2__RII_2, kf_RIIp_2__RII_2, kf_cAMP, IBMX_time, Rp_cAMPS_pAB_incubation_time, kf_H89, kf_RII_C_2__RII_2, kf_RII_C_2__RIIp_C_2, xi_i_Rp8_pCPT_cAMPS_pAB, xi_pAC, fourABnOH_level, ki_Sp8_Br_cAMPS_AM, xi_b_Sp8_Br_cAMPS, xi_b_Rp8_Br_cAMPS, H89_time, kdeg_cAMP, xi_AC_cAMP_Fsk, xi_b_Rp8_pCPT_cAMPS, ki_Rp_cAMPS_pAB, xi_i_Rp_cAMPS_pAB, KD_PDE_Csub, ki_IBMX, Fsk_time, Rp8_pCPT_cAMPS_pAB_incubation_time, PDE_total, ki_Rp8_Br_cAMPS_pAB, xi_kf_RII_2__RII_C_2, default, xi_i_Rp8_Br_cAMPS_pAB, Sp8_Br_cAMPS_AM_time, xi_KD_Rp8_pCPT_cAMPS, kp_AC, xi_pPDE, xi_KD_Rp_cAMPS, Rp8_Br_cAMPS_pAB_incubation_time, AC_total, kf_RIIp_C_2__RII_C_2, xi_KD_Sp8_Br_cAMPS, nuc, kf_RIIp_cAMP_C_2__RIIp_2, KD_cAMP, KD_IBMX, kf_Fsk, xi_i_Sp8_Br_cAMPS_AM, cyt, KD_H89, ks_AC_cAMP, Rp8_pCPT_cAMPS_pAB_level, Rp_cAMPS_pAB_level, kf_RII_2__RII_C_2, Fsk_level, kdp_AC, IBMX_level = p 
	b_Calpha_global, b_pRII_global, rel_open, s_Calpha_global, s_pRII_Western, s_pRII_global, rho_Calpha_Microscopy, rho_pRII_Microscopy, rho_pRII_Western, xi_rel_open = nonDynParam 
	if observableId == "pRII_Microscopy" 
		return nothing
	end

	if observableId == "pRII_Western" 
		return nothing
	end

	if observableId == "Calpha_Microscopy" 
		return nothing
	end

end

function evalDSdDp!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	pAC, Rp8_Br_cAMPS, Rp8_pCPT_cAMPS, PDE, Rp_cAMPS, RII_2, RIIp_Rp8_Br_cAMPS_C_2, cAMP, RIIp_Sp8_Br_cAMPS_C_2, IBMX, AC_Fsk, RIIp_C_2, Sp8_Br_cAMPS, RII_C_2, RIIp_Sp8_Br_cAMPS_2, Csub, Csub_H89, AC, RIIp_cAMP_2, pAC_Fsk, RIIp_Rp8_pCPT_cAMPS_C_2, pPDE, RIIp_Rp_cAMPS_C_2, RIIp_2, RIIp_cAMP_C_2= u 
	ki_Rp8_pCPT_cAMPS_pAB, xi_b_Rp_cAMPS, RII2_total, H89_level, fourABnOH_incubation_time, KD_Fsk, kdeg_cAMP_free, Rp8_Br_cAMPS_pAB_level, xi_KD_Rp8_Br_cAMPS, kf_PDE_Csub, Sp8_Br_cAMPS_AM_level, xi_kf_RII_C_2__RII_2, kf_RIIp_2__RII_2, kf_cAMP, IBMX_time, Rp_cAMPS_pAB_incubation_time, kf_H89, kf_RII_C_2__RII_2, kf_RII_C_2__RIIp_C_2, xi_i_Rp8_pCPT_cAMPS_pAB, xi_pAC, fourABnOH_level, ki_Sp8_Br_cAMPS_AM, xi_b_Sp8_Br_cAMPS, xi_b_Rp8_Br_cAMPS, H89_time, kdeg_cAMP, xi_AC_cAMP_Fsk, xi_b_Rp8_pCPT_cAMPS, ki_Rp_cAMPS_pAB, xi_i_Rp_cAMPS_pAB, KD_PDE_Csub, ki_IBMX, Fsk_time, Rp8_pCPT_cAMPS_pAB_incubation_time, PDE_total, ki_Rp8_Br_cAMPS_pAB, xi_kf_RII_2__RII_C_2, default, xi_i_Rp8_Br_cAMPS_pAB, Sp8_Br_cAMPS_AM_time, xi_KD_Rp8_pCPT_cAMPS, kp_AC, xi_pPDE, xi_KD_Rp_cAMPS, Rp8_Br_cAMPS_pAB_incubation_time, AC_total, kf_RIIp_C_2__RII_C_2, xi_KD_Sp8_Br_cAMPS, nuc, kf_RIIp_cAMP_C_2__RIIp_2, KD_cAMP, KD_IBMX, kf_Fsk, xi_i_Sp8_Br_cAMPS_AM, cyt, KD_H89, ks_AC_cAMP, Rp8_pCPT_cAMPS_pAB_level, Rp_cAMPS_pAB_level, kf_RII_2__RII_C_2, Fsk_level, kdp_AC, IBMX_level = p 
	b_Calpha_global, b_pRII_global, rel_open, s_Calpha_global, s_pRII_Western, s_pRII_global, rho_Calpha_Microscopy, rho_pRII_Microscopy, rho_pRII_Western, xi_rel_open = nonDynParam 
	if observableId == "pRII_Microscopy" 
		return nothing
	end

	if observableId == "pRII_Western" 
		return nothing
	end

	if observableId == "Calpha_Microscopy" 
		return nothing
	end

end

