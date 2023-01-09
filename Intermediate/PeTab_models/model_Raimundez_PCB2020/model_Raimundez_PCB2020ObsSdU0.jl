function evalYmod(u, t, dynPar, obsPar, nonDynParam, paramData, observableId, mapObsParam) 

	EGFR, pPI3K, MMET_METinh_2, pMAPK, MMET_2, pEGFR_EGF_2, pEGFR_EGF_2_i, pMPI3K, EGFR_EGF_2, MMET, pMMET_2, RAS_GTP, pMMET_2_i, pMMET_pEGFR_i, MMET_MMET_METinh, pAKT, EGFR_EGF, MMET_METinh, MMET_EGFR, EGFR_CET, pMMET_pEGFR, EGFR_MMET_METinh= u 
	KD_EGFR_CET, KD_EGFR_EGF, KD_METinh, d_AKTtotal__MKN1_2_HS746T, d_AKTtotal__fm_2_hm, d_MAPKtotal__MKN1_2_HS746T, d_MAPKtotal__fm_2_hm, d_MPI3Ktotal__fm_2_hm, d_PI3Ktotal__fm_2_hm, d_RAStotal__MKN1_2_HS746T, d_RAStotal__fm_2_hm, d_kdeg_membran__MKN1_2_HS746T, d_kdeg_pEGFR_EGF_2_i__MKN1_2_HS746T, d_kexp_pEGFR_EGF_2_i__MKN1_2_HS746T, d_kimp_pEGFR_EGF_2__MKN1_2_HS746T, d_ksyn_EGFR__MKN1_2_HS746T, d_ksyn_EGFR__fm_2_hm, d_ksyn_MMET__fm_2_hm, ka_AKT__MKN1, ka_MAPK__MKN1, ka_PI3K__basal, ka_PI3K__pEGFR_EGF_2, ka_RAS__basal__MKN1, ka_RAS__pEGFR_EGF_2__MKN1, kbin_EGFR_CET, kbin_EGFR_EGF, kdeg_membran__MKN1, kdeg_pEGFR_EGF_2_i__MKN1, kdim_EGFR_EGF, kdim_MMET, kdim_MMET_EGFR, kdim_MMETinh, kexp_pEGFR_EGF_2_i__MKN1, ki_AKT__MKN1, ki_MAPK, ki_PI3K__MKN1, ki_RAS__MKN1, kimp_pEGFR_EGF_2__MKN1, kpho_EGFR_EGF, kpho_MMET, kpho_MMET_EGFR, ksyn_MMET__HS746T_fm, xi_ka_PI3K_pMMET_2, xi_ka_PI3K_pMMET_pEGFR, xi_ka_RAS_pMMET_2, xi_ka_RAS_pMMET_pEGFR, xi_kdeg_pMMET_2_i, xi_kdeg_pMMET_pEGFR_i, xi_kdim_MMET, xi_kdim_MMET_EGFR, xi_kexp_pMMET_2_i, xi_kexp_pMMET_pEGFR_i, xi_ki_MPI3K, xi_kimp_pMMET_2, xi_kimp_pMMET_pEGFR, xi_kpho_MMET, xi_kpho_MMET_EGFR = dynPar 
	AKTtotal__MKN1_fm_C = paramData.paramVal[1] 
	MAPKtotal__MKN1_fm_C = paramData.paramVal[5] 
	MPI3Ktotal__MKN1_fm_C = paramData.paramVal[6] 
	PI3Ktotal__HS746T_fm_C = paramData.paramVal[7] 
	RAStotal__MKN1_fm_C = paramData.paramVal[8] 
	d_ka_AKT__MKN1_2_HS746T_C = paramData.paramVal[17] 
	d_ka_MAPK__MKN1_2_HS746T_C = paramData.paramVal[18] 
	d_ka_RAS__basal__MKN1_2_HS746T_C = paramData.paramVal[19] 
	d_ka_RAS__pEGFR_EGF_2__MKN1_2_HS746T_C = paramData.paramVal[20] 
	d_ki_AKT__MKN1_2_HS746T_C = paramData.paramVal[24] 
	d_ki_PI3K__MKN1_2_HS746T_C = paramData.paramVal[25] 
	d_ki_RAS__MKN1_2_HS746T_C = paramData.paramVal[26] 
	ksyn_EGFR__MKN1_fm_C = paramData.paramVal[54] 

	if observableId == "EGFR_ID10_HS746T_FM_EGF_CET" 
		observableParameter1_EGFR_ID10_HS746T_FM_EGF_CET = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_EGFR_ID10_HS746T_FM_EGF_CET * ( EGFR + EGFR_EGF + EGFR_CET + 2 * EGFR_EGF_2 + EGFR_MMET_METinh + MMET_EGFR + 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "EGFR_ID10_HS746T_FM_EGF" 
		observableParameter1_EGFR_ID10_HS746T_FM_EGF = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_EGFR_ID10_HS746T_FM_EGF * ( EGFR + EGFR_EGF + EGFR_CET + 2 * EGFR_EGF_2 + EGFR_MMET_METinh + MMET_EGFR + 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "EGFR_ID12_MKN1_FM_HM" 
		observableParameter1_EGFR_ID12_MKN1_FM_HM = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_EGFR_ID12_MKN1_FM_HM * ( EGFR + EGFR_EGF + EGFR_CET + 2 * EGFR_EGF_2 + EGFR_MMET_METinh + MMET_EGFR + 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "EGFR_ID3_HS746T_FM" 
		observableParameter1_EGFR_ID3_HS746T_FM = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_EGFR_ID3_HS746T_FM * ( EGFR + EGFR_EGF + EGFR_CET + 2 * EGFR_EGF_2 + EGFR_MMET_METinh + MMET_EGFR + 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "EGFR_ID3_MKN1_FM" 
		observableParameter1_EGFR_ID3_MKN1_FM = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_EGFR_ID3_MKN1_FM * ( EGFR + EGFR_EGF + EGFR_CET + 2 * EGFR_EGF_2 + EGFR_MMET_METinh + MMET_EGFR + 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "EGFR_ID4_HS746T_HM_3min" 
		observableParameter1_EGFR_ID4_HS746T_HM_3min = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_EGFR_ID4_HS746T_HM_3min * ( EGFR + EGFR_EGF + EGFR_CET + 2 * EGFR_EGF_2 + EGFR_MMET_METinh + MMET_EGFR + 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "EGFR_ID4_MKN1_HM_3min" 
		observableParameter1_EGFR_ID4_MKN1_HM_3min = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_EGFR_ID4_MKN1_HM_3min * ( EGFR + EGFR_EGF + EGFR_CET + 2 * EGFR_EGF_2 + EGFR_MMET_METinh + MMET_EGFR + 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "MMET_ID10_HS746T_FM_EGF_CET" 
		observableParameter1_MMET_ID10_HS746T_FM_EGF_CET = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_MMET_ID10_HS746T_FM_EGF_CET * ( EGFR_MMET_METinh + MMET + 2 * MMET_2 + MMET_EGFR + MMET_METinh + 2 * MMET_METinh_2 + 2 * MMET_MMET_METinh + 2 * pMMET_2 + 2 * pMMET_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "MMET_ID10_HS746T_FM_EGF" 
		observableParameter1_MMET_ID10_HS746T_FM_EGF = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_MMET_ID10_HS746T_FM_EGF * ( EGFR_MMET_METinh + MMET + 2 * MMET_2 + MMET_EGFR + MMET_METinh + 2 * MMET_METinh_2 + 2 * MMET_MMET_METinh + 2 * pMMET_2 + 2 * pMMET_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pAKT_ID12_MKN1_FM_HM" 
		observableParameter1_pAKT_ID12_MKN1_FM_HM = getObsOrSdParam(obsPar, mapObsParam)
		return pAKT * observableParameter1_pAKT_ID12_MKN1_FM_HM 
	end

	if observableId == "pAKT_ID2_HS746T_HM" 
		observableParameter1_pAKT_ID2_HS746T_HM = getObsOrSdParam(obsPar, mapObsParam)
		return pAKT * observableParameter1_pAKT_ID2_HS746T_HM 
	end

	if observableId == "pAKT_ID2_MKN1_HM" 
		observableParameter1_pAKT_ID2_MKN1_HM = getObsOrSdParam(obsPar, mapObsParam)
		return pAKT * observableParameter1_pAKT_ID2_MKN1_HM 
	end

	if observableId == "pAKT_ID5_MKN1_FM_15min" 
		observableParameter1_pAKT_ID5_MKN1_FM_15min = getObsOrSdParam(obsPar, mapObsParam)
		return pAKT * observableParameter1_pAKT_ID5_MKN1_FM_15min 
	end

	if observableId == "pAKT_ID5_MKN1_FM_1min" 
		observableParameter1_pAKT_ID5_MKN1_FM_1min = getObsOrSdParam(obsPar, mapObsParam)
		return pAKT * observableParameter1_pAKT_ID5_MKN1_FM_1min 
	end

	if observableId == "pAKT_ID5_MKN1_FM_240min" 
		observableParameter1_pAKT_ID5_MKN1_FM_240min = getObsOrSdParam(obsPar, mapObsParam)
		return pAKT * observableParameter1_pAKT_ID5_MKN1_FM_240min 
	end

	if observableId == "pAKT_ID5_MKN1_FM_3min" 
		observableParameter1_pAKT_ID5_MKN1_FM_3min = getObsOrSdParam(obsPar, mapObsParam)
		return pAKT * observableParameter1_pAKT_ID5_MKN1_FM_3min 
	end

	if observableId == "pAKT_ID7_HS746T_FM_15min" 
		observableParameter1_pAKT_ID7_HS746T_FM_15min = getObsOrSdParam(obsPar, mapObsParam)
		return pAKT * observableParameter1_pAKT_ID7_HS746T_FM_15min 
	end

	if observableId == "pAKT_ID7_HS746T_FM_1min" 
		observableParameter1_pAKT_ID7_HS746T_FM_1min = getObsOrSdParam(obsPar, mapObsParam)
		return pAKT * observableParameter1_pAKT_ID7_HS746T_FM_1min 
	end

	if observableId == "pAKT_ID7_HS746T_FM_240min" 
		observableParameter1_pAKT_ID7_HS746T_FM_240min = getObsOrSdParam(obsPar, mapObsParam)
		return pAKT * observableParameter1_pAKT_ID7_HS746T_FM_240min 
	end

	if observableId == "pAKT_ID7_HS746T_FM_3min" 
		observableParameter1_pAKT_ID7_HS746T_FM_3min = getObsOrSdParam(obsPar, mapObsParam)
		return pAKT * observableParameter1_pAKT_ID7_HS746T_FM_3min 
	end

	if observableId == "pAKT_ID8_HS746T_HM_15min" 
		observableParameter1_pAKT_ID8_HS746T_HM_15min = getObsOrSdParam(obsPar, mapObsParam)
		return pAKT * observableParameter1_pAKT_ID8_HS746T_HM_15min 
	end

	if observableId == "pAKT_ID8_HS746T_HM_1min" 
		observableParameter1_pAKT_ID8_HS746T_HM_1min = getObsOrSdParam(obsPar, mapObsParam)
		return pAKT * observableParameter1_pAKT_ID8_HS746T_HM_1min 
	end

	if observableId == "pAKT_ID8_HS746T_HM_240min" 
		observableParameter1_pAKT_ID8_HS746T_HM_240min = getObsOrSdParam(obsPar, mapObsParam)
		return pAKT * observableParameter1_pAKT_ID8_HS746T_HM_240min 
	end

	if observableId == "pAKT_ID8_HS746T_HM_3min" 
		observableParameter1_pAKT_ID8_HS746T_HM_3min = getObsOrSdParam(obsPar, mapObsParam)
		return pAKT * observableParameter1_pAKT_ID8_HS746T_HM_3min 
	end

	if observableId == "pAKT_ID9_HS746T_FM" 
		observableParameter1_pAKT_ID9_HS746T_FM = getObsOrSdParam(obsPar, mapObsParam)
		return pAKT * observableParameter1_pAKT_ID9_HS746T_FM 
	end

	if observableId == "pAKT_ID9_HS746T_HM_30EGF" 
		observableParameter1_pAKT_ID9_HS746T_HM_30EGF = getObsOrSdParam(obsPar, mapObsParam)
		return pAKT * observableParameter1_pAKT_ID9_HS746T_HM_30EGF 
	end

	if observableId == "pAKT_ID9_HS746T_HM_5EGF" 
		observableParameter1_pAKT_ID9_HS746T_HM_5EGF = getObsOrSdParam(obsPar, mapObsParam)
		return pAKT * observableParameter1_pAKT_ID9_HS746T_HM_5EGF 
	end

	if observableId == "pEGFR_ID10_HS746T_FM_EGF_CET" 
		observableParameter1_pEGFR_ID10_HS746T_FM_EGF_CET = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_ID10_HS746T_FM_EGF_CET * ( 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pEGFR_ID10_HS746T_FM_EGF" 
		observableParameter1_pEGFR_ID10_HS746T_FM_EGF = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_ID10_HS746T_FM_EGF * ( 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pEGFR_ID12_MKN1_FM_HM" 
		observableParameter1_pEGFR_ID12_MKN1_FM_HM = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_ID12_MKN1_FM_HM * ( 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pEGFR_ID2_HS746T_HM" 
		observableParameter1_pEGFR_ID2_HS746T_HM = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_ID2_HS746T_HM * ( 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pEGFR_ID2_MKN1_HM" 
		observableParameter1_pEGFR_ID2_MKN1_HM = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_ID2_MKN1_HM * ( 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pEGFR_ID3_HS746T_FM" 
		observableParameter1_pEGFR_ID3_HS746T_FM = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_ID3_HS746T_FM * ( 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pEGFR_ID3_MKN1_FM" 
		observableParameter1_pEGFR_ID3_MKN1_FM = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_ID3_MKN1_FM * ( 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pEGFR_ID4_HS746T_HM_3min" 
		observableParameter1_pEGFR_ID4_HS746T_HM_3min = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_ID4_HS746T_HM_3min * ( 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pEGFR_ID4_MKN1_HM_3min" 
		observableParameter1_pEGFR_ID4_MKN1_HM_3min = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_ID4_MKN1_HM_3min * ( 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pEGFR_ID5_MKN1_FM_15min" 
		observableParameter1_pEGFR_ID5_MKN1_FM_15min = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_ID5_MKN1_FM_15min * ( 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pEGFR_ID5_MKN1_FM_1min" 
		observableParameter1_pEGFR_ID5_MKN1_FM_1min = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_ID5_MKN1_FM_1min * ( 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pEGFR_ID5_MKN1_FM_240min" 
		observableParameter1_pEGFR_ID5_MKN1_FM_240min = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_ID5_MKN1_FM_240min * ( 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pEGFR_ID5_MKN1_FM_3min" 
		observableParameter1_pEGFR_ID5_MKN1_FM_3min = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_ID5_MKN1_FM_3min * ( 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pEGFR_ID6_MKN1_HM_15min" 
		observableParameter1_pEGFR_ID6_MKN1_HM_15min = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_ID6_MKN1_HM_15min * ( 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pEGFR_ID6_MKN1_HM_1min" 
		observableParameter1_pEGFR_ID6_MKN1_HM_1min = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_ID6_MKN1_HM_1min * ( 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pEGFR_ID6_MKN1_HM_240min" 
		observableParameter1_pEGFR_ID6_MKN1_HM_240min = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_ID6_MKN1_HM_240min * ( 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pEGFR_ID6_MKN1_HM_3min" 
		observableParameter1_pEGFR_ID6_MKN1_HM_3min = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_ID6_MKN1_HM_3min * ( 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pEGFR_ID7_HS746T_FM_15min" 
		observableParameter1_pEGFR_ID7_HS746T_FM_15min = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_ID7_HS746T_FM_15min * ( 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pEGFR_ID7_HS746T_FM_1min" 
		observableParameter1_pEGFR_ID7_HS746T_FM_1min = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_ID7_HS746T_FM_1min * ( 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pEGFR_ID7_HS746T_FM_240min" 
		observableParameter1_pEGFR_ID7_HS746T_FM_240min = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_ID7_HS746T_FM_240min * ( 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pEGFR_ID7_HS746T_FM_3min" 
		observableParameter1_pEGFR_ID7_HS746T_FM_3min = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_ID7_HS746T_FM_3min * ( 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pEGFR_ID8_HS746T_HM_15min" 
		observableParameter1_pEGFR_ID8_HS746T_HM_15min = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_ID8_HS746T_HM_15min * ( 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pEGFR_ID8_HS746T_HM_1min" 
		observableParameter1_pEGFR_ID8_HS746T_HM_1min = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_ID8_HS746T_HM_1min * ( 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pEGFR_ID8_HS746T_HM_240min" 
		observableParameter1_pEGFR_ID8_HS746T_HM_240min = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_ID8_HS746T_HM_240min * ( 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pEGFR_ID8_HS746T_HM_3min" 
		observableParameter1_pEGFR_ID8_HS746T_HM_3min = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_ID8_HS746T_HM_3min * ( 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pEGFR_ID9_HS746T_FM" 
		observableParameter1_pEGFR_ID9_HS746T_FM = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_ID9_HS746T_FM * ( 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pEGFR_ID9_HS746T_HM_30EGF" 
		observableParameter1_pEGFR_ID9_HS746T_HM_30EGF = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_ID9_HS746T_HM_30EGF * ( 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pEGFR_ID9_HS746T_HM_5EGF" 
		observableParameter1_pEGFR_ID9_HS746T_HM_5EGF = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_ID9_HS746T_HM_5EGF * ( 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pEGFR_ID9b_MKN1_HM_1EGF" 
		observableParameter1_pEGFR_ID9b_MKN1_HM_1EGF = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_ID9b_MKN1_HM_1EGF * ( 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pEGFR_ID9b_MKN1_HM_5EGF" 
		observableParameter1_pEGFR_ID9b_MKN1_HM_5EGF = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_ID9b_MKN1_HM_5EGF * ( 2 * pEGFR_EGF_2 + 2 * pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i ) 
	end

	if observableId == "pMAPK_ID10_HS746T_FM_EGF_CET" 
		observableParameter1_pMAPK_ID10_HS746T_FM_EGF_CET = getObsOrSdParam(obsPar, mapObsParam)
		return pMAPK * observableParameter1_pMAPK_ID10_HS746T_FM_EGF_CET 
	end

	if observableId == "pMAPK_ID10_HS746T_FM_EGF" 
		observableParameter1_pMAPK_ID10_HS746T_FM_EGF = getObsOrSdParam(obsPar, mapObsParam)
		return pMAPK * observableParameter1_pMAPK_ID10_HS746T_FM_EGF 
	end

	if observableId == "pMAPK_ID12_MKN1_FM_HM" 
		observableParameter1_pMAPK_ID12_MKN1_FM_HM = getObsOrSdParam(obsPar, mapObsParam)
		return pMAPK * observableParameter1_pMAPK_ID12_MKN1_FM_HM 
	end

	if observableId == "pMAPK_ID2_HS746T_HM" 
		observableParameter1_pMAPK_ID2_HS746T_HM = getObsOrSdParam(obsPar, mapObsParam)
		return pMAPK * observableParameter1_pMAPK_ID2_HS746T_HM 
	end

	if observableId == "pMAPK_ID2_MKN1_HM" 
		observableParameter1_pMAPK_ID2_MKN1_HM = getObsOrSdParam(obsPar, mapObsParam)
		return pMAPK * observableParameter1_pMAPK_ID2_MKN1_HM 
	end

	if observableId == "pMAPK_ID5_MKN1_FM_15min" 
		observableParameter1_pMAPK_ID5_MKN1_FM_15min = getObsOrSdParam(obsPar, mapObsParam)
		return pMAPK * observableParameter1_pMAPK_ID5_MKN1_FM_15min 
	end

	if observableId == "pMAPK_ID5_MKN1_FM_1min" 
		observableParameter1_pMAPK_ID5_MKN1_FM_1min = getObsOrSdParam(obsPar, mapObsParam)
		return pMAPK * observableParameter1_pMAPK_ID5_MKN1_FM_1min 
	end

	if observableId == "pMAPK_ID5_MKN1_FM_240min" 
		observableParameter1_pMAPK_ID5_MKN1_FM_240min = getObsOrSdParam(obsPar, mapObsParam)
		return pMAPK * observableParameter1_pMAPK_ID5_MKN1_FM_240min 
	end

	if observableId == "pMAPK_ID5_MKN1_FM_3min" 
		observableParameter1_pMAPK_ID5_MKN1_FM_3min = getObsOrSdParam(obsPar, mapObsParam)
		return pMAPK * observableParameter1_pMAPK_ID5_MKN1_FM_3min 
	end

	if observableId == "pMAPK_ID7_HS746T_FM_15min" 
		observableParameter1_pMAPK_ID7_HS746T_FM_15min = getObsOrSdParam(obsPar, mapObsParam)
		return pMAPK * observableParameter1_pMAPK_ID7_HS746T_FM_15min 
	end

	if observableId == "pMAPK_ID7_HS746T_FM_1min" 
		observableParameter1_pMAPK_ID7_HS746T_FM_1min = getObsOrSdParam(obsPar, mapObsParam)
		return pMAPK * observableParameter1_pMAPK_ID7_HS746T_FM_1min 
	end

	if observableId == "pMAPK_ID7_HS746T_FM_240min" 
		observableParameter1_pMAPK_ID7_HS746T_FM_240min = getObsOrSdParam(obsPar, mapObsParam)
		return pMAPK * observableParameter1_pMAPK_ID7_HS746T_FM_240min 
	end

	if observableId == "pMAPK_ID7_HS746T_FM_3min" 
		observableParameter1_pMAPK_ID7_HS746T_FM_3min = getObsOrSdParam(obsPar, mapObsParam)
		return pMAPK * observableParameter1_pMAPK_ID7_HS746T_FM_3min 
	end

	if observableId == "pMAPK_ID8_HS746T_HM_15min" 
		observableParameter1_pMAPK_ID8_HS746T_HM_15min = getObsOrSdParam(obsPar, mapObsParam)
		return pMAPK * observableParameter1_pMAPK_ID8_HS746T_HM_15min 
	end

	if observableId == "pMAPK_ID8_HS746T_HM_1min" 
		observableParameter1_pMAPK_ID8_HS746T_HM_1min = getObsOrSdParam(obsPar, mapObsParam)
		return pMAPK * observableParameter1_pMAPK_ID8_HS746T_HM_1min 
	end

	if observableId == "pMAPK_ID8_HS746T_HM_240min" 
		observableParameter1_pMAPK_ID8_HS746T_HM_240min = getObsOrSdParam(obsPar, mapObsParam)
		return pMAPK * observableParameter1_pMAPK_ID8_HS746T_HM_240min 
	end

	if observableId == "pMAPK_ID8_HS746T_HM_3min" 
		observableParameter1_pMAPK_ID8_HS746T_HM_3min = getObsOrSdParam(obsPar, mapObsParam)
		return pMAPK * observableParameter1_pMAPK_ID8_HS746T_HM_3min 
	end

	if observableId == "pMAPK_ID9_HS746T_FM" 
		observableParameter1_pMAPK_ID9_HS746T_FM = getObsOrSdParam(obsPar, mapObsParam)
		return pMAPK * observableParameter1_pMAPK_ID9_HS746T_FM 
	end

	if observableId == "pMAPK_ID9_HS746T_HM_30EGF" 
		observableParameter1_pMAPK_ID9_HS746T_HM_30EGF = getObsOrSdParam(obsPar, mapObsParam)
		return pMAPK * observableParameter1_pMAPK_ID9_HS746T_HM_30EGF 
	end

	if observableId == "pMAPK_ID9_HS746T_HM_5EGF" 
		observableParameter1_pMAPK_ID9_HS746T_HM_5EGF = getObsOrSdParam(obsPar, mapObsParam)
		return pMAPK * observableParameter1_pMAPK_ID9_HS746T_HM_5EGF 
	end

	if observableId == "pMAPK_ID9b_MKN1_HM_1EGF" 
		observableParameter1_pMAPK_ID9b_MKN1_HM_1EGF = getObsOrSdParam(obsPar, mapObsParam)
		return pMAPK * observableParameter1_pMAPK_ID9b_MKN1_HM_1EGF 
	end

	if observableId == "pMAPK_ID9b_MKN1_HM_5EGF" 
		observableParameter1_pMAPK_ID9b_MKN1_HM_5EGF = getObsOrSdParam(obsPar, mapObsParam)
		return pMAPK * observableParameter1_pMAPK_ID9b_MKN1_HM_5EGF 
	end

end

function evalU0!(u0Vec, paramVec) 

	d_ksyn_EGFR__fm_2_hm, xi_kexp_pMMET_2_i, d_ka_RAS__pEGFR_EGF_2__MKN1_2_HS746T, CET_bool1, EGF_bool1, d_ka_MAPK__MKN1_2_HS746T, ki_RAS__MKN1, d_ka_AKT__MKN1_2_HS746T, KD_EGFR_CET, xi_ka_PI3K_pMMET_2, d_MAPKtotal__fm_2_hm, d_ksyn_MMET__fm_2_hm, ka_PI3K__basal, d_kimp_pEGFR_EGF_2__MKN1_2_HS746T, d_AKTtotal__fm_2_hm, KD_METinh, xi_ki_MPI3K, d_kexp_pEGFR_EGF_2_i__MKN1_2_HS746T, d_ksyn_EGFR__MKN1_2_HS746T, d_MPI3Ktotal__fm_2_hm, d_kdeg_membran__MKN1_2_HS746T, MPI3Ktotal__MKN1_fm, KD_EGFR_EGF, METinh_bool1, kexp_pEGFR_EGF_2_i__MKN1, CET_level, kbin_EGFR_EGF, ksyn_EGFR__MKN1_fm, xi_kdeg_pMMET_pEGFR_i, d_ki_PI3K__MKN1_2_HS746T, xi_kpho_MMET_EGFR, MAPKtotal__MKN1_fm, xi_kexp_pMMET_pEGFR_i, d_RAStotal__MKN1_2_HS746T, kpho_EGFR_EGF, ka_MAPK__MKN1, kpho_MMET_EGFR, full_medium, xi_ka_RAS_pMMET_2, xi_kimp_pMMET_pEGFR, xi_kdeg_pMMET_2_i, ka_RAS__basal__MKN1, d_RAStotal__fm_2_hm, kdim_MMET, MKN1, ksyn_MMET__HS746T_fm, kdim_MMET_EGFR, ka_PI3K__pEGFR_EGF_2, kdeg_membran__MKN1, HS746T, d_ki_RAS__MKN1_2_HS746T, EGF_level, kdim_EGFR_EGF, kbin_EGFR_CET, ka_AKT__MKN1, d_ki_AKT__MKN1_2_HS746T, PI3Ktotal__HS746T_fm, ki_PI3K__MKN1, kpho_MMET, AKTtotal__MKN1_fm, xi_kimp_pMMET_2, relative_ksyn_EGFR, ka_RAS__pEGFR_EGF_2__MKN1, ki_MAPK, cyt, kimp_pEGFR_EGF_2__MKN1, d_AKTtotal__MKN1_2_HS746T, ki_AKT__MKN1, RAStotal__MKN1_fm, d_PI3Ktotal__fm_2_hm, xi_ka_PI3K_pMMET_pEGFR, xi_kdim_MMET_EGFR, xi_kpho_MMET, METinh_level, kdim_MMETinh, d_MAPKtotal__MKN1_2_HS746T, xi_ka_RAS_pMMET_pEGFR, kdeg_pEGFR_EGF_2_i__MKN1, xi_kdim_MMET, d_kdeg_pEGFR_EGF_2_i__MKN1_2_HS746T, d_ka_RAS__basal__MKN1_2_HS746T, hunger_medium = paramVec 

	EGFR = 0.0 
	pPI3K = 0.0 
	MMET_METinh_2 = 0.0 
	pMAPK = 0.0 
	MMET_2 = 0.0 
	pEGFR_EGF_2 = 0.0 
	pEGFR_EGF_2_i = 0.0 
	pMPI3K = 0.0 
	EGFR_EGF_2 = 0.0 
	MMET = 0.0 
	pMMET_2 = 0.0 
	RAS_GTP = 0.0 
	pMMET_2_i = 0.0 
	pMMET_pEGFR_i = 0.0 
	MMET_MMET_METinh = 0.0 
	pAKT = 0.0 
	EGFR_EGF = 0.0 
	MMET_METinh = 0.0 
	MMET_EGFR = 0.0 
	EGFR_CET = 0.0 
	pMMET_pEGFR = 0.0 
	EGFR_MMET_METinh = 0.0 

	u0Vec .= EGFR, pPI3K, MMET_METinh_2, pMAPK, MMET_2, pEGFR_EGF_2, pEGFR_EGF_2_i, pMPI3K, EGFR_EGF_2, MMET, pMMET_2, RAS_GTP, pMMET_2_i, pMMET_pEGFR_i, MMET_MMET_METinh, pAKT, EGFR_EGF, MMET_METinh, MMET_EGFR, EGFR_CET, pMMET_pEGFR, EGFR_MMET_METinh
end

function evalU0(paramVec) 

	d_ksyn_EGFR__fm_2_hm, xi_kexp_pMMET_2_i, d_ka_RAS__pEGFR_EGF_2__MKN1_2_HS746T, CET_bool1, EGF_bool1, d_ka_MAPK__MKN1_2_HS746T, ki_RAS__MKN1, d_ka_AKT__MKN1_2_HS746T, KD_EGFR_CET, xi_ka_PI3K_pMMET_2, d_MAPKtotal__fm_2_hm, d_ksyn_MMET__fm_2_hm, ka_PI3K__basal, d_kimp_pEGFR_EGF_2__MKN1_2_HS746T, d_AKTtotal__fm_2_hm, KD_METinh, xi_ki_MPI3K, d_kexp_pEGFR_EGF_2_i__MKN1_2_HS746T, d_ksyn_EGFR__MKN1_2_HS746T, d_MPI3Ktotal__fm_2_hm, d_kdeg_membran__MKN1_2_HS746T, MPI3Ktotal__MKN1_fm, KD_EGFR_EGF, METinh_bool1, kexp_pEGFR_EGF_2_i__MKN1, CET_level, kbin_EGFR_EGF, ksyn_EGFR__MKN1_fm, xi_kdeg_pMMET_pEGFR_i, d_ki_PI3K__MKN1_2_HS746T, xi_kpho_MMET_EGFR, MAPKtotal__MKN1_fm, xi_kexp_pMMET_pEGFR_i, d_RAStotal__MKN1_2_HS746T, kpho_EGFR_EGF, ka_MAPK__MKN1, kpho_MMET_EGFR, full_medium, xi_ka_RAS_pMMET_2, xi_kimp_pMMET_pEGFR, xi_kdeg_pMMET_2_i, ka_RAS__basal__MKN1, d_RAStotal__fm_2_hm, kdim_MMET, MKN1, ksyn_MMET__HS746T_fm, kdim_MMET_EGFR, ka_PI3K__pEGFR_EGF_2, kdeg_membran__MKN1, HS746T, d_ki_RAS__MKN1_2_HS746T, EGF_level, kdim_EGFR_EGF, kbin_EGFR_CET, ka_AKT__MKN1, d_ki_AKT__MKN1_2_HS746T, PI3Ktotal__HS746T_fm, ki_PI3K__MKN1, kpho_MMET, AKTtotal__MKN1_fm, xi_kimp_pMMET_2, relative_ksyn_EGFR, ka_RAS__pEGFR_EGF_2__MKN1, ki_MAPK, cyt, kimp_pEGFR_EGF_2__MKN1, d_AKTtotal__MKN1_2_HS746T, ki_AKT__MKN1, RAStotal__MKN1_fm, d_PI3Ktotal__fm_2_hm, xi_ka_PI3K_pMMET_pEGFR, xi_kdim_MMET_EGFR, xi_kpho_MMET, METinh_level, kdim_MMETinh, d_MAPKtotal__MKN1_2_HS746T, xi_ka_RAS_pMMET_pEGFR, kdeg_pEGFR_EGF_2_i__MKN1, xi_kdim_MMET, d_kdeg_pEGFR_EGF_2_i__MKN1_2_HS746T, d_ka_RAS__basal__MKN1_2_HS746T, hunger_medium = paramVec 

	EGFR = 0.0 
	pPI3K = 0.0 
	MMET_METinh_2 = 0.0 
	pMAPK = 0.0 
	MMET_2 = 0.0 
	pEGFR_EGF_2 = 0.0 
	pEGFR_EGF_2_i = 0.0 
	pMPI3K = 0.0 
	EGFR_EGF_2 = 0.0 
	MMET = 0.0 
	pMMET_2 = 0.0 
	RAS_GTP = 0.0 
	pMMET_2_i = 0.0 
	pMMET_pEGFR_i = 0.0 
	MMET_MMET_METinh = 0.0 
	pAKT = 0.0 
	EGFR_EGF = 0.0 
	MMET_METinh = 0.0 
	MMET_EGFR = 0.0 
	EGFR_CET = 0.0 
	pMMET_pEGFR = 0.0 
	EGFR_MMET_METinh = 0.0 

	 return [EGFR, pPI3K, MMET_METinh_2, pMAPK, MMET_2, pEGFR_EGF_2, pEGFR_EGF_2_i, pMPI3K, EGFR_EGF_2, MMET, pMMET_2, RAS_GTP, pMMET_2_i, pMMET_pEGFR_i, MMET_MMET_METinh, pAKT, EGFR_EGF, MMET_METinh, MMET_EGFR, EGFR_CET, pMMET_pEGFR, EGFR_MMET_METinh]
end

function evalSd!(u, t, sdPar, dynPar, nonDynParam, paramData, observableId, mapSdParam) 

	EGFR, pPI3K, MMET_METinh_2, pMAPK, MMET_2, pEGFR_EGF_2, pEGFR_EGF_2_i, pMPI3K, EGFR_EGF_2, MMET, pMMET_2, RAS_GTP, pMMET_2_i, pMMET_pEGFR_i, MMET_MMET_METinh, pAKT, EGFR_EGF, MMET_METinh, MMET_EGFR, EGFR_CET, pMMET_pEGFR, EGFR_MMET_METinh= u 
	KD_EGFR_CET, KD_EGFR_EGF, KD_METinh, d_AKTtotal__MKN1_2_HS746T, d_AKTtotal__fm_2_hm, d_MAPKtotal__MKN1_2_HS746T, d_MAPKtotal__fm_2_hm, d_MPI3Ktotal__fm_2_hm, d_PI3Ktotal__fm_2_hm, d_RAStotal__MKN1_2_HS746T, d_RAStotal__fm_2_hm, d_kdeg_membran__MKN1_2_HS746T, d_kdeg_pEGFR_EGF_2_i__MKN1_2_HS746T, d_kexp_pEGFR_EGF_2_i__MKN1_2_HS746T, d_kimp_pEGFR_EGF_2__MKN1_2_HS746T, d_ksyn_EGFR__MKN1_2_HS746T, d_ksyn_EGFR__fm_2_hm, d_ksyn_MMET__fm_2_hm, ka_AKT__MKN1, ka_MAPK__MKN1, ka_PI3K__basal, ka_PI3K__pEGFR_EGF_2, ka_RAS__basal__MKN1, ka_RAS__pEGFR_EGF_2__MKN1, kbin_EGFR_CET, kbin_EGFR_EGF, kdeg_membran__MKN1, kdeg_pEGFR_EGF_2_i__MKN1, kdim_EGFR_EGF, kdim_MMET, kdim_MMET_EGFR, kdim_MMETinh, kexp_pEGFR_EGF_2_i__MKN1, ki_AKT__MKN1, ki_MAPK, ki_PI3K__MKN1, ki_RAS__MKN1, kimp_pEGFR_EGF_2__MKN1, kpho_EGFR_EGF, kpho_MMET, kpho_MMET_EGFR, ksyn_MMET__HS746T_fm, xi_ka_PI3K_pMMET_2, xi_ka_PI3K_pMMET_pEGFR, xi_ka_RAS_pMMET_2, xi_ka_RAS_pMMET_pEGFR, xi_kdeg_pMMET_2_i, xi_kdeg_pMMET_pEGFR_i, xi_kdim_MMET, xi_kdim_MMET_EGFR, xi_kexp_pMMET_2_i, xi_kexp_pMMET_pEGFR_i, xi_ki_MPI3K, xi_kimp_pMMET_2, xi_kimp_pMMET_pEGFR, xi_kpho_MMET, xi_kpho_MMET_EGFR = dynPar 
	AKTtotal__MKN1_fm_C = paramData.paramVal[1] 
	MAPKtotal__MKN1_fm_C = paramData.paramVal[5] 
	MPI3Ktotal__MKN1_fm_C = paramData.paramVal[6] 
	PI3Ktotal__HS746T_fm_C = paramData.paramVal[7] 
	RAStotal__MKN1_fm_C = paramData.paramVal[8] 
	d_ka_AKT__MKN1_2_HS746T_C = paramData.paramVal[17] 
	d_ka_MAPK__MKN1_2_HS746T_C = paramData.paramVal[18] 
	d_ka_RAS__basal__MKN1_2_HS746T_C = paramData.paramVal[19] 
	d_ka_RAS__pEGFR_EGF_2__MKN1_2_HS746T_C = paramData.paramVal[20] 
	d_ki_AKT__MKN1_2_HS746T_C = paramData.paramVal[24] 
	d_ki_PI3K__MKN1_2_HS746T_C = paramData.paramVal[25] 
	d_ki_RAS__MKN1_2_HS746T_C = paramData.paramVal[26] 
	ksyn_EGFR__MKN1_fm_C = paramData.paramVal[54] 

	if observableId == "EGFR_ID10_HS746T_FM_EGF_CET" 
		noiseParameter1_EGFR_ID10_HS746T_FM_EGF_CET = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_EGFR_ID10_HS746T_FM_EGF_CET 
	end

	if observableId == "EGFR_ID10_HS746T_FM_EGF" 
		noiseParameter1_EGFR_ID10_HS746T_FM_EGF = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_EGFR_ID10_HS746T_FM_EGF 
	end

	if observableId == "EGFR_ID12_MKN1_FM_HM" 
		noiseParameter1_EGFR_ID12_MKN1_FM_HM = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_EGFR_ID12_MKN1_FM_HM 
	end

	if observableId == "EGFR_ID3_HS746T_FM" 
		noiseParameter1_EGFR_ID3_HS746T_FM = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_EGFR_ID3_HS746T_FM 
	end

	if observableId == "EGFR_ID3_MKN1_FM" 
		noiseParameter1_EGFR_ID3_MKN1_FM = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_EGFR_ID3_MKN1_FM 
	end

	if observableId == "EGFR_ID4_HS746T_HM_3min" 
		noiseParameter1_EGFR_ID4_HS746T_HM_3min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_EGFR_ID4_HS746T_HM_3min 
	end

	if observableId == "EGFR_ID4_MKN1_HM_3min" 
		noiseParameter1_EGFR_ID4_MKN1_HM_3min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_EGFR_ID4_MKN1_HM_3min 
	end

	if observableId == "MMET_ID10_HS746T_FM_EGF_CET" 
		noiseParameter1_MMET_ID10_HS746T_FM_EGF_CET = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_MMET_ID10_HS746T_FM_EGF_CET 
	end

	if observableId == "MMET_ID10_HS746T_FM_EGF" 
		noiseParameter1_MMET_ID10_HS746T_FM_EGF = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_MMET_ID10_HS746T_FM_EGF 
	end

	if observableId == "pAKT_ID12_MKN1_FM_HM" 
		noiseParameter1_pAKT_ID12_MKN1_FM_HM = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pAKT_ID12_MKN1_FM_HM 
	end

	if observableId == "pAKT_ID2_HS746T_HM" 
		noiseParameter1_pAKT_ID2_HS746T_HM = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pAKT_ID2_HS746T_HM 
	end

	if observableId == "pAKT_ID2_MKN1_HM" 
		noiseParameter1_pAKT_ID2_MKN1_HM = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pAKT_ID2_MKN1_HM 
	end

	if observableId == "pAKT_ID5_MKN1_FM_15min" 
		noiseParameter1_pAKT_ID5_MKN1_FM_15min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pAKT_ID5_MKN1_FM_15min 
	end

	if observableId == "pAKT_ID5_MKN1_FM_1min" 
		noiseParameter1_pAKT_ID5_MKN1_FM_1min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pAKT_ID5_MKN1_FM_1min 
	end

	if observableId == "pAKT_ID5_MKN1_FM_240min" 
		noiseParameter1_pAKT_ID5_MKN1_FM_240min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pAKT_ID5_MKN1_FM_240min 
	end

	if observableId == "pAKT_ID5_MKN1_FM_3min" 
		noiseParameter1_pAKT_ID5_MKN1_FM_3min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pAKT_ID5_MKN1_FM_3min 
	end

	if observableId == "pAKT_ID7_HS746T_FM_15min" 
		noiseParameter1_pAKT_ID7_HS746T_FM_15min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pAKT_ID7_HS746T_FM_15min 
	end

	if observableId == "pAKT_ID7_HS746T_FM_1min" 
		noiseParameter1_pAKT_ID7_HS746T_FM_1min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pAKT_ID7_HS746T_FM_1min 
	end

	if observableId == "pAKT_ID7_HS746T_FM_240min" 
		noiseParameter1_pAKT_ID7_HS746T_FM_240min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pAKT_ID7_HS746T_FM_240min 
	end

	if observableId == "pAKT_ID7_HS746T_FM_3min" 
		noiseParameter1_pAKT_ID7_HS746T_FM_3min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pAKT_ID7_HS746T_FM_3min 
	end

	if observableId == "pAKT_ID8_HS746T_HM_15min" 
		noiseParameter1_pAKT_ID8_HS746T_HM_15min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pAKT_ID8_HS746T_HM_15min 
	end

	if observableId == "pAKT_ID8_HS746T_HM_1min" 
		noiseParameter1_pAKT_ID8_HS746T_HM_1min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pAKT_ID8_HS746T_HM_1min 
	end

	if observableId == "pAKT_ID8_HS746T_HM_240min" 
		noiseParameter1_pAKT_ID8_HS746T_HM_240min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pAKT_ID8_HS746T_HM_240min 
	end

	if observableId == "pAKT_ID8_HS746T_HM_3min" 
		noiseParameter1_pAKT_ID8_HS746T_HM_3min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pAKT_ID8_HS746T_HM_3min 
	end

	if observableId == "pAKT_ID9_HS746T_FM" 
		noiseParameter1_pAKT_ID9_HS746T_FM = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pAKT_ID9_HS746T_FM 
	end

	if observableId == "pAKT_ID9_HS746T_HM_30EGF" 
		noiseParameter1_pAKT_ID9_HS746T_HM_30EGF = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pAKT_ID9_HS746T_HM_30EGF 
	end

	if observableId == "pAKT_ID9_HS746T_HM_5EGF" 
		noiseParameter1_pAKT_ID9_HS746T_HM_5EGF = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pAKT_ID9_HS746T_HM_5EGF 
	end

	if observableId == "pEGFR_ID10_HS746T_FM_EGF_CET" 
		noiseParameter1_pEGFR_ID10_HS746T_FM_EGF_CET = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_ID10_HS746T_FM_EGF_CET 
	end

	if observableId == "pEGFR_ID10_HS746T_FM_EGF" 
		noiseParameter1_pEGFR_ID10_HS746T_FM_EGF = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_ID10_HS746T_FM_EGF 
	end

	if observableId == "pEGFR_ID12_MKN1_FM_HM" 
		noiseParameter1_pEGFR_ID12_MKN1_FM_HM = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_ID12_MKN1_FM_HM 
	end

	if observableId == "pEGFR_ID2_HS746T_HM" 
		noiseParameter1_pEGFR_ID2_HS746T_HM = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_ID2_HS746T_HM 
	end

	if observableId == "pEGFR_ID2_MKN1_HM" 
		noiseParameter1_pEGFR_ID2_MKN1_HM = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_ID2_MKN1_HM 
	end

	if observableId == "pEGFR_ID3_HS746T_FM" 
		noiseParameter1_pEGFR_ID3_HS746T_FM = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_ID3_HS746T_FM 
	end

	if observableId == "pEGFR_ID3_MKN1_FM" 
		noiseParameter1_pEGFR_ID3_MKN1_FM = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_ID3_MKN1_FM 
	end

	if observableId == "pEGFR_ID4_HS746T_HM_3min" 
		noiseParameter1_pEGFR_ID4_HS746T_HM_3min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_ID4_HS746T_HM_3min 
	end

	if observableId == "pEGFR_ID4_MKN1_HM_3min" 
		noiseParameter1_pEGFR_ID4_MKN1_HM_3min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_ID4_MKN1_HM_3min 
	end

	if observableId == "pEGFR_ID5_MKN1_FM_15min" 
		noiseParameter1_pEGFR_ID5_MKN1_FM_15min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_ID5_MKN1_FM_15min 
	end

	if observableId == "pEGFR_ID5_MKN1_FM_1min" 
		noiseParameter1_pEGFR_ID5_MKN1_FM_1min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_ID5_MKN1_FM_1min 
	end

	if observableId == "pEGFR_ID5_MKN1_FM_240min" 
		noiseParameter1_pEGFR_ID5_MKN1_FM_240min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_ID5_MKN1_FM_240min 
	end

	if observableId == "pEGFR_ID5_MKN1_FM_3min" 
		noiseParameter1_pEGFR_ID5_MKN1_FM_3min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_ID5_MKN1_FM_3min 
	end

	if observableId == "pEGFR_ID6_MKN1_HM_15min" 
		noiseParameter1_pEGFR_ID6_MKN1_HM_15min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_ID6_MKN1_HM_15min 
	end

	if observableId == "pEGFR_ID6_MKN1_HM_1min" 
		noiseParameter1_pEGFR_ID6_MKN1_HM_1min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_ID6_MKN1_HM_1min 
	end

	if observableId == "pEGFR_ID6_MKN1_HM_240min" 
		noiseParameter1_pEGFR_ID6_MKN1_HM_240min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_ID6_MKN1_HM_240min 
	end

	if observableId == "pEGFR_ID6_MKN1_HM_3min" 
		noiseParameter1_pEGFR_ID6_MKN1_HM_3min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_ID6_MKN1_HM_3min 
	end

	if observableId == "pEGFR_ID7_HS746T_FM_15min" 
		noiseParameter1_pEGFR_ID7_HS746T_FM_15min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_ID7_HS746T_FM_15min 
	end

	if observableId == "pEGFR_ID7_HS746T_FM_1min" 
		noiseParameter1_pEGFR_ID7_HS746T_FM_1min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_ID7_HS746T_FM_1min 
	end

	if observableId == "pEGFR_ID7_HS746T_FM_240min" 
		noiseParameter1_pEGFR_ID7_HS746T_FM_240min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_ID7_HS746T_FM_240min 
	end

	if observableId == "pEGFR_ID7_HS746T_FM_3min" 
		noiseParameter1_pEGFR_ID7_HS746T_FM_3min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_ID7_HS746T_FM_3min 
	end

	if observableId == "pEGFR_ID8_HS746T_HM_15min" 
		noiseParameter1_pEGFR_ID8_HS746T_HM_15min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_ID8_HS746T_HM_15min 
	end

	if observableId == "pEGFR_ID8_HS746T_HM_1min" 
		noiseParameter1_pEGFR_ID8_HS746T_HM_1min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_ID8_HS746T_HM_1min 
	end

	if observableId == "pEGFR_ID8_HS746T_HM_240min" 
		noiseParameter1_pEGFR_ID8_HS746T_HM_240min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_ID8_HS746T_HM_240min 
	end

	if observableId == "pEGFR_ID8_HS746T_HM_3min" 
		noiseParameter1_pEGFR_ID8_HS746T_HM_3min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_ID8_HS746T_HM_3min 
	end

	if observableId == "pEGFR_ID9_HS746T_FM" 
		noiseParameter1_pEGFR_ID9_HS746T_FM = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_ID9_HS746T_FM 
	end

	if observableId == "pEGFR_ID9_HS746T_HM_30EGF" 
		noiseParameter1_pEGFR_ID9_HS746T_HM_30EGF = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_ID9_HS746T_HM_30EGF 
	end

	if observableId == "pEGFR_ID9_HS746T_HM_5EGF" 
		noiseParameter1_pEGFR_ID9_HS746T_HM_5EGF = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_ID9_HS746T_HM_5EGF 
	end

	if observableId == "pEGFR_ID9b_MKN1_HM_1EGF" 
		noiseParameter1_pEGFR_ID9b_MKN1_HM_1EGF = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_ID9b_MKN1_HM_1EGF 
	end

	if observableId == "pEGFR_ID9b_MKN1_HM_5EGF" 
		noiseParameter1_pEGFR_ID9b_MKN1_HM_5EGF = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_ID9b_MKN1_HM_5EGF 
	end

	if observableId == "pMAPK_ID10_HS746T_FM_EGF_CET" 
		noiseParameter1_pMAPK_ID10_HS746T_FM_EGF_CET = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pMAPK_ID10_HS746T_FM_EGF_CET 
	end

	if observableId == "pMAPK_ID10_HS746T_FM_EGF" 
		noiseParameter1_pMAPK_ID10_HS746T_FM_EGF = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pMAPK_ID10_HS746T_FM_EGF 
	end

	if observableId == "pMAPK_ID12_MKN1_FM_HM" 
		noiseParameter1_pMAPK_ID12_MKN1_FM_HM = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pMAPK_ID12_MKN1_FM_HM 
	end

	if observableId == "pMAPK_ID2_HS746T_HM" 
		noiseParameter1_pMAPK_ID2_HS746T_HM = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pMAPK_ID2_HS746T_HM 
	end

	if observableId == "pMAPK_ID2_MKN1_HM" 
		noiseParameter1_pMAPK_ID2_MKN1_HM = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pMAPK_ID2_MKN1_HM 
	end

	if observableId == "pMAPK_ID5_MKN1_FM_15min" 
		noiseParameter1_pMAPK_ID5_MKN1_FM_15min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pMAPK_ID5_MKN1_FM_15min 
	end

	if observableId == "pMAPK_ID5_MKN1_FM_1min" 
		noiseParameter1_pMAPK_ID5_MKN1_FM_1min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pMAPK_ID5_MKN1_FM_1min 
	end

	if observableId == "pMAPK_ID5_MKN1_FM_240min" 
		noiseParameter1_pMAPK_ID5_MKN1_FM_240min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pMAPK_ID5_MKN1_FM_240min 
	end

	if observableId == "pMAPK_ID5_MKN1_FM_3min" 
		noiseParameter1_pMAPK_ID5_MKN1_FM_3min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pMAPK_ID5_MKN1_FM_3min 
	end

	if observableId == "pMAPK_ID7_HS746T_FM_15min" 
		noiseParameter1_pMAPK_ID7_HS746T_FM_15min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pMAPK_ID7_HS746T_FM_15min 
	end

	if observableId == "pMAPK_ID7_HS746T_FM_1min" 
		noiseParameter1_pMAPK_ID7_HS746T_FM_1min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pMAPK_ID7_HS746T_FM_1min 
	end

	if observableId == "pMAPK_ID7_HS746T_FM_240min" 
		noiseParameter1_pMAPK_ID7_HS746T_FM_240min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pMAPK_ID7_HS746T_FM_240min 
	end

	if observableId == "pMAPK_ID7_HS746T_FM_3min" 
		noiseParameter1_pMAPK_ID7_HS746T_FM_3min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pMAPK_ID7_HS746T_FM_3min 
	end

	if observableId == "pMAPK_ID8_HS746T_HM_15min" 
		noiseParameter1_pMAPK_ID8_HS746T_HM_15min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pMAPK_ID8_HS746T_HM_15min 
	end

	if observableId == "pMAPK_ID8_HS746T_HM_1min" 
		noiseParameter1_pMAPK_ID8_HS746T_HM_1min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pMAPK_ID8_HS746T_HM_1min 
	end

	if observableId == "pMAPK_ID8_HS746T_HM_240min" 
		noiseParameter1_pMAPK_ID8_HS746T_HM_240min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pMAPK_ID8_HS746T_HM_240min 
	end

	if observableId == "pMAPK_ID8_HS746T_HM_3min" 
		noiseParameter1_pMAPK_ID8_HS746T_HM_3min = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pMAPK_ID8_HS746T_HM_3min 
	end

	if observableId == "pMAPK_ID9_HS746T_FM" 
		noiseParameter1_pMAPK_ID9_HS746T_FM = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pMAPK_ID9_HS746T_FM 
	end

	if observableId == "pMAPK_ID9_HS746T_HM_30EGF" 
		noiseParameter1_pMAPK_ID9_HS746T_HM_30EGF = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pMAPK_ID9_HS746T_HM_30EGF 
	end

	if observableId == "pMAPK_ID9_HS746T_HM_5EGF" 
		noiseParameter1_pMAPK_ID9_HS746T_HM_5EGF = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pMAPK_ID9_HS746T_HM_5EGF 
	end

	if observableId == "pMAPK_ID9b_MKN1_HM_1EGF" 
		noiseParameter1_pMAPK_ID9b_MKN1_HM_1EGF = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pMAPK_ID9b_MKN1_HM_1EGF 
	end

	if observableId == "pMAPK_ID9b_MKN1_HM_5EGF" 
		noiseParameter1_pMAPK_ID9b_MKN1_HM_5EGF = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pMAPK_ID9b_MKN1_HM_5EGF 
	end

end