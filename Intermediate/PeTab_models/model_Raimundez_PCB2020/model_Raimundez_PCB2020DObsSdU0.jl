function evalDYmodDu(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	EGFR, pPI3K, MMET_METinh_2, pMAPK, MMET_2, pEGFR_EGF_2, pEGFR_EGF_2_i, pMPI3K, EGFR_EGF_2, MMET, pMMET_2, RAS_GTP, pMMET_2_i, pMMET_pEGFR_i, MMET_MMET_METinh, pAKT, EGFR_EGF, MMET_METinh, MMET_EGFR, EGFR_CET, pMMET_pEGFR, EGFR_MMET_METinh= u 
	d_ksyn_EGFR__fm_2_hm, xi_kexp_pMMET_2_i, d_ka_RAS__pEGFR_EGF_2__MKN1_2_HS746T, CET_bool1, EGF_bool1, d_ka_MAPK__MKN1_2_HS746T, ki_RAS__MKN1, d_ka_AKT__MKN1_2_HS746T, KD_EGFR_CET, xi_ka_PI3K_pMMET_2, d_MAPKtotal__fm_2_hm, d_ksyn_MMET__fm_2_hm, ka_PI3K__basal, d_kimp_pEGFR_EGF_2__MKN1_2_HS746T, d_AKTtotal__fm_2_hm, KD_METinh, xi_ki_MPI3K, d_kexp_pEGFR_EGF_2_i__MKN1_2_HS746T, d_ksyn_EGFR__MKN1_2_HS746T, d_MPI3Ktotal__fm_2_hm, d_kdeg_membran__MKN1_2_HS746T, MPI3Ktotal__MKN1_fm, KD_EGFR_EGF, METinh_bool1, kexp_pEGFR_EGF_2_i__MKN1, CET_level, kbin_EGFR_EGF, ksyn_EGFR__MKN1_fm, xi_kdeg_pMMET_pEGFR_i, d_ki_PI3K__MKN1_2_HS746T, xi_kpho_MMET_EGFR, MAPKtotal__MKN1_fm, xi_kexp_pMMET_pEGFR_i, d_RAStotal__MKN1_2_HS746T, kpho_EGFR_EGF, ka_MAPK__MKN1, kpho_MMET_EGFR, full_medium, xi_ka_RAS_pMMET_2, xi_kimp_pMMET_pEGFR, xi_kdeg_pMMET_2_i, ka_RAS__basal__MKN1, d_RAStotal__fm_2_hm, kdim_MMET, MKN1, ksyn_MMET__HS746T_fm, kdim_MMET_EGFR, ka_PI3K__pEGFR_EGF_2, kdeg_membran__MKN1, HS746T, d_ki_RAS__MKN1_2_HS746T, EGF_level, kdim_EGFR_EGF, kbin_EGFR_CET, ka_AKT__MKN1, d_ki_AKT__MKN1_2_HS746T, PI3Ktotal__HS746T_fm, ki_PI3K__MKN1, kpho_MMET, AKTtotal__MKN1_fm, xi_kimp_pMMET_2, relative_ksyn_EGFR, ka_RAS__pEGFR_EGF_2__MKN1, ki_MAPK, cyt, kimp_pEGFR_EGF_2__MKN1, d_AKTtotal__MKN1_2_HS746T, ki_AKT__MKN1, RAStotal__MKN1_fm, d_PI3Ktotal__fm_2_hm, xi_ka_PI3K_pMMET_pEGFR, xi_kdim_MMET_EGFR, xi_kpho_MMET, METinh_level, kdim_MMETinh, d_MAPKtotal__MKN1_2_HS746T, xi_ka_RAS_pMMET_pEGFR, kdeg_pEGFR_EGF_2_i__MKN1, xi_kdim_MMET, d_kdeg_pEGFR_EGF_2_i__MKN1_2_HS746T, d_ka_RAS__basal__MKN1_2_HS746T, hunger_medium = p 
	if observableId == "EGFR_ID10_HS746T_FM_EGF_CET" 
		observableParameter1_EGFR_ID10_HS746T_FM_EGF_CET = getObsOrSdParam(obsPar, mapObsParam)
		out[1] = observableParameter1_EGFR_ID10_HS746T_FM_EGF_CET
		out[6] = 2observableParameter1_EGFR_ID10_HS746T_FM_EGF_CET
		out[7] = 2observableParameter1_EGFR_ID10_HS746T_FM_EGF_CET
		out[9] = 2observableParameter1_EGFR_ID10_HS746T_FM_EGF_CET
		out[14] = observableParameter1_EGFR_ID10_HS746T_FM_EGF_CET
		out[17] = observableParameter1_EGFR_ID10_HS746T_FM_EGF_CET
		out[19] = observableParameter1_EGFR_ID10_HS746T_FM_EGF_CET
		out[20] = observableParameter1_EGFR_ID10_HS746T_FM_EGF_CET
		out[21] = observableParameter1_EGFR_ID10_HS746T_FM_EGF_CET
		out[22] = observableParameter1_EGFR_ID10_HS746T_FM_EGF_CET
		return nothing
	end

	if observableId == "EGFR_ID10_HS746T_FM_EGF" 
		observableParameter1_EGFR_ID10_HS746T_FM_EGF = getObsOrSdParam(obsPar, mapObsParam)
		out[1] = observableParameter1_EGFR_ID10_HS746T_FM_EGF
		out[6] = 2observableParameter1_EGFR_ID10_HS746T_FM_EGF
		out[7] = 2observableParameter1_EGFR_ID10_HS746T_FM_EGF
		out[9] = 2observableParameter1_EGFR_ID10_HS746T_FM_EGF
		out[14] = observableParameter1_EGFR_ID10_HS746T_FM_EGF
		out[17] = observableParameter1_EGFR_ID10_HS746T_FM_EGF
		out[19] = observableParameter1_EGFR_ID10_HS746T_FM_EGF
		out[20] = observableParameter1_EGFR_ID10_HS746T_FM_EGF
		out[21] = observableParameter1_EGFR_ID10_HS746T_FM_EGF
		out[22] = observableParameter1_EGFR_ID10_HS746T_FM_EGF
		return nothing
	end

	if observableId == "EGFR_ID12_MKN1_FM_HM" 
		observableParameter1_EGFR_ID12_MKN1_FM_HM = getObsOrSdParam(obsPar, mapObsParam)
		out[1] = observableParameter1_EGFR_ID12_MKN1_FM_HM
		out[6] = 2observableParameter1_EGFR_ID12_MKN1_FM_HM
		out[7] = 2observableParameter1_EGFR_ID12_MKN1_FM_HM
		out[9] = 2observableParameter1_EGFR_ID12_MKN1_FM_HM
		out[14] = observableParameter1_EGFR_ID12_MKN1_FM_HM
		out[17] = observableParameter1_EGFR_ID12_MKN1_FM_HM
		out[19] = observableParameter1_EGFR_ID12_MKN1_FM_HM
		out[20] = observableParameter1_EGFR_ID12_MKN1_FM_HM
		out[21] = observableParameter1_EGFR_ID12_MKN1_FM_HM
		out[22] = observableParameter1_EGFR_ID12_MKN1_FM_HM
		return nothing
	end

	if observableId == "EGFR_ID3_HS746T_FM" 
		observableParameter1_EGFR_ID3_HS746T_FM = getObsOrSdParam(obsPar, mapObsParam)
		out[1] = observableParameter1_EGFR_ID3_HS746T_FM
		out[6] = 2observableParameter1_EGFR_ID3_HS746T_FM
		out[7] = 2observableParameter1_EGFR_ID3_HS746T_FM
		out[9] = 2observableParameter1_EGFR_ID3_HS746T_FM
		out[14] = observableParameter1_EGFR_ID3_HS746T_FM
		out[17] = observableParameter1_EGFR_ID3_HS746T_FM
		out[19] = observableParameter1_EGFR_ID3_HS746T_FM
		out[20] = observableParameter1_EGFR_ID3_HS746T_FM
		out[21] = observableParameter1_EGFR_ID3_HS746T_FM
		out[22] = observableParameter1_EGFR_ID3_HS746T_FM
		return nothing
	end

	if observableId == "EGFR_ID3_MKN1_FM" 
		observableParameter1_EGFR_ID3_MKN1_FM = getObsOrSdParam(obsPar, mapObsParam)
		out[1] = observableParameter1_EGFR_ID3_MKN1_FM
		out[6] = 2observableParameter1_EGFR_ID3_MKN1_FM
		out[7] = 2observableParameter1_EGFR_ID3_MKN1_FM
		out[9] = 2observableParameter1_EGFR_ID3_MKN1_FM
		out[14] = observableParameter1_EGFR_ID3_MKN1_FM
		out[17] = observableParameter1_EGFR_ID3_MKN1_FM
		out[19] = observableParameter1_EGFR_ID3_MKN1_FM
		out[20] = observableParameter1_EGFR_ID3_MKN1_FM
		out[21] = observableParameter1_EGFR_ID3_MKN1_FM
		out[22] = observableParameter1_EGFR_ID3_MKN1_FM
		return nothing
	end

	if observableId == "EGFR_ID4_HS746T_HM_3min" 
		observableParameter1_EGFR_ID4_HS746T_HM_3min = getObsOrSdParam(obsPar, mapObsParam)
		out[1] = observableParameter1_EGFR_ID4_HS746T_HM_3min
		out[6] = 2observableParameter1_EGFR_ID4_HS746T_HM_3min
		out[7] = 2observableParameter1_EGFR_ID4_HS746T_HM_3min
		out[9] = 2observableParameter1_EGFR_ID4_HS746T_HM_3min
		out[14] = observableParameter1_EGFR_ID4_HS746T_HM_3min
		out[17] = observableParameter1_EGFR_ID4_HS746T_HM_3min
		out[19] = observableParameter1_EGFR_ID4_HS746T_HM_3min
		out[20] = observableParameter1_EGFR_ID4_HS746T_HM_3min
		out[21] = observableParameter1_EGFR_ID4_HS746T_HM_3min
		out[22] = observableParameter1_EGFR_ID4_HS746T_HM_3min
		return nothing
	end

	if observableId == "EGFR_ID4_MKN1_HM_3min" 
		observableParameter1_EGFR_ID4_MKN1_HM_3min = getObsOrSdParam(obsPar, mapObsParam)
		out[1] = observableParameter1_EGFR_ID4_MKN1_HM_3min
		out[6] = 2observableParameter1_EGFR_ID4_MKN1_HM_3min
		out[7] = 2observableParameter1_EGFR_ID4_MKN1_HM_3min
		out[9] = 2observableParameter1_EGFR_ID4_MKN1_HM_3min
		out[14] = observableParameter1_EGFR_ID4_MKN1_HM_3min
		out[17] = observableParameter1_EGFR_ID4_MKN1_HM_3min
		out[19] = observableParameter1_EGFR_ID4_MKN1_HM_3min
		out[20] = observableParameter1_EGFR_ID4_MKN1_HM_3min
		out[21] = observableParameter1_EGFR_ID4_MKN1_HM_3min
		out[22] = observableParameter1_EGFR_ID4_MKN1_HM_3min
		return nothing
	end

	if observableId == "MMET_ID10_HS746T_FM_EGF_CET" 
		observableParameter1_MMET_ID10_HS746T_FM_EGF_CET = getObsOrSdParam(obsPar, mapObsParam)
		out[3] = 2observableParameter1_MMET_ID10_HS746T_FM_EGF_CET
		out[5] = 2observableParameter1_MMET_ID10_HS746T_FM_EGF_CET
		out[10] = observableParameter1_MMET_ID10_HS746T_FM_EGF_CET
		out[11] = 2observableParameter1_MMET_ID10_HS746T_FM_EGF_CET
		out[13] = 2observableParameter1_MMET_ID10_HS746T_FM_EGF_CET
		out[14] = observableParameter1_MMET_ID10_HS746T_FM_EGF_CET
		out[15] = 2observableParameter1_MMET_ID10_HS746T_FM_EGF_CET
		out[18] = observableParameter1_MMET_ID10_HS746T_FM_EGF_CET
		out[19] = observableParameter1_MMET_ID10_HS746T_FM_EGF_CET
		out[21] = observableParameter1_MMET_ID10_HS746T_FM_EGF_CET
		out[22] = observableParameter1_MMET_ID10_HS746T_FM_EGF_CET
		return nothing
	end

	if observableId == "MMET_ID10_HS746T_FM_EGF" 
		observableParameter1_MMET_ID10_HS746T_FM_EGF = getObsOrSdParam(obsPar, mapObsParam)
		out[3] = 2observableParameter1_MMET_ID10_HS746T_FM_EGF
		out[5] = 2observableParameter1_MMET_ID10_HS746T_FM_EGF
		out[10] = observableParameter1_MMET_ID10_HS746T_FM_EGF
		out[11] = 2observableParameter1_MMET_ID10_HS746T_FM_EGF
		out[13] = 2observableParameter1_MMET_ID10_HS746T_FM_EGF
		out[14] = observableParameter1_MMET_ID10_HS746T_FM_EGF
		out[15] = 2observableParameter1_MMET_ID10_HS746T_FM_EGF
		out[18] = observableParameter1_MMET_ID10_HS746T_FM_EGF
		out[19] = observableParameter1_MMET_ID10_HS746T_FM_EGF
		out[21] = observableParameter1_MMET_ID10_HS746T_FM_EGF
		out[22] = observableParameter1_MMET_ID10_HS746T_FM_EGF
		return nothing
	end

	if observableId == "pAKT_ID12_MKN1_FM_HM" 
		observableParameter1_pAKT_ID12_MKN1_FM_HM = getObsOrSdParam(obsPar, mapObsParam)
		out[16] = observableParameter1_pAKT_ID12_MKN1_FM_HM
		return nothing
	end

	if observableId == "pAKT_ID2_HS746T_HM" 
		observableParameter1_pAKT_ID2_HS746T_HM = getObsOrSdParam(obsPar, mapObsParam)
		out[16] = observableParameter1_pAKT_ID2_HS746T_HM
		return nothing
	end

	if observableId == "pAKT_ID2_MKN1_HM" 
		observableParameter1_pAKT_ID2_MKN1_HM = getObsOrSdParam(obsPar, mapObsParam)
		out[16] = observableParameter1_pAKT_ID2_MKN1_HM
		return nothing
	end

	if observableId == "pAKT_ID5_MKN1_FM_15min" 
		observableParameter1_pAKT_ID5_MKN1_FM_15min = getObsOrSdParam(obsPar, mapObsParam)
		out[16] = observableParameter1_pAKT_ID5_MKN1_FM_15min
		return nothing
	end

	if observableId == "pAKT_ID5_MKN1_FM_1min" 
		observableParameter1_pAKT_ID5_MKN1_FM_1min = getObsOrSdParam(obsPar, mapObsParam)
		out[16] = observableParameter1_pAKT_ID5_MKN1_FM_1min
		return nothing
	end

	if observableId == "pAKT_ID5_MKN1_FM_240min" 
		observableParameter1_pAKT_ID5_MKN1_FM_240min = getObsOrSdParam(obsPar, mapObsParam)
		out[16] = observableParameter1_pAKT_ID5_MKN1_FM_240min
		return nothing
	end

	if observableId == "pAKT_ID5_MKN1_FM_3min" 
		observableParameter1_pAKT_ID5_MKN1_FM_3min = getObsOrSdParam(obsPar, mapObsParam)
		out[16] = observableParameter1_pAKT_ID5_MKN1_FM_3min
		return nothing
	end

	if observableId == "pAKT_ID7_HS746T_FM_15min" 
		observableParameter1_pAKT_ID7_HS746T_FM_15min = getObsOrSdParam(obsPar, mapObsParam)
		out[16] = observableParameter1_pAKT_ID7_HS746T_FM_15min
		return nothing
	end

	if observableId == "pAKT_ID7_HS746T_FM_1min" 
		observableParameter1_pAKT_ID7_HS746T_FM_1min = getObsOrSdParam(obsPar, mapObsParam)
		out[16] = observableParameter1_pAKT_ID7_HS746T_FM_1min
		return nothing
	end

	if observableId == "pAKT_ID7_HS746T_FM_240min" 
		observableParameter1_pAKT_ID7_HS746T_FM_240min = getObsOrSdParam(obsPar, mapObsParam)
		out[16] = observableParameter1_pAKT_ID7_HS746T_FM_240min
		return nothing
	end

	if observableId == "pAKT_ID7_HS746T_FM_3min" 
		observableParameter1_pAKT_ID7_HS746T_FM_3min = getObsOrSdParam(obsPar, mapObsParam)
		out[16] = observableParameter1_pAKT_ID7_HS746T_FM_3min
		return nothing
	end

	if observableId == "pAKT_ID8_HS746T_HM_15min" 
		observableParameter1_pAKT_ID8_HS746T_HM_15min = getObsOrSdParam(obsPar, mapObsParam)
		out[16] = observableParameter1_pAKT_ID8_HS746T_HM_15min
		return nothing
	end

	if observableId == "pAKT_ID8_HS746T_HM_1min" 
		observableParameter1_pAKT_ID8_HS746T_HM_1min = getObsOrSdParam(obsPar, mapObsParam)
		out[16] = observableParameter1_pAKT_ID8_HS746T_HM_1min
		return nothing
	end

	if observableId == "pAKT_ID8_HS746T_HM_240min" 
		observableParameter1_pAKT_ID8_HS746T_HM_240min = getObsOrSdParam(obsPar, mapObsParam)
		out[16] = observableParameter1_pAKT_ID8_HS746T_HM_240min
		return nothing
	end

	if observableId == "pAKT_ID8_HS746T_HM_3min" 
		observableParameter1_pAKT_ID8_HS746T_HM_3min = getObsOrSdParam(obsPar, mapObsParam)
		out[16] = observableParameter1_pAKT_ID8_HS746T_HM_3min
		return nothing
	end

	if observableId == "pAKT_ID9_HS746T_FM" 
		observableParameter1_pAKT_ID9_HS746T_FM = getObsOrSdParam(obsPar, mapObsParam)
		out[16] = observableParameter1_pAKT_ID9_HS746T_FM
		return nothing
	end

	if observableId == "pAKT_ID9_HS746T_HM_30EGF" 
		observableParameter1_pAKT_ID9_HS746T_HM_30EGF = getObsOrSdParam(obsPar, mapObsParam)
		out[16] = observableParameter1_pAKT_ID9_HS746T_HM_30EGF
		return nothing
	end

	if observableId == "pAKT_ID9_HS746T_HM_5EGF" 
		observableParameter1_pAKT_ID9_HS746T_HM_5EGF = getObsOrSdParam(obsPar, mapObsParam)
		out[16] = observableParameter1_pAKT_ID9_HS746T_HM_5EGF
		return nothing
	end

	if observableId == "pEGFR_ID10_HS746T_FM_EGF_CET" 
		observableParameter1_pEGFR_ID10_HS746T_FM_EGF_CET = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = 2observableParameter1_pEGFR_ID10_HS746T_FM_EGF_CET
		out[7] = 2observableParameter1_pEGFR_ID10_HS746T_FM_EGF_CET
		out[14] = observableParameter1_pEGFR_ID10_HS746T_FM_EGF_CET
		out[21] = observableParameter1_pEGFR_ID10_HS746T_FM_EGF_CET
		return nothing
	end

	if observableId == "pEGFR_ID10_HS746T_FM_EGF" 
		observableParameter1_pEGFR_ID10_HS746T_FM_EGF = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = 2observableParameter1_pEGFR_ID10_HS746T_FM_EGF
		out[7] = 2observableParameter1_pEGFR_ID10_HS746T_FM_EGF
		out[14] = observableParameter1_pEGFR_ID10_HS746T_FM_EGF
		out[21] = observableParameter1_pEGFR_ID10_HS746T_FM_EGF
		return nothing
	end

	if observableId == "pEGFR_ID12_MKN1_FM_HM" 
		observableParameter1_pEGFR_ID12_MKN1_FM_HM = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = 2observableParameter1_pEGFR_ID12_MKN1_FM_HM
		out[7] = 2observableParameter1_pEGFR_ID12_MKN1_FM_HM
		out[14] = observableParameter1_pEGFR_ID12_MKN1_FM_HM
		out[21] = observableParameter1_pEGFR_ID12_MKN1_FM_HM
		return nothing
	end

	if observableId == "pEGFR_ID2_HS746T_HM" 
		observableParameter1_pEGFR_ID2_HS746T_HM = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = 2observableParameter1_pEGFR_ID2_HS746T_HM
		out[7] = 2observableParameter1_pEGFR_ID2_HS746T_HM
		out[14] = observableParameter1_pEGFR_ID2_HS746T_HM
		out[21] = observableParameter1_pEGFR_ID2_HS746T_HM
		return nothing
	end

	if observableId == "pEGFR_ID2_MKN1_HM" 
		observableParameter1_pEGFR_ID2_MKN1_HM = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = 2observableParameter1_pEGFR_ID2_MKN1_HM
		out[7] = 2observableParameter1_pEGFR_ID2_MKN1_HM
		out[14] = observableParameter1_pEGFR_ID2_MKN1_HM
		out[21] = observableParameter1_pEGFR_ID2_MKN1_HM
		return nothing
	end

	if observableId == "pEGFR_ID3_HS746T_FM" 
		observableParameter1_pEGFR_ID3_HS746T_FM = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = 2observableParameter1_pEGFR_ID3_HS746T_FM
		out[7] = 2observableParameter1_pEGFR_ID3_HS746T_FM
		out[14] = observableParameter1_pEGFR_ID3_HS746T_FM
		out[21] = observableParameter1_pEGFR_ID3_HS746T_FM
		return nothing
	end

	if observableId == "pEGFR_ID3_MKN1_FM" 
		observableParameter1_pEGFR_ID3_MKN1_FM = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = 2observableParameter1_pEGFR_ID3_MKN1_FM
		out[7] = 2observableParameter1_pEGFR_ID3_MKN1_FM
		out[14] = observableParameter1_pEGFR_ID3_MKN1_FM
		out[21] = observableParameter1_pEGFR_ID3_MKN1_FM
		return nothing
	end

	if observableId == "pEGFR_ID4_HS746T_HM_3min" 
		observableParameter1_pEGFR_ID4_HS746T_HM_3min = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = 2observableParameter1_pEGFR_ID4_HS746T_HM_3min
		out[7] = 2observableParameter1_pEGFR_ID4_HS746T_HM_3min
		out[14] = observableParameter1_pEGFR_ID4_HS746T_HM_3min
		out[21] = observableParameter1_pEGFR_ID4_HS746T_HM_3min
		return nothing
	end

	if observableId == "pEGFR_ID4_MKN1_HM_3min" 
		observableParameter1_pEGFR_ID4_MKN1_HM_3min = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = 2observableParameter1_pEGFR_ID4_MKN1_HM_3min
		out[7] = 2observableParameter1_pEGFR_ID4_MKN1_HM_3min
		out[14] = observableParameter1_pEGFR_ID4_MKN1_HM_3min
		out[21] = observableParameter1_pEGFR_ID4_MKN1_HM_3min
		return nothing
	end

	if observableId == "pEGFR_ID5_MKN1_FM_15min" 
		observableParameter1_pEGFR_ID5_MKN1_FM_15min = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = 2observableParameter1_pEGFR_ID5_MKN1_FM_15min
		out[7] = 2observableParameter1_pEGFR_ID5_MKN1_FM_15min
		out[14] = observableParameter1_pEGFR_ID5_MKN1_FM_15min
		out[21] = observableParameter1_pEGFR_ID5_MKN1_FM_15min
		return nothing
	end

	if observableId == "pEGFR_ID5_MKN1_FM_1min" 
		observableParameter1_pEGFR_ID5_MKN1_FM_1min = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = 2observableParameter1_pEGFR_ID5_MKN1_FM_1min
		out[7] = 2observableParameter1_pEGFR_ID5_MKN1_FM_1min
		out[14] = observableParameter1_pEGFR_ID5_MKN1_FM_1min
		out[21] = observableParameter1_pEGFR_ID5_MKN1_FM_1min
		return nothing
	end

	if observableId == "pEGFR_ID5_MKN1_FM_240min" 
		observableParameter1_pEGFR_ID5_MKN1_FM_240min = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = 2observableParameter1_pEGFR_ID5_MKN1_FM_240min
		out[7] = 2observableParameter1_pEGFR_ID5_MKN1_FM_240min
		out[14] = observableParameter1_pEGFR_ID5_MKN1_FM_240min
		out[21] = observableParameter1_pEGFR_ID5_MKN1_FM_240min
		return nothing
	end

	if observableId == "pEGFR_ID5_MKN1_FM_3min" 
		observableParameter1_pEGFR_ID5_MKN1_FM_3min = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = 2observableParameter1_pEGFR_ID5_MKN1_FM_3min
		out[7] = 2observableParameter1_pEGFR_ID5_MKN1_FM_3min
		out[14] = observableParameter1_pEGFR_ID5_MKN1_FM_3min
		out[21] = observableParameter1_pEGFR_ID5_MKN1_FM_3min
		return nothing
	end

	if observableId == "pEGFR_ID6_MKN1_HM_15min" 
		observableParameter1_pEGFR_ID6_MKN1_HM_15min = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = 2observableParameter1_pEGFR_ID6_MKN1_HM_15min
		out[7] = 2observableParameter1_pEGFR_ID6_MKN1_HM_15min
		out[14] = observableParameter1_pEGFR_ID6_MKN1_HM_15min
		out[21] = observableParameter1_pEGFR_ID6_MKN1_HM_15min
		return nothing
	end

	if observableId == "pEGFR_ID6_MKN1_HM_1min" 
		observableParameter1_pEGFR_ID6_MKN1_HM_1min = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = 2observableParameter1_pEGFR_ID6_MKN1_HM_1min
		out[7] = 2observableParameter1_pEGFR_ID6_MKN1_HM_1min
		out[14] = observableParameter1_pEGFR_ID6_MKN1_HM_1min
		out[21] = observableParameter1_pEGFR_ID6_MKN1_HM_1min
		return nothing
	end

	if observableId == "pEGFR_ID6_MKN1_HM_240min" 
		observableParameter1_pEGFR_ID6_MKN1_HM_240min = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = 2observableParameter1_pEGFR_ID6_MKN1_HM_240min
		out[7] = 2observableParameter1_pEGFR_ID6_MKN1_HM_240min
		out[14] = observableParameter1_pEGFR_ID6_MKN1_HM_240min
		out[21] = observableParameter1_pEGFR_ID6_MKN1_HM_240min
		return nothing
	end

	if observableId == "pEGFR_ID6_MKN1_HM_3min" 
		observableParameter1_pEGFR_ID6_MKN1_HM_3min = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = 2observableParameter1_pEGFR_ID6_MKN1_HM_3min
		out[7] = 2observableParameter1_pEGFR_ID6_MKN1_HM_3min
		out[14] = observableParameter1_pEGFR_ID6_MKN1_HM_3min
		out[21] = observableParameter1_pEGFR_ID6_MKN1_HM_3min
		return nothing
	end

	if observableId == "pEGFR_ID7_HS746T_FM_15min" 
		observableParameter1_pEGFR_ID7_HS746T_FM_15min = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = 2observableParameter1_pEGFR_ID7_HS746T_FM_15min
		out[7] = 2observableParameter1_pEGFR_ID7_HS746T_FM_15min
		out[14] = observableParameter1_pEGFR_ID7_HS746T_FM_15min
		out[21] = observableParameter1_pEGFR_ID7_HS746T_FM_15min
		return nothing
	end

	if observableId == "pEGFR_ID7_HS746T_FM_1min" 
		observableParameter1_pEGFR_ID7_HS746T_FM_1min = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = 2observableParameter1_pEGFR_ID7_HS746T_FM_1min
		out[7] = 2observableParameter1_pEGFR_ID7_HS746T_FM_1min
		out[14] = observableParameter1_pEGFR_ID7_HS746T_FM_1min
		out[21] = observableParameter1_pEGFR_ID7_HS746T_FM_1min
		return nothing
	end

	if observableId == "pEGFR_ID7_HS746T_FM_240min" 
		observableParameter1_pEGFR_ID7_HS746T_FM_240min = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = 2observableParameter1_pEGFR_ID7_HS746T_FM_240min
		out[7] = 2observableParameter1_pEGFR_ID7_HS746T_FM_240min
		out[14] = observableParameter1_pEGFR_ID7_HS746T_FM_240min
		out[21] = observableParameter1_pEGFR_ID7_HS746T_FM_240min
		return nothing
	end

	if observableId == "pEGFR_ID7_HS746T_FM_3min" 
		observableParameter1_pEGFR_ID7_HS746T_FM_3min = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = 2observableParameter1_pEGFR_ID7_HS746T_FM_3min
		out[7] = 2observableParameter1_pEGFR_ID7_HS746T_FM_3min
		out[14] = observableParameter1_pEGFR_ID7_HS746T_FM_3min
		out[21] = observableParameter1_pEGFR_ID7_HS746T_FM_3min
		return nothing
	end

	if observableId == "pEGFR_ID8_HS746T_HM_15min" 
		observableParameter1_pEGFR_ID8_HS746T_HM_15min = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = 2observableParameter1_pEGFR_ID8_HS746T_HM_15min
		out[7] = 2observableParameter1_pEGFR_ID8_HS746T_HM_15min
		out[14] = observableParameter1_pEGFR_ID8_HS746T_HM_15min
		out[21] = observableParameter1_pEGFR_ID8_HS746T_HM_15min
		return nothing
	end

	if observableId == "pEGFR_ID8_HS746T_HM_1min" 
		observableParameter1_pEGFR_ID8_HS746T_HM_1min = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = 2observableParameter1_pEGFR_ID8_HS746T_HM_1min
		out[7] = 2observableParameter1_pEGFR_ID8_HS746T_HM_1min
		out[14] = observableParameter1_pEGFR_ID8_HS746T_HM_1min
		out[21] = observableParameter1_pEGFR_ID8_HS746T_HM_1min
		return nothing
	end

	if observableId == "pEGFR_ID8_HS746T_HM_240min" 
		observableParameter1_pEGFR_ID8_HS746T_HM_240min = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = 2observableParameter1_pEGFR_ID8_HS746T_HM_240min
		out[7] = 2observableParameter1_pEGFR_ID8_HS746T_HM_240min
		out[14] = observableParameter1_pEGFR_ID8_HS746T_HM_240min
		out[21] = observableParameter1_pEGFR_ID8_HS746T_HM_240min
		return nothing
	end

	if observableId == "pEGFR_ID8_HS746T_HM_3min" 
		observableParameter1_pEGFR_ID8_HS746T_HM_3min = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = 2observableParameter1_pEGFR_ID8_HS746T_HM_3min
		out[7] = 2observableParameter1_pEGFR_ID8_HS746T_HM_3min
		out[14] = observableParameter1_pEGFR_ID8_HS746T_HM_3min
		out[21] = observableParameter1_pEGFR_ID8_HS746T_HM_3min
		return nothing
	end

	if observableId == "pEGFR_ID9_HS746T_FM" 
		observableParameter1_pEGFR_ID9_HS746T_FM = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = 2observableParameter1_pEGFR_ID9_HS746T_FM
		out[7] = 2observableParameter1_pEGFR_ID9_HS746T_FM
		out[14] = observableParameter1_pEGFR_ID9_HS746T_FM
		out[21] = observableParameter1_pEGFR_ID9_HS746T_FM
		return nothing
	end

	if observableId == "pEGFR_ID9_HS746T_HM_30EGF" 
		observableParameter1_pEGFR_ID9_HS746T_HM_30EGF = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = 2observableParameter1_pEGFR_ID9_HS746T_HM_30EGF
		out[7] = 2observableParameter1_pEGFR_ID9_HS746T_HM_30EGF
		out[14] = observableParameter1_pEGFR_ID9_HS746T_HM_30EGF
		out[21] = observableParameter1_pEGFR_ID9_HS746T_HM_30EGF
		return nothing
	end

	if observableId == "pEGFR_ID9_HS746T_HM_5EGF" 
		observableParameter1_pEGFR_ID9_HS746T_HM_5EGF = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = 2observableParameter1_pEGFR_ID9_HS746T_HM_5EGF
		out[7] = 2observableParameter1_pEGFR_ID9_HS746T_HM_5EGF
		out[14] = observableParameter1_pEGFR_ID9_HS746T_HM_5EGF
		out[21] = observableParameter1_pEGFR_ID9_HS746T_HM_5EGF
		return nothing
	end

	if observableId == "pEGFR_ID9b_MKN1_HM_1EGF" 
		observableParameter1_pEGFR_ID9b_MKN1_HM_1EGF = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = 2observableParameter1_pEGFR_ID9b_MKN1_HM_1EGF
		out[7] = 2observableParameter1_pEGFR_ID9b_MKN1_HM_1EGF
		out[14] = observableParameter1_pEGFR_ID9b_MKN1_HM_1EGF
		out[21] = observableParameter1_pEGFR_ID9b_MKN1_HM_1EGF
		return nothing
	end

	if observableId == "pEGFR_ID9b_MKN1_HM_5EGF" 
		observableParameter1_pEGFR_ID9b_MKN1_HM_5EGF = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = 2observableParameter1_pEGFR_ID9b_MKN1_HM_5EGF
		out[7] = 2observableParameter1_pEGFR_ID9b_MKN1_HM_5EGF
		out[14] = observableParameter1_pEGFR_ID9b_MKN1_HM_5EGF
		out[21] = observableParameter1_pEGFR_ID9b_MKN1_HM_5EGF
		return nothing
	end

	if observableId == "pMAPK_ID10_HS746T_FM_EGF_CET" 
		observableParameter1_pMAPK_ID10_HS746T_FM_EGF_CET = getObsOrSdParam(obsPar, mapObsParam)
		out[4] = observableParameter1_pMAPK_ID10_HS746T_FM_EGF_CET
		return nothing
	end

	if observableId == "pMAPK_ID10_HS746T_FM_EGF" 
		observableParameter1_pMAPK_ID10_HS746T_FM_EGF = getObsOrSdParam(obsPar, mapObsParam)
		out[4] = observableParameter1_pMAPK_ID10_HS746T_FM_EGF
		return nothing
	end

	if observableId == "pMAPK_ID12_MKN1_FM_HM" 
		observableParameter1_pMAPK_ID12_MKN1_FM_HM = getObsOrSdParam(obsPar, mapObsParam)
		out[4] = observableParameter1_pMAPK_ID12_MKN1_FM_HM
		return nothing
	end

	if observableId == "pMAPK_ID2_HS746T_HM" 
		observableParameter1_pMAPK_ID2_HS746T_HM = getObsOrSdParam(obsPar, mapObsParam)
		out[4] = observableParameter1_pMAPK_ID2_HS746T_HM
		return nothing
	end

	if observableId == "pMAPK_ID2_MKN1_HM" 
		observableParameter1_pMAPK_ID2_MKN1_HM = getObsOrSdParam(obsPar, mapObsParam)
		out[4] = observableParameter1_pMAPK_ID2_MKN1_HM
		return nothing
	end

	if observableId == "pMAPK_ID5_MKN1_FM_15min" 
		observableParameter1_pMAPK_ID5_MKN1_FM_15min = getObsOrSdParam(obsPar, mapObsParam)
		out[4] = observableParameter1_pMAPK_ID5_MKN1_FM_15min
		return nothing
	end

	if observableId == "pMAPK_ID5_MKN1_FM_1min" 
		observableParameter1_pMAPK_ID5_MKN1_FM_1min = getObsOrSdParam(obsPar, mapObsParam)
		out[4] = observableParameter1_pMAPK_ID5_MKN1_FM_1min
		return nothing
	end

	if observableId == "pMAPK_ID5_MKN1_FM_240min" 
		observableParameter1_pMAPK_ID5_MKN1_FM_240min = getObsOrSdParam(obsPar, mapObsParam)
		out[4] = observableParameter1_pMAPK_ID5_MKN1_FM_240min
		return nothing
	end

	if observableId == "pMAPK_ID5_MKN1_FM_3min" 
		observableParameter1_pMAPK_ID5_MKN1_FM_3min = getObsOrSdParam(obsPar, mapObsParam)
		out[4] = observableParameter1_pMAPK_ID5_MKN1_FM_3min
		return nothing
	end

	if observableId == "pMAPK_ID7_HS746T_FM_15min" 
		observableParameter1_pMAPK_ID7_HS746T_FM_15min = getObsOrSdParam(obsPar, mapObsParam)
		out[4] = observableParameter1_pMAPK_ID7_HS746T_FM_15min
		return nothing
	end

	if observableId == "pMAPK_ID7_HS746T_FM_1min" 
		observableParameter1_pMAPK_ID7_HS746T_FM_1min = getObsOrSdParam(obsPar, mapObsParam)
		out[4] = observableParameter1_pMAPK_ID7_HS746T_FM_1min
		return nothing
	end

	if observableId == "pMAPK_ID7_HS746T_FM_240min" 
		observableParameter1_pMAPK_ID7_HS746T_FM_240min = getObsOrSdParam(obsPar, mapObsParam)
		out[4] = observableParameter1_pMAPK_ID7_HS746T_FM_240min
		return nothing
	end

	if observableId == "pMAPK_ID7_HS746T_FM_3min" 
		observableParameter1_pMAPK_ID7_HS746T_FM_3min = getObsOrSdParam(obsPar, mapObsParam)
		out[4] = observableParameter1_pMAPK_ID7_HS746T_FM_3min
		return nothing
	end

	if observableId == "pMAPK_ID8_HS746T_HM_15min" 
		observableParameter1_pMAPK_ID8_HS746T_HM_15min = getObsOrSdParam(obsPar, mapObsParam)
		out[4] = observableParameter1_pMAPK_ID8_HS746T_HM_15min
		return nothing
	end

	if observableId == "pMAPK_ID8_HS746T_HM_1min" 
		observableParameter1_pMAPK_ID8_HS746T_HM_1min = getObsOrSdParam(obsPar, mapObsParam)
		out[4] = observableParameter1_pMAPK_ID8_HS746T_HM_1min
		return nothing
	end

	if observableId == "pMAPK_ID8_HS746T_HM_240min" 
		observableParameter1_pMAPK_ID8_HS746T_HM_240min = getObsOrSdParam(obsPar, mapObsParam)
		out[4] = observableParameter1_pMAPK_ID8_HS746T_HM_240min
		return nothing
	end

	if observableId == "pMAPK_ID8_HS746T_HM_3min" 
		observableParameter1_pMAPK_ID8_HS746T_HM_3min = getObsOrSdParam(obsPar, mapObsParam)
		out[4] = observableParameter1_pMAPK_ID8_HS746T_HM_3min
		return nothing
	end

	if observableId == "pMAPK_ID9_HS746T_FM" 
		observableParameter1_pMAPK_ID9_HS746T_FM = getObsOrSdParam(obsPar, mapObsParam)
		out[4] = observableParameter1_pMAPK_ID9_HS746T_FM
		return nothing
	end

	if observableId == "pMAPK_ID9_HS746T_HM_30EGF" 
		observableParameter1_pMAPK_ID9_HS746T_HM_30EGF = getObsOrSdParam(obsPar, mapObsParam)
		out[4] = observableParameter1_pMAPK_ID9_HS746T_HM_30EGF
		return nothing
	end

	if observableId == "pMAPK_ID9_HS746T_HM_5EGF" 
		observableParameter1_pMAPK_ID9_HS746T_HM_5EGF = getObsOrSdParam(obsPar, mapObsParam)
		out[4] = observableParameter1_pMAPK_ID9_HS746T_HM_5EGF
		return nothing
	end

	if observableId == "pMAPK_ID9b_MKN1_HM_1EGF" 
		observableParameter1_pMAPK_ID9b_MKN1_HM_1EGF = getObsOrSdParam(obsPar, mapObsParam)
		out[4] = observableParameter1_pMAPK_ID9b_MKN1_HM_1EGF
		return nothing
	end

	if observableId == "pMAPK_ID9b_MKN1_HM_5EGF" 
		observableParameter1_pMAPK_ID9b_MKN1_HM_5EGF = getObsOrSdParam(obsPar, mapObsParam)
		out[4] = observableParameter1_pMAPK_ID9b_MKN1_HM_5EGF
		return nothing
	end

end

function evalDYmodDp(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	EGFR, pPI3K, MMET_METinh_2, pMAPK, MMET_2, pEGFR_EGF_2, pEGFR_EGF_2_i, pMPI3K, EGFR_EGF_2, MMET, pMMET_2, RAS_GTP, pMMET_2_i, pMMET_pEGFR_i, MMET_MMET_METinh, pAKT, EGFR_EGF, MMET_METinh, MMET_EGFR, EGFR_CET, pMMET_pEGFR, EGFR_MMET_METinh= u 
	d_ksyn_EGFR__fm_2_hm, xi_kexp_pMMET_2_i, d_ka_RAS__pEGFR_EGF_2__MKN1_2_HS746T, CET_bool1, EGF_bool1, d_ka_MAPK__MKN1_2_HS746T, ki_RAS__MKN1, d_ka_AKT__MKN1_2_HS746T, KD_EGFR_CET, xi_ka_PI3K_pMMET_2, d_MAPKtotal__fm_2_hm, d_ksyn_MMET__fm_2_hm, ka_PI3K__basal, d_kimp_pEGFR_EGF_2__MKN1_2_HS746T, d_AKTtotal__fm_2_hm, KD_METinh, xi_ki_MPI3K, d_kexp_pEGFR_EGF_2_i__MKN1_2_HS746T, d_ksyn_EGFR__MKN1_2_HS746T, d_MPI3Ktotal__fm_2_hm, d_kdeg_membran__MKN1_2_HS746T, MPI3Ktotal__MKN1_fm, KD_EGFR_EGF, METinh_bool1, kexp_pEGFR_EGF_2_i__MKN1, CET_level, kbin_EGFR_EGF, ksyn_EGFR__MKN1_fm, xi_kdeg_pMMET_pEGFR_i, d_ki_PI3K__MKN1_2_HS746T, xi_kpho_MMET_EGFR, MAPKtotal__MKN1_fm, xi_kexp_pMMET_pEGFR_i, d_RAStotal__MKN1_2_HS746T, kpho_EGFR_EGF, ka_MAPK__MKN1, kpho_MMET_EGFR, full_medium, xi_ka_RAS_pMMET_2, xi_kimp_pMMET_pEGFR, xi_kdeg_pMMET_2_i, ka_RAS__basal__MKN1, d_RAStotal__fm_2_hm, kdim_MMET, MKN1, ksyn_MMET__HS746T_fm, kdim_MMET_EGFR, ka_PI3K__pEGFR_EGF_2, kdeg_membran__MKN1, HS746T, d_ki_RAS__MKN1_2_HS746T, EGF_level, kdim_EGFR_EGF, kbin_EGFR_CET, ka_AKT__MKN1, d_ki_AKT__MKN1_2_HS746T, PI3Ktotal__HS746T_fm, ki_PI3K__MKN1, kpho_MMET, AKTtotal__MKN1_fm, xi_kimp_pMMET_2, relative_ksyn_EGFR, ka_RAS__pEGFR_EGF_2__MKN1, ki_MAPK, cyt, kimp_pEGFR_EGF_2__MKN1, d_AKTtotal__MKN1_2_HS746T, ki_AKT__MKN1, RAStotal__MKN1_fm, d_PI3Ktotal__fm_2_hm, xi_ka_PI3K_pMMET_pEGFR, xi_kdim_MMET_EGFR, xi_kpho_MMET, METinh_level, kdim_MMETinh, d_MAPKtotal__MKN1_2_HS746T, xi_ka_RAS_pMMET_pEGFR, kdeg_pEGFR_EGF_2_i__MKN1, xi_kdim_MMET, d_kdeg_pEGFR_EGF_2_i__MKN1_2_HS746T, d_ka_RAS__basal__MKN1_2_HS746T, hunger_medium = p 
	if observableId == "EGFR_ID10_HS746T_FM_EGF_CET" 
		return nothing
	end

	if observableId == "EGFR_ID10_HS746T_FM_EGF" 
		return nothing
	end

	if observableId == "EGFR_ID12_MKN1_FM_HM" 
		return nothing
	end

	if observableId == "EGFR_ID3_HS746T_FM" 
		return nothing
	end

	if observableId == "EGFR_ID3_MKN1_FM" 
		return nothing
	end

	if observableId == "EGFR_ID4_HS746T_HM_3min" 
		return nothing
	end

	if observableId == "EGFR_ID4_MKN1_HM_3min" 
		return nothing
	end

	if observableId == "MMET_ID10_HS746T_FM_EGF_CET" 
		return nothing
	end

	if observableId == "MMET_ID10_HS746T_FM_EGF" 
		return nothing
	end

	if observableId == "pAKT_ID12_MKN1_FM_HM" 
		return nothing
	end

	if observableId == "pAKT_ID2_HS746T_HM" 
		return nothing
	end

	if observableId == "pAKT_ID2_MKN1_HM" 
		return nothing
	end

	if observableId == "pAKT_ID5_MKN1_FM_15min" 
		return nothing
	end

	if observableId == "pAKT_ID5_MKN1_FM_1min" 
		return nothing
	end

	if observableId == "pAKT_ID5_MKN1_FM_240min" 
		return nothing
	end

	if observableId == "pAKT_ID5_MKN1_FM_3min" 
		return nothing
	end

	if observableId == "pAKT_ID7_HS746T_FM_15min" 
		return nothing
	end

	if observableId == "pAKT_ID7_HS746T_FM_1min" 
		return nothing
	end

	if observableId == "pAKT_ID7_HS746T_FM_240min" 
		return nothing
	end

	if observableId == "pAKT_ID7_HS746T_FM_3min" 
		return nothing
	end

	if observableId == "pAKT_ID8_HS746T_HM_15min" 
		return nothing
	end

	if observableId == "pAKT_ID8_HS746T_HM_1min" 
		return nothing
	end

	if observableId == "pAKT_ID8_HS746T_HM_240min" 
		return nothing
	end

	if observableId == "pAKT_ID8_HS746T_HM_3min" 
		return nothing
	end

	if observableId == "pAKT_ID9_HS746T_FM" 
		return nothing
	end

	if observableId == "pAKT_ID9_HS746T_HM_30EGF" 
		return nothing
	end

	if observableId == "pAKT_ID9_HS746T_HM_5EGF" 
		return nothing
	end

	if observableId == "pEGFR_ID10_HS746T_FM_EGF_CET" 
		return nothing
	end

	if observableId == "pEGFR_ID10_HS746T_FM_EGF" 
		return nothing
	end

	if observableId == "pEGFR_ID12_MKN1_FM_HM" 
		return nothing
	end

	if observableId == "pEGFR_ID2_HS746T_HM" 
		return nothing
	end

	if observableId == "pEGFR_ID2_MKN1_HM" 
		return nothing
	end

	if observableId == "pEGFR_ID3_HS746T_FM" 
		return nothing
	end

	if observableId == "pEGFR_ID3_MKN1_FM" 
		return nothing
	end

	if observableId == "pEGFR_ID4_HS746T_HM_3min" 
		return nothing
	end

	if observableId == "pEGFR_ID4_MKN1_HM_3min" 
		return nothing
	end

	if observableId == "pEGFR_ID5_MKN1_FM_15min" 
		return nothing
	end

	if observableId == "pEGFR_ID5_MKN1_FM_1min" 
		return nothing
	end

	if observableId == "pEGFR_ID5_MKN1_FM_240min" 
		return nothing
	end

	if observableId == "pEGFR_ID5_MKN1_FM_3min" 
		return nothing
	end

	if observableId == "pEGFR_ID6_MKN1_HM_15min" 
		return nothing
	end

	if observableId == "pEGFR_ID6_MKN1_HM_1min" 
		return nothing
	end

	if observableId == "pEGFR_ID6_MKN1_HM_240min" 
		return nothing
	end

	if observableId == "pEGFR_ID6_MKN1_HM_3min" 
		return nothing
	end

	if observableId == "pEGFR_ID7_HS746T_FM_15min" 
		return nothing
	end

	if observableId == "pEGFR_ID7_HS746T_FM_1min" 
		return nothing
	end

	if observableId == "pEGFR_ID7_HS746T_FM_240min" 
		return nothing
	end

	if observableId == "pEGFR_ID7_HS746T_FM_3min" 
		return nothing
	end

	if observableId == "pEGFR_ID8_HS746T_HM_15min" 
		return nothing
	end

	if observableId == "pEGFR_ID8_HS746T_HM_1min" 
		return nothing
	end

	if observableId == "pEGFR_ID8_HS746T_HM_240min" 
		return nothing
	end

	if observableId == "pEGFR_ID8_HS746T_HM_3min" 
		return nothing
	end

	if observableId == "pEGFR_ID9_HS746T_FM" 
		return nothing
	end

	if observableId == "pEGFR_ID9_HS746T_HM_30EGF" 
		return nothing
	end

	if observableId == "pEGFR_ID9_HS746T_HM_5EGF" 
		return nothing
	end

	if observableId == "pEGFR_ID9b_MKN1_HM_1EGF" 
		return nothing
	end

	if observableId == "pEGFR_ID9b_MKN1_HM_5EGF" 
		return nothing
	end

	if observableId == "pMAPK_ID10_HS746T_FM_EGF_CET" 
		return nothing
	end

	if observableId == "pMAPK_ID10_HS746T_FM_EGF" 
		return nothing
	end

	if observableId == "pMAPK_ID12_MKN1_FM_HM" 
		return nothing
	end

	if observableId == "pMAPK_ID2_HS746T_HM" 
		return nothing
	end

	if observableId == "pMAPK_ID2_MKN1_HM" 
		return nothing
	end

	if observableId == "pMAPK_ID5_MKN1_FM_15min" 
		return nothing
	end

	if observableId == "pMAPK_ID5_MKN1_FM_1min" 
		return nothing
	end

	if observableId == "pMAPK_ID5_MKN1_FM_240min" 
		return nothing
	end

	if observableId == "pMAPK_ID5_MKN1_FM_3min" 
		return nothing
	end

	if observableId == "pMAPK_ID7_HS746T_FM_15min" 
		return nothing
	end

	if observableId == "pMAPK_ID7_HS746T_FM_1min" 
		return nothing
	end

	if observableId == "pMAPK_ID7_HS746T_FM_240min" 
		return nothing
	end

	if observableId == "pMAPK_ID7_HS746T_FM_3min" 
		return nothing
	end

	if observableId == "pMAPK_ID8_HS746T_HM_15min" 
		return nothing
	end

	if observableId == "pMAPK_ID8_HS746T_HM_1min" 
		return nothing
	end

	if observableId == "pMAPK_ID8_HS746T_HM_240min" 
		return nothing
	end

	if observableId == "pMAPK_ID8_HS746T_HM_3min" 
		return nothing
	end

	if observableId == "pMAPK_ID9_HS746T_FM" 
		return nothing
	end

	if observableId == "pMAPK_ID9_HS746T_HM_30EGF" 
		return nothing
	end

	if observableId == "pMAPK_ID9_HS746T_HM_5EGF" 
		return nothing
	end

	if observableId == "pMAPK_ID9b_MKN1_HM_1EGF" 
		return nothing
	end

	if observableId == "pMAPK_ID9b_MKN1_HM_5EGF" 
		return nothing
	end

end

function evalDSdDu!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	EGFR, pPI3K, MMET_METinh_2, pMAPK, MMET_2, pEGFR_EGF_2, pEGFR_EGF_2_i, pMPI3K, EGFR_EGF_2, MMET, pMMET_2, RAS_GTP, pMMET_2_i, pMMET_pEGFR_i, MMET_MMET_METinh, pAKT, EGFR_EGF, MMET_METinh, MMET_EGFR, EGFR_CET, pMMET_pEGFR, EGFR_MMET_METinh= u 
	d_ksyn_EGFR__fm_2_hm, xi_kexp_pMMET_2_i, d_ka_RAS__pEGFR_EGF_2__MKN1_2_HS746T, CET_bool1, EGF_bool1, d_ka_MAPK__MKN1_2_HS746T, ki_RAS__MKN1, d_ka_AKT__MKN1_2_HS746T, KD_EGFR_CET, xi_ka_PI3K_pMMET_2, d_MAPKtotal__fm_2_hm, d_ksyn_MMET__fm_2_hm, ka_PI3K__basal, d_kimp_pEGFR_EGF_2__MKN1_2_HS746T, d_AKTtotal__fm_2_hm, KD_METinh, xi_ki_MPI3K, d_kexp_pEGFR_EGF_2_i__MKN1_2_HS746T, d_ksyn_EGFR__MKN1_2_HS746T, d_MPI3Ktotal__fm_2_hm, d_kdeg_membran__MKN1_2_HS746T, MPI3Ktotal__MKN1_fm, KD_EGFR_EGF, METinh_bool1, kexp_pEGFR_EGF_2_i__MKN1, CET_level, kbin_EGFR_EGF, ksyn_EGFR__MKN1_fm, xi_kdeg_pMMET_pEGFR_i, d_ki_PI3K__MKN1_2_HS746T, xi_kpho_MMET_EGFR, MAPKtotal__MKN1_fm, xi_kexp_pMMET_pEGFR_i, d_RAStotal__MKN1_2_HS746T, kpho_EGFR_EGF, ka_MAPK__MKN1, kpho_MMET_EGFR, full_medium, xi_ka_RAS_pMMET_2, xi_kimp_pMMET_pEGFR, xi_kdeg_pMMET_2_i, ka_RAS__basal__MKN1, d_RAStotal__fm_2_hm, kdim_MMET, MKN1, ksyn_MMET__HS746T_fm, kdim_MMET_EGFR, ka_PI3K__pEGFR_EGF_2, kdeg_membran__MKN1, HS746T, d_ki_RAS__MKN1_2_HS746T, EGF_level, kdim_EGFR_EGF, kbin_EGFR_CET, ka_AKT__MKN1, d_ki_AKT__MKN1_2_HS746T, PI3Ktotal__HS746T_fm, ki_PI3K__MKN1, kpho_MMET, AKTtotal__MKN1_fm, xi_kimp_pMMET_2, relative_ksyn_EGFR, ka_RAS__pEGFR_EGF_2__MKN1, ki_MAPK, cyt, kimp_pEGFR_EGF_2__MKN1, d_AKTtotal__MKN1_2_HS746T, ki_AKT__MKN1, RAStotal__MKN1_fm, d_PI3Ktotal__fm_2_hm, xi_ka_PI3K_pMMET_pEGFR, xi_kdim_MMET_EGFR, xi_kpho_MMET, METinh_level, kdim_MMETinh, d_MAPKtotal__MKN1_2_HS746T, xi_ka_RAS_pMMET_pEGFR, kdeg_pEGFR_EGF_2_i__MKN1, xi_kdim_MMET, d_kdeg_pEGFR_EGF_2_i__MKN1_2_HS746T, d_ka_RAS__basal__MKN1_2_HS746T, hunger_medium = p 
	if observableId == "EGFR_ID10_HS746T_FM_EGF_CET" 
		return nothing
	end

	if observableId == "EGFR_ID10_HS746T_FM_EGF" 
		return nothing
	end

	if observableId == "EGFR_ID12_MKN1_FM_HM" 
		return nothing
	end

	if observableId == "EGFR_ID3_HS746T_FM" 
		return nothing
	end

	if observableId == "EGFR_ID3_MKN1_FM" 
		return nothing
	end

	if observableId == "EGFR_ID4_HS746T_HM_3min" 
		return nothing
	end

	if observableId == "EGFR_ID4_MKN1_HM_3min" 
		return nothing
	end

	if observableId == "MMET_ID10_HS746T_FM_EGF_CET" 
		return nothing
	end

	if observableId == "MMET_ID10_HS746T_FM_EGF" 
		return nothing
	end

	if observableId == "pAKT_ID12_MKN1_FM_HM" 
		return nothing
	end

	if observableId == "pAKT_ID2_HS746T_HM" 
		return nothing
	end

	if observableId == "pAKT_ID2_MKN1_HM" 
		return nothing
	end

	if observableId == "pAKT_ID5_MKN1_FM_15min" 
		return nothing
	end

	if observableId == "pAKT_ID5_MKN1_FM_1min" 
		return nothing
	end

	if observableId == "pAKT_ID5_MKN1_FM_240min" 
		return nothing
	end

	if observableId == "pAKT_ID5_MKN1_FM_3min" 
		return nothing
	end

	if observableId == "pAKT_ID7_HS746T_FM_15min" 
		return nothing
	end

	if observableId == "pAKT_ID7_HS746T_FM_1min" 
		return nothing
	end

	if observableId == "pAKT_ID7_HS746T_FM_240min" 
		return nothing
	end

	if observableId == "pAKT_ID7_HS746T_FM_3min" 
		return nothing
	end

	if observableId == "pAKT_ID8_HS746T_HM_15min" 
		return nothing
	end

	if observableId == "pAKT_ID8_HS746T_HM_1min" 
		return nothing
	end

	if observableId == "pAKT_ID8_HS746T_HM_240min" 
		return nothing
	end

	if observableId == "pAKT_ID8_HS746T_HM_3min" 
		return nothing
	end

	if observableId == "pAKT_ID9_HS746T_FM" 
		return nothing
	end

	if observableId == "pAKT_ID9_HS746T_HM_30EGF" 
		return nothing
	end

	if observableId == "pAKT_ID9_HS746T_HM_5EGF" 
		return nothing
	end

	if observableId == "pEGFR_ID10_HS746T_FM_EGF_CET" 
		return nothing
	end

	if observableId == "pEGFR_ID10_HS746T_FM_EGF" 
		return nothing
	end

	if observableId == "pEGFR_ID12_MKN1_FM_HM" 
		return nothing
	end

	if observableId == "pEGFR_ID2_HS746T_HM" 
		return nothing
	end

	if observableId == "pEGFR_ID2_MKN1_HM" 
		return nothing
	end

	if observableId == "pEGFR_ID3_HS746T_FM" 
		return nothing
	end

	if observableId == "pEGFR_ID3_MKN1_FM" 
		return nothing
	end

	if observableId == "pEGFR_ID4_HS746T_HM_3min" 
		return nothing
	end

	if observableId == "pEGFR_ID4_MKN1_HM_3min" 
		return nothing
	end

	if observableId == "pEGFR_ID5_MKN1_FM_15min" 
		return nothing
	end

	if observableId == "pEGFR_ID5_MKN1_FM_1min" 
		return nothing
	end

	if observableId == "pEGFR_ID5_MKN1_FM_240min" 
		return nothing
	end

	if observableId == "pEGFR_ID5_MKN1_FM_3min" 
		return nothing
	end

	if observableId == "pEGFR_ID6_MKN1_HM_15min" 
		return nothing
	end

	if observableId == "pEGFR_ID6_MKN1_HM_1min" 
		return nothing
	end

	if observableId == "pEGFR_ID6_MKN1_HM_240min" 
		return nothing
	end

	if observableId == "pEGFR_ID6_MKN1_HM_3min" 
		return nothing
	end

	if observableId == "pEGFR_ID7_HS746T_FM_15min" 
		return nothing
	end

	if observableId == "pEGFR_ID7_HS746T_FM_1min" 
		return nothing
	end

	if observableId == "pEGFR_ID7_HS746T_FM_240min" 
		return nothing
	end

	if observableId == "pEGFR_ID7_HS746T_FM_3min" 
		return nothing
	end

	if observableId == "pEGFR_ID8_HS746T_HM_15min" 
		return nothing
	end

	if observableId == "pEGFR_ID8_HS746T_HM_1min" 
		return nothing
	end

	if observableId == "pEGFR_ID8_HS746T_HM_240min" 
		return nothing
	end

	if observableId == "pEGFR_ID8_HS746T_HM_3min" 
		return nothing
	end

	if observableId == "pEGFR_ID9_HS746T_FM" 
		return nothing
	end

	if observableId == "pEGFR_ID9_HS746T_HM_30EGF" 
		return nothing
	end

	if observableId == "pEGFR_ID9_HS746T_HM_5EGF" 
		return nothing
	end

	if observableId == "pEGFR_ID9b_MKN1_HM_1EGF" 
		return nothing
	end

	if observableId == "pEGFR_ID9b_MKN1_HM_5EGF" 
		return nothing
	end

	if observableId == "pMAPK_ID10_HS746T_FM_EGF_CET" 
		return nothing
	end

	if observableId == "pMAPK_ID10_HS746T_FM_EGF" 
		return nothing
	end

	if observableId == "pMAPK_ID12_MKN1_FM_HM" 
		return nothing
	end

	if observableId == "pMAPK_ID2_HS746T_HM" 
		return nothing
	end

	if observableId == "pMAPK_ID2_MKN1_HM" 
		return nothing
	end

	if observableId == "pMAPK_ID5_MKN1_FM_15min" 
		return nothing
	end

	if observableId == "pMAPK_ID5_MKN1_FM_1min" 
		return nothing
	end

	if observableId == "pMAPK_ID5_MKN1_FM_240min" 
		return nothing
	end

	if observableId == "pMAPK_ID5_MKN1_FM_3min" 
		return nothing
	end

	if observableId == "pMAPK_ID7_HS746T_FM_15min" 
		return nothing
	end

	if observableId == "pMAPK_ID7_HS746T_FM_1min" 
		return nothing
	end

	if observableId == "pMAPK_ID7_HS746T_FM_240min" 
		return nothing
	end

	if observableId == "pMAPK_ID7_HS746T_FM_3min" 
		return nothing
	end

	if observableId == "pMAPK_ID8_HS746T_HM_15min" 
		return nothing
	end

	if observableId == "pMAPK_ID8_HS746T_HM_1min" 
		return nothing
	end

	if observableId == "pMAPK_ID8_HS746T_HM_240min" 
		return nothing
	end

	if observableId == "pMAPK_ID8_HS746T_HM_3min" 
		return nothing
	end

	if observableId == "pMAPK_ID9_HS746T_FM" 
		return nothing
	end

	if observableId == "pMAPK_ID9_HS746T_HM_30EGF" 
		return nothing
	end

	if observableId == "pMAPK_ID9_HS746T_HM_5EGF" 
		return nothing
	end

	if observableId == "pMAPK_ID9b_MKN1_HM_1EGF" 
		return nothing
	end

	if observableId == "pMAPK_ID9b_MKN1_HM_5EGF" 
		return nothing
	end

end

function evalDSdDp!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	EGFR, pPI3K, MMET_METinh_2, pMAPK, MMET_2, pEGFR_EGF_2, pEGFR_EGF_2_i, pMPI3K, EGFR_EGF_2, MMET, pMMET_2, RAS_GTP, pMMET_2_i, pMMET_pEGFR_i, MMET_MMET_METinh, pAKT, EGFR_EGF, MMET_METinh, MMET_EGFR, EGFR_CET, pMMET_pEGFR, EGFR_MMET_METinh= u 
	d_ksyn_EGFR__fm_2_hm, xi_kexp_pMMET_2_i, d_ka_RAS__pEGFR_EGF_2__MKN1_2_HS746T, CET_bool1, EGF_bool1, d_ka_MAPK__MKN1_2_HS746T, ki_RAS__MKN1, d_ka_AKT__MKN1_2_HS746T, KD_EGFR_CET, xi_ka_PI3K_pMMET_2, d_MAPKtotal__fm_2_hm, d_ksyn_MMET__fm_2_hm, ka_PI3K__basal, d_kimp_pEGFR_EGF_2__MKN1_2_HS746T, d_AKTtotal__fm_2_hm, KD_METinh, xi_ki_MPI3K, d_kexp_pEGFR_EGF_2_i__MKN1_2_HS746T, d_ksyn_EGFR__MKN1_2_HS746T, d_MPI3Ktotal__fm_2_hm, d_kdeg_membran__MKN1_2_HS746T, MPI3Ktotal__MKN1_fm, KD_EGFR_EGF, METinh_bool1, kexp_pEGFR_EGF_2_i__MKN1, CET_level, kbin_EGFR_EGF, ksyn_EGFR__MKN1_fm, xi_kdeg_pMMET_pEGFR_i, d_ki_PI3K__MKN1_2_HS746T, xi_kpho_MMET_EGFR, MAPKtotal__MKN1_fm, xi_kexp_pMMET_pEGFR_i, d_RAStotal__MKN1_2_HS746T, kpho_EGFR_EGF, ka_MAPK__MKN1, kpho_MMET_EGFR, full_medium, xi_ka_RAS_pMMET_2, xi_kimp_pMMET_pEGFR, xi_kdeg_pMMET_2_i, ka_RAS__basal__MKN1, d_RAStotal__fm_2_hm, kdim_MMET, MKN1, ksyn_MMET__HS746T_fm, kdim_MMET_EGFR, ka_PI3K__pEGFR_EGF_2, kdeg_membran__MKN1, HS746T, d_ki_RAS__MKN1_2_HS746T, EGF_level, kdim_EGFR_EGF, kbin_EGFR_CET, ka_AKT__MKN1, d_ki_AKT__MKN1_2_HS746T, PI3Ktotal__HS746T_fm, ki_PI3K__MKN1, kpho_MMET, AKTtotal__MKN1_fm, xi_kimp_pMMET_2, relative_ksyn_EGFR, ka_RAS__pEGFR_EGF_2__MKN1, ki_MAPK, cyt, kimp_pEGFR_EGF_2__MKN1, d_AKTtotal__MKN1_2_HS746T, ki_AKT__MKN1, RAStotal__MKN1_fm, d_PI3Ktotal__fm_2_hm, xi_ka_PI3K_pMMET_pEGFR, xi_kdim_MMET_EGFR, xi_kpho_MMET, METinh_level, kdim_MMETinh, d_MAPKtotal__MKN1_2_HS746T, xi_ka_RAS_pMMET_pEGFR, kdeg_pEGFR_EGF_2_i__MKN1, xi_kdim_MMET, d_kdeg_pEGFR_EGF_2_i__MKN1_2_HS746T, d_ka_RAS__basal__MKN1_2_HS746T, hunger_medium = p 
	if observableId == "EGFR_ID10_HS746T_FM_EGF_CET" 
		return nothing
	end

	if observableId == "EGFR_ID10_HS746T_FM_EGF" 
		return nothing
	end

	if observableId == "EGFR_ID12_MKN1_FM_HM" 
		return nothing
	end

	if observableId == "EGFR_ID3_HS746T_FM" 
		return nothing
	end

	if observableId == "EGFR_ID3_MKN1_FM" 
		return nothing
	end

	if observableId == "EGFR_ID4_HS746T_HM_3min" 
		return nothing
	end

	if observableId == "EGFR_ID4_MKN1_HM_3min" 
		return nothing
	end

	if observableId == "MMET_ID10_HS746T_FM_EGF_CET" 
		return nothing
	end

	if observableId == "MMET_ID10_HS746T_FM_EGF" 
		return nothing
	end

	if observableId == "pAKT_ID12_MKN1_FM_HM" 
		return nothing
	end

	if observableId == "pAKT_ID2_HS746T_HM" 
		return nothing
	end

	if observableId == "pAKT_ID2_MKN1_HM" 
		return nothing
	end

	if observableId == "pAKT_ID5_MKN1_FM_15min" 
		return nothing
	end

	if observableId == "pAKT_ID5_MKN1_FM_1min" 
		return nothing
	end

	if observableId == "pAKT_ID5_MKN1_FM_240min" 
		return nothing
	end

	if observableId == "pAKT_ID5_MKN1_FM_3min" 
		return nothing
	end

	if observableId == "pAKT_ID7_HS746T_FM_15min" 
		return nothing
	end

	if observableId == "pAKT_ID7_HS746T_FM_1min" 
		return nothing
	end

	if observableId == "pAKT_ID7_HS746T_FM_240min" 
		return nothing
	end

	if observableId == "pAKT_ID7_HS746T_FM_3min" 
		return nothing
	end

	if observableId == "pAKT_ID8_HS746T_HM_15min" 
		return nothing
	end

	if observableId == "pAKT_ID8_HS746T_HM_1min" 
		return nothing
	end

	if observableId == "pAKT_ID8_HS746T_HM_240min" 
		return nothing
	end

	if observableId == "pAKT_ID8_HS746T_HM_3min" 
		return nothing
	end

	if observableId == "pAKT_ID9_HS746T_FM" 
		return nothing
	end

	if observableId == "pAKT_ID9_HS746T_HM_30EGF" 
		return nothing
	end

	if observableId == "pAKT_ID9_HS746T_HM_5EGF" 
		return nothing
	end

	if observableId == "pEGFR_ID10_HS746T_FM_EGF_CET" 
		return nothing
	end

	if observableId == "pEGFR_ID10_HS746T_FM_EGF" 
		return nothing
	end

	if observableId == "pEGFR_ID12_MKN1_FM_HM" 
		return nothing
	end

	if observableId == "pEGFR_ID2_HS746T_HM" 
		return nothing
	end

	if observableId == "pEGFR_ID2_MKN1_HM" 
		return nothing
	end

	if observableId == "pEGFR_ID3_HS746T_FM" 
		return nothing
	end

	if observableId == "pEGFR_ID3_MKN1_FM" 
		return nothing
	end

	if observableId == "pEGFR_ID4_HS746T_HM_3min" 
		return nothing
	end

	if observableId == "pEGFR_ID4_MKN1_HM_3min" 
		return nothing
	end

	if observableId == "pEGFR_ID5_MKN1_FM_15min" 
		return nothing
	end

	if observableId == "pEGFR_ID5_MKN1_FM_1min" 
		return nothing
	end

	if observableId == "pEGFR_ID5_MKN1_FM_240min" 
		return nothing
	end

	if observableId == "pEGFR_ID5_MKN1_FM_3min" 
		return nothing
	end

	if observableId == "pEGFR_ID6_MKN1_HM_15min" 
		return nothing
	end

	if observableId == "pEGFR_ID6_MKN1_HM_1min" 
		return nothing
	end

	if observableId == "pEGFR_ID6_MKN1_HM_240min" 
		return nothing
	end

	if observableId == "pEGFR_ID6_MKN1_HM_3min" 
		return nothing
	end

	if observableId == "pEGFR_ID7_HS746T_FM_15min" 
		return nothing
	end

	if observableId == "pEGFR_ID7_HS746T_FM_1min" 
		return nothing
	end

	if observableId == "pEGFR_ID7_HS746T_FM_240min" 
		return nothing
	end

	if observableId == "pEGFR_ID7_HS746T_FM_3min" 
		return nothing
	end

	if observableId == "pEGFR_ID8_HS746T_HM_15min" 
		return nothing
	end

	if observableId == "pEGFR_ID8_HS746T_HM_1min" 
		return nothing
	end

	if observableId == "pEGFR_ID8_HS746T_HM_240min" 
		return nothing
	end

	if observableId == "pEGFR_ID8_HS746T_HM_3min" 
		return nothing
	end

	if observableId == "pEGFR_ID9_HS746T_FM" 
		return nothing
	end

	if observableId == "pEGFR_ID9_HS746T_HM_30EGF" 
		return nothing
	end

	if observableId == "pEGFR_ID9_HS746T_HM_5EGF" 
		return nothing
	end

	if observableId == "pEGFR_ID9b_MKN1_HM_1EGF" 
		return nothing
	end

	if observableId == "pEGFR_ID9b_MKN1_HM_5EGF" 
		return nothing
	end

	if observableId == "pMAPK_ID10_HS746T_FM_EGF_CET" 
		return nothing
	end

	if observableId == "pMAPK_ID10_HS746T_FM_EGF" 
		return nothing
	end

	if observableId == "pMAPK_ID12_MKN1_FM_HM" 
		return nothing
	end

	if observableId == "pMAPK_ID2_HS746T_HM" 
		return nothing
	end

	if observableId == "pMAPK_ID2_MKN1_HM" 
		return nothing
	end

	if observableId == "pMAPK_ID5_MKN1_FM_15min" 
		return nothing
	end

	if observableId == "pMAPK_ID5_MKN1_FM_1min" 
		return nothing
	end

	if observableId == "pMAPK_ID5_MKN1_FM_240min" 
		return nothing
	end

	if observableId == "pMAPK_ID5_MKN1_FM_3min" 
		return nothing
	end

	if observableId == "pMAPK_ID7_HS746T_FM_15min" 
		return nothing
	end

	if observableId == "pMAPK_ID7_HS746T_FM_1min" 
		return nothing
	end

	if observableId == "pMAPK_ID7_HS746T_FM_240min" 
		return nothing
	end

	if observableId == "pMAPK_ID7_HS746T_FM_3min" 
		return nothing
	end

	if observableId == "pMAPK_ID8_HS746T_HM_15min" 
		return nothing
	end

	if observableId == "pMAPK_ID8_HS746T_HM_1min" 
		return nothing
	end

	if observableId == "pMAPK_ID8_HS746T_HM_240min" 
		return nothing
	end

	if observableId == "pMAPK_ID8_HS746T_HM_3min" 
		return nothing
	end

	if observableId == "pMAPK_ID9_HS746T_FM" 
		return nothing
	end

	if observableId == "pMAPK_ID9_HS746T_HM_30EGF" 
		return nothing
	end

	if observableId == "pMAPK_ID9_HS746T_HM_5EGF" 
		return nothing
	end

	if observableId == "pMAPK_ID9b_MKN1_HM_1EGF" 
		return nothing
	end

	if observableId == "pMAPK_ID9b_MKN1_HM_5EGF" 
		return nothing
	end

end

