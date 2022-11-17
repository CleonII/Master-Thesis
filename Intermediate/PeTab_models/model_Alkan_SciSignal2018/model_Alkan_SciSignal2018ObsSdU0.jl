function evalYmod(u, t, dynPar, obsPar, nonDynParam, paramData, observableId, mapObsParam) 

	Cells_Dead, Cells_Apo1, pATR_G2, yH2AX_G2, Wip1_mRNA_S, pATM_S, p21_G2, Cells_SSBDamage_S, pp53_S, pATM_G2, Cells_Apo2, Cells_Cycle_G2, Cells_Apo_ReOx, yH2AX_S, p21_mRNA_S, Wip1_S, Wip1_mRNA_G2, Cells_Apo4, Cells_DSBDamage_G2, Cells_Apo, pChk2_G2, Wip1_G2, pATR_S, p21_S, Cells_DSBDamage_S, p21_mRNA_G2, Space, pChk2_S, Cells_Cycle_S, Cells_Apo3, pChk1_G2, pDNAPK_S, pChk1_S, pDNAPK_G2, pp53_G2, Cells= u 
	k_damage_dox_dsb, k_damage_gem_ssb, k_damage_sn38_ssb, k_death_reox, k_dox_apo, k_rep_hr, k_rep_nhej, k_rep_nhej_sat, k_rep_s, k_ssb_to_dsb, k_ssb_to_dsb_sn38, kt, p_atm_act_dsb, p_atr_act_atm, p_atr_act_ssb, p_chk1_act, p_chk1_dea_wip1, p_chk2_act, p_chk2_dea_wip1, p_dnapk_act, p_dnapk_dea_wip1, p_h2ax_act_atm, p_h2ax_act_atr, p_h2ax_act_dnapk, p_h2ax_dea, p_mrna_exp_inh_dox, p_p21_mrna_turn, p_p21_turn, p_p53_act_atm, p_p53_act_atr, p_p53_act_chk1, p_p53_act_chk2, p_wip1_mrna_turn, p_wip1_turn = dynPar 
	init_Cells_C = paramData.paramVal[1] 
	init_Cells_Cycle_G2_rel_C = paramData.paramVal[2] 
	init_Cells_Cycle_S_rel_C = paramData.paramVal[3] 
	init_Space_C = paramData.paramVal[4] 
	k_apo_dsb_g2_C = paramData.paramVal[5] 
	k_apo_dsb_s_C = paramData.paramVal[6] 
	k_apo_ssb_C = paramData.paramVal[7] 
	k_cyc_arr_chk1_C = paramData.paramVal[8] 
	k_cyc_arr_chk2_C = paramData.paramVal[9] 
	k_death_C = paramData.paramVal[13] 
	k_death_delay_C = paramData.paramVal[14] 
	k_dox_kd_C = paramData.paramVal[17] 
	k_lyse_C = paramData.paramVal[18] 
	kt_apo_C = paramData.paramVal[26] 
	p_mrna_exp_inh_dox_kd_C = paramData.paramVal[41] 

	if observableId == "CellsCasp_count" 
		return Cells_Apo4 + Cells_Apo_ReOx 
	end

	if observableId == "Cells_count" 
		return Cells + Cells_Cycle_S + Cells_Cycle_G2 + Cells_SSBDamage_S + Cells_DSBDamage_S + Cells_DSBDamage_G2 + Cells_Apo + Cells_Apo1 + Cells_Apo2 + Cells_Apo3 + Cells_Apo4 + Cells_Apo_ReOx 
	end

	if observableId == "Wip1_mRNA_fold" 
		observableParameter1_Wip1_mRNA_fold = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_Wip1_mRNA_fold * ( Wip1_mRNA_G2 + Wip1_mRNA_S ) + 1 
	end

	if observableId == "p21_mRNA_fold" 
		observableParameter1_p21_mRNA_fold = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_p21_mRNA_fold * ( p21_mRNA_S + p21_mRNA_G2 ) + 1 
	end

	if observableId == "pATM_au" 
		observableParameter1_pATM_au = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pATM_au * ( pATM_G2 + pATM_S ) 
	end

	if observableId == "pChk1_au" 
		observableParameter1_pChk1_au = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pChk1_au * ( pChk1_G2 + pChk1_S ) 
	end

	if observableId == "pChk2_au" 
		observableParameter1_pChk2_au = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pChk2_au * ( pChk2_S + pChk2_G2 ) 
	end

	if observableId == "pDNAPK_au" 
		observableParameter1_pDNAPK_au = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pDNAPK_au * ( pDNAPK_G2 + pDNAPK_S ) 
	end

	if observableId == "pp53_au" 
		observableParameter1_pp53_au = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pp53_au * ( pp53_S + pp53_G2 ) 
	end

	if observableId == "tp21_au" 
		observableParameter1_tp21_au = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_tp21_au * ( p21_S + p21_G2 ) 
	end

	if observableId == "tp53_au" 
		observableParameter1_tp53_au = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_tp53_au * ( pp53_S + pp53_G2 ) 
	end

	if observableId == "yH2AX_au" 
		observableParameter1_yH2AX_au = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_yH2AX_au * ( yH2AX_S + yH2AX_G2 ) 
	end

end

function evalU0!(u0Vec, paramVec) 

	k_dox_apo, p_p53_act_atm, p_chk1_dea_wip1, k_rep_hr, p_p53_act_chk1, k_rep_nhej, k_death, k_dox_kd, init_Cells_Cycle_G2_rel, k_apo_dsb_g2, k_rep_s, p_p21_turn, p_h2ax_act_atm, p_p53_act_chk2, p_atr_act_ssb, p_mrna_exp_inh_dox_kd, p_chk2_dea_wip1, k_cyc_arr_chk2, p_wip1_mrna_turn, p_mrna_exp_inh_dox, Gem_level, k_ssb_to_dsb, k_rep_nhej_sat, Dox_level, k_lyse, k_damage_sn38_ssb, p_chk1_act, p_atr_act_atm, default, k_apo_ssb, k_apo_dsb_s, k_cyc_arr_chk1, p_h2ax_act_atr, k_death_delay, p_h2ax_act_dnapk, k_damage_dox_dsb, p_chk2_act, p_atm_act_dsb, kt, p_dnapk_act, k_damage_gem_ssb, p_p53_act_atr, init_Cells_Cycle_S_rel, k_death_reox, kt_apo, k_ssb_to_dsb_sn38, init_Space, init_Cells, p_dnapk_dea_wip1, SN38_level, p_h2ax_dea, p_p21_mrna_turn, p_wip1_turn = paramVec 

	Cells_Dead = 0.0 
	Cells_Apo1 = 0.0 
	pATR_G2 = 0.0 
	yH2AX_G2 = 0.0 
	Wip1_mRNA_S = 0.0 
	pATM_S = 0.0 
	p21_G2 = 0.0 
	Cells_SSBDamage_S = 0.0 
	pp53_S = 0.0 
	pATM_G2 = 0.0 
	Cells_Apo2 = 0.0 
	Cells_Cycle_G2 = - init_Cells * init_Cells_Cycle_G2_rel * ( init_Cells_Cycle_S_rel - 1 ) 
	Cells_Apo_ReOx = 0.0 
	yH2AX_S = 0.0 
	p21_mRNA_S = 0.0 
	Wip1_S = 0.0 
	Wip1_mRNA_G2 = 0.0 
	Cells_Apo4 = 0.0 
	Cells_DSBDamage_G2 = 0.0 
	Cells_Apo = 0.0 
	pChk2_G2 = 0.0 
	Wip1_G2 = 0.0 
	pATR_S = 0.0 
	p21_S = 0.0 
	Cells_DSBDamage_S = 0.0 
	p21_mRNA_G2 = 0.0 
	Space = init_Space 
	pChk2_S = 0.0 
	Cells_Cycle_S = init_Cells * init_Cells_Cycle_S_rel 
	Cells_Apo3 = 0.0 
	pChk1_G2 = 0.0 
	pDNAPK_S = 0.0 
	pChk1_S = 0.0 
	pDNAPK_G2 = 0.0 
	pp53_G2 = 0.0 
	Cells = init_Cells * ( init_Cells_Cycle_G2_rel - 1 ) * ( init_Cells_Cycle_S_rel - 1 ) 

	u0Vec .= Cells_Dead, Cells_Apo1, pATR_G2, yH2AX_G2, Wip1_mRNA_S, pATM_S, p21_G2, Cells_SSBDamage_S, pp53_S, pATM_G2, Cells_Apo2, Cells_Cycle_G2, Cells_Apo_ReOx, yH2AX_S, p21_mRNA_S, Wip1_S, Wip1_mRNA_G2, Cells_Apo4, Cells_DSBDamage_G2, Cells_Apo, pChk2_G2, Wip1_G2, pATR_S, p21_S, Cells_DSBDamage_S, p21_mRNA_G2, Space, pChk2_S, Cells_Cycle_S, Cells_Apo3, pChk1_G2, pDNAPK_S, pChk1_S, pDNAPK_G2, pp53_G2, Cells
end

function evalSd!(u, t, sdPar, dynPar, nonDynParam, paramData, observableId, mapSdParam) 

	Cells_Dead, Cells_Apo1, pATR_G2, yH2AX_G2, Wip1_mRNA_S, pATM_S, p21_G2, Cells_SSBDamage_S, pp53_S, pATM_G2, Cells_Apo2, Cells_Cycle_G2, Cells_Apo_ReOx, yH2AX_S, p21_mRNA_S, Wip1_S, Wip1_mRNA_G2, Cells_Apo4, Cells_DSBDamage_G2, Cells_Apo, pChk2_G2, Wip1_G2, pATR_S, p21_S, Cells_DSBDamage_S, p21_mRNA_G2, Space, pChk2_S, Cells_Cycle_S, Cells_Apo3, pChk1_G2, pDNAPK_S, pChk1_S, pDNAPK_G2, pp53_G2, Cells= u 
	k_damage_dox_dsb, k_damage_gem_ssb, k_damage_sn38_ssb, k_death_reox, k_dox_apo, k_rep_hr, k_rep_nhej, k_rep_nhej_sat, k_rep_s, k_ssb_to_dsb, k_ssb_to_dsb_sn38, kt, p_atm_act_dsb, p_atr_act_atm, p_atr_act_ssb, p_chk1_act, p_chk1_dea_wip1, p_chk2_act, p_chk2_dea_wip1, p_dnapk_act, p_dnapk_dea_wip1, p_h2ax_act_atm, p_h2ax_act_atr, p_h2ax_act_dnapk, p_h2ax_dea, p_mrna_exp_inh_dox, p_p21_mrna_turn, p_p21_turn, p_p53_act_atm, p_p53_act_atr, p_p53_act_chk1, p_p53_act_chk2, p_wip1_mrna_turn, p_wip1_turn = dynPar 
	init_Cells_C = paramData.paramVal[1] 
	init_Cells_Cycle_G2_rel_C = paramData.paramVal[2] 
	init_Cells_Cycle_S_rel_C = paramData.paramVal[3] 
	init_Space_C = paramData.paramVal[4] 
	k_apo_dsb_g2_C = paramData.paramVal[5] 
	k_apo_dsb_s_C = paramData.paramVal[6] 
	k_apo_ssb_C = paramData.paramVal[7] 
	k_cyc_arr_chk1_C = paramData.paramVal[8] 
	k_cyc_arr_chk2_C = paramData.paramVal[9] 
	k_death_C = paramData.paramVal[13] 
	k_death_delay_C = paramData.paramVal[14] 
	k_dox_kd_C = paramData.paramVal[17] 
	k_lyse_C = paramData.paramVal[18] 
	kt_apo_C = paramData.paramVal[26] 
	p_mrna_exp_inh_dox_kd_C = paramData.paramVal[41] 

	if observableId == "CellsCasp_count" 
		noiseParameter1_CellsCasp_count = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_CellsCasp_count 
	end

	if observableId == "Cells_count" 
		noiseParameter1_Cells_count = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_Cells_count 
	end

	if observableId == "Wip1_mRNA_fold" 
		noiseParameter1_Wip1_mRNA_fold = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_Wip1_mRNA_fold 
	end

	if observableId == "p21_mRNA_fold" 
		noiseParameter1_p21_mRNA_fold = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_p21_mRNA_fold 
	end

	if observableId == "pATM_au" 
		noiseParameter1_pATM_au = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pATM_au 
	end

	if observableId == "pChk1_au" 
		noiseParameter1_pChk1_au = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pChk1_au 
	end

	if observableId == "pChk2_au" 
		noiseParameter1_pChk2_au = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pChk2_au 
	end

	if observableId == "pDNAPK_au" 
		noiseParameter1_pDNAPK_au = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pDNAPK_au 
	end

	if observableId == "pp53_au" 
		noiseParameter1_pp53_au = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pp53_au 
	end

	if observableId == "tp21_au" 
		noiseParameter1_tp21_au = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_tp21_au 
	end

	if observableId == "tp53_au" 
		noiseParameter1_tp53_au = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_tp53_au 
	end

	if observableId == "yH2AX_au" 
		noiseParameter1_yH2AX_au = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_yH2AX_au 
	end

end