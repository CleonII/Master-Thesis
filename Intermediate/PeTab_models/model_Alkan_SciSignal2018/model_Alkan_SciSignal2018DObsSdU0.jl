function evalDYmodDu(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	Cells_Dead, Cells_Apo1, pATR_G2, yH2AX_G2, Wip1_mRNA_S, pATM_S, p21_G2, Cells_SSBDamage_S, pp53_S, pATM_G2, Cells_Apo2, Cells_Cycle_G2, Cells_Apo_ReOx, yH2AX_S, p21_mRNA_S, Wip1_S, Wip1_mRNA_G2, Cells_Apo4, Cells_DSBDamage_G2, Cells_Apo, pChk2_G2, Wip1_G2, pATR_S, p21_S, Cells_DSBDamage_S, p21_mRNA_G2, Space, pChk2_S, Cells_Cycle_S, Cells_Apo3, pChk1_G2, pDNAPK_S, pChk1_S, pDNAPK_G2, pp53_G2, Cells= u 
	k_dox_apo, p_p53_act_atm, p_chk1_dea_wip1, k_rep_hr, p_p53_act_chk1, k_rep_nhej, k_death, k_dox_kd, init_Cells_Cycle_G2_rel, k_apo_dsb_g2, k_rep_s, p_p21_turn, p_h2ax_act_atm, p_p53_act_chk2, p_atr_act_ssb, p_mrna_exp_inh_dox_kd, p_chk2_dea_wip1, k_cyc_arr_chk2, p_wip1_mrna_turn, p_mrna_exp_inh_dox, Gem_level, k_ssb_to_dsb, k_rep_nhej_sat, Dox_level, k_lyse, k_damage_sn38_ssb, p_chk1_act, p_atr_act_atm, default, k_apo_ssb, k_apo_dsb_s, k_cyc_arr_chk1, p_h2ax_act_atr, k_death_delay, p_h2ax_act_dnapk, k_damage_dox_dsb, p_chk2_act, p_atm_act_dsb, kt, p_dnapk_act, k_damage_gem_ssb, p_p53_act_atr, init_Cells_Cycle_S_rel, k_death_reox, kt_apo, k_ssb_to_dsb_sn38, init_Space, init_Cells, p_dnapk_dea_wip1, SN38_level, p_h2ax_dea, p_p21_mrna_turn, p_wip1_turn = p 
	if observableId == "CellsCasp_count" 
		out[13] = 1
		out[18] = 1
		return nothing
	end

	if observableId == "Cells_count" 
		out[2] = 1
		out[8] = 1
		out[11] = 1
		out[12] = 1
		out[13] = 1
		out[18] = 1
		out[19] = 1
		out[20] = 1
		out[25] = 1
		out[29] = 1
		out[30] = 1
		out[36] = 1
		return nothing
	end

	if observableId == "Wip1_mRNA_fold" 
		observableParameter1_Wip1_mRNA_fold = getObsOrSdParam(obsPar, mapObsParam)
		out[5] = observableParameter1_Wip1_mRNA_fold
		out[17] = observableParameter1_Wip1_mRNA_fold
		return nothing
	end

	if observableId == "p21_mRNA_fold" 
		observableParameter1_p21_mRNA_fold = getObsOrSdParam(obsPar, mapObsParam)
		out[15] = observableParameter1_p21_mRNA_fold
		out[26] = observableParameter1_p21_mRNA_fold
		return nothing
	end

	if observableId == "pATM_au" 
		observableParameter1_pATM_au = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = observableParameter1_pATM_au
		out[10] = observableParameter1_pATM_au
		return nothing
	end

	if observableId == "pChk1_au" 
		observableParameter1_pChk1_au = getObsOrSdParam(obsPar, mapObsParam)
		out[31] = observableParameter1_pChk1_au
		out[33] = observableParameter1_pChk1_au
		return nothing
	end

	if observableId == "pChk2_au" 
		observableParameter1_pChk2_au = getObsOrSdParam(obsPar, mapObsParam)
		out[21] = observableParameter1_pChk2_au
		out[28] = observableParameter1_pChk2_au
		return nothing
	end

	if observableId == "pDNAPK_au" 
		observableParameter1_pDNAPK_au = getObsOrSdParam(obsPar, mapObsParam)
		out[32] = observableParameter1_pDNAPK_au
		out[34] = observableParameter1_pDNAPK_au
		return nothing
	end

	if observableId == "pp53_au" 
		observableParameter1_pp53_au = getObsOrSdParam(obsPar, mapObsParam)
		out[9] = observableParameter1_pp53_au
		out[35] = observableParameter1_pp53_au
		return nothing
	end

	if observableId == "tp21_au" 
		observableParameter1_tp21_au = getObsOrSdParam(obsPar, mapObsParam)
		out[7] = observableParameter1_tp21_au
		out[24] = observableParameter1_tp21_au
		return nothing
	end

	if observableId == "tp53_au" 
		observableParameter1_tp53_au = getObsOrSdParam(obsPar, mapObsParam)
		out[9] = observableParameter1_tp53_au
		out[35] = observableParameter1_tp53_au
		return nothing
	end

	if observableId == "yH2AX_au" 
		observableParameter1_yH2AX_au = getObsOrSdParam(obsPar, mapObsParam)
		out[4] = observableParameter1_yH2AX_au
		out[14] = observableParameter1_yH2AX_au
		return nothing
	end

end

function evalDYmodDp(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	Cells_Dead, Cells_Apo1, pATR_G2, yH2AX_G2, Wip1_mRNA_S, pATM_S, p21_G2, Cells_SSBDamage_S, pp53_S, pATM_G2, Cells_Apo2, Cells_Cycle_G2, Cells_Apo_ReOx, yH2AX_S, p21_mRNA_S, Wip1_S, Wip1_mRNA_G2, Cells_Apo4, Cells_DSBDamage_G2, Cells_Apo, pChk2_G2, Wip1_G2, pATR_S, p21_S, Cells_DSBDamage_S, p21_mRNA_G2, Space, pChk2_S, Cells_Cycle_S, Cells_Apo3, pChk1_G2, pDNAPK_S, pChk1_S, pDNAPK_G2, pp53_G2, Cells= u 
	k_dox_apo, p_p53_act_atm, p_chk1_dea_wip1, k_rep_hr, p_p53_act_chk1, k_rep_nhej, k_death, k_dox_kd, init_Cells_Cycle_G2_rel, k_apo_dsb_g2, k_rep_s, p_p21_turn, p_h2ax_act_atm, p_p53_act_chk2, p_atr_act_ssb, p_mrna_exp_inh_dox_kd, p_chk2_dea_wip1, k_cyc_arr_chk2, p_wip1_mrna_turn, p_mrna_exp_inh_dox, Gem_level, k_ssb_to_dsb, k_rep_nhej_sat, Dox_level, k_lyse, k_damage_sn38_ssb, p_chk1_act, p_atr_act_atm, default, k_apo_ssb, k_apo_dsb_s, k_cyc_arr_chk1, p_h2ax_act_atr, k_death_delay, p_h2ax_act_dnapk, k_damage_dox_dsb, p_chk2_act, p_atm_act_dsb, kt, p_dnapk_act, k_damage_gem_ssb, p_p53_act_atr, init_Cells_Cycle_S_rel, k_death_reox, kt_apo, k_ssb_to_dsb_sn38, init_Space, init_Cells, p_dnapk_dea_wip1, SN38_level, p_h2ax_dea, p_p21_mrna_turn, p_wip1_turn = p 
	if observableId == "CellsCasp_count" 
		return nothing
	end

	if observableId == "Cells_count" 
		return nothing
	end

	if observableId == "Wip1_mRNA_fold" 
		return nothing
	end

	if observableId == "p21_mRNA_fold" 
		return nothing
	end

	if observableId == "pATM_au" 
		return nothing
	end

	if observableId == "pChk1_au" 
		return nothing
	end

	if observableId == "pChk2_au" 
		return nothing
	end

	if observableId == "pDNAPK_au" 
		return nothing
	end

	if observableId == "pp53_au" 
		return nothing
	end

	if observableId == "tp21_au" 
		return nothing
	end

	if observableId == "tp53_au" 
		return nothing
	end

	if observableId == "yH2AX_au" 
		return nothing
	end

end

function evalDSdDu!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	Cells_Dead, Cells_Apo1, pATR_G2, yH2AX_G2, Wip1_mRNA_S, pATM_S, p21_G2, Cells_SSBDamage_S, pp53_S, pATM_G2, Cells_Apo2, Cells_Cycle_G2, Cells_Apo_ReOx, yH2AX_S, p21_mRNA_S, Wip1_S, Wip1_mRNA_G2, Cells_Apo4, Cells_DSBDamage_G2, Cells_Apo, pChk2_G2, Wip1_G2, pATR_S, p21_S, Cells_DSBDamage_S, p21_mRNA_G2, Space, pChk2_S, Cells_Cycle_S, Cells_Apo3, pChk1_G2, pDNAPK_S, pChk1_S, pDNAPK_G2, pp53_G2, Cells= u 
	k_dox_apo, p_p53_act_atm, p_chk1_dea_wip1, k_rep_hr, p_p53_act_chk1, k_rep_nhej, k_death, k_dox_kd, init_Cells_Cycle_G2_rel, k_apo_dsb_g2, k_rep_s, p_p21_turn, p_h2ax_act_atm, p_p53_act_chk2, p_atr_act_ssb, p_mrna_exp_inh_dox_kd, p_chk2_dea_wip1, k_cyc_arr_chk2, p_wip1_mrna_turn, p_mrna_exp_inh_dox, Gem_level, k_ssb_to_dsb, k_rep_nhej_sat, Dox_level, k_lyse, k_damage_sn38_ssb, p_chk1_act, p_atr_act_atm, default, k_apo_ssb, k_apo_dsb_s, k_cyc_arr_chk1, p_h2ax_act_atr, k_death_delay, p_h2ax_act_dnapk, k_damage_dox_dsb, p_chk2_act, p_atm_act_dsb, kt, p_dnapk_act, k_damage_gem_ssb, p_p53_act_atr, init_Cells_Cycle_S_rel, k_death_reox, kt_apo, k_ssb_to_dsb_sn38, init_Space, init_Cells, p_dnapk_dea_wip1, SN38_level, p_h2ax_dea, p_p21_mrna_turn, p_wip1_turn = p 
	if observableId == "CellsCasp_count" 
		return nothing
	end

	if observableId == "Cells_count" 
		return nothing
	end

	if observableId == "Wip1_mRNA_fold" 
		return nothing
	end

	if observableId == "p21_mRNA_fold" 
		return nothing
	end

	if observableId == "pATM_au" 
		return nothing
	end

	if observableId == "pChk1_au" 
		return nothing
	end

	if observableId == "pChk2_au" 
		return nothing
	end

	if observableId == "pDNAPK_au" 
		return nothing
	end

	if observableId == "pp53_au" 
		return nothing
	end

	if observableId == "tp21_au" 
		return nothing
	end

	if observableId == "tp53_au" 
		return nothing
	end

	if observableId == "yH2AX_au" 
		return nothing
	end

end

function evalDSdDp!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	Cells_Dead, Cells_Apo1, pATR_G2, yH2AX_G2, Wip1_mRNA_S, pATM_S, p21_G2, Cells_SSBDamage_S, pp53_S, pATM_G2, Cells_Apo2, Cells_Cycle_G2, Cells_Apo_ReOx, yH2AX_S, p21_mRNA_S, Wip1_S, Wip1_mRNA_G2, Cells_Apo4, Cells_DSBDamage_G2, Cells_Apo, pChk2_G2, Wip1_G2, pATR_S, p21_S, Cells_DSBDamage_S, p21_mRNA_G2, Space, pChk2_S, Cells_Cycle_S, Cells_Apo3, pChk1_G2, pDNAPK_S, pChk1_S, pDNAPK_G2, pp53_G2, Cells= u 
	k_dox_apo, p_p53_act_atm, p_chk1_dea_wip1, k_rep_hr, p_p53_act_chk1, k_rep_nhej, k_death, k_dox_kd, init_Cells_Cycle_G2_rel, k_apo_dsb_g2, k_rep_s, p_p21_turn, p_h2ax_act_atm, p_p53_act_chk2, p_atr_act_ssb, p_mrna_exp_inh_dox_kd, p_chk2_dea_wip1, k_cyc_arr_chk2, p_wip1_mrna_turn, p_mrna_exp_inh_dox, Gem_level, k_ssb_to_dsb, k_rep_nhej_sat, Dox_level, k_lyse, k_damage_sn38_ssb, p_chk1_act, p_atr_act_atm, default, k_apo_ssb, k_apo_dsb_s, k_cyc_arr_chk1, p_h2ax_act_atr, k_death_delay, p_h2ax_act_dnapk, k_damage_dox_dsb, p_chk2_act, p_atm_act_dsb, kt, p_dnapk_act, k_damage_gem_ssb, p_p53_act_atr, init_Cells_Cycle_S_rel, k_death_reox, kt_apo, k_ssb_to_dsb_sn38, init_Space, init_Cells, p_dnapk_dea_wip1, SN38_level, p_h2ax_dea, p_p21_mrna_turn, p_wip1_turn = p 
	if observableId == "CellsCasp_count" 
		return nothing
	end

	if observableId == "Cells_count" 
		return nothing
	end

	if observableId == "Wip1_mRNA_fold" 
		return nothing
	end

	if observableId == "p21_mRNA_fold" 
		return nothing
	end

	if observableId == "pATM_au" 
		return nothing
	end

	if observableId == "pChk1_au" 
		return nothing
	end

	if observableId == "pChk2_au" 
		return nothing
	end

	if observableId == "pDNAPK_au" 
		return nothing
	end

	if observableId == "pp53_au" 
		return nothing
	end

	if observableId == "tp21_au" 
		return nothing
	end

	if observableId == "tp53_au" 
		return nothing
	end

	if observableId == "yH2AX_au" 
		return nothing
	end

end

