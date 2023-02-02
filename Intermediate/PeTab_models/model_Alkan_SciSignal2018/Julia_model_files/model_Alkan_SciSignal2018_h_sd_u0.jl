#u[1] = Cells_Dead, u[2] = Cells_Apo1, u[3] = pATR_G2, u[4] = yH2AX_G2, u[5] = Wip1_mRNA_S, u[6] = pATM_S, u[7] = p21_G2, u[8] = Cells_SSBDamage_S, u[9] = pp53_S, u[10] = pATM_G2, u[11] = Cells_Apo2, u[12] = Cells_Cycle_G2, u[13] = Cells_Apo_ReOx, u[14] = yH2AX_S, u[15] = p21_mRNA_S, u[16] = Wip1_S, u[17] = Wip1_mRNA_G2, u[18] = Cells_Apo4, u[19] = Cells_DSBDamage_G2, u[20] = Cells_Apo, u[21] = pChk2_G2, u[22] = Wip1_G2, u[23] = pATR_S, u[24] = p21_S, u[25] = Cells_DSBDamage_S, u[26] = p21_mRNA_G2, u[27] = Space, u[28] = pChk2_S, u[29] = Cells_Cycle_S, u[30] = Cells_Apo3, u[31] = pChk1_G2, u[32] = pDNAPK_S, u[33] = pChk1_S, u[34] = pDNAPK_G2, u[35] = pp53_G2, u[36] = Cells
#θ_dynamicNames[1] = k_damage_dox_dsb, θ_dynamicNames[2] = k_damage_gem_ssb, θ_dynamicNames[3] = k_damage_sn38_ssb, θ_dynamicNames[4] = k_death_reox, θ_dynamicNames[5] = k_dox_apo, θ_dynamicNames[6] = k_rep_hr, θ_dynamicNames[7] = k_rep_nhej, θ_dynamicNames[8] = k_rep_nhej_sat, θ_dynamicNames[9] = k_rep_s, θ_dynamicNames[10] = k_ssb_to_dsb, θ_dynamicNames[11] = k_ssb_to_dsb_sn38, θ_dynamicNames[12] = kt, θ_dynamicNames[13] = p_atm_act_dsb, θ_dynamicNames[14] = p_atr_act_atm, θ_dynamicNames[15] = p_atr_act_ssb, θ_dynamicNames[16] = p_chk1_act, θ_dynamicNames[17] = p_chk1_dea_wip1, θ_dynamicNames[18] = p_chk2_act, θ_dynamicNames[19] = p_chk2_dea_wip1, θ_dynamicNames[20] = p_dnapk_act, θ_dynamicNames[21] = p_dnapk_dea_wip1, θ_dynamicNames[22] = p_h2ax_act_atm, θ_dynamicNames[23] = p_h2ax_act_atr, θ_dynamicNames[24] = p_h2ax_act_dnapk, θ_dynamicNames[25] = p_h2ax_dea, θ_dynamicNames[26] = p_mrna_exp_inh_dox, θ_dynamicNames[27] = p_p21_mrna_turn, θ_dynamicNames[28] = p_p21_turn, θ_dynamicNames[29] = p_p53_act_atm, θ_dynamicNames[30] = p_p53_act_atr, θ_dynamicNames[31] = p_p53_act_chk1, θ_dynamicNames[32] = p_p53_act_chk2, θ_dynamicNames[33] = p_wip1_mrna_turn, θ_dynamicNames[34] = p_wip1_turn
##parameterInfo.nominalValue[1] = init_Cells_C 
#parameterInfo.nominalValue[2] = init_Cells_Cycle_G2_rel_C 
#parameterInfo.nominalValue[3] = init_Cells_Cycle_S_rel_C 
#parameterInfo.nominalValue[4] = init_Space_C 
#parameterInfo.nominalValue[5] = k_apo_dsb_g2_C 
#parameterInfo.nominalValue[6] = k_apo_dsb_s_C 
#parameterInfo.nominalValue[7] = k_apo_ssb_C 
#parameterInfo.nominalValue[8] = k_cyc_arr_chk1_C 
#parameterInfo.nominalValue[9] = k_cyc_arr_chk2_C 
#parameterInfo.nominalValue[13] = k_death_C 
#parameterInfo.nominalValue[14] = k_death_delay_C 
#parameterInfo.nominalValue[17] = k_dox_kd_C 
#parameterInfo.nominalValue[18] = k_lyse_C 
#parameterInfo.nominalValue[26] = kt_apo_C 
#parameterInfo.nominalValue[41] = p_mrna_exp_inh_dox_kd_C 


function compute_h(u::AbstractVector, t::Real, θ_dynamic::AbstractVector, θ_observable::AbstractVector,
                   θ_nonDynamic::AbstractVector, parameterInfo::ParametersInfo, observableId::Symbol, 
                      parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :CellsCasp_count 
		return u[18] + u[13] 
	end

	if observableId == :Cells_count 
		return u[36] + u[29] + u[12] + u[8] + u[25] + u[19] + u[20] + u[2] + u[11] + u[30] + u[18] + u[13] 
	end

	if observableId == :Wip1_mRNA_fold 
		observableParameter1_Wip1_mRNA_fold = getObsOrSdParam(θ_observable, parameterMap)
		return observableParameter1_Wip1_mRNA_fold * ( u[17] + u[5] ) + 1 
	end

	if observableId == :p21_mRNA_fold 
		observableParameter1_p21_mRNA_fold = getObsOrSdParam(θ_observable, parameterMap)
		return observableParameter1_p21_mRNA_fold * ( u[15] + u[26] ) + 1 
	end

	if observableId == :pATM_au 
		observableParameter1_pATM_au = getObsOrSdParam(θ_observable, parameterMap)
		return observableParameter1_pATM_au * ( u[10] + u[6] ) 
	end

	if observableId == :pChk1_au 
		observableParameter1_pChk1_au = getObsOrSdParam(θ_observable, parameterMap)
		return observableParameter1_pChk1_au * ( u[31] + u[33] ) 
	end

	if observableId == :pChk2_au 
		observableParameter1_pChk2_au = getObsOrSdParam(θ_observable, parameterMap)
		return observableParameter1_pChk2_au * ( u[28] + u[21] ) 
	end

	if observableId == :pDNAPK_au 
		observableParameter1_pDNAPK_au = getObsOrSdParam(θ_observable, parameterMap)
		return observableParameter1_pDNAPK_au * ( u[34] + u[32] ) 
	end

	if observableId == :pp53_au 
		observableParameter1_pp53_au = getObsOrSdParam(θ_observable, parameterMap)
		return observableParameter1_pp53_au * ( u[9] + u[35] ) 
	end

	if observableId == :tp21_au 
		observableParameter1_tp21_au = getObsOrSdParam(θ_observable, parameterMap)
		return observableParameter1_tp21_au * ( u[24] + u[7] ) 
	end

	if observableId == :tp53_au 
		observableParameter1_tp53_au = getObsOrSdParam(θ_observable, parameterMap)
		return observableParameter1_tp53_au * ( u[9] + u[35] ) 
	end

	if observableId == :yH2AX_au 
		observableParameter1_yH2AX_au = getObsOrSdParam(θ_observable, parameterMap)
		return observableParameter1_yH2AX_au * ( u[14] + u[4] ) 
	end

end

function compute_u0!(u0::AbstractVector, pODEProblem::AbstractVector) 

	#pODEProblem[1] = k_dox_apo, pODEProblem[2] = p_p53_act_atm, pODEProblem[3] = p_chk1_dea_wip1, pODEProblem[4] = k_rep_hr, pODEProblem[5] = p_p53_act_chk1, pODEProblem[6] = k_rep_nhej, pODEProblem[7] = k_death, pODEProblem[8] = k_dox_kd, pODEProblem[9] = init_Cells_Cycle_G2_rel, pODEProblem[10] = k_apo_dsb_g2, pODEProblem[11] = k_rep_s, pODEProblem[12] = p_p21_turn, pODEProblem[13] = p_h2ax_act_atm, pODEProblem[14] = p_p53_act_chk2, pODEProblem[15] = p_atr_act_ssb, pODEProblem[16] = p_mrna_exp_inh_dox_kd, pODEProblem[17] = p_chk2_dea_wip1, pODEProblem[18] = k_cyc_arr_chk2, pODEProblem[19] = p_wip1_mrna_turn, pODEProblem[20] = p_mrna_exp_inh_dox, pODEProblem[21] = Gem_level, pODEProblem[22] = k_ssb_to_dsb, pODEProblem[23] = k_rep_nhej_sat, pODEProblem[24] = Dox_level, pODEProblem[25] = k_lyse, pODEProblem[26] = k_damage_sn38_ssb, pODEProblem[27] = p_chk1_act, pODEProblem[28] = p_atr_act_atm, pODEProblem[29] = default, pODEProblem[30] = k_apo_ssb, pODEProblem[31] = k_apo_dsb_s, pODEProblem[32] = k_cyc_arr_chk1, pODEProblem[33] = p_h2ax_act_atr, pODEProblem[34] = k_death_delay, pODEProblem[35] = p_h2ax_act_dnapk, pODEProblem[36] = k_damage_dox_dsb, pODEProblem[37] = p_chk2_act, pODEProblem[38] = p_atm_act_dsb, pODEProblem[39] = kt, pODEProblem[40] = p_dnapk_act, pODEProblem[41] = k_damage_gem_ssb, pODEProblem[42] = p_p53_act_atr, pODEProblem[43] = init_Cells_Cycle_S_rel, pODEProblem[44] = k_death_reox, pODEProblem[45] = kt_apo, pODEProblem[46] = k_ssb_to_dsb_sn38, pODEProblem[47] = init_Space, pODEProblem[48] = init_Cells, pODEProblem[49] = p_dnapk_dea_wip1, pODEProblem[50] = SN38_level, pODEProblem[51] = p_h2ax_dea, pODEProblem[52] = p_p21_mrna_turn, pODEProblem[53] = p_wip1_turn

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
	Cells_Cycle_G2 = - pODEProblem[48] * pODEProblem[9] * ( pODEProblem[43] - 1 ) 
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
	Space = pODEProblem[47] 
	pChk2_S = 0.0 
	Cells_Cycle_S = pODEProblem[48] * pODEProblem[43] 
	Cells_Apo3 = 0.0 
	pChk1_G2 = 0.0 
	pDNAPK_S = 0.0 
	pChk1_S = 0.0 
	pDNAPK_G2 = 0.0 
	pp53_G2 = 0.0 
	Cells = pODEProblem[48] * ( pODEProblem[9] - 1 ) * ( pODEProblem[43] - 1 ) 

	u0 .= Cells_Dead, Cells_Apo1, pATR_G2, yH2AX_G2, Wip1_mRNA_S, pATM_S, p21_G2, Cells_SSBDamage_S, pp53_S, pATM_G2, Cells_Apo2, Cells_Cycle_G2, Cells_Apo_ReOx, yH2AX_S, p21_mRNA_S, Wip1_S, Wip1_mRNA_G2, Cells_Apo4, Cells_DSBDamage_G2, Cells_Apo, pChk2_G2, Wip1_G2, pATR_S, p21_S, Cells_DSBDamage_S, p21_mRNA_G2, Space, pChk2_S, Cells_Cycle_S, Cells_Apo3, pChk1_G2, pDNAPK_S, pChk1_S, pDNAPK_G2, pp53_G2, Cells
end

function compute_u0(pODEProblem::AbstractVector)::AbstractVector 

	#pODEProblem[1] = k_dox_apo, pODEProblem[2] = p_p53_act_atm, pODEProblem[3] = p_chk1_dea_wip1, pODEProblem[4] = k_rep_hr, pODEProblem[5] = p_p53_act_chk1, pODEProblem[6] = k_rep_nhej, pODEProblem[7] = k_death, pODEProblem[8] = k_dox_kd, pODEProblem[9] = init_Cells_Cycle_G2_rel, pODEProblem[10] = k_apo_dsb_g2, pODEProblem[11] = k_rep_s, pODEProblem[12] = p_p21_turn, pODEProblem[13] = p_h2ax_act_atm, pODEProblem[14] = p_p53_act_chk2, pODEProblem[15] = p_atr_act_ssb, pODEProblem[16] = p_mrna_exp_inh_dox_kd, pODEProblem[17] = p_chk2_dea_wip1, pODEProblem[18] = k_cyc_arr_chk2, pODEProblem[19] = p_wip1_mrna_turn, pODEProblem[20] = p_mrna_exp_inh_dox, pODEProblem[21] = Gem_level, pODEProblem[22] = k_ssb_to_dsb, pODEProblem[23] = k_rep_nhej_sat, pODEProblem[24] = Dox_level, pODEProblem[25] = k_lyse, pODEProblem[26] = k_damage_sn38_ssb, pODEProblem[27] = p_chk1_act, pODEProblem[28] = p_atr_act_atm, pODEProblem[29] = default, pODEProblem[30] = k_apo_ssb, pODEProblem[31] = k_apo_dsb_s, pODEProblem[32] = k_cyc_arr_chk1, pODEProblem[33] = p_h2ax_act_atr, pODEProblem[34] = k_death_delay, pODEProblem[35] = p_h2ax_act_dnapk, pODEProblem[36] = k_damage_dox_dsb, pODEProblem[37] = p_chk2_act, pODEProblem[38] = p_atm_act_dsb, pODEProblem[39] = kt, pODEProblem[40] = p_dnapk_act, pODEProblem[41] = k_damage_gem_ssb, pODEProblem[42] = p_p53_act_atr, pODEProblem[43] = init_Cells_Cycle_S_rel, pODEProblem[44] = k_death_reox, pODEProblem[45] = kt_apo, pODEProblem[46] = k_ssb_to_dsb_sn38, pODEProblem[47] = init_Space, pODEProblem[48] = init_Cells, pODEProblem[49] = p_dnapk_dea_wip1, pODEProblem[50] = SN38_level, pODEProblem[51] = p_h2ax_dea, pODEProblem[52] = p_p21_mrna_turn, pODEProblem[53] = p_wip1_turn

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
	Cells_Cycle_G2 = - pODEProblem[48] * pODEProblem[9] * ( pODEProblem[43] - 1 ) 
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
	Space = pODEProblem[47] 
	pChk2_S = 0.0 
	Cells_Cycle_S = pODEProblem[48] * pODEProblem[43] 
	Cells_Apo3 = 0.0 
	pChk1_G2 = 0.0 
	pDNAPK_S = 0.0 
	pChk1_S = 0.0 
	pDNAPK_G2 = 0.0 
	pp53_G2 = 0.0 
	Cells = pODEProblem[48] * ( pODEProblem[9] - 1 ) * ( pODEProblem[43] - 1 ) 

	 return [Cells_Dead, Cells_Apo1, pATR_G2, yH2AX_G2, Wip1_mRNA_S, pATM_S, p21_G2, Cells_SSBDamage_S, pp53_S, pATM_G2, Cells_Apo2, Cells_Cycle_G2, Cells_Apo_ReOx, yH2AX_S, p21_mRNA_S, Wip1_S, Wip1_mRNA_G2, Cells_Apo4, Cells_DSBDamage_G2, Cells_Apo, pChk2_G2, Wip1_G2, pATR_S, p21_S, Cells_DSBDamage_S, p21_mRNA_G2, Space, pChk2_S, Cells_Cycle_S, Cells_Apo3, pChk1_G2, pDNAPK_S, pChk1_S, pDNAPK_G2, pp53_G2, Cells]
end

function compute_σ(u::AbstractVector, t::Real, θ_sd::AbstractVector, θ_dynamic::AbstractVector, θ_nonDynamic::AbstractVector, 
                   parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :CellsCasp_count 
		noiseParameter1_CellsCasp_count = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_CellsCasp_count 
	end

	if observableId == :Cells_count 
		noiseParameter1_Cells_count = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_Cells_count 
	end

	if observableId == :Wip1_mRNA_fold 
		noiseParameter1_Wip1_mRNA_fold = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_Wip1_mRNA_fold 
	end

	if observableId == :p21_mRNA_fold 
		noiseParameter1_p21_mRNA_fold = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_p21_mRNA_fold 
	end

	if observableId == :pATM_au 
		noiseParameter1_pATM_au = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_pATM_au 
	end

	if observableId == :pChk1_au 
		noiseParameter1_pChk1_au = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_pChk1_au 
	end

	if observableId == :pChk2_au 
		noiseParameter1_pChk2_au = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_pChk2_au 
	end

	if observableId == :pDNAPK_au 
		noiseParameter1_pDNAPK_au = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_pDNAPK_au 
	end

	if observableId == :pp53_au 
		noiseParameter1_pp53_au = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_pp53_au 
	end

	if observableId == :tp21_au 
		noiseParameter1_tp21_au = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_tp21_au 
	end

	if observableId == :tp53_au 
		noiseParameter1_tp53_au = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_tp53_au 
	end

	if observableId == :yH2AX_au 
		noiseParameter1_yH2AX_au = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_yH2AX_au 
	end

end