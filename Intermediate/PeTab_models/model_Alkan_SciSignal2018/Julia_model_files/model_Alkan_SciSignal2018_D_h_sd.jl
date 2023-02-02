#u[1] = Cells_Dead, u[2] = Cells_Apo1, u[3] = pATR_G2, u[4] = yH2AX_G2, u[5] = Wip1_mRNA_S, u[6] = pATM_S, u[7] = p21_G2, u[8] = Cells_SSBDamage_S, u[9] = pp53_S, u[10] = pATM_G2, u[11] = Cells_Apo2, u[12] = Cells_Cycle_G2, u[13] = Cells_Apo_ReOx, u[14] = yH2AX_S, u[15] = p21_mRNA_S, u[16] = Wip1_S, u[17] = Wip1_mRNA_G2, u[18] = Cells_Apo4, u[19] = Cells_DSBDamage_G2, u[20] = Cells_Apo, u[21] = pChk2_G2, u[22] = Wip1_G2, u[23] = pATR_S, u[24] = p21_S, u[25] = Cells_DSBDamage_S, u[26] = p21_mRNA_G2, u[27] = Space, u[28] = pChk2_S, u[29] = Cells_Cycle_S, u[30] = Cells_Apo3, u[31] = pChk1_G2, u[32] = pDNAPK_S, u[33] = pChk1_S, u[34] = pDNAPK_G2, u[35] = pp53_G2, u[36] = Cells
#pODEProblem[1] = k_dox_apo, pODEProblem[2] = p_p53_act_atm, pODEProblem[3] = p_chk1_dea_wip1, pODEProblem[4] = k_rep_hr, pODEProblem[5] = p_p53_act_chk1, pODEProblem[6] = k_rep_nhej, pODEProblem[7] = k_death, pODEProblem[8] = k_dox_kd, pODEProblem[9] = init_Cells_Cycle_G2_rel, pODEProblem[10] = k_apo_dsb_g2, pODEProblem[11] = k_rep_s, pODEProblem[12] = p_p21_turn, pODEProblem[13] = p_h2ax_act_atm, pODEProblem[14] = p_p53_act_chk2, pODEProblem[15] = p_atr_act_ssb, pODEProblem[16] = p_mrna_exp_inh_dox_kd, pODEProblem[17] = p_chk2_dea_wip1, pODEProblem[18] = k_cyc_arr_chk2, pODEProblem[19] = p_wip1_mrna_turn, pODEProblem[20] = p_mrna_exp_inh_dox, pODEProblem[21] = Gem_level, pODEProblem[22] = k_ssb_to_dsb, pODEProblem[23] = k_rep_nhej_sat, pODEProblem[24] = Dox_level, pODEProblem[25] = k_lyse, pODEProblem[26] = k_damage_sn38_ssb, pODEProblem[27] = p_chk1_act, pODEProblem[28] = p_atr_act_atm, pODEProblem[29] = default, pODEProblem[30] = k_apo_ssb, pODEProblem[31] = k_apo_dsb_s, pODEProblem[32] = k_cyc_arr_chk1, pODEProblem[33] = p_h2ax_act_atr, pODEProblem[34] = k_death_delay, pODEProblem[35] = p_h2ax_act_dnapk, pODEProblem[36] = k_damage_dox_dsb, pODEProblem[37] = p_chk2_act, pODEProblem[38] = p_atm_act_dsb, pODEProblem[39] = kt, pODEProblem[40] = p_dnapk_act, pODEProblem[41] = k_damage_gem_ssb, pODEProblem[42] = p_p53_act_atr, pODEProblem[43] = init_Cells_Cycle_S_rel, pODEProblem[44] = k_death_reox, pODEProblem[45] = kt_apo, pODEProblem[46] = k_ssb_to_dsb_sn38, pODEProblem[47] = init_Space, pODEProblem[48] = init_Cells, pODEProblem[49] = p_dnapk_dea_wip1, pODEProblem[50] = SN38_level, pODEProblem[51] = p_h2ax_dea, pODEProblem[52] = p_p21_mrna_turn, pODEProblem[53] = p_wip1_turn
#
function compute_∂h∂u!(u, t::Real, pODEProblem::AbstractVector, θ_observable::AbstractVector,                    
                       θ_nonDynamic::AbstractVector, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :CellsCasp_count 
		out[13] = 1
		out[18] = 1
		return nothing
	end

	if observableId == :Cells_count 
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

	if observableId == :Wip1_mRNA_fold 
		observableParameter1_Wip1_mRNA_fold = getObsOrSdParam(θ_observable, parameterMap)
		out[5] = observableParameter1_Wip1_mRNA_fold
		out[17] = observableParameter1_Wip1_mRNA_fold
		return nothing
	end

	if observableId == :p21_mRNA_fold 
		observableParameter1_p21_mRNA_fold = getObsOrSdParam(θ_observable, parameterMap)
		out[15] = observableParameter1_p21_mRNA_fold
		out[26] = observableParameter1_p21_mRNA_fold
		return nothing
	end

	if observableId == :pATM_au 
		observableParameter1_pATM_au = getObsOrSdParam(θ_observable, parameterMap)
		out[6] = observableParameter1_pATM_au
		out[10] = observableParameter1_pATM_au
		return nothing
	end

	if observableId == :pChk1_au 
		observableParameter1_pChk1_au = getObsOrSdParam(θ_observable, parameterMap)
		out[31] = observableParameter1_pChk1_au
		out[33] = observableParameter1_pChk1_au
		return nothing
	end

	if observableId == :pChk2_au 
		observableParameter1_pChk2_au = getObsOrSdParam(θ_observable, parameterMap)
		out[21] = observableParameter1_pChk2_au
		out[28] = observableParameter1_pChk2_au
		return nothing
	end

	if observableId == :pDNAPK_au 
		observableParameter1_pDNAPK_au = getObsOrSdParam(θ_observable, parameterMap)
		out[32] = observableParameter1_pDNAPK_au
		out[34] = observableParameter1_pDNAPK_au
		return nothing
	end

	if observableId == :pp53_au 
		observableParameter1_pp53_au = getObsOrSdParam(θ_observable, parameterMap)
		out[9] = observableParameter1_pp53_au
		out[35] = observableParameter1_pp53_au
		return nothing
	end

	if observableId == :tp21_au 
		observableParameter1_tp21_au = getObsOrSdParam(θ_observable, parameterMap)
		out[7] = observableParameter1_tp21_au
		out[24] = observableParameter1_tp21_au
		return nothing
	end

	if observableId == :tp53_au 
		observableParameter1_tp53_au = getObsOrSdParam(θ_observable, parameterMap)
		out[9] = observableParameter1_tp53_au
		out[35] = observableParameter1_tp53_au
		return nothing
	end

	if observableId == :yH2AX_au 
		observableParameter1_yH2AX_au = getObsOrSdParam(θ_observable, parameterMap)
		out[4] = observableParameter1_yH2AX_au
		out[14] = observableParameter1_yH2AX_au
		return nothing
	end

end

function compute_∂h∂p!(u, t::Real, pODEProblem::AbstractVector, θ_observable::AbstractVector,
                       θ_nonDynamic::AbstractVector, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :CellsCasp_count 
		return nothing
	end

	if observableId == :Cells_count 
		return nothing
	end

	if observableId == :Wip1_mRNA_fold 
		return nothing
	end

	if observableId == :p21_mRNA_fold 
		return nothing
	end

	if observableId == :pATM_au 
		return nothing
	end

	if observableId == :pChk1_au 
		return nothing
	end

	if observableId == :pChk2_au 
		return nothing
	end

	if observableId == :pDNAPK_au 
		return nothing
	end

	if observableId == :pp53_au 
		return nothing
	end

	if observableId == :tp21_au 
		return nothing
	end

	if observableId == :tp53_au 
		return nothing
	end

	if observableId == :yH2AX_au 
		return nothing
	end

end

function compute_∂σ∂σu!(u, t::Real, θ_sd::AbstractVector, pODEProblem::AbstractVector, θ_nonDynamic::AbstractVector, 
                        parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :CellsCasp_count 
		return nothing
	end

	if observableId == :Cells_count 
		return nothing
	end

	if observableId == :Wip1_mRNA_fold 
		return nothing
	end

	if observableId == :p21_mRNA_fold 
		return nothing
	end

	if observableId == :pATM_au 
		return nothing
	end

	if observableId == :pChk1_au 
		return nothing
	end

	if observableId == :pChk2_au 
		return nothing
	end

	if observableId == :pDNAPK_au 
		return nothing
	end

	if observableId == :pp53_au 
		return nothing
	end

	if observableId == :tp21_au 
		return nothing
	end

	if observableId == :tp53_au 
		return nothing
	end

	if observableId == :yH2AX_au 
		return nothing
	end

end

function compute_∂σ∂σp!(u, t::Real, θ_sd::AbstractVector, pODEProblem::AbstractVector, θ_nonDynamic::AbstractVector, 
                        parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :CellsCasp_count 
		return nothing
	end

	if observableId == :Cells_count 
		return nothing
	end

	if observableId == :Wip1_mRNA_fold 
		return nothing
	end

	if observableId == :p21_mRNA_fold 
		return nothing
	end

	if observableId == :pATM_au 
		return nothing
	end

	if observableId == :pChk1_au 
		return nothing
	end

	if observableId == :pChk2_au 
		return nothing
	end

	if observableId == :pDNAPK_au 
		return nothing
	end

	if observableId == :pp53_au 
		return nothing
	end

	if observableId == :tp21_au 
		return nothing
	end

	if observableId == :tp53_au 
		return nothing
	end

	if observableId == :yH2AX_au 
		return nothing
	end

end

