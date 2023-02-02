# Model name: model_Alkan_SciSignal2018
# Number of parameters: 52
# Number of species: 36
function getODEModel_model_Alkan_SciSignal2018()

    ### Define independent and dependent variables
    ModelingToolkit.@variables t Cells_Dead(t) Cells_Apo1(t) pATR_G2(t) yH2AX_G2(t) Wip1_mRNA_S(t) pATM_S(t) p21_G2(t) Cells_SSBDamage_S(t) pp53_S(t) pATM_G2(t) Cells_Apo2(t) Cells_Cycle_G2(t) Cells_Apo_ReOx(t) yH2AX_S(t) p21_mRNA_S(t) Wip1_S(t) Wip1_mRNA_G2(t) Cells_Apo4(t) Cells_DSBDamage_G2(t) Cells_Apo(t) pChk2_G2(t) Wip1_G2(t) pATR_S(t) p21_S(t) Cells_DSBDamage_S(t) p21_mRNA_G2(t) Space(t) pChk2_S(t) Cells_Cycle_S(t) Cells_Apo3(t) pChk1_G2(t) pDNAPK_S(t) pChk1_S(t) pDNAPK_G2(t) pp53_G2(t) Cells(t)

    ### Store dependent variables in array for ODESystem command
    stateArray = [Cells_Dead, Cells_Apo1, pATR_G2, yH2AX_G2, Wip1_mRNA_S, pATM_S, p21_G2, Cells_SSBDamage_S, pp53_S, pATM_G2, Cells_Apo2, Cells_Cycle_G2, Cells_Apo_ReOx, yH2AX_S, p21_mRNA_S, Wip1_S, Wip1_mRNA_G2, Cells_Apo4, Cells_DSBDamage_G2, Cells_Apo, pChk2_G2, Wip1_G2, pATR_S, p21_S, Cells_DSBDamage_S, p21_mRNA_G2, Space, pChk2_S, Cells_Cycle_S, Cells_Apo3, pChk1_G2, pDNAPK_S, pChk1_S, pDNAPK_G2, pp53_G2, Cells]

    ### Define variable parameters

    ### Define potential algebraic variables

    ### Define parameters
    ModelingToolkit.@parameters k_dox_apo p_p53_act_atm p_chk1_dea_wip1 k_rep_hr p_p53_act_chk1 k_rep_nhej k_death k_dox_kd init_Cells_Cycle_G2_rel k_apo_dsb_g2 k_rep_s p_p21_turn p_h2ax_act_atm p_p53_act_chk2 p_atr_act_ssb p_mrna_exp_inh_dox_kd p_chk2_dea_wip1 k_cyc_arr_chk2 p_wip1_mrna_turn p_mrna_exp_inh_dox Gem_level k_ssb_to_dsb k_rep_nhej_sat Dox_level k_lyse k_damage_sn38_ssb p_chk1_act p_atr_act_atm default k_apo_ssb k_apo_dsb_s k_cyc_arr_chk1 p_h2ax_act_atr k_death_delay p_h2ax_act_dnapk k_damage_dox_dsb p_chk2_act p_atm_act_dsb kt p_dnapk_act k_damage_gem_ssb p_p53_act_atr init_Cells_Cycle_S_rel k_death_reox kt_apo k_ssb_to_dsb_sn38 init_Space init_Cells p_dnapk_dea_wip1 SN38_level p_h2ax_dea p_p21_mrna_turn p_wip1_turn

    ### Store parameters in array for ODESystem command
    parameterArray = [k_dox_apo, p_p53_act_atm, p_chk1_dea_wip1, k_rep_hr, p_p53_act_chk1, k_rep_nhej, k_death, k_dox_kd, init_Cells_Cycle_G2_rel, k_apo_dsb_g2, k_rep_s, p_p21_turn, p_h2ax_act_atm, p_p53_act_chk2, p_atr_act_ssb, p_mrna_exp_inh_dox_kd, p_chk2_dea_wip1, k_cyc_arr_chk2, p_wip1_mrna_turn, p_mrna_exp_inh_dox, Gem_level, k_ssb_to_dsb, k_rep_nhej_sat, Dox_level, k_lyse, k_damage_sn38_ssb, p_chk1_act, p_atr_act_atm, default, k_apo_ssb, k_apo_dsb_s, k_cyc_arr_chk1, p_h2ax_act_atr, k_death_delay, p_h2ax_act_dnapk, k_damage_dox_dsb, p_chk2_act, p_atm_act_dsb, kt, p_dnapk_act, k_damage_gem_ssb, p_p53_act_atr, init_Cells_Cycle_S_rel, k_death_reox, kt_apo, k_ssb_to_dsb_sn38, init_Space, init_Cells, p_dnapk_dea_wip1, SN38_level, p_h2ax_dea, p_p21_mrna_turn, p_wip1_turn]

    ### Define an operator for the differentiation w.r.t. time
    D = Differential(t)

    ### Continious events ###

    ### Discrete events ###

    ### Derivatives ###
    eqs = [
    D(Cells_Dead) ~ +1.0 * ( 1 /default ) * (Cells_Apo4 * k_death)+1.0 * ( 1 /default ) * (Cells_Apo_ReOx * k_death_reox)-1.0 * ( 1 /default ) * (Cells_Dead * k_lyse),
    D(Cells_Apo1) ~ +1.0 * ( 1 /default ) * (Cells_Apo * k_death_delay)-1.0 * ( 1 /default ) * (Cells_Apo1 * k_death_delay)-1.0 * ( 1 /default ) * (Cells_Apo1 * k_dox_apo * ((Dox_level)^(6) / ((Dox_level)^(6) + (k_dox_kd)^(6)))),
    D(pATR_G2) ~ +1.0 * ( 1 /default ) * (-pATM_G2 * p_atr_act_atm * (pATR_G2 - 1))-1.0 * ( 1 /default ) * (pATR_G2 / 24),
    D(yH2AX_G2) ~ +1.0 * ( 1 /default ) * (-pATM_G2 * p_h2ax_act_atm * (yH2AX_G2 - 1))+1.0 * ( 1 /default ) * (-pATR_G2 * p_h2ax_act_atr * (yH2AX_G2 - 1))+1.0 * ( 1 /default ) * (-pDNAPK_G2 * p_h2ax_act_dnapk * (yH2AX_G2 - 1))-1.0 * ( 1 /default ) * (p_h2ax_dea * yH2AX_G2),
    D(Wip1_mRNA_S) ~ +1.0 * ( 1 /default ) * (p_wip1_mrna_turn * pp53_S / (p_mrna_exp_inh_dox * ((Dox_level)^(6) / ((Dox_level)^(6) + (p_mrna_exp_inh_dox_kd)^(6))) + 1))-1.0 * ( 1 /default ) * (Wip1_mRNA_S * p_wip1_mrna_turn),
    D(pATM_S) ~ +1.0 * ( 1 /default ) * (-(Cells_DSBDamage_S * p_atm_act_dsb * (pATM_S - 1)) / (Cells + Cells_Apo + Cells_Apo1 + Cells_Apo2 + Cells_Apo3 + Cells_Apo4 + Cells_Cycle_G2 + Cells_Cycle_S + Cells_Apo_ReOx + Cells_DSBDamage_G2 + Cells_DSBDamage_S + Cells_SSBDamage_S))-1.0 * ( 1 /default ) * (pATM_S / 24),
    D(p21_G2) ~ +1.0 * ( 1 /default ) * (p21_mRNA_G2 * p_p21_turn)-1.0 * ( 1 /default ) * (p21_G2 * p_p21_turn),
    D(Cells_SSBDamage_S) ~ -1.0 * ( 1 /default ) * (Cells_SSBDamage_S * k_dox_apo * ((Dox_level)^(6) / ((Dox_level)^(6) + (k_dox_kd)^(6))))+1.0 * ( 1 /default ) * (Gem_level * Cells_Cycle_S * k_damage_gem_ssb)+1.0 * ( 1 /default ) * (SN38_level * Cells_Cycle_S * k_damage_sn38_ssb)-1.0 * ( 1 /default ) * (Cells_SSBDamage_S * (k_ssb_to_dsb + k_ssb_to_dsb_sn38 * ifelse(SN38_level > 0, 1, 0)))-1.0 * ( 1 /default ) * (Cells_SSBDamage_S * k_rep_s * pATR_S)-1.0 * ( 1 /default ) * (Cells_SSBDamage_S * k_apo_ssb * pp53_S)-1.0 * ( 1 /default ) * (Cells_SSBDamage_S * kt_apo / (k_cyc_arr_chk1 * pChk1_S + 1)),
    D(pp53_S) ~ +1.0 * ( 1 /default ) * (-pATR_S * p_p53_act_atr * (pp53_S - 1))+1.0 * ( 1 /default ) * (-pChk1_S * p_p53_act_chk1 * (pp53_S - 1))+1.0 * ( 1 /default ) * (-pATM_S * p_p53_act_atm * (pp53_S - 1))+1.0 * ( 1 /default ) * (-pChk2_S * p_p53_act_chk2 * (pp53_S - 1))-1.0 * ( 1 /default ) * (pp53_S / 24),
    D(pATM_G2) ~ +1.0 * ( 1 /default ) * (-(Cells_DSBDamage_G2 * p_atm_act_dsb * (pATM_G2 - 1)) / (Cells + Cells_Apo + Cells_Apo1 + Cells_Apo2 + Cells_Apo3 + Cells_Apo4 + 2 * Cells_Cycle_G2 + Cells_Apo_ReOx + Cells_DSBDamage_G2 + Cells_DSBDamage_S + Cells_SSBDamage_S))-1.0 * ( 1 /default ) * (pATM_G2 / 24),
    D(Cells_Apo2) ~ +1.0 * ( 1 /default ) * (Cells_Apo1 * k_death_delay)-1.0 * ( 1 /default ) * (Cells_Apo2 * k_death_delay)-1.0 * ( 1 /default ) * (Cells_Apo2 * k_dox_apo * ((Dox_level)^(6) / ((Dox_level)^(6) + (k_dox_kd)^(6)))),
    D(Cells_Cycle_G2) ~ +1.0 * ( 1 /default ) * (Cells_Cycle_S * kt * (init_Cells_Cycle_G2_rel - 1) * (init_Cells_Cycle_S_rel - 1) * (init_Cells * init_Cells_Cycle_S_rel + init_Cells * (init_Cells_Cycle_G2_rel - 1) * (init_Cells_Cycle_S_rel - 1) - 2 * init_Cells * init_Cells_Cycle_G2_rel * (init_Cells_Cycle_S_rel - 1)) / (init_Cells_Cycle_S_rel * (2 * init_Cells * init_Cells_Cycle_S_rel + init_Cells * (init_Cells_Cycle_G2_rel - 1) * (init_Cells_Cycle_S_rel - 1) - 2 * init_Cells * init_Cells_Cycle_G2_rel * (init_Cells_Cycle_S_rel - 1))))-1.0 * ( 1 /default ) * (-(Cells_Cycle_G2 * init_Cells * kt * (init_Cells_Cycle_G2_rel - 1) * (init_Cells_Cycle_S_rel - 1) * (init_Cells * init_Cells_Cycle_S_rel + init_Cells * (init_Cells_Cycle_G2_rel - 1) * (init_Cells_Cycle_S_rel - 1) - init_Cells * init_Cells_Cycle_G2_rel * (init_Cells_Cycle_S_rel - 1))) / ((init_Cells)^(2) * init_Cells_Cycle_G2_rel * (init_Cells_Cycle_G2_rel - 1) * (init_Cells_Cycle_S_rel - 1)^(2) - 2 * (init_Cells)^(2) * (init_Cells_Cycle_G2_rel)^(2) * (init_Cells_Cycle_S_rel - 1)^(2) + 2 * (init_Cells)^(2) * init_Cells_Cycle_G2_rel * init_Cells_Cycle_S_rel * (init_Cells_Cycle_S_rel - 1)))-1.0 * ( 1 /default ) * (Cells_Cycle_G2 * k_dox_apo * ((Dox_level)^(6) / ((Dox_level)^(6) + (k_dox_kd)^(6))))-1.0 * ( 1 /default ) * (Cells_Cycle_G2 * k_damage_dox_dsb * Dox_level)+1.0 * ( 1 /default ) * (Cells_DSBDamage_G2 * k_rep_nhej * pDNAPK_G2 / (Cells_DSBDamage_G2 + k_rep_nhej_sat))+1.0 * ( 1 /default ) * (Cells_DSBDamage_G2 * k_rep_hr * pATM_G2 * pATR_G2),
    D(Cells_Apo_ReOx) ~ +1.0 * ( 1 /default ) * (Cells * k_dox_apo * ((Dox_level)^(6) / ((Dox_level)^(6) + (k_dox_kd)^(6))))+1.0 * ( 1 /default ) * (Cells_Cycle_S * k_dox_apo * ((Dox_level)^(6) / ((Dox_level)^(6) + (k_dox_kd)^(6))))+1.0 * ( 1 /default ) * (Cells_Cycle_G2 * k_dox_apo * ((Dox_level)^(6) / ((Dox_level)^(6) + (k_dox_kd)^(6))))+1.0 * ( 1 /default ) * (Cells_SSBDamage_S * k_dox_apo * ((Dox_level)^(6) / ((Dox_level)^(6) + (k_dox_kd)^(6))))+1.0 * ( 1 /default ) * (Cells_DSBDamage_S * k_dox_apo * ((Dox_level)^(6) / ((Dox_level)^(6) + (k_dox_kd)^(6))))+1.0 * ( 1 /default ) * (Cells_DSBDamage_G2 * k_dox_apo * ((Dox_level)^(6) / ((Dox_level)^(6) + (k_dox_kd)^(6))))+1.0 * ( 1 /default ) * (Cells_Apo * k_dox_apo * ((Dox_level)^(6) / ((Dox_level)^(6) + (k_dox_kd)^(6))))+1.0 * ( 1 /default ) * (Cells_Apo1 * k_dox_apo * ((Dox_level)^(6) / ((Dox_level)^(6) + (k_dox_kd)^(6))))+1.0 * ( 1 /default ) * (Cells_Apo2 * k_dox_apo * ((Dox_level)^(6) / ((Dox_level)^(6) + (k_dox_kd)^(6))))+1.0 * ( 1 /default ) * (Cells_Apo3 * k_dox_apo * ((Dox_level)^(6) / ((Dox_level)^(6) + (k_dox_kd)^(6))))+1.0 * ( 1 /default ) * (Cells_Apo4 * k_dox_apo * ((Dox_level)^(6) / ((Dox_level)^(6) + (k_dox_kd)^(6))))-1.0 * ( 1 /default ) * (Cells_Apo_ReOx * k_death_reox),
    D(yH2AX_S) ~ +1.0 * ( 1 /default ) * (-pATM_S * p_h2ax_act_atm * (yH2AX_S - 1))+1.0 * ( 1 /default ) * (-pATR_S * p_h2ax_act_atr * (yH2AX_S - 1))+1.0 * ( 1 /default ) * (-pDNAPK_S * p_h2ax_act_dnapk * (yH2AX_S - 1))-1.0 * ( 1 /default ) * (p_h2ax_dea * yH2AX_S),
    D(p21_mRNA_S) ~ +1.0 * ( 1 /default ) * (p_p21_mrna_turn * pp53_S / (p_mrna_exp_inh_dox * ((Dox_level)^(6) / ((Dox_level)^(6) + (p_mrna_exp_inh_dox_kd)^(6))) + 1))-1.0 * ( 1 /default ) * (p21_mRNA_S * p_p21_mrna_turn),
    D(Wip1_S) ~ +1.0 * ( 1 /default ) * (Wip1_mRNA_S * p_wip1_turn)-1.0 * ( 1 /default ) * (Wip1_S * p_wip1_turn),
    D(Wip1_mRNA_G2) ~ +1.0 * ( 1 /default ) * (p_wip1_mrna_turn * pp53_G2 / (p_mrna_exp_inh_dox * ((Dox_level)^(6) / ((Dox_level)^(6) + (p_mrna_exp_inh_dox_kd)^(6))) + 1))-1.0 * ( 1 /default ) * (Wip1_mRNA_G2 * p_wip1_mrna_turn),
    D(Cells_Apo4) ~ +1.0 * ( 1 /default ) * (Cells_Apo3 * k_death_delay)-1.0 * ( 1 /default ) * (Cells_Apo4 * k_death)-1.0 * ( 1 /default ) * (Cells_Apo4 * k_dox_apo * ((Dox_level)^(6) / ((Dox_level)^(6) + (k_dox_kd)^(6)))),
    D(Cells_DSBDamage_G2) ~ -1.0 * ( 1 /default ) * (Cells_DSBDamage_G2 * k_dox_apo * ((Dox_level)^(6) / ((Dox_level)^(6) + (k_dox_kd)^(6))))+1.0 * ( 1 /default ) * (Cells_Cycle_G2 * k_damage_dox_dsb * Dox_level)-1.0 * ( 1 /default ) * (Cells_DSBDamage_G2 * k_rep_nhej * pDNAPK_G2 / (Cells_DSBDamage_G2 + k_rep_nhej_sat))-1.0 * ( 1 /default ) * (Cells_DSBDamage_G2 * k_rep_hr * pATM_G2 * pATR_G2)-1.0 * ( 1 /default ) * (Cells_DSBDamage_G2 * k_apo_dsb_g2 * pp53_G2)-1.0 * ( 1 /default ) * (Cells_DSBDamage_G2 * kt_apo / (k_cyc_arr_chk2 * pChk2_G2 + 1)),
    D(Cells_Apo) ~ -1.0 * ( 1 /default ) * (Cells_Apo * k_death_delay)-1.0 * ( 1 /default ) * (Cells_Apo * k_dox_apo * ((Dox_level)^(6) / ((Dox_level)^(6) + (k_dox_kd)^(6))))+1.0 * ( 1 /default ) * (Cells_SSBDamage_S * k_apo_ssb * pp53_S)+1.0 * ( 1 /default ) * (Cells_DSBDamage_S * k_apo_dsb_s * pp53_S)+1.0 * ( 1 /default ) * (Cells_DSBDamage_G2 * k_apo_dsb_g2 * pp53_G2)+2.0 * ( 1 /default ) * (Cells_SSBDamage_S * kt_apo / (k_cyc_arr_chk1 * pChk1_S + 1))+2.0 * ( 1 /default ) * (Cells_DSBDamage_S * kt_apo / (k_cyc_arr_chk1 * pChk1_S + 1))+2.0 * ( 1 /default ) * (Cells_DSBDamage_G2 * kt_apo / (k_cyc_arr_chk2 * pChk2_G2 + 1)),
    D(pChk2_G2) ~ +1.0 * ( 1 /default ) * (-pATM_G2 * p_chk2_act * (pChk2_G2 - 1))-1.0 * ( 1 /default ) * (pChk2_G2 / 24)-1.0 * ( 1 /default ) * (Wip1_G2 * pChk2_G2 * p_chk2_dea_wip1),
    D(Wip1_G2) ~ +1.0 * ( 1 /default ) * (Wip1_mRNA_G2 * p_wip1_turn)-1.0 * ( 1 /default ) * (Wip1_G2 * p_wip1_turn),
    D(pATR_S) ~ +1.0 * ( 1 /default ) * (-(Cells_SSBDamage_S * p_atr_act_ssb * (pATR_S - 1)) / (Cells + Cells_Apo + Cells_Apo1 + Cells_Apo2 + Cells_Apo3 + Cells_Apo4 + Cells_Cycle_G2 + Cells_Cycle_S + Cells_Apo_ReOx + Cells_DSBDamage_G2 + Cells_DSBDamage_S + Cells_SSBDamage_S))+1.0 * ( 1 /default ) * (-pATM_S * p_atr_act_atm * (pATR_S - 1))-1.0 * ( 1 /default ) * (pATR_S / 24),
    D(p21_S) ~ +1.0 * ( 1 /default ) * (p21_mRNA_S * p_p21_turn)-1.0 * ( 1 /default ) * (p21_S * p_p21_turn),
    D(Cells_DSBDamage_S) ~ -1.0 * ( 1 /default ) * (Cells_DSBDamage_S * k_dox_apo * ((Dox_level)^(6) / ((Dox_level)^(6) + (k_dox_kd)^(6))))+1.0 * ( 1 /default ) * (Cells_SSBDamage_S * (k_ssb_to_dsb + k_ssb_to_dsb_sn38 * ifelse(SN38_level > 0, 1, 0)))-1.0 * ( 1 /default ) * (Cells_DSBDamage_S * k_rep_nhej * pDNAPK_S / (Cells_DSBDamage_S + k_rep_nhej_sat))-1.0 * ( 1 /default ) * (Cells_DSBDamage_S * k_apo_dsb_s * pp53_S)-1.0 * ( 1 /default ) * (Cells_DSBDamage_S * kt_apo / (k_cyc_arr_chk1 * pChk1_S + 1)),
    D(p21_mRNA_G2) ~ +1.0 * ( 1 /default ) * (p_p21_mrna_turn * pp53_G2 / (p_mrna_exp_inh_dox * ((Dox_level)^(6) / ((Dox_level)^(6) + (p_mrna_exp_inh_dox_kd)^(6))) + 1))-1.0 * ( 1 /default ) * (p21_mRNA_G2 * p_p21_mrna_turn),
    D(Space) ~ -1.0 * ( 1 /default ) * (Cells * Space * kt / init_Space)+1.0 * ( 1 /default ) * (Cells_Dead * k_lyse),
    D(pChk2_S) ~ +1.0 * ( 1 /default ) * (-pATM_S * p_chk2_act * (pChk2_S - 1))-1.0 * ( 1 /default ) * (pChk2_S / 24)-1.0 * ( 1 /default ) * (Wip1_S * pChk2_S * p_chk2_dea_wip1),
    D(Cells_Cycle_S) ~ +1.0 * ( 1 /default ) * (Cells * Space * kt / init_Space)-1.0 * ( 1 /default ) * (Cells_Cycle_S * kt * (init_Cells_Cycle_G2_rel - 1) * (init_Cells_Cycle_S_rel - 1) * (init_Cells * init_Cells_Cycle_S_rel + init_Cells * (init_Cells_Cycle_G2_rel - 1) * (init_Cells_Cycle_S_rel - 1) - 2 * init_Cells * init_Cells_Cycle_G2_rel * (init_Cells_Cycle_S_rel - 1)) / (init_Cells_Cycle_S_rel * (2 * init_Cells * init_Cells_Cycle_S_rel + init_Cells * (init_Cells_Cycle_G2_rel - 1) * (init_Cells_Cycle_S_rel - 1) - 2 * init_Cells * init_Cells_Cycle_G2_rel * (init_Cells_Cycle_S_rel - 1))))-1.0 * ( 1 /default ) * (Cells_Cycle_S * k_dox_apo * ((Dox_level)^(6) / ((Dox_level)^(6) + (k_dox_kd)^(6))))-1.0 * ( 1 /default ) * (Gem_level * Cells_Cycle_S * k_damage_gem_ssb)-1.0 * ( 1 /default ) * (SN38_level * Cells_Cycle_S * k_damage_sn38_ssb)+1.0 * ( 1 /default ) * (Cells_SSBDamage_S * k_rep_s * pATR_S)+1.0 * ( 1 /default ) * (Cells_DSBDamage_S * k_rep_nhej * pDNAPK_S / (Cells_DSBDamage_S + k_rep_nhej_sat)),
    D(Cells_Apo3) ~ +1.0 * ( 1 /default ) * (Cells_Apo2 * k_death_delay)-1.0 * ( 1 /default ) * (Cells_Apo3 * k_death_delay)-1.0 * ( 1 /default ) * (Cells_Apo3 * k_dox_apo * ((Dox_level)^(6) / ((Dox_level)^(6) + (k_dox_kd)^(6)))),
    D(pChk1_G2) ~ +1.0 * ( 1 /default ) * (-pATR_G2 * p_chk1_act * (pChk1_G2 - 1))-1.0 * ( 1 /default ) * (pChk1_G2 / 24)-1.0 * ( 1 /default ) * (Wip1_G2 * pChk1_G2 * p_chk1_dea_wip1),
    D(pDNAPK_S) ~ +1.0 * ( 1 /default ) * (-(Cells_DSBDamage_S * p_dnapk_act * (pDNAPK_S - 1)) / (Cells + Cells_Apo + Cells_Apo1 + Cells_Apo2 + Cells_Apo3 + Cells_Apo4 + Cells_Cycle_G2 + Cells_Cycle_S + Cells_Apo_ReOx + Cells_DSBDamage_G2 + Cells_DSBDamage_S + Cells_SSBDamage_S))-1.0 * ( 1 /default ) * (pDNAPK_S / 24)-1.0 * ( 1 /default ) * (Wip1_S * pDNAPK_S * p_dnapk_dea_wip1),
    D(pChk1_S) ~ +1.0 * ( 1 /default ) * (-pATR_S * p_chk1_act * (pChk1_S - 1))-1.0 * ( 1 /default ) * (pChk1_S / 24)-1.0 * ( 1 /default ) * (Wip1_S * pChk1_S * p_chk1_dea_wip1),
    D(pDNAPK_G2) ~ +1.0 * ( 1 /default ) * (-(Cells_DSBDamage_G2 * p_dnapk_act * (pDNAPK_G2 - 1)) / (Cells + Cells_Apo + Cells_Apo1 + Cells_Apo2 + Cells_Apo3 + Cells_Apo4 + 2 * Cells_Cycle_G2 + Cells_Apo_ReOx + Cells_DSBDamage_G2 + Cells_DSBDamage_S + Cells_SSBDamage_S))-1.0 * ( 1 /default ) * (pDNAPK_G2 / 24)-1.0 * ( 1 /default ) * (Wip1_G2 * pDNAPK_G2 * p_dnapk_dea_wip1),
    D(pp53_G2) ~ +1.0 * ( 1 /default ) * (-pATR_G2 * p_p53_act_atr * (pp53_G2 - 1))+1.0 * ( 1 /default ) * (-pChk1_G2 * p_p53_act_chk1 * (pp53_G2 - 1))+1.0 * ( 1 /default ) * (-pATM_G2 * p_p53_act_atm * (pp53_G2 - 1))+1.0 * ( 1 /default ) * (-pChk2_G2 * p_p53_act_chk2 * (pp53_G2 - 1))-1.0 * ( 1 /default ) * (pp53_G2 / 24),
    D(Cells) ~ -1.0 * ( 1 /default ) * (Cells * Space * kt / init_Space)+2.0 * ( 1 /default ) * (-(Cells_Cycle_G2 * init_Cells * kt * (init_Cells_Cycle_G2_rel - 1) * (init_Cells_Cycle_S_rel - 1) * (init_Cells * init_Cells_Cycle_S_rel + init_Cells * (init_Cells_Cycle_G2_rel - 1) * (init_Cells_Cycle_S_rel - 1) - init_Cells * init_Cells_Cycle_G2_rel * (init_Cells_Cycle_S_rel - 1))) / ((init_Cells)^(2) * init_Cells_Cycle_G2_rel * (init_Cells_Cycle_G2_rel - 1) * (init_Cells_Cycle_S_rel - 1)^(2) - 2 * (init_Cells)^(2) * (init_Cells_Cycle_G2_rel)^(2) * (init_Cells_Cycle_S_rel - 1)^(2) + 2 * (init_Cells)^(2) * init_Cells_Cycle_G2_rel * init_Cells_Cycle_S_rel * (init_Cells_Cycle_S_rel - 1)))-1.0 * ( 1 /default ) * (Cells * k_dox_apo * ((Dox_level)^(6) / ((Dox_level)^(6) + (k_dox_kd)^(6))))
    ]

    @named sys = ODESystem(eqs, t, stateArray, parameterArray)

    ### Initial species concentrations ###
    initialSpeciesValues = [
    Cells_Dead => 0.0,
    Cells_Apo1 => 0.0,
    pATR_G2 => 0.0,
    yH2AX_G2 => 0.0,
    Wip1_mRNA_S => 0.0,
    pATM_S => 0.0,
    p21_G2 => 0.0,
    Cells_SSBDamage_S => 0.0,
    pp53_S => 0.0,
    pATM_G2 => 0.0,
    Cells_Apo2 => 0.0,
    Cells_Cycle_G2 => -init_Cells * init_Cells_Cycle_G2_rel * (init_Cells_Cycle_S_rel - 1),
    Cells_Apo_ReOx => 0.0,
    yH2AX_S => 0.0,
    p21_mRNA_S => 0.0,
    Wip1_S => 0.0,
    Wip1_mRNA_G2 => 0.0,
    Cells_Apo4 => 0.0,
    Cells_DSBDamage_G2 => 0.0,
    Cells_Apo => 0.0,
    pChk2_G2 => 0.0,
    Wip1_G2 => 0.0,
    pATR_S => 0.0,
    p21_S => 0.0,
    Cells_DSBDamage_S => 0.0,
    p21_mRNA_G2 => 0.0,
    Space => init_Space,
    pChk2_S => 0.0,
    Cells_Cycle_S => init_Cells * init_Cells_Cycle_S_rel,
    Cells_Apo3 => 0.0,
    pChk1_G2 => 0.0,
    pDNAPK_S => 0.0,
    pChk1_S => 0.0,
    pDNAPK_G2 => 0.0,
    pp53_G2 => 0.0,
    Cells => init_Cells * (init_Cells_Cycle_G2_rel - 1) * (init_Cells_Cycle_S_rel - 1)
    ]

    ### SBML file parameter values ###
    trueParameterValues = [
    k_dox_apo => 0.0488009964455022,
    p_p53_act_atm => 0.172667542315577,
    p_chk1_dea_wip1 => 9.99999998885131,
    k_rep_hr => 1.00000000000009e-5,
    p_p53_act_chk1 => 1.00000000000855e-5,
    k_rep_nhej => 0.324197628164391,
    k_death => 0.103177842663178,
    k_dox_kd => 0.356451133426244,
    init_Cells_Cycle_G2_rel => 0.4,
    k_apo_dsb_g2 => 1.0e-5,
    k_rep_s => 1.00000000000009e-5,
    p_p21_turn => 0.999999999999921,
    p_h2ax_act_atm => 0.00135444548009859,
    p_p53_act_chk2 => 0.0512391341519817,
    p_atr_act_ssb => 0.0791834194872235,
    p_mrna_exp_inh_dox_kd => 0.356451179497663,
    p_chk2_dea_wip1 => 9.9999999999992,
    k_cyc_arr_chk2 => 1000.0,
    p_wip1_mrna_turn => 3.40101630069606,
    p_mrna_exp_inh_dox => 999.999999999919,
    Gem_level => 0.1,
    k_ssb_to_dsb => 0.0304996268234734,
    k_rep_nhej_sat => 0.00686700298036127,
    Dox_level => 0.0,
    k_lyse => 10.0,
    k_damage_sn38_ssb => 174.434737098025,
    p_chk1_act => 8.43249925096952,
    p_atr_act_atm => 0.0902182570058027,
    default => 1.0,
    k_apo_ssb => 1.0e-5,
    k_apo_dsb_s => 1.0e-5,
    k_cyc_arr_chk1 => 999.99999999992,
    p_h2ax_act_atr => 0.00537764161229014,
    k_death_delay => 0.141041612206608,
    p_h2ax_act_dnapk => 1.00000000000962e-5,
    k_damage_dox_dsb => 214.760281583402,
    p_chk2_act => 6.89308197813853,
    p_atm_act_dsb => 0.999999999999924,
    kt => 0.187392093037445,
    p_dnapk_act => 9.99999764924192,
    k_damage_gem_ssb => 49.8700470112818,
    p_p53_act_atr => 1.00000000000009e-5,
    init_Cells_Cycle_S_rel => 0.5,
    k_death_reox => 0.034124623752161,
    kt_apo => 0.123501569444653,
    k_ssb_to_dsb_sn38 => 999.999674332206,
    init_Space => 100.0,
    init_Cells => 1.0,
    p_dnapk_dea_wip1 => 6.58806922708483,
    SN38_level => 0.0,
    p_h2ax_dea => 1.00003466144413e-5,
    p_p21_mrna_turn => 0.109354404673231,
    p_wip1_turn => 3.41164292406599
    ]

    return sys, initialSpeciesValues, trueParameterValues

end
