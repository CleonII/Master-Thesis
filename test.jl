# Number of parameters: 52
# Number of species: 36

### Extra functions
pow(a,b) = a^b

### True parameter values ###
trueParameterValues = [0.1, 0.0, 0.0, 1.0, 0.4, 0.5, 100.0, 174.434737098025, 999.999674332206, 214.760281583402, 0.0488009964455022, 0.356451133426244, 1.0e-5, 1.0e-5, 1.0e-5, 999.99999999992, 1000.0, 49.8700470112818, 0.103177842663178, 0.141041612206608, 0.034124623752161, 10.0, 1.00000000000009e-5, 0.324197628164391, 0.00686700298036127, 1.00000000000009e-5, 0.0304996268234734, 0.187392093037445, 0.123501569444653, 999.999999999919, 0.356451179497663, 0.999999999999924, 0.0902182570058027, 0.0791834194872235, 8.43249925096952, 9.99999998885131, 6.89308197813853, 9.9999999999992, 9.99999764924192, 6.58806922708483, 0.00135444548009859, 0.00537764161229014, 1.00000000000962e-5, 1.00003466144413e-5, 0.109354404673231, 0.999999999999921, 0.172667542315577, 1.00000000000009e-5, 1.00000000000855e-5, 0.0512391341519817, 3.40101630069606, 3.41164292406599]

### Initial species concentrations ###
initialSpeciesValue = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

function ODE_model!(du, u, p, t)

### Parameters ###
Gem_level, Dox_level, SN38_level, init_Cells, init_Cells_Cycle_G2_rel, init_Cells_Cycle_S_rel, init_Space, k_damage_sn38_ssb, k_ssb_to_dsb_sn38, k_damage_dox_dsb, k_dox_apo, k_dox_kd, k_apo_dsb_g2, k_apo_dsb_s, k_apo_ssb, k_cyc_arr_chk1, k_cyc_arr_chk2, k_damage_gem_ssb, k_death, k_death_delay, k_death_reox, k_lyse, k_rep_hr, k_rep_nhej, k_rep_nhej_sat, k_rep_s, k_ssb_to_dsb, kt, kt_apo, p_mrna_exp_inh_dox, p_mrna_exp_inh_dox_kd, p_atm_act_dsb, p_atr_act_atm, p_atr_act_ssb, p_chk1_act, p_chk1_dea_wip1, p_chk2_act, p_chk2_dea_wip1, p_dnapk_act, p_dnapk_dea_wip1, p_h2ax_act_atm, p_h2ax_act_atr, p_h2ax_act_dnapk, p_h2ax_dea, p_p21_mrna_turn, p_p21_turn, p_p53_act_atm, p_p53_act_atr, p_p53_act_chk1, p_p53_act_chk2, p_wip1_mrna_turn, p_wip1_turn = p

### Compartments ###
default = 1.0

### Species ###
Space, Cells, Cells_Cycle_S, Cells_Cycle_G2, Cells_SSBDamage_S, Cells_DSBDamage_S, Cells_DSBDamage_G2, Cells_Apo, Cells_Apo1, Cells_Apo2, Cells_Apo3, Cells_Apo4, Cells_Apo_ReOx, Cells_Dead, pATR_S, pATR_G2, pChk1_S, pChk1_G2, pATM_S, pATM_G2, pChk2_S, pChk2_G2, pDNAPK_S, pDNAPK_G2, pp53_S, pp53_G2, yH2AX_S, yH2AX_G2, p21_mRNA_S, p21_S, p21_mRNA_G2, p21_G2, Wip1_mRNA_S, Wip1_S, Wip1_mRNA_G2, Wip1_G2 = u

### Function definitions ###

### Events ###

### Rules ###

### Derivatives ###
du[1] = -1.0 * (Cells * Space * kt / init_Space)+1.0 * (Cells_Dead * k_lyse)
du[2] = -1.0 * (Cells * Space * kt / init_Space)+2.0 * (-(Cells_Cycle_G2 * init_Cells * kt * (init_Cells_Cycle_G2_rel - 1) * (init_Cells_Cycle_S_rel - 1) * (init_Cells * init_Cells_Cycle_S_rel + init_Cells * (init_Cells_Cycle_G2_rel - 1) * (init_Cells_Cycle_S_rel - 1) - init_Cells * init_Cells_Cycle_G2_rel * (init_Cells_Cycle_S_rel - 1))) / (pow(init_Cells, 2) * init_Cells_Cycle_G2_rel * (init_Cells_Cycle_G2_rel - 1) * pow(init_Cells_Cycle_S_rel - 1, 2) - 2 * pow(init_Cells, 2) * pow(init_Cells_Cycle_G2_rel, 2) * pow(init_Cells_Cycle_S_rel - 1, 2) + 2 * pow(init_Cells, 2) * init_Cells_Cycle_G2_rel * init_Cells_Cycle_S_rel * (init_Cells_Cycle_S_rel - 1)))-1.0 * (Cells * k_dox_apo * (pow(Dox_level, 6) / (pow(Dox_level, 6) + pow(k_dox_kd, 6))))
du[3] = +1.0 * (Cells * Space * kt / init_Space)-1.0 * (Cells_Cycle_S * kt * (init_Cells_Cycle_G2_rel - 1) * (init_Cells_Cycle_S_rel - 1) * (init_Cells * init_Cells_Cycle_S_rel + init_Cells * (init_Cells_Cycle_G2_rel - 1) * (init_Cells_Cycle_S_rel - 1) - 2 * init_Cells * init_Cells_Cycle_G2_rel * (init_Cells_Cycle_S_rel - 1)) / (init_Cells_Cycle_S_rel * (2 * init_Cells * init_Cells_Cycle_S_rel + init_Cells * (init_Cells_Cycle_G2_rel - 1) * (init_Cells_Cycle_S_rel - 1) - 2 * init_Cells * init_Cells_Cycle_G2_rel * (init_Cells_Cycle_S_rel - 1))))-1.0 * (Cells_Cycle_S * k_dox_apo * (pow(Dox_level, 6) / (pow(Dox_level, 6) + pow(k_dox_kd, 6))))-1.0 * (Gem_level * Cells_Cycle_S * k_damage_gem_ssb)-1.0 * (SN38_level * Cells_Cycle_S * k_damage_sn38_ssb)+1.0 * (Cells_SSBDamage_S * k_rep_s * pATR_S)+1.0 * (Cells_DSBDamage_S * k_rep_nhej * pDNAPK_S / (Cells_DSBDamage_S + k_rep_nhej_sat))
du[4] = +1.0 * (Cells_Cycle_S * kt * (init_Cells_Cycle_G2_rel - 1) * (init_Cells_Cycle_S_rel - 1) * (init_Cells * init_Cells_Cycle_S_rel + init_Cells * (init_Cells_Cycle_G2_rel - 1) * (init_Cells_Cycle_S_rel - 1) - 2 * init_Cells * init_Cells_Cycle_G2_rel * (init_Cells_Cycle_S_rel - 1)) / (init_Cells_Cycle_S_rel * (2 * init_Cells * init_Cells_Cycle_S_rel + init_Cells * (init_Cells_Cycle_G2_rel - 1) * (init_Cells_Cycle_S_rel - 1) - 2 * init_Cells * init_Cells_Cycle_G2_rel * (init_Cells_Cycle_S_rel - 1))))-1.0 * (-(Cells_Cycle_G2 * init_Cells * kt * (init_Cells_Cycle_G2_rel - 1) * (init_Cells_Cycle_S_rel - 1) * (init_Cells * init_Cells_Cycle_S_rel + init_Cells * (init_Cells_Cycle_G2_rel - 1) * (init_Cells_Cycle_S_rel - 1) - init_Cells * init_Cells_Cycle_G2_rel * (init_Cells_Cycle_S_rel - 1))) / (pow(init_Cells, 2) * init_Cells_Cycle_G2_rel * (init_Cells_Cycle_G2_rel - 1) * pow(init_Cells_Cycle_S_rel - 1, 2) - 2 * pow(init_Cells, 2) * pow(init_Cells_Cycle_G2_rel, 2) * pow(init_Cells_Cycle_S_rel - 1, 2) + 2 * pow(init_Cells, 2) * init_Cells_Cycle_G2_rel * init_Cells_Cycle_S_rel * (init_Cells_Cycle_S_rel - 1)))-1.0 * (Cells_Cycle_G2 * k_dox_apo * (pow(Dox_level, 6) / (pow(Dox_level, 6) + pow(k_dox_kd, 6))))-1.0 * (Cells_Cycle_G2 * k_damage_dox_dsb * Dox_level)+1.0 * (Cells_DSBDamage_G2 * k_rep_nhej * pDNAPK_G2 / (Cells_DSBDamage_G2 + k_rep_nhej_sat))+1.0 * (Cells_DSBDamage_G2 * k_rep_hr * pATM_G2 * pATR_G2)
du[5] = -1.0 * (Cells_SSBDamage_S * k_dox_apo * (pow(Dox_level, 6) / (pow(Dox_level, 6) + pow(k_dox_kd, 6))))+1.0 * (Gem_level * Cells_Cycle_S * k_damage_gem_ssb)+1.0 * (SN38_level * Cells_Cycle_S * k_damage_sn38_ssb)-1.0 * (Cells_SSBDamage_S * (k_ssb_to_dsb + k_ssb_to_dsb_sn38 * begin; if SN38_level > 0; 1; elseif SN38_level <= 0; 0; end; end))-1.0 * (Cells_SSBDamage_S * k_rep_s * pATR_S)-1.0 * (Cells_SSBDamage_S * k_apo_ssb * pp53_S)-1.0 * (Cells_SSBDamage_S * kt_apo / (k_cyc_arr_chk1 * pChk1_S + 1))
du[6] = -1.0 * (Cells_DSBDamage_S * k_dox_apo * (pow(Dox_level, 6) / (pow(Dox_level, 6) + pow(k_dox_kd, 6))))+1.0 * (Cells_SSBDamage_S * (k_ssb_to_dsb + k_ssb_to_dsb_sn38 * begin; if SN38_level > 0; 1; elseif SN38_level <= 0; 0; end; end))-1.0 * (Cells_DSBDamage_S * k_rep_nhej * pDNAPK_S / (Cells_DSBDamage_S + k_rep_nhej_sat))-1.0 * (Cells_DSBDamage_S * k_apo_dsb_s * pp53_S)-1.0 * (Cells_DSBDamage_S * kt_apo / (k_cyc_arr_chk1 * pChk1_S + 1))
du[7] = -1.0 * (Cells_DSBDamage_G2 * k_dox_apo * (pow(Dox_level, 6) / (pow(Dox_level, 6) + pow(k_dox_kd, 6))))+1.0 * (Cells_Cycle_G2 * k_damage_dox_dsb * Dox_level)-1.0 * (Cells_DSBDamage_G2 * k_rep_nhej * pDNAPK_G2 / (Cells_DSBDamage_G2 + k_rep_nhej_sat))-1.0 * (Cells_DSBDamage_G2 * k_rep_hr * pATM_G2 * pATR_G2)-1.0 * (Cells_DSBDamage_G2 * k_apo_dsb_g2 * pp53_G2)-1.0 * (Cells_DSBDamage_G2 * kt_apo / (k_cyc_arr_chk2 * pChk2_G2 + 1))
du[8] = -1.0 * (Cells_Apo * k_death_delay)-1.0 * (Cells_Apo * k_dox_apo * (pow(Dox_level, 6) / (pow(Dox_level, 6) + pow(k_dox_kd, 6))))+1.0 * (Cells_SSBDamage_S * k_apo_ssb * pp53_S)+1.0 * (Cells_DSBDamage_S * k_apo_dsb_s * pp53_S)+1.0 * (Cells_DSBDamage_G2 * k_apo_dsb_g2 * pp53_G2)+2.0 * (Cells_SSBDamage_S * kt_apo / (k_cyc_arr_chk1 * pChk1_S + 1))+2.0 * (Cells_DSBDamage_S * kt_apo / (k_cyc_arr_chk1 * pChk1_S + 1))+2.0 * (Cells_DSBDamage_G2 * kt_apo / (k_cyc_arr_chk2 * pChk2_G2 + 1))
du[9] = +1.0 * (Cells_Apo * k_death_delay)-1.0 * (Cells_Apo1 * k_death_delay)-1.0 * (Cells_Apo1 * k_dox_apo * (pow(Dox_level, 6) / (pow(Dox_level, 6) + pow(k_dox_kd, 6))))
du[10] = +1.0 * (Cells_Apo1 * k_death_delay)-1.0 * (Cells_Apo2 * k_death_delay)-1.0 * (Cells_Apo2 * k_dox_apo * (pow(Dox_level, 6) / (pow(Dox_level, 6) + pow(k_dox_kd, 6))))
du[11] = +1.0 * (Cells_Apo2 * k_death_delay)-1.0 * (Cells_Apo3 * k_death_delay)-1.0 * (Cells_Apo3 * k_dox_apo * (pow(Dox_level, 6) / (pow(Dox_level, 6) + pow(k_dox_kd, 6))))
du[12] = +1.0 * (Cells_Apo3 * k_death_delay)-1.0 * (Cells_Apo4 * k_death)-1.0 * (Cells_Apo4 * k_dox_apo * (pow(Dox_level, 6) / (pow(Dox_level, 6) + pow(k_dox_kd, 6))))
du[13] = +1.0 * (Cells * k_dox_apo * (pow(Dox_level, 6) / (pow(Dox_level, 6) + pow(k_dox_kd, 6))))+1.0 * (Cells_Cycle_S * k_dox_apo * (pow(Dox_level, 6) / (pow(Dox_level, 6) + pow(k_dox_kd, 6))))+1.0 * (Cells_Cycle_G2 * k_dox_apo * (pow(Dox_level, 6) / (pow(Dox_level, 6) + pow(k_dox_kd, 6))))+1.0 * (Cells_SSBDamage_S * k_dox_apo * (pow(Dox_level, 6) / (pow(Dox_level, 6) + pow(k_dox_kd, 6))))+1.0 * (Cells_DSBDamage_S * k_dox_apo * (pow(Dox_level, 6) / (pow(Dox_level, 6) + pow(k_dox_kd, 6))))+1.0 * (Cells_DSBDamage_G2 * k_dox_apo * (pow(Dox_level, 6) / (pow(Dox_level, 6) + pow(k_dox_kd, 6))))+1.0 * (Cells_Apo * k_dox_apo * (pow(Dox_level, 6) / (pow(Dox_level, 6) + pow(k_dox_kd, 6))))+1.0 * (Cells_Apo1 * k_dox_apo * (pow(Dox_level, 6) / (pow(Dox_level, 6) + pow(k_dox_kd, 6))))+1.0 * (Cells_Apo2 * k_dox_apo * (pow(Dox_level, 6) / (pow(Dox_level, 6) + pow(k_dox_kd, 6))))+1.0 * (Cells_Apo3 * k_dox_apo * (pow(Dox_level, 6) / (pow(Dox_level, 6) + pow(k_dox_kd, 6))))+1.0 * (Cells_Apo4 * k_dox_apo * (pow(Dox_level, 6) / (pow(Dox_level, 6) + pow(k_dox_kd, 6))))-1.0 * (Cells_Apo_ReOx * k_death_reox)
du[14] = +1.0 * (Cells_Apo4 * k_death)+1.0 * (Cells_Apo_ReOx * k_death_reox)-1.0 * (Cells_Dead * k_lyse)
du[15] = +1.0 * (-(Cells_SSBDamage_S * p_atr_act_ssb * (pATR_S - 1)) / (Cells + Cells_Apo + Cells_Apo1 + Cells_Apo2 + Cells_Apo3 + Cells_Apo4 + Cells_Cycle_G2 + Cells_Cycle_S + Cells_Apo_ReOx + Cells_DSBDamage_G2 + Cells_DSBDamage_S + Cells_SSBDamage_S))+1.0 * (-pATM_S * p_atr_act_atm * (pATR_S - 1))-1.0 * (pATR_S / 24)
du[16] = +1.0 * (-pATM_G2 * p_atr_act_atm * (pATR_G2 - 1))-1.0 * (pATR_G2 / 24)
du[17] = +1.0 * (-pATR_S * p_chk1_act * (pChk1_S - 1))-1.0 * (pChk1_S / 24)-1.0 * (Wip1_S * pChk1_S * p_chk1_dea_wip1)
du[18] = +1.0 * (-pATR_G2 * p_chk1_act * (pChk1_G2 - 1))-1.0 * (pChk1_G2 / 24)-1.0 * (Wip1_G2 * pChk1_G2 * p_chk1_dea_wip1)
du[19] = +1.0 * (-(Cells_DSBDamage_S * p_atm_act_dsb * (pATM_S - 1)) / (Cells + Cells_Apo + Cells_Apo1 + Cells_Apo2 + Cells_Apo3 + Cells_Apo4 + Cells_Cycle_G2 + Cells_Cycle_S + Cells_Apo_ReOx + Cells_DSBDamage_G2 + Cells_DSBDamage_S + Cells_SSBDamage_S))-1.0 * (pATM_S / 24)
du[20] = +1.0 * (-(Cells_DSBDamage_G2 * p_atm_act_dsb * (pATM_G2 - 1)) / (Cells + Cells_Apo + Cells_Apo1 + Cells_Apo2 + Cells_Apo3 + Cells_Apo4 + 2 * Cells_Cycle_G2 + Cells_Apo_ReOx + Cells_DSBDamage_G2 + Cells_DSBDamage_S + Cells_SSBDamage_S))-1.0 * (pATM_G2 / 24)
du[21] = +1.0 * (-pATM_S * p_chk2_act * (pChk2_S - 1))-1.0 * (pChk2_S / 24)-1.0 * (Wip1_S * pChk2_S * p_chk2_dea_wip1)
du[22] = +1.0 * (-pATM_G2 * p_chk2_act * (pChk2_G2 - 1))-1.0 * (pChk2_G2 / 24)-1.0 * (Wip1_G2 * pChk2_G2 * p_chk2_dea_wip1)
du[23] = +1.0 * (-(Cells_DSBDamage_S * p_dnapk_act * (pDNAPK_S - 1)) / (Cells + Cells_Apo + Cells_Apo1 + Cells_Apo2 + Cells_Apo3 + Cells_Apo4 + Cells_Cycle_G2 + Cells_Cycle_S + Cells_Apo_ReOx + Cells_DSBDamage_G2 + Cells_DSBDamage_S + Cells_SSBDamage_S))-1.0 * (pDNAPK_S / 24)-1.0 * (Wip1_S * pDNAPK_S * p_dnapk_dea_wip1)
du[24] = +1.0 * (-(Cells_DSBDamage_G2 * p_dnapk_act * (pDNAPK_G2 - 1)) / (Cells + Cells_Apo + Cells_Apo1 + Cells_Apo2 + Cells_Apo3 + Cells_Apo4 + 2 * Cells_Cycle_G2 + Cells_Apo_ReOx + Cells_DSBDamage_G2 + Cells_DSBDamage_S + Cells_SSBDamage_S))-1.0 * (pDNAPK_G2 / 24)-1.0 * (Wip1_G2 * pDNAPK_G2 * p_dnapk_dea_wip1)
du[25] = +1.0 * (-pATR_S * p_p53_act_atr * (pp53_S - 1))+1.0 * (-pChk1_S * p_p53_act_chk1 * (pp53_S - 1))+1.0 * (-pATM_S * p_p53_act_atm * (pp53_S - 1))+1.0 * (-pChk2_S * p_p53_act_chk2 * (pp53_S - 1))-1.0 * (pp53_S / 24)
du[26] = +1.0 * (-pATR_G2 * p_p53_act_atr * (pp53_G2 - 1))+1.0 * (-pChk1_G2 * p_p53_act_chk1 * (pp53_G2 - 1))+1.0 * (-pATM_G2 * p_p53_act_atm * (pp53_G2 - 1))+1.0 * (-pChk2_G2 * p_p53_act_chk2 * (pp53_G2 - 1))-1.0 * (pp53_G2 / 24)
du[27] = +1.0 * (-pATM_S * p_h2ax_act_atm * (yH2AX_S - 1))+1.0 * (-pATR_S * p_h2ax_act_atr * (yH2AX_S - 1))+1.0 * (-pDNAPK_S * p_h2ax_act_dnapk * (yH2AX_S - 1))-1.0 * (p_h2ax_dea * yH2AX_S)
du[28] = +1.0 * (-pATM_G2 * p_h2ax_act_atm * (yH2AX_G2 - 1))+1.0 * (-pATR_G2 * p_h2ax_act_atr * (yH2AX_G2 - 1))+1.0 * (-pDNAPK_G2 * p_h2ax_act_dnapk * (yH2AX_G2 - 1))-1.0 * (p_h2ax_dea * yH2AX_G2)
du[29] = +1.0 * (p_p21_mrna_turn * pp53_S / (p_mrna_exp_inh_dox * (pow(Dox_level, 6) / (pow(Dox_level, 6) + pow(p_mrna_exp_inh_dox_kd, 6))) + 1))-1.0 * (p21_mRNA_S * p_p21_mrna_turn)
du[30] = +1.0 * (p21_mRNA_S * p_p21_turn)-1.0 * (p21_S * p_p21_turn)
du[31] = +1.0 * (p_p21_mrna_turn * pp53_G2 / (p_mrna_exp_inh_dox * (pow(Dox_level, 6) / (pow(Dox_level, 6) + pow(p_mrna_exp_inh_dox_kd, 6))) + 1))-1.0 * (p21_mRNA_G2 * p_p21_mrna_turn)
du[32] = +1.0 * (p21_mRNA_G2 * p_p21_turn)-1.0 * (p21_G2 * p_p21_turn)
du[33] = +1.0 * (p_wip1_mrna_turn * pp53_S / (p_mrna_exp_inh_dox * (pow(Dox_level, 6) / (pow(Dox_level, 6) + pow(p_mrna_exp_inh_dox_kd, 6))) + 1))-1.0 * (Wip1_mRNA_S * p_wip1_mrna_turn)
du[34] = +1.0 * (Wip1_mRNA_S * p_wip1_turn)-1.0 * (Wip1_S * p_wip1_turn)
du[35] = +1.0 * (p_wip1_mrna_turn * pp53_G2 / (p_mrna_exp_inh_dox * (pow(Dox_level, 6) / (pow(Dox_level, 6) + pow(p_mrna_exp_inh_dox_kd, 6))) + 1))-1.0 * (Wip1_mRNA_G2 * p_wip1_mrna_turn)
du[36] = +1.0 * (Wip1_mRNA_G2 * p_wip1_turn)-1.0 * (Wip1_G2 * p_wip1_turn)

nothing
end
