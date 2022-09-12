# Model name: model_Lucarelli_CellSystems2018
# Number of parameters: 106
# Number of species: 33
function getODEModel_model_Lucarelli_CellSystems2018()

    ### Define constant parameters
    cell = 13.17

    ### Define independent and dependent variables
    ModelingToolkit.@variables t ppS3_ppS3_ppS3(t) ppS3_S4_S4(t) geneH(t) geneI(t) geneJ(t) TGFb_pRec(t) geneE(t) S4(t) ppS2_ppS2_ppS2(t) geneK(t) pS2(t) pS3(t) geneC(t) ppS2_ppS2_ppS3(t) geneF(t) S4_S4_S4(t) TGFb(t) S3(t) S2(t) S2_S4_S4(t) ppS3_ppS3_S4(t) geneB(t) ppS3(t) geneA(t) geneD(t) ppS2_S4_S4(t) geneL(t) ppS2_ppS3_S4(t) geneG(t) Rec(t) ppS2_ppS2_S4(t) ppS2_ppS3_ppS3(t) ppS2(t)

    ### Define variable parameters

    ### Define dummy variable
    ModelingToolkit.@variables dummyVariable(t)

    ### Define parameters
    ModelingToolkit.@parameters geneC_inh3 geneJ_act1 geneG_inh1 init_Rec geneC_act3 geneA_inh2 geneK_turn geneA_act3 geneI_act2 geneA_act2 geneD_inh2 k_223 k_233 geneB_act1 S_dephosphos geneG_act1 geneH_inh2 geneI_act1 geneB_inh3 geneI_inh2 geneC_turn geneJ_inh1 geneK_act3 geneJ_inh3 geneG_inh3 geneC_act2 init_TGFb geneL_inh1 geneE_inh3 geneF_act3 geneF_inh1 geneD_act3 geneE_act2 geneD_inh3 k_234 geneH_act2 geneA_turn geneL_act1 geneH_act1 geneL_inh2 geneB_turn init_S4 geneL_inh3 khomo2 geneK_act2 geneA_inh3 geneJ_act2 geneI_inh3 geneD_act1 geneJ_turn geneG_act3 geneL_act2 Rec_act geneI_act3 k_224 pRec_degind geneE_inh2 geneF_inh2 geneC_act1 geneD_inh1 k_344 init_S3 geneI_inh1 geneL_turn geneF_act2 k_on_u geneE_inh1 geneH_act3 geneH_inh1 geneK_inh3 geneE_act1 geneA_inh1 geneB_inh1 geneH_inh3 geneD_turn S_dephos geneL_act3 S_phos geneG_turn geneA_act1 geneI_turn khomo3 geneC_inh1 geneG_act2 init_S2 geneK_inh2 k_334 geneB_inh2 geneH_turn geneJ_inh2 khomo4 geneB_act2 geneF_turn geneJ_act3 geneK_inh1 geneB_act3 geneK_act1 geneE_turn geneE_act3 geneG_inh2 kdiss_SS k_244 geneD_act2 geneF_inh3 geneF_act1 geneC_inh2

    ### Define an operator for the differentiation w.r.t. time
    D = Differential(t)

    ### Events ###

    ### Derivatives ###
    eqs = [
    D(ppS3_ppS3_ppS3) ~ +1.0 * ( 1 /cell ) * (cell * khomo3 * (ppS3)^(3))-1.0 * ( 1 /cell ) * (cell * 3 * S_dephosphos * ppS3_ppS3_ppS3),
    D(ppS3_S4_S4) ~ +1.0 * ( 1 /cell ) * (cell * (S4)^(2) * k_344 * ppS3)-1.0 * ( 1 /cell ) * (cell * S_dephosphos * ppS3_S4_S4),
    D(geneH) ~ +1.0 * ( 1 /cell ) * (cell * ((geneH_turn + geneH_act2 * ppS2_S4_S4 + geneH_act1 * ppS2_ppS3_S4 + geneH_act3 * ppS2_ppS3_ppS3) / (geneH_inh2 * ppS2_S4_S4 + geneH_inh1 * ppS2_ppS3_S4 + geneH_inh3 * ppS2_ppS3_ppS3 + 1)))-1.0 * ( 1 /cell ) * (cell * geneH * geneH_turn),
    D(geneI) ~ +1.0 * ( 1 /cell ) * (cell * ((geneI_turn + geneI_act2 * ppS2_S4_S4 + geneI_act1 * ppS2_ppS3_S4 + geneI_act3 * ppS2_ppS3_ppS3) / (geneI_inh2 * ppS2_S4_S4 + geneI_inh1 * ppS2_ppS3_S4 + geneI_inh3 * ppS2_ppS3_ppS3 + 1)))-1.0 * ( 1 /cell ) * (cell * geneI * geneI_turn),
    D(geneJ) ~ +1.0 * ( 1 /cell ) * (cell * ((geneJ_turn + geneJ_act2 * ppS2_S4_S4 + geneJ_act1 * ppS2_ppS3_S4 + geneJ_act3 * ppS2_ppS3_ppS3) / (geneJ_inh2 * ppS2_S4_S4 + geneJ_inh1 * ppS2_ppS3_S4 + geneJ_inh3 * ppS2_ppS3_ppS3 + 1)))-1.0 * ( 1 /cell ) * (cell * geneJ * geneJ_turn),
    D(TGFb_pRec) ~ +1.0 * ( 1 /cell ) * (cell * Rec * Rec_act * TGFb)-1.0 * ( 1 /cell ) * (cell * TGFb_pRec * pRec_degind),
    D(geneE) ~ +1.0 * ( 1 /cell ) * (cell * ((geneE_turn + geneE_act2 * ppS2_S4_S4 + geneE_act1 * ppS2_ppS3_S4 + geneE_act3 * ppS2_ppS3_ppS3) / (geneE_inh2 * ppS2_S4_S4 + geneE_inh1 * ppS2_ppS3_S4 + geneE_inh3 * ppS2_ppS3_ppS3 + 1)))-1.0 * ( 1 /cell ) * (cell * geneE * geneE_turn),
    D(S4) ~ -2.0 * ( 1 /cell ) * (cell * S2 * (S4)^(2) * k_on_u)+2.0 * ( 1 /cell ) * (cell * S2_S4_S4 * kdiss_SS)-3.0 * ( 1 /cell ) * (cell * (S4)^(3) * khomo4)+3.0 * ( 1 /cell ) * (cell * S4_S4_S4 * kdiss_SS)-1.0 * ( 1 /cell ) * (cell * S4 * k_224 * (ppS2)^(2))+1.0 * ( 1 /cell ) * (cell * 2 * S_dephosphos * ppS2_ppS2_S4)-1.0 * ( 1 /cell ) * (cell * S4 * k_334 * (ppS3)^(2))+1.0 * ( 1 /cell ) * (cell * 2 * S_dephosphos * ppS3_ppS3_S4)-2.0 * ( 1 /cell ) * (cell * (S4)^(2) * k_244 * ppS2)+2.0 * ( 1 /cell ) * (cell * S_dephosphos * ppS2_S4_S4)-2.0 * ( 1 /cell ) * (cell * (S4)^(2) * k_344 * ppS3)+2.0 * ( 1 /cell ) * (cell * S_dephosphos * ppS3_S4_S4)-1.0 * ( 1 /cell ) * (cell * S4 * k_234 * ppS2 * ppS3)+1.0 * ( 1 /cell ) * (cell * S_dephosphos * ppS2_ppS3_S4)+1.0 * ( 1 /cell ) * (cell * S_dephosphos * ppS2_ppS3_S4),
    D(ppS2_ppS2_ppS2) ~ +1.0 * ( 1 /cell ) * (cell * khomo2 * (ppS2)^(3))-1.0 * ( 1 /cell ) * (cell * 3 * S_dephosphos * ppS2_ppS2_ppS2),
    D(geneK) ~ +1.0 * ( 1 /cell ) * (cell * ((geneK_turn + geneK_act2 * ppS2_S4_S4 + geneK_act1 * ppS2_ppS3_S4 + geneK_act3 * ppS2_ppS3_ppS3) / (geneK_inh2 * ppS2_S4_S4 + geneK_inh1 * ppS2_ppS3_S4 + geneK_inh3 * ppS2_ppS3_ppS3 + 1)))-1.0 * ( 1 /cell ) * (cell * geneK * geneK_turn),
    D(pS2) ~ +1.0 * ( 1 /cell ) * (cell * 3 * S_dephosphos * ppS2_ppS2_ppS2)+1.0 * ( 1 /cell ) * (cell * S_dephosphos * ppS2)-1.0 * ( 1 /cell ) * (cell * S_dephos * pS2)+1.0 * ( 1 /cell ) * (cell * 2 * S_dephosphos * ppS2_ppS2_ppS3)+1.0 * ( 1 /cell ) * (cell * 2 * S_dephosphos * ppS2_ppS2_S4)+1.0 * ( 1 /cell ) * (cell * S_dephosphos * ppS2_ppS3_ppS3)+1.0 * ( 1 /cell ) * (cell * S_dephosphos * ppS2_S4_S4)+1.0 * ( 1 /cell ) * (cell * S_dephosphos * ppS2_ppS3_S4),
    D(pS3) ~ +1.0 * ( 1 /cell ) * (cell * 3 * S_dephosphos * ppS3_ppS3_ppS3)+1.0 * ( 1 /cell ) * (cell * S_dephosphos * ppS3)-1.0 * ( 1 /cell ) * (cell * S_dephos * pS3)+1.0 * ( 1 /cell ) * (cell * S_dephosphos * ppS2_ppS2_ppS3)+1.0 * ( 1 /cell ) * (cell * 2 * S_dephosphos * ppS2_ppS3_ppS3)+1.0 * ( 1 /cell ) * (cell * 2 * S_dephosphos * ppS3_ppS3_S4)+1.0 * ( 1 /cell ) * (cell * S_dephosphos * ppS3_S4_S4)+1.0 * ( 1 /cell ) * (cell * S_dephosphos * ppS2_ppS3_S4),
    D(geneC) ~ +1.0 * ( 1 /cell ) * (cell * ((geneC_turn + geneC_act2 * ppS2_S4_S4 + geneC_act1 * ppS2_ppS3_S4 + geneC_act3 * ppS2_ppS3_ppS3) / (geneC_inh2 * ppS2_S4_S4 + geneC_inh1 * ppS2_ppS3_S4 + geneC_inh3 * ppS2_ppS3_ppS3 + 1)))-1.0 * ( 1 /cell ) * (cell * geneC * geneC_turn),
    D(ppS2_ppS2_ppS3) ~ +1.0 * ( 1 /cell ) * (cell * k_223 * (ppS2)^(2) * ppS3)-1.0 * ( 1 /cell ) * (cell * 2 * S_dephosphos * ppS2_ppS2_ppS3)-1.0 * ( 1 /cell ) * (cell * S_dephosphos * ppS2_ppS2_ppS3),
    D(geneF) ~ +1.0 * ( 1 /cell ) * (cell * ((geneF_turn + geneF_act2 * ppS2_S4_S4 + geneF_act1 * ppS2_ppS3_S4 + geneF_act3 * ppS2_ppS3_ppS3) / (geneF_inh2 * ppS2_S4_S4 + geneF_inh1 * ppS2_ppS3_S4 + geneF_inh3 * ppS2_ppS3_ppS3 + 1)))-1.0 * ( 1 /cell ) * (cell * geneF * geneF_turn),
    D(S4_S4_S4) ~ +1.0 * ( 1 /cell ) * (cell * (S4)^(3) * khomo4)-1.0 * ( 1 /cell ) * (cell * S4_S4_S4 * kdiss_SS),
    D(TGFb) ~ -1.0 * ( 1 /cell ) * (cell * Rec * Rec_act * TGFb),
    D(S3) ~ -1.0 * ( 1 /cell ) * (cell * S3 * S_phos * TGFb_pRec)+1.0 * ( 1 /cell ) * (cell * S_dephos * pS3),
    D(S2) ~ -1.0 * ( 1 /cell ) * (cell * S2 * (S4)^(2) * k_on_u)+1.0 * ( 1 /cell ) * (cell * S2_S4_S4 * kdiss_SS)-1.0 * ( 1 /cell ) * (cell * S2 * S_phos * TGFb_pRec)+1.0 * ( 1 /cell ) * (cell * S_dephos * pS2),
    D(S2_S4_S4) ~ +1.0 * ( 1 /cell ) * (cell * S2 * (S4)^(2) * k_on_u)-1.0 * ( 1 /cell ) * (cell * S2_S4_S4 * kdiss_SS),
    D(ppS3_ppS3_S4) ~ +1.0 * ( 1 /cell ) * (cell * S4 * k_334 * (ppS3)^(2))-1.0 * ( 1 /cell ) * (cell * 2 * S_dephosphos * ppS3_ppS3_S4),
    D(geneB) ~ +1.0 * ( 1 /cell ) * (cell * ((geneB_turn + geneB_act2 * ppS2_S4_S4 + geneB_act1 * ppS2_ppS3_S4 + geneB_act3 * ppS2_ppS3_ppS3) / (geneB_inh2 * ppS2_S4_S4 + geneB_inh1 * ppS2_ppS3_S4 + geneB_inh3 * ppS2_ppS3_ppS3 + 1)))-1.0 * ( 1 /cell ) * (cell * geneB * geneB_turn),
    D(ppS3) ~ -3.0 * ( 1 /cell ) * (cell * khomo3 * (ppS3)^(3))+2.0 * ( 1 /cell ) * (cell * 3 * S_dephosphos * ppS3_ppS3_ppS3)+1.0 * ( 1 /cell ) * (cell * S3 * S_phos * TGFb_pRec)-1.0 * ( 1 /cell ) * (cell * S_dephosphos * ppS3)-1.0 * ( 1 /cell ) * (cell * k_223 * (ppS2)^(2) * ppS3)+1.0 * ( 1 /cell ) * (cell * 2 * S_dephosphos * ppS2_ppS2_ppS3)-2.0 * ( 1 /cell ) * (cell * k_233 * ppS2 * (ppS3)^(2))+2.0 * ( 1 /cell ) * (cell * S_dephosphos * ppS2_ppS3_ppS3)+1.0 * ( 1 /cell ) * (cell * 2 * S_dephosphos * ppS2_ppS3_ppS3)-2.0 * ( 1 /cell ) * (cell * S4 * k_334 * (ppS3)^(2))+1.0 * ( 1 /cell ) * (cell * 2 * S_dephosphos * ppS3_ppS3_S4)-1.0 * ( 1 /cell ) * (cell * (S4)^(2) * k_344 * ppS3)-1.0 * ( 1 /cell ) * (cell * S4 * k_234 * ppS2 * ppS3)+1.0 * ( 1 /cell ) * (cell * S_dephosphos * ppS2_ppS3_S4),
    D(geneA) ~ +1.0 * ( 1 /cell ) * (cell * ((geneA_turn + geneA_act2 * ppS2_S4_S4 + geneA_act1 * ppS2_ppS3_S4 + geneA_act3 * ppS2_ppS3_ppS3) / (geneA_inh2 * ppS2_S4_S4 + geneA_inh1 * ppS2_ppS3_S4 + geneA_inh3 * ppS2_ppS3_ppS3 + 1)))-1.0 * ( 1 /cell ) * (cell * geneA * geneA_turn),
    D(geneD) ~ +1.0 * ( 1 /cell ) * (cell * ((geneD_turn + geneD_act2 * ppS2_S4_S4 + geneD_act1 * ppS2_ppS3_S4 + geneD_act3 * ppS2_ppS3_ppS3) / (geneD_inh2 * ppS2_S4_S4 + geneD_inh1 * ppS2_ppS3_S4 + geneD_inh3 * ppS2_ppS3_ppS3 + 1)))-1.0 * ( 1 /cell ) * (cell * geneD * geneD_turn),
    D(ppS2_S4_S4) ~ +1.0 * ( 1 /cell ) * (cell * (S4)^(2) * k_244 * ppS2)-1.0 * ( 1 /cell ) * (cell * S_dephosphos * ppS2_S4_S4),
    D(geneL) ~ +1.0 * ( 1 /cell ) * (cell * ((geneL_turn + geneL_act2 * ppS2_S4_S4 + geneL_act1 * ppS2_ppS3_S4 + geneL_act3 * ppS2_ppS3_ppS3) / (geneL_inh2 * ppS2_S4_S4 + geneL_inh1 * ppS2_ppS3_S4 + geneL_inh3 * ppS2_ppS3_ppS3 + 1)))-1.0 * ( 1 /cell ) * (cell * geneL * geneL_turn),
    D(ppS2_ppS3_S4) ~ +1.0 * ( 1 /cell ) * (cell * S4 * k_234 * ppS2 * ppS3)-1.0 * ( 1 /cell ) * (cell * S_dephosphos * ppS2_ppS3_S4)-1.0 * ( 1 /cell ) * (cell * S_dephosphos * ppS2_ppS3_S4),
    D(geneG) ~ +1.0 * ( 1 /cell ) * (cell * ((geneG_turn + geneG_act2 * ppS2_S4_S4 + geneG_act1 * ppS2_ppS3_S4 + geneG_act3 * ppS2_ppS3_ppS3) / (geneG_inh2 * ppS2_S4_S4 + geneG_inh1 * ppS2_ppS3_S4 + geneG_inh3 * ppS2_ppS3_ppS3 + 1)))-1.0 * ( 1 /cell ) * (cell * geneG * geneG_turn),
    D(Rec) ~ -1.0 * ( 1 /cell ) * (cell * Rec * Rec_act * TGFb),
    D(ppS2_ppS2_S4) ~ +1.0 * ( 1 /cell ) * (cell * S4 * k_224 * (ppS2)^(2))-1.0 * ( 1 /cell ) * (cell * 2 * S_dephosphos * ppS2_ppS2_S4),
    D(ppS2_ppS3_ppS3) ~ +1.0 * ( 1 /cell ) * (cell * k_233 * ppS2 * (ppS3)^(2))-1.0 * ( 1 /cell ) * (cell * S_dephosphos * ppS2_ppS3_ppS3)-1.0 * ( 1 /cell ) * (cell * 2 * S_dephosphos * ppS2_ppS3_ppS3),
    D(ppS2) ~ -3.0 * ( 1 /cell ) * (cell * khomo2 * (ppS2)^(3))+2.0 * ( 1 /cell ) * (cell * 3 * S_dephosphos * ppS2_ppS2_ppS2)+1.0 * ( 1 /cell ) * (cell * S2 * S_phos * TGFb_pRec)-1.0 * ( 1 /cell ) * (cell * S_dephosphos * ppS2)-2.0 * ( 1 /cell ) * (cell * k_223 * (ppS2)^(2) * ppS3)+1.0 * ( 1 /cell ) * (cell * 2 * S_dephosphos * ppS2_ppS2_ppS3)+2.0 * ( 1 /cell ) * (cell * S_dephosphos * ppS2_ppS2_ppS3)-2.0 * ( 1 /cell ) * (cell * S4 * k_224 * (ppS2)^(2))+1.0 * ( 1 /cell ) * (cell * 2 * S_dephosphos * ppS2_ppS2_S4)-1.0 * ( 1 /cell ) * (cell * k_233 * ppS2 * (ppS3)^(2))+1.0 * ( 1 /cell ) * (cell * 2 * S_dephosphos * ppS2_ppS3_ppS3)-1.0 * ( 1 /cell ) * (cell * (S4)^(2) * k_244 * ppS2)-1.0 * ( 1 /cell ) * (cell * S4 * k_234 * ppS2 * ppS3)+1.0 * ( 1 /cell ) * (cell * S_dephosphos * ppS2_ppS3_S4),
    D(dummyVariable) ~ +init_S4+init_Rec+init_TGFb+init_S2+init_S3
    ]

    @named sys = ODESystem(eqs)

    ### Initial species concentrations ###
    initialSpeciesValues = [
    ppS3_ppS3_ppS3 => 0.0,
    ppS3_S4_S4 => 0.0,
    geneH => 1.0,
    geneI => 1.0,
    geneJ => 1.0,
    TGFb_pRec => 0.0,
    geneE => 1.0,
    S4 => init_S4,
    ppS2_ppS2_ppS2 => 0.0,
    geneK => 1.0,
    pS2 => 0.0,
    pS3 => 0.0,
    geneC => 1.0,
    ppS2_ppS2_ppS3 => 0.0,
    geneF => 1.0,
    S4_S4_S4 => 0.0,
    TGFb => init_TGFb,
    S3 => init_S3,
    S2 => init_S2,
    S2_S4_S4 => 0.0,
    ppS3_ppS3_S4 => 0.0,
    geneB => 1.0,
    ppS3 => 0.0,
    geneA => 1.0,
    geneD => 1.0,
    ppS2_S4_S4 => 0.0,
    geneL => 1.0,
    ppS2_ppS3_S4 => 0.0,
    geneG => 1.0,
    Rec => init_Rec,
    ppS2_ppS2_S4 => 0.0,
    ppS2_ppS3_ppS3 => 0.0,
    ppS2 => 0.0,
    dummyVariable => 0.0]

    ### True parameter values ###
    trueParameterValues = [
    geneC_inh3 => 0.035379693264711,
    geneJ_act1 => 0.0,
    geneG_inh1 => 0.0,
    init_Rec => 1.83614829885072,
    geneC_act3 => 0.0,
    geneA_inh2 => 0.0462523981984203,
    geneK_turn => 0.00359915570352007,
    geneA_act3 => 0.0,
    geneI_act2 => 0.0,
    geneA_act2 => 0.000797280887701997,
    geneD_inh2 => 0.0868803385899104,
    k_223 => 0.0,
    k_233 => 0.148655616115831,
    geneB_act1 => 0.461699170803527,
    S_dephosphos => 0.0489949146416608,
    geneG_act1 => 626.768522107231,
    geneH_inh2 => 19.0317907364409,
    geneI_act1 => 9.30152789301652,
    geneB_inh3 => 0.407329851468312,
    geneI_inh2 => 320.763519415643,
    geneC_turn => 0.00378630404521753,
    geneJ_inh1 => 0.47739150796461,
    geneK_act3 => 0.0,
    geneJ_inh3 => 0.00851371884219428,
    geneG_inh3 => 0.000815621975867545,
    geneC_act2 => 0.0,
    init_TGFb => 0.0,
    geneL_inh1 => 1.29976503414042,
    geneE_inh3 => 9.67035015330535,
    geneF_act3 => 20.6442242166861,
    geneF_inh1 => 3.63157061548523,
    geneD_act3 => 0.0,
    geneE_act2 => 1.03661055568903,
    geneD_inh3 => 0.101419105775767,
    k_234 => 0.000403100938582095,
    geneH_act2 => 89.9866382703511,
    geneA_turn => 0.00455788978022206,
    geneL_act1 => 0.0,
    geneH_act1 => 0.0,
    geneL_inh2 => 0.062323573525662,
    geneB_turn => 0.0112459404708192,
    init_S4 => 0.0,
    geneL_inh3 => 0.13414856027346,
    khomo2 => 0.0,
    geneK_act2 => 0.000284224503418874,
    geneA_inh3 => 0.026381058989095,
    geneJ_act2 => 0.0,
    geneI_inh3 => 999.981624754706,
    geneD_act1 => 0.0111673720417551,
    geneJ_turn => 0.0140190729281921,
    geneG_act3 => 0.0,
    geneL_act2 => 0.000474849473374749,
    Rec_act => 0.00598343088790046,
    geneI_act3 => 0.0,
    k_224 => 0.0,
    pRec_degind => 0.0404498417896043,
    geneE_inh2 => 1.4392548717489,
    geneF_inh2 => 0.0,
    geneC_act1 => 0.00714636599139809,
    geneD_inh1 => 0.0,
    k_344 => 0.0,
    init_S3 => 0.0,
    geneI_inh1 => 0.0,
    geneL_turn => 0.0149504509077553,
    geneF_act2 => 0.135655043663672,
    k_on_u => 0.0,
    geneE_inh1 => 8.20035989079744,
    geneH_act3 => 999.921986303898,
    geneH_inh1 => 999.998960700206,
    geneK_inh3 => 2.43594567888689,
    geneE_act1 => 0.0,
    geneA_inh1 => 0.0,
    geneB_inh1 => 0.0,
    geneH_inh3 => 218.356690800299,
    geneD_turn => 0.000802160834341515,
    S_dephos => 0.286397320264704,
    geneL_act3 => 0.0,
    S_phos => 0.379720142751521,
    geneG_turn => 295.971427530835,
    geneA_act1 => 0.0141389661164018,
    geneI_turn => 0.00328197841814657,
    khomo3 => 0.0,
    geneC_inh1 => 0.0,
    geneG_act2 => 55.9534537604082,
    init_S2 => 0.0,
    geneK_inh2 => 1.36603583936842,
    k_334 => 0.0,
    geneB_inh2 => 0.747200749545747,
    geneH_turn => 0.992827361229694,
    geneJ_inh2 => 0.00985305306597475,
    khomo4 => 0.0,
    geneB_act2 => 0.073852775177442,
    geneF_turn => 36.6724104315148,
    geneJ_act3 => 0.0,
    geneK_inh1 => 0.0,
    geneB_act3 => 0.0304255473876269,
    geneK_act1 => 0.0125133830064198,
    geneE_turn => 0.124380538819967,
    geneE_act3 => 6.10224669858575,
    geneG_inh2 => 0.0175854489108708,
    kdiss_SS => 0.0,
    k_244 => 8.56307475687022e-6,
    geneD_act2 => 0.000268942132758932,
    geneF_inh3 => 0.0,
    geneF_act1 => 998.299919973209,
    geneC_inh2 => 0.018441582221635]

    return sys, initialSpeciesValues, trueParameterValues

end
