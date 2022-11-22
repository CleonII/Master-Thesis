function evalYmod(u, t, dynPar, obsPar, nonDynParam, paramData, observableId, mapObsParam) 

	ppS3_ppS3_ppS3, ppS3_S4_S4, geneH, geneI, geneJ, TGFb_pRec, geneE, S4, ppS2_ppS2_ppS2, geneK, pS2, pS3, geneC, ppS2_ppS2_ppS3, geneF, S4_S4_S4, TGFb, S3, S2, S2_S4_S4, ppS3_ppS3_S4, geneB, ppS3, geneA, geneD, ppS2_S4_S4, geneL, ppS2_ppS3_S4, geneG, Rec, ppS2_ppS2_S4, ppS2_ppS3_ppS3, ppS2= u 
	Rec_act, S2tot, S3tot, S4tot, S_dephos, S_dephosphos, S_phos, geneA_act1, geneA_act2, geneA_inh2, geneA_inh3, geneA_turn, geneB_act1, geneB_act2, geneB_act3, geneB_inh2, geneB_inh3, geneB_turn, geneC_act1, geneC_inh2, geneC_inh3, geneC_turn, geneD_act1, geneD_act2, geneD_inh2, geneD_inh3, geneD_turn, geneE_act2, geneE_act3, geneE_inh1, geneE_inh2, geneE_inh3, geneE_turn, geneF_act1, geneF_act2, geneF_act3, geneF_inh1, geneF_turn, geneG_act1, geneG_act2, geneG_inh2, geneG_inh3, geneG_turn, geneH_act2, geneH_act3, geneH_inh1, geneH_inh2, geneH_inh3, geneH_turn, geneI_act1, geneI_inh2, geneI_inh3, geneI_turn, geneJ_inh1, geneJ_inh2, geneJ_inh3, geneJ_turn, geneK_act1, geneK_act2, geneK_inh2, geneK_inh3, geneK_turn, geneL_act2, geneL_inh1, geneL_inh2, geneL_inh3, geneL_turn, init_Rec, k_233, k_234, k_244, pRec_degind = dynPar 
	sd_Bmp4_nExpID100, sd_Cxcl15_nExpID100, sd_Dnmt3a_nExpID100, sd_Dusp5_nExpID100, sd_Jun_nExpID100, sd_Klf10_nExpID100, sd_Pdk4_nExpID100, sd_Ski_nExpID100, sd_Skil_nExpID100, sd_Smad7_nExpID100, sd_Sox4_nExpID100, sd_Tgfa_nExpID100 = nonDynParam 
	geneA_act3_C = paramData.paramVal[10] 
	geneA_inh1_C = paramData.paramVal[11] 
	geneB_inh1_C = paramData.paramVal[18] 
	geneC_act2_C = paramData.paramVal[23] 
	geneC_act3_C = paramData.paramVal[24] 
	geneC_inh1_C = paramData.paramVal[25] 
	geneD_act3_C = paramData.paramVal[31] 
	geneD_inh1_C = paramData.paramVal[32] 
	geneE_act1_C = paramData.paramVal[36] 
	geneF_inh2_C = paramData.paramVal[47] 
	geneF_inh3_C = paramData.paramVal[48] 
	geneG_act3_C = paramData.paramVal[52] 
	geneG_inh1_C = paramData.paramVal[53] 
	geneH_act1_C = paramData.paramVal[57] 
	geneI_act2_C = paramData.paramVal[65] 
	geneI_act3_C = paramData.paramVal[66] 
	geneI_inh1_C = paramData.paramVal[67] 
	geneJ_act1_C = paramData.paramVal[71] 
	geneJ_act2_C = paramData.paramVal[72] 
	geneJ_act3_C = paramData.paramVal[73] 
	geneK_act3_C = paramData.paramVal[80] 
	geneK_inh1_C = paramData.paramVal[81] 
	geneL_act1_C = paramData.paramVal[85] 
	geneL_act3_C = paramData.paramVal[87] 
	k_223_C = paramData.paramVal[93] 
	k_224_C = paramData.paramVal[94] 
	k_334_C = paramData.paramVal[98] 
	k_344_C = paramData.paramVal[99] 
	k_on_u_C = paramData.paramVal[100] 
	kdiss_SS_C = paramData.paramVal[101] 
	khomo2_C = paramData.paramVal[102] 
	khomo3_C = paramData.paramVal[103] 
	khomo4_C = paramData.paramVal[104] 

	if observableId == "Bmp4_gene_expression_OE_nExpID100" 
		return geneH + 0.000001 
	end

	if observableId == "Cxcl15_gene_expression_OE_nExpID100" 
		return geneI + 0.000001 
	end

	if observableId == "Dnmt3a_gene_expression_OE_nExpID100" 
		return geneC + 0.000001 
	end

	if observableId == "Dusp5_gene_expression_OE_nExpID100" 
		return geneJ + 0.000001 
	end

	if observableId == "Jun_gene_expression_OE_nExpID100" 
		return geneE + 0.000001 
	end

	if observableId == "Klf10_gene_expression_OE_nExpID100" 
		return geneG + 0.000001 
	end

	if observableId == "Pdk4_gene_expression_OE_nExpID100" 
		return geneL + 0.000001 
	end

	if observableId == "S23IP_S2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5" 
		return ( 100 * ( S2 + S2_S4_S4 ) ) / ( S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 ) 
	end

	if observableId == "S23IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return ( 100 * ( S2 + S2_S4_S4 ) ) / ( S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 ) 
	end

	if observableId == "S23IP_S3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5" 
		return ( 100 * S3 ) / ( S3 + pS3 + ppS3 + ppS3_S4_S4 + ppS2_ppS3_S4 + 2 * ppS3_ppS3_S4 + ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 3 * ppS3_ppS3_ppS3 ) 
	end

	if observableId == "S23IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return ( 100 * S3 ) / ( S3 + pS3 + ppS3 + ppS3_S4_S4 + ppS2_ppS3_S4 + 2 * ppS3_ppS3_S4 + ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 3 * ppS3_ppS3_ppS3 ) 
	end

	if observableId == "S23IP_pS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5" 
		return ( 100 * pS2 ) / ( S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 ) 
	end

	if observableId == "S23IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return ( 100 * pS2 ) / ( S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 ) 
	end

	if observableId == "S23IP_pS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5" 
		return ( 100 * pS3 ) / ( S3 + pS3 + ppS3 + ppS3_S4_S4 + ppS2_ppS3_S4 + 2 * ppS3_ppS3_S4 + ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 3 * ppS3_ppS3_ppS3 ) 
	end

	if observableId == "S23IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return ( 100 * pS3 ) / ( S3 + pS3 + ppS3 + ppS3_S4_S4 + ppS2_ppS3_S4 + 2 * ppS3_ppS3_S4 + ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 3 * ppS3_ppS3_ppS3 ) 
	end

	if observableId == "S23IP_ppS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5" 
		return ( 100 * ( ppS2 + ppS2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 ) ) / ( S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 ) 
	end

	if observableId == "S23IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return ( 100 * ( ppS2 + ppS2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 ) ) / ( S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 ) 
	end

	if observableId == "S23IP_ppS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5" 
		return ( 100 * ( ppS3 + ppS3_S4_S4 + ppS2_ppS3_S4 + 2 * ppS3_ppS3_S4 + ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 3 * ppS3_ppS3_ppS3 ) ) / ( S3 + pS3 + ppS3 + ppS3_S4_S4 + ppS2_ppS3_S4 + 2 * ppS3_ppS3_S4 + ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 3 * ppS3_ppS3_ppS3 ) 
	end

	if observableId == "S23IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return ( 100 * ( ppS3 + ppS3_S4_S4 + ppS2_ppS3_S4 + 2 * ppS3_ppS3_S4 + ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 3 * ppS3_ppS3_ppS3 ) ) / ( S3 + pS3 + ppS3 + ppS3_S4_S4 + ppS2_ppS3_S4 + 2 * ppS3_ppS3_S4 + ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 3 * ppS3_ppS3_ppS3 ) 
	end

	if observableId == "S23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1" 
		return S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 
	end

	if observableId == "S23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10" 
		return ( 100 * ( S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 ) ) / ( S2 + S3 + S2_S4_S4 + pS2 + pS3 + ppS2 + ppS3 + ppS2_S4_S4 + ppS3_S4_S4 + 2 * ppS2_ppS2_S4 + 2 * ppS2_ppS3_S4 + 2 * ppS3_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 3 * ppS2_ppS2_ppS3 + 3 * ppS2_ppS3_ppS3 + 3 * ppS3_ppS3_ppS3 ) 
	end

	if observableId == "S23IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return ( 100 * ( S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 ) ) / ( S2 + S3 + 3 * S2_S4_S4 + pS2 + pS3 + ppS2 + ppS3 + 3 * ppS2_S4_S4 + 3 * ppS3_S4_S4 + 3 * ppS2_ppS2_S4 + 4 * ppS2_ppS3_S4 + 3 * ppS3_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 3 * ppS2_ppS2_ppS3 + 3 * ppS2_ppS3_ppS3 + 3 * ppS3_ppS3_ppS3 ) 
	end

	if observableId == "S23IP_tS3_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10" 
		return ( 100 * ( S3 + pS3 + ppS3 + ppS3_S4_S4 + ppS2_ppS3_S4 + 2 * ppS3_ppS3_S4 + ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 3 * ppS3_ppS3_ppS3 ) ) / ( S2 + S3 + S2_S4_S4 + pS2 + pS3 + ppS2 + ppS3 + ppS2_S4_S4 + ppS3_S4_S4 + 2 * ppS2_ppS2_S4 + 2 * ppS2_ppS3_S4 + 2 * ppS3_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 3 * ppS2_ppS2_ppS3 + 3 * ppS2_ppS3_ppS3 + 3 * ppS3_ppS3_ppS3 ) 
	end

	if observableId == "S23IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return ( 100 * ( S3 + pS3 + ppS3 + ppS3_S4_S4 + ppS2_ppS3_S4 + 2 * ppS3_ppS3_S4 + ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 3 * ppS3_ppS3_ppS3 ) ) / ( S2 + S3 + 3 * S2_S4_S4 + pS2 + pS3 + ppS2 + ppS3 + 3 * ppS2_S4_S4 + 3 * ppS3_S4_S4 + 3 * ppS2_ppS2_S4 + 4 * ppS2_ppS3_S4 + 3 * ppS3_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 3 * ppS2_ppS2_ppS3 + 3 * ppS2_ppS3_ppS3 + 3 * ppS3_ppS3_ppS3 ) 
	end

	if observableId == "S23IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return ( 100 * ( 2 * S2_S4_S4 + 2 * ppS2_S4_S4 + 2 * ppS3_S4_S4 + ppS2_ppS2_S4 + 2 * ppS2_ppS3_S4 + ppS3_ppS3_S4 ) ) / ( S2 + S3 + 3 * S2_S4_S4 + pS2 + pS3 + ppS2 + ppS3 + 3 * ppS2_S4_S4 + 3 * ppS3_S4_S4 + 3 * ppS2_ppS2_S4 + 4 * ppS2_ppS3_S4 + 3 * ppS3_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 3 * ppS2_ppS2_ppS3 + 3 * ppS2_ppS3_ppS3 + 3 * ppS3_ppS3_ppS3 ) 
	end

	if observableId == "S2IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return ( 100 * ( S2 + S2_S4_S4 ) ) / ( S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 ) 
	end

	if observableId == "S2IP_S2_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return ( 100 * ( S2 + S2_S4_S4 ) ) / ( S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 ) 
	end

	if observableId == "S2IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 0 
	end

	if observableId == "S2IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return ( 100 * pS2 ) / ( S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 ) 
	end

	if observableId == "S2IP_pS2_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return ( 100 * pS2 ) / ( S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 ) 
	end

	if observableId == "S2IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 0 
	end

	if observableId == "S2IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return ( 100 * ( ppS2 + ppS2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 ) ) / ( S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 ) 
	end

	if observableId == "S2IP_ppS2_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return ( 100 * ( ppS2 + ppS2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 ) ) / ( S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 ) 
	end

	if observableId == "S2IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 100 
	end

	if observableId == "S2IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return ( 100 * ( S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 ) ) / ( S2 + 3 * S2_S4_S4 + pS2 + ppS2 + 3 * ppS2_S4_S4 + 3 * ppS2_ppS2_S4 + 3 * ppS2_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 3 * ppS2_ppS2_ppS3 + 3 * ppS2_ppS3_ppS3 ) 
	end

	if observableId == "S2IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return ( 100 * ( S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 ) ) / ( S2 + 3 * S2_S4_S4 + pS2 + ppS2 + 3 * ppS2_S4_S4 + 3 * ppS2_ppS2_S4 + 3 * ppS2_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 3 * ppS2_ppS2_ppS3 + 3 * ppS2_ppS3_ppS3 ) 
	end

	if observableId == "S2IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return ( 100 * ( ppS2_ppS3_S4 + ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 ) ) / ( S2 + 3 * S2_S4_S4 + pS2 + ppS2 + 3 * ppS2_S4_S4 + 3 * ppS2_ppS2_S4 + 3 * ppS2_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 3 * ppS2_ppS2_ppS3 + 3 * ppS2_ppS3_ppS3 ) 
	end

	if observableId == "S2IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return ( 100 * ( ppS2_ppS3_S4 + ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 ) ) / ( S2 + 3 * S2_S4_S4 + pS2 + ppS2 + 3 * ppS2_S4_S4 + 3 * ppS2_ppS2_S4 + 3 * ppS2_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 3 * ppS2_ppS2_ppS3 + 3 * ppS2_ppS3_ppS3 ) 
	end

	if observableId == "S2IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return ( 100 * ( 2 * ppS2_S4_S4 + ppS2_ppS2_S4 + ppS2_ppS3_S4 ) ) / ( S2 + 3 * S2_S4_S4 + pS2 + ppS2 + 3 * ppS2_S4_S4 + 3 * ppS2_ppS2_S4 + 3 * ppS2_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 3 * ppS2_ppS2_ppS3 + 3 * ppS2_ppS3_ppS3 ) 
	end

	if observableId == "S2IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return ( 100 * ( 2 * ppS2_S4_S4 + ppS2_ppS2_S4 + ppS2_ppS3_S4 ) ) / ( S2 + 3 * S2_S4_S4 + pS2 + ppS2 + 3 * ppS2_S4_S4 + 3 * ppS2_ppS2_S4 + 3 * ppS2_ppS3_S4 + 3 * ppS2_ppS2_ppS2 + 3 * ppS2_ppS2_ppS3 + 3 * ppS2_ppS3_ppS3 ) 
	end

	if observableId == "S3IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 0 
	end

	if observableId == "S3IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return ( 100 * S3 ) / ( S3 + pS3 + ppS3 + ppS3_S4_S4 + ppS2_ppS3_S4 + 2 * ppS3_ppS3_S4 + ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 3 * ppS3_ppS3_ppS3 ) 
	end

	if observableId == "S3IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 0 
	end

	if observableId == "S3IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return ( 100 * pS3 ) / ( S3 + pS3 + ppS3 + ppS3_S4_S4 + ppS2_ppS3_S4 + 2 * ppS3_ppS3_S4 + ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 3 * ppS3_ppS3_ppS3 ) 
	end

	if observableId == "S3IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 100 
	end

	if observableId == "S3IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return ( 100 * ( ppS3 + ppS3_S4_S4 + ppS2_ppS3_S4 + 2 * ppS3_ppS3_S4 + ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 3 * ppS3_ppS3_ppS3 ) ) / ( S3 + pS3 + ppS3 + ppS3_S4_S4 + ppS2_ppS3_S4 + 2 * ppS3_ppS3_S4 + ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 3 * ppS3_ppS3_ppS3 ) 
	end

	if observableId == "S3IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return ( 100 * ( ppS2_ppS3_S4 + 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 ) ) / ( S3 + pS3 + ppS3 + 3 * ppS3_S4_S4 + 3 * ppS2_ppS3_S4 + 3 * ppS3_ppS3_S4 + 3 * ppS2_ppS2_ppS3 + 3 * ppS2_ppS3_ppS3 + 3 * ppS3_ppS3_ppS3 ) 
	end

	if observableId == "S3IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return ( 100 * ( ppS2_ppS3_S4 + 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 ) ) / ( S3 + pS3 + ppS3 + 3 * ppS3_S4_S4 + 3 * ppS2_ppS3_S4 + 3 * ppS3_ppS3_S4 + 3 * ppS2_ppS2_ppS3 + 3 * ppS2_ppS3_ppS3 + 3 * ppS3_ppS3_ppS3 ) 
	end

	if observableId == "S3IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return ( 100 * ( S3 + pS3 + ppS3 + ppS3_S4_S4 + ppS2_ppS3_S4 + 2 * ppS3_ppS3_S4 + ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 3 * ppS3_ppS3_ppS3 ) ) / ( S3 + pS3 + ppS3 + 3 * ppS3_S4_S4 + 3 * ppS2_ppS3_S4 + 3 * ppS3_ppS3_S4 + 3 * ppS2_ppS2_ppS3 + 3 * ppS2_ppS3_ppS3 + 3 * ppS3_ppS3_ppS3 ) 
	end

	if observableId == "S3IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return ( 100 * ( S3 + pS3 + ppS3 + ppS3_S4_S4 + ppS2_ppS3_S4 + 2 * ppS3_ppS3_S4 + ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 3 * ppS3_ppS3_ppS3 ) ) / ( S3 + pS3 + ppS3 + 3 * ppS3_S4_S4 + 3 * ppS2_ppS3_S4 + 3 * ppS3_ppS3_S4 + 3 * ppS2_ppS2_ppS3 + 3 * ppS2_ppS3_ppS3 + 3 * ppS3_ppS3_ppS3 ) 
	end

	if observableId == "S3IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return ( 100 * ( 2 * ppS3_S4_S4 + ppS2_ppS3_S4 + ppS3_ppS3_S4 ) ) / ( S3 + pS3 + ppS3 + 3 * ppS3_S4_S4 + 3 * ppS2_ppS3_S4 + 3 * ppS3_ppS3_S4 + 3 * ppS2_ppS2_ppS3 + 3 * ppS2_ppS3_ppS3 + 3 * ppS3_ppS3_ppS3 ) 
	end

	if observableId == "S3IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return ( 100 * ( 2 * ppS3_S4_S4 + ppS2_ppS3_S4 + ppS3_ppS3_S4 ) ) / ( S3 + pS3 + ppS3 + 3 * ppS3_S4_S4 + 3 * ppS2_ppS3_S4 + 3 * ppS3_ppS3_S4 + 3 * ppS2_ppS2_ppS3 + 3 * ppS2_ppS3_ppS3 + 3 * ppS3_ppS3_ppS3 ) 
	end

	if observableId == "S4IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return ( 100 * ( S2_S4_S4 + ppS2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_S4 ) ) / ( S4 + 3 * S2_S4_S4 + 3 * S4_S4_S4 + 3 * ppS2_S4_S4 + 3 * ppS3_S4_S4 + 3 * ppS2_ppS2_S4 + 3 * ppS2_ppS3_S4 + 3 * ppS3_ppS3_S4 ) 
	end

	if observableId == "S4IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return ( 100 * ( S2_S4_S4 + ppS2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_S4 ) ) / ( S4 + 3 * S2_S4_S4 + 3 * S4_S4_S4 + 3 * ppS2_S4_S4 + 3 * ppS3_S4_S4 + 3 * ppS2_ppS2_S4 + 3 * ppS2_ppS3_S4 + 3 * ppS3_ppS3_S4 ) 
	end

	if observableId == "S4IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return ( 100 * ( ppS3_S4_S4 + ppS2_ppS3_S4 + 2 * ppS3_ppS3_S4 ) ) / ( S4 + 3 * S2_S4_S4 + 3 * S4_S4_S4 + 3 * ppS2_S4_S4 + 3 * ppS3_S4_S4 + 3 * ppS2_ppS2_S4 + 3 * ppS2_ppS3_S4 + 3 * ppS3_ppS3_S4 ) 
	end

	if observableId == "S4IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return ( 100 * ( ppS3_S4_S4 + ppS2_ppS3_S4 + 2 * ppS3_ppS3_S4 ) ) / ( S4 + 3 * S2_S4_S4 + 3 * S4_S4_S4 + 3 * ppS2_S4_S4 + 3 * ppS3_S4_S4 + 3 * ppS2_ppS2_S4 + 3 * ppS2_ppS3_S4 + 3 * ppS3_ppS3_S4 ) 
	end

	if observableId == "S4IP_tS4_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1" 
		return S4 + 2 * S2_S4_S4 + 3 * S4_S4_S4 + 2 * ppS2_S4_S4 + 2 * ppS3_S4_S4 + ppS2_ppS2_S4 + ppS2_ppS3_S4 + ppS3_ppS3_S4 
	end

	if observableId == "S4IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return ( 100 * ( S4 + 2 * S2_S4_S4 + 3 * S4_S4_S4 + 2 * ppS2_S4_S4 + 2 * ppS3_S4_S4 + ppS2_ppS2_S4 + ppS2_ppS3_S4 + ppS3_ppS3_S4 ) ) / ( S4 + 3 * S2_S4_S4 + 3 * S4_S4_S4 + 3 * ppS2_S4_S4 + 3 * ppS3_S4_S4 + 3 * ppS2_ppS2_S4 + 3 * ppS2_ppS3_S4 + 3 * ppS3_ppS3_S4 ) 
	end

	if observableId == "S4IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return ( 100 * ( S4 + 2 * S2_S4_S4 + 3 * S4_S4_S4 + 2 * ppS2_S4_S4 + 2 * ppS3_S4_S4 + ppS2_ppS2_S4 + ppS2_ppS3_S4 + ppS3_ppS3_S4 ) ) / ( S4 + 3 * S2_S4_S4 + 3 * S4_S4_S4 + 3 * ppS2_S4_S4 + 3 * ppS3_S4_S4 + 3 * ppS2_ppS2_S4 + 3 * ppS2_ppS3_S4 + 3 * ppS3_ppS3_S4 ) 
	end

	if observableId == "Ski_gene_expression_OE_nExpID100" 
		return geneA + 0.000001 
	end

	if observableId == "Skil_gene_expression_OE_nExpID100" 
		return geneB + 0.000001 
	end

	if observableId == "Smad7_gene_expression_OE_nExpID100" 
		return geneF + 0.000001 
	end

	if observableId == "Sox4_gene_expression_OE_nExpID100" 
		return geneD + 0.000001 
	end

	if observableId == "TGFbR_PL_PN_TGFbR_hepa16_nExpID11" 
		return Rec + TGFb_pRec 
	end

	if observableId == "Tgfa_gene_expression_OE_nExpID100" 
		return geneK + 0.000001 
	end

end

function evalU0!(u0Vec, paramVec) 

	geneC_inh3, geneJ_act1, geneG_inh1, init_Rec, geneC_act3, geneA_inh2, geneK_turn, geneA_act3, geneI_act2, geneA_act2, geneD_inh2, k_223, k_233, geneB_act1, S_dephosphos, geneG_act1, geneH_inh2, geneI_act1, geneB_inh3, geneI_inh2, geneC_turn, geneJ_inh1, geneK_act3, geneJ_inh3, geneG_inh3, geneC_act2, init_TGFb, geneL_inh1, geneE_inh3, geneF_act3, cell, geneF_inh1, geneD_act3, geneE_act2, geneD_inh3, k_234, geneH_act2, geneA_turn, geneL_act1, geneH_act1, geneL_inh2, geneB_turn, init_S4, geneL_inh3, khomo2, geneK_act2, geneA_inh3, geneJ_act2, geneI_inh3, geneD_act1, geneJ_turn, geneG_act3, geneL_act2, Rec_act, geneI_act3, k_224, pRec_degind, geneE_inh2, geneF_inh2, geneC_act1, geneD_inh1, k_344, init_S3, geneI_inh1, geneL_turn, geneF_act2, k_on_u, geneE_inh1, geneH_act3, geneH_inh1, geneK_inh3, geneE_act1, geneA_inh1, geneB_inh1, geneH_inh3, geneD_turn, S_dephos, geneL_act3, S_phos, geneG_turn, geneA_act1, geneI_turn, khomo3, geneC_inh1, geneG_act2, init_S2, geneK_inh2, k_334, geneB_inh2, geneH_turn, geneJ_inh2, khomo4, geneB_act2, geneF_turn, geneJ_act3, geneK_inh1, geneB_act3, geneK_act1, geneE_turn, geneE_act3, geneG_inh2, kdiss_SS, k_244, geneD_act2, geneF_inh3, geneF_act1, geneC_inh2 = paramVec 

	ppS3_ppS3_ppS3 = 0.0 
	ppS3_S4_S4 = 0.0 
	geneH = 1.0 
	geneI = 1.0 
	geneJ = 1.0 
	TGFb_pRec = 0.0 
	geneE = 1.0 
	S4 = init_S4 
	ppS2_ppS2_ppS2 = 0.0 
	geneK = 1.0 
	pS2 = 0.0 
	pS3 = 0.0 
	geneC = 1.0 
	ppS2_ppS2_ppS3 = 0.0 
	geneF = 1.0 
	S4_S4_S4 = 0.0 
	TGFb = init_TGFb 
	S3 = init_S3 
	S2 = init_S2 
	S2_S4_S4 = 0.0 
	ppS3_ppS3_S4 = 0.0 
	geneB = 1.0 
	ppS3 = 0.0 
	geneA = 1.0 
	geneD = 1.0 
	ppS2_S4_S4 = 0.0 
	geneL = 1.0 
	ppS2_ppS3_S4 = 0.0 
	geneG = 1.0 
	Rec = init_Rec 
	ppS2_ppS2_S4 = 0.0 
	ppS2_ppS3_ppS3 = 0.0 
	ppS2 = 0.0 

	u0Vec .= ppS3_ppS3_ppS3, ppS3_S4_S4, geneH, geneI, geneJ, TGFb_pRec, geneE, S4, ppS2_ppS2_ppS2, geneK, pS2, pS3, geneC, ppS2_ppS2_ppS3, geneF, S4_S4_S4, TGFb, S3, S2, S2_S4_S4, ppS3_ppS3_S4, geneB, ppS3, geneA, geneD, ppS2_S4_S4, geneL, ppS2_ppS3_S4, geneG, Rec, ppS2_ppS2_S4, ppS2_ppS3_ppS3, ppS2
end

function evalU0(paramVec) 

	geneC_inh3, geneJ_act1, geneG_inh1, init_Rec, geneC_act3, geneA_inh2, geneK_turn, geneA_act3, geneI_act2, geneA_act2, geneD_inh2, k_223, k_233, geneB_act1, S_dephosphos, geneG_act1, geneH_inh2, geneI_act1, geneB_inh3, geneI_inh2, geneC_turn, geneJ_inh1, geneK_act3, geneJ_inh3, geneG_inh3, geneC_act2, init_TGFb, geneL_inh1, geneE_inh3, geneF_act3, cell, geneF_inh1, geneD_act3, geneE_act2, geneD_inh3, k_234, geneH_act2, geneA_turn, geneL_act1, geneH_act1, geneL_inh2, geneB_turn, init_S4, geneL_inh3, khomo2, geneK_act2, geneA_inh3, geneJ_act2, geneI_inh3, geneD_act1, geneJ_turn, geneG_act3, geneL_act2, Rec_act, geneI_act3, k_224, pRec_degind, geneE_inh2, geneF_inh2, geneC_act1, geneD_inh1, k_344, init_S3, geneI_inh1, geneL_turn, geneF_act2, k_on_u, geneE_inh1, geneH_act3, geneH_inh1, geneK_inh3, geneE_act1, geneA_inh1, geneB_inh1, geneH_inh3, geneD_turn, S_dephos, geneL_act3, S_phos, geneG_turn, geneA_act1, geneI_turn, khomo3, geneC_inh1, geneG_act2, init_S2, geneK_inh2, k_334, geneB_inh2, geneH_turn, geneJ_inh2, khomo4, geneB_act2, geneF_turn, geneJ_act3, geneK_inh1, geneB_act3, geneK_act1, geneE_turn, geneE_act3, geneG_inh2, kdiss_SS, k_244, geneD_act2, geneF_inh3, geneF_act1, geneC_inh2 = paramVec 

	ppS3_ppS3_ppS3 = 0.0 
	ppS3_S4_S4 = 0.0 
	geneH = 1.0 
	geneI = 1.0 
	geneJ = 1.0 
	TGFb_pRec = 0.0 
	geneE = 1.0 
	S4 = init_S4 
	ppS2_ppS2_ppS2 = 0.0 
	geneK = 1.0 
	pS2 = 0.0 
	pS3 = 0.0 
	geneC = 1.0 
	ppS2_ppS2_ppS3 = 0.0 
	geneF = 1.0 
	S4_S4_S4 = 0.0 
	TGFb = init_TGFb 
	S3 = init_S3 
	S2 = init_S2 
	S2_S4_S4 = 0.0 
	ppS3_ppS3_S4 = 0.0 
	geneB = 1.0 
	ppS3 = 0.0 
	geneA = 1.0 
	geneD = 1.0 
	ppS2_S4_S4 = 0.0 
	geneL = 1.0 
	ppS2_ppS3_S4 = 0.0 
	geneG = 1.0 
	Rec = init_Rec 
	ppS2_ppS2_S4 = 0.0 
	ppS2_ppS3_ppS3 = 0.0 
	ppS2 = 0.0 

	 return [ppS3_ppS3_ppS3, ppS3_S4_S4, geneH, geneI, geneJ, TGFb_pRec, geneE, S4, ppS2_ppS2_ppS2, geneK, pS2, pS3, geneC, ppS2_ppS2_ppS3, geneF, S4_S4_S4, TGFb, S3, S2, S2_S4_S4, ppS3_ppS3_S4, geneB, ppS3, geneA, geneD, ppS2_S4_S4, geneL, ppS2_ppS3_S4, geneG, Rec, ppS2_ppS2_S4, ppS2_ppS3_ppS3, ppS2]
end

function evalSd!(u, t, sdPar, dynPar, nonDynParam, paramData, observableId, mapSdParam) 

	ppS3_ppS3_ppS3, ppS3_S4_S4, geneH, geneI, geneJ, TGFb_pRec, geneE, S4, ppS2_ppS2_ppS2, geneK, pS2, pS3, geneC, ppS2_ppS2_ppS3, geneF, S4_S4_S4, TGFb, S3, S2, S2_S4_S4, ppS3_ppS3_S4, geneB, ppS3, geneA, geneD, ppS2_S4_S4, geneL, ppS2_ppS3_S4, geneG, Rec, ppS2_ppS2_S4, ppS2_ppS3_ppS3, ppS2= u 
	Rec_act, S2tot, S3tot, S4tot, S_dephos, S_dephosphos, S_phos, geneA_act1, geneA_act2, geneA_inh2, geneA_inh3, geneA_turn, geneB_act1, geneB_act2, geneB_act3, geneB_inh2, geneB_inh3, geneB_turn, geneC_act1, geneC_inh2, geneC_inh3, geneC_turn, geneD_act1, geneD_act2, geneD_inh2, geneD_inh3, geneD_turn, geneE_act2, geneE_act3, geneE_inh1, geneE_inh2, geneE_inh3, geneE_turn, geneF_act1, geneF_act2, geneF_act3, geneF_inh1, geneF_turn, geneG_act1, geneG_act2, geneG_inh2, geneG_inh3, geneG_turn, geneH_act2, geneH_act3, geneH_inh1, geneH_inh2, geneH_inh3, geneH_turn, geneI_act1, geneI_inh2, geneI_inh3, geneI_turn, geneJ_inh1, geneJ_inh2, geneJ_inh3, geneJ_turn, geneK_act1, geneK_act2, geneK_inh2, geneK_inh3, geneK_turn, geneL_act2, geneL_inh1, geneL_inh2, geneL_inh3, geneL_turn, init_Rec, k_233, k_234, k_244, pRec_degind = dynPar 
	sd_Bmp4_nExpID100, sd_Cxcl15_nExpID100, sd_Dnmt3a_nExpID100, sd_Dusp5_nExpID100, sd_Jun_nExpID100, sd_Klf10_nExpID100, sd_Pdk4_nExpID100, sd_Ski_nExpID100, sd_Skil_nExpID100, sd_Smad7_nExpID100, sd_Sox4_nExpID100, sd_Tgfa_nExpID100 = nonDynParam 
	geneA_act3_C = paramData.paramVal[10] 
	geneA_inh1_C = paramData.paramVal[11] 
	geneB_inh1_C = paramData.paramVal[18] 
	geneC_act2_C = paramData.paramVal[23] 
	geneC_act3_C = paramData.paramVal[24] 
	geneC_inh1_C = paramData.paramVal[25] 
	geneD_act3_C = paramData.paramVal[31] 
	geneD_inh1_C = paramData.paramVal[32] 
	geneE_act1_C = paramData.paramVal[36] 
	geneF_inh2_C = paramData.paramVal[47] 
	geneF_inh3_C = paramData.paramVal[48] 
	geneG_act3_C = paramData.paramVal[52] 
	geneG_inh1_C = paramData.paramVal[53] 
	geneH_act1_C = paramData.paramVal[57] 
	geneI_act2_C = paramData.paramVal[65] 
	geneI_act3_C = paramData.paramVal[66] 
	geneI_inh1_C = paramData.paramVal[67] 
	geneJ_act1_C = paramData.paramVal[71] 
	geneJ_act2_C = paramData.paramVal[72] 
	geneJ_act3_C = paramData.paramVal[73] 
	geneK_act3_C = paramData.paramVal[80] 
	geneK_inh1_C = paramData.paramVal[81] 
	geneL_act1_C = paramData.paramVal[85] 
	geneL_act3_C = paramData.paramVal[87] 
	k_223_C = paramData.paramVal[93] 
	k_224_C = paramData.paramVal[94] 
	k_334_C = paramData.paramVal[98] 
	k_344_C = paramData.paramVal[99] 
	k_on_u_C = paramData.paramVal[100] 
	kdiss_SS_C = paramData.paramVal[101] 
	khomo2_C = paramData.paramVal[102] 
	khomo3_C = paramData.paramVal[103] 
	khomo4_C = paramData.paramVal[104] 

	if observableId == "Bmp4_gene_expression_OE_nExpID100" 
		return sd_Bmp4_nExpID100 
	end

	if observableId == "Cxcl15_gene_expression_OE_nExpID100" 
		return sd_Cxcl15_nExpID100 
	end

	if observableId == "Dnmt3a_gene_expression_OE_nExpID100" 
		return sd_Dnmt3a_nExpID100 
	end

	if observableId == "Dusp5_gene_expression_OE_nExpID100" 
		return sd_Dusp5_nExpID100 
	end

	if observableId == "Jun_gene_expression_OE_nExpID100" 
		return sd_Jun_nExpID100 
	end

	if observableId == "Klf10_gene_expression_OE_nExpID100" 
		return sd_Klf10_nExpID100 
	end

	if observableId == "Pdk4_gene_expression_OE_nExpID100" 
		return sd_Pdk4_nExpID100 
	end

	if observableId == "S23IP_S2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5" 
		return 5 
	end

	if observableId == "S23IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 5 
	end

	if observableId == "S23IP_S3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5" 
		return 5 
	end

	if observableId == "S23IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 5 
	end

	if observableId == "S23IP_pS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5" 
		return 5 
	end

	if observableId == "S23IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 5 
	end

	if observableId == "S23IP_pS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5" 
		return 5 
	end

	if observableId == "S23IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 5 
	end

	if observableId == "S23IP_ppS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5" 
		return 5 
	end

	if observableId == "S23IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 5 
	end

	if observableId == "S23IP_ppS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5" 
		return 5 
	end

	if observableId == "S23IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 5 
	end

	if observableId == "S23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1" 
		return 15 
	end

	if observableId == "S23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10" 
		return 1 
	end

	if observableId == "S23IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 5 
	end

	if observableId == "S23IP_tS3_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10" 
		return 1 
	end

	if observableId == "S23IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 5 
	end

	if observableId == "S23IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 5 
	end

	if observableId == "S2IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 5 
	end

	if observableId == "S2IP_S2_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return 5 
	end

	if observableId == "S2IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 5 
	end

	if observableId == "S2IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 5 
	end

	if observableId == "S2IP_pS2_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return 5 
	end

	if observableId == "S2IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 5 
	end

	if observableId == "S2IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 5 
	end

	if observableId == "S2IP_ppS2_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return 5 
	end

	if observableId == "S2IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 5 
	end

	if observableId == "S2IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 5 
	end

	if observableId == "S2IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return 5 
	end

	if observableId == "S2IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 5 
	end

	if observableId == "S2IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return 5 
	end

	if observableId == "S2IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 5 
	end

	if observableId == "S2IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return 5 
	end

	if observableId == "S3IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 5 
	end

	if observableId == "S3IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 5 
	end

	if observableId == "S3IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 5 
	end

	if observableId == "S3IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 5 
	end

	if observableId == "S3IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 5 
	end

	if observableId == "S3IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 5 
	end

	if observableId == "S3IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 5 
	end

	if observableId == "S3IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return 5 
	end

	if observableId == "S3IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 5 
	end

	if observableId == "S3IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return 5 
	end

	if observableId == "S3IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 5 
	end

	if observableId == "S3IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return 5 
	end

	if observableId == "S4IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 5 
	end

	if observableId == "S4IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return 5 
	end

	if observableId == "S4IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 5 
	end

	if observableId == "S4IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return 5 
	end

	if observableId == "S4IP_tS4_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1" 
		return 15 
	end

	if observableId == "S4IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return 5 
	end

	if observableId == "S4IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return 5 
	end

	if observableId == "Ski_gene_expression_OE_nExpID100" 
		return sd_Ski_nExpID100 
	end

	if observableId == "Skil_gene_expression_OE_nExpID100" 
		return sd_Skil_nExpID100 
	end

	if observableId == "Smad7_gene_expression_OE_nExpID100" 
		return sd_Smad7_nExpID100 
	end

	if observableId == "Sox4_gene_expression_OE_nExpID100" 
		return sd_Sox4_nExpID100 
	end

	if observableId == "TGFbR_PL_PN_TGFbR_hepa16_nExpID11" 
		return 0.092 
	end

	if observableId == "Tgfa_gene_expression_OE_nExpID100" 
		return sd_Tgfa_nExpID100 
	end

end