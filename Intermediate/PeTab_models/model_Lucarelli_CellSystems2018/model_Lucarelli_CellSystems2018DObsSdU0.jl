function evalDYmodDu(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	ppS3_ppS3_ppS3, ppS3_S4_S4, geneH, geneI, geneJ, TGFb_pRec, geneE, S4, ppS2_ppS2_ppS2, geneK, pS2, pS3, geneC, ppS2_ppS2_ppS3, geneF, S4_S4_S4, TGFb, S3, S2, S2_S4_S4, ppS3_ppS3_S4, geneB, ppS3, geneA, geneD, ppS2_S4_S4, geneL, ppS2_ppS3_S4, geneG, Rec, ppS2_ppS2_S4, ppS2_ppS3_ppS3, ppS2= u 
	geneC_inh3, geneJ_act1, geneG_inh1, init_Rec, geneC_act3, geneA_inh2, geneK_turn, geneA_act3, geneI_act2, geneA_act2, geneD_inh2, k_223, k_233, geneB_act1, S_dephosphos, geneG_act1, geneH_inh2, geneI_act1, geneB_inh3, geneI_inh2, geneC_turn, geneJ_inh1, geneK_act3, geneJ_inh3, geneG_inh3, geneC_act2, init_TGFb, geneL_inh1, geneE_inh3, geneF_act3, cell, geneF_inh1, geneD_act3, geneE_act2, geneD_inh3, k_234, geneH_act2, geneA_turn, geneL_act1, geneH_act1, geneL_inh2, geneB_turn, init_S4, geneL_inh3, khomo2, geneK_act2, geneA_inh3, geneJ_act2, geneI_inh3, geneD_act1, geneJ_turn, geneG_act3, geneL_act2, Rec_act, geneI_act3, k_224, pRec_degind, geneE_inh2, geneF_inh2, geneC_act1, geneD_inh1, k_344, init_S3, geneI_inh1, geneL_turn, geneF_act2, k_on_u, geneE_inh1, geneH_act3, geneH_inh1, geneK_inh3, geneE_act1, geneA_inh1, geneB_inh1, geneH_inh3, geneD_turn, S_dephos, geneL_act3, S_phos, geneG_turn, geneA_act1, geneI_turn, khomo3, geneC_inh1, geneG_act2, init_S2, geneK_inh2, k_334, geneB_inh2, geneH_turn, geneJ_inh2, khomo4, geneB_act2, geneF_turn, geneJ_act3, geneK_inh1, geneB_act3, geneK_act1, geneE_turn, geneE_act3, geneG_inh2, kdiss_SS, k_244, geneD_act2, geneF_inh3, geneF_act1, geneC_inh2 = p 
	sd_Bmp4_nExpID100, sd_Cxcl15_nExpID100, sd_Dnmt3a_nExpID100, sd_Dusp5_nExpID100, sd_Jun_nExpID100, sd_Klf10_nExpID100, sd_Pdk4_nExpID100, sd_Ski_nExpID100, sd_Skil_nExpID100, sd_Smad7_nExpID100, sd_Sox4_nExpID100, sd_Tgfa_nExpID100 = nonDynParam 
	if observableId == :Bmp4_gene_expression_OE_nExpID100 
		out[3] = 1
		return nothing
	end

	if observableId == :Cxcl15_gene_expression_OE_nExpID100 
		out[4] = 1
		return nothing
	end

	if observableId == :Dnmt3a_gene_expression_OE_nExpID100 
		out[13] = 1
		return nothing
	end

	if observableId == :Dusp5_gene_expression_OE_nExpID100 
		out[5] = 1
		return nothing
	end

	if observableId == :Jun_gene_expression_OE_nExpID100 
		out[7] = 1
		return nothing
	end

	if observableId == :Klf10_gene_expression_OE_nExpID100 
		out[29] = 1
		return nothing
	end

	if observableId == :Pdk4_gene_expression_OE_nExpID100 
		out[27] = 1
		return nothing
	end

	if observableId == :S23IP_S2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 
		out[9] = (-300S2 - 300S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[11] = (-100S2 - 100S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[14] = (-200S2 - 200S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[19] = (100.0pS2 + 100.0ppS2 + 100.0ppS2_S4_S4 + 100.0ppS2_ppS3_S4 + 100.0ppS2_ppS3_ppS3 + 200.0ppS2_ppS2_S4 + 200.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[20] = (100.0pS2 + 100.0ppS2 + 100.0ppS2_S4_S4 + 100.0ppS2_ppS3_S4 + 100.0ppS2_ppS3_ppS3 + 200.0ppS2_ppS2_S4 + 200.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[26] = (-100S2 - 100S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[28] = (-100S2 - 100S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[31] = (-200S2 - 200S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[32] = (-100S2 - 100S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[33] = (-100S2 - 100S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		return nothing
	end

	if observableId == :S23IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		out[9] = (-300S2 - 300S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[11] = (-100S2 - 100S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[14] = (-200S2 - 200S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[19] = (100.0pS2 + 100.0ppS2 + 100.0ppS2_S4_S4 + 100.0ppS2_ppS3_S4 + 100.0ppS2_ppS3_ppS3 + 200.0ppS2_ppS2_S4 + 200.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[20] = (100.0pS2 + 100.0ppS2 + 100.0ppS2_S4_S4 + 100.0ppS2_ppS3_S4 + 100.0ppS2_ppS3_ppS3 + 200.0ppS2_ppS2_S4 + 200.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[26] = (-100S2 - 100S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[28] = (-100S2 - 100S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[31] = (-200S2 - 200S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[32] = (-100S2 - 100S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[33] = (-100S2 - 100S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		return nothing
	end

	if observableId == :S23IP_S3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 
		out[1] = (-300S3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[2] = (-100S3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[12] = (-100S3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[14] = (-100S3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[18] = (100.0pS3 + 100.0ppS2_ppS2_ppS3 + 100.0ppS2_ppS3_S4 + 200.0ppS2_ppS3_ppS3 + 100.0ppS3 + 100.0ppS3_S4_S4 + 200.0ppS3_ppS3_S4 + 300.0ppS3_ppS3_ppS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[21] = (-200S3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[23] = (-100S3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[28] = (-100S3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[32] = (-200S3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		return nothing
	end

	if observableId == :S23IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		out[1] = (-300S3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[2] = (-100S3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[12] = (-100S3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[14] = (-100S3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[18] = (100.0pS3 + 100.0ppS2_ppS2_ppS3 + 100.0ppS2_ppS3_S4 + 200.0ppS2_ppS3_ppS3 + 100.0ppS3 + 100.0ppS3_S4_S4 + 200.0ppS3_ppS3_S4 + 300.0ppS3_ppS3_ppS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[21] = (-200S3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[23] = (-100S3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[28] = (-100S3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[32] = (-200S3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		return nothing
	end

	if observableId == :S23IP_pS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 
		out[9] = (-300pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[11] = (100.0S2 + 100.0S2_S4_S4 + 100.0ppS2 + 100.0ppS2_S4_S4 + 100.0ppS2_ppS3_S4 + 100.0ppS2_ppS3_ppS3 + 200.0ppS2_ppS2_S4 + 200.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[14] = (-200pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[19] = (-100pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[20] = (-100pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[26] = (-100pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[28] = (-100pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[31] = (-200pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[32] = (-100pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[33] = (-100pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		return nothing
	end

	if observableId == :S23IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		out[9] = (-300pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[11] = (100.0S2 + 100.0S2_S4_S4 + 100.0ppS2 + 100.0ppS2_S4_S4 + 100.0ppS2_ppS3_S4 + 100.0ppS2_ppS3_ppS3 + 200.0ppS2_ppS2_S4 + 200.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[14] = (-200pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[19] = (-100pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[20] = (-100pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[26] = (-100pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[28] = (-100pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[31] = (-200pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[32] = (-100pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[33] = (-100pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		return nothing
	end

	if observableId == :S23IP_pS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 
		out[1] = (-300pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[2] = (-100pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[12] = (100.0S3 + 100.0ppS2_ppS2_ppS3 + 100.0ppS2_ppS3_S4 + 200.0ppS2_ppS3_ppS3 + 100.0ppS3 + 100.0ppS3_S4_S4 + 200.0ppS3_ppS3_S4 + 300.0ppS3_ppS3_ppS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[14] = (-100pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[18] = (-100pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[21] = (-200pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[23] = (-100pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[28] = (-100pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[32] = (-200pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		return nothing
	end

	if observableId == :S23IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		out[1] = (-300pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[2] = (-100pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[12] = (100.0S3 + 100.0ppS2_ppS2_ppS3 + 100.0ppS2_ppS3_S4 + 200.0ppS2_ppS3_ppS3 + 100.0ppS3 + 100.0ppS3_S4_S4 + 200.0ppS3_ppS3_S4 + 300.0ppS3_ppS3_ppS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[14] = (-100pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[18] = (-100pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[21] = (-200pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[23] = (-100pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[28] = (-100pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[32] = (-200pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		return nothing
	end

	if observableId == :S23IP_ppS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 
		out[9] = (300.0S2 + 300.0S2_S4_S4 + 300.0pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[11] = (-100ppS2 - 100ppS2_S4_S4 - 100ppS2_ppS3_S4 - 100ppS2_ppS3_ppS3 - 200ppS2_ppS2_S4 - 200ppS2_ppS2_ppS3 - 300ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[14] = (200.0S2 + 200.0S2_S4_S4 + 200.0pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[19] = (-100ppS2 - 100ppS2_S4_S4 - 100ppS2_ppS3_S4 - 100ppS2_ppS3_ppS3 - 200ppS2_ppS2_S4 - 200ppS2_ppS2_ppS3 - 300ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[20] = (-100ppS2 - 100ppS2_S4_S4 - 100ppS2_ppS3_S4 - 100ppS2_ppS3_ppS3 - 200ppS2_ppS2_S4 - 200ppS2_ppS2_ppS3 - 300ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[26] = (100.0S2 + 100.0S2_S4_S4 + 100.0pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[28] = (100.0S2 + 100.0S2_S4_S4 + 100.0pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[31] = (200.0S2 + 200.0S2_S4_S4 + 200.0pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[32] = (100.0S2 + 100.0S2_S4_S4 + 100.0pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[33] = (100.0S2 + 100.0S2_S4_S4 + 100.0pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		return nothing
	end

	if observableId == :S23IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		out[9] = (300.0S2 + 300.0S2_S4_S4 + 300.0pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[11] = (-100ppS2 - 100ppS2_S4_S4 - 100ppS2_ppS3_S4 - 100ppS2_ppS3_ppS3 - 200ppS2_ppS2_S4 - 200ppS2_ppS2_ppS3 - 300ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[14] = (200.0S2 + 200.0S2_S4_S4 + 200.0pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[19] = (-100ppS2 - 100ppS2_S4_S4 - 100ppS2_ppS3_S4 - 100ppS2_ppS3_ppS3 - 200ppS2_ppS2_S4 - 200ppS2_ppS2_ppS3 - 300ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[20] = (-100ppS2 - 100ppS2_S4_S4 - 100ppS2_ppS3_S4 - 100ppS2_ppS3_ppS3 - 200ppS2_ppS2_S4 - 200ppS2_ppS2_ppS3 - 300ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[26] = (100.0S2 + 100.0S2_S4_S4 + 100.0pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[28] = (100.0S2 + 100.0S2_S4_S4 + 100.0pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[31] = (200.0S2 + 200.0S2_S4_S4 + 200.0pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[32] = (100.0S2 + 100.0S2_S4_S4 + 100.0pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[33] = (100.0S2 + 100.0S2_S4_S4 + 100.0pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		return nothing
	end

	if observableId == :S23IP_ppS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 
		out[1] = (300.0S3 + 300.0pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[2] = (100.0S3 + 100.0pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[12] = (-100ppS2_ppS2_ppS3 - 100ppS2_ppS3_S4 - 200ppS2_ppS3_ppS3 - 100ppS3 - 100ppS3_S4_S4 - 200ppS3_ppS3_S4 - 300ppS3_ppS3_ppS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[14] = (100.0S3 + 100.0pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[18] = (-100ppS2_ppS2_ppS3 - 100ppS2_ppS3_S4 - 200ppS2_ppS3_ppS3 - 100ppS3 - 100ppS3_S4_S4 - 200ppS3_ppS3_S4 - 300ppS3_ppS3_ppS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[21] = (200.0S3 + 200.0pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[23] = (100.0S3 + 100.0pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[28] = (100.0S3 + 100.0pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[32] = (200.0S3 + 200.0pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		return nothing
	end

	if observableId == :S23IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		out[1] = (300.0S3 + 300.0pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[2] = (100.0S3 + 100.0pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[12] = (-100ppS2_ppS2_ppS3 - 100ppS2_ppS3_S4 - 200ppS2_ppS3_ppS3 - 100ppS3 - 100ppS3_S4_S4 - 200ppS3_ppS3_S4 - 300ppS3_ppS3_ppS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[14] = (100.0S3 + 100.0pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[18] = (-100ppS2_ppS2_ppS3 - 100ppS2_ppS3_S4 - 200ppS2_ppS3_ppS3 - 100ppS3 - 100ppS3_S4_S4 - 200ppS3_ppS3_S4 - 300ppS3_ppS3_ppS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[21] = (200.0S3 + 200.0pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[23] = (100.0S3 + 100.0pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[28] = (100.0S3 + 100.0pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[32] = (200.0S3 + 200.0pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		return nothing
	end

	if observableId == :S23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1 
		out[9] = 3
		out[11] = 1
		out[14] = 2
		out[19] = 1
		out[20] = 1
		out[26] = 1
		out[28] = 1
		out[31] = 2
		out[32] = 1
		out[33] = 1
		return nothing
	end

	if observableId == :S23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10 
		out[1] = (-300S2 - 300S2_S4_S4 - 300pS2 - 300ppS2 - 300ppS2_S4_S4 - 300ppS2_ppS3_S4 - 300ppS2_ppS3_ppS3 - 600ppS2_ppS2_S4 - 600ppS2_ppS2_ppS3 - 900ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		out[2] = (-100S2 - 100S2_S4_S4 - 100pS2 - 100ppS2 - 100ppS2_S4_S4 - 100ppS2_ppS3_S4 - 100ppS2_ppS3_ppS3 - 200ppS2_ppS2_S4 - 200ppS2_ppS2_ppS3 - 300ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		out[9] = (300.0S3 + 300.0pS3 + 300.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS3_S4 + 600.0ppS2_ppS3_ppS3 + 300.0ppS3 + 300.0ppS3_S4_S4 + 600.0ppS3_ppS3_S4 + 900.0ppS3_ppS3_ppS3) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		out[11] = (100.0S3 + 100.0pS3 + 100.0ppS2_ppS2_ppS3 + 100.0ppS2_ppS3_S4 + 200.0ppS2_ppS3_ppS3 + 100.0ppS3 + 100.0ppS3_S4_S4 + 200.0ppS3_ppS3_S4 + 300.0ppS3_ppS3_ppS3) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		out[12] = (-100S2 - 100S2_S4_S4 - 100pS2 - 100ppS2 - 100ppS2_S4_S4 - 100ppS2_ppS3_S4 - 100ppS2_ppS3_ppS3 - 200ppS2_ppS2_S4 - 200ppS2_ppS2_ppS3 - 300ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		out[14] = (200.0S3 + 200.0pS3 + 100.0ppS2_ppS3_S4 + 200.0ppS3 + 200.0ppS3_S4_S4 + 300.0ppS2_ppS3_ppS3 + 400.0ppS3_ppS3_S4 + 600.0ppS3_ppS3_ppS3 - 100.0S2 - 100.0S2_S4_S4 - 100.0pS2 - 100.0ppS2 - 100.0ppS2_S4_S4 - 200.0ppS2_ppS2_S4 - 300.0ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		out[18] = (-100S2 - 100S2_S4_S4 - 100pS2 - 100ppS2 - 100ppS2_S4_S4 - 100ppS2_ppS3_S4 - 100ppS2_ppS3_ppS3 - 200ppS2_ppS2_S4 - 200ppS2_ppS2_ppS3 - 300ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		out[19] = (100.0S3 + 100.0pS3 + 100.0ppS2_ppS2_ppS3 + 100.0ppS2_ppS3_S4 + 200.0ppS2_ppS3_ppS3 + 100.0ppS3 + 100.0ppS3_S4_S4 + 200.0ppS3_ppS3_S4 + 300.0ppS3_ppS3_ppS3) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		out[20] = (100.0S3 + 100.0pS3 + 100.0ppS2_ppS2_ppS3 + 100.0ppS2_ppS3_S4 + 200.0ppS2_ppS3_ppS3 + 100.0ppS3 + 100.0ppS3_S4_S4 + 200.0ppS3_ppS3_S4 + 300.0ppS3_ppS3_ppS3) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		out[21] = (-200S2 - 200S2_S4_S4 - 200pS2 - 200ppS2 - 200ppS2_S4_S4 - 200ppS2_ppS3_S4 - 200ppS2_ppS3_ppS3 - 400ppS2_ppS2_S4 - 400ppS2_ppS2_ppS3 - 600ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		out[23] = (-100S2 - 100S2_S4_S4 - 100pS2 - 100ppS2 - 100ppS2_S4_S4 - 100ppS2_ppS3_S4 - 100ppS2_ppS3_ppS3 - 200ppS2_ppS2_S4 - 200ppS2_ppS2_ppS3 - 300ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		out[26] = (100.0S3 + 100.0pS3 + 100.0ppS2_ppS2_ppS3 + 100.0ppS2_ppS3_S4 + 200.0ppS2_ppS3_ppS3 + 100.0ppS3 + 100.0ppS3_S4_S4 + 200.0ppS3_ppS3_S4 + 300.0ppS3_ppS3_ppS3) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		out[28] = (100.0S3 + 100.0pS3 + 100.0ppS2_ppS3_ppS3 + 100.0ppS3 + 100.0ppS3_S4_S4 + 200.0ppS3_ppS3_S4 + 300.0ppS3_ppS3_ppS3 - 100.0S2 - 100.0S2_S4_S4 - 100.0pS2 - 100.0ppS2 - 100.0ppS2_S4_S4 - 200.0ppS2_ppS2_S4 - 100.0ppS2_ppS2_ppS3 - 300.0ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		out[31] = (200.0S3 + 200.0pS3 + 200.0ppS2_ppS2_ppS3 + 200.0ppS2_ppS3_S4 + 400.0ppS2_ppS3_ppS3 + 200.0ppS3 + 200.0ppS3_S4_S4 + 400.0ppS3_ppS3_S4 + 600.0ppS3_ppS3_ppS3) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		out[32] = (100.0S3 + 100.0pS3 + 100.0ppS3 + 100.0ppS3_S4_S4 + 200.0ppS3_ppS3_S4 + 300.0ppS3_ppS3_ppS3 - 200.0S2 - 200.0S2_S4_S4 - 200.0pS2 - 200.0ppS2 - 200.0ppS2_S4_S4 - 400.0ppS2_ppS2_S4 - 300.0ppS2_ppS2_ppS3 - 100.0ppS2_ppS3_S4 - 600.0ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		out[33] = (100.0S3 + 100.0pS3 + 100.0ppS2_ppS2_ppS3 + 100.0ppS2_ppS3_S4 + 200.0ppS2_ppS3_ppS3 + 100.0ppS3 + 100.0ppS3_S4_S4 + 200.0ppS3_ppS3_S4 + 300.0ppS3_ppS3_ppS3) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		return nothing
	end

	if observableId == :S23IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		out[1] = (-300S2 - 300S2_S4_S4 - 300pS2 - 300ppS2 - 300ppS2_S4_S4 - 300ppS2_ppS3_S4 - 300ppS2_ppS3_ppS3 - 600ppS2_ppS2_S4 - 600ppS2_ppS2_ppS3 - 900ppS2_ppS2_ppS2) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[2] = (-300S2 - 300S2_S4_S4 - 300pS2 - 300ppS2 - 300ppS2_S4_S4 - 300ppS2_ppS3_S4 - 300ppS2_ppS3_ppS3 - 600ppS2_ppS2_S4 - 600ppS2_ppS2_ppS3 - 900ppS2_ppS2_ppS2) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[9] = (600.0S2_S4_S4 + 300.0S3 + 300.0pS3 + 600.0ppS2_S4_S4 + 300.0ppS2_ppS2_S4 + 300.0ppS2_ppS2_ppS3 + 300.0ppS3 + 600.0ppS2_ppS3_ppS3 + 900.0ppS2_ppS3_S4 + 900.0ppS3_S4_S4 + 900.0ppS3_ppS3_S4 + 900.0ppS3_ppS3_ppS3) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[11] = (200.0S2_S4_S4 + 100.0S3 + 100.0pS3 + 200.0ppS2_S4_S4 + 100.0ppS2_ppS2_S4 + 100.0ppS2_ppS2_ppS3 + 100.0ppS3 + 200.0ppS2_ppS3_ppS3 + 300.0ppS2_ppS3_S4 + 300.0ppS3_S4_S4 + 300.0ppS3_ppS3_S4 + 300.0ppS3_ppS3_ppS3) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[12] = (-100S2 - 100S2_S4_S4 - 100pS2 - 100ppS2 - 100ppS2_S4_S4 - 100ppS2_ppS3_S4 - 100ppS2_ppS3_ppS3 - 200ppS2_ppS2_S4 - 200ppS2_ppS2_ppS3 - 300ppS2_ppS2_ppS2) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[14] = (300.0S2_S4_S4 + 200.0S3 + 200.0pS3 + 300.0ppS2_S4_S4 + 200.0ppS3 + 300.0ppS2_ppS3_ppS3 + 500.0ppS2_ppS3_S4 + 600.0ppS3_S4_S4 + 600.0ppS3_ppS3_S4 + 600.0ppS3_ppS3_ppS3 - 100.0S2 - 100.0pS2 - 100.0ppS2 - 300.0ppS2_ppS2_ppS2) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[18] = (-100S2 - 100S2_S4_S4 - 100pS2 - 100ppS2 - 100ppS2_S4_S4 - 100ppS2_ppS3_S4 - 100ppS2_ppS3_ppS3 - 200ppS2_ppS2_S4 - 200ppS2_ppS2_ppS3 - 300ppS2_ppS2_ppS2) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[19] = (200.0S2_S4_S4 + 100.0S3 + 100.0pS3 + 200.0ppS2_S4_S4 + 100.0ppS2_ppS2_S4 + 100.0ppS2_ppS2_ppS3 + 100.0ppS3 + 200.0ppS2_ppS3_ppS3 + 300.0ppS2_ppS3_S4 + 300.0ppS3_S4_S4 + 300.0ppS3_ppS3_S4 + 300.0ppS3_ppS3_ppS3) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[20] = (100.0S3 + 100.0pS3 + 100.0ppS2_ppS3_S4 + 100.0ppS3 + 300.0ppS3_S4_S4 + 300.0ppS3_ppS3_S4 + 300.0ppS3_ppS3_ppS3 - 200.0S2 - 200.0pS2 - 200.0ppS2 - 300.0ppS2_ppS2_S4 - 300.0ppS2_ppS2_ppS3 - 600.0ppS2_ppS2_ppS2) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[21] = (-300S2 - 300S2_S4_S4 - 300pS2 - 300ppS2 - 300ppS2_S4_S4 - 300ppS2_ppS3_S4 - 300ppS2_ppS3_ppS3 - 600ppS2_ppS2_S4 - 600ppS2_ppS2_ppS3 - 900ppS2_ppS2_ppS2) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[23] = (-100S2 - 100S2_S4_S4 - 100pS2 - 100ppS2 - 100ppS2_S4_S4 - 100ppS2_ppS3_S4 - 100ppS2_ppS3_ppS3 - 200ppS2_ppS2_S4 - 200ppS2_ppS2_ppS3 - 300ppS2_ppS2_ppS2) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[26] = (100.0S3 + 100.0pS3 + 100.0ppS2_ppS3_S4 + 100.0ppS3 + 300.0ppS3_S4_S4 + 300.0ppS3_ppS3_S4 + 300.0ppS3_ppS3_ppS3 - 200.0S2 - 200.0pS2 - 200.0ppS2 - 300.0ppS2_ppS2_S4 - 300.0ppS2_ppS2_ppS3 - 600.0ppS2_ppS2_ppS2) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[28] = (100.0S3 + 100.0pS3 + 100.0ppS3 + 300.0ppS3_S4_S4 + 300.0ppS3_ppS3_S4 + 300.0ppS3_ppS3_ppS3 - 100.0S2_S4_S4 - 100.0ppS2_S4_S4 - 300.0S2 - 300.0pS2 - 300.0ppS2 - 500.0ppS2_ppS2_S4 - 500.0ppS2_ppS2_ppS3 - 100.0ppS2_ppS3_ppS3 - 900.0ppS2_ppS2_ppS2) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[31] = (300.0S2_S4_S4 + 200.0S3 + 200.0pS3 + 300.0ppS2_S4_S4 + 200.0ppS3 + 300.0ppS2_ppS3_ppS3 + 500.0ppS2_ppS3_S4 + 600.0ppS3_S4_S4 + 600.0ppS3_ppS3_S4 + 600.0ppS3_ppS3_ppS3 - 100.0S2 - 100.0pS2 - 100.0ppS2 - 300.0ppS2_ppS2_ppS2) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[32] = (100.0S3 + 100.0pS3 + 100.0ppS2_ppS3_S4 + 100.0ppS3 + 300.0ppS3_S4_S4 + 300.0ppS3_ppS3_S4 + 300.0ppS3_ppS3_ppS3 - 200.0S2 - 200.0pS2 - 200.0ppS2 - 300.0ppS2_ppS2_S4 - 300.0ppS2_ppS2_ppS3 - 600.0ppS2_ppS2_ppS2) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[33] = (200.0S2_S4_S4 + 100.0S3 + 100.0pS3 + 200.0ppS2_S4_S4 + 100.0ppS2_ppS2_S4 + 100.0ppS2_ppS2_ppS3 + 100.0ppS3 + 200.0ppS2_ppS3_ppS3 + 300.0ppS2_ppS3_S4 + 300.0ppS3_S4_S4 + 300.0ppS3_ppS3_S4 + 300.0ppS3_ppS3_ppS3) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		return nothing
	end

	if observableId == :S23IP_tS3_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10 
		out[1] = (300.0S2 + 300.0S2_S4_S4 + 300.0pS2 + 300.0ppS2 + 300.0ppS2_S4_S4 + 300.0ppS2_ppS3_S4 + 300.0ppS2_ppS3_ppS3 + 600.0ppS2_ppS2_S4 + 600.0ppS2_ppS2_ppS3 + 900.0ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		out[2] = (100.0S2 + 100.0S2_S4_S4 + 100.0pS2 + 100.0ppS2 + 100.0ppS2_S4_S4 + 100.0ppS2_ppS3_S4 + 100.0ppS2_ppS3_ppS3 + 200.0ppS2_ppS2_S4 + 200.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		out[9] = (-300S3 - 300pS3 - 300ppS2_ppS2_ppS3 - 300ppS2_ppS3_S4 - 600ppS2_ppS3_ppS3 - 300ppS3 - 300ppS3_S4_S4 - 600ppS3_ppS3_S4 - 900ppS3_ppS3_ppS3) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		out[11] = (-100S3 - 100pS3 - 100ppS2_ppS2_ppS3 - 100ppS2_ppS3_S4 - 200ppS2_ppS3_ppS3 - 100ppS3 - 100ppS3_S4_S4 - 200ppS3_ppS3_S4 - 300ppS3_ppS3_ppS3) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		out[12] = (100.0S2 + 100.0S2_S4_S4 + 100.0pS2 + 100.0ppS2 + 100.0ppS2_S4_S4 + 100.0ppS2_ppS3_S4 + 100.0ppS2_ppS3_ppS3 + 200.0ppS2_ppS2_S4 + 200.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		out[14] = (100.0S2 + 100.0S2_S4_S4 + 100.0pS2 + 100.0ppS2 + 100.0ppS2_S4_S4 + 200.0ppS2_ppS2_S4 + 300.0ppS2_ppS2_ppS2 - 200.0S3 - 200.0pS3 - 100.0ppS2_ppS3_S4 - 200.0ppS3 - 200.0ppS3_S4_S4 - 300.0ppS2_ppS3_ppS3 - 400.0ppS3_ppS3_S4 - 600.0ppS3_ppS3_ppS3) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		out[18] = (100.0S2 + 100.0S2_S4_S4 + 100.0pS2 + 100.0ppS2 + 100.0ppS2_S4_S4 + 100.0ppS2_ppS3_S4 + 100.0ppS2_ppS3_ppS3 + 200.0ppS2_ppS2_S4 + 200.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		out[19] = (-100S3 - 100pS3 - 100ppS2_ppS2_ppS3 - 100ppS2_ppS3_S4 - 200ppS2_ppS3_ppS3 - 100ppS3 - 100ppS3_S4_S4 - 200ppS3_ppS3_S4 - 300ppS3_ppS3_ppS3) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		out[20] = (-100S3 - 100pS3 - 100ppS2_ppS2_ppS3 - 100ppS2_ppS3_S4 - 200ppS2_ppS3_ppS3 - 100ppS3 - 100ppS3_S4_S4 - 200ppS3_ppS3_S4 - 300ppS3_ppS3_ppS3) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		out[21] = (200.0S2 + 200.0S2_S4_S4 + 200.0pS2 + 200.0ppS2 + 200.0ppS2_S4_S4 + 200.0ppS2_ppS3_S4 + 200.0ppS2_ppS3_ppS3 + 400.0ppS2_ppS2_S4 + 400.0ppS2_ppS2_ppS3 + 600.0ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		out[23] = (100.0S2 + 100.0S2_S4_S4 + 100.0pS2 + 100.0ppS2 + 100.0ppS2_S4_S4 + 100.0ppS2_ppS3_S4 + 100.0ppS2_ppS3_ppS3 + 200.0ppS2_ppS2_S4 + 200.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		out[26] = (-100S3 - 100pS3 - 100ppS2_ppS2_ppS3 - 100ppS2_ppS3_S4 - 200ppS2_ppS3_ppS3 - 100ppS3 - 100ppS3_S4_S4 - 200ppS3_ppS3_S4 - 300ppS3_ppS3_ppS3) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		out[28] = (100.0S2 + 100.0S2_S4_S4 + 100.0pS2 + 100.0ppS2 + 100.0ppS2_S4_S4 + 200.0ppS2_ppS2_S4 + 100.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS2_ppS2 - 100.0S3 - 100.0pS3 - 100.0ppS2_ppS3_ppS3 - 100.0ppS3 - 100.0ppS3_S4_S4 - 200.0ppS3_ppS3_S4 - 300.0ppS3_ppS3_ppS3) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		out[31] = (-200S3 - 200pS3 - 200ppS2_ppS2_ppS3 - 200ppS2_ppS3_S4 - 400ppS2_ppS3_ppS3 - 200ppS3 - 200ppS3_S4_S4 - 400ppS3_ppS3_S4 - 600ppS3_ppS3_ppS3) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		out[32] = (200.0S2 + 200.0S2_S4_S4 + 200.0pS2 + 200.0ppS2 + 200.0ppS2_S4_S4 + 400.0ppS2_ppS2_S4 + 300.0ppS2_ppS2_ppS3 + 100.0ppS2_ppS3_S4 + 600.0ppS2_ppS2_ppS2 - 100.0S3 - 100.0pS3 - 100.0ppS3 - 100.0ppS3_S4_S4 - 200.0ppS3_ppS3_S4 - 300.0ppS3_ppS3_ppS3) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		out[33] = (-100S3 - 100pS3 - 100ppS2_ppS2_ppS3 - 100ppS2_ppS3_S4 - 200ppS2_ppS3_ppS3 - 100ppS3 - 100ppS3_S4_S4 - 200ppS3_ppS3_S4 - 300ppS3_ppS3_ppS3) / ((S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS2_S4 + 2ppS2_ppS3_S4 + 2ppS3_ppS3_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_ppS3_ppS3)^2)
		return nothing
	end

	if observableId == :S23IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		out[1] = (300.0S2 + 300.0pS2 + 300.0ppS2 + 900.0S2_S4_S4 + 900.0ppS2_S4_S4 + 900.0ppS2_ppS2_S4 + 600.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS3_ppS3 + 600.0ppS3_S4_S4 + 300.0ppS3_ppS3_S4 + 900.0ppS2_ppS2_ppS2 + 900.0ppS2_ppS3_S4) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[2] = (100.0S2 + 100.0pS2 + 100.0ppS2 + 300.0S2_S4_S4 + 300.0ppS2_S4_S4 + 300.0ppS2_ppS2_S4 + 100.0ppS2_ppS3_S4 + 300.0ppS2_ppS2_ppS2 - 200.0S3 - 200.0pS3 - 200.0ppS3 - 300.0ppS2_ppS3_ppS3 - 300.0ppS3_ppS3_S4 - 600.0ppS3_ppS3_ppS3) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[9] = (-300S3 - 300pS3 - 300ppS2_ppS2_ppS3 - 300ppS2_ppS3_S4 - 600ppS2_ppS3_ppS3 - 300ppS3 - 300ppS3_S4_S4 - 600ppS3_ppS3_S4 - 900ppS3_ppS3_ppS3) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[11] = (-100S3 - 100pS3 - 100ppS2_ppS2_ppS3 - 100ppS2_ppS3_S4 - 200ppS2_ppS3_ppS3 - 100ppS3 - 100ppS3_S4_S4 - 200ppS3_ppS3_S4 - 300ppS3_ppS3_ppS3) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[12] = (100.0S2 + 100.0pS2 + 100.0ppS2 + 300.0S2_S4_S4 + 300.0ppS2_S4_S4 + 300.0ppS2_ppS2_S4 + 200.0ppS2_ppS2_ppS3 + 100.0ppS2_ppS3_ppS3 + 200.0ppS3_S4_S4 + 100.0ppS3_ppS3_S4 + 300.0ppS2_ppS2_ppS2 + 300.0ppS2_ppS3_S4) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[14] = (100.0S2 + 100.0pS2 + 100.0ppS2 + 300.0S2_S4_S4 + 300.0ppS2_S4_S4 + 300.0ppS2_ppS2_S4 + 100.0ppS2_ppS3_S4 + 300.0ppS2_ppS2_ppS2 - 200.0S3 - 200.0pS3 - 200.0ppS3 - 300.0ppS2_ppS3_ppS3 - 300.0ppS3_ppS3_S4 - 600.0ppS3_ppS3_ppS3) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[18] = (100.0S2 + 100.0pS2 + 100.0ppS2 + 300.0S2_S4_S4 + 300.0ppS2_S4_S4 + 300.0ppS2_ppS2_S4 + 200.0ppS2_ppS2_ppS3 + 100.0ppS2_ppS3_ppS3 + 200.0ppS3_S4_S4 + 100.0ppS3_ppS3_S4 + 300.0ppS2_ppS2_ppS2 + 300.0ppS2_ppS3_S4) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[19] = (-100S3 - 100pS3 - 100ppS2_ppS2_ppS3 - 100ppS2_ppS3_S4 - 200ppS2_ppS3_ppS3 - 100ppS3 - 100ppS3_S4_S4 - 200ppS3_ppS3_S4 - 300ppS3_ppS3_ppS3) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[20] = (-300S3 - 300pS3 - 300ppS2_ppS2_ppS3 - 300ppS2_ppS3_S4 - 600ppS2_ppS3_ppS3 - 300ppS3 - 300ppS3_S4_S4 - 600ppS3_ppS3_S4 - 900ppS3_ppS3_ppS3) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[21] = (200.0S2 + 200.0pS2 + 200.0ppS2 + 600.0S2_S4_S4 + 600.0ppS2_S4_S4 + 600.0ppS2_ppS2_S4 + 300.0ppS2_ppS2_ppS3 + 300.0ppS3_S4_S4 + 500.0ppS2_ppS3_S4 + 600.0ppS2_ppS2_ppS2 - 100.0S3 - 100.0pS3 - 100.0ppS3 - 300.0ppS3_ppS3_ppS3) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[23] = (100.0S2 + 100.0pS2 + 100.0ppS2 + 300.0S2_S4_S4 + 300.0ppS2_S4_S4 + 300.0ppS2_ppS2_S4 + 200.0ppS2_ppS2_ppS3 + 100.0ppS2_ppS3_ppS3 + 200.0ppS3_S4_S4 + 100.0ppS3_ppS3_S4 + 300.0ppS2_ppS2_ppS2 + 300.0ppS2_ppS3_S4) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[26] = (-300S3 - 300pS3 - 300ppS2_ppS2_ppS3 - 300ppS2_ppS3_S4 - 600ppS2_ppS3_ppS3 - 300ppS3 - 300ppS3_S4_S4 - 600ppS3_ppS3_S4 - 900ppS3_ppS3_ppS3) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[28] = (100.0S2 + 100.0pS2 + 100.0ppS2 + 300.0S2_S4_S4 + 300.0ppS2_S4_S4 + 300.0ppS2_ppS2_S4 + 300.0ppS2_ppS2_ppS2 - 300.0S3 - 300.0pS3 - 100.0ppS2_ppS2_ppS3 - 100.0ppS3_S4_S4 - 300.0ppS3 - 500.0ppS2_ppS3_ppS3 - 500.0ppS3_ppS3_S4 - 900.0ppS3_ppS3_ppS3) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[31] = (-300S3 - 300pS3 - 300ppS2_ppS2_ppS3 - 300ppS2_ppS3_S4 - 600ppS2_ppS3_ppS3 - 300ppS3 - 300ppS3_S4_S4 - 600ppS3_ppS3_S4 - 900ppS3_ppS3_ppS3) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[32] = (200.0S2 + 200.0pS2 + 200.0ppS2 + 600.0S2_S4_S4 + 600.0ppS2_S4_S4 + 600.0ppS2_ppS2_S4 + 300.0ppS2_ppS2_ppS3 + 300.0ppS3_S4_S4 + 500.0ppS2_ppS3_S4 + 600.0ppS2_ppS2_ppS2 - 100.0S3 - 100.0pS3 - 100.0ppS3 - 300.0ppS3_ppS3_ppS3) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[33] = (-100S3 - 100pS3 - 100ppS2_ppS2_ppS3 - 100ppS2_ppS3_S4 - 200ppS2_ppS3_ppS3 - 100ppS3 - 100ppS3_S4_S4 - 200ppS3_ppS3_S4 - 300ppS3_ppS3_ppS3) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		return nothing
	end

	if observableId == :S23IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		out[1] = (-600S2_S4_S4 - 600ppS2_S4_S4 - 300ppS2_ppS2_S4 - 600ppS2_ppS3_S4 - 600ppS3_S4_S4 - 300ppS3_ppS3_S4) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[2] = (200.0S2 + 200.0S3 + 200.0pS2 + 200.0pS3 + 200.0ppS2 + 300.0ppS2_ppS2_S4 + 200.0ppS2_ppS3_S4 + 200.0ppS3 + 300.0ppS3_ppS3_S4 + 600.0ppS2_ppS2_ppS2 + 600.0ppS2_ppS2_ppS3 + 600.0ppS2_ppS3_ppS3 + 600.0ppS3_ppS3_ppS3) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[9] = (-600S2_S4_S4 - 600ppS2_S4_S4 - 300ppS2_ppS2_S4 - 600ppS2_ppS3_S4 - 600ppS3_S4_S4 - 300ppS3_ppS3_S4) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[11] = (-200S2_S4_S4 - 200ppS2_S4_S4 - 100ppS2_ppS2_S4 - 200ppS2_ppS3_S4 - 200ppS3_S4_S4 - 100ppS3_ppS3_S4) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[12] = (-200S2_S4_S4 - 200ppS2_S4_S4 - 100ppS2_ppS2_S4 - 200ppS2_ppS3_S4 - 200ppS3_S4_S4 - 100ppS3_ppS3_S4) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[14] = (-600S2_S4_S4 - 600ppS2_S4_S4 - 300ppS2_ppS2_S4 - 600ppS2_ppS3_S4 - 600ppS3_S4_S4 - 300ppS3_ppS3_S4) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[18] = (-200S2_S4_S4 - 200ppS2_S4_S4 - 100ppS2_ppS2_S4 - 200ppS2_ppS3_S4 - 200ppS3_S4_S4 - 100ppS3_ppS3_S4) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[19] = (-200S2_S4_S4 - 200ppS2_S4_S4 - 100ppS2_ppS2_S4 - 200ppS2_ppS3_S4 - 200ppS3_S4_S4 - 100ppS3_ppS3_S4) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[20] = (200.0S2 + 200.0S3 + 200.0pS2 + 200.0pS3 + 200.0ppS2 + 300.0ppS2_ppS2_S4 + 200.0ppS2_ppS3_S4 + 200.0ppS3 + 300.0ppS3_ppS3_S4 + 600.0ppS2_ppS2_ppS2 + 600.0ppS2_ppS2_ppS3 + 600.0ppS2_ppS3_ppS3 + 600.0ppS3_ppS3_ppS3) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[21] = (100.0S2 + 100.0S3 + 100.0pS2 + 100.0pS3 + 100.0ppS2 + 100.0ppS3 + 300.0ppS2_ppS2_ppS2 + 300.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS3_ppS3 + 300.0ppS3_ppS3_ppS3 - 300.0S2_S4_S4 - 300.0ppS2_S4_S4 - 200.0ppS2_ppS3_S4 - 300.0ppS3_S4_S4) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[23] = (-200S2_S4_S4 - 200ppS2_S4_S4 - 100ppS2_ppS2_S4 - 200ppS2_ppS3_S4 - 200ppS3_S4_S4 - 100ppS3_ppS3_S4) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[26] = (200.0S2 + 200.0S3 + 200.0pS2 + 200.0pS3 + 200.0ppS2 + 300.0ppS2_ppS2_S4 + 200.0ppS2_ppS3_S4 + 200.0ppS3 + 300.0ppS3_ppS3_S4 + 600.0ppS2_ppS2_ppS2 + 600.0ppS2_ppS2_ppS3 + 600.0ppS2_ppS3_ppS3 + 600.0ppS3_ppS3_ppS3) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[28] = (200.0S2 + 200.0S3 + 200.0pS2 + 200.0pS3 + 200.0ppS2 + 200.0ppS2_ppS2_S4 + 200.0ppS3 + 200.0ppS3_ppS3_S4 + 600.0ppS2_ppS2_ppS2 + 600.0ppS2_ppS2_ppS3 + 600.0ppS2_ppS3_ppS3 + 600.0ppS3_ppS3_ppS3 - 200.0S2_S4_S4 - 200.0ppS2_S4_S4 - 200.0ppS3_S4_S4) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[31] = (100.0S2 + 100.0S3 + 100.0pS2 + 100.0pS3 + 100.0ppS2 + 100.0ppS3 + 300.0ppS2_ppS2_ppS2 + 300.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS3_ppS3 + 300.0ppS3_ppS3_ppS3 - 300.0S2_S4_S4 - 300.0ppS2_S4_S4 - 200.0ppS2_ppS3_S4 - 300.0ppS3_S4_S4) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[32] = (-600S2_S4_S4 - 600ppS2_S4_S4 - 300ppS2_ppS2_S4 - 600ppS2_ppS3_S4 - 600ppS3_S4_S4 - 300ppS3_ppS3_S4) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		out[33] = (-200S2_S4_S4 - 200ppS2_S4_S4 - 100ppS2_ppS2_S4 - 200ppS2_ppS3_S4 - 200ppS3_S4_S4 - 100ppS3_ppS3_S4) / ((S2 + S3 + pS2 + pS3 + ppS2 + ppS3 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3 + 4ppS2_ppS3_S4)^2)
		return nothing
	end

	if observableId == :S2IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		out[9] = (-300S2 - 300S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[11] = (-100S2 - 100S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[14] = (-200S2 - 200S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[19] = (100.0pS2 + 100.0ppS2 + 100.0ppS2_S4_S4 + 100.0ppS2_ppS3_S4 + 100.0ppS2_ppS3_ppS3 + 200.0ppS2_ppS2_S4 + 200.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[20] = (100.0pS2 + 100.0ppS2 + 100.0ppS2_S4_S4 + 100.0ppS2_ppS3_S4 + 100.0ppS2_ppS3_ppS3 + 200.0ppS2_ppS2_S4 + 200.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[26] = (-100S2 - 100S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[28] = (-100S2 - 100S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[31] = (-200S2 - 200S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[32] = (-100S2 - 100S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[33] = (-100S2 - 100S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		return nothing
	end

	if observableId == :S2IP_S2_MSPL2_complexes_all_hepa16_time_course_nExpID6 
		out[9] = (-300S2 - 300S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[11] = (-100S2 - 100S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[14] = (-200S2 - 200S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[19] = (100.0pS2 + 100.0ppS2 + 100.0ppS2_S4_S4 + 100.0ppS2_ppS3_S4 + 100.0ppS2_ppS3_ppS3 + 200.0ppS2_ppS2_S4 + 200.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[20] = (100.0pS2 + 100.0ppS2 + 100.0ppS2_S4_S4 + 100.0ppS2_ppS3_S4 + 100.0ppS2_ppS3_ppS3 + 200.0ppS2_ppS2_S4 + 200.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[26] = (-100S2 - 100S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[28] = (-100S2 - 100S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[31] = (-200S2 - 200S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[32] = (-100S2 - 100S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[33] = (-100S2 - 100S2_S4_S4) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		return nothing
	end

	if observableId == :S2IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S2IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		out[9] = (-300pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[11] = (100.0S2 + 100.0S2_S4_S4 + 100.0ppS2 + 100.0ppS2_S4_S4 + 100.0ppS2_ppS3_S4 + 100.0ppS2_ppS3_ppS3 + 200.0ppS2_ppS2_S4 + 200.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[14] = (-200pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[19] = (-100pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[20] = (-100pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[26] = (-100pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[28] = (-100pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[31] = (-200pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[32] = (-100pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[33] = (-100pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		return nothing
	end

	if observableId == :S2IP_pS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 
		out[9] = (-300pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[11] = (100.0S2 + 100.0S2_S4_S4 + 100.0ppS2 + 100.0ppS2_S4_S4 + 100.0ppS2_ppS3_S4 + 100.0ppS2_ppS3_ppS3 + 200.0ppS2_ppS2_S4 + 200.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[14] = (-200pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[19] = (-100pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[20] = (-100pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[26] = (-100pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[28] = (-100pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[31] = (-200pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[32] = (-100pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[33] = (-100pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		return nothing
	end

	if observableId == :S2IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S2IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		out[9] = (300.0S2 + 300.0S2_S4_S4 + 300.0pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[11] = (-100ppS2 - 100ppS2_S4_S4 - 100ppS2_ppS3_S4 - 100ppS2_ppS3_ppS3 - 200ppS2_ppS2_S4 - 200ppS2_ppS2_ppS3 - 300ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[14] = (200.0S2 + 200.0S2_S4_S4 + 200.0pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[19] = (-100ppS2 - 100ppS2_S4_S4 - 100ppS2_ppS3_S4 - 100ppS2_ppS3_ppS3 - 200ppS2_ppS2_S4 - 200ppS2_ppS2_ppS3 - 300ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[20] = (-100ppS2 - 100ppS2_S4_S4 - 100ppS2_ppS3_S4 - 100ppS2_ppS3_ppS3 - 200ppS2_ppS2_S4 - 200ppS2_ppS2_ppS3 - 300ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[26] = (100.0S2 + 100.0S2_S4_S4 + 100.0pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[28] = (100.0S2 + 100.0S2_S4_S4 + 100.0pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[31] = (200.0S2 + 200.0S2_S4_S4 + 200.0pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[32] = (100.0S2 + 100.0S2_S4_S4 + 100.0pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[33] = (100.0S2 + 100.0S2_S4_S4 + 100.0pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		return nothing
	end

	if observableId == :S2IP_ppS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 
		out[9] = (300.0S2 + 300.0S2_S4_S4 + 300.0pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[11] = (-100ppS2 - 100ppS2_S4_S4 - 100ppS2_ppS3_S4 - 100ppS2_ppS3_ppS3 - 200ppS2_ppS2_S4 - 200ppS2_ppS2_ppS3 - 300ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[14] = (200.0S2 + 200.0S2_S4_S4 + 200.0pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[19] = (-100ppS2 - 100ppS2_S4_S4 - 100ppS2_ppS3_S4 - 100ppS2_ppS3_ppS3 - 200ppS2_ppS2_S4 - 200ppS2_ppS2_ppS3 - 300ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[20] = (-100ppS2 - 100ppS2_S4_S4 - 100ppS2_ppS3_S4 - 100ppS2_ppS3_ppS3 - 200ppS2_ppS2_S4 - 200ppS2_ppS2_ppS3 - 300ppS2_ppS2_ppS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[26] = (100.0S2 + 100.0S2_S4_S4 + 100.0pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[28] = (100.0S2 + 100.0S2_S4_S4 + 100.0pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[31] = (200.0S2 + 200.0S2_S4_S4 + 200.0pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[32] = (100.0S2 + 100.0S2_S4_S4 + 100.0pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		out[33] = (100.0S2 + 100.0S2_S4_S4 + 100.0pS2) / ((S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3 + 2ppS2_ppS2_S4 + 2ppS2_ppS2_ppS3 + 3ppS2_ppS2_ppS2)^2)
		return nothing
	end

	if observableId == :S2IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S2IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		out[9] = (600.0S2_S4_S4 + 600.0ppS2_S4_S4 + 300.0ppS2_ppS2_S4 + 300.0ppS2_ppS2_ppS3 + 600.0ppS2_ppS3_S4 + 600.0ppS2_ppS3_ppS3) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[11] = (200.0S2_S4_S4 + 200.0ppS2_S4_S4 + 100.0ppS2_ppS2_S4 + 100.0ppS2_ppS2_ppS3 + 200.0ppS2_ppS3_S4 + 200.0ppS2_ppS3_ppS3) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[14] = (300.0S2_S4_S4 + 300.0ppS2_S4_S4 + 300.0ppS2_ppS3_S4 + 300.0ppS2_ppS3_ppS3 - 100.0S2 - 100.0pS2 - 100.0ppS2 - 300.0ppS2_ppS2_ppS2) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[19] = (200.0S2_S4_S4 + 200.0ppS2_S4_S4 + 100.0ppS2_ppS2_S4 + 100.0ppS2_ppS2_ppS3 + 200.0ppS2_ppS3_S4 + 200.0ppS2_ppS3_ppS3) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[20] = (-200.0S2 - 200.0pS2 - 200.0ppS2 - 300.0ppS2_ppS2_S4 - 300.0ppS2_ppS2_ppS3 - 600.0ppS2_ppS2_ppS2) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[26] = (-200.0S2 - 200.0pS2 - 200.0ppS2 - 300.0ppS2_ppS2_S4 - 300.0ppS2_ppS2_ppS3 - 600.0ppS2_ppS2_ppS2) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[28] = (-200.0S2 - 200.0pS2 - 200.0ppS2 - 300.0ppS2_ppS2_S4 - 300.0ppS2_ppS2_ppS3 - 600.0ppS2_ppS2_ppS2) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[31] = (300.0S2_S4_S4 + 300.0ppS2_S4_S4 + 300.0ppS2_ppS3_S4 + 300.0ppS2_ppS3_ppS3 - 100.0S2 - 100.0pS2 - 100.0ppS2 - 300.0ppS2_ppS2_ppS2) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[32] = (-200.0S2 - 200.0pS2 - 200.0ppS2 - 300.0ppS2_ppS2_S4 - 300.0ppS2_ppS2_ppS3 - 600.0ppS2_ppS2_ppS2) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[33] = (200.0S2_S4_S4 + 200.0ppS2_S4_S4 + 100.0ppS2_ppS2_S4 + 100.0ppS2_ppS2_ppS3 + 200.0ppS2_ppS3_S4 + 200.0ppS2_ppS3_ppS3) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		return nothing
	end

	if observableId == :S2IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 
		out[9] = (600.0S2_S4_S4 + 600.0ppS2_S4_S4 + 300.0ppS2_ppS2_S4 + 300.0ppS2_ppS2_ppS3 + 600.0ppS2_ppS3_S4 + 600.0ppS2_ppS3_ppS3) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[11] = (200.0S2_S4_S4 + 200.0ppS2_S4_S4 + 100.0ppS2_ppS2_S4 + 100.0ppS2_ppS2_ppS3 + 200.0ppS2_ppS3_S4 + 200.0ppS2_ppS3_ppS3) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[14] = (300.0S2_S4_S4 + 300.0ppS2_S4_S4 + 300.0ppS2_ppS3_S4 + 300.0ppS2_ppS3_ppS3 - 100.0S2 - 100.0pS2 - 100.0ppS2 - 300.0ppS2_ppS2_ppS2) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[19] = (200.0S2_S4_S4 + 200.0ppS2_S4_S4 + 100.0ppS2_ppS2_S4 + 100.0ppS2_ppS2_ppS3 + 200.0ppS2_ppS3_S4 + 200.0ppS2_ppS3_ppS3) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[20] = (-200.0S2 - 200.0pS2 - 200.0ppS2 - 300.0ppS2_ppS2_S4 - 300.0ppS2_ppS2_ppS3 - 600.0ppS2_ppS2_ppS2) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[26] = (-200.0S2 - 200.0pS2 - 200.0ppS2 - 300.0ppS2_ppS2_S4 - 300.0ppS2_ppS2_ppS3 - 600.0ppS2_ppS2_ppS2) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[28] = (-200.0S2 - 200.0pS2 - 200.0ppS2 - 300.0ppS2_ppS2_S4 - 300.0ppS2_ppS2_ppS3 - 600.0ppS2_ppS2_ppS2) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[31] = (300.0S2_S4_S4 + 300.0ppS2_S4_S4 + 300.0ppS2_ppS3_S4 + 300.0ppS2_ppS3_ppS3 - 100.0S2 - 100.0pS2 - 100.0ppS2 - 300.0ppS2_ppS2_ppS2) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[32] = (-200.0S2 - 200.0pS2 - 200.0ppS2 - 300.0ppS2_ppS2_S4 - 300.0ppS2_ppS2_ppS3 - 600.0ppS2_ppS2_ppS2) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[33] = (200.0S2_S4_S4 + 200.0ppS2_S4_S4 + 100.0ppS2_ppS2_S4 + 100.0ppS2_ppS2_ppS3 + 200.0ppS2_ppS3_S4 + 200.0ppS2_ppS3_ppS3) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		return nothing
	end

	if observableId == :S2IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		out[9] = (-300ppS2_ppS2_ppS3 - 300ppS2_ppS3_S4 - 600ppS2_ppS3_ppS3) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[11] = (-100ppS2_ppS2_ppS3 - 100ppS2_ppS3_S4 - 200ppS2_ppS3_ppS3) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[14] = (100.0S2 + 100.0pS2 + 300.0S2_S4_S4 + 100.0ppS2 + 300.0ppS2_S4_S4 + 300.0ppS2_ppS2_S4 + 300.0ppS2_ppS2_ppS2 - 300.0ppS2_ppS3_ppS3) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[19] = (-100ppS2_ppS2_ppS3 - 100ppS2_ppS3_S4 - 200ppS2_ppS3_ppS3) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[20] = (-300ppS2_ppS2_ppS3 - 300ppS2_ppS3_S4 - 600ppS2_ppS3_ppS3) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[26] = (-300ppS2_ppS2_ppS3 - 300ppS2_ppS3_S4 - 600ppS2_ppS3_ppS3) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[28] = (100.0S2 + 100.0pS2 + 300.0S2_S4_S4 + 100.0ppS2 + 300.0ppS2_S4_S4 + 300.0ppS2_ppS2_S4 + 300.0ppS2_ppS2_ppS2 - 300.0ppS2_ppS3_ppS3) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[31] = (-300ppS2_ppS2_ppS3 - 300ppS2_ppS3_S4 - 600ppS2_ppS3_ppS3) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[32] = (200.0S2 + 200.0pS2 + 600.0S2_S4_S4 + 200.0ppS2 + 300.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS3_S4 + 600.0ppS2_S4_S4 + 600.0ppS2_ppS2_S4 + 600.0ppS2_ppS2_ppS2) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[33] = (-100ppS2_ppS2_ppS3 - 100ppS2_ppS3_S4 - 200ppS2_ppS3_ppS3) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		return nothing
	end

	if observableId == :S2IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6 
		out[9] = (-300ppS2_ppS2_ppS3 - 300ppS2_ppS3_S4 - 600ppS2_ppS3_ppS3) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[11] = (-100ppS2_ppS2_ppS3 - 100ppS2_ppS3_S4 - 200ppS2_ppS3_ppS3) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[14] = (100.0S2 + 100.0pS2 + 300.0S2_S4_S4 + 100.0ppS2 + 300.0ppS2_S4_S4 + 300.0ppS2_ppS2_S4 + 300.0ppS2_ppS2_ppS2 - 300.0ppS2_ppS3_ppS3) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[19] = (-100ppS2_ppS2_ppS3 - 100ppS2_ppS3_S4 - 200ppS2_ppS3_ppS3) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[20] = (-300ppS2_ppS2_ppS3 - 300ppS2_ppS3_S4 - 600ppS2_ppS3_ppS3) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[26] = (-300ppS2_ppS2_ppS3 - 300ppS2_ppS3_S4 - 600ppS2_ppS3_ppS3) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[28] = (100.0S2 + 100.0pS2 + 300.0S2_S4_S4 + 100.0ppS2 + 300.0ppS2_S4_S4 + 300.0ppS2_ppS2_S4 + 300.0ppS2_ppS2_ppS2 - 300.0ppS2_ppS3_ppS3) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[31] = (-300ppS2_ppS2_ppS3 - 300ppS2_ppS3_S4 - 600ppS2_ppS3_ppS3) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[32] = (200.0S2 + 200.0pS2 + 600.0S2_S4_S4 + 200.0ppS2 + 300.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS3_S4 + 600.0ppS2_S4_S4 + 600.0ppS2_ppS2_S4 + 600.0ppS2_ppS2_ppS2) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[33] = (-100ppS2_ppS2_ppS3 - 100ppS2_ppS3_S4 - 200ppS2_ppS3_ppS3) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		return nothing
	end

	if observableId == :S2IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		out[9] = (-600ppS2_S4_S4 - 300ppS2_ppS2_S4 - 300ppS2_ppS3_S4) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[11] = (-200ppS2_S4_S4 - 100ppS2_ppS2_S4 - 100ppS2_ppS3_S4) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[14] = (-600ppS2_S4_S4 - 300ppS2_ppS2_S4 - 300ppS2_ppS3_S4) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[19] = (-200ppS2_S4_S4 - 100ppS2_ppS2_S4 - 100ppS2_ppS3_S4) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[20] = (-600ppS2_S4_S4 - 300ppS2_ppS2_S4 - 300ppS2_ppS3_S4) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[26] = (200.0S2 + 200.0pS2 + 600.0S2_S4_S4 + 200.0ppS2 + 300.0ppS2_ppS2_S4 + 300.0ppS2_ppS3_S4 + 600.0ppS2_ppS2_ppS2 + 600.0ppS2_ppS2_ppS3 + 600.0ppS2_ppS3_ppS3) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[28] = (100.0S2 + 100.0pS2 + 300.0S2_S4_S4 + 100.0ppS2 + 300.0ppS2_ppS2_ppS2 + 300.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS3_ppS3 - 300.0ppS2_S4_S4) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[31] = (100.0S2 + 100.0pS2 + 300.0S2_S4_S4 + 100.0ppS2 + 300.0ppS2_ppS2_ppS2 + 300.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS3_ppS3 - 300.0ppS2_S4_S4) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[32] = (-600ppS2_S4_S4 - 300ppS2_ppS2_S4 - 300ppS2_ppS3_S4) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[33] = (-200ppS2_S4_S4 - 100ppS2_ppS2_S4 - 100ppS2_ppS3_S4) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		return nothing
	end

	if observableId == :S2IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6 
		out[9] = (-600ppS2_S4_S4 - 300ppS2_ppS2_S4 - 300ppS2_ppS3_S4) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[11] = (-200ppS2_S4_S4 - 100ppS2_ppS2_S4 - 100ppS2_ppS3_S4) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[14] = (-600ppS2_S4_S4 - 300ppS2_ppS2_S4 - 300ppS2_ppS3_S4) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[19] = (-200ppS2_S4_S4 - 100ppS2_ppS2_S4 - 100ppS2_ppS3_S4) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[20] = (-600ppS2_S4_S4 - 300ppS2_ppS2_S4 - 300ppS2_ppS3_S4) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[26] = (200.0S2 + 200.0pS2 + 600.0S2_S4_S4 + 200.0ppS2 + 300.0ppS2_ppS2_S4 + 300.0ppS2_ppS3_S4 + 600.0ppS2_ppS2_ppS2 + 600.0ppS2_ppS2_ppS3 + 600.0ppS2_ppS3_ppS3) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[28] = (100.0S2 + 100.0pS2 + 300.0S2_S4_S4 + 100.0ppS2 + 300.0ppS2_ppS2_ppS2 + 300.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS3_ppS3 - 300.0ppS2_S4_S4) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[31] = (100.0S2 + 100.0pS2 + 300.0S2_S4_S4 + 100.0ppS2 + 300.0ppS2_ppS2_ppS2 + 300.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS3_ppS3 - 300.0ppS2_S4_S4) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[32] = (-600ppS2_S4_S4 - 300ppS2_ppS2_S4 - 300ppS2_ppS3_S4) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		out[33] = (-200ppS2_S4_S4 - 100ppS2_ppS2_S4 - 100ppS2_ppS3_S4) / ((S2 + pS2 + ppS2 + 3S2_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS2_ppS2 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3)^2)
		return nothing
	end

	if observableId == :S3IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S3IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		out[1] = (-300S3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[2] = (-100S3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[12] = (-100S3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[14] = (-100S3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[18] = (100.0pS3 + 100.0ppS2_ppS2_ppS3 + 100.0ppS2_ppS3_S4 + 200.0ppS2_ppS3_ppS3 + 100.0ppS3 + 100.0ppS3_S4_S4 + 200.0ppS3_ppS3_S4 + 300.0ppS3_ppS3_ppS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[21] = (-200S3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[23] = (-100S3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[28] = (-100S3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[32] = (-200S3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		return nothing
	end

	if observableId == :S3IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S3IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		out[1] = (-300pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[2] = (-100pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[12] = (100.0S3 + 100.0ppS2_ppS2_ppS3 + 100.0ppS2_ppS3_S4 + 200.0ppS2_ppS3_ppS3 + 100.0ppS3 + 100.0ppS3_S4_S4 + 200.0ppS3_ppS3_S4 + 300.0ppS3_ppS3_ppS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[14] = (-100pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[18] = (-100pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[21] = (-200pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[23] = (-100pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[28] = (-100pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[32] = (-200pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		return nothing
	end

	if observableId == :S3IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S3IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		out[1] = (300.0S3 + 300.0pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[2] = (100.0S3 + 100.0pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[12] = (-100ppS2_ppS2_ppS3 - 100ppS2_ppS3_S4 - 200ppS2_ppS3_ppS3 - 100ppS3 - 100ppS3_S4_S4 - 200ppS3_ppS3_S4 - 300ppS3_ppS3_ppS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[14] = (100.0S3 + 100.0pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[18] = (-100ppS2_ppS2_ppS3 - 100ppS2_ppS3_S4 - 200ppS2_ppS3_ppS3 - 100ppS3 - 100ppS3_S4_S4 - 200ppS3_ppS3_S4 - 300ppS3_ppS3_ppS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[21] = (200.0S3 + 200.0pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[23] = (100.0S3 + 100.0pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[28] = (100.0S3 + 100.0pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[32] = (200.0S3 + 200.0pS3) / ((S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS3 + ppS3_S4_S4 + 2ppS2_ppS3_ppS3 + 2ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		return nothing
	end

	if observableId == :S3IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		out[1] = (-600ppS2_ppS2_ppS3 - 300ppS2_ppS3_S4 - 300ppS2_ppS3_ppS3) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[2] = (-600ppS2_ppS2_ppS3 - 300ppS2_ppS3_S4 - 300ppS2_ppS3_ppS3) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[12] = (-200ppS2_ppS2_ppS3 - 100ppS2_ppS3_S4 - 100ppS2_ppS3_ppS3) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[14] = (200.0S3 + 200.0pS3 + 200.0ppS3 + 300.0ppS2_ppS3_S4 + 300.0ppS2_ppS3_ppS3 + 600.0ppS3_S4_S4 + 600.0ppS3_ppS3_S4 + 600.0ppS3_ppS3_ppS3) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[18] = (-200ppS2_ppS2_ppS3 - 100ppS2_ppS3_S4 - 100ppS2_ppS3_ppS3) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[21] = (-600ppS2_ppS2_ppS3 - 300ppS2_ppS3_S4 - 300ppS2_ppS3_ppS3) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[23] = (-200ppS2_ppS2_ppS3 - 100ppS2_ppS3_S4 - 100ppS2_ppS3_ppS3) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[28] = (100.0S3 + 100.0pS3 + 100.0ppS3 + 300.0ppS3_S4_S4 + 300.0ppS3_ppS3_S4 + 300.0ppS3_ppS3_ppS3 - 300.0ppS2_ppS2_ppS3) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[32] = (100.0S3 + 100.0pS3 + 100.0ppS3 + 300.0ppS3_S4_S4 + 300.0ppS3_ppS3_S4 + 300.0ppS3_ppS3_ppS3 - 300.0ppS2_ppS2_ppS3) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		return nothing
	end

	if observableId == :S3IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 
		out[1] = (-600ppS2_ppS2_ppS3 - 300ppS2_ppS3_S4 - 300ppS2_ppS3_ppS3) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[2] = (-600ppS2_ppS2_ppS3 - 300ppS2_ppS3_S4 - 300ppS2_ppS3_ppS3) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[12] = (-200ppS2_ppS2_ppS3 - 100ppS2_ppS3_S4 - 100ppS2_ppS3_ppS3) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[14] = (200.0S3 + 200.0pS3 + 200.0ppS3 + 300.0ppS2_ppS3_S4 + 300.0ppS2_ppS3_ppS3 + 600.0ppS3_S4_S4 + 600.0ppS3_ppS3_S4 + 600.0ppS3_ppS3_ppS3) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[18] = (-200ppS2_ppS2_ppS3 - 100ppS2_ppS3_S4 - 100ppS2_ppS3_ppS3) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[21] = (-600ppS2_ppS2_ppS3 - 300ppS2_ppS3_S4 - 300ppS2_ppS3_ppS3) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[23] = (-200ppS2_ppS2_ppS3 - 100ppS2_ppS3_S4 - 100ppS2_ppS3_ppS3) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[28] = (100.0S3 + 100.0pS3 + 100.0ppS3 + 300.0ppS3_S4_S4 + 300.0ppS3_ppS3_S4 + 300.0ppS3_ppS3_ppS3 - 300.0ppS2_ppS2_ppS3) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[32] = (100.0S3 + 100.0pS3 + 100.0ppS3 + 300.0ppS3_S4_S4 + 300.0ppS3_ppS3_S4 + 300.0ppS3_ppS3_ppS3 - 300.0ppS2_ppS2_ppS3) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		return nothing
	end

	if observableId == :S3IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		out[1] = (600.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS3_ppS3 + 600.0ppS2_ppS3_S4 + 600.0ppS3_S4_S4 + 300.0ppS3_ppS3_S4) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[2] = (-200.0S3 - 200.0pS3 - 300.0ppS2_ppS3_ppS3 - 200.0ppS3 - 300.0ppS3_ppS3_S4 - 600.0ppS3_ppS3_ppS3) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[12] = (200.0ppS2_ppS2_ppS3 + 100.0ppS2_ppS3_ppS3 + 200.0ppS2_ppS3_S4 + 200.0ppS3_S4_S4 + 100.0ppS3_ppS3_S4) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[14] = (-200.0S3 - 200.0pS3 - 300.0ppS2_ppS3_ppS3 - 200.0ppS3 - 300.0ppS3_ppS3_S4 - 600.0ppS3_ppS3_ppS3) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[18] = (200.0ppS2_ppS2_ppS3 + 100.0ppS2_ppS3_ppS3 + 200.0ppS2_ppS3_S4 + 200.0ppS3_S4_S4 + 100.0ppS3_ppS3_S4) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[21] = (300.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS3_S4 + 300.0ppS3_S4_S4 - 100.0S3 - 100.0pS3 - 100.0ppS3 - 300.0ppS3_ppS3_ppS3) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[23] = (200.0ppS2_ppS2_ppS3 + 100.0ppS2_ppS3_ppS3 + 200.0ppS2_ppS3_S4 + 200.0ppS3_S4_S4 + 100.0ppS3_ppS3_S4) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[28] = (-200.0S3 - 200.0pS3 - 300.0ppS2_ppS3_ppS3 - 200.0ppS3 - 300.0ppS3_ppS3_S4 - 600.0ppS3_ppS3_ppS3) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[32] = (300.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS3_S4 + 300.0ppS3_S4_S4 - 100.0S3 - 100.0pS3 - 100.0ppS3 - 300.0ppS3_ppS3_ppS3) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		return nothing
	end

	if observableId == :S3IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6 
		out[1] = (600.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS3_ppS3 + 600.0ppS2_ppS3_S4 + 600.0ppS3_S4_S4 + 300.0ppS3_ppS3_S4) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[2] = (-200.0S3 - 200.0pS3 - 300.0ppS2_ppS3_ppS3 - 200.0ppS3 - 300.0ppS3_ppS3_S4 - 600.0ppS3_ppS3_ppS3) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[12] = (200.0ppS2_ppS2_ppS3 + 100.0ppS2_ppS3_ppS3 + 200.0ppS2_ppS3_S4 + 200.0ppS3_S4_S4 + 100.0ppS3_ppS3_S4) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[14] = (-200.0S3 - 200.0pS3 - 300.0ppS2_ppS3_ppS3 - 200.0ppS3 - 300.0ppS3_ppS3_S4 - 600.0ppS3_ppS3_ppS3) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[18] = (200.0ppS2_ppS2_ppS3 + 100.0ppS2_ppS3_ppS3 + 200.0ppS2_ppS3_S4 + 200.0ppS3_S4_S4 + 100.0ppS3_ppS3_S4) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[21] = (300.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS3_S4 + 300.0ppS3_S4_S4 - 100.0S3 - 100.0pS3 - 100.0ppS3 - 300.0ppS3_ppS3_ppS3) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[23] = (200.0ppS2_ppS2_ppS3 + 100.0ppS2_ppS3_ppS3 + 200.0ppS2_ppS3_S4 + 200.0ppS3_S4_S4 + 100.0ppS3_ppS3_S4) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[28] = (-200.0S3 - 200.0pS3 - 300.0ppS2_ppS3_ppS3 - 200.0ppS3 - 300.0ppS3_ppS3_S4 - 600.0ppS3_ppS3_ppS3) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[32] = (300.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS3_S4 + 300.0ppS3_S4_S4 - 100.0S3 - 100.0pS3 - 100.0ppS3 - 300.0ppS3_ppS3_ppS3) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		return nothing
	end

	if observableId == :S3IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		out[1] = (-300ppS2_ppS3_S4 - 600ppS3_S4_S4 - 300ppS3_ppS3_S4) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[2] = (200.0S3 + 200.0pS3 + 600.0ppS2_ppS2_ppS3 + 200.0ppS3 + 300.0ppS2_ppS3_S4 + 600.0ppS2_ppS3_ppS3 + 300.0ppS3_ppS3_S4 + 600.0ppS3_ppS3_ppS3) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[12] = (-100ppS2_ppS3_S4 - 200ppS3_S4_S4 - 100ppS3_ppS3_S4) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[14] = (-300ppS2_ppS3_S4 - 600ppS3_S4_S4 - 300ppS3_ppS3_S4) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[18] = (-100ppS2_ppS3_S4 - 200ppS3_S4_S4 - 100ppS3_ppS3_S4) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[21] = (100.0S3 + 100.0pS3 + 300.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS3_ppS3 + 100.0ppS3 + 300.0ppS3_ppS3_ppS3 - 300.0ppS3_S4_S4) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[23] = (-100ppS2_ppS3_S4 - 200ppS3_S4_S4 - 100ppS3_ppS3_S4) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[28] = (100.0S3 + 100.0pS3 + 300.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS3_ppS3 + 100.0ppS3 + 300.0ppS3_ppS3_ppS3 - 300.0ppS3_S4_S4) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[32] = (-300ppS2_ppS3_S4 - 600ppS3_S4_S4 - 300ppS3_ppS3_S4) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		return nothing
	end

	if observableId == :S3IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6 
		out[1] = (-300ppS2_ppS3_S4 - 600ppS3_S4_S4 - 300ppS3_ppS3_S4) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[2] = (200.0S3 + 200.0pS3 + 600.0ppS2_ppS2_ppS3 + 200.0ppS3 + 300.0ppS2_ppS3_S4 + 600.0ppS2_ppS3_ppS3 + 300.0ppS3_ppS3_S4 + 600.0ppS3_ppS3_ppS3) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[12] = (-100ppS2_ppS3_S4 - 200ppS3_S4_S4 - 100ppS3_ppS3_S4) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[14] = (-300ppS2_ppS3_S4 - 600ppS3_S4_S4 - 300ppS3_ppS3_S4) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[18] = (-100ppS2_ppS3_S4 - 200ppS3_S4_S4 - 100ppS3_ppS3_S4) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[21] = (100.0S3 + 100.0pS3 + 300.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS3_ppS3 + 100.0ppS3 + 300.0ppS3_ppS3_ppS3 - 300.0ppS3_S4_S4) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[23] = (-100ppS2_ppS3_S4 - 200ppS3_S4_S4 - 100ppS3_ppS3_S4) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[28] = (100.0S3 + 100.0pS3 + 300.0ppS2_ppS2_ppS3 + 300.0ppS2_ppS3_ppS3 + 100.0ppS3 + 300.0ppS3_ppS3_ppS3 - 300.0ppS3_S4_S4) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		out[32] = (-300ppS2_ppS3_S4 - 600ppS3_S4_S4 - 300ppS3_ppS3_S4) / ((S3 + pS3 + ppS3 + 3ppS2_ppS2_ppS3 + 3ppS2_ppS3_S4 + 3ppS2_ppS3_ppS3 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4 + 3ppS3_ppS3_ppS3)^2)
		return nothing
	end

	if observableId == :S4IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		out[2] = (-300S2_S4_S4 - 300ppS2_S4_S4 - 300ppS2_ppS3_S4 - 600ppS2_ppS2_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[8] = (-100S2_S4_S4 - 100ppS2_S4_S4 - 100ppS2_ppS3_S4 - 200ppS2_ppS2_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[16] = (-300S2_S4_S4 - 300ppS2_S4_S4 - 300ppS2_ppS3_S4 - 600ppS2_ppS2_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[20] = (100.0S4 + 300.0S4_S4_S4 + 300.0ppS3_S4_S4 + 300.0ppS3_ppS3_S4 - 300.0ppS2_ppS2_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[21] = (-300S2_S4_S4 - 300ppS2_S4_S4 - 300ppS2_ppS3_S4 - 600ppS2_ppS2_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[26] = (100.0S4 + 300.0S4_S4_S4 + 300.0ppS3_S4_S4 + 300.0ppS3_ppS3_S4 - 300.0ppS2_ppS2_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[28] = (100.0S4 + 300.0S4_S4_S4 + 300.0ppS3_S4_S4 + 300.0ppS3_ppS3_S4 - 300.0ppS2_ppS2_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[31] = (300.0S2_S4_S4 + 200.0S4 + 600.0S4_S4_S4 + 300.0ppS2_S4_S4 + 300.0ppS2_ppS3_S4 + 600.0ppS3_S4_S4 + 600.0ppS3_ppS3_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		return nothing
	end

	if observableId == :S4IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 
		out[2] = (-300S2_S4_S4 - 300ppS2_S4_S4 - 300ppS2_ppS3_S4 - 600ppS2_ppS2_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[8] = (-100S2_S4_S4 - 100ppS2_S4_S4 - 100ppS2_ppS3_S4 - 200ppS2_ppS2_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[16] = (-300S2_S4_S4 - 300ppS2_S4_S4 - 300ppS2_ppS3_S4 - 600ppS2_ppS2_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[20] = (100.0S4 + 300.0S4_S4_S4 + 300.0ppS3_S4_S4 + 300.0ppS3_ppS3_S4 - 300.0ppS2_ppS2_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[21] = (-300S2_S4_S4 - 300ppS2_S4_S4 - 300ppS2_ppS3_S4 - 600ppS2_ppS2_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[26] = (100.0S4 + 300.0S4_S4_S4 + 300.0ppS3_S4_S4 + 300.0ppS3_ppS3_S4 - 300.0ppS2_ppS2_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[28] = (100.0S4 + 300.0S4_S4_S4 + 300.0ppS3_S4_S4 + 300.0ppS3_ppS3_S4 - 300.0ppS2_ppS2_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[31] = (300.0S2_S4_S4 + 200.0S4 + 600.0S4_S4_S4 + 300.0ppS2_S4_S4 + 300.0ppS2_ppS3_S4 + 600.0ppS3_S4_S4 + 600.0ppS3_ppS3_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		return nothing
	end

	if observableId == :S4IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		out[2] = (300.0S2_S4_S4 + 100.0S4 + 300.0S4_S4_S4 + 300.0ppS2_S4_S4 + 300.0ppS2_ppS2_S4 - 300.0ppS3_ppS3_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[8] = (-100ppS2_ppS3_S4 - 100ppS3_S4_S4 - 200ppS3_ppS3_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[16] = (-300ppS2_ppS3_S4 - 300ppS3_S4_S4 - 600ppS3_ppS3_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[20] = (-300ppS2_ppS3_S4 - 300ppS3_S4_S4 - 600ppS3_ppS3_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[21] = (600.0S2_S4_S4 + 200.0S4 + 300.0ppS2_ppS3_S4 + 300.0ppS3_S4_S4 + 600.0S4_S4_S4 + 600.0ppS2_S4_S4 + 600.0ppS2_ppS2_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[26] = (-300ppS2_ppS3_S4 - 300ppS3_S4_S4 - 600ppS3_ppS3_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[28] = (300.0S2_S4_S4 + 100.0S4 + 300.0S4_S4_S4 + 300.0ppS2_S4_S4 + 300.0ppS2_ppS2_S4 - 300.0ppS3_ppS3_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[31] = (-300ppS2_ppS3_S4 - 300ppS3_S4_S4 - 600ppS3_ppS3_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		return nothing
	end

	if observableId == :S4IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6 
		out[2] = (300.0S2_S4_S4 + 100.0S4 + 300.0S4_S4_S4 + 300.0ppS2_S4_S4 + 300.0ppS2_ppS2_S4 - 300.0ppS3_ppS3_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[8] = (-100ppS2_ppS3_S4 - 100ppS3_S4_S4 - 200ppS3_ppS3_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[16] = (-300ppS2_ppS3_S4 - 300ppS3_S4_S4 - 600ppS3_ppS3_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[20] = (-300ppS2_ppS3_S4 - 300ppS3_S4_S4 - 600ppS3_ppS3_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[21] = (600.0S2_S4_S4 + 200.0S4 + 300.0ppS2_ppS3_S4 + 300.0ppS3_S4_S4 + 600.0S4_S4_S4 + 600.0ppS2_S4_S4 + 600.0ppS2_ppS2_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[26] = (-300ppS2_ppS3_S4 - 300ppS3_S4_S4 - 600ppS3_ppS3_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[28] = (300.0S2_S4_S4 + 100.0S4 + 300.0S4_S4_S4 + 300.0ppS2_S4_S4 + 300.0ppS2_ppS2_S4 - 300.0ppS3_ppS3_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[31] = (-300ppS2_ppS3_S4 - 300ppS3_S4_S4 - 600ppS3_ppS3_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		return nothing
	end

	if observableId == :S4IP_tS4_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1 
		out[2] = 2
		out[8] = 1
		out[16] = 3
		out[20] = 2
		out[21] = 1
		out[26] = 2
		out[28] = 1
		out[31] = 1
		return nothing
	end

	if observableId == :S4IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		out[2] = (300.0ppS2_ppS2_S4 + 300.0ppS2_ppS3_S4 + 300.0ppS3_ppS3_S4 - 100.0S4 - 300.0S4_S4_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[8] = (100.0S2_S4_S4 + 100.0ppS2_S4_S4 + 100.0ppS3_S4_S4 + 200.0ppS2_ppS2_S4 + 200.0ppS2_ppS3_S4 + 200.0ppS3_ppS3_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[16] = (300.0S2_S4_S4 + 300.0ppS2_S4_S4 + 300.0ppS3_S4_S4 + 600.0ppS2_ppS2_S4 + 600.0ppS2_ppS3_S4 + 600.0ppS3_ppS3_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[20] = (300.0ppS2_ppS2_S4 + 300.0ppS2_ppS3_S4 + 300.0ppS3_ppS3_S4 - 100.0S4 - 300.0S4_S4_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[21] = (-300.0S2_S4_S4 - 200.0S4 - 600.0S4_S4_S4 - 300.0ppS2_S4_S4 - 300.0ppS3_S4_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[26] = (300.0ppS2_ppS2_S4 + 300.0ppS2_ppS3_S4 + 300.0ppS3_ppS3_S4 - 100.0S4 - 300.0S4_S4_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[28] = (-300.0S2_S4_S4 - 200.0S4 - 600.0S4_S4_S4 - 300.0ppS2_S4_S4 - 300.0ppS3_S4_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[31] = (-300.0S2_S4_S4 - 200.0S4 - 600.0S4_S4_S4 - 300.0ppS2_S4_S4 - 300.0ppS3_S4_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		return nothing
	end

	if observableId == :S4IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6 
		out[2] = (300.0ppS2_ppS2_S4 + 300.0ppS2_ppS3_S4 + 300.0ppS3_ppS3_S4 - 100.0S4 - 300.0S4_S4_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[8] = (100.0S2_S4_S4 + 100.0ppS2_S4_S4 + 100.0ppS3_S4_S4 + 200.0ppS2_ppS2_S4 + 200.0ppS2_ppS3_S4 + 200.0ppS3_ppS3_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[16] = (300.0S2_S4_S4 + 300.0ppS2_S4_S4 + 300.0ppS3_S4_S4 + 600.0ppS2_ppS2_S4 + 600.0ppS2_ppS3_S4 + 600.0ppS3_ppS3_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[20] = (300.0ppS2_ppS2_S4 + 300.0ppS2_ppS3_S4 + 300.0ppS3_ppS3_S4 - 100.0S4 - 300.0S4_S4_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[21] = (-300.0S2_S4_S4 - 200.0S4 - 600.0S4_S4_S4 - 300.0ppS2_S4_S4 - 300.0ppS3_S4_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[26] = (300.0ppS2_ppS2_S4 + 300.0ppS2_ppS3_S4 + 300.0ppS3_ppS3_S4 - 100.0S4 - 300.0S4_S4_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[28] = (-300.0S2_S4_S4 - 200.0S4 - 600.0S4_S4_S4 - 300.0ppS2_S4_S4 - 300.0ppS3_S4_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		out[31] = (-300.0S2_S4_S4 - 200.0S4 - 600.0S4_S4_S4 - 300.0ppS2_S4_S4 - 300.0ppS3_S4_S4) / ((S4 + 3S2_S4_S4 + 3S4_S4_S4 + 3ppS2_S4_S4 + 3ppS2_ppS2_S4 + 3ppS2_ppS3_S4 + 3ppS3_S4_S4 + 3ppS3_ppS3_S4)^2)
		return nothing
	end

	if observableId == :Ski_gene_expression_OE_nExpID100 
		out[24] = 1
		return nothing
	end

	if observableId == :Skil_gene_expression_OE_nExpID100 
		out[22] = 1
		return nothing
	end

	if observableId == :Smad7_gene_expression_OE_nExpID100 
		out[15] = 1
		return nothing
	end

	if observableId == :Sox4_gene_expression_OE_nExpID100 
		out[25] = 1
		return nothing
	end

	if observableId == :TGFbR_PL_PN_TGFbR_hepa16_nExpID11 
		out[6] = 1
		out[30] = 1
		return nothing
	end

	if observableId == :Tgfa_gene_expression_OE_nExpID100 
		out[10] = 1
		return nothing
	end

end

function evalDYmodDp(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	ppS3_ppS3_ppS3, ppS3_S4_S4, geneH, geneI, geneJ, TGFb_pRec, geneE, S4, ppS2_ppS2_ppS2, geneK, pS2, pS3, geneC, ppS2_ppS2_ppS3, geneF, S4_S4_S4, TGFb, S3, S2, S2_S4_S4, ppS3_ppS3_S4, geneB, ppS3, geneA, geneD, ppS2_S4_S4, geneL, ppS2_ppS3_S4, geneG, Rec, ppS2_ppS2_S4, ppS2_ppS3_ppS3, ppS2= u 
	geneC_inh3, geneJ_act1, geneG_inh1, init_Rec, geneC_act3, geneA_inh2, geneK_turn, geneA_act3, geneI_act2, geneA_act2, geneD_inh2, k_223, k_233, geneB_act1, S_dephosphos, geneG_act1, geneH_inh2, geneI_act1, geneB_inh3, geneI_inh2, geneC_turn, geneJ_inh1, geneK_act3, geneJ_inh3, geneG_inh3, geneC_act2, init_TGFb, geneL_inh1, geneE_inh3, geneF_act3, cell, geneF_inh1, geneD_act3, geneE_act2, geneD_inh3, k_234, geneH_act2, geneA_turn, geneL_act1, geneH_act1, geneL_inh2, geneB_turn, init_S4, geneL_inh3, khomo2, geneK_act2, geneA_inh3, geneJ_act2, geneI_inh3, geneD_act1, geneJ_turn, geneG_act3, geneL_act2, Rec_act, geneI_act3, k_224, pRec_degind, geneE_inh2, geneF_inh2, geneC_act1, geneD_inh1, k_344, init_S3, geneI_inh1, geneL_turn, geneF_act2, k_on_u, geneE_inh1, geneH_act3, geneH_inh1, geneK_inh3, geneE_act1, geneA_inh1, geneB_inh1, geneH_inh3, geneD_turn, S_dephos, geneL_act3, S_phos, geneG_turn, geneA_act1, geneI_turn, khomo3, geneC_inh1, geneG_act2, init_S2, geneK_inh2, k_334, geneB_inh2, geneH_turn, geneJ_inh2, khomo4, geneB_act2, geneF_turn, geneJ_act3, geneK_inh1, geneB_act3, geneK_act1, geneE_turn, geneE_act3, geneG_inh2, kdiss_SS, k_244, geneD_act2, geneF_inh3, geneF_act1, geneC_inh2 = p 
	sd_Bmp4_nExpID100, sd_Cxcl15_nExpID100, sd_Dnmt3a_nExpID100, sd_Dusp5_nExpID100, sd_Jun_nExpID100, sd_Klf10_nExpID100, sd_Pdk4_nExpID100, sd_Ski_nExpID100, sd_Skil_nExpID100, sd_Smad7_nExpID100, sd_Sox4_nExpID100, sd_Tgfa_nExpID100 = nonDynParam 
	if observableId == :Bmp4_gene_expression_OE_nExpID100 
		return nothing
	end

	if observableId == :Cxcl15_gene_expression_OE_nExpID100 
		return nothing
	end

	if observableId == :Dnmt3a_gene_expression_OE_nExpID100 
		return nothing
	end

	if observableId == :Dusp5_gene_expression_OE_nExpID100 
		return nothing
	end

	if observableId == :Jun_gene_expression_OE_nExpID100 
		return nothing
	end

	if observableId == :Klf10_gene_expression_OE_nExpID100 
		return nothing
	end

	if observableId == :Pdk4_gene_expression_OE_nExpID100 
		return nothing
	end

	if observableId == :S23IP_S2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 
		return nothing
	end

	if observableId == :S23IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S23IP_S3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 
		return nothing
	end

	if observableId == :S23IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S23IP_pS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 
		return nothing
	end

	if observableId == :S23IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S23IP_pS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 
		return nothing
	end

	if observableId == :S23IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S23IP_ppS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 
		return nothing
	end

	if observableId == :S23IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S23IP_ppS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 
		return nothing
	end

	if observableId == :S23IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1 
		return nothing
	end

	if observableId == :S23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10 
		return nothing
	end

	if observableId == :S23IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S23IP_tS3_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10 
		return nothing
	end

	if observableId == :S23IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S23IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S2IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S2IP_S2_MSPL2_complexes_all_hepa16_time_course_nExpID6 
		return nothing
	end

	if observableId == :S2IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S2IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S2IP_pS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 
		return nothing
	end

	if observableId == :S2IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S2IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S2IP_ppS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 
		return nothing
	end

	if observableId == :S2IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S2IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S2IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 
		return nothing
	end

	if observableId == :S2IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S2IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6 
		return nothing
	end

	if observableId == :S2IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S2IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6 
		return nothing
	end

	if observableId == :S3IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S3IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S3IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S3IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S3IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S3IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S3IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S3IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 
		return nothing
	end

	if observableId == :S3IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S3IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6 
		return nothing
	end

	if observableId == :S3IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S3IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6 
		return nothing
	end

	if observableId == :S4IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S4IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 
		return nothing
	end

	if observableId == :S4IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S4IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6 
		return nothing
	end

	if observableId == :S4IP_tS4_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1 
		return nothing
	end

	if observableId == :S4IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 
		return nothing
	end

	if observableId == :S4IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6 
		return nothing
	end

	if observableId == :Ski_gene_expression_OE_nExpID100 
		return nothing
	end

	if observableId == :Skil_gene_expression_OE_nExpID100 
		return nothing
	end

	if observableId == :Smad7_gene_expression_OE_nExpID100 
		return nothing
	end

	if observableId == :Sox4_gene_expression_OE_nExpID100 
		return nothing
	end

	if observableId == :TGFbR_PL_PN_TGFbR_hepa16_nExpID11 
		return nothing
	end

	if observableId == :Tgfa_gene_expression_OE_nExpID100 
		return nothing
	end

end

function evalDSdDu!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	ppS3_ppS3_ppS3, ppS3_S4_S4, geneH, geneI, geneJ, TGFb_pRec, geneE, S4, ppS2_ppS2_ppS2, geneK, pS2, pS3, geneC, ppS2_ppS2_ppS3, geneF, S4_S4_S4, TGFb, S3, S2, S2_S4_S4, ppS3_ppS3_S4, geneB, ppS3, geneA, geneD, ppS2_S4_S4, geneL, ppS2_ppS3_S4, geneG, Rec, ppS2_ppS2_S4, ppS2_ppS3_ppS3, ppS2= u 
	geneC_inh3, geneJ_act1, geneG_inh1, init_Rec, geneC_act3, geneA_inh2, geneK_turn, geneA_act3, geneI_act2, geneA_act2, geneD_inh2, k_223, k_233, geneB_act1, S_dephosphos, geneG_act1, geneH_inh2, geneI_act1, geneB_inh3, geneI_inh2, geneC_turn, geneJ_inh1, geneK_act3, geneJ_inh3, geneG_inh3, geneC_act2, init_TGFb, geneL_inh1, geneE_inh3, geneF_act3, cell, geneF_inh1, geneD_act3, geneE_act2, geneD_inh3, k_234, geneH_act2, geneA_turn, geneL_act1, geneH_act1, geneL_inh2, geneB_turn, init_S4, geneL_inh3, khomo2, geneK_act2, geneA_inh3, geneJ_act2, geneI_inh3, geneD_act1, geneJ_turn, geneG_act3, geneL_act2, Rec_act, geneI_act3, k_224, pRec_degind, geneE_inh2, geneF_inh2, geneC_act1, geneD_inh1, k_344, init_S3, geneI_inh1, geneL_turn, geneF_act2, k_on_u, geneE_inh1, geneH_act3, geneH_inh1, geneK_inh3, geneE_act1, geneA_inh1, geneB_inh1, geneH_inh3, geneD_turn, S_dephos, geneL_act3, S_phos, geneG_turn, geneA_act1, geneI_turn, khomo3, geneC_inh1, geneG_act2, init_S2, geneK_inh2, k_334, geneB_inh2, geneH_turn, geneJ_inh2, khomo4, geneB_act2, geneF_turn, geneJ_act3, geneK_inh1, geneB_act3, geneK_act1, geneE_turn, geneE_act3, geneG_inh2, kdiss_SS, k_244, geneD_act2, geneF_inh3, geneF_act1, geneC_inh2 = p 
	sd_Bmp4_nExpID100, sd_Cxcl15_nExpID100, sd_Dnmt3a_nExpID100, sd_Dusp5_nExpID100, sd_Jun_nExpID100, sd_Klf10_nExpID100, sd_Pdk4_nExpID100, sd_Ski_nExpID100, sd_Skil_nExpID100, sd_Smad7_nExpID100, sd_Sox4_nExpID100, sd_Tgfa_nExpID100 = nonDynParam 
	if observableId == "Bmp4_gene_expression_OE_nExpID100" 
		return nothing
	end

	if observableId == "Cxcl15_gene_expression_OE_nExpID100" 
		return nothing
	end

	if observableId == "Dnmt3a_gene_expression_OE_nExpID100" 
		return nothing
	end

	if observableId == "Dusp5_gene_expression_OE_nExpID100" 
		return nothing
	end

	if observableId == "Jun_gene_expression_OE_nExpID100" 
		return nothing
	end

	if observableId == "Klf10_gene_expression_OE_nExpID100" 
		return nothing
	end

	if observableId == "Pdk4_gene_expression_OE_nExpID100" 
		return nothing
	end

	if observableId == "S23IP_S2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5" 
		return nothing
	end

	if observableId == "S23IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S23IP_S3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5" 
		return nothing
	end

	if observableId == "S23IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S23IP_pS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5" 
		return nothing
	end

	if observableId == "S23IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S23IP_pS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5" 
		return nothing
	end

	if observableId == "S23IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S23IP_ppS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5" 
		return nothing
	end

	if observableId == "S23IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S23IP_ppS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5" 
		return nothing
	end

	if observableId == "S23IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1" 
		return nothing
	end

	if observableId == "S23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10" 
		return nothing
	end

	if observableId == "S23IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S23IP_tS3_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10" 
		return nothing
	end

	if observableId == "S23IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S23IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S2IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S2IP_S2_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return nothing
	end

	if observableId == "S2IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S2IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S2IP_pS2_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return nothing
	end

	if observableId == "S2IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S2IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S2IP_ppS2_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return nothing
	end

	if observableId == "S2IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S2IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S2IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return nothing
	end

	if observableId == "S2IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S2IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return nothing
	end

	if observableId == "S2IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S2IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return nothing
	end

	if observableId == "S3IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S3IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S3IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S3IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S3IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S3IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S3IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S3IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return nothing
	end

	if observableId == "S3IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S3IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return nothing
	end

	if observableId == "S3IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S3IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return nothing
	end

	if observableId == "S4IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S4IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return nothing
	end

	if observableId == "S4IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S4IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return nothing
	end

	if observableId == "S4IP_tS4_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1" 
		return nothing
	end

	if observableId == "S4IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S4IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return nothing
	end

	if observableId == "Ski_gene_expression_OE_nExpID100" 
		return nothing
	end

	if observableId == "Skil_gene_expression_OE_nExpID100" 
		return nothing
	end

	if observableId == "Smad7_gene_expression_OE_nExpID100" 
		return nothing
	end

	if observableId == "Sox4_gene_expression_OE_nExpID100" 
		return nothing
	end

	if observableId == "TGFbR_PL_PN_TGFbR_hepa16_nExpID11" 
		return nothing
	end

	if observableId == "Tgfa_gene_expression_OE_nExpID100" 
		return nothing
	end

end

function evalDSdDp!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	ppS3_ppS3_ppS3, ppS3_S4_S4, geneH, geneI, geneJ, TGFb_pRec, geneE, S4, ppS2_ppS2_ppS2, geneK, pS2, pS3, geneC, ppS2_ppS2_ppS3, geneF, S4_S4_S4, TGFb, S3, S2, S2_S4_S4, ppS3_ppS3_S4, geneB, ppS3, geneA, geneD, ppS2_S4_S4, geneL, ppS2_ppS3_S4, geneG, Rec, ppS2_ppS2_S4, ppS2_ppS3_ppS3, ppS2= u 
	geneC_inh3, geneJ_act1, geneG_inh1, init_Rec, geneC_act3, geneA_inh2, geneK_turn, geneA_act3, geneI_act2, geneA_act2, geneD_inh2, k_223, k_233, geneB_act1, S_dephosphos, geneG_act1, geneH_inh2, geneI_act1, geneB_inh3, geneI_inh2, geneC_turn, geneJ_inh1, geneK_act3, geneJ_inh3, geneG_inh3, geneC_act2, init_TGFb, geneL_inh1, geneE_inh3, geneF_act3, cell, geneF_inh1, geneD_act3, geneE_act2, geneD_inh3, k_234, geneH_act2, geneA_turn, geneL_act1, geneH_act1, geneL_inh2, geneB_turn, init_S4, geneL_inh3, khomo2, geneK_act2, geneA_inh3, geneJ_act2, geneI_inh3, geneD_act1, geneJ_turn, geneG_act3, geneL_act2, Rec_act, geneI_act3, k_224, pRec_degind, geneE_inh2, geneF_inh2, geneC_act1, geneD_inh1, k_344, init_S3, geneI_inh1, geneL_turn, geneF_act2, k_on_u, geneE_inh1, geneH_act3, geneH_inh1, geneK_inh3, geneE_act1, geneA_inh1, geneB_inh1, geneH_inh3, geneD_turn, S_dephos, geneL_act3, S_phos, geneG_turn, geneA_act1, geneI_turn, khomo3, geneC_inh1, geneG_act2, init_S2, geneK_inh2, k_334, geneB_inh2, geneH_turn, geneJ_inh2, khomo4, geneB_act2, geneF_turn, geneJ_act3, geneK_inh1, geneB_act3, geneK_act1, geneE_turn, geneE_act3, geneG_inh2, kdiss_SS, k_244, geneD_act2, geneF_inh3, geneF_act1, geneC_inh2 = p 
	sd_Bmp4_nExpID100, sd_Cxcl15_nExpID100, sd_Dnmt3a_nExpID100, sd_Dusp5_nExpID100, sd_Jun_nExpID100, sd_Klf10_nExpID100, sd_Pdk4_nExpID100, sd_Ski_nExpID100, sd_Skil_nExpID100, sd_Smad7_nExpID100, sd_Sox4_nExpID100, sd_Tgfa_nExpID100 = nonDynParam 
	if observableId == "Bmp4_gene_expression_OE_nExpID100" 
		return nothing
	end

	if observableId == "Cxcl15_gene_expression_OE_nExpID100" 
		return nothing
	end

	if observableId == "Dnmt3a_gene_expression_OE_nExpID100" 
		return nothing
	end

	if observableId == "Dusp5_gene_expression_OE_nExpID100" 
		return nothing
	end

	if observableId == "Jun_gene_expression_OE_nExpID100" 
		return nothing
	end

	if observableId == "Klf10_gene_expression_OE_nExpID100" 
		return nothing
	end

	if observableId == "Pdk4_gene_expression_OE_nExpID100" 
		return nothing
	end

	if observableId == "S23IP_S2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5" 
		return nothing
	end

	if observableId == "S23IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S23IP_S3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5" 
		return nothing
	end

	if observableId == "S23IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S23IP_pS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5" 
		return nothing
	end

	if observableId == "S23IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S23IP_pS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5" 
		return nothing
	end

	if observableId == "S23IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S23IP_ppS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5" 
		return nothing
	end

	if observableId == "S23IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S23IP_ppS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5" 
		return nothing
	end

	if observableId == "S23IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1" 
		return nothing
	end

	if observableId == "S23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10" 
		return nothing
	end

	if observableId == "S23IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S23IP_tS3_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10" 
		return nothing
	end

	if observableId == "S23IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S23IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S2IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S2IP_S2_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return nothing
	end

	if observableId == "S2IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S2IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S2IP_pS2_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return nothing
	end

	if observableId == "S2IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S2IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S2IP_ppS2_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return nothing
	end

	if observableId == "S2IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S2IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S2IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return nothing
	end

	if observableId == "S2IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S2IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return nothing
	end

	if observableId == "S2IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S2IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return nothing
	end

	if observableId == "S3IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S3IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S3IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S3IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S3IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S3IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S3IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S3IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return nothing
	end

	if observableId == "S3IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S3IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return nothing
	end

	if observableId == "S3IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S3IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return nothing
	end

	if observableId == "S4IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S4IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return nothing
	end

	if observableId == "S4IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S4IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return nothing
	end

	if observableId == "S4IP_tS4_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1" 
		return nothing
	end

	if observableId == "S4IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2" 
		return nothing
	end

	if observableId == "S4IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6" 
		return nothing
	end

	if observableId == "Ski_gene_expression_OE_nExpID100" 
		return nothing
	end

	if observableId == "Skil_gene_expression_OE_nExpID100" 
		return nothing
	end

	if observableId == "Smad7_gene_expression_OE_nExpID100" 
		return nothing
	end

	if observableId == "Sox4_gene_expression_OE_nExpID100" 
		return nothing
	end

	if observableId == "TGFbR_PL_PN_TGFbR_hepa16_nExpID11" 
		return nothing
	end

	if observableId == "Tgfa_gene_expression_OE_nExpID100" 
		return nothing
	end

end

