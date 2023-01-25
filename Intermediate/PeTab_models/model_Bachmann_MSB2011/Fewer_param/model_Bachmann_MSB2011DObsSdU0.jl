function evalDYmodDu(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	p1EpoRpJAK2, pSTAT5, EpoRJAK2_CIS, SOCS3nRNA4, SOCS3RNA, SHP1, STAT5, EpoRJAK2, CISnRNA1, SOCS3nRNA1, SOCS3nRNA2, CISnRNA3, CISnRNA4, SOCS3, CISnRNA5, SOCS3nRNA5, SOCS3nRNA3, SHP1Act, npSTAT5, p12EpoRpJAK2, p2EpoRpJAK2, CIS, EpoRpJAK2, CISnRNA2, CISRNA= u 
	STAT5Exp, STAT5Imp, init_SOCS3_multiplier, EpoRCISRemove, STAT5ActEpoR, SHP1ActEpoR, JAK2EpoRDeaSHP1, CISTurn, SOCS3Turn, init_EpoRJAK2_CIS, SOCS3Inh, ActD, init_CIS_multiplier, cyt, CISRNAEqc, JAK2ActEpo, Epo, SOCS3oe, CISInh, SHP1Dea, SOCS3EqcOE, CISRNADelay, init_SHP1, CISEqcOE, EpoRActJAK2, SOCS3RNAEqc, CISEqc, SHP1ProOE, SOCS3RNADelay, init_STAT5, CISoe, CISRNATurn, init_SHP1_multiplier, init_EpoRJAK2, nuc, EpoRCISInh, STAT5ActJAK2, SOCS3RNATurn, SOCS3Eqc = p 
	if observableId == "observable_CISRNA_foldA" 
		observableParameter1_observable_CISRNA_foldA = getObsOrSdParam(obsPar, mapObsParam)
		out[25] = observableParameter1_observable_CISRNA_foldA / CISRNAEqc
		return nothing
	end

	if observableId == "observable_CISRNA_foldB" 
		observableParameter1_observable_CISRNA_foldB = getObsOrSdParam(obsPar, mapObsParam)
		out[25] = observableParameter1_observable_CISRNA_foldB / CISRNAEqc
		return nothing
	end

	if observableId == "observable_CISRNA_foldC" 
		observableParameter1_observable_CISRNA_foldC = getObsOrSdParam(obsPar, mapObsParam)
		out[25] = observableParameter1_observable_CISRNA_foldC / CISRNAEqc
		return nothing
	end

	if observableId == "observable_CIS_abs" 
		out[22] = 1
		return nothing
	end

	if observableId == "observable_CIS_au" 
		observableParameter1_observable_CIS_au, observableParameter2_observable_CIS_au = getObsOrSdParam(obsPar, mapObsParam)
		out[22] = observableParameter2_observable_CIS_au / CISEqc
		return nothing
	end

	if observableId == "observable_CIS_au1" 
		observableParameter1_observable_CIS_au1 = getObsOrSdParam(obsPar, mapObsParam)
		out[22] = observableParameter1_observable_CIS_au1 / CISEqc
		return nothing
	end

	if observableId == "observable_CIS_au2" 
		observableParameter1_observable_CIS_au2 = getObsOrSdParam(obsPar, mapObsParam)
		out[22] = observableParameter1_observable_CIS_au2 / CISEqc
		return nothing
	end

	if observableId == "observable_SHP1_abs" 
		out[6] = 1
		out[18] = 1
		return nothing
	end

	if observableId == "observable_SOCS3RNA_foldA" 
		observableParameter1_observable_SOCS3RNA_foldA = getObsOrSdParam(obsPar, mapObsParam)
		out[5] = observableParameter1_observable_SOCS3RNA_foldA / SOCS3RNAEqc
		return nothing
	end

	if observableId == "observable_SOCS3RNA_foldB" 
		observableParameter1_observable_SOCS3RNA_foldB = getObsOrSdParam(obsPar, mapObsParam)
		out[5] = observableParameter1_observable_SOCS3RNA_foldB / SOCS3RNAEqc
		return nothing
	end

	if observableId == "observable_SOCS3RNA_foldC" 
		observableParameter1_observable_SOCS3RNA_foldC = getObsOrSdParam(obsPar, mapObsParam)
		out[5] = observableParameter1_observable_SOCS3RNA_foldC / SOCS3RNAEqc
		return nothing
	end

	if observableId == "observable_SOCS3_abs" 
		out[14] = 1
		return nothing
	end

	if observableId == "observable_SOCS3_au" 
		observableParameter1_observable_SOCS3_au, observableParameter2_observable_SOCS3_au = getObsOrSdParam(obsPar, mapObsParam)
		out[14] = observableParameter2_observable_SOCS3_au / SOCS3Eqc
		return nothing
	end

	if observableId == "observable_STAT5_abs" 
		out[7] = 1
		return nothing
	end

	if observableId == "observable_pEpoR_au" 
		observableParameter1_observable_pEpoR_au, observableParameter2_observable_pEpoR_au = getObsOrSdParam(obsPar, mapObsParam)
		out[1] = (16observableParameter2_observable_pEpoR_au) / init_EpoRJAK2
		out[20] = (16observableParameter2_observable_pEpoR_au) / init_EpoRJAK2
		out[21] = (16observableParameter2_observable_pEpoR_au) / init_EpoRJAK2
		return nothing
	end

	if observableId == "observable_pJAK2_au" 
		observableParameter1_observable_pJAK2_au, observableParameter2_observable_pJAK2_au = getObsOrSdParam(obsPar, mapObsParam)
		out[1] = (2observableParameter2_observable_pJAK2_au) / init_EpoRJAK2
		out[20] = (2observableParameter2_observable_pJAK2_au) / init_EpoRJAK2
		out[21] = (2observableParameter2_observable_pJAK2_au) / init_EpoRJAK2
		out[23] = (2observableParameter2_observable_pJAK2_au) / init_EpoRJAK2
		return nothing
	end

	if observableId == "observable_pSTAT5B_rel" 
		observableParameter1_observable_pSTAT5B_rel = getObsOrSdParam(obsPar, mapObsParam)
		out[2] = (100.0STAT5) / ((STAT5 + pSTAT5)^2)
		out[7] = (-100pSTAT5) / ((STAT5 + pSTAT5)^2)
		return nothing
	end

	if observableId == "observable_pSTAT5_au" 
		observableParameter1_observable_pSTAT5_au, observableParameter2_observable_pSTAT5_au = getObsOrSdParam(obsPar, mapObsParam)
		out[2] = observableParameter2_observable_pSTAT5_au / init_STAT5
		return nothing
	end

	if observableId == "observable_tSHP1_au" 
		observableParameter1_observable_tSHP1_au = getObsOrSdParam(obsPar, mapObsParam)
		out[6] = observableParameter1_observable_tSHP1_au / init_SHP1
		out[18] = observableParameter1_observable_tSHP1_au / init_SHP1
		return nothing
	end

	if observableId == "observable_tSTAT5_au" 
		observableParameter1_observable_tSTAT5_au = getObsOrSdParam(obsPar, mapObsParam)
		out[2] = observableParameter1_observable_tSTAT5_au / init_STAT5
		out[7] = observableParameter1_observable_tSTAT5_au / init_STAT5
		return nothing
	end

end

function evalDYmodDp(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	p1EpoRpJAK2, pSTAT5, EpoRJAK2_CIS, SOCS3nRNA4, SOCS3RNA, SHP1, STAT5, EpoRJAK2, CISnRNA1, SOCS3nRNA1, SOCS3nRNA2, CISnRNA3, CISnRNA4, SOCS3, CISnRNA5, SOCS3nRNA5, SOCS3nRNA3, SHP1Act, npSTAT5, p12EpoRpJAK2, p2EpoRpJAK2, CIS, EpoRpJAK2, CISnRNA2, CISRNA= u 
	STAT5Exp, STAT5Imp, init_SOCS3_multiplier, EpoRCISRemove, STAT5ActEpoR, SHP1ActEpoR, JAK2EpoRDeaSHP1, CISTurn, SOCS3Turn, init_EpoRJAK2_CIS, SOCS3Inh, ActD, init_CIS_multiplier, cyt, CISRNAEqc, JAK2ActEpo, Epo, SOCS3oe, CISInh, SHP1Dea, SOCS3EqcOE, CISRNADelay, init_SHP1, CISEqcOE, EpoRActJAK2, SOCS3RNAEqc, CISEqc, SHP1ProOE, SOCS3RNADelay, init_STAT5, CISoe, CISRNATurn, init_SHP1_multiplier, init_EpoRJAK2, nuc, EpoRCISInh, STAT5ActJAK2, SOCS3RNATurn, SOCS3Eqc = p 
	if observableId == "observable_CISRNA_foldA" 
		observableParameter1_observable_CISRNA_foldA = getObsOrSdParam(obsPar, mapObsParam)
		out[15] = (-CISRNA*observableParameter1_observable_CISRNA_foldA) / (CISRNAEqc^2)
		return nothing
	end

	if observableId == "observable_CISRNA_foldB" 
		observableParameter1_observable_CISRNA_foldB = getObsOrSdParam(obsPar, mapObsParam)
		out[15] = (-CISRNA*observableParameter1_observable_CISRNA_foldB) / (CISRNAEqc^2)
		return nothing
	end

	if observableId == "observable_CISRNA_foldC" 
		observableParameter1_observable_CISRNA_foldC = getObsOrSdParam(obsPar, mapObsParam)
		out[15] = (-CISRNA*observableParameter1_observable_CISRNA_foldC) / (CISRNAEqc^2)
		return nothing
	end

	if observableId == "observable_CIS_abs" 
		return nothing
	end

	if observableId == "observable_CIS_au" 
		observableParameter1_observable_CIS_au, observableParameter2_observable_CIS_au = getObsOrSdParam(obsPar, mapObsParam)
		out[27] = (-CIS*observableParameter2_observable_CIS_au) / (CISEqc^2)
		return nothing
	end

	if observableId == "observable_CIS_au1" 
		observableParameter1_observable_CIS_au1 = getObsOrSdParam(obsPar, mapObsParam)
		out[27] = (-CIS*observableParameter1_observable_CIS_au1) / (CISEqc^2)
		return nothing
	end

	if observableId == "observable_CIS_au2" 
		observableParameter1_observable_CIS_au2 = getObsOrSdParam(obsPar, mapObsParam)
		out[27] = (-CIS*observableParameter1_observable_CIS_au2) / (CISEqc^2)
		return nothing
	end

	if observableId == "observable_SHP1_abs" 
		return nothing
	end

	if observableId == "observable_SOCS3RNA_foldA" 
		observableParameter1_observable_SOCS3RNA_foldA = getObsOrSdParam(obsPar, mapObsParam)
		out[26] = (-SOCS3RNA*observableParameter1_observable_SOCS3RNA_foldA) / (SOCS3RNAEqc^2)
		return nothing
	end

	if observableId == "observable_SOCS3RNA_foldB" 
		observableParameter1_observable_SOCS3RNA_foldB = getObsOrSdParam(obsPar, mapObsParam)
		out[26] = (-SOCS3RNA*observableParameter1_observable_SOCS3RNA_foldB) / (SOCS3RNAEqc^2)
		return nothing
	end

	if observableId == "observable_SOCS3RNA_foldC" 
		observableParameter1_observable_SOCS3RNA_foldC = getObsOrSdParam(obsPar, mapObsParam)
		out[26] = (-SOCS3RNA*observableParameter1_observable_SOCS3RNA_foldC) / (SOCS3RNAEqc^2)
		return nothing
	end

	if observableId == "observable_SOCS3_abs" 
		return nothing
	end

	if observableId == "observable_SOCS3_au" 
		observableParameter1_observable_SOCS3_au, observableParameter2_observable_SOCS3_au = getObsOrSdParam(obsPar, mapObsParam)
		out[39] = (-SOCS3*observableParameter2_observable_SOCS3_au) / (SOCS3Eqc^2)
		return nothing
	end

	if observableId == "observable_STAT5_abs" 
		return nothing
	end

	if observableId == "observable_pEpoR_au" 
		observableParameter1_observable_pEpoR_au, observableParameter2_observable_pEpoR_au = getObsOrSdParam(obsPar, mapObsParam)
		out[34] = (observableParameter2_observable_pEpoR_au*(-16p12EpoRpJAK2 - 16p1EpoRpJAK2 - 16p2EpoRpJAK2)) / (init_EpoRJAK2^2)
		return nothing
	end

	if observableId == "observable_pJAK2_au" 
		observableParameter1_observable_pJAK2_au, observableParameter2_observable_pJAK2_au = getObsOrSdParam(obsPar, mapObsParam)
		out[34] = (-observableParameter2_observable_pJAK2_au*(2EpoRpJAK2 + 2p12EpoRpJAK2 + 2p1EpoRpJAK2 + 2p2EpoRpJAK2)) / (init_EpoRJAK2^2)
		return nothing
	end

	if observableId == "observable_pSTAT5B_rel" 
		return nothing
	end

	if observableId == "observable_pSTAT5_au" 
		observableParameter1_observable_pSTAT5_au, observableParameter2_observable_pSTAT5_au = getObsOrSdParam(obsPar, mapObsParam)
		out[30] = (-observableParameter2_observable_pSTAT5_au*pSTAT5) / (init_STAT5^2)
		return nothing
	end

	if observableId == "observable_tSHP1_au" 
		observableParameter1_observable_tSHP1_au = getObsOrSdParam(obsPar, mapObsParam)
		out[23] = (-observableParameter1_observable_tSHP1_au*(SHP1 + SHP1Act)) / (init_SHP1^2)
		return nothing
	end

	if observableId == "observable_tSTAT5_au" 
		observableParameter1_observable_tSTAT5_au = getObsOrSdParam(obsPar, mapObsParam)
		out[30] = (observableParameter1_observable_tSTAT5_au*(-STAT5 - pSTAT5)) / (init_STAT5^2)
		return nothing
	end

end

function evalDSdDu!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	p1EpoRpJAK2, pSTAT5, EpoRJAK2_CIS, SOCS3nRNA4, SOCS3RNA, SHP1, STAT5, EpoRJAK2, CISnRNA1, SOCS3nRNA1, SOCS3nRNA2, CISnRNA3, CISnRNA4, SOCS3, CISnRNA5, SOCS3nRNA5, SOCS3nRNA3, SHP1Act, npSTAT5, p12EpoRpJAK2, p2EpoRpJAK2, CIS, EpoRpJAK2, CISnRNA2, CISRNA= u 
	STAT5Exp, STAT5Imp, init_SOCS3_multiplier, EpoRCISRemove, STAT5ActEpoR, SHP1ActEpoR, JAK2EpoRDeaSHP1, CISTurn, SOCS3Turn, init_EpoRJAK2_CIS, SOCS3Inh, ActD, init_CIS_multiplier, cyt, CISRNAEqc, JAK2ActEpo, Epo, SOCS3oe, CISInh, SHP1Dea, SOCS3EqcOE, CISRNADelay, init_SHP1, CISEqcOE, EpoRActJAK2, SOCS3RNAEqc, CISEqc, SHP1ProOE, SOCS3RNADelay, init_STAT5, CISoe, CISRNATurn, init_SHP1_multiplier, init_EpoRJAK2, nuc, EpoRCISInh, STAT5ActJAK2, SOCS3RNATurn, SOCS3Eqc = p 
	if observableId == "observable_CISRNA_foldA" 
		return nothing
	end

	if observableId == "observable_CISRNA_foldB" 
		return nothing
	end

	if observableId == "observable_CISRNA_foldC" 
		return nothing
	end

	if observableId == "observable_CIS_abs" 
		return nothing
	end

	if observableId == "observable_CIS_au" 
		return nothing
	end

	if observableId == "observable_CIS_au1" 
		return nothing
	end

	if observableId == "observable_CIS_au2" 
		return nothing
	end

	if observableId == "observable_SHP1_abs" 
		return nothing
	end

	if observableId == "observable_SOCS3RNA_foldA" 
		return nothing
	end

	if observableId == "observable_SOCS3RNA_foldB" 
		return nothing
	end

	if observableId == "observable_SOCS3RNA_foldC" 
		return nothing
	end

	if observableId == "observable_SOCS3_abs" 
		return nothing
	end

	if observableId == "observable_SOCS3_au" 
		return nothing
	end

	if observableId == "observable_STAT5_abs" 
		return nothing
	end

	if observableId == "observable_pEpoR_au" 
		return nothing
	end

	if observableId == "observable_pJAK2_au" 
		return nothing
	end

	if observableId == "observable_pSTAT5B_rel" 
		return nothing
	end

	if observableId == "observable_pSTAT5_au" 
		return nothing
	end

	if observableId == "observable_tSHP1_au" 
		return nothing
	end

	if observableId == "observable_tSTAT5_au" 
		return nothing
	end

end

function evalDSdDp!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	p1EpoRpJAK2, pSTAT5, EpoRJAK2_CIS, SOCS3nRNA4, SOCS3RNA, SHP1, STAT5, EpoRJAK2, CISnRNA1, SOCS3nRNA1, SOCS3nRNA2, CISnRNA3, CISnRNA4, SOCS3, CISnRNA5, SOCS3nRNA5, SOCS3nRNA3, SHP1Act, npSTAT5, p12EpoRpJAK2, p2EpoRpJAK2, CIS, EpoRpJAK2, CISnRNA2, CISRNA= u 
	STAT5Exp, STAT5Imp, init_SOCS3_multiplier, EpoRCISRemove, STAT5ActEpoR, SHP1ActEpoR, JAK2EpoRDeaSHP1, CISTurn, SOCS3Turn, init_EpoRJAK2_CIS, SOCS3Inh, ActD, init_CIS_multiplier, cyt, CISRNAEqc, JAK2ActEpo, Epo, SOCS3oe, CISInh, SHP1Dea, SOCS3EqcOE, CISRNADelay, init_SHP1, CISEqcOE, EpoRActJAK2, SOCS3RNAEqc, CISEqc, SHP1ProOE, SOCS3RNADelay, init_STAT5, CISoe, CISRNATurn, init_SHP1_multiplier, init_EpoRJAK2, nuc, EpoRCISInh, STAT5ActJAK2, SOCS3RNATurn, SOCS3Eqc = p 
	if observableId == "observable_CISRNA_foldA" 
		return nothing
	end

	if observableId == "observable_CISRNA_foldB" 
		return nothing
	end

	if observableId == "observable_CISRNA_foldC" 
		return nothing
	end

	if observableId == "observable_CIS_abs" 
		return nothing
	end

	if observableId == "observable_CIS_au" 
		return nothing
	end

	if observableId == "observable_CIS_au1" 
		return nothing
	end

	if observableId == "observable_CIS_au2" 
		return nothing
	end

	if observableId == "observable_SHP1_abs" 
		return nothing
	end

	if observableId == "observable_SOCS3RNA_foldA" 
		return nothing
	end

	if observableId == "observable_SOCS3RNA_foldB" 
		return nothing
	end

	if observableId == "observable_SOCS3RNA_foldC" 
		return nothing
	end

	if observableId == "observable_SOCS3_abs" 
		return nothing
	end

	if observableId == "observable_SOCS3_au" 
		return nothing
	end

	if observableId == "observable_STAT5_abs" 
		return nothing
	end

	if observableId == "observable_pEpoR_au" 
		return nothing
	end

	if observableId == "observable_pJAK2_au" 
		return nothing
	end

	if observableId == "observable_pSTAT5B_rel" 
		return nothing
	end

	if observableId == "observable_pSTAT5_au" 
		return nothing
	end

	if observableId == "observable_tSHP1_au" 
		return nothing
	end

	if observableId == "observable_tSTAT5_au" 
		return nothing
	end

end

