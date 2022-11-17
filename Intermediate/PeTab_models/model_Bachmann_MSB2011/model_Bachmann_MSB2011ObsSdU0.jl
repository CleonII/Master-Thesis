function evalYmod(u, t, dynPar, obsPar, nonDynParam, paramData, observableId, mapObsParam) 

	p1EpoRpJAK2, pSTAT5, EpoRJAK2_CIS, SOCS3nRNA4, SOCS3RNA, SHP1, STAT5, EpoRJAK2, CISnRNA1, SOCS3nRNA1, SOCS3nRNA2, CISnRNA3, CISnRNA4, SOCS3, CISnRNA5, SOCS3nRNA5, SOCS3nRNA3, SHP1Act, npSTAT5, p12EpoRpJAK2, p2EpoRpJAK2, CIS, EpoRpJAK2, CISnRNA2, CISRNA= u 
	CISEqc, CISEqcOE, CISInh, CISRNADelay, CISRNATurn, CISTurn, EpoRActJAK2, EpoRCISInh, EpoRCISRemove, JAK2ActEpo, JAK2EpoRDeaSHP1, SHP1ActEpoR, SHP1Dea, SHP1ProOE, SOCS3Eqc, SOCS3EqcOE, SOCS3Inh, SOCS3RNADelay, SOCS3RNATurn, SOCS3Turn, STAT5ActEpoR, STAT5ActJAK2, STAT5Exp, STAT5Imp, init_EpoRJAK2, init_SHP1, init_STAT5 = dynPar 
	CISRNAEqc_C = paramData.paramVal[5] 
	SOCS3RNAEqc_C = paramData.paramVal[20] 

	if observableId == "observable_CISRNA_foldA" 
		observableParameter1_observable_CISRNA_foldA = getObsOrSdParam(obsPar, mapObsParam)
		return CISRNA * observableParameter1_observable_CISRNA_foldA / CISRNAEqc_C + 1 
	end

	if observableId == "observable_CISRNA_foldB" 
		observableParameter1_observable_CISRNA_foldB = getObsOrSdParam(obsPar, mapObsParam)
		return CISRNA * observableParameter1_observable_CISRNA_foldB / CISRNAEqc_C + 1 
	end

	if observableId == "observable_CISRNA_foldC" 
		observableParameter1_observable_CISRNA_foldC = getObsOrSdParam(obsPar, mapObsParam)
		return CISRNA * observableParameter1_observable_CISRNA_foldC / CISRNAEqc_C + 1 
	end

	if observableId == "observable_CIS_abs" 
		return CIS 
	end

	if observableId == "observable_CIS_au" 
		observableParameter1_observable_CIS_au, observableParameter2_observable_CIS_au = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_observable_CIS_au + CIS * observableParameter2_observable_CIS_au / CISEqc 
	end

	if observableId == "observable_CIS_au1" 
		observableParameter1_observable_CIS_au1 = getObsOrSdParam(obsPar, mapObsParam)
		return CIS * observableParameter1_observable_CIS_au1 / CISEqc 
	end

	if observableId == "observable_CIS_au2" 
		observableParameter1_observable_CIS_au2 = getObsOrSdParam(obsPar, mapObsParam)
		return CIS * observableParameter1_observable_CIS_au2 / CISEqc 
	end

	if observableId == "observable_SHP1_abs" 
		return SHP1 + SHP1Act 
	end

	if observableId == "observable_SOCS3RNA_foldA" 
		observableParameter1_observable_SOCS3RNA_foldA = getObsOrSdParam(obsPar, mapObsParam)
		return SOCS3RNA * observableParameter1_observable_SOCS3RNA_foldA / SOCS3RNAEqc_C + 1 
	end

	if observableId == "observable_SOCS3RNA_foldB" 
		observableParameter1_observable_SOCS3RNA_foldB = getObsOrSdParam(obsPar, mapObsParam)
		return SOCS3RNA * observableParameter1_observable_SOCS3RNA_foldB / SOCS3RNAEqc_C + 1 
	end

	if observableId == "observable_SOCS3RNA_foldC" 
		observableParameter1_observable_SOCS3RNA_foldC = getObsOrSdParam(obsPar, mapObsParam)
		return SOCS3RNA * observableParameter1_observable_SOCS3RNA_foldC / SOCS3RNAEqc_C + 1 
	end

	if observableId == "observable_SOCS3_abs" 
		return SOCS3 
	end

	if observableId == "observable_SOCS3_au" 
		observableParameter1_observable_SOCS3_au, observableParameter2_observable_SOCS3_au = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_observable_SOCS3_au + SOCS3 * observableParameter2_observable_SOCS3_au / SOCS3Eqc 
	end

	if observableId == "observable_STAT5_abs" 
		return STAT5 
	end

	if observableId == "observable_pEpoR_au" 
		observableParameter1_observable_pEpoR_au, observableParameter2_observable_pEpoR_au = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_observable_pEpoR_au + observableParameter2_observable_pEpoR_au * ( 16 * p12EpoRpJAK2 + 16 * p1EpoRpJAK2 + 16 * p2EpoRpJAK2 ) / init_EpoRJAK2 
	end

	if observableId == "observable_pJAK2_au" 
		observableParameter1_observable_pJAK2_au, observableParameter2_observable_pJAK2_au = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_observable_pJAK2_au + observableParameter2_observable_pJAK2_au * ( 2 * EpoRpJAK2 + 2 * p12EpoRpJAK2 + 2 * p1EpoRpJAK2 + 2 * p2EpoRpJAK2 ) / init_EpoRJAK2 
	end

	if observableId == "observable_pSTAT5B_rel" 
		observableParameter1_observable_pSTAT5B_rel = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_observable_pSTAT5B_rel + 100 * pSTAT5 / ( STAT5 + pSTAT5 ) 
	end

	if observableId == "observable_pSTAT5_au" 
		observableParameter1_observable_pSTAT5_au, observableParameter2_observable_pSTAT5_au = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_observable_pSTAT5_au + pSTAT5 * observableParameter2_observable_pSTAT5_au / init_STAT5 
	end

	if observableId == "observable_tSHP1_au" 
		observableParameter1_observable_tSHP1_au = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_observable_tSHP1_au * ( SHP1 + SHP1Act ) / init_SHP1 
	end

	if observableId == "observable_tSTAT5_au" 
		observableParameter1_observable_tSTAT5_au = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_observable_tSTAT5_au * ( STAT5 + pSTAT5 ) / init_STAT5 
	end

end

function evalU0!(u0Vec, paramVec) 

	STAT5Exp, STAT5Imp, init_SOCS3_multiplier, EpoRCISRemove, STAT5ActEpoR, SHP1ActEpoR, JAK2EpoRDeaSHP1, CISTurn, SOCS3Turn, init_EpoRJAK2_CIS, SOCS3Inh, ActD, init_CIS_multiplier, cyt, CISRNAEqc, JAK2ActEpo, Epo, SOCS3oe, CISInh, SHP1Dea, SOCS3EqcOE, CISRNADelay, init_SHP1, CISEqcOE, EpoRActJAK2, SOCS3RNAEqc, CISEqc, SHP1ProOE, SOCS3RNADelay, init_STAT5, CISoe, CISRNATurn, init_SHP1_multiplier, init_EpoRJAK2, nuc, EpoRCISInh, STAT5ActJAK2, SOCS3RNATurn, SOCS3Eqc = paramVec 

	p1EpoRpJAK2 = 0.0 
	pSTAT5 = 0.0 
	EpoRJAK2_CIS = init_EpoRJAK2_CIS 
	SOCS3nRNA4 = 0.0 
	SOCS3RNA = 0.0 
	SHP1 = init_SHP1 * ( 1 + SHP1ProOE * init_SHP1_multiplier ) 
	STAT5 = init_STAT5 
	EpoRJAK2 = init_EpoRJAK2 
	CISnRNA1 = 0.0 
	SOCS3nRNA1 = 0.0 
	SOCS3nRNA2 = 0.0 
	CISnRNA3 = 0.0 
	CISnRNA4 = 0.0 
	SOCS3 = SOCS3Eqc * SOCS3EqcOE * init_SOCS3_multiplier 
	CISnRNA5 = 0.0 
	SOCS3nRNA5 = 0.0 
	SOCS3nRNA3 = 0.0 
	SHP1Act = 0.0 
	npSTAT5 = 0.0 
	p12EpoRpJAK2 = 0.0 
	p2EpoRpJAK2 = 0.0 
	CIS = CISEqc * CISEqcOE * init_CIS_multiplier 
	EpoRpJAK2 = 0.0 
	CISnRNA2 = 0.0 
	CISRNA = 0.0 

	u0Vec .= p1EpoRpJAK2, pSTAT5, EpoRJAK2_CIS, SOCS3nRNA4, SOCS3RNA, SHP1, STAT5, EpoRJAK2, CISnRNA1, SOCS3nRNA1, SOCS3nRNA2, CISnRNA3, CISnRNA4, SOCS3, CISnRNA5, SOCS3nRNA5, SOCS3nRNA3, SHP1Act, npSTAT5, p12EpoRpJAK2, p2EpoRpJAK2, CIS, EpoRpJAK2, CISnRNA2, CISRNA
end

function evalSd!(u, t, sdPar, dynPar, nonDynParam, paramData, observableId, mapSdParam) 

	p1EpoRpJAK2, pSTAT5, EpoRJAK2_CIS, SOCS3nRNA4, SOCS3RNA, SHP1, STAT5, EpoRJAK2, CISnRNA1, SOCS3nRNA1, SOCS3nRNA2, CISnRNA3, CISnRNA4, SOCS3, CISnRNA5, SOCS3nRNA5, SOCS3nRNA3, SHP1Act, npSTAT5, p12EpoRpJAK2, p2EpoRpJAK2, CIS, EpoRpJAK2, CISnRNA2, CISRNA= u 
	CISEqc, CISEqcOE, CISInh, CISRNADelay, CISRNATurn, CISTurn, EpoRActJAK2, EpoRCISInh, EpoRCISRemove, JAK2ActEpo, JAK2EpoRDeaSHP1, SHP1ActEpoR, SHP1Dea, SHP1ProOE, SOCS3Eqc, SOCS3EqcOE, SOCS3Inh, SOCS3RNADelay, SOCS3RNATurn, SOCS3Turn, STAT5ActEpoR, STAT5ActJAK2, STAT5Exp, STAT5Imp, init_EpoRJAK2, init_SHP1, init_STAT5 = dynPar 
	CISRNAEqc_C = paramData.paramVal[5] 
	SOCS3RNAEqc_C = paramData.paramVal[20] 

	if observableId == "observable_CISRNA_foldA" 
		noiseParameter1_observable_CISRNA_foldA = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_CISRNA_foldA 
	end

	if observableId == "observable_CISRNA_foldB" 
		noiseParameter1_observable_CISRNA_foldB = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_CISRNA_foldB 
	end

	if observableId == "observable_CISRNA_foldC" 
		noiseParameter1_observable_CISRNA_foldC = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_CISRNA_foldC 
	end

	if observableId == "observable_CIS_abs" 
		noiseParameter1_observable_CIS_abs = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_CIS_abs 
	end

	if observableId == "observable_CIS_au" 
		noiseParameter1_observable_CIS_au = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_CIS_au 
	end

	if observableId == "observable_CIS_au1" 
		noiseParameter1_observable_CIS_au1 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_CIS_au1 
	end

	if observableId == "observable_CIS_au2" 
		noiseParameter1_observable_CIS_au2 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_CIS_au2 
	end

	if observableId == "observable_SHP1_abs" 
		noiseParameter1_observable_SHP1_abs = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_SHP1_abs 
	end

	if observableId == "observable_SOCS3RNA_foldA" 
		noiseParameter1_observable_SOCS3RNA_foldA = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_SOCS3RNA_foldA 
	end

	if observableId == "observable_SOCS3RNA_foldB" 
		noiseParameter1_observable_SOCS3RNA_foldB = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_SOCS3RNA_foldB 
	end

	if observableId == "observable_SOCS3RNA_foldC" 
		noiseParameter1_observable_SOCS3RNA_foldC = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_SOCS3RNA_foldC 
	end

	if observableId == "observable_SOCS3_abs" 
		noiseParameter1_observable_SOCS3_abs = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_SOCS3_abs 
	end

	if observableId == "observable_SOCS3_au" 
		noiseParameter1_observable_SOCS3_au = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_SOCS3_au 
	end

	if observableId == "observable_STAT5_abs" 
		noiseParameter1_observable_STAT5_abs = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_STAT5_abs 
	end

	if observableId == "observable_pEpoR_au" 
		noiseParameter1_observable_pEpoR_au = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_pEpoR_au 
	end

	if observableId == "observable_pJAK2_au" 
		noiseParameter1_observable_pJAK2_au = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_pJAK2_au 
	end

	if observableId == "observable_pSTAT5B_rel" 
		noiseParameter1_observable_pSTAT5B_rel = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_pSTAT5B_rel 
	end

	if observableId == "observable_pSTAT5_au" 
		noiseParameter1_observable_pSTAT5_au, noiseParameter2_observable_pSTAT5_au = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_pSTAT5_au + noiseParameter2_observable_pSTAT5_au 
	end

	if observableId == "observable_tSHP1_au" 
		noiseParameter1_observable_tSHP1_au = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_tSHP1_au 
	end

	if observableId == "observable_tSTAT5_au" 
		noiseParameter1_observable_tSTAT5_au = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_tSTAT5_au 
	end

end