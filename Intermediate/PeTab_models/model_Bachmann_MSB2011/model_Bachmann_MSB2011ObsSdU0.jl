function evalYmod(u, t, dynPar, obsPar, paramData, obsData, observableId, simulationId) 

	p1EpoRpJAK2, pSTAT5, EpoRJAK2_CIS, SOCS3nRNA4, SOCS3RNA, SHP1, STAT5, EpoRJAK2, CISnRNA1, SOCS3nRNA1, SOCS3nRNA2, CISnRNA3, CISnRNA4, SOCS3, CISnRNA5, SOCS3nRNA5, SOCS3nRNA3, SHP1Act, npSTAT5, p12EpoRpJAK2, p2EpoRpJAK2, CIS, EpoRpJAK2, CISnRNA2, CISRNA, dummyVariable= u 
	CISEqc, CISEqcOE, CISInh, CISRNADelay, CISRNATurn, CISTurn, EpoRActJAK2, EpoRCISInh, EpoRCISRemove, JAK2ActEpo, JAK2EpoRDeaSHP1, SHP1ActEpoR, SHP1Dea, SHP1ProOE, SOCS3Eqc, SOCS3EqcOE, SOCS3Inh, SOCS3RNADelay, SOCS3RNATurn, SOCS3Turn, STAT5ActEpoR, STAT5ActJAK2, STAT5Exp, STAT5Imp, init_EpoRJAK2, init_SHP1, init_STAT5 = dynPar 
	CISRNAEqc_C = paramData.paramVal[5] 
	SOCS3RNAEqc_C = paramData.paramVal[20] 

	if observableId == "observable_CISRNA_foldA" 
		observableParameter1_observable_CISRNA_foldA = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return CISRNA * observableParameter1_observable_CISRNA_foldA / CISRNAEqc_C + 1 
	end

	if observableId == "observable_CISRNA_foldB" 
		observableParameter1_observable_CISRNA_foldB = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return CISRNA * observableParameter1_observable_CISRNA_foldB / CISRNAEqc_C + 1 
	end

	if observableId == "observable_CISRNA_foldC" 
		observableParameter1_observable_CISRNA_foldC = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return CISRNA * observableParameter1_observable_CISRNA_foldC / CISRNAEqc_C + 1 
	end

	if observableId == "observable_CIS_abs" 
		return CIS 
	end

	if observableId == "observable_CIS_au" 
		observableParameter1_observable_CIS_au, observableParameter2_observable_CIS_au = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return observableParameter1_observable_CIS_au + CIS * observableParameter2_observable_CIS_au / CISEqc 
	end

	if observableId == "observable_CIS_au1" 
		observableParameter1_observable_CIS_au1 = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return CIS * observableParameter1_observable_CIS_au1 / CISEqc 
	end

	if observableId == "observable_CIS_au2" 
		observableParameter1_observable_CIS_au2 = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return CIS * observableParameter1_observable_CIS_au2 / CISEqc 
	end

	if observableId == "observable_SHP1_abs" 
		return SHP1 + SHP1Act 
	end

	if observableId == "observable_SOCS3RNA_foldA" 
		observableParameter1_observable_SOCS3RNA_foldA = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return SOCS3RNA * observableParameter1_observable_SOCS3RNA_foldA / SOCS3RNAEqc_C + 1 
	end

	if observableId == "observable_SOCS3RNA_foldB" 
		observableParameter1_observable_SOCS3RNA_foldB = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return SOCS3RNA * observableParameter1_observable_SOCS3RNA_foldB / SOCS3RNAEqc_C + 1 
	end

	if observableId == "observable_SOCS3RNA_foldC" 
		observableParameter1_observable_SOCS3RNA_foldC = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return SOCS3RNA * observableParameter1_observable_SOCS3RNA_foldC / SOCS3RNAEqc_C + 1 
	end

	if observableId == "observable_SOCS3_abs" 
		return SOCS3 
	end

	if observableId == "observable_SOCS3_au" 
		observableParameter1_observable_SOCS3_au, observableParameter2_observable_SOCS3_au = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return observableParameter1_observable_SOCS3_au + SOCS3 * observableParameter2_observable_SOCS3_au / SOCS3Eqc 
	end

	if observableId == "observable_STAT5_abs" 
		return STAT5 
	end

	if observableId == "observable_pEpoR_au" 
		observableParameter1_observable_pEpoR_au, observableParameter2_observable_pEpoR_au = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return observableParameter1_observable_pEpoR_au + observableParameter2_observable_pEpoR_au * ( 16 * p12EpoRpJAK2 + 16 * p1EpoRpJAK2 + 16 * p2EpoRpJAK2 ) / init_EpoRJAK2 
	end

	if observableId == "observable_pJAK2_au" 
		observableParameter1_observable_pJAK2_au, observableParameter2_observable_pJAK2_au = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return observableParameter1_observable_pJAK2_au + observableParameter2_observable_pJAK2_au * ( 2 * EpoRpJAK2 + 2 * p12EpoRpJAK2 + 2 * p1EpoRpJAK2 + 2 * p2EpoRpJAK2 ) / init_EpoRJAK2 
	end

	if observableId == "observable_pSTAT5B_rel" 
		observableParameter1_observable_pSTAT5B_rel = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return observableParameter1_observable_pSTAT5B_rel + 100 * pSTAT5 / ( STAT5 + pSTAT5 ) 
	end

	if observableId == "observable_pSTAT5_au" 
		observableParameter1_observable_pSTAT5_au, observableParameter2_observable_pSTAT5_au = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return observableParameter1_observable_pSTAT5_au + pSTAT5 * observableParameter2_observable_pSTAT5_au / init_STAT5 
	end

	if observableId == "observable_tSHP1_au" 
		observableParameter1_observable_tSHP1_au = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return observableParameter1_observable_tSHP1_au * ( SHP1 + SHP1Act ) / init_SHP1 
	end

	if observableId == "observable_tSTAT5_au" 
		observableParameter1_observable_tSTAT5_au = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return observableParameter1_observable_tSTAT5_au * ( STAT5 + pSTAT5 ) / init_STAT5 
	end

end

function evalU0!(u0Vec, paramVec) 

	init_SHP1, EpoRCISInh, EpoRActJAK2, SOCS3Eqc, SOCS3Inh, JAK2EpoRDeaSHP1, init_EpoRJAK2, STAT5Imp, STAT5ActJAK2, CISEqc, STAT5ActEpoR, CISInh, EpoRCISRemove, SOCS3RNADelay, nuc, cyt, SOCS3RNATurn, SHP1Dea, SHP1ActEpoR, STAT5Exp, JAK2ActEpo, Epo, CISRNAEqc, CISRNATurn, CISRNADelay, init_STAT5, ActD, SOCS3RNAEqc, SOCS3Turn, SOCS3EqcOE, SOCS3oe, CISEqcOE, CISTurn, CISoe, init_SOCS3_multiplier, SHP1ProOE, init_EpoRJAK2_CIS, init_SHP1_multiplier, init_CIS_multiplier = paramVec 

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
	dummyVariable = 0.0 

	u0Vec .= p1EpoRpJAK2, pSTAT5, EpoRJAK2_CIS, SOCS3nRNA4, SOCS3RNA, SHP1, STAT5, EpoRJAK2, CISnRNA1, SOCS3nRNA1, SOCS3nRNA2, CISnRNA3, CISnRNA4, SOCS3, CISnRNA5, SOCS3nRNA5, SOCS3nRNA3, SHP1Act, npSTAT5, p12EpoRpJAK2, p2EpoRpJAK2, CIS, EpoRpJAK2, CISnRNA2, CISRNA, dummyVariable
end

function evalSd!(u, t, sdPar, dynPar, paramData, obsData, observableId, simulationId) 

	p1EpoRpJAK2, pSTAT5, EpoRJAK2_CIS, SOCS3nRNA4, SOCS3RNA, SHP1, STAT5, EpoRJAK2, CISnRNA1, SOCS3nRNA1, SOCS3nRNA2, CISnRNA3, CISnRNA4, SOCS3, CISnRNA5, SOCS3nRNA5, SOCS3nRNA3, SHP1Act, npSTAT5, p12EpoRpJAK2, p2EpoRpJAK2, CIS, EpoRpJAK2, CISnRNA2, CISRNA, dummyVariable= u 
	CISEqc, CISEqcOE, CISInh, CISRNADelay, CISRNATurn, CISTurn, EpoRActJAK2, EpoRCISInh, EpoRCISRemove, JAK2ActEpo, JAK2EpoRDeaSHP1, SHP1ActEpoR, SHP1Dea, SHP1ProOE, SOCS3Eqc, SOCS3EqcOE, SOCS3Inh, SOCS3RNADelay, SOCS3RNATurn, SOCS3Turn, STAT5ActEpoR, STAT5ActJAK2, STAT5Exp, STAT5Imp, init_EpoRJAK2, init_SHP1, init_STAT5 = dynPar 
	CISRNAEqc_C = paramData.paramVal[5] 
	SOCS3RNAEqc_C = paramData.paramVal[20] 

	if observableId == "observable_CISRNA_foldA" 
		noiseParameter1_observable_CISRNA_foldA = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_CISRNA_foldA 
	end

	if observableId == "observable_CISRNA_foldB" 
		noiseParameter1_observable_CISRNA_foldB = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_CISRNA_foldB 
	end

	if observableId == "observable_CISRNA_foldC" 
		noiseParameter1_observable_CISRNA_foldC = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_CISRNA_foldC 
	end

	if observableId == "observable_CIS_abs" 
		noiseParameter1_observable_CIS_abs = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_CIS_abs 
	end

	if observableId == "observable_CIS_au" 
		noiseParameter1_observable_CIS_au = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_CIS_au 
	end

	if observableId == "observable_CIS_au1" 
		noiseParameter1_observable_CIS_au1 = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_CIS_au1 
	end

	if observableId == "observable_CIS_au2" 
		noiseParameter1_observable_CIS_au2 = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_CIS_au2 
	end

	if observableId == "observable_SHP1_abs" 
		noiseParameter1_observable_SHP1_abs = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_SHP1_abs 
	end

	if observableId == "observable_SOCS3RNA_foldA" 
		noiseParameter1_observable_SOCS3RNA_foldA = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_SOCS3RNA_foldA 
	end

	if observableId == "observable_SOCS3RNA_foldB" 
		noiseParameter1_observable_SOCS3RNA_foldB = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_SOCS3RNA_foldB 
	end

	if observableId == "observable_SOCS3RNA_foldC" 
		noiseParameter1_observable_SOCS3RNA_foldC = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_SOCS3RNA_foldC 
	end

	if observableId == "observable_SOCS3_abs" 
		noiseParameter1_observable_SOCS3_abs = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_SOCS3_abs 
	end

	if observableId == "observable_SOCS3_au" 
		noiseParameter1_observable_SOCS3_au = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_SOCS3_au 
	end

	if observableId == "observable_STAT5_abs" 
		noiseParameter1_observable_STAT5_abs = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_STAT5_abs 
	end

	if observableId == "observable_pEpoR_au" 
		noiseParameter1_observable_pEpoR_au = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_pEpoR_au 
	end

	if observableId == "observable_pJAK2_au" 
		noiseParameter1_observable_pJAK2_au = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_pJAK2_au 
	end

	if observableId == "observable_pSTAT5B_rel" 
		noiseParameter1_observable_pSTAT5B_rel = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_pSTAT5B_rel 
	end

	if observableId == "observable_pSTAT5_au" 
		noiseParameter1_observable_pSTAT5_au, noiseParameter2_observable_pSTAT5_au = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_pSTAT5_au + noiseParameter2_observable_pSTAT5_au 
	end

	if observableId == "observable_tSHP1_au" 
		noiseParameter1_observable_tSHP1_au = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_tSHP1_au 
	end

	if observableId == "observable_tSTAT5_au" 
		noiseParameter1_observable_tSTAT5_au = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_tSTAT5_au 
	end

end