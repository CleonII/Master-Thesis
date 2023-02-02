#u[1] = p1EpoRpJAK2, u[2] = pSTAT5, u[3] = EpoRJAK2_CIS, u[4] = SOCS3nRNA4, u[5] = SOCS3RNA, u[6] = SHP1, u[7] = STAT5, u[8] = EpoRJAK2, u[9] = CISnRNA1, u[10] = SOCS3nRNA1, u[11] = SOCS3nRNA2, u[12] = CISnRNA3, u[13] = CISnRNA4, u[14] = SOCS3, u[15] = CISnRNA5, u[16] = SOCS3nRNA5, u[17] = SOCS3nRNA3, u[18] = SHP1Act, u[19] = npSTAT5, u[20] = p12EpoRpJAK2, u[21] = p2EpoRpJAK2, u[22] = CIS, u[23] = EpoRpJAK2, u[24] = CISnRNA2, u[25] = CISRNA
#θ_dynamicNames[1] = CISEqc, θ_dynamicNames[2] = CISEqcOE, θ_dynamicNames[3] = CISInh, θ_dynamicNames[4] = CISRNADelay, θ_dynamicNames[5] = CISRNATurn, θ_dynamicNames[6] = CISTurn, θ_dynamicNames[7] = EpoRActJAK2, θ_dynamicNames[8] = EpoRCISInh, θ_dynamicNames[9] = EpoRCISRemove, θ_dynamicNames[10] = JAK2ActEpo, θ_dynamicNames[11] = JAK2EpoRDeaSHP1, θ_dynamicNames[12] = SHP1ActEpoR, θ_dynamicNames[13] = SHP1Dea, θ_dynamicNames[14] = SHP1ProOE, θ_dynamicNames[15] = SOCS3Eqc, θ_dynamicNames[16] = SOCS3EqcOE, θ_dynamicNames[17] = SOCS3Inh, θ_dynamicNames[18] = SOCS3RNADelay, θ_dynamicNames[19] = SOCS3RNATurn, θ_dynamicNames[20] = SOCS3Turn, θ_dynamicNames[21] = STAT5ActEpoR, θ_dynamicNames[22] = STAT5ActJAK2, θ_dynamicNames[23] = STAT5Exp, θ_dynamicNames[24] = STAT5Imp, θ_dynamicNames[25] = init_EpoRJAK2, θ_dynamicNames[26] = init_SHP1, θ_dynamicNames[27] = init_STAT5
##parameterInfo.nominalValue[5] = CISRNAEqc_C 
#parameterInfo.nominalValue[20] = SOCS3RNAEqc_C 


function compute_h(u::AbstractVector, t::Real, θ_dynamic::AbstractVector, θ_observable::AbstractVector,
                   θ_nonDynamic::AbstractVector, parameterInfo::ParametersInfo, observableId::Symbol, 
                      parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :observable_CISRNA_foldA 
		observableParameter1_observable_CISRNA_foldA = getObsOrSdParam(θ_observable, parameterMap)
		return u[25] * observableParameter1_observable_CISRNA_foldA / parameterInfo.nominalValue[5] + 1 
	end

	if observableId == :observable_CISRNA_foldB 
		observableParameter1_observable_CISRNA_foldB = getObsOrSdParam(θ_observable, parameterMap)
		return u[25] * observableParameter1_observable_CISRNA_foldB / parameterInfo.nominalValue[5] + 1 
	end

	if observableId == :observable_CISRNA_foldC 
		observableParameter1_observable_CISRNA_foldC = getObsOrSdParam(θ_observable, parameterMap)
		return u[25] * observableParameter1_observable_CISRNA_foldC / parameterInfo.nominalValue[5] + 1 
	end

	if observableId == :observable_CIS_abs 
		return u[22] 
	end

	if observableId == :observable_CIS_au 
		observableParameter1_observable_CIS_au, observableParameter2_observable_CIS_au = getObsOrSdParam(θ_observable, parameterMap)
		return observableParameter1_observable_CIS_au + u[22] * observableParameter2_observable_CIS_au / θ_dynamic[1] 
	end

	if observableId == :observable_CIS_au1 
		observableParameter1_observable_CIS_au1 = getObsOrSdParam(θ_observable, parameterMap)
		return u[22] * observableParameter1_observable_CIS_au1 / θ_dynamic[1] 
	end

	if observableId == :observable_CIS_au2 
		observableParameter1_observable_CIS_au2 = getObsOrSdParam(θ_observable, parameterMap)
		return u[22] * observableParameter1_observable_CIS_au2 / θ_dynamic[1] 
	end

	if observableId == :observable_SHP1_abs 
		return u[6] + u[18] 
	end

	if observableId == :observable_SOCS3RNA_foldA 
		observableParameter1_observable_SOCS3RNA_foldA = getObsOrSdParam(θ_observable, parameterMap)
		return u[5] * observableParameter1_observable_SOCS3RNA_foldA / parameterInfo.nominalValue[20] + 1 
	end

	if observableId == :observable_SOCS3RNA_foldB 
		observableParameter1_observable_SOCS3RNA_foldB = getObsOrSdParam(θ_observable, parameterMap)
		return u[5] * observableParameter1_observable_SOCS3RNA_foldB / parameterInfo.nominalValue[20] + 1 
	end

	if observableId == :observable_SOCS3RNA_foldC 
		observableParameter1_observable_SOCS3RNA_foldC = getObsOrSdParam(θ_observable, parameterMap)
		return u[5] * observableParameter1_observable_SOCS3RNA_foldC / parameterInfo.nominalValue[20] + 1 
	end

	if observableId == :observable_SOCS3_abs 
		return u[14] 
	end

	if observableId == :observable_SOCS3_au 
		observableParameter1_observable_SOCS3_au, observableParameter2_observable_SOCS3_au = getObsOrSdParam(θ_observable, parameterMap)
		return observableParameter1_observable_SOCS3_au + u[14] * observableParameter2_observable_SOCS3_au / θ_dynamic[15] 
	end

	if observableId == :observable_STAT5_abs 
		return u[7] 
	end

	if observableId == :observable_pEpoR_au 
		observableParameter1_observable_pEpoR_au, observableParameter2_observable_pEpoR_au = getObsOrSdParam(θ_observable, parameterMap)
		return observableParameter1_observable_pEpoR_au + observableParameter2_observable_pEpoR_au * ( 16 * u[20] + 16 * u[1] + 16 * u[21] ) / θ_dynamic[25] 
	end

	if observableId == :observable_pJAK2_au 
		observableParameter1_observable_pJAK2_au, observableParameter2_observable_pJAK2_au = getObsOrSdParam(θ_observable, parameterMap)
		return observableParameter1_observable_pJAK2_au + observableParameter2_observable_pJAK2_au * ( 2 * u[23] + 2 * u[20] + 2 * u[1] + 2 * u[21] ) / θ_dynamic[25] 
	end

	if observableId == :observable_pSTAT5B_rel 
		observableParameter1_observable_pSTAT5B_rel = getObsOrSdParam(θ_observable, parameterMap)
		return observableParameter1_observable_pSTAT5B_rel + 100 * u[2] / ( u[7] + u[2] ) 
	end

	if observableId == :observable_pSTAT5_au 
		observableParameter1_observable_pSTAT5_au, observableParameter2_observable_pSTAT5_au = getObsOrSdParam(θ_observable, parameterMap)
		return observableParameter1_observable_pSTAT5_au + u[2] * observableParameter2_observable_pSTAT5_au / θ_dynamic[27] 
	end

	if observableId == :observable_tSHP1_au 
		observableParameter1_observable_tSHP1_au = getObsOrSdParam(θ_observable, parameterMap)
		return observableParameter1_observable_tSHP1_au * ( u[6] + u[18] ) / θ_dynamic[26] 
	end

	if observableId == :observable_tSTAT5_au 
		observableParameter1_observable_tSTAT5_au = getObsOrSdParam(θ_observable, parameterMap)
		return observableParameter1_observable_tSTAT5_au * ( u[7] + u[2] ) / θ_dynamic[27] 
	end

end

function compute_u0!(u0::AbstractVector, pODEProblem::AbstractVector) 

	#pODEProblem[1] = STAT5Exp, pODEProblem[2] = STAT5Imp, pODEProblem[3] = init_SOCS3_multiplier, pODEProblem[4] = EpoRCISRemove, pODEProblem[5] = STAT5ActEpoR, pODEProblem[6] = SHP1ActEpoR, pODEProblem[7] = JAK2EpoRDeaSHP1, pODEProblem[8] = CISTurn, pODEProblem[9] = SOCS3Turn, pODEProblem[10] = init_EpoRJAK2_CIS, pODEProblem[11] = SOCS3Inh, pODEProblem[12] = ActD, pODEProblem[13] = init_CIS_multiplier, pODEProblem[14] = cyt, pODEProblem[15] = CISRNAEqc, pODEProblem[16] = JAK2ActEpo, pODEProblem[17] = Epo, pODEProblem[18] = SOCS3oe, pODEProblem[19] = CISInh, pODEProblem[20] = SHP1Dea, pODEProblem[21] = SOCS3EqcOE, pODEProblem[22] = CISRNADelay, pODEProblem[23] = init_SHP1, pODEProblem[24] = CISEqcOE, pODEProblem[25] = EpoRActJAK2, pODEProblem[26] = SOCS3RNAEqc, pODEProblem[27] = CISEqc, pODEProblem[28] = SHP1ProOE, pODEProblem[29] = SOCS3RNADelay, pODEProblem[30] = init_STAT5, pODEProblem[31] = CISoe, pODEProblem[32] = CISRNATurn, pODEProblem[33] = init_SHP1_multiplier, pODEProblem[34] = init_EpoRJAK2, pODEProblem[35] = nuc, pODEProblem[36] = EpoRCISInh, pODEProblem[37] = STAT5ActJAK2, pODEProblem[38] = SOCS3RNATurn, pODEProblem[39] = SOCS3Eqc

	p1EpoRpJAK2 = 0.0 
	pSTAT5 = 0.0 
	EpoRJAK2_CIS = pODEProblem[10] 
	SOCS3nRNA4 = 0.0 
	SOCS3RNA = 0.0 
	SHP1 = pODEProblem[23] * ( 1 + pODEProblem[28] * pODEProblem[33] ) 
	STAT5 = pODEProblem[30] 
	EpoRJAK2 = pODEProblem[34] 
	CISnRNA1 = 0.0 
	SOCS3nRNA1 = 0.0 
	SOCS3nRNA2 = 0.0 
	CISnRNA3 = 0.0 
	CISnRNA4 = 0.0 
	SOCS3 = pODEProblem[39] * pODEProblem[21] * pODEProblem[3] 
	CISnRNA5 = 0.0 
	SOCS3nRNA5 = 0.0 
	SOCS3nRNA3 = 0.0 
	SHP1Act = 0.0 
	npSTAT5 = 0.0 
	p12EpoRpJAK2 = 0.0 
	p2EpoRpJAK2 = 0.0 
	CIS = pODEProblem[27] * pODEProblem[24] * pODEProblem[13] 
	EpoRpJAK2 = 0.0 
	CISnRNA2 = 0.0 
	CISRNA = 0.0 

	u0 .= p1EpoRpJAK2, pSTAT5, EpoRJAK2_CIS, SOCS3nRNA4, SOCS3RNA, SHP1, STAT5, EpoRJAK2, CISnRNA1, SOCS3nRNA1, SOCS3nRNA2, CISnRNA3, CISnRNA4, SOCS3, CISnRNA5, SOCS3nRNA5, SOCS3nRNA3, SHP1Act, npSTAT5, p12EpoRpJAK2, p2EpoRpJAK2, CIS, EpoRpJAK2, CISnRNA2, CISRNA
end

function compute_u0(pODEProblem::AbstractVector)::AbstractVector 

	#pODEProblem[1] = STAT5Exp, pODEProblem[2] = STAT5Imp, pODEProblem[3] = init_SOCS3_multiplier, pODEProblem[4] = EpoRCISRemove, pODEProblem[5] = STAT5ActEpoR, pODEProblem[6] = SHP1ActEpoR, pODEProblem[7] = JAK2EpoRDeaSHP1, pODEProblem[8] = CISTurn, pODEProblem[9] = SOCS3Turn, pODEProblem[10] = init_EpoRJAK2_CIS, pODEProblem[11] = SOCS3Inh, pODEProblem[12] = ActD, pODEProblem[13] = init_CIS_multiplier, pODEProblem[14] = cyt, pODEProblem[15] = CISRNAEqc, pODEProblem[16] = JAK2ActEpo, pODEProblem[17] = Epo, pODEProblem[18] = SOCS3oe, pODEProblem[19] = CISInh, pODEProblem[20] = SHP1Dea, pODEProblem[21] = SOCS3EqcOE, pODEProblem[22] = CISRNADelay, pODEProblem[23] = init_SHP1, pODEProblem[24] = CISEqcOE, pODEProblem[25] = EpoRActJAK2, pODEProblem[26] = SOCS3RNAEqc, pODEProblem[27] = CISEqc, pODEProblem[28] = SHP1ProOE, pODEProblem[29] = SOCS3RNADelay, pODEProblem[30] = init_STAT5, pODEProblem[31] = CISoe, pODEProblem[32] = CISRNATurn, pODEProblem[33] = init_SHP1_multiplier, pODEProblem[34] = init_EpoRJAK2, pODEProblem[35] = nuc, pODEProblem[36] = EpoRCISInh, pODEProblem[37] = STAT5ActJAK2, pODEProblem[38] = SOCS3RNATurn, pODEProblem[39] = SOCS3Eqc

	p1EpoRpJAK2 = 0.0 
	pSTAT5 = 0.0 
	EpoRJAK2_CIS = pODEProblem[10] 
	SOCS3nRNA4 = 0.0 
	SOCS3RNA = 0.0 
	SHP1 = pODEProblem[23] * ( 1 + pODEProblem[28] * pODEProblem[33] ) 
	STAT5 = pODEProblem[30] 
	EpoRJAK2 = pODEProblem[34] 
	CISnRNA1 = 0.0 
	SOCS3nRNA1 = 0.0 
	SOCS3nRNA2 = 0.0 
	CISnRNA3 = 0.0 
	CISnRNA4 = 0.0 
	SOCS3 = pODEProblem[39] * pODEProblem[21] * pODEProblem[3] 
	CISnRNA5 = 0.0 
	SOCS3nRNA5 = 0.0 
	SOCS3nRNA3 = 0.0 
	SHP1Act = 0.0 
	npSTAT5 = 0.0 
	p12EpoRpJAK2 = 0.0 
	p2EpoRpJAK2 = 0.0 
	CIS = pODEProblem[27] * pODEProblem[24] * pODEProblem[13] 
	EpoRpJAK2 = 0.0 
	CISnRNA2 = 0.0 
	CISRNA = 0.0 

	 return [p1EpoRpJAK2, pSTAT5, EpoRJAK2_CIS, SOCS3nRNA4, SOCS3RNA, SHP1, STAT5, EpoRJAK2, CISnRNA1, SOCS3nRNA1, SOCS3nRNA2, CISnRNA3, CISnRNA4, SOCS3, CISnRNA5, SOCS3nRNA5, SOCS3nRNA3, SHP1Act, npSTAT5, p12EpoRpJAK2, p2EpoRpJAK2, CIS, EpoRpJAK2, CISnRNA2, CISRNA]
end

function compute_σ(u::AbstractVector, t::Real, θ_sd::AbstractVector, θ_dynamic::AbstractVector, θ_nonDynamic::AbstractVector, 
                   parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :observable_CISRNA_foldA 
		noiseParameter1_observable_CISRNA_foldA = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_CISRNA_foldA 
	end

	if observableId == :observable_CISRNA_foldB 
		noiseParameter1_observable_CISRNA_foldB = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_CISRNA_foldB 
	end

	if observableId == :observable_CISRNA_foldC 
		noiseParameter1_observable_CISRNA_foldC = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_CISRNA_foldC 
	end

	if observableId == :observable_CIS_abs 
		noiseParameter1_observable_CIS_abs = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_CIS_abs 
	end

	if observableId == :observable_CIS_au 
		noiseParameter1_observable_CIS_au = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_CIS_au 
	end

	if observableId == :observable_CIS_au1 
		noiseParameter1_observable_CIS_au1 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_CIS_au1 
	end

	if observableId == :observable_CIS_au2 
		noiseParameter1_observable_CIS_au2 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_CIS_au2 
	end

	if observableId == :observable_SHP1_abs 
		noiseParameter1_observable_SHP1_abs = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_SHP1_abs 
	end

	if observableId == :observable_SOCS3RNA_foldA 
		noiseParameter1_observable_SOCS3RNA_foldA = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_SOCS3RNA_foldA 
	end

	if observableId == :observable_SOCS3RNA_foldB 
		noiseParameter1_observable_SOCS3RNA_foldB = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_SOCS3RNA_foldB 
	end

	if observableId == :observable_SOCS3RNA_foldC 
		noiseParameter1_observable_SOCS3RNA_foldC = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_SOCS3RNA_foldC 
	end

	if observableId == :observable_SOCS3_abs 
		noiseParameter1_observable_SOCS3_abs = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_SOCS3_abs 
	end

	if observableId == :observable_SOCS3_au 
		noiseParameter1_observable_SOCS3_au = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_SOCS3_au 
	end

	if observableId == :observable_STAT5_abs 
		noiseParameter1_observable_STAT5_abs = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_STAT5_abs 
	end

	if observableId == :observable_pEpoR_au 
		noiseParameter1_observable_pEpoR_au = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_pEpoR_au 
	end

	if observableId == :observable_pJAK2_au 
		noiseParameter1_observable_pJAK2_au = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_pJAK2_au 
	end

	if observableId == :observable_pSTAT5B_rel 
		noiseParameter1_observable_pSTAT5B_rel = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_pSTAT5B_rel 
	end

	if observableId == :observable_pSTAT5_au 
		noiseParameter1_observable_pSTAT5_au, noiseParameter2_observable_pSTAT5_au = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_pSTAT5_au + noiseParameter2_observable_pSTAT5_au 
	end

	if observableId == :observable_tSHP1_au 
		noiseParameter1_observable_tSHP1_au = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_tSHP1_au 
	end

	if observableId == :observable_tSTAT5_au 
		noiseParameter1_observable_tSTAT5_au = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_tSTAT5_au 
	end

end