function evalYmod(u, t, dynPar, obsPar, paramData, observableId, mapObsParam) 

	pAkt_S6, pAkt, pS6, EGFR, pEGFR_Akt, pEGFR, Akt, S6, EGF_EGFR, EGF, dummyVariable= u 
	EGFR_turnover, init_AKT, init_EGFR, init_S6, reaction_1_k1, reaction_1_k2, reaction_2_k1, reaction_2_k2, reaction_3_k1, reaction_4_k1, reaction_5_k1, reaction_5_k2, reaction_6_k1, reaction_7_k1, reaction_8_k1, reaction_9_k1 = dynPar 

	if observableId == "pAkt_tot" 
		observableParameter1_pAkt_tot = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pAkt_tot * ( pAkt + pAkt_S6 ) 
	end

	if observableId == "pEGFR_tot" 
		observableParameter1_pEGFR_tot = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_pEGFR_tot * ( pEGFR + pEGFR_Akt ) 
	end

	if observableId == "pS6_tot" 
		observableParameter1_pS6_tot = getObsOrSdParam(obsPar, mapObsParam)
		return pS6 * observableParameter1_pS6_tot 
	end

end

function evalU0!(u0Vec, paramVec) 

	reaction_5_k2, reaction_5_k1, reaction_6_k1, reaction_3_k1, reaction_7_k1, reaction_8_k1, EGFR_turnover, reaction_1_k2, reaction_1_k1, reaction_2_k2, reaction_2_k1, reaction_4_k1, reaction_9_k1, EGF_rate, EGF_end, init_AKT, init_S6, EGF_0, init_EGFR = paramVec 

	pAkt_S6 = 0.0 
	pAkt = 0.0 
	pS6 = 0.0 
	EGFR = init_EGFR 
	pEGFR_Akt = 0.0 
	pEGFR = 0.0 
	Akt = init_AKT 
	S6 = init_S6 
	EGF_EGFR = 0.0 
	EGF = EGF_0 
	dummyVariable = 0.0 

	u0Vec .= pAkt_S6, pAkt, pS6, EGFR, pEGFR_Akt, pEGFR, Akt, S6, EGF_EGFR, EGF, dummyVariable
end

function evalSd!(u, t, sdPar, dynPar, paramData, observableId, mapSdParam) 

	pAkt_S6, pAkt, pS6, EGFR, pEGFR_Akt, pEGFR, Akt, S6, EGF_EGFR, EGF, dummyVariable= u 
	EGFR_turnover, init_AKT, init_EGFR, init_S6, reaction_1_k1, reaction_1_k2, reaction_2_k1, reaction_2_k2, reaction_3_k1, reaction_4_k1, reaction_5_k1, reaction_5_k2, reaction_6_k1, reaction_7_k1, reaction_8_k1, reaction_9_k1 = dynPar 

	if observableId == "pAkt_tot" 
		noiseParameter1_pAkt_tot = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pAkt_tot 
	end

	if observableId == "pEGFR_tot" 
		noiseParameter1_pEGFR_tot = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pEGFR_tot 
	end

	if observableId == "pS6_tot" 
		noiseParameter1_pS6_tot = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pS6_tot 
	end

end