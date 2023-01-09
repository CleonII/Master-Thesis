function evalYmod(u, t, dynPar, obsPar, nonDynParam, paramData, observableId, mapObsParam) 

	pAkt_S6, pAkt, pS6, EGFR, pEGFR_Akt, pEGFR, Akt, S6, EGF_EGFR= u 
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

	EGF_end, reaction_5_k1, reaction_2_k2, init_AKT, init_EGFR, EGF_bool1, EGF_rate, EGFR_turnover, reaction_1_k1, reaction_1_k2, reaction_8_k1, reaction_4_k1, reaction_6_k1, reaction_2_k1, init_S6, reaction_7_k1, reaction_9_k1, reaction_3_k1, reaction_5_k2, Cell, EGF_0 = paramVec 

	pAkt_S6 = 0.0 
	pAkt = 0.0 
	pS6 = 0.0 
	EGFR = init_EGFR 
	pEGFR_Akt = 0.0 
	pEGFR = 0.0 
	Akt = init_AKT 
	S6 = init_S6 
	EGF_EGFR = 0.0 

	u0Vec .= pAkt_S6, pAkt, pS6, EGFR, pEGFR_Akt, pEGFR, Akt, S6, EGF_EGFR
end

function evalU0(paramVec) 

	EGF_end, reaction_5_k1, reaction_2_k2, init_AKT, init_EGFR, EGF_bool1, EGF_rate, EGFR_turnover, reaction_1_k1, reaction_1_k2, reaction_8_k1, reaction_4_k1, reaction_6_k1, reaction_2_k1, init_S6, reaction_7_k1, reaction_9_k1, reaction_3_k1, reaction_5_k2, Cell, EGF_0 = paramVec 

	pAkt_S6 = 0.0 
	pAkt = 0.0 
	pS6 = 0.0 
	EGFR = init_EGFR 
	pEGFR_Akt = 0.0 
	pEGFR = 0.0 
	Akt = init_AKT 
	S6 = init_S6 
	EGF_EGFR = 0.0 

	 return [pAkt_S6, pAkt, pS6, EGFR, pEGFR_Akt, pEGFR, Akt, S6, EGF_EGFR]
end

function evalSd!(u, t, sdPar, dynPar, nonDynParam, paramData, observableId, mapSdParam) 

	pAkt_S6, pAkt, pS6, EGFR, pEGFR_Akt, pEGFR, Akt, S6, EGF_EGFR= u 
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