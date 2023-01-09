function evalDYmodDu(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	pAkt_S6, pAkt, pS6, EGFR, pEGFR_Akt, pEGFR, Akt, S6, EGF_EGFR= u 
	EGF_end, reaction_5_k1, reaction_2_k2, init_AKT, init_EGFR, EGF_bool1, EGF_rate, EGFR_turnover, reaction_1_k1, reaction_1_k2, reaction_8_k1, reaction_4_k1, reaction_6_k1, reaction_2_k1, init_S6, reaction_7_k1, reaction_9_k1, reaction_3_k1, reaction_5_k2, Cell, EGF_0 = p 
	if observableId == "pAkt_tot" 
		observableParameter1_pAkt_tot = getObsOrSdParam(obsPar, mapObsParam)
		out[1] = observableParameter1_pAkt_tot
		out[2] = observableParameter1_pAkt_tot
		return nothing
	end

	if observableId == "pEGFR_tot" 
		observableParameter1_pEGFR_tot = getObsOrSdParam(obsPar, mapObsParam)
		out[5] = observableParameter1_pEGFR_tot
		out[6] = observableParameter1_pEGFR_tot
		return nothing
	end

	if observableId == "pS6_tot" 
		observableParameter1_pS6_tot = getObsOrSdParam(obsPar, mapObsParam)
		out[3] = observableParameter1_pS6_tot
		return nothing
	end

end

function evalDYmodDp(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	pAkt_S6, pAkt, pS6, EGFR, pEGFR_Akt, pEGFR, Akt, S6, EGF_EGFR= u 
	EGF_end, reaction_5_k1, reaction_2_k2, init_AKT, init_EGFR, EGF_bool1, EGF_rate, EGFR_turnover, reaction_1_k1, reaction_1_k2, reaction_8_k1, reaction_4_k1, reaction_6_k1, reaction_2_k1, init_S6, reaction_7_k1, reaction_9_k1, reaction_3_k1, reaction_5_k2, Cell, EGF_0 = p 
	if observableId == "pAkt_tot" 
		return nothing
	end

	if observableId == "pEGFR_tot" 
		return nothing
	end

	if observableId == "pS6_tot" 
		return nothing
	end

end

function evalDSdDu!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	pAkt_S6, pAkt, pS6, EGFR, pEGFR_Akt, pEGFR, Akt, S6, EGF_EGFR= u 
	EGF_end, reaction_5_k1, reaction_2_k2, init_AKT, init_EGFR, EGF_bool1, EGF_rate, EGFR_turnover, reaction_1_k1, reaction_1_k2, reaction_8_k1, reaction_4_k1, reaction_6_k1, reaction_2_k1, init_S6, reaction_7_k1, reaction_9_k1, reaction_3_k1, reaction_5_k2, Cell, EGF_0 = p 
	if observableId == "pAkt_tot" 
		return nothing
	end

	if observableId == "pEGFR_tot" 
		return nothing
	end

	if observableId == "pS6_tot" 
		return nothing
	end

end

function evalDSdDp!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	pAkt_S6, pAkt, pS6, EGFR, pEGFR_Akt, pEGFR, Akt, S6, EGF_EGFR= u 
	EGF_end, reaction_5_k1, reaction_2_k2, init_AKT, init_EGFR, EGF_bool1, EGF_rate, EGFR_turnover, reaction_1_k1, reaction_1_k2, reaction_8_k1, reaction_4_k1, reaction_6_k1, reaction_2_k1, init_S6, reaction_7_k1, reaction_9_k1, reaction_3_k1, reaction_5_k2, Cell, EGF_0 = p 
	if observableId == "pAkt_tot" 
		return nothing
	end

	if observableId == "pEGFR_tot" 
		return nothing
	end

	if observableId == "pS6_tot" 
		return nothing
	end

end

