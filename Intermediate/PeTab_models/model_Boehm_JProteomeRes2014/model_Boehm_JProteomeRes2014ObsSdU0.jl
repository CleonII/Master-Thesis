function evalYmod(u, t, dynPar, obsPar, nonDynParam, paramData, observableId, mapObsParam) 

	STAT5A, pApA, nucpApB, nucpBpB, STAT5B, pApB, nucpApA, pBpB= u 
	Epo_degradation_BaF3, k_exp_hetero, k_exp_homo, k_imp_hetero, k_imp_homo, k_phos = dynPar 
	ratio_C = paramData.paramVal[7] 
	specC17_C = paramData.paramVal[11] 

	if observableId == "pSTAT5A_rel" 
		return ( 100 * pApB + 200 * pApA * specC17_C ) / ( pApB + STAT5A * specC17_C + 2 * pApA * specC17_C ) 
	end

	if observableId == "pSTAT5B_rel" 
		return - ( 100 * pApB - 200 * pBpB * ( specC17_C - 1 ) ) / ( ( STAT5B * ( specC17_C - 1 ) - pApB ) + 2 * pBpB * ( specC17_C - 1 ) ) 
	end

	if observableId == "rSTAT5A_rel" 
		return ( 100 * pApB + 100 * STAT5A * specC17_C + 200 * pApA * specC17_C ) / ( 2 * pApB + STAT5A * specC17_C + 2 * pApA * specC17_C - STAT5B * ( specC17_C - 1 ) - 2 * pBpB * ( specC17_C - 1 ) ) 
	end

end

function evalU0!(u0Vec, paramVec) 

	ratio, k_imp_homo, k_exp_hetero, cyt, k_phos, specC17, Epo_degradation_BaF3, k_exp_homo, nuc, k_imp_hetero = paramVec 

	STAT5A = 207.6 * ratio 
	pApA = 0.0 
	nucpApB = 0.0 
	nucpBpB = 0.0 
	STAT5B = 207.6 - 207.6 * ratio 
	pApB = 0.0 
	nucpApA = 0.0 
	pBpB = 0.0 

	u0Vec .= STAT5A, pApA, nucpApB, nucpBpB, STAT5B, pApB, nucpApA, pBpB
end

function evalU0(paramVec) 

	ratio, k_imp_homo, k_exp_hetero, cyt, k_phos, specC17, Epo_degradation_BaF3, k_exp_homo, nuc, k_imp_hetero = paramVec 

	STAT5A = 207.6 * ratio 
	pApA = 0.0 
	nucpApB = 0.0 
	nucpBpB = 0.0 
	STAT5B = 207.6 - 207.6 * ratio 
	pApB = 0.0 
	nucpApA = 0.0 
	pBpB = 0.0 

	 return [STAT5A, pApA, nucpApB, nucpBpB, STAT5B, pApB, nucpApA, pBpB]
end

function evalSd!(u, t, sdPar, dynPar, nonDynParam, paramData, observableId, mapSdParam) 

	STAT5A, pApA, nucpApB, nucpBpB, STAT5B, pApB, nucpApA, pBpB= u 
	Epo_degradation_BaF3, k_exp_hetero, k_exp_homo, k_imp_hetero, k_imp_homo, k_phos = dynPar 
	ratio_C = paramData.paramVal[7] 
	specC17_C = paramData.paramVal[11] 

	if observableId == "pSTAT5A_rel" 
		noiseParameter1_pSTAT5A_rel = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pSTAT5A_rel 
	end

	if observableId == "pSTAT5B_rel" 
		noiseParameter1_pSTAT5B_rel = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_pSTAT5B_rel 
	end

	if observableId == "rSTAT5A_rel" 
		noiseParameter1_rSTAT5A_rel = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_rSTAT5A_rel 
	end

end