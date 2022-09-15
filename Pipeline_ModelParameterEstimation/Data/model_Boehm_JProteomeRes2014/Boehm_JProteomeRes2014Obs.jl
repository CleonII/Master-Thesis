function Boehm_JProteomeRes2014(u, t, dynPar, obsPar, paramData, obsData, observableId, simulationId) 

	STAT5A, pApA, nucpApB, nucpBpB, STAT5B, pApB, nucpApA, pBpB, dummyVariable= u 
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

function Boehm_JProteomeRes2014_t0!(u0Vec, paramVec) 

	k_exp_hetero, Epo_degradation_BaF3, k_exp_homo, k_phos, k_imp_homo, k_imp_hetero, specC17, ratio = paramVec 

	STAT5A = 207.6 * ratio 
	pApA = 0.0 
	nucpApB = 0.0 
	nucpBpB = 0.0 
	STAT5B = 207.6 - 207.6 * ratio 
	pApB = 0.0 
	nucpApA = 0.0 
	pBpB = 0.0 
	dummyVariable = 0.0 

	u0Vec .= STAT5A, pApA, nucpApB, nucpBpB, STAT5B, pApB, nucpApA, pBpB, dummyVariable
end

function Boehm_JProteomeRes2014_sd!(u, t, sdPar, dynPar, paramData, obsData, observableId, simulationId) 

	STAT5A, pApA, nucpApB, nucpBpB, STAT5B, pApB, nucpApA, pBpB, dummyVariable= u 
	Epo_degradation_BaF3, k_exp_hetero, k_exp_homo, k_imp_hetero, k_imp_homo, k_phos = dynPar 
	ratio_C = paramData.paramVal[7] 
	specC17_C = paramData.paramVal[11] 

	if observableId == "pSTAT5A_rel" 
		noiseParameter1_pSTAT5A_rel = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_pSTAT5A_rel 
	end

	if observableId == "pSTAT5B_rel" 
		noiseParameter1_pSTAT5B_rel = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_pSTAT5B_rel 
	end

	if observableId == "rSTAT5A_rel" 
		noiseParameter1_rSTAT5A_rel = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_rSTAT5A_rel 
	end

end