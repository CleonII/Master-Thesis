function evalDYmodDu(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	STAT5A, pApA, nucpApB, nucpBpB, STAT5B, pApB, nucpApA, pBpB= u 
	ratio, k_imp_homo, k_exp_hetero, cyt, k_phos, specC17, Epo_degradation_BaF3, k_exp_homo, nuc, k_imp_hetero = p 
	if observableId == "pSTAT5A_rel" 
		out[1] = (specC17*(-100pApB - 200pApA*specC17)) / ((pApB + STAT5A*specC17 + 2pApA*specC17)^2)
		out[2] = (200.0STAT5A*(specC17^2)) / ((pApB + STAT5A*specC17 + 2pApA*specC17)^2)
		out[6] = (100.0STAT5A*specC17) / ((pApB + STAT5A*specC17 + 2pApA*specC17)^2)
		return nothing
	end

	if observableId == "pSTAT5B_rel" 
		out[5] = ((specC17 - 1)*(100pApB + 200pBpB - 200pBpB*specC17)) / ((STAT5B*specC17 + 2pBpB*specC17 - STAT5B - pApB - 2pBpB)^2)
		out[6] = (100.0STAT5B - 100.0STAT5B*specC17) / ((STAT5B*specC17 + 2pBpB*specC17 - STAT5B - pApB - 2pBpB)^2)
		out[8] = (200.0STAT5B + 200.0STAT5B*(specC17^2) - 400.0STAT5B*specC17) / ((STAT5B*specC17 + 2pBpB*specC17 - STAT5B - pApB - 2pBpB)^2)
		return nothing
	end

	if observableId == "rSTAT5A_rel" 
		out[1] = (100.0STAT5B*specC17 + 100.0pApB*specC17 + 200.0pBpB*specC17 - 100.0STAT5B*(specC17^2) - 200.0pBpB*(specC17^2)) / ((STAT5B + 2pApB + 2pBpB + STAT5A*specC17 + 2pApA*specC17 - STAT5B*specC17 - 2pBpB*specC17)^2)
		out[2] = (200.0STAT5B*specC17 + 200.0pApB*specC17 + 400.0pBpB*specC17 - 200.0STAT5B*(specC17^2) - 400.0pBpB*(specC17^2)) / ((STAT5B + 2pApB + 2pBpB + STAT5A*specC17 + 2pApA*specC17 - STAT5B*specC17 - 2pBpB*specC17)^2)
		out[5] = ((1 - specC17)*(-100pApB - 100STAT5A*specC17 - 200pApA*specC17)) / ((STAT5B + 2pApB + 2pBpB + STAT5A*specC17 + 2pApA*specC17 - STAT5B*specC17 - 2pBpB*specC17)^2)
		out[6] = (100.0STAT5B + 200.0pBpB - 100.0STAT5A*specC17 - 100.0STAT5B*specC17 - 200.0pApA*specC17 - 200.0pBpB*specC17) / ((STAT5B + 2pApB + 2pBpB + STAT5A*specC17 + 2pApA*specC17 - STAT5B*specC17 - 2pBpB*specC17)^2)
		out[8] = ((2specC17 - 2)*(100pApB + 100STAT5A*specC17 + 200pApA*specC17)) / ((STAT5B + 2pApB + 2pBpB + STAT5A*specC17 + 2pApA*specC17 - STAT5B*specC17 - 2pBpB*specC17)^2)
		return nothing
	end

end

function evalDYmodDp(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	STAT5A, pApA, nucpApB, nucpBpB, STAT5B, pApB, nucpApA, pBpB= u 
	ratio, k_imp_homo, k_exp_hetero, cyt, k_phos, specC17, Epo_degradation_BaF3, k_exp_homo, nuc, k_imp_hetero = p 
	if observableId == "pSTAT5A_rel" 
		out[6] = (-100.0STAT5A*pApB) / ((pApB + STAT5A*specC17 + 2pApA*specC17)^2)
		return nothing
	end

	if observableId == "pSTAT5B_rel" 
		out[6] = (100.0STAT5B*pApB) / ((STAT5B*specC17 + 2pBpB*specC17 - STAT5B - pApB - 2pBpB)^2)
		return nothing
	end

	if observableId == "rSTAT5A_rel" 
		out[6] = (100.0STAT5A*STAT5B + 100.0STAT5A*pApB + 100.0STAT5B*pApB + 200.0STAT5A*pBpB + 200.0STAT5B*pApA + 200.0pApA*pApB + 200.0pApB*pBpB + 400.0pApA*pBpB) / ((STAT5B + 2pApB + 2pBpB + STAT5A*specC17 + 2pApA*specC17 - STAT5B*specC17 - 2pBpB*specC17)^2)
		return nothing
	end

end

function evalDSdDu!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	STAT5A, pApA, nucpApB, nucpBpB, STAT5B, pApB, nucpApA, pBpB= u 
	ratio, k_imp_homo, k_exp_hetero, cyt, k_phos, specC17, Epo_degradation_BaF3, k_exp_homo, nuc, k_imp_hetero = p 
	if observableId == "pSTAT5A_rel" 
		return nothing
	end

	if observableId == "pSTAT5B_rel" 
		return nothing
	end

	if observableId == "rSTAT5A_rel" 
		return nothing
	end

end

function evalDSdDp!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	STAT5A, pApA, nucpApB, nucpBpB, STAT5B, pApB, nucpApA, pBpB= u 
	ratio, k_imp_homo, k_exp_hetero, cyt, k_phos, specC17, Epo_degradation_BaF3, k_exp_homo, nuc, k_imp_hetero = p 
	if observableId == "pSTAT5A_rel" 
		return nothing
	end

	if observableId == "pSTAT5B_rel" 
		return nothing
	end

	if observableId == "rSTAT5A_rel" 
		return nothing
	end

end

