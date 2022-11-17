function evalYmod(u, t, dynPar, obsPar, nonDynParam, paramData, observableId, mapObsParam) 

	Naive, Pathogen, LateEffector, EarlyEffector, Memory= u 
	delta_EL, delta_LM, delta_NE, mu_EE, mu_LE, mu_LL, mu_N, mu_P, mu_PE, mu_PL, rho_E, rho_P = dynPar 

	if observableId == "observable_EarlyEffector" 
		return EarlyEffector 
	end

	if observableId == "observable_LateEffector" 
		return LateEffector 
	end

	if observableId == "observable_Memory" 
		return Memory 
	end

	if observableId == "observable_Naive" 
		return Naive 
	end

end

function evalU0!(u0Vec, paramVec) 

	mu_LL, delta_NE, mu_PE, mu_P, mu_PL, delta_EL, mu_EE, default, mu_N, rho_E, delta_LM, rho_P, mu_LE = paramVec 

	Naive = 8090.0 
	Pathogen = 1.0 
	LateEffector = 0.0 
	EarlyEffector = 0.0 
	Memory = 0.0 

	u0Vec .= Naive, Pathogen, LateEffector, EarlyEffector, Memory
end

function evalSd!(u, t, sdPar, dynPar, nonDynParam, paramData, observableId, mapSdParam) 

	Naive, Pathogen, LateEffector, EarlyEffector, Memory= u 
	delta_EL, delta_LM, delta_NE, mu_EE, mu_LE, mu_LL, mu_N, mu_P, mu_PE, mu_PL, rho_E, rho_P = dynPar 

	if observableId == "observable_EarlyEffector" 
		noiseParameter1_observable_EarlyEffector = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_EarlyEffector 
	end

	if observableId == "observable_LateEffector" 
		noiseParameter1_observable_LateEffector = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_LateEffector 
	end

	if observableId == "observable_Memory" 
		noiseParameter1_observable_Memory = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_Memory 
	end

	if observableId == "observable_Naive" 
		noiseParameter1_observable_Naive = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_Naive 
	end

end