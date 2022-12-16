function evalDYmodDu(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	Naive, Pathogen, LateEffector, EarlyEffector, Memory= u 
	mu_LL, delta_NE, mu_PE, mu_P, mu_PL, delta_EL, mu_EE, default, mu_N, rho_E, delta_LM, rho_P, mu_LE = p 
	if observableId == "observable_EarlyEffector" 
		out[4] = 1
		return nothing
	end

	if observableId == "observable_LateEffector" 
		out[3] = 1
		return nothing
	end

	if observableId == "observable_Memory" 
		out[5] = 1
		return nothing
	end

	if observableId == "observable_Naive" 
		out[1] = 1
		return nothing
	end

end

function evalDYmodDp(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	Naive, Pathogen, LateEffector, EarlyEffector, Memory= u 
	mu_LL, delta_NE, mu_PE, mu_P, mu_PL, delta_EL, mu_EE, default, mu_N, rho_E, delta_LM, rho_P, mu_LE = p 
	if observableId == "observable_EarlyEffector" 
		return nothing
	end

	if observableId == "observable_LateEffector" 
		return nothing
	end

	if observableId == "observable_Memory" 
		return nothing
	end

	if observableId == "observable_Naive" 
		return nothing
	end

end

function evalDSdDu!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	Naive, Pathogen, LateEffector, EarlyEffector, Memory= u 
	mu_LL, delta_NE, mu_PE, mu_P, mu_PL, delta_EL, mu_EE, default, mu_N, rho_E, delta_LM, rho_P, mu_LE = p 
	if observableId == "observable_EarlyEffector" 
		return nothing
	end

	if observableId == "observable_LateEffector" 
		return nothing
	end

	if observableId == "observable_Memory" 
		return nothing
	end

	if observableId == "observable_Naive" 
		return nothing
	end

end

function evalDSdDp!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	Naive, Pathogen, LateEffector, EarlyEffector, Memory= u 
	mu_LL, delta_NE, mu_PE, mu_P, mu_PL, delta_EL, mu_EE, default, mu_N, rho_E, delta_LM, rho_P, mu_LE = p 
	if observableId == "observable_EarlyEffector" 
		return nothing
	end

	if observableId == "observable_LateEffector" 
		return nothing
	end

	if observableId == "observable_Memory" 
		return nothing
	end

	if observableId == "observable_Naive" 
		return nothing
	end

end

