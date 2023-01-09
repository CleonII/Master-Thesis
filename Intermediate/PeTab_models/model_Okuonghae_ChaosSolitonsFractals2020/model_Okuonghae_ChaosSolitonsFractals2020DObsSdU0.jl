function evalDYmodDu(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	detected_cumulative, symptomatic, asymptomatic, detected, exposed, deceased, recovered, susceptible= u 
	psi, nu, gamma_i, asymptomatic_start, delta, symptomatic_start, sigma, theta, Lagos, alpha, d_0, transmission_rate_effective, eps, exposed_start, gamma_0, d_D, gamma_a = p 
	if observableId == "observable_Cumulative" 
		out[1] = 1
		return nothing
	end

	if observableId == "observable_ActiveCases" 
		out[4] = 1
		return nothing
	end

end

function evalDYmodDp(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	detected_cumulative, symptomatic, asymptomatic, detected, exposed, deceased, recovered, susceptible= u 
	psi, nu, gamma_i, asymptomatic_start, delta, symptomatic_start, sigma, theta, Lagos, alpha, d_0, transmission_rate_effective, eps, exposed_start, gamma_0, d_D, gamma_a = p 
	if observableId == "observable_Cumulative" 
		return nothing
	end

	if observableId == "observable_ActiveCases" 
		return nothing
	end

end

function evalDSdDu!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	detected_cumulative, symptomatic, asymptomatic, detected, exposed, deceased, recovered, susceptible= u 
	psi, nu, gamma_i, asymptomatic_start, delta, symptomatic_start, sigma, theta, Lagos, alpha, d_0, transmission_rate_effective, eps, exposed_start, gamma_0, d_D, gamma_a = p 
	if observableId == "observable_Cumulative" 
		return nothing
	end

	if observableId == "observable_ActiveCases" 
		return nothing
	end

end

function evalDSdDp!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	detected_cumulative, symptomatic, asymptomatic, detected, exposed, deceased, recovered, susceptible= u 
	psi, nu, gamma_i, asymptomatic_start, delta, symptomatic_start, sigma, theta, Lagos, alpha, d_0, transmission_rate_effective, eps, exposed_start, gamma_0, d_D, gamma_a = p 
	if observableId == "observable_Cumulative" 
		return nothing
	end

	if observableId == "observable_ActiveCases" 
		return nothing
	end

end

