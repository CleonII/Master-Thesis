function evalYmod(u, t, dynPar, obsPar, paramData, obsData, observableId, simulationId) 

	detected_cumulative, symptomatic, asymptomatic, detected, exposed, deceased, recovered, susceptible, dummyVariable= u 
	alpha, d_0, d_D, gamma_0, gamma_a, gamma_i, nu, psi, sigma, theta, asymptomatic_start, symptomatic_start, exposed_start, transmission_rate_effective = dynPar 

	if observableId == "observable_Cumulative" 
		return detected_cumulative 
	end

	if observableId == "observable_ActiveCases" 
		return detected 
	end

end

function evalU0!(u0Vec, paramVec) 

	psi, theta, sigma, d_0, nu, gamma_0, gamma_a, d_D, gamma_i, alpha, eps, transmission_rate_effective, delta, exposed_start, asymptomatic_start, symptomatic_start = paramVec 

	detected_cumulative = 1.0 
	symptomatic = symptomatic_start 
	asymptomatic = asymptomatic_start 
	detected = 1.0 
	exposed = exposed_start 
	deceased = 0.0 
	recovered = 0.0 
	susceptible = 1.4367982e7 
	dummyVariable = 0.0 

	u0Vec .= detected_cumulative, symptomatic, asymptomatic, detected, exposed, deceased, recovered, susceptible, dummyVariable
end

function evalSd!(u, t, sdPar, dynPar, paramData, obsData, observableId, simulationId) 

	detected_cumulative, symptomatic, asymptomatic, detected, exposed, deceased, recovered, susceptible, dummyVariable= u 
	alpha, d_0, d_D, gamma_0, gamma_a, gamma_i, nu, psi, sigma, theta, asymptomatic_start, symptomatic_start, exposed_start, transmission_rate_effective = dynPar 

	if observableId == "observable_Cumulative" 
		noiseParameter1_observable_Cumulative = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_Cumulative 
	end

	if observableId == "observable_ActiveCases" 
		noiseParameter1_observable_ActiveCases = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_ActiveCases 
	end

end