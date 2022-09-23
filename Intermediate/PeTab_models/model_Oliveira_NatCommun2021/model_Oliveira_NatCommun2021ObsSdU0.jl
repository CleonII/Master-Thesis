function evalYmod(u, t, dynPar, obsPar, paramData, obsData, observableId, simulationId) 

	Hospital, Symptomatic, Cumulative_cases, Asymptomatic, Exposed, ICU, Recovered, Deaths, Susceptible, dummyVariable= u 
	beta_0, beta_1, beta_2, t_1, t_2, delta_, h_hosp_rate, gamma_h, gamma_u, exposed_init_concentration, asymptomatic_init_concentration, symptomatic_init_concentration = dynPar 
	p_symp_rate_C = paramData.paramVal[1] 
	kappa_C = paramData.paramVal[2] 
	gamma_a_C = paramData.paramVal[3] 
	gamma_s_C = paramData.paramVal[4] 
	xi_C = paramData.paramVal[5] 
	mu_h_C = paramData.paramVal[6] 
	mu_u_C = paramData.paramVal[7] 
	omega_h_C = paramData.paramVal[8] 
	omega_u_C = paramData.paramVal[9] 
	sd_cumulative_deaths_C = paramData.paramVal[19] 
	sd_cumulative_cases_C = paramData.paramVal[20] 
	population_C = paramData.paramVal[21] 

	if observableId == "cumulative_deaths" 
		return Deaths 
	end

	if observableId == "cumulative_cases" 
		return Cumulative_cases 
	end

end

function evalU0!(u0Vec, paramVec) 

	gamma_u, xi, gamma_s, omega_u, mu_h, gamma_h, omega_h, h_hosp_rate, mu_u, kappa, p_symp_rate, Interior, gamma_a, population, beta_2, delta_, beta_0, beta_1, beta_2_multiplier, t_1, t_2, exposed_init_concentration, asymptomatic_init_concentration, symptomatic_init_concentration = paramVec 

	Hospital = 0.0 
	Symptomatic = population * symptomatic_init_concentration 
	Cumulative_cases = population * symptomatic_init_concentration 
	Asymptomatic = asymptomatic_init_concentration * population 
	Exposed = exposed_init_concentration * population 
	ICU = 0.0 
	Recovered = 0.0 
	Deaths = 0.0 
	Susceptible = population * ( 1 - asymptomatic_init_concentration - exposed_init_concentration - symptomatic_init_concentration ) 
	dummyVariable = 0.0 

	u0Vec .= Hospital, Symptomatic, Cumulative_cases, Asymptomatic, Exposed, ICU, Recovered, Deaths, Susceptible, dummyVariable
end

function evalSd!(u, t, sdPar, dynPar, paramData, obsData, observableId, simulationId) 

	Hospital, Symptomatic, Cumulative_cases, Asymptomatic, Exposed, ICU, Recovered, Deaths, Susceptible, dummyVariable= u 
	beta_0, beta_1, beta_2, t_1, t_2, delta_, h_hosp_rate, gamma_h, gamma_u, exposed_init_concentration, asymptomatic_init_concentration, symptomatic_init_concentration = dynPar 
	p_symp_rate_C = paramData.paramVal[1] 
	kappa_C = paramData.paramVal[2] 
	gamma_a_C = paramData.paramVal[3] 
	gamma_s_C = paramData.paramVal[4] 
	xi_C = paramData.paramVal[5] 
	mu_h_C = paramData.paramVal[6] 
	mu_u_C = paramData.paramVal[7] 
	omega_h_C = paramData.paramVal[8] 
	omega_u_C = paramData.paramVal[9] 
	sd_cumulative_deaths_C = paramData.paramVal[19] 
	sd_cumulative_cases_C = paramData.paramVal[20] 
	population_C = paramData.paramVal[21] 

	if observableId == "cumulative_deaths" 
		noiseParameter1_cumulative_deaths = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_cumulative_deaths 
	end

	if observableId == "cumulative_cases" 
		noiseParameter1_cumulative_cases = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_cumulative_cases 
	end

end