function evalYmod(u, t, dynPar, obsPar, nonDynParam, paramData, observableId, mapObsParam) 

	Hospital, Symptomatic, Cumulative_cases, Asymptomatic, Exposed, ICU, Recovered, Deaths, Susceptible= u 
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

	asymptomatic_init_concentration, beta_2_multiplier, t_2, gamma_u, exposed_init_concentration, omega_u, kappa, h_hosp_rate, xi, delta_, t_1, beta_0, symptomatic_init_concentration, Interior, mu_u, omega_h, mu_h, beta_2, beta_1, population, p_symp_rate, gamma_s, gamma_h, gamma_a = paramVec 

	Hospital = 0.0 
	Symptomatic = population * symptomatic_init_concentration 
	Cumulative_cases = population * symptomatic_init_concentration 
	Asymptomatic = asymptomatic_init_concentration * population 
	Exposed = exposed_init_concentration * population 
	ICU = 0.0 
	Recovered = 0.0 
	Deaths = 0.0 
	Susceptible = population * ( 1 - asymptomatic_init_concentration - exposed_init_concentration - symptomatic_init_concentration ) 

	u0Vec .= Hospital, Symptomatic, Cumulative_cases, Asymptomatic, Exposed, ICU, Recovered, Deaths, Susceptible
end

function evalSd!(u, t, sdPar, dynPar, nonDynParam, paramData, observableId, mapSdParam) 

	Hospital, Symptomatic, Cumulative_cases, Asymptomatic, Exposed, ICU, Recovered, Deaths, Susceptible= u 
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
		noiseParameter1_cumulative_deaths = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_cumulative_deaths 
	end

	if observableId == "cumulative_cases" 
		noiseParameter1_cumulative_cases = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_cumulative_cases 
	end

end