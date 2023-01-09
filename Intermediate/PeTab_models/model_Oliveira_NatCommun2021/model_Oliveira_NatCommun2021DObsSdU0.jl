function evalDYmodDu(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	Hospital, Symptomatic, Cumulative_cases, Asymptomatic, Exposed, ICU, Recovered, Deaths, Susceptible= u 
	asymptomatic_init_concentration, beta_2_multiplier, t_2, gamma_u, exposed_init_concentration, omega_u, kappa, h_hosp_rate, xi, delta_, t_1, beta_0, symptomatic_init_concentration, Interior, mu_u, omega_h, beta_bool1, mu_h, beta_2, beta_1, beta_bool2, population, p_symp_rate, gamma_s, gamma_h, gamma_a = p 
	if observableId == "cumulative_deaths" 
		out[8] = 1
		return nothing
	end

	if observableId == "cumulative_cases" 
		out[3] = 1
		return nothing
	end

end

function evalDYmodDp(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	Hospital, Symptomatic, Cumulative_cases, Asymptomatic, Exposed, ICU, Recovered, Deaths, Susceptible= u 
	asymptomatic_init_concentration, beta_2_multiplier, t_2, gamma_u, exposed_init_concentration, omega_u, kappa, h_hosp_rate, xi, delta_, t_1, beta_0, symptomatic_init_concentration, Interior, mu_u, omega_h, beta_bool1, mu_h, beta_2, beta_1, beta_bool2, population, p_symp_rate, gamma_s, gamma_h, gamma_a = p 
	if observableId == "cumulative_deaths" 
		return nothing
	end

	if observableId == "cumulative_cases" 
		return nothing
	end

end

function evalDSdDu!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	Hospital, Symptomatic, Cumulative_cases, Asymptomatic, Exposed, ICU, Recovered, Deaths, Susceptible= u 
	asymptomatic_init_concentration, beta_2_multiplier, t_2, gamma_u, exposed_init_concentration, omega_u, kappa, h_hosp_rate, xi, delta_, t_1, beta_0, symptomatic_init_concentration, Interior, mu_u, omega_h, beta_bool1, mu_h, beta_2, beta_1, beta_bool2, population, p_symp_rate, gamma_s, gamma_h, gamma_a = p 
	if observableId == "cumulative_deaths" 
		return nothing
	end

	if observableId == "cumulative_cases" 
		return nothing
	end

end

function evalDSdDp!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	Hospital, Symptomatic, Cumulative_cases, Asymptomatic, Exposed, ICU, Recovered, Deaths, Susceptible= u 
	asymptomatic_init_concentration, beta_2_multiplier, t_2, gamma_u, exposed_init_concentration, omega_u, kappa, h_hosp_rate, xi, delta_, t_1, beta_0, symptomatic_init_concentration, Interior, mu_u, omega_h, beta_bool1, mu_h, beta_2, beta_1, beta_bool2, population, p_symp_rate, gamma_s, gamma_h, gamma_a = p 
	if observableId == "cumulative_deaths" 
		return nothing
	end

	if observableId == "cumulative_cases" 
		return nothing
	end

end

