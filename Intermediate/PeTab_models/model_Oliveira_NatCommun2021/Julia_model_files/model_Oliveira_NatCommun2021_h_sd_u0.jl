#u[1] = Hospital, u[2] = Symptomatic, u[3] = Cumulative_cases, u[4] = Asymptomatic, u[5] = Exposed, u[6] = ICU, u[7] = Recovered, u[8] = Deaths, u[9] = Susceptible
#θ_dynamicNames[1] = beta_0, θ_dynamicNames[2] = beta_1, θ_dynamicNames[3] = beta_2, θ_dynamicNames[4] = t_1, θ_dynamicNames[5] = t_2, θ_dynamicNames[6] = delta_, θ_dynamicNames[7] = h_hosp_rate, θ_dynamicNames[8] = gamma_h, θ_dynamicNames[9] = gamma_u, θ_dynamicNames[10] = exposed_init_concentration, θ_dynamicNames[11] = asymptomatic_init_concentration, θ_dynamicNames[12] = symptomatic_init_concentration
##parameterInfo.nominalValue[1] = p_symp_rate_C 
#parameterInfo.nominalValue[2] = kappa_C 
#parameterInfo.nominalValue[3] = gamma_a_C 
#parameterInfo.nominalValue[4] = gamma_s_C 
#parameterInfo.nominalValue[5] = xi_C 
#parameterInfo.nominalValue[6] = mu_h_C 
#parameterInfo.nominalValue[7] = mu_u_C 
#parameterInfo.nominalValue[8] = omega_h_C 
#parameterInfo.nominalValue[9] = omega_u_C 
#parameterInfo.nominalValue[19] = sd_cumulative_deaths_C 
#parameterInfo.nominalValue[20] = sd_cumulative_cases_C 
#parameterInfo.nominalValue[21] = population_C 


function compute_h(u::AbstractVector, t::Real, θ_dynamic::AbstractVector, θ_observable::AbstractVector,
                   θ_nonDynamic::AbstractVector, parameterInfo::ParametersInfo, observableId::Symbol, 
                      parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :cumulative_deaths 
		return u[8] 
	end

	if observableId == :cumulative_cases 
		return u[3] 
	end

end

function compute_u0!(u0::AbstractVector, pODEProblem::AbstractVector) 

	#pODEProblem[1] = asymptomatic_init_concentration, pODEProblem[2] = beta_2_multiplier, pODEProblem[3] = t_2, pODEProblem[4] = gamma_u, pODEProblem[5] = exposed_init_concentration, pODEProblem[6] = omega_u, pODEProblem[7] = kappa, pODEProblem[8] = h_hosp_rate, pODEProblem[9] = xi, pODEProblem[10] = delta_, pODEProblem[11] = t_1, pODEProblem[12] = beta_0, pODEProblem[13] = symptomatic_init_concentration, pODEProblem[14] = Interior, pODEProblem[15] = mu_u, pODEProblem[16] = omega_h, pODEProblem[17] = beta_bool1, pODEProblem[18] = mu_h, pODEProblem[19] = beta_2, pODEProblem[20] = beta_1, pODEProblem[21] = beta_bool2, pODEProblem[22] = population, pODEProblem[23] = p_symp_rate, pODEProblem[24] = gamma_s, pODEProblem[25] = gamma_h, pODEProblem[26] = gamma_a

	Hospital = 0.0 
	Symptomatic = pODEProblem[22] * pODEProblem[13] 
	Cumulative_cases = pODEProblem[22] * pODEProblem[13] 
	Asymptomatic = pODEProblem[1] * pODEProblem[22] 
	Exposed = pODEProblem[5] * pODEProblem[22] 
	ICU = 0.0 
	Recovered = 0.0 
	Deaths = 0.0 
	Susceptible = pODEProblem[22] * ( 1 - pODEProblem[1] - pODEProblem[5] - pODEProblem[13] ) 

	u0 .= Hospital, Symptomatic, Cumulative_cases, Asymptomatic, Exposed, ICU, Recovered, Deaths, Susceptible
end

function compute_u0(pODEProblem::AbstractVector)::AbstractVector 

	#pODEProblem[1] = asymptomatic_init_concentration, pODEProblem[2] = beta_2_multiplier, pODEProblem[3] = t_2, pODEProblem[4] = gamma_u, pODEProblem[5] = exposed_init_concentration, pODEProblem[6] = omega_u, pODEProblem[7] = kappa, pODEProblem[8] = h_hosp_rate, pODEProblem[9] = xi, pODEProblem[10] = delta_, pODEProblem[11] = t_1, pODEProblem[12] = beta_0, pODEProblem[13] = symptomatic_init_concentration, pODEProblem[14] = Interior, pODEProblem[15] = mu_u, pODEProblem[16] = omega_h, pODEProblem[17] = beta_bool1, pODEProblem[18] = mu_h, pODEProblem[19] = beta_2, pODEProblem[20] = beta_1, pODEProblem[21] = beta_bool2, pODEProblem[22] = population, pODEProblem[23] = p_symp_rate, pODEProblem[24] = gamma_s, pODEProblem[25] = gamma_h, pODEProblem[26] = gamma_a

	Hospital = 0.0 
	Symptomatic = pODEProblem[22] * pODEProblem[13] 
	Cumulative_cases = pODEProblem[22] * pODEProblem[13] 
	Asymptomatic = pODEProblem[1] * pODEProblem[22] 
	Exposed = pODEProblem[5] * pODEProblem[22] 
	ICU = 0.0 
	Recovered = 0.0 
	Deaths = 0.0 
	Susceptible = pODEProblem[22] * ( 1 - pODEProblem[1] - pODEProblem[5] - pODEProblem[13] ) 

	 return [Hospital, Symptomatic, Cumulative_cases, Asymptomatic, Exposed, ICU, Recovered, Deaths, Susceptible]
end

function compute_σ(u::AbstractVector, t::Real, θ_sd::AbstractVector, θ_dynamic::AbstractVector, θ_nonDynamic::AbstractVector, 
                   parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :cumulative_deaths 
		noiseParameter1_cumulative_deaths = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_cumulative_deaths 
	end

	if observableId == :cumulative_cases 
		noiseParameter1_cumulative_cases = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_cumulative_cases 
	end

end