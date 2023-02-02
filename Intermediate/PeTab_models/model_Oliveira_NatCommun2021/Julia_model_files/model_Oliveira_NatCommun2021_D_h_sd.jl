#u[1] = Hospital, u[2] = Symptomatic, u[3] = Cumulative_cases, u[4] = Asymptomatic, u[5] = Exposed, u[6] = ICU, u[7] = Recovered, u[8] = Deaths, u[9] = Susceptible
#pODEProblem[1] = asymptomatic_init_concentration, pODEProblem[2] = beta_2_multiplier, pODEProblem[3] = t_2, pODEProblem[4] = gamma_u, pODEProblem[5] = exposed_init_concentration, pODEProblem[6] = omega_u, pODEProblem[7] = kappa, pODEProblem[8] = h_hosp_rate, pODEProblem[9] = xi, pODEProblem[10] = delta_, pODEProblem[11] = t_1, pODEProblem[12] = beta_0, pODEProblem[13] = symptomatic_init_concentration, pODEProblem[14] = Interior, pODEProblem[15] = mu_u, pODEProblem[16] = omega_h, pODEProblem[17] = beta_bool1, pODEProblem[18] = mu_h, pODEProblem[19] = beta_2, pODEProblem[20] = beta_1, pODEProblem[21] = beta_bool2, pODEProblem[22] = population, pODEProblem[23] = p_symp_rate, pODEProblem[24] = gamma_s, pODEProblem[25] = gamma_h, pODEProblem[26] = gamma_a
#
function compute_∂h∂u!(u, t::Real, pODEProblem::AbstractVector, θ_observable::AbstractVector,                    
                       θ_nonDynamic::AbstractVector, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :cumulative_deaths 
		out[8] = 1
		return nothing
	end

	if observableId == :cumulative_cases 
		out[3] = 1
		return nothing
	end

end

function compute_∂h∂p!(u, t::Real, pODEProblem::AbstractVector, θ_observable::AbstractVector,
                       θ_nonDynamic::AbstractVector, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :cumulative_deaths 
		return nothing
	end

	if observableId == :cumulative_cases 
		return nothing
	end

end

function compute_∂σ∂σu!(u, t::Real, θ_sd::AbstractVector, pODEProblem::AbstractVector, θ_nonDynamic::AbstractVector, 
                        parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :cumulative_deaths 
		return nothing
	end

	if observableId == :cumulative_cases 
		return nothing
	end

end

function compute_∂σ∂σp!(u, t::Real, θ_sd::AbstractVector, pODEProblem::AbstractVector, θ_nonDynamic::AbstractVector, 
                        parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :cumulative_deaths 
		return nothing
	end

	if observableId == :cumulative_cases 
		return nothing
	end

end

