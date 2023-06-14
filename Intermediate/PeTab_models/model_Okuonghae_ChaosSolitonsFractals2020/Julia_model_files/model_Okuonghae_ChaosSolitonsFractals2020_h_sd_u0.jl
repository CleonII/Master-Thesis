#u[1] = detected_cumulative, u[2] = symptomatic, u[3] = asymptomatic, u[4] = detected, u[5] = exposed, u[6] = deceased, u[7] = recovered, u[8] = susceptible
#pODEProblemNames[1] = psi, pODEProblemNames[2] = nu, pODEProblemNames[3] = gamma_i, pODEProblemNames[4] = asymptomatic_start, pODEProblemNames[5] = delta, pODEProblemNames[6] = symptomatic_start, pODEProblemNames[7] = sigma, pODEProblemNames[8] = theta, pODEProblemNames[9] = Lagos, pODEProblemNames[10] = alpha, pODEProblemNames[11] = d_0, pODEProblemNames[12] = transmission_rate_effective, pODEProblemNames[13] = eps, pODEProblemNames[14] = exposed_start, pODEProblemNames[15] = gamma_0, pODEProblemNames[16] = d_D, pODEProblemNames[17] = gamma_a
#

function compute_h(u::AbstractVector, t::Real, pODEProblem::AbstractVector, θ_observable::AbstractVector,
                   θ_nonDynamic::AbstractVector, parameterInfo::ParametersInfo, observableId::Symbol, 
                      parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :observable_Cumulative 
		return u[1] 
	end

	if observableId == :observable_ActiveCases 
		return u[4] 
	end

end

function compute_u0!(u0::AbstractVector, pODEProblem::AbstractVector) 

	#pODEProblem[1] = psi, pODEProblem[2] = nu, pODEProblem[3] = gamma_i, pODEProblem[4] = asymptomatic_start, pODEProblem[5] = delta, pODEProblem[6] = symptomatic_start, pODEProblem[7] = sigma, pODEProblem[8] = theta, pODEProblem[9] = Lagos, pODEProblem[10] = alpha, pODEProblem[11] = d_0, pODEProblem[12] = transmission_rate_effective, pODEProblem[13] = eps, pODEProblem[14] = exposed_start, pODEProblem[15] = gamma_0, pODEProblem[16] = d_D, pODEProblem[17] = gamma_a

	detected_cumulative = 1.0 
	symptomatic = pODEProblem[6] 
	asymptomatic = pODEProblem[4] 
	detected = 1.0 
	exposed = pODEProblem[14] 
	deceased = 0.0 
	recovered = 0.0 
	susceptible = 1.4367982e7 

	u0 .= detected_cumulative, symptomatic, asymptomatic, detected, exposed, deceased, recovered, susceptible
end

function compute_u0(pODEProblem::AbstractVector)::AbstractVector 

	#pODEProblem[1] = psi, pODEProblem[2] = nu, pODEProblem[3] = gamma_i, pODEProblem[4] = asymptomatic_start, pODEProblem[5] = delta, pODEProblem[6] = symptomatic_start, pODEProblem[7] = sigma, pODEProblem[8] = theta, pODEProblem[9] = Lagos, pODEProblem[10] = alpha, pODEProblem[11] = d_0, pODEProblem[12] = transmission_rate_effective, pODEProblem[13] = eps, pODEProblem[14] = exposed_start, pODEProblem[15] = gamma_0, pODEProblem[16] = d_D, pODEProblem[17] = gamma_a

	detected_cumulative = 1.0 
	symptomatic = pODEProblem[6] 
	asymptomatic = pODEProblem[4] 
	detected = 1.0 
	exposed = pODEProblem[14] 
	deceased = 0.0 
	recovered = 0.0 
	susceptible = 1.4367982e7 

	 return [detected_cumulative, symptomatic, asymptomatic, detected, exposed, deceased, recovered, susceptible]
end

function compute_σ(u::AbstractVector, t::Real, θ_sd::AbstractVector, pODEProblem::AbstractVector, θ_nonDynamic::AbstractVector, 
                   parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :observable_Cumulative 
		noiseParameter1_observable_Cumulative = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_Cumulative 
	end

	if observableId == :observable_ActiveCases 
		noiseParameter1_observable_ActiveCases = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_ActiveCases 
	end

end