#u[1] = detected_cumulative, u[2] = symptomatic, u[3] = asymptomatic, u[4] = detected, u[5] = exposed, u[6] = deceased, u[7] = recovered, u[8] = susceptible
#pODEProblem[1] = psi, pODEProblem[2] = nu, pODEProblem[3] = gamma_i, pODEProblem[4] = asymptomatic_start, pODEProblem[5] = delta, pODEProblem[6] = symptomatic_start, pODEProblem[7] = sigma, pODEProblem[8] = theta, pODEProblem[9] = Lagos, pODEProblem[10] = alpha, pODEProblem[11] = d_0, pODEProblem[12] = transmission_rate_effective, pODEProblem[13] = eps, pODEProblem[14] = exposed_start, pODEProblem[15] = gamma_0, pODEProblem[16] = d_D, pODEProblem[17] = gamma_a
#
function compute_∂h∂u!(u, t::Real, pODEProblem::AbstractVector, θ_observable::AbstractVector,                    
                       θ_nonDynamic::AbstractVector, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :observable_Cumulative 
		out[1] = 1
		return nothing
	end

	if observableId == :observable_ActiveCases 
		out[4] = 1
		return nothing
	end

end

function compute_∂h∂p!(u, t::Real, pODEProblem::AbstractVector, θ_observable::AbstractVector,
                       θ_nonDynamic::AbstractVector, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :observable_Cumulative 
		return nothing
	end

	if observableId == :observable_ActiveCases 
		return nothing
	end

end

function compute_∂σ∂σu!(u, t::Real, θ_sd::AbstractVector, pODEProblem::AbstractVector, θ_nonDynamic::AbstractVector, 
                        parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :observable_Cumulative 
		return nothing
	end

	if observableId == :observable_ActiveCases 
		return nothing
	end

end

function compute_∂σ∂σp!(u, t::Real, θ_sd::AbstractVector, pODEProblem::AbstractVector, θ_nonDynamic::AbstractVector, 
                        parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :observable_Cumulative 
		return nothing
	end

	if observableId == :observable_ActiveCases 
		return nothing
	end

end

