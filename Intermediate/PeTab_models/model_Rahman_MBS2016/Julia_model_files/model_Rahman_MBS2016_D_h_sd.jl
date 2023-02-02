#u[1] = treated_weak, u[2] = infected_moderate, u[3] = treated_normal, u[4] = infected_normal, u[5] = infected_weak, u[6] = treated_moderate, u[7] = susceptible
#pODEProblem[1] = c1, pODEProblem[2] = treated_normal_death_rate, pODEProblem[3] = infected_normal_transmission_rate_relative, pODEProblem[4] = infected_normal_treatment_rate, pODEProblem[5] = infected_moderate_treatment_rate, pODEProblem[6] = infected_weak_transmission_rate_relative, pODEProblem[7] = infected_moderate_transmission_rate, pODEProblem[8] = infected_normal_worsen_rate, pODEProblem[9] = susceptible_death_rate, pODEProblem[10] = treated_moderate_improve_rate, pODEProblem[11] = infected_weak_treatment_rate, pODEProblem[12] = infected_normal_death_rate, pODEProblem[13] = treated_weak_death_rate, pODEProblem[14] = treated_moderate_death_rate, pODEProblem[15] = recruitment_rate, pODEProblem[16] = infected_moderate_death_rate, pODEProblem[17] = behavioural_change_rate, pODEProblem[18] = infected_moderate_worsen_rate, pODEProblem[19] = treated_weak_improve_rate, pODEProblem[20] = infected_weak_death_rate
#
function compute_∂h∂u!(u, t::Real, pODEProblem::AbstractVector, θ_observable::AbstractVector,                    
                       θ_nonDynamic::AbstractVector, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :observable_prevalence 
		out[1] = (100u[7]) / ((u[2] + u[4] + u[5] + u[7] + u[6] + u[3] + u[1])^2)
		out[2] = (100u[7]) / ((u[2] + u[4] + u[5] + u[7] + u[6] + u[3] + u[1])^2)
		out[3] = (100u[7]) / ((u[2] + u[4] + u[5] + u[7] + u[6] + u[3] + u[1])^2)
		out[4] = (100u[7]) / ((u[2] + u[4] + u[5] + u[7] + u[6] + u[3] + u[1])^2)
		out[5] = (100u[7]) / ((u[2] + u[4] + u[5] + u[7] + u[6] + u[3] + u[1])^2)
		out[6] = (100u[7]) / ((u[2] + u[4] + u[5] + u[7] + u[6] + u[3] + u[1])^2)
		out[7] = (-100.0u[2] - 100.0u[4] - 100.0u[5] - 100.0u[6] - 100.0u[3] - 100.0u[1]) / ((u[2] + u[4] + u[5] + u[7] + u[6] + u[3] + u[1])^2)
		return nothing
	end

end

function compute_∂h∂p!(u, t::Real, pODEProblem::AbstractVector, θ_observable::AbstractVector,
                       θ_nonDynamic::AbstractVector, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :observable_prevalence 
		return nothing
	end

end

function compute_∂σ∂σu!(u, t::Real, θ_sd::AbstractVector, pODEProblem::AbstractVector, θ_nonDynamic::AbstractVector, 
                        parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :observable_prevalence 
		return nothing
	end

end

function compute_∂σ∂σp!(u, t::Real, θ_sd::AbstractVector, pODEProblem::AbstractVector, θ_nonDynamic::AbstractVector, 
                        parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :observable_prevalence 
		return nothing
	end

end

