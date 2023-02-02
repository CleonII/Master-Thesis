#u[1] = treated_weak, u[2] = infected_moderate, u[3] = treated_normal, u[4] = infected_normal, u[5] = infected_weak, u[6] = treated_moderate, u[7] = susceptible
#θ_dynamicNames[1] = infected_normal_transmission_rate_relative, θ_dynamicNames[2] = infected_moderate_transmission_rate, θ_dynamicNames[3] = infected_weak_transmission_rate_relative, θ_dynamicNames[4] = infected_weak_treatment_rate, θ_dynamicNames[5] = infected_normal_worsen_rate, θ_dynamicNames[6] = infected_moderate_worsen_rate, θ_dynamicNames[7] = treated_moderate_improve_rate, θ_dynamicNames[8] = treated_weak_improve_rate, θ_dynamicNames[9] = behavioural_change_rate
#

function compute_h(u::AbstractVector, t::Real, θ_dynamic::AbstractVector, θ_observable::AbstractVector,
                   θ_nonDynamic::AbstractVector, parameterInfo::ParametersInfo, observableId::Symbol, 
                      parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :observable_prevalence 
		return ( ( 1  - u[7] /  ( u[7] + u[4] + u[2] + u[5] + u[3] + u[6] + u[1] ) )  * 100 ) 
	end

end

function compute_u0!(u0::AbstractVector, pODEProblem::AbstractVector) 

	#pODEProblem[1] = c1, pODEProblem[2] = treated_normal_death_rate, pODEProblem[3] = infected_normal_transmission_rate_relative, pODEProblem[4] = infected_normal_treatment_rate, pODEProblem[5] = infected_moderate_treatment_rate, pODEProblem[6] = infected_weak_transmission_rate_relative, pODEProblem[7] = infected_moderate_transmission_rate, pODEProblem[8] = infected_normal_worsen_rate, pODEProblem[9] = susceptible_death_rate, pODEProblem[10] = treated_moderate_improve_rate, pODEProblem[11] = infected_weak_treatment_rate, pODEProblem[12] = infected_normal_death_rate, pODEProblem[13] = treated_weak_death_rate, pODEProblem[14] = treated_moderate_death_rate, pODEProblem[15] = recruitment_rate, pODEProblem[16] = infected_moderate_death_rate, pODEProblem[17] = behavioural_change_rate, pODEProblem[18] = infected_moderate_worsen_rate, pODEProblem[19] = treated_weak_improve_rate, pODEProblem[20] = infected_weak_death_rate

	treated_weak = 0.0 
	infected_moderate = 9000.0 
	treated_normal = 0.0 
	infected_normal = 16300.0 
	infected_weak = 11000.0 
	treated_moderate = 0.0 
	susceptible = 1.794e7 

	u0 .= treated_weak, infected_moderate, treated_normal, infected_normal, infected_weak, treated_moderate, susceptible
end

function compute_u0(pODEProblem::AbstractVector)::AbstractVector 

	#pODEProblem[1] = c1, pODEProblem[2] = treated_normal_death_rate, pODEProblem[3] = infected_normal_transmission_rate_relative, pODEProblem[4] = infected_normal_treatment_rate, pODEProblem[5] = infected_moderate_treatment_rate, pODEProblem[6] = infected_weak_transmission_rate_relative, pODEProblem[7] = infected_moderate_transmission_rate, pODEProblem[8] = infected_normal_worsen_rate, pODEProblem[9] = susceptible_death_rate, pODEProblem[10] = treated_moderate_improve_rate, pODEProblem[11] = infected_weak_treatment_rate, pODEProblem[12] = infected_normal_death_rate, pODEProblem[13] = treated_weak_death_rate, pODEProblem[14] = treated_moderate_death_rate, pODEProblem[15] = recruitment_rate, pODEProblem[16] = infected_moderate_death_rate, pODEProblem[17] = behavioural_change_rate, pODEProblem[18] = infected_moderate_worsen_rate, pODEProblem[19] = treated_weak_improve_rate, pODEProblem[20] = infected_weak_death_rate

	treated_weak = 0.0 
	infected_moderate = 9000.0 
	treated_normal = 0.0 
	infected_normal = 16300.0 
	infected_weak = 11000.0 
	treated_moderate = 0.0 
	susceptible = 1.794e7 

	 return [treated_weak, infected_moderate, treated_normal, infected_normal, infected_weak, treated_moderate, susceptible]
end

function compute_σ(u::AbstractVector, t::Real, θ_sd::AbstractVector, θ_dynamic::AbstractVector, θ_nonDynamic::AbstractVector, 
                   parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :observable_prevalence 
		noiseParameter1_observable_prevalence = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_prevalence 
	end

end