function evalYmod(u, t, dynPar, obsPar, nonDynParam, paramData, observableId, mapObsParam) 

	treated_weak, infected_moderate, treated_normal, infected_normal, infected_weak, treated_moderate, susceptible= u 
	infected_normal_transmission_rate_relative, infected_moderate_transmission_rate, infected_weak_transmission_rate_relative, infected_weak_treatment_rate, infected_normal_worsen_rate, infected_moderate_worsen_rate, treated_moderate_improve_rate, treated_weak_improve_rate, behavioural_change_rate = dynPar 

	prevalence = (1 - susceptible / (susceptible + infected_normal + infected_moderate + infected_weak + treated_normal + treated_moderate + treated_weak)) * 100

	if observableId == "observable_prevalence" 
		return prevalence 
	end

end

function evalU0!(u0Vec, paramVec) 

	c1, treated_normal_death_rate, infected_normal_transmission_rate_relative, infected_normal_treatment_rate, infected_moderate_treatment_rate, infected_weak_transmission_rate_relative, infected_moderate_transmission_rate, infected_normal_worsen_rate, susceptible_death_rate, treated_moderate_improve_rate, infected_weak_treatment_rate, infected_normal_death_rate, treated_weak_death_rate, treated_moderate_death_rate, recruitment_rate, infected_moderate_death_rate, behavioural_change_rate, infected_moderate_worsen_rate, treated_weak_improve_rate, infected_weak_death_rate = paramVec 

	treated_weak = 0.0 
	infected_moderate = 9000.0 
	treated_normal = 0.0 
	infected_normal = 16300.0 
	infected_weak = 11000.0 
	treated_moderate = 0.0 
	susceptible = 1.794e7 

	u0Vec .= treated_weak, infected_moderate, treated_normal, infected_normal, infected_weak, treated_moderate, susceptible
end

function evalU0(paramVec) 

	c1, treated_normal_death_rate, infected_normal_transmission_rate_relative, infected_normal_treatment_rate, infected_moderate_treatment_rate, infected_weak_transmission_rate_relative, infected_moderate_transmission_rate, infected_normal_worsen_rate, susceptible_death_rate, treated_moderate_improve_rate, infected_weak_treatment_rate, infected_normal_death_rate, treated_weak_death_rate, treated_moderate_death_rate, recruitment_rate, infected_moderate_death_rate, behavioural_change_rate, infected_moderate_worsen_rate, treated_weak_improve_rate, infected_weak_death_rate = paramVec 

	treated_weak = 0.0 
	infected_moderate = 9000.0 
	treated_normal = 0.0 
	infected_normal = 16300.0 
	infected_weak = 11000.0 
	treated_moderate = 0.0 
	susceptible = 1.794e7 

	 return [treated_weak, infected_moderate, treated_normal, infected_normal, infected_weak, treated_moderate, susceptible]
end

function evalSd!(u, t, sdPar, dynPar, nonDynParam, paramData, observableId, mapSdParam) 

	treated_weak, infected_moderate, treated_normal, infected_normal, infected_weak, treated_moderate, susceptible= u 
	infected_normal_transmission_rate_relative, infected_moderate_transmission_rate, infected_weak_transmission_rate_relative, infected_weak_treatment_rate, infected_normal_worsen_rate, infected_moderate_worsen_rate, treated_moderate_improve_rate, treated_weak_improve_rate, behavioural_change_rate = dynPar 

	if observableId == "observable_prevalence" 
		noiseParameter1_observable_prevalence = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_prevalence 
	end

end