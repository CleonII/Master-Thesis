function evalDYmodDu(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	treated_weak, infected_moderate, treated_normal, infected_normal, infected_weak, treated_moderate, susceptible= u 
	c1, treated_normal_death_rate, infected_normal_transmission_rate_relative, infected_normal_treatment_rate, infected_moderate_treatment_rate, infected_weak_transmission_rate_relative, infected_moderate_transmission_rate, infected_normal_worsen_rate, susceptible_death_rate, treated_moderate_improve_rate, infected_weak_treatment_rate, infected_normal_death_rate, treated_weak_death_rate, treated_moderate_death_rate, recruitment_rate, infected_moderate_death_rate, behavioural_change_rate, infected_moderate_worsen_rate, treated_weak_improve_rate, infected_weak_death_rate = p 
	prevalence = (1 - susceptible / (susceptible + infected_normal + infected_moderate + infected_weak + treated_normal + treated_moderate + treated_weak)) * 100

	if observableId == "observable_prevalence" 
		return nothing
	end

end

function evalDYmodDp(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	treated_weak, infected_moderate, treated_normal, infected_normal, infected_weak, treated_moderate, susceptible= u 
	c1, treated_normal_death_rate, infected_normal_transmission_rate_relative, infected_normal_treatment_rate, infected_moderate_treatment_rate, infected_weak_transmission_rate_relative, infected_moderate_transmission_rate, infected_normal_worsen_rate, susceptible_death_rate, treated_moderate_improve_rate, infected_weak_treatment_rate, infected_normal_death_rate, treated_weak_death_rate, treated_moderate_death_rate, recruitment_rate, infected_moderate_death_rate, behavioural_change_rate, infected_moderate_worsen_rate, treated_weak_improve_rate, infected_weak_death_rate = p 
	prevalence = (1 - susceptible / (susceptible + infected_normal + infected_moderate + infected_weak + treated_normal + treated_moderate + treated_weak)) * 100

	if observableId == "observable_prevalence" 
		return nothing
	end

end

function evalDSdDu!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	treated_weak, infected_moderate, treated_normal, infected_normal, infected_weak, treated_moderate, susceptible= u 
	c1, treated_normal_death_rate, infected_normal_transmission_rate_relative, infected_normal_treatment_rate, infected_moderate_treatment_rate, infected_weak_transmission_rate_relative, infected_moderate_transmission_rate, infected_normal_worsen_rate, susceptible_death_rate, treated_moderate_improve_rate, infected_weak_treatment_rate, infected_normal_death_rate, treated_weak_death_rate, treated_moderate_death_rate, recruitment_rate, infected_moderate_death_rate, behavioural_change_rate, infected_moderate_worsen_rate, treated_weak_improve_rate, infected_weak_death_rate = p 
	if observableId == "observable_prevalence" 
		return nothing
	end

end

function evalDSdDp!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	treated_weak, infected_moderate, treated_normal, infected_normal, infected_weak, treated_moderate, susceptible= u 
	c1, treated_normal_death_rate, infected_normal_transmission_rate_relative, infected_normal_treatment_rate, infected_moderate_treatment_rate, infected_weak_transmission_rate_relative, infected_moderate_transmission_rate, infected_normal_worsen_rate, susceptible_death_rate, treated_moderate_improve_rate, infected_weak_treatment_rate, infected_normal_death_rate, treated_weak_death_rate, treated_moderate_death_rate, recruitment_rate, infected_moderate_death_rate, behavioural_change_rate, infected_moderate_worsen_rate, treated_weak_improve_rate, infected_weak_death_rate = p 
	if observableId == "observable_prevalence" 
		return nothing
	end

end

