function evalDYmodDu(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	Quarantined_Infected, Confirmed_Infected, Susceptible, Unquarantined_Infected= u 
	Total_Pop_Hubei, R_Stage_I_Wuhan, gamma_2_Stage_III_Wuhan, Total_Pop_China, Trigger_Wuhan, gamma_1_Stage_I_China, gamma_1_Stage_I_Hubei, sigma, R_Stage_II_Wuhan, R_Stage_II_China, gamma_1_Stage_I_Wuhan, gamma_2_Stage_II_Wuhan, Trigger_China, gamma_2_Stage_II_China, R_Stage_III_Wuhan, Trigger_Stage_I, gamma_2_Stage_I_Wuhan, R_Stage_I_China, Trigger_Stage_III, gamma_1_Stage_II_China, gamma_1_Stage_III_Wuhan, gamma_2_Stage_I_China, Trigger_Stage_II, R_Stage_II_Hubei, Total_Pop_Wuhan, gamma_1_Stage_II_Wuhan, gamma_2_Stage_I_Hubei, R_Stage_I_Hubei, gamma_1_Stage_II_Hubei, Trigger_Hubei, gamma_2_Stage_II_Hubei, China = p 
	Cumulative_Infected = Confirmed_Infected + Quarantined_Infected + Unquarantined_Infected

	if observableId == "observable_confirmed_infected" 
		out[2] = 1
		return nothing
	end

	if observableId == "observable_susceptible" 
		out[3] = 1
		return nothing
	end

	if observableId == "observable_unquarantined_infected" 
		out[4] = 1
		return nothing
	end

	if observableId == "observable_quarantined_infected" 
		out[1] = 1
		return nothing
	end

	if observableId == "observable_cumulative_infected" 
		return nothing
	end

end

function evalDYmodDp(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	Quarantined_Infected, Confirmed_Infected, Susceptible, Unquarantined_Infected= u 
	Total_Pop_Hubei, R_Stage_I_Wuhan, gamma_2_Stage_III_Wuhan, Total_Pop_China, Trigger_Wuhan, gamma_1_Stage_I_China, gamma_1_Stage_I_Hubei, sigma, R_Stage_II_Wuhan, R_Stage_II_China, gamma_1_Stage_I_Wuhan, gamma_2_Stage_II_Wuhan, Trigger_China, gamma_2_Stage_II_China, R_Stage_III_Wuhan, Trigger_Stage_I, gamma_2_Stage_I_Wuhan, R_Stage_I_China, Trigger_Stage_III, gamma_1_Stage_II_China, gamma_1_Stage_III_Wuhan, gamma_2_Stage_I_China, Trigger_Stage_II, R_Stage_II_Hubei, Total_Pop_Wuhan, gamma_1_Stage_II_Wuhan, gamma_2_Stage_I_Hubei, R_Stage_I_Hubei, gamma_1_Stage_II_Hubei, Trigger_Hubei, gamma_2_Stage_II_Hubei, China = p 
	Cumulative_Infected = Confirmed_Infected + Quarantined_Infected + Unquarantined_Infected

	if observableId == "observable_confirmed_infected" 
		return nothing
	end

	if observableId == "observable_susceptible" 
		return nothing
	end

	if observableId == "observable_unquarantined_infected" 
		return nothing
	end

	if observableId == "observable_quarantined_infected" 
		return nothing
	end

	if observableId == "observable_cumulative_infected" 
		return nothing
	end

end

function evalDSdDu!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	Quarantined_Infected, Confirmed_Infected, Susceptible, Unquarantined_Infected= u 
	Total_Pop_Hubei, R_Stage_I_Wuhan, gamma_2_Stage_III_Wuhan, Total_Pop_China, Trigger_Wuhan, gamma_1_Stage_I_China, gamma_1_Stage_I_Hubei, sigma, R_Stage_II_Wuhan, R_Stage_II_China, gamma_1_Stage_I_Wuhan, gamma_2_Stage_II_Wuhan, Trigger_China, gamma_2_Stage_II_China, R_Stage_III_Wuhan, Trigger_Stage_I, gamma_2_Stage_I_Wuhan, R_Stage_I_China, Trigger_Stage_III, gamma_1_Stage_II_China, gamma_1_Stage_III_Wuhan, gamma_2_Stage_I_China, Trigger_Stage_II, R_Stage_II_Hubei, Total_Pop_Wuhan, gamma_1_Stage_II_Wuhan, gamma_2_Stage_I_Hubei, R_Stage_I_Hubei, gamma_1_Stage_II_Hubei, Trigger_Hubei, gamma_2_Stage_II_Hubei, China = p 
	if observableId == "observable_confirmed_infected" 
		return nothing
	end

	if observableId == "observable_susceptible" 
		return nothing
	end

	if observableId == "observable_unquarantined_infected" 
		return nothing
	end

	if observableId == "observable_quarantined_infected" 
		return nothing
	end

	if observableId == "observable_cumulative_infected" 
		return nothing
	end

end

function evalDSdDp!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	Quarantined_Infected, Confirmed_Infected, Susceptible, Unquarantined_Infected= u 
	Total_Pop_Hubei, R_Stage_I_Wuhan, gamma_2_Stage_III_Wuhan, Total_Pop_China, Trigger_Wuhan, gamma_1_Stage_I_China, gamma_1_Stage_I_Hubei, sigma, R_Stage_II_Wuhan, R_Stage_II_China, gamma_1_Stage_I_Wuhan, gamma_2_Stage_II_Wuhan, Trigger_China, gamma_2_Stage_II_China, R_Stage_III_Wuhan, Trigger_Stage_I, gamma_2_Stage_I_Wuhan, R_Stage_I_China, Trigger_Stage_III, gamma_1_Stage_II_China, gamma_1_Stage_III_Wuhan, gamma_2_Stage_I_China, Trigger_Stage_II, R_Stage_II_Hubei, Total_Pop_Wuhan, gamma_1_Stage_II_Wuhan, gamma_2_Stage_I_Hubei, R_Stage_I_Hubei, gamma_1_Stage_II_Hubei, Trigger_Hubei, gamma_2_Stage_II_Hubei, China = p 
	if observableId == "observable_confirmed_infected" 
		return nothing
	end

	if observableId == "observable_susceptible" 
		return nothing
	end

	if observableId == "observable_unquarantined_infected" 
		return nothing
	end

	if observableId == "observable_quarantined_infected" 
		return nothing
	end

	if observableId == "observable_cumulative_infected" 
		return nothing
	end

end

