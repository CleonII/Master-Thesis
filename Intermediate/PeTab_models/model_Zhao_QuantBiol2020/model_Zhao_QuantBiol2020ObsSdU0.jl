function evalYmod(u, t, dynPar, obsPar, paramData, observableId, mapObsParam) 

	Quarantined_Infected, Confirmed_Infected, Susceptible, Unquarantined_Infected, dummyVariable= u 
	R_Stage_I_Wuhan, gamma_1_Stage_I_Wuhan, gamma_2_Stage_I_Wuhan, R_Stage_II_Wuhan, gamma_1_Stage_II_Wuhan, gamma_2_Stage_II_Wuhan, R_Stage_III_Wuhan, gamma_1_Stage_III_Wuhan, gamma_2_Stage_III_Wuhan, R_Stage_I_Hubei, gamma_1_Stage_I_Hubei, gamma_2_Stage_I_Hubei, R_Stage_II_Hubei, gamma_1_Stage_II_Hubei, gamma_2_Stage_II_Hubei, R_Stage_I_China, gamma_1_Stage_I_China, gamma_2_Stage_I_China, R_Stage_II_China, gamma_1_Stage_II_China, gamma_2_Stage_II_China = dynPar 
	sigma_C = paramData.paramVal[1] 

	Cumulative_Infected = Confirmed_Infected+Quarantined_Infected+Unquarantined_Infected

	if observableId == "observable_confirmed_infected" 
		return Confirmed_Infected 
	end

	if observableId == "observable_susceptible" 
		return Susceptible 
	end

	if observableId == "observable_unquarantined_infected" 
		return Unquarantined_Infected 
	end

	if observableId == "observable_quarantined_infected" 
		return Quarantined_Infected 
	end

	if observableId == "observable_cumulative_infected" 
		return  Cumulative_Infected
	end

end

function evalU0!(u0Vec, paramVec) 

	gamma_1_Stage_I_Wuhan, gamma_1_Stage_II_China, gamma_1_Stage_II_Wuhan, gamma_1_Stage_III_Wuhan, gamma_2_Stage_II_China, Trigger_Stage_III, Trigger_Stage_I, Trigger_Hubei, gamma_2_Stage_II_Wuhan, gamma_2_Stage_I_China, Trigger_China, Trigger_Wuhan, gamma_1_Stage_II_Hubei, gamma_2_Stage_I_Wuhan, gamma_2_Stage_III_Wuhan, gamma_2_Stage_I_Hubei, gamma_2_Stage_II_Hubei, sigma, Trigger_Stage_II, gamma_1_Stage_I_China, gamma_1_Stage_I_Hubei, R_Stage_II_Hubei, R_Stage_III_Wuhan, Total_Pop_Hubei, R_Stage_I_Hubei, R_Stage_II_China, R_Stage_I_Wuhan, Total_Pop_China, Total_Pop_Wuhan, R_Stage_I_China, R_Stage_II_Wuhan = paramVec 

	Quarantined_Infected = 2000 * Trigger_China * Trigger_Stage_II + 1500 * Trigger_Hubei * Trigger_Stage_II + 5000 * Trigger_Stage_III * Trigger_Wuhan 
	Confirmed_Infected = Trigger_Wuhan * ( 258 * Trigger_Stage_I + 2000 * Trigger_Stage_II + 36000 * Trigger_Stage_III ) + 4000 * Trigger_China * Trigger_Stage_II + 1600 * Trigger_Hubei * Trigger_Stage_II 
	Susceptible = Total_Pop_China * Trigger_China + Total_Pop_Hubei * Trigger_Hubei + Total_Pop_Wuhan * Trigger_Wuhan 
	Unquarantined_Infected = Trigger_China * ( 291 * Trigger_Stage_I + 2800 * Trigger_Stage_II ) + Trigger_Hubei * ( 270 * Trigger_Stage_I + 5700 * Trigger_Stage_II ) + Trigger_Wuhan * ( 258 * Trigger_Stage_I + 4000 * Trigger_Stage_III + 15270 * Trigger_Stage_II ) 
	dummyVariable = 0.0 

	u0Vec .= Quarantined_Infected, Confirmed_Infected, Susceptible, Unquarantined_Infected, dummyVariable
end

function evalSd!(u, t, sdPar, dynPar, paramData, observableId, mapSdParam) 

	Quarantined_Infected, Confirmed_Infected, Susceptible, Unquarantined_Infected, dummyVariable= u 
	R_Stage_I_Wuhan, gamma_1_Stage_I_Wuhan, gamma_2_Stage_I_Wuhan, R_Stage_II_Wuhan, gamma_1_Stage_II_Wuhan, gamma_2_Stage_II_Wuhan, R_Stage_III_Wuhan, gamma_1_Stage_III_Wuhan, gamma_2_Stage_III_Wuhan, R_Stage_I_Hubei, gamma_1_Stage_I_Hubei, gamma_2_Stage_I_Hubei, R_Stage_II_Hubei, gamma_1_Stage_II_Hubei, gamma_2_Stage_II_Hubei, R_Stage_I_China, gamma_1_Stage_I_China, gamma_2_Stage_I_China, R_Stage_II_China, gamma_1_Stage_II_China, gamma_2_Stage_II_China = dynPar 
	sigma_C = paramData.paramVal[1] 

	if observableId == "observable_confirmed_infected" 
		noiseParameter1_observable_confirmed_infected = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_confirmed_infected 
	end

	if observableId == "observable_susceptible" 
		noiseParameter1_observable_susceptible = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_susceptible 
	end

	if observableId == "observable_unquarantined_infected" 
		noiseParameter1_observable_confirmed_infected = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_confirmed_infected 
	end

	if observableId == "observable_quarantined_infected" 
		noiseParameter1_observable_quarantined_infected = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_quarantined_infected 
	end

	if observableId == "observable_cumulative_infected" 
		noiseParameter1_observable_cumulative_infected = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_cumulative_infected 
	end

end