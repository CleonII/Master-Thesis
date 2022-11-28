function evalYmod(u, t, dynPar, obsPar, nonDynParam, paramData, observableId, mapObsParam) 

	Infected, Recovered, Susceptible= u 
	gamma_NY, Ro_NY, Io_NY = dynPar 
	gamma_CA_C = paramData.paramVal[2] 
	Ro_CA_C = paramData.paramVal[4] 
	Io_CA_C = paramData.paramVal[6] 
	sd_observable_infected_NY_C = paramData.paramVal[7] 
	sd_observable_infected_CA_C = paramData.paramVal[8] 

	if observableId == "observable_susceptible" 
		return Susceptible 
	end

	if observableId == "observable_infected" 
		return Infected 
	end

	if observableId == "observable_resistant" 
		return Recovered 
	end

end

function evalU0!(u0Vec, paramVec) 

	Lockdown_NY_end, Pop_CA, Io_CA, Pop_NY, Io_NY, Trigger_NY, USA___CA__NY, Lockdown_CA_start, gamma_NY, Trigger_Lockdown, Ro_NY, Lockdown_CA_end, gamma_CA, Ro_CA, Lockdown_NY_start = paramVec 

	Infected = ( Io_CA * ( 1 - Trigger_NY ) + Io_NY * Trigger_NY ) / ( Pop_NY * Trigger_NY + Pop_CA * ( 1 - Trigger_NY ) ) 
	Recovered = 0.0 
	Susceptible = ( Pop_NY * Trigger_NY + Pop_CA * ( 1 - Trigger_NY ) + ( - Io_CA * ( 1 - Trigger_NY ) - Io_NY * Trigger_NY ) / ( Pop_NY * Trigger_NY + Pop_CA * ( 1 - Trigger_NY ) ) ) / ( Pop_NY * Trigger_NY + Pop_CA * ( 1 - Trigger_NY ) ) 

	u0Vec .= Infected, Recovered, Susceptible
end

function evalU0(paramVec) 

	Lockdown_NY_end, Pop_CA, Io_CA, Pop_NY, Io_NY, Trigger_NY, USA___CA__NY, Lockdown_CA_start, gamma_NY, Trigger_Lockdown, Ro_NY, Lockdown_CA_end, gamma_CA, Ro_CA, Lockdown_NY_start = paramVec 

	Infected = ( Io_CA * ( 1 - Trigger_NY ) + Io_NY * Trigger_NY ) / ( Pop_NY * Trigger_NY + Pop_CA * ( 1 - Trigger_NY ) ) 
	Recovered = 0.0 
	Susceptible = ( Pop_NY * Trigger_NY + Pop_CA * ( 1 - Trigger_NY ) + ( - Io_CA * ( 1 - Trigger_NY ) - Io_NY * Trigger_NY ) / ( Pop_NY * Trigger_NY + Pop_CA * ( 1 - Trigger_NY ) ) ) / ( Pop_NY * Trigger_NY + Pop_CA * ( 1 - Trigger_NY ) ) 

	 return [Infected, Recovered, Susceptible]
end

function evalSd!(u, t, sdPar, dynPar, nonDynParam, paramData, observableId, mapSdParam) 

	Infected, Recovered, Susceptible= u 
	gamma_NY, Ro_NY, Io_NY = dynPar 
	gamma_CA_C = paramData.paramVal[2] 
	Ro_CA_C = paramData.paramVal[4] 
	Io_CA_C = paramData.paramVal[6] 
	sd_observable_infected_NY_C = paramData.paramVal[7] 
	sd_observable_infected_CA_C = paramData.paramVal[8] 

	if observableId == "observable_susceptible" 
		noiseParameter1_observable_susceptible = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_susceptible 
	end

	if observableId == "observable_infected" 
		noiseParameter1_observable_infected = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_infected 
	end

	if observableId == "observable_resistant" 
		noiseParameter1_observable_resistant = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_resistant 
	end

end