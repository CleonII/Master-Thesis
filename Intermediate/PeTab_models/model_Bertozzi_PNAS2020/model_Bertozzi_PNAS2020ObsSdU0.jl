function evalYmod(u, t, dynPar, obsPar, paramData, obsData, observableId, simulationId) 

	Infected, Recovered, Susceptible, Ro, dummyVariable= u 
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

	gamma_CA, Trigger_NY, gamma_NY, Io_CA, Lockdown_CA_end, Pop_CA, Trigger_Lockdown, Pop_NY, Lockdown_NY_end, Lockdown_NY_start, Ro_NY, Io_NY, Lockdown_CA_start, Ro_CA = paramVec 

	Infected = ( Io_CA * ( 1 - Trigger_NY ) + Io_NY * Trigger_NY ) / ( Pop_NY * Trigger_NY + Pop_CA * ( 1 - Trigger_NY ) ) 
	Recovered = 0.0 
	Susceptible = ( Pop_NY * Trigger_NY + Pop_CA * ( 1 - Trigger_NY ) + ( - Io_CA * ( 1 - Trigger_NY ) - Io_NY * Trigger_NY ) / ( Pop_NY * Trigger_NY + Pop_CA * ( 1 - Trigger_NY ) ) ) / ( Pop_NY * Trigger_NY + Pop_CA * ( 1 - Trigger_NY ) ) 
	Ro = 2.7 
	dummyVariable = 0.0 

	u0Vec .= Infected, Recovered, Susceptible, Ro, dummyVariable
end

function evalSd!(u, t, sdPar, dynPar, paramData, obsData, observableId, simulationId) 

	Infected, Recovered, Susceptible, Ro, dummyVariable= u 
	gamma_NY, Ro_NY, Io_NY = dynPar 
	gamma_CA_C = paramData.paramVal[2] 
	Ro_CA_C = paramData.paramVal[4] 
	Io_CA_C = paramData.paramVal[6] 
	sd_observable_infected_NY_C = paramData.paramVal[7] 
	sd_observable_infected_CA_C = paramData.paramVal[8] 

	if observableId == "observable_susceptible" 
		noiseParameter1_observable_susceptible = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_susceptible 
	end

	if observableId == "observable_infected" 
		noiseParameter1_observable_infected = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_infected 
	end

	if observableId == "observable_resistant" 
		noiseParameter1_observable_resistant = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_resistant 
	end

end