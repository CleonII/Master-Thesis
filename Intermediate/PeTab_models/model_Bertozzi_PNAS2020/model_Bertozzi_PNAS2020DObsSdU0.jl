function evalDYmodDu(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	Infected, Recovered, Susceptible= u 
	Lockdown_NY_end, Pop_CA, Ro_bool1, Ro_bool2, Io_CA, Ro_bool4, Pop_NY, Io_NY, Trigger_NY, USA___CA__NY, Lockdown_CA_start, gamma_NY, Trigger_Lockdown, Ro_NY, Lockdown_CA_end, Ro_bool5, gamma_CA, Ro_CA, Lockdown_NY_start = p 
	if observableId == "observable_susceptible" 
		out[3] = 1
		return nothing
	end

	if observableId == "observable_infected" 
		out[1] = 1
		return nothing
	end

	if observableId == "observable_resistant" 
		out[2] = 1
		return nothing
	end

end

function evalDYmodDp(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	Infected, Recovered, Susceptible= u 
	Lockdown_NY_end, Pop_CA, Ro_bool1, Ro_bool2, Io_CA, Ro_bool4, Pop_NY, Io_NY, Trigger_NY, USA___CA__NY, Lockdown_CA_start, gamma_NY, Trigger_Lockdown, Ro_NY, Lockdown_CA_end, Ro_bool5, gamma_CA, Ro_CA, Lockdown_NY_start = p 
	if observableId == "observable_susceptible" 
		return nothing
	end

	if observableId == "observable_infected" 
		return nothing
	end

	if observableId == "observable_resistant" 
		return nothing
	end

end

function evalDSdDu!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	Infected, Recovered, Susceptible= u 
	Lockdown_NY_end, Pop_CA, Ro_bool1, Ro_bool2, Io_CA, Ro_bool4, Pop_NY, Io_NY, Trigger_NY, USA___CA__NY, Lockdown_CA_start, gamma_NY, Trigger_Lockdown, Ro_NY, Lockdown_CA_end, Ro_bool5, gamma_CA, Ro_CA, Lockdown_NY_start = p 
	if observableId == "observable_susceptible" 
		return nothing
	end

	if observableId == "observable_infected" 
		return nothing
	end

	if observableId == "observable_resistant" 
		return nothing
	end

end

function evalDSdDp!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	Infected, Recovered, Susceptible= u 
	Lockdown_NY_end, Pop_CA, Ro_bool1, Ro_bool2, Io_CA, Ro_bool4, Pop_NY, Io_NY, Trigger_NY, USA___CA__NY, Lockdown_CA_start, gamma_NY, Trigger_Lockdown, Ro_NY, Lockdown_CA_end, Ro_bool5, gamma_CA, Ro_CA, Lockdown_NY_start = p 
	if observableId == "observable_susceptible" 
		return nothing
	end

	if observableId == "observable_infected" 
		return nothing
	end

	if observableId == "observable_resistant" 
		return nothing
	end

end

