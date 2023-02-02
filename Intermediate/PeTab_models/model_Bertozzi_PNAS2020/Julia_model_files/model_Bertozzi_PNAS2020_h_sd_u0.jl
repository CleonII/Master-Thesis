#u[1] = Infected, u[2] = Recovered, u[3] = Susceptible
#θ_dynamicNames[1] = gamma_NY, θ_dynamicNames[2] = Ro_NY, θ_dynamicNames[3] = Io_NY
##parameterInfo.nominalValue[2] = gamma_CA_C 
#parameterInfo.nominalValue[4] = Ro_CA_C 
#parameterInfo.nominalValue[6] = Io_CA_C 
#parameterInfo.nominalValue[7] = sd_observable_infected_NY_C 
#parameterInfo.nominalValue[8] = sd_observable_infected_CA_C 


function compute_h(u::AbstractVector, t::Real, θ_dynamic::AbstractVector, θ_observable::AbstractVector,
                   θ_nonDynamic::AbstractVector, parameterInfo::ParametersInfo, observableId::Symbol, 
                      parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :observable_susceptible 
		return u[3] 
	end

	if observableId == :observable_infected 
		return u[1] 
	end

	if observableId == :observable_resistant 
		return u[2] 
	end

end

function compute_u0!(u0::AbstractVector, pODEProblem::AbstractVector) 

	#pODEProblem[1] = Lockdown_NY_end, pODEProblem[2] = Pop_CA, pODEProblem[3] = Ro_bool1, pODEProblem[4] = Ro_bool2, pODEProblem[5] = Io_CA, pODEProblem[6] = Ro_bool3, pODEProblem[7] = Ro_bool4, pODEProblem[8] = Pop_NY, pODEProblem[9] = Io_NY, pODEProblem[10] = Trigger_NY, pODEProblem[11] = USA___CA__NY, pODEProblem[12] = Lockdown_CA_start, pODEProblem[13] = gamma_NY, pODEProblem[14] = Trigger_Lockdown, pODEProblem[15] = Ro_NY, pODEProblem[16] = Lockdown_CA_end, pODEProblem[17] = gamma_CA, pODEProblem[18] = Ro_CA, pODEProblem[19] = Lockdown_NY_start

	Infected = ( pODEProblem[5] * ( 1 - pODEProblem[10] ) + pODEProblem[9] * pODEProblem[10] ) / ( pODEProblem[8] * pODEProblem[10] + pODEProblem[2] * ( 1 - pODEProblem[10] ) ) 
	Recovered = 0.0 
	Susceptible = ( pODEProblem[8] * pODEProblem[10] + pODEProblem[2] * ( 1 - pODEProblem[10] ) + ( - pODEProblem[5] * ( 1 - pODEProblem[10] ) - pODEProblem[9] * pODEProblem[10] ) / ( pODEProblem[8] * pODEProblem[10] + pODEProblem[2] * ( 1 - pODEProblem[10] ) ) ) / ( pODEProblem[8] * pODEProblem[10] + pODEProblem[2] * ( 1 - pODEProblem[10] ) ) 

	u0 .= Infected, Recovered, Susceptible
end

function compute_u0(pODEProblem::AbstractVector)::AbstractVector 

	#pODEProblem[1] = Lockdown_NY_end, pODEProblem[2] = Pop_CA, pODEProblem[3] = Ro_bool1, pODEProblem[4] = Ro_bool2, pODEProblem[5] = Io_CA, pODEProblem[6] = Ro_bool3, pODEProblem[7] = Ro_bool4, pODEProblem[8] = Pop_NY, pODEProblem[9] = Io_NY, pODEProblem[10] = Trigger_NY, pODEProblem[11] = USA___CA__NY, pODEProblem[12] = Lockdown_CA_start, pODEProblem[13] = gamma_NY, pODEProblem[14] = Trigger_Lockdown, pODEProblem[15] = Ro_NY, pODEProblem[16] = Lockdown_CA_end, pODEProblem[17] = gamma_CA, pODEProblem[18] = Ro_CA, pODEProblem[19] = Lockdown_NY_start

	Infected = ( pODEProblem[5] * ( 1 - pODEProblem[10] ) + pODEProblem[9] * pODEProblem[10] ) / ( pODEProblem[8] * pODEProblem[10] + pODEProblem[2] * ( 1 - pODEProblem[10] ) ) 
	Recovered = 0.0 
	Susceptible = ( pODEProblem[8] * pODEProblem[10] + pODEProblem[2] * ( 1 - pODEProblem[10] ) + ( - pODEProblem[5] * ( 1 - pODEProblem[10] ) - pODEProblem[9] * pODEProblem[10] ) / ( pODEProblem[8] * pODEProblem[10] + pODEProblem[2] * ( 1 - pODEProblem[10] ) ) ) / ( pODEProblem[8] * pODEProblem[10] + pODEProblem[2] * ( 1 - pODEProblem[10] ) ) 

	 return [Infected, Recovered, Susceptible]
end

function compute_σ(u::AbstractVector, t::Real, θ_sd::AbstractVector, θ_dynamic::AbstractVector, θ_nonDynamic::AbstractVector, 
                   parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :observable_susceptible 
		noiseParameter1_observable_susceptible = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_susceptible 
	end

	if observableId == :observable_infected 
		noiseParameter1_observable_infected = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_infected 
	end

	if observableId == :observable_resistant 
		noiseParameter1_observable_resistant = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_resistant 
	end

end