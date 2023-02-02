#u[1] = Quarantined_Infected, u[2] = Confirmed_Infected, u[3] = Susceptible, u[4] = Unquarantined_Infected
#θ_dynamicNames[1] = R_Stage_I_Wuhan, θ_dynamicNames[2] = gamma_1_Stage_I_Wuhan, θ_dynamicNames[3] = gamma_2_Stage_I_Wuhan, θ_dynamicNames[4] = R_Stage_II_Wuhan, θ_dynamicNames[5] = gamma_1_Stage_II_Wuhan, θ_dynamicNames[6] = gamma_2_Stage_II_Wuhan, θ_dynamicNames[7] = R_Stage_III_Wuhan, θ_dynamicNames[8] = gamma_1_Stage_III_Wuhan, θ_dynamicNames[9] = gamma_2_Stage_III_Wuhan, θ_dynamicNames[10] = R_Stage_I_Hubei, θ_dynamicNames[11] = gamma_1_Stage_I_Hubei, θ_dynamicNames[12] = gamma_2_Stage_I_Hubei, θ_dynamicNames[13] = R_Stage_II_Hubei, θ_dynamicNames[14] = gamma_1_Stage_II_Hubei, θ_dynamicNames[15] = gamma_2_Stage_II_Hubei, θ_dynamicNames[16] = R_Stage_I_China, θ_dynamicNames[17] = gamma_1_Stage_I_China, θ_dynamicNames[18] = gamma_2_Stage_I_China, θ_dynamicNames[19] = R_Stage_II_China, θ_dynamicNames[20] = gamma_1_Stage_II_China, θ_dynamicNames[21] = gamma_2_Stage_II_China
##parameterInfo.nominalValue[1] = sigma_C 


function compute_h(u::AbstractVector, t::Real, θ_dynamic::AbstractVector, θ_observable::AbstractVector,
                   θ_nonDynamic::AbstractVector, parameterInfo::ParametersInfo, observableId::Symbol, 
                      parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :observable_confirmed_infected 
		return u[2] 
	end

	if observableId == :observable_susceptible 
		return u[3] 
	end

	if observableId == :observable_unquarantined_infected 
		return u[4] 
	end

	if observableId == :observable_quarantined_infected 
		return u[1] 
	end

	if observableId == :observable_cumulative_infected 
		return ( u[2] + u[1] + u[4] ) 
	end

end

function compute_u0!(u0::AbstractVector, pODEProblem::AbstractVector) 

	#pODEProblem[1] = Total_Pop_Hubei, pODEProblem[2] = R_Stage_I_Wuhan, pODEProblem[3] = gamma_2_Stage_III_Wuhan, pODEProblem[4] = Total_Pop_China, pODEProblem[5] = Trigger_Wuhan, pODEProblem[6] = gamma_1_Stage_I_China, pODEProblem[7] = gamma_1_Stage_I_Hubei, pODEProblem[8] = sigma, pODEProblem[9] = R_Stage_II_Wuhan, pODEProblem[10] = R_Stage_II_China, pODEProblem[11] = gamma_1_Stage_I_Wuhan, pODEProblem[12] = gamma_2_Stage_II_Wuhan, pODEProblem[13] = Trigger_China, pODEProblem[14] = gamma_2_Stage_II_China, pODEProblem[15] = R_Stage_III_Wuhan, pODEProblem[16] = Trigger_Stage_I, pODEProblem[17] = gamma_2_Stage_I_Wuhan, pODEProblem[18] = R_Stage_I_China, pODEProblem[19] = Trigger_Stage_III, pODEProblem[20] = gamma_1_Stage_II_China, pODEProblem[21] = gamma_1_Stage_III_Wuhan, pODEProblem[22] = gamma_2_Stage_I_China, pODEProblem[23] = Trigger_Stage_II, pODEProblem[24] = R_Stage_II_Hubei, pODEProblem[25] = Total_Pop_Wuhan, pODEProblem[26] = gamma_1_Stage_II_Wuhan, pODEProblem[27] = gamma_2_Stage_I_Hubei, pODEProblem[28] = R_Stage_I_Hubei, pODEProblem[29] = gamma_1_Stage_II_Hubei, pODEProblem[30] = Trigger_Hubei, pODEProblem[31] = gamma_2_Stage_II_Hubei, pODEProblem[32] = China

	Quarantined_Infected = 2000 * pODEProblem[13] * pODEProblem[23] + 1500 * pODEProblem[30] * pODEProblem[23] + 5000 * pODEProblem[19] * pODEProblem[5] 
	Confirmed_Infected = pODEProblem[5] * ( 258 * pODEProblem[16] + 2000 * pODEProblem[23] + 36000 * pODEProblem[19] ) + 4000 * pODEProblem[13] * pODEProblem[23] + 1600 * pODEProblem[30] * pODEProblem[23] 
	Susceptible = pODEProblem[4] * pODEProblem[13] + pODEProblem[1] * pODEProblem[30] + pODEProblem[25] * pODEProblem[5] 
	Unquarantined_Infected = pODEProblem[13] * ( 291 * pODEProblem[16] + 2800 * pODEProblem[23] ) + pODEProblem[30] * ( 270 * pODEProblem[16] + 5700 * pODEProblem[23] ) + pODEProblem[5] * ( 258 * pODEProblem[16] + 4000 * pODEProblem[19] + 15270 * pODEProblem[23] ) 

	u0 .= Quarantined_Infected, Confirmed_Infected, Susceptible, Unquarantined_Infected
end

function compute_u0(pODEProblem::AbstractVector)::AbstractVector 

	#pODEProblem[1] = Total_Pop_Hubei, pODEProblem[2] = R_Stage_I_Wuhan, pODEProblem[3] = gamma_2_Stage_III_Wuhan, pODEProblem[4] = Total_Pop_China, pODEProblem[5] = Trigger_Wuhan, pODEProblem[6] = gamma_1_Stage_I_China, pODEProblem[7] = gamma_1_Stage_I_Hubei, pODEProblem[8] = sigma, pODEProblem[9] = R_Stage_II_Wuhan, pODEProblem[10] = R_Stage_II_China, pODEProblem[11] = gamma_1_Stage_I_Wuhan, pODEProblem[12] = gamma_2_Stage_II_Wuhan, pODEProblem[13] = Trigger_China, pODEProblem[14] = gamma_2_Stage_II_China, pODEProblem[15] = R_Stage_III_Wuhan, pODEProblem[16] = Trigger_Stage_I, pODEProblem[17] = gamma_2_Stage_I_Wuhan, pODEProblem[18] = R_Stage_I_China, pODEProblem[19] = Trigger_Stage_III, pODEProblem[20] = gamma_1_Stage_II_China, pODEProblem[21] = gamma_1_Stage_III_Wuhan, pODEProblem[22] = gamma_2_Stage_I_China, pODEProblem[23] = Trigger_Stage_II, pODEProblem[24] = R_Stage_II_Hubei, pODEProblem[25] = Total_Pop_Wuhan, pODEProblem[26] = gamma_1_Stage_II_Wuhan, pODEProblem[27] = gamma_2_Stage_I_Hubei, pODEProblem[28] = R_Stage_I_Hubei, pODEProblem[29] = gamma_1_Stage_II_Hubei, pODEProblem[30] = Trigger_Hubei, pODEProblem[31] = gamma_2_Stage_II_Hubei, pODEProblem[32] = China

	Quarantined_Infected = 2000 * pODEProblem[13] * pODEProblem[23] + 1500 * pODEProblem[30] * pODEProblem[23] + 5000 * pODEProblem[19] * pODEProblem[5] 
	Confirmed_Infected = pODEProblem[5] * ( 258 * pODEProblem[16] + 2000 * pODEProblem[23] + 36000 * pODEProblem[19] ) + 4000 * pODEProblem[13] * pODEProblem[23] + 1600 * pODEProblem[30] * pODEProblem[23] 
	Susceptible = pODEProblem[4] * pODEProblem[13] + pODEProblem[1] * pODEProblem[30] + pODEProblem[25] * pODEProblem[5] 
	Unquarantined_Infected = pODEProblem[13] * ( 291 * pODEProblem[16] + 2800 * pODEProblem[23] ) + pODEProblem[30] * ( 270 * pODEProblem[16] + 5700 * pODEProblem[23] ) + pODEProblem[5] * ( 258 * pODEProblem[16] + 4000 * pODEProblem[19] + 15270 * pODEProblem[23] ) 

	 return [Quarantined_Infected, Confirmed_Infected, Susceptible, Unquarantined_Infected]
end

function compute_σ(u::AbstractVector, t::Real, θ_sd::AbstractVector, θ_dynamic::AbstractVector, θ_nonDynamic::AbstractVector, 
                   parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :observable_confirmed_infected 
		noiseParameter1_observable_confirmed_infected = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_confirmed_infected 
	end

	if observableId == :observable_susceptible 
		noiseParameter1_observable_susceptible = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_susceptible 
	end

	if observableId == :observable_unquarantined_infected 
		noiseParameter1_observable_confirmed_infected = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_confirmed_infected 
	end

	if observableId == :observable_quarantined_infected 
		noiseParameter1_observable_quarantined_infected = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_quarantined_infected 
	end

	if observableId == :observable_cumulative_infected 
		noiseParameter1_observable_cumulative_infected = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_cumulative_infected 
	end

end