#u[1] = Quarantined_Infected, u[2] = Confirmed_Infected, u[3] = Susceptible, u[4] = Unquarantined_Infected
#pODEProblem[1] = Total_Pop_Hubei, pODEProblem[2] = R_Stage_I_Wuhan, pODEProblem[3] = gamma_2_Stage_III_Wuhan, pODEProblem[4] = Total_Pop_China, pODEProblem[5] = Trigger_Wuhan, pODEProblem[6] = gamma_1_Stage_I_China, pODEProblem[7] = gamma_1_Stage_I_Hubei, pODEProblem[8] = sigma, pODEProblem[9] = R_Stage_II_Wuhan, pODEProblem[10] = R_Stage_II_China, pODEProblem[11] = gamma_1_Stage_I_Wuhan, pODEProblem[12] = gamma_2_Stage_II_Wuhan, pODEProblem[13] = Trigger_China, pODEProblem[14] = gamma_2_Stage_II_China, pODEProblem[15] = R_Stage_III_Wuhan, pODEProblem[16] = Trigger_Stage_I, pODEProblem[17] = gamma_2_Stage_I_Wuhan, pODEProblem[18] = R_Stage_I_China, pODEProblem[19] = Trigger_Stage_III, pODEProblem[20] = gamma_1_Stage_II_China, pODEProblem[21] = gamma_1_Stage_III_Wuhan, pODEProblem[22] = gamma_2_Stage_I_China, pODEProblem[23] = Trigger_Stage_II, pODEProblem[24] = R_Stage_II_Hubei, pODEProblem[25] = Total_Pop_Wuhan, pODEProblem[26] = gamma_1_Stage_II_Wuhan, pODEProblem[27] = gamma_2_Stage_I_Hubei, pODEProblem[28] = R_Stage_I_Hubei, pODEProblem[29] = gamma_1_Stage_II_Hubei, pODEProblem[30] = Trigger_Hubei, pODEProblem[31] = gamma_2_Stage_II_Hubei, pODEProblem[32] = China
#
function compute_∂h∂u!(u, t::Real, pODEProblem::AbstractVector, θ_observable::AbstractVector,                    
                       θ_nonDynamic::AbstractVector, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :observable_confirmed_infected 
		out[2] = 1
		return nothing
	end

	if observableId == :observable_susceptible 
		out[3] = 1
		return nothing
	end

	if observableId == :observable_unquarantined_infected 
		out[4] = 1
		return nothing
	end

	if observableId == :observable_quarantined_infected 
		out[1] = 1
		return nothing
	end

	if observableId == :observable_cumulative_infected 
		out[1] = 1
		out[2] = 1
		out[4] = 1
		return nothing
	end

end

function compute_∂h∂p!(u, t::Real, pODEProblem::AbstractVector, θ_observable::AbstractVector,
                       θ_nonDynamic::AbstractVector, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :observable_confirmed_infected 
		return nothing
	end

	if observableId == :observable_susceptible 
		return nothing
	end

	if observableId == :observable_unquarantined_infected 
		return nothing
	end

	if observableId == :observable_quarantined_infected 
		return nothing
	end

	if observableId == :observable_cumulative_infected 
		return nothing
	end

end

function compute_∂σ∂σu!(u, t::Real, θ_sd::AbstractVector, pODEProblem::AbstractVector, θ_nonDynamic::AbstractVector, 
                        parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :observable_confirmed_infected 
		return nothing
	end

	if observableId == :observable_susceptible 
		return nothing
	end

	if observableId == :observable_unquarantined_infected 
		return nothing
	end

	if observableId == :observable_quarantined_infected 
		return nothing
	end

	if observableId == :observable_cumulative_infected 
		return nothing
	end

end

function compute_∂σ∂σp!(u, t::Real, θ_sd::AbstractVector, pODEProblem::AbstractVector, θ_nonDynamic::AbstractVector, 
                        parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :observable_confirmed_infected 
		return nothing
	end

	if observableId == :observable_susceptible 
		return nothing
	end

	if observableId == :observable_unquarantined_infected 
		return nothing
	end

	if observableId == :observable_quarantined_infected 
		return nothing
	end

	if observableId == :observable_cumulative_infected 
		return nothing
	end

end

