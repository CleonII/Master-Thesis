#u[1] = Infected, u[2] = Recovered, u[3] = Susceptible
#pODEProblem[1] = Lockdown_NY_end, pODEProblem[2] = Pop_CA, pODEProblem[3] = Ro_bool1, pODEProblem[4] = Ro_bool2, pODEProblem[5] = Io_CA, pODEProblem[6] = Ro_bool3, pODEProblem[7] = Ro_bool4, pODEProblem[8] = Pop_NY, pODEProblem[9] = Io_NY, pODEProblem[10] = Trigger_NY, pODEProblem[11] = USA___CA__NY, pODEProblem[12] = Lockdown_CA_start, pODEProblem[13] = gamma_NY, pODEProblem[14] = Trigger_Lockdown, pODEProblem[15] = Ro_NY, pODEProblem[16] = Lockdown_CA_end, pODEProblem[17] = gamma_CA, pODEProblem[18] = Ro_CA, pODEProblem[19] = Lockdown_NY_start
#
function compute_∂h∂u!(u, t::Real, pODEProblem::AbstractVector, θ_observable::AbstractVector,                    
                       θ_nonDynamic::AbstractVector, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :observable_susceptible 
		out[3] = 1
		return nothing
	end

	if observableId == :observable_infected 
		out[1] = 1
		return nothing
	end

	if observableId == :observable_resistant 
		out[2] = 1
		return nothing
	end

end

function compute_∂h∂p!(u, t::Real, pODEProblem::AbstractVector, θ_observable::AbstractVector,
                       θ_nonDynamic::AbstractVector, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :observable_susceptible 
		return nothing
	end

	if observableId == :observable_infected 
		return nothing
	end

	if observableId == :observable_resistant 
		return nothing
	end

end

function compute_∂σ∂σu!(u, t::Real, θ_sd::AbstractVector, pODEProblem::AbstractVector, θ_nonDynamic::AbstractVector, 
                        parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :observable_susceptible 
		return nothing
	end

	if observableId == :observable_infected 
		return nothing
	end

	if observableId == :observable_resistant 
		return nothing
	end

end

function compute_∂σ∂σp!(u, t::Real, θ_sd::AbstractVector, pODEProblem::AbstractVector, θ_nonDynamic::AbstractVector, 
                        parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :observable_susceptible 
		return nothing
	end

	if observableId == :observable_infected 
		return nothing
	end

	if observableId == :observable_resistant 
		return nothing
	end

end

