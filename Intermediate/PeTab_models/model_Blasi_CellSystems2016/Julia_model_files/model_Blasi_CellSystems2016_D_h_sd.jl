#u[1] = x_k5k12k16, u[2] = x_k8, u[3] = x_k16, u[4] = x_0ac, u[5] = x_k12, u[6] = x_k5k8, u[7] = x_k5k12, u[8] = x_k12k16, u[9] = x_k8k12k16, u[10] = x_k5, u[11] = x_k5k16, u[12] = x_k5k8k12, u[13] = x_k8k12, u[14] = x_4ac, u[15] = x_k8k16, u[16] = x_k5k8k16
#pODEProblem[1] = a_k5_k5k12, pODEProblem[2] = a_k8, pODEProblem[3] = d, pODEProblem[4] = default, pODEProblem[5] = a_k12k16_k8k12k16, pODEProblem[6] = a_basal, pODEProblem[7] = a_k5k12_k5k8k12, pODEProblem[8] = a_k8k12k16_4ac, pODEProblem[9] = a_k12_k5k12, pODEProblem[10] = a_k16_k12k16
#
function compute_∂h∂u!(u, t::Real, pODEProblem::AbstractVector, θ_observable::AbstractVector,                    
                       θ_nonDynamic::AbstractVector, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :observable_0ac 
		out[4] = 1
		return nothing
	end

	if observableId == :observable_4ac 
		out[14] = 1
		return nothing
	end

	if observableId == :observable_k12 
		out[5] = 1
		return nothing
	end

	if observableId == :observable_k12k16 
		out[8] = 1
		return nothing
	end

	if observableId == :observable_k16 
		out[3] = 1
		return nothing
	end

	if observableId == :observable_k5 
		out[10] = 1
		return nothing
	end

	if observableId == :observable_k5k12 
		out[7] = 1
		return nothing
	end

	if observableId == :observable_k5k12k16 
		out[1] = 1
		return nothing
	end

	if observableId == :observable_k5k16 
		out[11] = 1
		return nothing
	end

	if observableId == :observable_k5k8 
		out[6] = 1
		return nothing
	end

	if observableId == :observable_k5k8k12 
		out[12] = 1
		return nothing
	end

	if observableId == :observable_k5k8k16 
		out[16] = 1
		return nothing
	end

	if observableId == :observable_k8 
		out[2] = 1
		return nothing
	end

	if observableId == :observable_k8k12 
		out[13] = 1
		return nothing
	end

	if observableId == :observable_k8k12k16 
		out[9] = 1
		return nothing
	end

	if observableId == :observable_k8k16 
		out[15] = 1
		return nothing
	end

end

function compute_∂h∂p!(u, t::Real, pODEProblem::AbstractVector, θ_observable::AbstractVector,
                       θ_nonDynamic::AbstractVector, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :observable_0ac 
		return nothing
	end

	if observableId == :observable_4ac 
		return nothing
	end

	if observableId == :observable_k12 
		return nothing
	end

	if observableId == :observable_k12k16 
		return nothing
	end

	if observableId == :observable_k16 
		return nothing
	end

	if observableId == :observable_k5 
		return nothing
	end

	if observableId == :observable_k5k12 
		return nothing
	end

	if observableId == :observable_k5k12k16 
		return nothing
	end

	if observableId == :observable_k5k16 
		return nothing
	end

	if observableId == :observable_k5k8 
		return nothing
	end

	if observableId == :observable_k5k8k12 
		return nothing
	end

	if observableId == :observable_k5k8k16 
		return nothing
	end

	if observableId == :observable_k8 
		return nothing
	end

	if observableId == :observable_k8k12 
		return nothing
	end

	if observableId == :observable_k8k12k16 
		return nothing
	end

	if observableId == :observable_k8k16 
		return nothing
	end

end

function compute_∂σ∂σu!(u, t::Real, θ_sd::AbstractVector, pODEProblem::AbstractVector, θ_nonDynamic::AbstractVector, 
                        parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :observable_0ac 
		return nothing
	end

	if observableId == :observable_4ac 
		return nothing
	end

	if observableId == :observable_k12 
		return nothing
	end

	if observableId == :observable_k12k16 
		return nothing
	end

	if observableId == :observable_k16 
		return nothing
	end

	if observableId == :observable_k5 
		return nothing
	end

	if observableId == :observable_k5k12 
		return nothing
	end

	if observableId == :observable_k5k12k16 
		return nothing
	end

	if observableId == :observable_k5k16 
		return nothing
	end

	if observableId == :observable_k5k8 
		return nothing
	end

	if observableId == :observable_k5k8k12 
		return nothing
	end

	if observableId == :observable_k5k8k16 
		return nothing
	end

	if observableId == :observable_k8 
		return nothing
	end

	if observableId == :observable_k8k12 
		return nothing
	end

	if observableId == :observable_k8k12k16 
		return nothing
	end

	if observableId == :observable_k8k16 
		return nothing
	end

end

function compute_∂σ∂σp!(u, t::Real, θ_sd::AbstractVector, pODEProblem::AbstractVector, θ_nonDynamic::AbstractVector, 
                        parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :observable_0ac 
		return nothing
	end

	if observableId == :observable_4ac 
		return nothing
	end

	if observableId == :observable_k12 
		return nothing
	end

	if observableId == :observable_k12k16 
		return nothing
	end

	if observableId == :observable_k16 
		return nothing
	end

	if observableId == :observable_k5 
		return nothing
	end

	if observableId == :observable_k5k12 
		return nothing
	end

	if observableId == :observable_k5k12k16 
		return nothing
	end

	if observableId == :observable_k5k16 
		return nothing
	end

	if observableId == :observable_k5k8 
		return nothing
	end

	if observableId == :observable_k5k8k12 
		return nothing
	end

	if observableId == :observable_k5k8k16 
		return nothing
	end

	if observableId == :observable_k8 
		return nothing
	end

	if observableId == :observable_k8k12 
		return nothing
	end

	if observableId == :observable_k8k12k16 
		return nothing
	end

	if observableId == :observable_k8k16 
		return nothing
	end

end

