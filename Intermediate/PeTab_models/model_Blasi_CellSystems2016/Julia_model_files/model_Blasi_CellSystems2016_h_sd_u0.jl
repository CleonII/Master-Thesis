#u[1] = x_k5k12k16, u[2] = x_k8, u[3] = x_k16, u[4] = x_0ac, u[5] = x_k12, u[6] = x_k5k8, u[7] = x_k5k12, u[8] = x_k12k16, u[9] = x_k8k12k16, u[10] = x_k5, u[11] = x_k5k16, u[12] = x_k5k8k12, u[13] = x_k8k12, u[14] = x_4ac, u[15] = x_k8k16, u[16] = x_k5k8k16
#θ_dynamicNames[1] = a_basal, θ_dynamicNames[2] = a_k8, θ_dynamicNames[3] = a_k5_k5k12, θ_dynamicNames[4] = a_k12_k5k12, θ_dynamicNames[5] = a_k16_k12k16, θ_dynamicNames[6] = a_k5k12_k5k8k12, θ_dynamicNames[7] = a_k12k16_k8k12k16, θ_dynamicNames[8] = a_k8k12k16_4ac
##parameterInfo.nominalValue[9] = d_C 


function compute_h(u::AbstractVector, t::Real, θ_dynamic::AbstractVector, θ_observable::AbstractVector,
                   θ_nonDynamic::AbstractVector, parameterInfo::ParametersInfo, observableId::Symbol, 
                      parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :observable_0ac 
		return u[4] 
	end

	if observableId == :observable_4ac 
		return u[14] 
	end

	if observableId == :observable_k12 
		return u[5] 
	end

	if observableId == :observable_k12k16 
		return u[8] 
	end

	if observableId == :observable_k16 
		return u[3] 
	end

	if observableId == :observable_k5 
		return u[10] 
	end

	if observableId == :observable_k5k12 
		return u[7] 
	end

	if observableId == :observable_k5k12k16 
		return u[1] 
	end

	if observableId == :observable_k5k16 
		return u[11] 
	end

	if observableId == :observable_k5k8 
		return u[6] 
	end

	if observableId == :observable_k5k8k12 
		return u[12] 
	end

	if observableId == :observable_k5k8k16 
		return u[16] 
	end

	if observableId == :observable_k8 
		return u[2] 
	end

	if observableId == :observable_k8k12 
		return u[13] 
	end

	if observableId == :observable_k8k12k16 
		return u[9] 
	end

	if observableId == :observable_k8k16 
		return u[15] 
	end

end

function compute_u0!(u0::AbstractVector, pODEProblem::AbstractVector) 

	#pODEProblem[1] = a_k5_k5k12, pODEProblem[2] = a_k8, pODEProblem[3] = d, pODEProblem[4] = default, pODEProblem[5] = a_k12k16_k8k12k16, pODEProblem[6] = a_basal, pODEProblem[7] = a_k5k12_k5k8k12, pODEProblem[8] = a_k8k12k16_4ac, pODEProblem[9] = a_k12_k5k12, pODEProblem[10] = a_k16_k12k16

	x_k5k12k16 = 0.0 
	x_k8 = 0.0 
	x_k16 = 0.0 
	x_0ac = 1.0 
	x_k12 = 0.0 
	x_k5k8 = 0.0 
	x_k5k12 = 0.0 
	x_k12k16 = 0.0 
	x_k8k12k16 = 0.0 
	x_k5 = 0.0 
	x_k5k16 = 0.0 
	x_k5k8k12 = 0.0 
	x_k8k12 = 0.0 
	x_4ac = 0.0 
	x_k8k16 = 0.0 
	x_k5k8k16 = 0.0 

	u0 .= x_k5k12k16, x_k8, x_k16, x_0ac, x_k12, x_k5k8, x_k5k12, x_k12k16, x_k8k12k16, x_k5, x_k5k16, x_k5k8k12, x_k8k12, x_4ac, x_k8k16, x_k5k8k16
end

function compute_u0(pODEProblem::AbstractVector)::AbstractVector 

	#pODEProblem[1] = a_k5_k5k12, pODEProblem[2] = a_k8, pODEProblem[3] = d, pODEProblem[4] = default, pODEProblem[5] = a_k12k16_k8k12k16, pODEProblem[6] = a_basal, pODEProblem[7] = a_k5k12_k5k8k12, pODEProblem[8] = a_k8k12k16_4ac, pODEProblem[9] = a_k12_k5k12, pODEProblem[10] = a_k16_k12k16

	x_k5k12k16 = 0.0 
	x_k8 = 0.0 
	x_k16 = 0.0 
	x_0ac = 1.0 
	x_k12 = 0.0 
	x_k5k8 = 0.0 
	x_k5k12 = 0.0 
	x_k12k16 = 0.0 
	x_k8k12k16 = 0.0 
	x_k5 = 0.0 
	x_k5k16 = 0.0 
	x_k5k8k12 = 0.0 
	x_k8k12 = 0.0 
	x_4ac = 0.0 
	x_k8k16 = 0.0 
	x_k5k8k16 = 0.0 

	 return [x_k5k12k16, x_k8, x_k16, x_0ac, x_k12, x_k5k8, x_k5k12, x_k12k16, x_k8k12k16, x_k5, x_k5k16, x_k5k8k12, x_k8k12, x_4ac, x_k8k16, x_k5k8k16]
end

function compute_σ(u::AbstractVector, t::Real, θ_sd::AbstractVector, θ_dynamic::AbstractVector, θ_nonDynamic::AbstractVector, 
                   parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :observable_0ac 
		noiseParameter1_observable_0ac = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_0ac 
	end

	if observableId == :observable_4ac 
		noiseParameter1_observable_4ac = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_4ac 
	end

	if observableId == :observable_k12 
		noiseParameter1_observable_k12 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_k12 
	end

	if observableId == :observable_k12k16 
		noiseParameter1_observable_k12k16 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_k12k16 
	end

	if observableId == :observable_k16 
		noiseParameter1_observable_k16 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_k16 
	end

	if observableId == :observable_k5 
		noiseParameter1_observable_k5 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_k5 
	end

	if observableId == :observable_k5k12 
		noiseParameter1_observable_k5k12 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_k5k12 
	end

	if observableId == :observable_k5k12k16 
		noiseParameter1_observable_k5k12k16 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_k5k12k16 
	end

	if observableId == :observable_k5k16 
		noiseParameter1_observable_k5k16 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_k5k16 
	end

	if observableId == :observable_k5k8 
		noiseParameter1_observable_k5k8 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_k5k8 
	end

	if observableId == :observable_k5k8k12 
		noiseParameter1_observable_k5k8k12 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_k5k8k12 
	end

	if observableId == :observable_k5k8k16 
		noiseParameter1_observable_k5k8k16 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_k5k8k16 
	end

	if observableId == :observable_k8 
		noiseParameter1_observable_k8 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_k8 
	end

	if observableId == :observable_k8k12 
		noiseParameter1_observable_k8k12 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_k8k12 
	end

	if observableId == :observable_k8k12k16 
		noiseParameter1_observable_k8k12k16 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_k8k12k16 
	end

	if observableId == :observable_k8k16 
		noiseParameter1_observable_k8k16 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_k8k16 
	end

end