#u[1] = b10, u[2] = bio, u[3] = ohbio, u[4] = zea, u[5] = bcry, u[6] = ohb10, u[7] = bcar
#θ_dynamicNames[1] = init_b10_1, θ_dynamicNames[2] = init_bcar1, θ_dynamicNames[3] = init_bcar2, θ_dynamicNames[4] = init_bcry_1, θ_dynamicNames[5] = init_ohb10_1, θ_dynamicNames[6] = init_zea_1, θ_dynamicNames[7] = k5, θ_dynamicNames[8] = kb1, θ_dynamicNames[9] = kb2, θ_dynamicNames[10] = kc1, θ_dynamicNames[11] = kc2, θ_dynamicNames[12] = kc4, θ_dynamicNames[13] = szea
#

function compute_h(u::AbstractVector, t::Real, θ_dynamic::AbstractVector, θ_observable::AbstractVector,
                   θ_nonDynamic::AbstractVector, parameterInfo::ParametersInfo, observableId::Symbol, 
                      parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :ob10 
		return u[1] 
	end

	if observableId == :obcar 
		return u[7] 
	end

	if observableId == :obcry 
		return u[5] 
	end

	if observableId == :obio 
		return u[2] 
	end

	if observableId == :oohb10 
		return u[6] 
	end

	if observableId == :ozea 
		return u[4] 
	end

end

function compute_u0!(u0::AbstractVector, pODEProblem::AbstractVector) 

	#pODEProblem[1] = kc2_multiplier, pODEProblem[2] = init_zea, pODEProblem[3] = kc4_multiplier, pODEProblem[4] = cyt, pODEProblem[5] = k5_multiplier, pODEProblem[6] = kc1_multiplier, pODEProblem[7] = init_b10, pODEProblem[8] = init_bcry, pODEProblem[9] = kb1_multiplier, pODEProblem[10] = kb2_multiplier, pODEProblem[11] = kc1, pODEProblem[12] = kc4, pODEProblem[13] = init_ohb10, pODEProblem[14] = init_bcar, pODEProblem[15] = kc2, pODEProblem[16] = kb2, pODEProblem[17] = k5, pODEProblem[18] = kb1

	b10 = pODEProblem[7] 
	bio = 0.0 
	ohbio = 0.0 
	zea = pODEProblem[2] 
	bcry = pODEProblem[8] 
	ohb10 = pODEProblem[13] 
	bcar = pODEProblem[14] 

	u0 .= b10, bio, ohbio, zea, bcry, ohb10, bcar
end

function compute_u0(pODEProblem::AbstractVector)::AbstractVector 

	#pODEProblem[1] = kc2_multiplier, pODEProblem[2] = init_zea, pODEProblem[3] = kc4_multiplier, pODEProblem[4] = cyt, pODEProblem[5] = k5_multiplier, pODEProblem[6] = kc1_multiplier, pODEProblem[7] = init_b10, pODEProblem[8] = init_bcry, pODEProblem[9] = kb1_multiplier, pODEProblem[10] = kb2_multiplier, pODEProblem[11] = kc1, pODEProblem[12] = kc4, pODEProblem[13] = init_ohb10, pODEProblem[14] = init_bcar, pODEProblem[15] = kc2, pODEProblem[16] = kb2, pODEProblem[17] = k5, pODEProblem[18] = kb1

	b10 = pODEProblem[7] 
	bio = 0.0 
	ohbio = 0.0 
	zea = pODEProblem[2] 
	bcry = pODEProblem[8] 
	ohb10 = pODEProblem[13] 
	bcar = pODEProblem[14] 

	 return [b10, bio, ohbio, zea, bcry, ohb10, bcar]
end

function compute_σ(u::AbstractVector, t::Real, θ_sd::AbstractVector, θ_dynamic::AbstractVector, θ_nonDynamic::AbstractVector, 
                   parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :ob10 
		noiseParameter1_ob10 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_ob10 
	end

	if observableId == :obcar 
		noiseParameter1_obcar = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_obcar 
	end

	if observableId == :obcry 
		noiseParameter1_obcry = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_obcry 
	end

	if observableId == :obio 
		noiseParameter1_obio = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_obio 
	end

	if observableId == :oohb10 
		noiseParameter1_oohb10 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_oohb10 
	end

	if observableId == :ozea 
		noiseParameter1_ozea = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_ozea 
	end

end