#u[1] = X_protein, u[2] = GFP_mRNA, u[3] = Y_mRNA, u[4] = X_mRNA, u[5] = Z_mRNA, u[6] = Z_protein, u[7] = Y_protein, u[8] = GFP
#θ_dynamicNames[1] = KM, θ_dynamicNames[2] = eff, θ_dynamicNames[3] = eff_GFP, θ_dynamicNames[4] = init_GFP, θ_dynamicNames[5] = init_GFP_mRNA, θ_dynamicNames[6] = init_X_mRNA, θ_dynamicNames[7] = init_X_protein, θ_dynamicNames[8] = init_Y_mRNA, θ_dynamicNames[9] = init_Y_protein, θ_dynamicNames[10] = init_Z_mRNA, θ_dynamicNames[11] = init_Z_protein, θ_dynamicNames[12] = n_Hill, θ_dynamicNames[13] = tau_mRNA, θ_dynamicNames[14] = tau_mRNA_GFP, θ_dynamicNames[15] = tau_prot, θ_dynamicNames[16] = tau_prot_GFP, θ_dynamicNames[17] = tps_active, θ_dynamicNames[18] = tps_repr
#

function compute_h(u::AbstractVector, t::Real, θ_dynamic::AbstractVector, θ_observable::AbstractVector,
                   θ_nonDynamic::AbstractVector, parameterInfo::ParametersInfo, observableId::Symbol, 
                      parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :fluorescence 
		observableParameter1_fluorescence, observableParameter2_fluorescence = getObsOrSdParam(θ_observable, parameterMap)
		return observableParameter1_fluorescence + u[8] * observableParameter2_fluorescence 
	end

end

function compute_u0!(u0::AbstractVector, pODEProblem::AbstractVector) 

	#pODEProblem[1] = tau_mRNA, pODEProblem[2] = tps_repr, pODEProblem[3] = init_GFP, pODEProblem[4] = n_Hill, pODEProblem[5] = init_Y_mRNA, pODEProblem[6] = init_Z_mRNA, pODEProblem[7] = init_X_protein, pODEProblem[8] = init_Y_protein, pODEProblem[9] = tau_mRNA_GFP, pODEProblem[10] = eff, pODEProblem[11] = init_Z_protein, pODEProblem[12] = tau_prot_GFP, pODEProblem[13] = eff_GFP, pODEProblem[14] = tps_active, pODEProblem[15] = init_GFP_mRNA, pODEProblem[16] = KM, pODEProblem[17] = init_X_mRNA, pODEProblem[18] = tau_prot, pODEProblem[19] = cell

	X_protein = pODEProblem[7] 
	GFP_mRNA = pODEProblem[15] 
	Y_mRNA = pODEProblem[5] 
	X_mRNA = pODEProblem[17] 
	Z_mRNA = pODEProblem[6] 
	Z_protein = pODEProblem[11] 
	Y_protein = pODEProblem[8] 
	GFP = pODEProblem[3] 

	u0 .= X_protein, GFP_mRNA, Y_mRNA, X_mRNA, Z_mRNA, Z_protein, Y_protein, GFP
end

function compute_u0(pODEProblem::AbstractVector)::AbstractVector 

	#pODEProblem[1] = tau_mRNA, pODEProblem[2] = tps_repr, pODEProblem[3] = init_GFP, pODEProblem[4] = n_Hill, pODEProblem[5] = init_Y_mRNA, pODEProblem[6] = init_Z_mRNA, pODEProblem[7] = init_X_protein, pODEProblem[8] = init_Y_protein, pODEProblem[9] = tau_mRNA_GFP, pODEProblem[10] = eff, pODEProblem[11] = init_Z_protein, pODEProblem[12] = tau_prot_GFP, pODEProblem[13] = eff_GFP, pODEProblem[14] = tps_active, pODEProblem[15] = init_GFP_mRNA, pODEProblem[16] = KM, pODEProblem[17] = init_X_mRNA, pODEProblem[18] = tau_prot, pODEProblem[19] = cell

	X_protein = pODEProblem[7] 
	GFP_mRNA = pODEProblem[15] 
	Y_mRNA = pODEProblem[5] 
	X_mRNA = pODEProblem[17] 
	Z_mRNA = pODEProblem[6] 
	Z_protein = pODEProblem[11] 
	Y_protein = pODEProblem[8] 
	GFP = pODEProblem[3] 

	 return [X_protein, GFP_mRNA, Y_mRNA, X_mRNA, Z_mRNA, Z_protein, Y_protein, GFP]
end

function compute_σ(u::AbstractVector, t::Real, θ_sd::AbstractVector, θ_dynamic::AbstractVector, θ_nonDynamic::AbstractVector, 
                   parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :fluorescence 
		noiseParameter1_fluorescence = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_fluorescence 
	end

end