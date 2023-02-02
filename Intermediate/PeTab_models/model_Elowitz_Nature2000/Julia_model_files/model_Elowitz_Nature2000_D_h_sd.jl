#u[1] = X_protein, u[2] = GFP_mRNA, u[3] = Y_mRNA, u[4] = X_mRNA, u[5] = Z_mRNA, u[6] = Z_protein, u[7] = Y_protein, u[8] = GFP
#pODEProblem[1] = tau_mRNA, pODEProblem[2] = tps_repr, pODEProblem[3] = init_GFP, pODEProblem[4] = n_Hill, pODEProblem[5] = init_Y_mRNA, pODEProblem[6] = init_Z_mRNA, pODEProblem[7] = init_X_protein, pODEProblem[8] = init_Y_protein, pODEProblem[9] = tau_mRNA_GFP, pODEProblem[10] = eff, pODEProblem[11] = init_Z_protein, pODEProblem[12] = tau_prot_GFP, pODEProblem[13] = eff_GFP, pODEProblem[14] = tps_active, pODEProblem[15] = init_GFP_mRNA, pODEProblem[16] = KM, pODEProblem[17] = init_X_mRNA, pODEProblem[18] = tau_prot, pODEProblem[19] = cell
#
function compute_∂h∂u!(u, t::Real, pODEProblem::AbstractVector, θ_observable::AbstractVector,                    
                       θ_nonDynamic::AbstractVector, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :fluorescence 
		observableParameter1_fluorescence, observableParameter2_fluorescence = getObsOrSdParam(θ_observable, parameterMap)
		out[8] = observableParameter2_fluorescence
		return nothing
	end

end

function compute_∂h∂p!(u, t::Real, pODEProblem::AbstractVector, θ_observable::AbstractVector,
                       θ_nonDynamic::AbstractVector, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :fluorescence 
		return nothing
	end

end

function compute_∂σ∂σu!(u, t::Real, θ_sd::AbstractVector, pODEProblem::AbstractVector, θ_nonDynamic::AbstractVector, 
                        parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :fluorescence 
		return nothing
	end

end

function compute_∂σ∂σp!(u, t::Real, θ_sd::AbstractVector, pODEProblem::AbstractVector, θ_nonDynamic::AbstractVector, 
                        parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :fluorescence 
		return nothing
	end

end

