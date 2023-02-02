#u[1] = A_state, u[2] = Y_state, u[3] = Z_state
#pODEProblem[1] = v0, pODEProblem[2] = Ky, pODEProblem[3] = Vm3, pODEProblem[4] = K2, pODEProblem[5] = Kz, pODEProblem[6] = v1, pODEProblem[7] = Vm2, pODEProblem[8] = beta_par, pODEProblem[9] = init_Y_state, pODEProblem[10] = extracellular, pODEProblem[11] = n_par, pODEProblem[12] = K_par, pODEProblem[13] = Kd, pODEProblem[14] = cytosol, pODEProblem[15] = epsilon_par, pODEProblem[16] = intravesicular, pODEProblem[17] = Kp, pODEProblem[18] = Kf, pODEProblem[19] = init_A_state, pODEProblem[20] = Vd, pODEProblem[21] = init_Z_state, pODEProblem[22] = Vp, pODEProblem[23] = Ka
#
function compute_∂h∂u!(u, t::Real, pODEProblem::AbstractVector, θ_observable::AbstractVector,                    
                       θ_nonDynamic::AbstractVector, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :Ca 
		observableParameter1_Ca, observableParameter2_Ca = getObsOrSdParam(θ_observable, parameterMap)
		out[3] = observableParameter2_Ca
		return nothing
	end

end

function compute_∂h∂p!(u, t::Real, pODEProblem::AbstractVector, θ_observable::AbstractVector,
                       θ_nonDynamic::AbstractVector, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :Ca 
		return nothing
	end

end

function compute_∂σ∂σu!(u, t::Real, θ_sd::AbstractVector, pODEProblem::AbstractVector, θ_nonDynamic::AbstractVector, 
                        parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :Ca 
		return nothing
	end

end

function compute_∂σ∂σp!(u, t::Real, θ_sd::AbstractVector, pODEProblem::AbstractVector, θ_nonDynamic::AbstractVector, 
                        parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :Ca 
		return nothing
	end

end

