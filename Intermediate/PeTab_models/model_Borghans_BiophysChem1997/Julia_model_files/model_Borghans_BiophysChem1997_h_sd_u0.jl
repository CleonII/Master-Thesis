#u[1] = A_state, u[2] = Y_state, u[3] = Z_state
#θ_dynamicNames[1] = K2, θ_dynamicNames[2] = K_par, θ_dynamicNames[3] = Ka, θ_dynamicNames[4] = Kd, θ_dynamicNames[5] = Kf, θ_dynamicNames[6] = Kp, θ_dynamicNames[7] = Ky, θ_dynamicNames[8] = Kz, θ_dynamicNames[9] = Vd, θ_dynamicNames[10] = Vm2, θ_dynamicNames[11] = Vm3, θ_dynamicNames[12] = Vp, θ_dynamicNames[13] = beta_par, θ_dynamicNames[14] = epsilon_par, θ_dynamicNames[15] = init_A_state, θ_dynamicNames[16] = init_Y_state, θ_dynamicNames[17] = init_Z_state, θ_dynamicNames[18] = n_par, θ_dynamicNames[19] = v0, θ_dynamicNames[20] = v1
#

function compute_h(u::AbstractVector, t::Real, θ_dynamic::AbstractVector, θ_observable::AbstractVector,
                   θ_nonDynamic::AbstractVector, parameterInfo::ParametersInfo, observableId::Symbol, 
                      parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :Ca 
		observableParameter1_Ca, observableParameter2_Ca = getObsOrSdParam(θ_observable, parameterMap)
		return observableParameter1_Ca + u[3] * observableParameter2_Ca 
	end

end

function compute_u0!(u0::AbstractVector, pODEProblem::AbstractVector) 

	#pODEProblem[1] = v0, pODEProblem[2] = Ky, pODEProblem[3] = Vm3, pODEProblem[4] = K2, pODEProblem[5] = Kz, pODEProblem[6] = v1, pODEProblem[7] = Vm2, pODEProblem[8] = beta_par, pODEProblem[9] = init_Y_state, pODEProblem[10] = extracellular, pODEProblem[11] = n_par, pODEProblem[12] = K_par, pODEProblem[13] = Kd, pODEProblem[14] = cytosol, pODEProblem[15] = epsilon_par, pODEProblem[16] = intravesicular, pODEProblem[17] = Kp, pODEProblem[18] = Kf, pODEProblem[19] = init_A_state, pODEProblem[20] = Vd, pODEProblem[21] = init_Z_state, pODEProblem[22] = Vp, pODEProblem[23] = Ka

	A_state = pODEProblem[19] 
	Y_state = pODEProblem[9] 
	Z_state = pODEProblem[21] 

	u0 .= A_state, Y_state, Z_state
end

function compute_u0(pODEProblem::AbstractVector)::AbstractVector 

	#pODEProblem[1] = v0, pODEProblem[2] = Ky, pODEProblem[3] = Vm3, pODEProblem[4] = K2, pODEProblem[5] = Kz, pODEProblem[6] = v1, pODEProblem[7] = Vm2, pODEProblem[8] = beta_par, pODEProblem[9] = init_Y_state, pODEProblem[10] = extracellular, pODEProblem[11] = n_par, pODEProblem[12] = K_par, pODEProblem[13] = Kd, pODEProblem[14] = cytosol, pODEProblem[15] = epsilon_par, pODEProblem[16] = intravesicular, pODEProblem[17] = Kp, pODEProblem[18] = Kf, pODEProblem[19] = init_A_state, pODEProblem[20] = Vd, pODEProblem[21] = init_Z_state, pODEProblem[22] = Vp, pODEProblem[23] = Ka

	A_state = pODEProblem[19] 
	Y_state = pODEProblem[9] 
	Z_state = pODEProblem[21] 

	 return [A_state, Y_state, Z_state]
end

function compute_σ(u::AbstractVector, t::Real, θ_sd::AbstractVector, θ_dynamic::AbstractVector, θ_nonDynamic::AbstractVector, 
                   parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :Ca 
		noiseParameter1_Ca = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_Ca 
	end

end