function evalYmod(u, t, dynPar, obsPar, nonDynParam, paramData, observableId, mapObsParam) 

	A_state, Y_state, Z_state= u 
	K2, K_par, Ka, Kd, Kf, Kp, Ky, Kz, Vd, Vm2, Vm3, Vp, beta_par, epsilon_par, init_A_state, init_Y_state, init_Z_state, n_par, v0, v1 = dynPar 

	if observableId == "Ca" 
		observableParameter1_Ca, observableParameter2_Ca = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_Ca + Z_state * observableParameter2_Ca 
	end

end

function evalU0!(u0Vec, paramVec) 

	v0, Ky, Vm3, K2, Kz, v1, Vm2, beta_par, init_Y_state, extracellular, n_par, K_par, Kd, cytosol, epsilon_par, intravesicular, Kp, Kf, init_A_state, Vd, init_Z_state, Vp, Ka = paramVec 

	A_state = init_A_state 
	Y_state = init_Y_state 
	Z_state = init_Z_state 

	u0Vec .= A_state, Y_state, Z_state
end

function evalU0(paramVec) 

	v0, Ky, Vm3, K2, Kz, v1, Vm2, beta_par, init_Y_state, extracellular, n_par, K_par, Kd, cytosol, epsilon_par, intravesicular, Kp, Kf, init_A_state, Vd, init_Z_state, Vp, Ka = paramVec 

	A_state = init_A_state 
	Y_state = init_Y_state 
	Z_state = init_Z_state 

	 return [A_state, Y_state, Z_state]
end

function evalSd!(u, t, sdPar, dynPar, nonDynParam, paramData, observableId, mapSdParam) 

	A_state, Y_state, Z_state= u 
	K2, K_par, Ka, Kd, Kf, Kp, Ky, Kz, Vd, Vm2, Vm3, Vp, beta_par, epsilon_par, init_A_state, init_Y_state, init_Z_state, n_par, v0, v1 = dynPar 

	if observableId == "Ca" 
		noiseParameter1_Ca = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_Ca 
	end

end