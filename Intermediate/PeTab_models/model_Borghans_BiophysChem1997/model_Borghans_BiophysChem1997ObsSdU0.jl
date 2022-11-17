function evalYmod(u, t, dynPar, obsPar, nonDynParam, paramData, observableId, mapObsParam) 

	A_state, Y_state, Z_state, dummyVariable= u 
	K2, K_par, Ka, Kd, Kf, Kp, Ky, Kz, Vd, Vm2, Vm3, Vp, beta_par, epsilon_par, init_A_state, init_Y_state, init_Z_state, n_par, v0, v1 = dynPar 

	if observableId == "Ca" 
		observableParameter1_Ca, observableParameter2_Ca = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_Ca + Z_state * observableParameter2_Ca 
	end

end

function evalU0!(u0Vec, paramVec) 

	Vp, Kd, beta_par, Kp, n_par, Vd, epsilon_par, cytosol, Kf, K2, Vm2, Vm3, Ky, Kz, Ka, intravesicular, v0, v1, K_par, extracellular, init_Y_state, init_Z_state, init_A_state = paramVec 

	A_state = init_A_state 
	Y_state = init_Y_state 
	Z_state = init_Z_state 
	dummyVariable = 0.0 

	u0Vec .= A_state, Y_state, Z_state, dummyVariable
end

function evalSd!(u, t, sdPar, dynPar, nonDynParam, paramData, observableId, mapSdParam) 

	A_state, Y_state, Z_state, dummyVariable= u 
	K2, K_par, Ka, Kd, Kf, Kp, Ky, Kz, Vd, Vm2, Vm3, Vp, beta_par, epsilon_par, init_A_state, init_Y_state, init_Z_state, n_par, v0, v1 = dynPar 

	if observableId == "Ca" 
		noiseParameter1_Ca = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_Ca 
	end

end