function evalDYmodDu(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	A_state, Y_state, Z_state= u 
	v0, Ky, Vm3, K2, Kz, v1, Vm2, beta_par, init_Y_state, extracellular, n_par, K_par, Kd, cytosol, epsilon_par, intravesicular, Kp, Kf, init_A_state, Vd, init_Z_state, Vp, Ka = p 
	if observableId == "Ca" 
		observableParameter1_Ca, observableParameter2_Ca = getObsOrSdParam(obsPar, mapObsParam)
		out[3] = observableParameter2_Ca
		return nothing
	end

end

function evalDYmodDp(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	A_state, Y_state, Z_state= u 
	v0, Ky, Vm3, K2, Kz, v1, Vm2, beta_par, init_Y_state, extracellular, n_par, K_par, Kd, cytosol, epsilon_par, intravesicular, Kp, Kf, init_A_state, Vd, init_Z_state, Vp, Ka = p 
	if observableId == "Ca" 
		return nothing
	end

end

function evalDSdDu!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	A_state, Y_state, Z_state= u 
	v0, Ky, Vm3, K2, Kz, v1, Vm2, beta_par, init_Y_state, extracellular, n_par, K_par, Kd, cytosol, epsilon_par, intravesicular, Kp, Kf, init_A_state, Vd, init_Z_state, Vp, Ka = p 
	if observableId == "Ca" 
		return nothing
	end

end

function evalDSdDp!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	A_state, Y_state, Z_state= u 
	v0, Ky, Vm3, K2, Kz, v1, Vm2, beta_par, init_Y_state, extracellular, n_par, K_par, Kd, cytosol, epsilon_par, intravesicular, Kp, Kf, init_A_state, Vd, init_Z_state, Vp, Ka = p 
	if observableId == "Ca" 
		return nothing
	end

end

