function evalYmod(u, t, dynPar, obsPar, nonDynParam, paramData, observableId, mapObsParam) 

	IPR_S, IPR_I2, IPR_R, IPR_O, IPR_I1, IPR_A= u 
	k1, k2, k3, k4, k_1, k_2, k_3, k_4, l2, l4, l6, l_2, l_4, l_6 = dynPar 

	if observableId == "open_probability" 
		return ( 0.9 * IPR_A + 0.1 * IPR_O ) ^ 4 
	end

end

function evalU0!(u0Vec, paramVec) 

	l_4, k_4, IP3, k4, k_2, l2, l_2, default, l_6, k1, k_3, l6, membrane, k3, l4, Ca, k2, k_1 = paramVec 

	IPR_S = 0.0 
	IPR_I2 = 0.0 
	IPR_R = 1.0 
	IPR_O = 0.0 
	IPR_I1 = 0.0 
	IPR_A = 0.0 

	u0Vec .= IPR_S, IPR_I2, IPR_R, IPR_O, IPR_I1, IPR_A
end

function evalSd!(u, t, sdPar, dynPar, nonDynParam, paramData, observableId, mapSdParam) 

	IPR_S, IPR_I2, IPR_R, IPR_O, IPR_I1, IPR_A= u 
	k1, k2, k3, k4, k_1, k_2, k_3, k_4, l2, l4, l6, l_2, l_4, l_6 = dynPar 

	if observableId == "open_probability" 
		noiseParameter1_open_probability = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_open_probability 
	end

end