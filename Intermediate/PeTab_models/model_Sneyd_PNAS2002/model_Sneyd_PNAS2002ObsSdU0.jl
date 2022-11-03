function evalYmod(u, t, dynPar, obsPar, paramData, observableId, mapObsParam) 

	IPR_S, IPR_I2, IPR_R, IPR_O, IPR_I1, IPR_A, dummyVariable= u 
	k1, k2, k3, k4, k_1, k_2, k_3, k_4, l2, l4, l6, l_2, l_4, l_6 = dynPar 

	if observableId == "open_probability" 
		return ( 0.9 * IPR_A + 0.1 * IPR_O ) ^ 4 
	end

end

function evalU0!(u0Vec, paramVec) 

	Ca, k_4, l6, membrane, k_3, k4, k3, l_6, k_1, l2, k1, l_2, l4, l_4, k2, IP3, k_2, default = paramVec 

	IPR_S = 0.0 
	IPR_I2 = 0.0 
	IPR_R = 1.0 
	IPR_O = 0.0 
	IPR_I1 = 0.0 
	IPR_A = 0.0 
	dummyVariable = 0.0 

	u0Vec .= IPR_S, IPR_I2, IPR_R, IPR_O, IPR_I1, IPR_A, dummyVariable
end

function evalSd!(u, t, sdPar, dynPar, paramData, observableId, mapSdParam) 

	IPR_S, IPR_I2, IPR_R, IPR_O, IPR_I1, IPR_A, dummyVariable= u 
	k1, k2, k3, k4, k_1, k_2, k_3, k_4, l2, l4, l6, l_2, l_4, l_6 = dynPar 

	if observableId == "open_probability" 
		noiseParameter1_open_probability = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_open_probability 
	end

end