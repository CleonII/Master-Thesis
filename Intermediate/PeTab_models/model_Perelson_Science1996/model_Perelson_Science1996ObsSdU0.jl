function evalYmod(u, t, dynPar, obsPar, nonDynParam, paramData, observableId, mapObsParam) 

	Vni, V, Vin, Tstar= u 
	c, delta = dynPar 
	NN_C = paramData.paramVal[1] 
	T0_C = paramData.paramVal[2] 
	K0_C = paramData.paramVal[5] 

	if observableId == "task0_model0_perelson1_V" 
		return V 
	end

end

function evalU0!(u0Vec, paramVec) 

	c, T0, default, K0, NN, delta = paramVec 

	Vni = 0.0 
	V = 1.86e6 
	Vin = 1.86e6 
	Tstar = 15061.32075 

	u0Vec .= Vni, V, Vin, Tstar
end

function evalSd!(u, t, sdPar, dynPar, nonDynParam, paramData, observableId, mapSdParam) 

	Vni, V, Vin, Tstar= u 
	c, delta = dynPar 
	NN_C = paramData.paramVal[1] 
	T0_C = paramData.paramVal[2] 
	K0_C = paramData.paramVal[5] 

	if observableId == "task0_model0_perelson1_V" 
		noiseParameter1_task0_model0_perelson1_V = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_task0_model0_perelson1_V 
	end

end