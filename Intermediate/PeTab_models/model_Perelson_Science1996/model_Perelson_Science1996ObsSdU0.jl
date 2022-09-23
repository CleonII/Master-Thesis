function evalYmod(u, t, dynPar, obsPar, paramData, obsData, observableId, simulationId) 

	Vni, V, Vin, Tstar, dummyVariable= u 
	c, delta = dynPar 
	NN_C = paramData.paramVal[1] 
	T0_C = paramData.paramVal[2] 
	K0_C = paramData.paramVal[5] 

	if observableId == "task0_model0_perelson1_V" 
		return V 
	end

end

function evalU0!(u0Vec, paramVec) 

	NN, delta, default, c, T0, K0 = paramVec 

	Vni = 0.0 
	V = 1.86e6 
	Vin = 1.86e6 
	Tstar = 15061.32075 
	dummyVariable = 0.0 

	u0Vec .= Vni, V, Vin, Tstar, dummyVariable
end

function evalSd!(u, t, sdPar, dynPar, paramData, obsData, observableId, simulationId) 

	Vni, V, Vin, Tstar, dummyVariable= u 
	c, delta = dynPar 
	NN_C = paramData.paramVal[1] 
	T0_C = paramData.paramVal[2] 
	K0_C = paramData.paramVal[5] 

	if observableId == "task0_model0_perelson1_V" 
		noiseParameter1_task0_model0_perelson1_V = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_task0_model0_perelson1_V 
	end

end