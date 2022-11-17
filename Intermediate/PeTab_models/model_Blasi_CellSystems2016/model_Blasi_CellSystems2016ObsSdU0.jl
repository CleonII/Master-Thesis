function evalYmod(u, t, dynPar, obsPar, nonDynParam, paramData, observableId, mapObsParam) 

	x_k5k12k16, x_k8, x_k16, x_0ac, x_k12, x_k5k8, x_k5k12, x_k12k16, x_k8k12k16, x_k5, x_k5k16, x_k5k8k12, x_k8k12, x_4ac, x_k8k16, x_k5k8k16= u 
	a_basal, a_k8, a_k5_k5k12, a_k12_k5k12, a_k16_k12k16, a_k5k12_k5k8k12, a_k12k16_k8k12k16, a_k8k12k16_4ac = dynPar 
	d_C = paramData.paramVal[9] 

	if observableId == "observable_0ac" 
		return x_0ac 
	end

	if observableId == "observable_4ac" 
		return x_4ac 
	end

	if observableId == "observable_k12" 
		return x_k12 
	end

	if observableId == "observable_k12k16" 
		return x_k12k16 
	end

	if observableId == "observable_k16" 
		return x_k16 
	end

	if observableId == "observable_k5" 
		return x_k5 
	end

	if observableId == "observable_k5k12" 
		return x_k5k12 
	end

	if observableId == "observable_k5k12k16" 
		return x_k5k12k16 
	end

	if observableId == "observable_k5k16" 
		return x_k5k16 
	end

	if observableId == "observable_k5k8" 
		return x_k5k8 
	end

	if observableId == "observable_k5k8k12" 
		return x_k5k8k12 
	end

	if observableId == "observable_k5k8k16" 
		return x_k5k8k16 
	end

	if observableId == "observable_k8" 
		return x_k8 
	end

	if observableId == "observable_k8k12" 
		return x_k8k12 
	end

	if observableId == "observable_k8k12k16" 
		return x_k8k12k16 
	end

	if observableId == "observable_k8k16" 
		return x_k8k16 
	end

end

function evalU0!(u0Vec, paramVec) 

	a_k5_k5k12, a_k8, d, default, a_k12k16_k8k12k16, a_basal, a_k5k12_k5k8k12, a_k8k12k16_4ac, a_k12_k5k12, a_k16_k12k16 = paramVec 

	x_k5k12k16 = 0.0 
	x_k8 = 0.0 
	x_k16 = 0.0 
	x_0ac = 1.0 
	x_k12 = 0.0 
	x_k5k8 = 0.0 
	x_k5k12 = 0.0 
	x_k12k16 = 0.0 
	x_k8k12k16 = 0.0 
	x_k5 = 0.0 
	x_k5k16 = 0.0 
	x_k5k8k12 = 0.0 
	x_k8k12 = 0.0 
	x_4ac = 0.0 
	x_k8k16 = 0.0 
	x_k5k8k16 = 0.0 

	u0Vec .= x_k5k12k16, x_k8, x_k16, x_0ac, x_k12, x_k5k8, x_k5k12, x_k12k16, x_k8k12k16, x_k5, x_k5k16, x_k5k8k12, x_k8k12, x_4ac, x_k8k16, x_k5k8k16
end

function evalSd!(u, t, sdPar, dynPar, nonDynParam, paramData, observableId, mapSdParam) 

	x_k5k12k16, x_k8, x_k16, x_0ac, x_k12, x_k5k8, x_k5k12, x_k12k16, x_k8k12k16, x_k5, x_k5k16, x_k5k8k12, x_k8k12, x_4ac, x_k8k16, x_k5k8k16= u 
	a_basal, a_k8, a_k5_k5k12, a_k12_k5k12, a_k16_k12k16, a_k5k12_k5k8k12, a_k12k16_k8k12k16, a_k8k12k16_4ac = dynPar 
	d_C = paramData.paramVal[9] 

	if observableId == "observable_0ac" 
		noiseParameter1_observable_0ac = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_0ac 
	end

	if observableId == "observable_4ac" 
		noiseParameter1_observable_4ac = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_4ac 
	end

	if observableId == "observable_k12" 
		noiseParameter1_observable_k12 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_k12 
	end

	if observableId == "observable_k12k16" 
		noiseParameter1_observable_k12k16 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_k12k16 
	end

	if observableId == "observable_k16" 
		noiseParameter1_observable_k16 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_k16 
	end

	if observableId == "observable_k5" 
		noiseParameter1_observable_k5 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_k5 
	end

	if observableId == "observable_k5k12" 
		noiseParameter1_observable_k5k12 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_k5k12 
	end

	if observableId == "observable_k5k12k16" 
		noiseParameter1_observable_k5k12k16 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_k5k12k16 
	end

	if observableId == "observable_k5k16" 
		noiseParameter1_observable_k5k16 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_k5k16 
	end

	if observableId == "observable_k5k8" 
		noiseParameter1_observable_k5k8 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_k5k8 
	end

	if observableId == "observable_k5k8k12" 
		noiseParameter1_observable_k5k8k12 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_k5k8k12 
	end

	if observableId == "observable_k5k8k16" 
		noiseParameter1_observable_k5k8k16 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_k5k8k16 
	end

	if observableId == "observable_k8" 
		noiseParameter1_observable_k8 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_k8 
	end

	if observableId == "observable_k8k12" 
		noiseParameter1_observable_k8k12 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_k8k12 
	end

	if observableId == "observable_k8k12k16" 
		noiseParameter1_observable_k8k12k16 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_k8k12k16 
	end

	if observableId == "observable_k8k16" 
		noiseParameter1_observable_k8k16 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_k8k16 
	end

end