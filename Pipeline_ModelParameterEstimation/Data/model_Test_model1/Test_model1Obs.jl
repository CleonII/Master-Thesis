function Test_model1(u, t, dynPar, obsPar, paramData, obsData, observableId, simulationId) 

	sebastian, damiano, dummyVariable= u 
	alpha, beta, gamma, delta = dynPar 
	alpha_scale_C = paramData.paramVal[7] 

	if observableId == "sebastian_measurement" 
		return sebastian 
	end

	if observableId == "damiano_measurement" 
		return damiano 
	end

end

function Test_model1_t0!(u0Vec, paramVec) 

	beta, alpha, alpha_scale, delta, gamma = paramVec 

	sebastian = 8.0 
	damiano = 4.0 
	dummyVariable = 0.0 

	u0Vec .= sebastian, damiano, dummyVariable
end

function Test_model1_sd!(u, t, sdPar, dynPar, paramData, obsData, observableId, simulationId) 

	sebastian, damiano, dummyVariable= u 
	alpha, beta, gamma, delta = dynPar 
	alpha_scale_C = paramData.paramVal[7] 

	if observableId == "sebastian_measurement" 
		noiseParameter1_sebastian_measurement = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_sebastian_measurement 
	end

	if observableId == "damiano_measurement" 
		noiseParameter1_damiano_measurement = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_damiano_measurement 
	end

end