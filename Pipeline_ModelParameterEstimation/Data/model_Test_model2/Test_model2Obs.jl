function Test_model2(u, t, dynPar, obsPar, paramData, obsData, observableId, simulationId) 

	sebastian, damiano= u 
	alpha, beta = dynPar 

	if observableId == "sebastian_measurement" 
		return sebastian 
	end

	if observableId == "damiano_measurement" 
		return damiano 
	end

end

function Test_model2_t0!(u0Vec, paramVec) 

	alpha, beta = paramVec 

	sebastian = 8 
	damiano = 4 

	u0Vec .= sebastian, damiano
end

function Test_model2_sd!(u, t, sdPar, dynPar, paramData, obsData, observableId, simulationId) 

	sebastian, damiano= u 
	alpha, beta = dynPar 

	if observableId == "sebastian_measurement" 
		noiseParameter1_sebastian_measurement = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_sebastian_measurement 
	end

	if observableId == "damiano_measurement" 
		noiseParameter1_damiano_measurement = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_damiano_measurement 
	end

end