function evalYmod(u, t, dynPar, obsPar, paramData, observableId, mapObsParam) 

	IR2, IR2in, Rec2, IR1in, Uptake1, Uptake2, InsulinFragments, IR1, Rec1, Ins, BoundUnspec, dummyVariable= u 
	ini_R1, ini_R2fold, ka1, ka2fold, kd1, kd2fold, kin, kin2, koff_unspec, kon_unspec, kout, kout2, kout_frag = dynPar 

	if observableId == "observable_IR1" 
		observableParameter2_observable_IR1, observableParameter1_observable_IR1 = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter2_observable_IR1 * ( IR1 + IR1in + observableParameter1_observable_IR1 ) 
	end

	if observableId == "observable_IR2" 
		observableParameter2_observable_IR2, observableParameter1_observable_IR2 = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter2_observable_IR2 * ( IR2 + IR2in + observableParameter1_observable_IR2 ) 
	end

	if observableId == "observable_IRsum" 
		observableParameter2_observable_IRsum, observableParameter1_observable_IRsum = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter2_observable_IRsum * ( 0.605 * IR1 + 0.395 * IR2 + 0.605 * IR1in + 0.395 * IR2in + observableParameter1_observable_IRsum ) 
	end

	if observableId == "observable_Insulin" 
		observableParameter1_observable_Insulin, observableParameter2_observable_Insulin, observableParameter3_observable_Insulin, observableParameter3_observable_Insulin, observableParameter4_observable_Insulin = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_observable_Insulin + observableParameter2_observable_Insulin * ( Ins + InsulinFragments * observableParameter3_observable_Insulin ) / ( ( Ins + InsulinFragments * observableParameter3_observable_Insulin ) / observableParameter4_observable_Insulin + 1 ) 
	end

end

function evalU0!(u0Vec, paramVec) 

	kout2, kd2fold, ka2fold, kd1, kin2, ka1, kout_frag, kin, kout, koff_unspec, kon_unspec, ini_R2fold, ini_R1, init_Ins = paramVec 

	IR2 = 0.0 
	IR2in = 0.0 
	Rec2 = ini_R1 * ini_R2fold 
	IR1in = 0.0 
	Uptake1 = 0.0 
	Uptake2 = 0.0 
	InsulinFragments = 0.0 
	IR1 = 0.0 
	Rec1 = ini_R1 
	Ins = init_Ins 
	BoundUnspec = 0.0 
	dummyVariable = 0.0 

	u0Vec .= IR2, IR2in, Rec2, IR1in, Uptake1, Uptake2, InsulinFragments, IR1, Rec1, Ins, BoundUnspec, dummyVariable
end

function evalSd!(u, t, sdPar, dynPar, paramData, observableId, mapSdParam) 

	IR2, IR2in, Rec2, IR1in, Uptake1, Uptake2, InsulinFragments, IR1, Rec1, Ins, BoundUnspec, dummyVariable= u 
	ini_R1, ini_R2fold, ka1, ka2fold, kd1, kd2fold, kin, kin2, koff_unspec, kon_unspec, kout, kout2, kout_frag = dynPar 

	if observableId == "observable_IR1" 
		noiseParameter1_observable_IR1 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_IR1 
	end

	if observableId == "observable_IR2" 
		noiseParameter1_observable_IR2 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_IR2 
	end

	if observableId == "observable_IRsum" 
		noiseParameter1_observable_IRsum = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_IRsum 
	end

	if observableId == "observable_Insulin" 
		noiseParameter1_observable_Insulin = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_Insulin 
	end

end