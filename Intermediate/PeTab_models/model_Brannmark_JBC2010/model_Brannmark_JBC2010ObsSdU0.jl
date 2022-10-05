function evalYmod(u, t, dynPar, obsPar, paramData, observableId, mapObsParam) 

	IRp, IR, IRins, IRiP, IRS, X, IRi, IRSiP, Xp, dummyVariable= u 
	k1a, k1aBasic, k1b, k1c, k1d, k1e, k1f, k1g, k1r, k21, k22, k3, km2, km3 = dynPar 

	if observableId == "IR1_P" 
		observableParameter1_IR1_P = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_IR1_P * ( IRp + IRiP ) 
	end

	if observableId == "IRS1_P" 
		observableParameter1_IRS1_P = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_IRS1_P * IRSiP 
	end

	if observableId == "IRS1_P_DosR" 
		observableParameter1_IRS1_P_DosR = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_IRS1_P_DosR * IRSiP 
	end

end

function evalU0!(u0Vec, paramVec) 

	k1g, k1d, k1c, k1b, insulin_dose_2, insulin_time_1, k1a, insulin_time_2, k1aBasic, k1r, insulin_dose_1, k1e, k1f, k21, km2, k22, km3, k3, k_IRSiP_DosR, default = paramVec 

	IRp = 1.7629010620181e-9 
	IR = 9.94957642787569 
	IRins = 0.0173972221725393 
	IRiP = 1.11590026152296e-5 
	IRS = 9.86699348701367 
	X = 9.99984199487351 
	IRi = 0.0330151891862681 
	IRSiP = 0.133006512986336 
	Xp = 0.000158005126497888 
	dummyVariable = 0.0 

	u0Vec .= IRp, IR, IRins, IRiP, IRS, X, IRi, IRSiP, Xp, dummyVariable
end

function evalSd!(u, t, sdPar, dynPar, paramData, observableId, mapSdParam) 

	IRp, IR, IRins, IRiP, IRS, X, IRi, IRSiP, Xp, dummyVariable= u 
	k1a, k1aBasic, k1b, k1c, k1d, k1e, k1f, k1g, k1r, k21, k22, k3, km2, km3 = dynPar 

	if observableId == "IR1_P" 
		noiseParameter1_IR1_P = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_IR1_P 
	end

	if observableId == "IRS1_P" 
		noiseParameter1_IRS1_P = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_IRS1_P 
	end

	if observableId == "IRS1_P_DosR" 
		noiseParameter1_IRS1_P_DosR = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_IRS1_P_DosR 
	end

end