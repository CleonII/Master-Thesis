function evalDYmodDu(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	IRp, IR, IRins, IRiP, IRS, X, IRi, IRSiP, Xp= u 
	k1c, k21, insulin_bool1, k1g, insulin_dose_2, k1a, insulin_dose_1, k1aBasic, k1d, insulin_time_1, insulin_time_2, cyt, k22, insulin_bool2, default, k1r, k1f, k1b, k3, km2, k1e, k_IRSiP_DosR, km3 = p 
	if observableId == :IR1_P 
		observableParameter1_IR1_P = getObsOrSdParam(obsPar, mapObsParam)
		out[1] = observableParameter1_IR1_P
		out[4] = observableParameter1_IR1_P
		return nothing
	end

	if observableId == :IRS1_P 
		observableParameter1_IRS1_P = getObsOrSdParam(obsPar, mapObsParam)
		out[8] = observableParameter1_IRS1_P
		return nothing
	end

	if observableId == :IRS1_P_DosR 
		observableParameter1_IRS1_P_DosR = getObsOrSdParam(obsPar, mapObsParam)
		out[8] = observableParameter1_IRS1_P_DosR
		return nothing
	end

end

function evalDYmodDp(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	IRp, IR, IRins, IRiP, IRS, X, IRi, IRSiP, Xp= u 
	k1c, k21, insulin_bool1, k1g, insulin_dose_2, k1a, insulin_dose_1, k1aBasic, k1d, insulin_time_1, insulin_time_2, cyt, k22, insulin_bool2, default, k1r, k1f, k1b, k3, km2, k1e, k_IRSiP_DosR, km3 = p 
	if observableId == :IR1_P 
		return nothing
	end

	if observableId == :IRS1_P 
		return nothing
	end

	if observableId == :IRS1_P_DosR 
		return nothing
	end

end

function evalDSdDu!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	IRp, IR, IRins, IRiP, IRS, X, IRi, IRSiP, Xp= u 
	k1c, k21, insulin_bool1, k1g, insulin_dose_2, k1a, insulin_dose_1, k1aBasic, k1d, insulin_time_1, insulin_time_2, cyt, k22, insulin_bool2, default, k1r, k1f, k1b, k3, km2, k1e, k_IRSiP_DosR, km3 = p 
	if observableId == "IR1_P" 
		return nothing
	end

	if observableId == "IRS1_P" 
		return nothing
	end

	if observableId == "IRS1_P_DosR" 
		return nothing
	end

end

function evalDSdDp!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	IRp, IR, IRins, IRiP, IRS, X, IRi, IRSiP, Xp= u 
	k1c, k21, insulin_bool1, k1g, insulin_dose_2, k1a, insulin_dose_1, k1aBasic, k1d, insulin_time_1, insulin_time_2, cyt, k22, insulin_bool2, default, k1r, k1f, k1b, k3, km2, k1e, k_IRSiP_DosR, km3 = p 
	if observableId == "IR1_P" 
		return nothing
	end

	if observableId == "IRS1_P" 
		return nothing
	end

	if observableId == "IRS1_P_DosR" 
		return nothing
	end

end

