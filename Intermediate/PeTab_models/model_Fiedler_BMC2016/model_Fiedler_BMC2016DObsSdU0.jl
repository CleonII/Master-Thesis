function evalDYmodDu(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	RAF, MEK, pMEK, pERK, pRAF, ERK= u 
	ERK_total, UO126, k10, RAF_total, K_3, cyt, k4, Sorafenib, K_2, k6, k11, tau1, MEK_total, K_1, k3, tau2, k5, k2 = p 
	if observableId == :pErk 
		observableParameter1_pErk = getObsOrSdParam(obsPar, mapObsParam)
		out[4] = observableParameter1_pErk
		return nothing
	end

	if observableId == :pMek 
		observableParameter1_pMek = getObsOrSdParam(obsPar, mapObsParam)
		out[3] = observableParameter1_pMek
		return nothing
	end

end

function evalDYmodDp(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	RAF, MEK, pMEK, pERK, pRAF, ERK= u 
	ERK_total, UO126, k10, RAF_total, K_3, cyt, k4, Sorafenib, K_2, k6, k11, tau1, MEK_total, K_1, k3, tau2, k5, k2 = p 
	if observableId == :pErk 
		return nothing
	end

	if observableId == :pMek 
		return nothing
	end

end

function evalDSdDu!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	RAF, MEK, pMEK, pERK, pRAF, ERK= u 
	ERK_total, UO126, k10, RAF_total, K_3, cyt, k4, Sorafenib, K_2, k6, k11, tau1, MEK_total, K_1, k3, tau2, k5, k2 = p 
	if observableId == "pErk" 
		return nothing
	end

	if observableId == "pMek" 
		return nothing
	end

end

function evalDSdDp!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	RAF, MEK, pMEK, pERK, pRAF, ERK= u 
	ERK_total, UO126, k10, RAF_total, K_3, cyt, k4, Sorafenib, K_2, k6, k11, tau1, MEK_total, K_1, k3, tau2, k5, k2 = p 
	if observableId == "pErk" 
		return nothing
	end

	if observableId == "pMek" 
		return nothing
	end

end

