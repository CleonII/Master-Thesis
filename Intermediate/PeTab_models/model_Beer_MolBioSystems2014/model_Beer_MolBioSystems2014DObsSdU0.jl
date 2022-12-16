function evalDYmodDu(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	Glu, cGlu, Ind, Bac= u 
	kdegi, medium, Bacmax, ksyn, kdim, tau, init_Bac, beta = p 
	if observableId == "Bacnorm" 
		out[4] = 1
		return nothing
	end

	if observableId == "IndconcNormRange" 
		out[3] = 1
		return nothing
	end

end

function evalDYmodDp(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	Glu, cGlu, Ind, Bac= u 
	kdegi, medium, Bacmax, ksyn, kdim, tau, init_Bac, beta = p 
	if observableId == "Bacnorm" 
		return nothing
	end

	if observableId == "IndconcNormRange" 
		return nothing
	end

end

function evalDSdDu!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	Glu, cGlu, Ind, Bac= u 
	kdegi, medium, Bacmax, ksyn, kdim, tau, init_Bac, beta = p 
	if observableId == "Bacnorm" 
		return nothing
	end

	if observableId == "IndconcNormRange" 
		return nothing
	end

end

function evalDSdDp!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	Glu, cGlu, Ind, Bac= u 
	kdegi, medium, Bacmax, ksyn, kdim, tau, init_Bac, beta = p 
	if observableId == "Bacnorm" 
		return nothing
	end

	if observableId == "IndconcNormRange" 
		return nothing
	end

end

