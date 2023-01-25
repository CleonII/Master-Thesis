function evalYmod(u, t, dynPar, obsPar, nonDynParam, paramData, observableId, mapObsParam) 

	X_protein, GFP_mRNA, Y_mRNA, X_mRNA, Z_mRNA, Z_protein, Y_protein, GFP= u 
	KM, eff, eff_GFP, init_GFP, init_GFP_mRNA, init_X_mRNA, init_X_protein, init_Y_mRNA, init_Y_protein, init_Z_mRNA, init_Z_protein, n_Hill, tau_mRNA, tau_mRNA_GFP, tau_prot, tau_prot_GFP, tps_active, tps_repr = dynPar 

	if observableId == :fluorescence 
		observableParameter1_fluorescence, observableParameter2_fluorescence = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_fluorescence + GFP * observableParameter2_fluorescence 
	end

end

function evalU0!(u0Vec, paramVec) 

	tau_mRNA, tps_repr, init_GFP, n_Hill, init_Y_mRNA, init_Z_mRNA, init_X_protein, init_Y_protein, tau_mRNA_GFP, eff, init_Z_protein, tau_prot_GFP, eff_GFP, tps_active, init_GFP_mRNA, KM, init_X_mRNA, tau_prot, cell = paramVec 

	X_protein = init_X_protein 
	GFP_mRNA = init_GFP_mRNA 
	Y_mRNA = init_Y_mRNA 
	X_mRNA = init_X_mRNA 
	Z_mRNA = init_Z_mRNA 
	Z_protein = init_Z_protein 
	Y_protein = init_Y_protein 
	GFP = init_GFP 

	u0Vec .= X_protein, GFP_mRNA, Y_mRNA, X_mRNA, Z_mRNA, Z_protein, Y_protein, GFP
end

function evalU0(paramVec) 

	tau_mRNA, tps_repr, init_GFP, n_Hill, init_Y_mRNA, init_Z_mRNA, init_X_protein, init_Y_protein, tau_mRNA_GFP, eff, init_Z_protein, tau_prot_GFP, eff_GFP, tps_active, init_GFP_mRNA, KM, init_X_mRNA, tau_prot, cell = paramVec 

	X_protein = init_X_protein 
	GFP_mRNA = init_GFP_mRNA 
	Y_mRNA = init_Y_mRNA 
	X_mRNA = init_X_mRNA 
	Z_mRNA = init_Z_mRNA 
	Z_protein = init_Z_protein 
	Y_protein = init_Y_protein 
	GFP = init_GFP 

	 return [X_protein, GFP_mRNA, Y_mRNA, X_mRNA, Z_mRNA, Z_protein, Y_protein, GFP]
end

function evalSd!(u, t, sdPar, dynPar, nonDynParam, paramData, observableId, mapSdParam) 

	X_protein, GFP_mRNA, Y_mRNA, X_mRNA, Z_mRNA, Z_protein, Y_protein, GFP= u 
	KM, eff, eff_GFP, init_GFP, init_GFP_mRNA, init_X_mRNA, init_X_protein, init_Y_mRNA, init_Y_protein, init_Z_mRNA, init_Z_protein, n_Hill, tau_mRNA, tau_mRNA_GFP, tau_prot, tau_prot_GFP, tps_active, tps_repr = dynPar 

	if observableId == :fluorescence 
		noiseParameter1_fluorescence = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_fluorescence 
	end

end