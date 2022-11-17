function evalYmod(u, t, dynPar, obsPar, nonDynParam, paramData, observableId, mapObsParam) 

	X_protein, GFP_mRNA, Y_mRNA, X_mRNA, Z_mRNA, Z_protein, Y_protein, GFP, dummyVariable= u 
	KM, eff, eff_GFP, init_GFP, init_GFP_mRNA, init_X_mRNA, init_X_protein, init_Y_mRNA, init_Y_protein, init_Z_mRNA, init_Z_protein, n_Hill, tau_mRNA, tau_mRNA_GFP, tau_prot, tau_prot_GFP, tps_active, tps_repr = dynPar 

	if observableId == "fluorescence" 
		observableParameter1_fluorescence, observableParameter2_fluorescence = getObsOrSdParam(obsPar, mapObsParam)
		return observableParameter1_fluorescence + GFP * observableParameter2_fluorescence 
	end

end

function evalU0!(u0Vec, paramVec) 

	tau_mRNA, tau_prot, eff, tau_mRNA_GFP, KM, tps_repr, n_Hill, tps_active, tau_prot_GFP, eff_GFP, init_X_protein, init_GFP_mRNA, init_GFP, init_Y_mRNA, init_Z_protein, init_Z_mRNA, init_X_mRNA, init_Y_protein = paramVec 

	X_protein = init_X_protein 
	GFP_mRNA = init_GFP_mRNA 
	Y_mRNA = init_Y_mRNA 
	X_mRNA = init_X_mRNA 
	Z_mRNA = init_Z_mRNA 
	Z_protein = init_Z_protein 
	Y_protein = init_Y_protein 
	GFP = init_GFP 
	dummyVariable = 0.0 

	u0Vec .= X_protein, GFP_mRNA, Y_mRNA, X_mRNA, Z_mRNA, Z_protein, Y_protein, GFP, dummyVariable
end

function evalSd!(u, t, sdPar, dynPar, nonDynPar, paramData, observableId, mapSdParam) 

	X_protein, GFP_mRNA, Y_mRNA, X_mRNA, Z_mRNA, Z_protein, Y_protein, GFP, dummyVariable= u 
	KM, eff, eff_GFP, init_GFP, init_GFP_mRNA, init_X_mRNA, init_X_protein, init_Y_mRNA, init_Y_protein, init_Z_mRNA, init_Z_protein, n_Hill, tau_mRNA, tau_mRNA_GFP, tau_prot, tau_prot_GFP, tps_active, tps_repr = dynPar 

	if observableId == "fluorescence" 
		noiseParameter1_fluorescence = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_fluorescence 
	end

end