function evalDYmodDu(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	X_protein, GFP_mRNA, Y_mRNA, X_mRNA, Z_mRNA, Z_protein, Y_protein, GFP= u 
	tau_mRNA, tps_repr, init_GFP, n_Hill, init_Y_mRNA, init_Z_mRNA, init_X_protein, init_Y_protein, tau_mRNA_GFP, eff, init_Z_protein, tau_prot_GFP, eff_GFP, tps_active, init_GFP_mRNA, KM, init_X_mRNA, tau_prot, cell = p 
	if observableId == :fluorescence 
		observableParameter1_fluorescence, observableParameter2_fluorescence = getObsOrSdParam(obsPar, mapObsParam)
		out[8] = observableParameter2_fluorescence
		return nothing
	end

end

function evalDYmodDp(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	X_protein, GFP_mRNA, Y_mRNA, X_mRNA, Z_mRNA, Z_protein, Y_protein, GFP= u 
	tau_mRNA, tps_repr, init_GFP, n_Hill, init_Y_mRNA, init_Z_mRNA, init_X_protein, init_Y_protein, tau_mRNA_GFP, eff, init_Z_protein, tau_prot_GFP, eff_GFP, tps_active, init_GFP_mRNA, KM, init_X_mRNA, tau_prot, cell = p 
	if observableId == :fluorescence 
		return nothing
	end

end

function evalDSdDu!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	X_protein, GFP_mRNA, Y_mRNA, X_mRNA, Z_mRNA, Z_protein, Y_protein, GFP= u 
	tau_mRNA, tps_repr, init_GFP, n_Hill, init_Y_mRNA, init_Z_mRNA, init_X_protein, init_Y_protein, tau_mRNA_GFP, eff, init_Z_protein, tau_prot_GFP, eff_GFP, tps_active, init_GFP_mRNA, KM, init_X_mRNA, tau_prot, cell = p 
	if observableId == "fluorescence" 
		return nothing
	end

end

function evalDSdDp!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	X_protein, GFP_mRNA, Y_mRNA, X_mRNA, Z_mRNA, Z_protein, Y_protein, GFP= u 
	tau_mRNA, tps_repr, init_GFP, n_Hill, init_Y_mRNA, init_Z_mRNA, init_X_protein, init_Y_protein, tau_mRNA_GFP, eff, init_Z_protein, tau_prot_GFP, eff_GFP, tps_active, init_GFP_mRNA, KM, init_X_mRNA, tau_prot, cell = p 
	if observableId == "fluorescence" 
		return nothing
	end

end

