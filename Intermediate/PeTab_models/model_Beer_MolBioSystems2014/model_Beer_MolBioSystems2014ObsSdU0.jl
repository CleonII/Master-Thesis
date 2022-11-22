function evalYmod(u, t, dynPar, obsPar, nonDynParam, paramData, observableId, mapObsParam) 

	Glu, cGlu, Ind, Bac= u 
	Bacmax_typeIDT1_ExpID1, Bacmax_typeIDT1_ExpID2, Bacmax_typeIDT1_ExpID3, Bacmax_typeIDT1_ExpID4, Bacmax_typeIDT1_ExpID5, Bacmax_typeIDT1_ExpID6, Bacmax_typeIDT3_ExpID1, Bacmax_typeIDT3_ExpID2, Bacmax_typeIDT3_ExpID3, Bacmax_typeIDT3_ExpID4, Bacmax_typeIDT3_ExpID5, Bacmax_typeIDT3_ExpID6, Bacmax_typeIDT5_ExpID1, Bacmax_typeIDT5_ExpID2, Bacmax_typeIDT5_ExpID3, Bacmax_typeIDT5_ExpID4, Bacmax_typeIDT5_ExpID5, Bacmax_typeIDT5_ExpID6, Bacmax_typeIDwt_ExpID4, beta_typeIDT1_ExpID1, beta_typeIDT1_ExpID2, beta_typeIDT1_ExpID3, beta_typeIDT1_ExpID4, beta_typeIDT1_ExpID5, beta_typeIDT1_ExpID6, beta_typeIDT3_ExpID1, beta_typeIDT3_ExpID2, beta_typeIDT3_ExpID3, beta_typeIDT3_ExpID4, beta_typeIDT3_ExpID5, beta_typeIDT3_ExpID6, beta_typeIDT5_ExpID1, beta_typeIDT5_ExpID2, beta_typeIDT5_ExpID3, beta_typeIDT5_ExpID4, beta_typeIDT5_ExpID5, beta_typeIDT5_ExpID6, beta_typeIDwt_ExpID4, init_Bac, kdegi_typeIDT1, kdegi_typeIDT3, kdegi_typeIDT5, kdegi_typeIDwt, kdim_typeIDT1, kdim_typeIDT3, kdim_typeIDT5, kdim_typeIDwt, ksyn_typeIDT1, ksyn_typeIDT3, ksyn_typeIDT5, ksyn_typeIDwt, tau_typeIDT1_ExpID1, tau_typeIDT1_ExpID2, tau_typeIDT1_ExpID3, tau_typeIDT1_ExpID4, tau_typeIDT1_ExpID5, tau_typeIDT1_ExpID6, tau_typeIDT3_ExpID1, tau_typeIDT3_ExpID2, tau_typeIDT3_ExpID3, tau_typeIDT3_ExpID4, tau_typeIDT3_ExpID5, tau_typeIDT3_ExpID6, tau_typeIDT5_ExpID1, tau_typeIDT5_ExpID2, tau_typeIDT5_ExpID3, tau_typeIDT5_ExpID4, tau_typeIDT5_ExpID5, tau_typeIDT5_ExpID6, tau_typeIDwt_ExpID4 = dynPar 

	if observableId == "Bacnorm" 
		return Bac 
	end

	if observableId == "IndconcNormRange" 
		return Ind 
	end

end

function evalU0!(u0Vec, paramVec) 

	kdegi, medium, Bacmax, ksyn, kdim, tau, init_Bac, beta = paramVec 

	Glu = 10.0 
	cGlu = 0.0 
	Ind = 0.0 
	Bac = init_Bac 

	u0Vec .= Glu, cGlu, Ind, Bac
end

function evalU0(paramVec) 

	kdegi, medium, Bacmax, ksyn, kdim, tau, init_Bac, beta = paramVec 

	Glu = 10.0 
	cGlu = 0.0 
	Ind = 0.0 
	Bac = init_Bac 

	 return [Glu, cGlu, Ind, Bac]
end

function evalSd!(u, t, sdPar, dynPar, nonDynParam, paramData, observableId, mapSdParam) 

	Glu, cGlu, Ind, Bac= u 
	Bacmax_typeIDT1_ExpID1, Bacmax_typeIDT1_ExpID2, Bacmax_typeIDT1_ExpID3, Bacmax_typeIDT1_ExpID4, Bacmax_typeIDT1_ExpID5, Bacmax_typeIDT1_ExpID6, Bacmax_typeIDT3_ExpID1, Bacmax_typeIDT3_ExpID2, Bacmax_typeIDT3_ExpID3, Bacmax_typeIDT3_ExpID4, Bacmax_typeIDT3_ExpID5, Bacmax_typeIDT3_ExpID6, Bacmax_typeIDT5_ExpID1, Bacmax_typeIDT5_ExpID2, Bacmax_typeIDT5_ExpID3, Bacmax_typeIDT5_ExpID4, Bacmax_typeIDT5_ExpID5, Bacmax_typeIDT5_ExpID6, Bacmax_typeIDwt_ExpID4, beta_typeIDT1_ExpID1, beta_typeIDT1_ExpID2, beta_typeIDT1_ExpID3, beta_typeIDT1_ExpID4, beta_typeIDT1_ExpID5, beta_typeIDT1_ExpID6, beta_typeIDT3_ExpID1, beta_typeIDT3_ExpID2, beta_typeIDT3_ExpID3, beta_typeIDT3_ExpID4, beta_typeIDT3_ExpID5, beta_typeIDT3_ExpID6, beta_typeIDT5_ExpID1, beta_typeIDT5_ExpID2, beta_typeIDT5_ExpID3, beta_typeIDT5_ExpID4, beta_typeIDT5_ExpID5, beta_typeIDT5_ExpID6, beta_typeIDwt_ExpID4, init_Bac, kdegi_typeIDT1, kdegi_typeIDT3, kdegi_typeIDT5, kdegi_typeIDwt, kdim_typeIDT1, kdim_typeIDT3, kdim_typeIDT5, kdim_typeIDwt, ksyn_typeIDT1, ksyn_typeIDT3, ksyn_typeIDT5, ksyn_typeIDwt, tau_typeIDT1_ExpID1, tau_typeIDT1_ExpID2, tau_typeIDT1_ExpID3, tau_typeIDT1_ExpID4, tau_typeIDT1_ExpID5, tau_typeIDT1_ExpID6, tau_typeIDT3_ExpID1, tau_typeIDT3_ExpID2, tau_typeIDT3_ExpID3, tau_typeIDT3_ExpID4, tau_typeIDT3_ExpID5, tau_typeIDT3_ExpID6, tau_typeIDT5_ExpID1, tau_typeIDT5_ExpID2, tau_typeIDT5_ExpID3, tau_typeIDT5_ExpID4, tau_typeIDT5_ExpID5, tau_typeIDT5_ExpID6, tau_typeIDwt_ExpID4 = dynPar 

	if observableId == "Bacnorm" 
		noiseParameter1_Bacnorm = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_Bacnorm 
	end

	if observableId == "IndconcNormRange" 
		noiseParameter1_IndconcNormRange = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_IndconcNormRange 
	end

end