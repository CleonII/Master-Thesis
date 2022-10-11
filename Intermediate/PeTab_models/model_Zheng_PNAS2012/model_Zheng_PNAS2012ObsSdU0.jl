function evalYmod(u, t, dynPar, obsPar, paramData, observableId, mapObsParam) 

	K27me0K36me0, K27me2K36me3, K27me2K36me0, K27me0K36me1, K27me2K36me1, K27me0K36me2, K27me1K36me2, K27me3K36me2, K27me1K36me1, K27me1K36me3, K27me2K36me2, K27me0K36me3, K27me3K36me1, K27me3K36me0, K27me1K36me0= u 
	inflowp, k00_01, k00_10, k01_00, k01_02, k01_11, k02_01, k02_03, k02_12, k03_02, k03_13, k10_00, k10_11, k10_20, k11_01, k11_10, k11_12, k11_21, k12_02, k12_11, k12_13, k12_22, k13_03, k13_12, k13_23, k20_10, k20_21, k20_30, k21_11, k21_20, k21_22, k21_31, k22_12, k22_21, k22_23, k22_32, k23_13, k23_22, k30_20, k30_31, k31_21, k31_30, k31_32, k32_22, k32_31 = dynPar 

	if observableId == "observable_K27me0K36me0" 
		return K27me0K36me0 / ( K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2 ) 
	end

	if observableId == "observable_K27me0K36me1" 
		return K27me0K36me1 / ( K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2 ) 
	end

	if observableId == "observable_K27me0K36me2" 
		return K27me0K36me2 / ( K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2 ) 
	end

	if observableId == "observable_K27me0K36me3" 
		return K27me0K36me3 / ( K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2 ) 
	end

	if observableId == "observable_K27me1K36me0" 
		return K27me1K36me0 / ( K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2 ) 
	end

	if observableId == "observable_K27me1K36me1" 
		return K27me1K36me1 / ( K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2 ) 
	end

	if observableId == "observable_K27me1K36me2" 
		return K27me1K36me2 / ( K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2 ) 
	end

	if observableId == "observable_K27me1K36me3" 
		return K27me1K36me3 / ( K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2 ) 
	end

	if observableId == "observable_K27me2K36me0" 
		return K27me2K36me0 / ( K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2 ) 
	end

	if observableId == "observable_K27me2K36me1" 
		return K27me2K36me1 / ( K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2 ) 
	end

	if observableId == "observable_K27me2K36me2" 
		return K27me2K36me2 / ( K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2 ) 
	end

	if observableId == "observable_K27me2K36me3" 
		return K27me2K36me3 / ( K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2 ) 
	end

	if observableId == "observable_K27me3K36me0" 
		return K27me3K36me0 / ( K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2 ) 
	end

	if observableId == "observable_K27me3K36me1" 
		return K27me3K36me1 / ( K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2 ) 
	end

	if observableId == "observable_K27me3K36me2" 
		return K27me3K36me2 / ( K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2 ) 
	end

end

function evalU0!(u0Vec, paramVec) 

	inflowp, k00_01, k10_00, default, dilution, k01_00, k00_10, k13_23, k22_23, k23_13, k23_22, k10_20, k20_30, k20_10, k21_20, k30_20, k20_21, k02_01, k01_11, k11_01, k01_02, k21_31, k11_21, k21_22, k21_11, k31_21, k22_21, k02_12, k03_02, k02_03, k12_02, k22_12, k12_22, k11_12, k13_12, k12_11, k12_13, k32_22, k22_32, k31_32, k32_31, k10_11, k11_10, k13_03, k03_13, k31_30, k30_31 = paramVec 

	K27me0K36me0 = 0.00417724976345759 
	K27me2K36me3 = 0.00471831436002134 
	K27me2K36me0 = 0.00632744816295157 
	K27me0K36me1 = 0.0102104668587641 
	K27me2K36me1 = 0.0143896310177379 
	K27me0K36me2 = 0.169690316239546 
	K27me1K36me2 = 0.594249755169037 
	K27me3K36me2 = 0.00136041631795562 
	K27me1K36me1 = 0.0078328187288069 
	K27me1K36me3 = 0.102748675077958 
	K27me2K36me2 = 0.0263372634996529 
	K27me0K36me3 = 0.0504935214807544 
	K27me3K36me1 = 0.00250831034920277 
	K27me3K36me0 = 0.00330168411604165 
	K27me1K36me0 = 0.00165412810279407 

	u0Vec .= K27me0K36me0, K27me2K36me3, K27me2K36me0, K27me0K36me1, K27me2K36me1, K27me0K36me2, K27me1K36me2, K27me3K36me2, K27me1K36me1, K27me1K36me3, K27me2K36me2, K27me0K36me3, K27me3K36me1, K27me3K36me0, K27me1K36me0
end

function evalSd!(u, t, sdPar, dynPar, paramData, observableId, mapSdParam) 

	K27me0K36me0, K27me2K36me3, K27me2K36me0, K27me0K36me1, K27me2K36me1, K27me0K36me2, K27me1K36me2, K27me3K36me2, K27me1K36me1, K27me1K36me3, K27me2K36me2, K27me0K36me3, K27me3K36me1, K27me3K36me0, K27me1K36me0= u 
	inflowp, k00_01, k00_10, k01_00, k01_02, k01_11, k02_01, k02_03, k02_12, k03_02, k03_13, k10_00, k10_11, k10_20, k11_01, k11_10, k11_12, k11_21, k12_02, k12_11, k12_13, k12_22, k13_03, k13_12, k13_23, k20_10, k20_21, k20_30, k21_11, k21_20, k21_22, k21_31, k22_12, k22_21, k22_23, k22_32, k23_13, k23_22, k30_20, k30_31, k31_21, k31_30, k31_32, k32_22, k32_31 = dynPar 

	if observableId == "observable_K27me0K36me0" 
		noiseParameter1_observable_K27me0K36me0 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_K27me0K36me0 
	end

	if observableId == "observable_K27me0K36me1" 
		noiseParameter1_observable_K27me0K36me1 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_K27me0K36me1 
	end

	if observableId == "observable_K27me0K36me2" 
		noiseParameter1_observable_K27me0K36me2 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_K27me0K36me2 
	end

	if observableId == "observable_K27me0K36me3" 
		noiseParameter1_observable_K27me0K36me3 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_K27me0K36me3 
	end

	if observableId == "observable_K27me1K36me0" 
		noiseParameter1_observable_K27me1K36me0 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_K27me1K36me0 
	end

	if observableId == "observable_K27me1K36me1" 
		noiseParameter1_observable_K27me1K36me1 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_K27me1K36me1 
	end

	if observableId == "observable_K27me1K36me2" 
		noiseParameter1_observable_K27me1K36me2 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_K27me1K36me2 
	end

	if observableId == "observable_K27me1K36me3" 
		noiseParameter1_observable_K27me1K36me3 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_K27me1K36me3 
	end

	if observableId == "observable_K27me2K36me0" 
		noiseParameter1_observable_K27me2K36me0 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_K27me2K36me0 
	end

	if observableId == "observable_K27me2K36me1" 
		noiseParameter1_observable_K27me2K36me1 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_K27me2K36me1 
	end

	if observableId == "observable_K27me2K36me2" 
		noiseParameter1_observable_K27me2K36me2 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_K27me2K36me2 
	end

	if observableId == "observable_K27me2K36me3" 
		noiseParameter1_observable_K27me2K36me3 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_K27me2K36me3 
	end

	if observableId == "observable_K27me3K36me0" 
		noiseParameter1_observable_K27me3K36me0 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_K27me3K36me0 
	end

	if observableId == "observable_K27me3K36me1" 
		noiseParameter1_observable_K27me3K36me1 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_K27me3K36me1 
	end

	if observableId == "observable_K27me3K36me2" 
		noiseParameter1_observable_K27me3K36me2 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_K27me3K36me2 
	end

end