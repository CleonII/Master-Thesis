#u[1] = K27me0K36me0, u[2] = K27me2K36me3, u[3] = K27me2K36me0, u[4] = K27me0K36me1, u[5] = K27me2K36me1, u[6] = K27me0K36me2, u[7] = K27me1K36me2, u[8] = K27me3K36me2, u[9] = K27me1K36me1, u[10] = K27me1K36me3, u[11] = K27me2K36me2, u[12] = K27me0K36me3, u[13] = K27me3K36me1, u[14] = K27me3K36me0, u[15] = K27me1K36me0
#θ_dynamicNames[1] = inflowp, θ_dynamicNames[2] = k00_01, θ_dynamicNames[3] = k00_10, θ_dynamicNames[4] = k01_00, θ_dynamicNames[5] = k01_02, θ_dynamicNames[6] = k01_11, θ_dynamicNames[7] = k02_01, θ_dynamicNames[8] = k02_03, θ_dynamicNames[9] = k02_12, θ_dynamicNames[10] = k03_02, θ_dynamicNames[11] = k03_13, θ_dynamicNames[12] = k10_00, θ_dynamicNames[13] = k10_11, θ_dynamicNames[14] = k10_20, θ_dynamicNames[15] = k11_01, θ_dynamicNames[16] = k11_10, θ_dynamicNames[17] = k11_12, θ_dynamicNames[18] = k11_21, θ_dynamicNames[19] = k12_02, θ_dynamicNames[20] = k12_11, θ_dynamicNames[21] = k12_13, θ_dynamicNames[22] = k12_22, θ_dynamicNames[23] = k13_03, θ_dynamicNames[24] = k13_12, θ_dynamicNames[25] = k13_23, θ_dynamicNames[26] = k20_10, θ_dynamicNames[27] = k20_21, θ_dynamicNames[28] = k20_30, θ_dynamicNames[29] = k21_11, θ_dynamicNames[30] = k21_20, θ_dynamicNames[31] = k21_22, θ_dynamicNames[32] = k21_31, θ_dynamicNames[33] = k22_12, θ_dynamicNames[34] = k22_21, θ_dynamicNames[35] = k22_23, θ_dynamicNames[36] = k22_32, θ_dynamicNames[37] = k23_13, θ_dynamicNames[38] = k23_22, θ_dynamicNames[39] = k30_20, θ_dynamicNames[40] = k30_31, θ_dynamicNames[41] = k31_21, θ_dynamicNames[42] = k31_30, θ_dynamicNames[43] = k31_32, θ_dynamicNames[44] = k32_22, θ_dynamicNames[45] = k32_31
#

function compute_h(u::AbstractVector, t::Real, θ_dynamic::AbstractVector, θ_observable::AbstractVector,
                   θ_nonDynamic::AbstractVector, parameterInfo::ParametersInfo, observableId::Symbol, 
                      parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :observable_K27me0K36me0 
		return u[1] / ( u[1] + u[4] + u[6] + u[12] + u[15] + u[9] + u[7] + u[10] + u[3] + u[5] + u[11] + u[2] + u[14] + u[13] + u[8] ) 
	end

	if observableId == :observable_K27me0K36me1 
		return u[4] / ( u[1] + u[4] + u[6] + u[12] + u[15] + u[9] + u[7] + u[10] + u[3] + u[5] + u[11] + u[2] + u[14] + u[13] + u[8] ) 
	end

	if observableId == :observable_K27me0K36me2 
		return u[6] / ( u[1] + u[4] + u[6] + u[12] + u[15] + u[9] + u[7] + u[10] + u[3] + u[5] + u[11] + u[2] + u[14] + u[13] + u[8] ) 
	end

	if observableId == :observable_K27me0K36me3 
		return u[12] / ( u[1] + u[4] + u[6] + u[12] + u[15] + u[9] + u[7] + u[10] + u[3] + u[5] + u[11] + u[2] + u[14] + u[13] + u[8] ) 
	end

	if observableId == :observable_K27me1K36me0 
		return u[15] / ( u[1] + u[4] + u[6] + u[12] + u[15] + u[9] + u[7] + u[10] + u[3] + u[5] + u[11] + u[2] + u[14] + u[13] + u[8] ) 
	end

	if observableId == :observable_K27me1K36me1 
		return u[9] / ( u[1] + u[4] + u[6] + u[12] + u[15] + u[9] + u[7] + u[10] + u[3] + u[5] + u[11] + u[2] + u[14] + u[13] + u[8] ) 
	end

	if observableId == :observable_K27me1K36me2 
		return u[7] / ( u[1] + u[4] + u[6] + u[12] + u[15] + u[9] + u[7] + u[10] + u[3] + u[5] + u[11] + u[2] + u[14] + u[13] + u[8] ) 
	end

	if observableId == :observable_K27me1K36me3 
		return u[10] / ( u[1] + u[4] + u[6] + u[12] + u[15] + u[9] + u[7] + u[10] + u[3] + u[5] + u[11] + u[2] + u[14] + u[13] + u[8] ) 
	end

	if observableId == :observable_K27me2K36me0 
		return u[3] / ( u[1] + u[4] + u[6] + u[12] + u[15] + u[9] + u[7] + u[10] + u[3] + u[5] + u[11] + u[2] + u[14] + u[13] + u[8] ) 
	end

	if observableId == :observable_K27me2K36me1 
		return u[5] / ( u[1] + u[4] + u[6] + u[12] + u[15] + u[9] + u[7] + u[10] + u[3] + u[5] + u[11] + u[2] + u[14] + u[13] + u[8] ) 
	end

	if observableId == :observable_K27me2K36me2 
		return u[11] / ( u[1] + u[4] + u[6] + u[12] + u[15] + u[9] + u[7] + u[10] + u[3] + u[5] + u[11] + u[2] + u[14] + u[13] + u[8] ) 
	end

	if observableId == :observable_K27me2K36me3 
		return u[2] / ( u[1] + u[4] + u[6] + u[12] + u[15] + u[9] + u[7] + u[10] + u[3] + u[5] + u[11] + u[2] + u[14] + u[13] + u[8] ) 
	end

	if observableId == :observable_K27me3K36me0 
		return u[14] / ( u[1] + u[4] + u[6] + u[12] + u[15] + u[9] + u[7] + u[10] + u[3] + u[5] + u[11] + u[2] + u[14] + u[13] + u[8] ) 
	end

	if observableId == :observable_K27me3K36me1 
		return u[13] / ( u[1] + u[4] + u[6] + u[12] + u[15] + u[9] + u[7] + u[10] + u[3] + u[5] + u[11] + u[2] + u[14] + u[13] + u[8] ) 
	end

	if observableId == :observable_K27me3K36me2 
		return u[8] / ( u[1] + u[4] + u[6] + u[12] + u[15] + u[9] + u[7] + u[10] + u[3] + u[5] + u[11] + u[2] + u[14] + u[13] + u[8] ) 
	end

end

function compute_u0!(u0::AbstractVector, pODEProblem::AbstractVector) 

	#pODEProblem[1] = k20_10, pODEProblem[2] = k11_12, pODEProblem[3] = k22_32, pODEProblem[4] = k21_11, pODEProblem[5] = k13_12, pODEProblem[6] = k01_00, pODEProblem[7] = k11_01, pODEProblem[8] = k22_21, pODEProblem[9] = k22_12, pODEProblem[10] = k11_10, pODEProblem[11] = k12_11, pODEProblem[12] = k21_31, pODEProblem[13] = k02_12, pODEProblem[14] = k30_20, pODEProblem[15] = dilution, pODEProblem[16] = k31_21, pODEProblem[17] = k23_22, pODEProblem[18] = k02_01, pODEProblem[19] = k03_13, pODEProblem[20] = k10_11, pODEProblem[21] = k21_20, pODEProblem[22] = k10_00, pODEProblem[23] = k30_31, pODEProblem[24] = k20_30, pODEProblem[25] = k20_21, pODEProblem[26] = k00_10, pODEProblem[27] = k12_13, pODEProblem[28] = k01_11, pODEProblem[29] = k02_03, pODEProblem[30] = k00_01, pODEProblem[31] = k03_02, pODEProblem[32] = k23_13, pODEProblem[33] = k32_31, pODEProblem[34] = default, pODEProblem[35] = k13_03, pODEProblem[36] = k31_30, pODEProblem[37] = k12_02, pODEProblem[38] = k31_32, pODEProblem[39] = k01_02, pODEProblem[40] = k13_23, pODEProblem[41] = k21_22, pODEProblem[42] = k12_22, pODEProblem[43] = k11_21, pODEProblem[44] = k22_23, pODEProblem[45] = k10_20, pODEProblem[46] = inflowp, pODEProblem[47] = k32_22

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

	u0 .= K27me0K36me0, K27me2K36me3, K27me2K36me0, K27me0K36me1, K27me2K36me1, K27me0K36me2, K27me1K36me2, K27me3K36me2, K27me1K36me1, K27me1K36me3, K27me2K36me2, K27me0K36me3, K27me3K36me1, K27me3K36me0, K27me1K36me0
end

function compute_u0(pODEProblem::AbstractVector)::AbstractVector 

	#pODEProblem[1] = k20_10, pODEProblem[2] = k11_12, pODEProblem[3] = k22_32, pODEProblem[4] = k21_11, pODEProblem[5] = k13_12, pODEProblem[6] = k01_00, pODEProblem[7] = k11_01, pODEProblem[8] = k22_21, pODEProblem[9] = k22_12, pODEProblem[10] = k11_10, pODEProblem[11] = k12_11, pODEProblem[12] = k21_31, pODEProblem[13] = k02_12, pODEProblem[14] = k30_20, pODEProblem[15] = dilution, pODEProblem[16] = k31_21, pODEProblem[17] = k23_22, pODEProblem[18] = k02_01, pODEProblem[19] = k03_13, pODEProblem[20] = k10_11, pODEProblem[21] = k21_20, pODEProblem[22] = k10_00, pODEProblem[23] = k30_31, pODEProblem[24] = k20_30, pODEProblem[25] = k20_21, pODEProblem[26] = k00_10, pODEProblem[27] = k12_13, pODEProblem[28] = k01_11, pODEProblem[29] = k02_03, pODEProblem[30] = k00_01, pODEProblem[31] = k03_02, pODEProblem[32] = k23_13, pODEProblem[33] = k32_31, pODEProblem[34] = default, pODEProblem[35] = k13_03, pODEProblem[36] = k31_30, pODEProblem[37] = k12_02, pODEProblem[38] = k31_32, pODEProblem[39] = k01_02, pODEProblem[40] = k13_23, pODEProblem[41] = k21_22, pODEProblem[42] = k12_22, pODEProblem[43] = k11_21, pODEProblem[44] = k22_23, pODEProblem[45] = k10_20, pODEProblem[46] = inflowp, pODEProblem[47] = k32_22

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

	 return [K27me0K36me0, K27me2K36me3, K27me2K36me0, K27me0K36me1, K27me2K36me1, K27me0K36me2, K27me1K36me2, K27me3K36me2, K27me1K36me1, K27me1K36me3, K27me2K36me2, K27me0K36me3, K27me3K36me1, K27me3K36me0, K27me1K36me0]
end

function compute_σ(u::AbstractVector, t::Real, θ_sd::AbstractVector, θ_dynamic::AbstractVector, θ_nonDynamic::AbstractVector, 
                   parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :observable_K27me0K36me0 
		noiseParameter1_observable_K27me0K36me0 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_K27me0K36me0 
	end

	if observableId == :observable_K27me0K36me1 
		noiseParameter1_observable_K27me0K36me1 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_K27me0K36me1 
	end

	if observableId == :observable_K27me0K36me2 
		noiseParameter1_observable_K27me0K36me2 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_K27me0K36me2 
	end

	if observableId == :observable_K27me0K36me3 
		noiseParameter1_observable_K27me0K36me3 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_K27me0K36me3 
	end

	if observableId == :observable_K27me1K36me0 
		noiseParameter1_observable_K27me1K36me0 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_K27me1K36me0 
	end

	if observableId == :observable_K27me1K36me1 
		noiseParameter1_observable_K27me1K36me1 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_K27me1K36me1 
	end

	if observableId == :observable_K27me1K36me2 
		noiseParameter1_observable_K27me1K36me2 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_K27me1K36me2 
	end

	if observableId == :observable_K27me1K36me3 
		noiseParameter1_observable_K27me1K36me3 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_K27me1K36me3 
	end

	if observableId == :observable_K27me2K36me0 
		noiseParameter1_observable_K27me2K36me0 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_K27me2K36me0 
	end

	if observableId == :observable_K27me2K36me1 
		noiseParameter1_observable_K27me2K36me1 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_K27me2K36me1 
	end

	if observableId == :observable_K27me2K36me2 
		noiseParameter1_observable_K27me2K36me2 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_K27me2K36me2 
	end

	if observableId == :observable_K27me2K36me3 
		noiseParameter1_observable_K27me2K36me3 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_K27me2K36me3 
	end

	if observableId == :observable_K27me3K36me0 
		noiseParameter1_observable_K27me3K36me0 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_K27me3K36me0 
	end

	if observableId == :observable_K27me3K36me1 
		noiseParameter1_observable_K27me3K36me1 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_K27me3K36me1 
	end

	if observableId == :observable_K27me3K36me2 
		noiseParameter1_observable_K27me3K36me2 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_K27me3K36me2 
	end

end