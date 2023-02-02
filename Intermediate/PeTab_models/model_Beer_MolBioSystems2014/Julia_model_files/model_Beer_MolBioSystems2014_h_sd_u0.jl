#u[1] = Glu, u[2] = cGlu, u[3] = Ind, u[4] = Bac
#θ_dynamicNames[1] = Bacmax_typeIDT1_ExpID1, θ_dynamicNames[2] = Bacmax_typeIDT1_ExpID2, θ_dynamicNames[3] = Bacmax_typeIDT1_ExpID3, θ_dynamicNames[4] = Bacmax_typeIDT1_ExpID4, θ_dynamicNames[5] = Bacmax_typeIDT1_ExpID5, θ_dynamicNames[6] = Bacmax_typeIDT1_ExpID6, θ_dynamicNames[7] = Bacmax_typeIDT3_ExpID1, θ_dynamicNames[8] = Bacmax_typeIDT3_ExpID2, θ_dynamicNames[9] = Bacmax_typeIDT3_ExpID3, θ_dynamicNames[10] = Bacmax_typeIDT3_ExpID4, θ_dynamicNames[11] = Bacmax_typeIDT3_ExpID5, θ_dynamicNames[12] = Bacmax_typeIDT3_ExpID6, θ_dynamicNames[13] = Bacmax_typeIDT5_ExpID1, θ_dynamicNames[14] = Bacmax_typeIDT5_ExpID2, θ_dynamicNames[15] = Bacmax_typeIDT5_ExpID3, θ_dynamicNames[16] = Bacmax_typeIDT5_ExpID4, θ_dynamicNames[17] = Bacmax_typeIDT5_ExpID5, θ_dynamicNames[18] = Bacmax_typeIDT5_ExpID6, θ_dynamicNames[19] = Bacmax_typeIDwt_ExpID4, θ_dynamicNames[20] = beta_typeIDT1_ExpID1, θ_dynamicNames[21] = beta_typeIDT1_ExpID2, θ_dynamicNames[22] = beta_typeIDT1_ExpID3, θ_dynamicNames[23] = beta_typeIDT1_ExpID4, θ_dynamicNames[24] = beta_typeIDT1_ExpID5, θ_dynamicNames[25] = beta_typeIDT1_ExpID6, θ_dynamicNames[26] = beta_typeIDT3_ExpID1, θ_dynamicNames[27] = beta_typeIDT3_ExpID2, θ_dynamicNames[28] = beta_typeIDT3_ExpID3, θ_dynamicNames[29] = beta_typeIDT3_ExpID4, θ_dynamicNames[30] = beta_typeIDT3_ExpID5, θ_dynamicNames[31] = beta_typeIDT3_ExpID6, θ_dynamicNames[32] = beta_typeIDT5_ExpID1, θ_dynamicNames[33] = beta_typeIDT5_ExpID2, θ_dynamicNames[34] = beta_typeIDT5_ExpID3, θ_dynamicNames[35] = beta_typeIDT5_ExpID4, θ_dynamicNames[36] = beta_typeIDT5_ExpID5, θ_dynamicNames[37] = beta_typeIDT5_ExpID6, θ_dynamicNames[38] = beta_typeIDwt_ExpID4, θ_dynamicNames[39] = init_Bac, θ_dynamicNames[40] = kdegi_typeIDT1, θ_dynamicNames[41] = kdegi_typeIDT3, θ_dynamicNames[42] = kdegi_typeIDT5, θ_dynamicNames[43] = kdegi_typeIDwt, θ_dynamicNames[44] = kdim_typeIDT1, θ_dynamicNames[45] = kdim_typeIDT3, θ_dynamicNames[46] = kdim_typeIDT5, θ_dynamicNames[47] = kdim_typeIDwt, θ_dynamicNames[48] = ksyn_typeIDT1, θ_dynamicNames[49] = ksyn_typeIDT3, θ_dynamicNames[50] = ksyn_typeIDT5, θ_dynamicNames[51] = ksyn_typeIDwt, θ_dynamicNames[52] = tau_typeIDT1_ExpID1, θ_dynamicNames[53] = tau_typeIDT1_ExpID2, θ_dynamicNames[54] = tau_typeIDT1_ExpID3, θ_dynamicNames[55] = tau_typeIDT1_ExpID4, θ_dynamicNames[56] = tau_typeIDT1_ExpID5, θ_dynamicNames[57] = tau_typeIDT1_ExpID6, θ_dynamicNames[58] = tau_typeIDT3_ExpID1, θ_dynamicNames[59] = tau_typeIDT3_ExpID2, θ_dynamicNames[60] = tau_typeIDT3_ExpID3, θ_dynamicNames[61] = tau_typeIDT3_ExpID4, θ_dynamicNames[62] = tau_typeIDT3_ExpID5, θ_dynamicNames[63] = tau_typeIDT3_ExpID6, θ_dynamicNames[64] = tau_typeIDT5_ExpID1, θ_dynamicNames[65] = tau_typeIDT5_ExpID2, θ_dynamicNames[66] = tau_typeIDT5_ExpID3, θ_dynamicNames[67] = tau_typeIDT5_ExpID4, θ_dynamicNames[68] = tau_typeIDT5_ExpID5, θ_dynamicNames[69] = tau_typeIDT5_ExpID6, θ_dynamicNames[70] = tau_typeIDwt_ExpID4
#

function compute_h(u::AbstractVector, t::Real, θ_dynamic::AbstractVector, θ_observable::AbstractVector,
                   θ_nonDynamic::AbstractVector, parameterInfo::ParametersInfo, observableId::Symbol, 
                      parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :Bacnorm 
		return u[4] 
	end

	if observableId == :IndconcNormRange 
		return u[3] 
	end

end

function compute_u0!(u0::AbstractVector, pODEProblem::AbstractVector) 

	#pODEProblem[1] = lag_bool1, pODEProblem[2] = kdegi, pODEProblem[3] = medium, pODEProblem[4] = Bacmax, pODEProblem[5] = ksyn, pODEProblem[6] = kdim, pODEProblem[7] = tau, pODEProblem[8] = init_Bac, pODEProblem[9] = beta

	Glu = 10.0 
	cGlu = 0.0 
	Ind = 0.0 
	Bac = pODEProblem[8] 

	u0 .= Glu, cGlu, Ind, Bac
end

function compute_u0(pODEProblem::AbstractVector)::AbstractVector 

	#pODEProblem[1] = lag_bool1, pODEProblem[2] = kdegi, pODEProblem[3] = medium, pODEProblem[4] = Bacmax, pODEProblem[5] = ksyn, pODEProblem[6] = kdim, pODEProblem[7] = tau, pODEProblem[8] = init_Bac, pODEProblem[9] = beta

	Glu = 10.0 
	cGlu = 0.0 
	Ind = 0.0 
	Bac = pODEProblem[8] 

	 return [Glu, cGlu, Ind, Bac]
end

function compute_σ(u::AbstractVector, t::Real, θ_sd::AbstractVector, θ_dynamic::AbstractVector, θ_nonDynamic::AbstractVector, 
                   parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :Bacnorm 
		noiseParameter1_Bacnorm = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_Bacnorm 
	end

	if observableId == :IndconcNormRange 
		noiseParameter1_IndconcNormRange = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_IndconcNormRange 
	end

end