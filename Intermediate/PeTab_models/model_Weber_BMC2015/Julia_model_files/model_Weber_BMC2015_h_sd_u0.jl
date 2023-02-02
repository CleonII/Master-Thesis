#u[1] = CERTERa, u[2] = PI4K3B, u[3] = CERT, u[4] = CERTTGNa, u[5] = PKDDAGa, u[6] = PKD, u[7] = PI4K3Ba
#θ_dynamicNames[1] = a11, θ_dynamicNames[2] = a12, θ_dynamicNames[3] = a21, θ_dynamicNames[4] = a22, θ_dynamicNames[5] = a31, θ_dynamicNames[6] = a32, θ_dynamicNames[7] = a33, θ_dynamicNames[8] = m11, θ_dynamicNames[9] = m22, θ_dynamicNames[10] = m31, θ_dynamicNames[11] = m33, θ_dynamicNames[12] = p11, θ_dynamicNames[13] = p12, θ_dynamicNames[14] = p13, θ_dynamicNames[15] = p21, θ_dynamicNames[16] = p22, θ_dynamicNames[17] = p31, θ_dynamicNames[18] = p32, θ_dynamicNames[19] = p33, θ_dynamicNames[20] = pu3, θ_dynamicNames[21] = pu4, θ_dynamicNames[22] = pu5, θ_dynamicNames[23] = pu6, θ_dynamicNames[24] = s12, θ_dynamicNames[25] = s21, θ_dynamicNames[26] = s31
##parameterInfo.nominalValue[20] = pu2_C 
#parameterInfo.nominalValue[33] = std_yPKDt_C 
#parameterInfo.nominalValue[34] = std_yPI4K3Bt_C 
#parameterInfo.nominalValue[35] = std_yCERTt_C 


function compute_h(u::AbstractVector, t::Real, θ_dynamic::AbstractVector, θ_observable::AbstractVector,
                   θ_nonDynamic::AbstractVector, parameterInfo::ParametersInfo, observableId::Symbol, 
                      parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :yCERTpRN24 
		observableParameter1_yCERTpRN24 = getObsOrSdParam(θ_observable, parameterMap)
		return u[3] * observableParameter1_yCERTpRN24 / ( u[3] + u[1] + u[4] ) 
	end

	if observableId == :yCERTt 
		return u[3] + u[1] + u[4] 
	end

	if observableId == :yPI4K3BpRN24 
		observableParameter1_yPI4K3BpRN24 = getObsOrSdParam(θ_observable, parameterMap)
		return u[7] * observableParameter1_yPI4K3BpRN24 / ( u[2] + u[7] ) 
	end

	if observableId == :yPI4K3Bt 
		return u[2] + u[7] 
	end

	if observableId == :yPKDpN0 
		observableParameter1_yPKDpN0 = getObsOrSdParam(θ_observable, parameterMap)
		return u[5] * observableParameter1_yPKDpN0 
	end

	if observableId == :yPKDpN24 
		observableParameter1_yPKDpN24 = getObsOrSdParam(θ_observable, parameterMap)
		return u[5] * observableParameter1_yPKDpN24 
	end

	if observableId == :yPKDpN25 
		observableParameter1_yPKDpN25 = getObsOrSdParam(θ_observable, parameterMap)
		return u[5] * observableParameter1_yPKDpN25 
	end

	if observableId == :yPKDt 
		return u[6] + u[5] 
	end

end

function compute_u0!(u0::AbstractVector, pODEProblem::AbstractVector) 

	#pODEProblem[1] = a21, pODEProblem[2] = m33, pODEProblem[3] = u4_bool1, pODEProblem[4] = u5_bool1, pODEProblem[5] = kb_NB142_70_dose, pODEProblem[6] = a22, pODEProblem[7] = kb_NB142_70_time, pODEProblem[8] = p12, pODEProblem[9] = p33, pODEProblem[10] = s31, pODEProblem[11] = PdBu_dose, pODEProblem[12] = a12, pODEProblem[13] = p22, pODEProblem[14] = a33, pODEProblem[15] = p13, pODEProblem[16] = cyt, pODEProblem[17] = pu5, pODEProblem[18] = u6_bool1, pODEProblem[19] = pu2, pODEProblem[20] = p31, pODEProblem[21] = pu4, pODEProblem[22] = PdBu_time, pODEProblem[23] = s12, pODEProblem[24] = m11, pODEProblem[25] = m31, pODEProblem[26] = u2, pODEProblem[27] = Ect_Expr_PI4K3beta_flag, pODEProblem[28] = a11, pODEProblem[29] = p32, pODEProblem[30] = p21, pODEProblem[31] = pu6, pODEProblem[32] = s21, pODEProblem[33] = pu3, pODEProblem[34] = m22, pODEProblem[35] = p11, pODEProblem[36] = Ect_Expr_CERT_flag, pODEProblem[37] = a31, pODEProblem[38] = u3_bool1, pODEProblem[39] = a32

	CERTERa = 3.19483885902e7 
	PI4K3B = 1.5775405394e6 
	CERT = 160797.7364 
	CERTTGNa = 4.20828286681e7 
	PKDDAGa = 123.8608 
	PKD = 466534.7994 
	PI4K3Ba = 332054.5041 

	u0 .= CERTERa, PI4K3B, CERT, CERTTGNa, PKDDAGa, PKD, PI4K3Ba
end

function compute_u0(pODEProblem::AbstractVector)::AbstractVector 

	#pODEProblem[1] = a21, pODEProblem[2] = m33, pODEProblem[3] = u4_bool1, pODEProblem[4] = u5_bool1, pODEProblem[5] = kb_NB142_70_dose, pODEProblem[6] = a22, pODEProblem[7] = kb_NB142_70_time, pODEProblem[8] = p12, pODEProblem[9] = p33, pODEProblem[10] = s31, pODEProblem[11] = PdBu_dose, pODEProblem[12] = a12, pODEProblem[13] = p22, pODEProblem[14] = a33, pODEProblem[15] = p13, pODEProblem[16] = cyt, pODEProblem[17] = pu5, pODEProblem[18] = u6_bool1, pODEProblem[19] = pu2, pODEProblem[20] = p31, pODEProblem[21] = pu4, pODEProblem[22] = PdBu_time, pODEProblem[23] = s12, pODEProblem[24] = m11, pODEProblem[25] = m31, pODEProblem[26] = u2, pODEProblem[27] = Ect_Expr_PI4K3beta_flag, pODEProblem[28] = a11, pODEProblem[29] = p32, pODEProblem[30] = p21, pODEProblem[31] = pu6, pODEProblem[32] = s21, pODEProblem[33] = pu3, pODEProblem[34] = m22, pODEProblem[35] = p11, pODEProblem[36] = Ect_Expr_CERT_flag, pODEProblem[37] = a31, pODEProblem[38] = u3_bool1, pODEProblem[39] = a32

	CERTERa = 3.19483885902e7 
	PI4K3B = 1.5775405394e6 
	CERT = 160797.7364 
	CERTTGNa = 4.20828286681e7 
	PKDDAGa = 123.8608 
	PKD = 466534.7994 
	PI4K3Ba = 332054.5041 

	 return [CERTERa, PI4K3B, CERT, CERTTGNa, PKDDAGa, PKD, PI4K3Ba]
end

function compute_σ(u::AbstractVector, t::Real, θ_sd::AbstractVector, θ_dynamic::AbstractVector, θ_nonDynamic::AbstractVector, 
                   parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :yCERTpRN24 
		noiseParameter1_yCERTpRN24 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_yCERTpRN24 
	end

	if observableId == :yCERTt 
		noiseParameter1_yCERTt = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_yCERTt 
	end

	if observableId == :yPI4K3BpRN24 
		noiseParameter1_yPI4K3BpRN24 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_yPI4K3BpRN24 
	end

	if observableId == :yPI4K3Bt 
		noiseParameter1_yPI4K3Bt = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_yPI4K3Bt 
	end

	if observableId == :yPKDpN0 
		noiseParameter1_yPKDpN0 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_yPKDpN0 
	end

	if observableId == :yPKDpN24 
		noiseParameter1_yPKDpN24 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_yPKDpN24 
	end

	if observableId == :yPKDpN25 
		noiseParameter1_yPKDpN25 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_yPKDpN25 
	end

	if observableId == :yPKDt 
		noiseParameter1_yPKDt = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_yPKDt 
	end

end