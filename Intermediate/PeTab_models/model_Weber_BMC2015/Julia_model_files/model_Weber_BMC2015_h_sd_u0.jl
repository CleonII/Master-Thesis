#u[1] = CERTERa, u[2] = PI4K3B, u[3] = CERT, u[4] = CERTTGNa, u[5] = PKDDAGa, u[6] = PKD, u[7] = PI4K3Ba
#pODEProblemNames[1] = a21, pODEProblemNames[2] = m33, pODEProblemNames[3] = u4_bool1, pODEProblemNames[4] = u5_bool1, pODEProblemNames[5] = kb_NB142_70_dose, pODEProblemNames[6] = a22, pODEProblemNames[7] = kb_NB142_70_time, pODEProblemNames[8] = p12, pODEProblemNames[9] = p33, pODEProblemNames[10] = s31, pODEProblemNames[11] = PdBu_dose, pODEProblemNames[12] = a12, pODEProblemNames[13] = p22, pODEProblemNames[14] = a33, pODEProblemNames[15] = p13, pODEProblemNames[16] = cyt, pODEProblemNames[17] = pu5, pODEProblemNames[18] = u6_bool1, pODEProblemNames[19] = pu2, pODEProblemNames[20] = p31, pODEProblemNames[21] = pu4, pODEProblemNames[22] = PdBu_time, pODEProblemNames[23] = s12, pODEProblemNames[24] = m11, pODEProblemNames[25] = m31, pODEProblemNames[26] = u2, pODEProblemNames[27] = Ect_Expr_PI4K3beta_flag, pODEProblemNames[28] = a11, pODEProblemNames[29] = p32, pODEProblemNames[30] = p21, pODEProblemNames[31] = pu6, pODEProblemNames[32] = s21, pODEProblemNames[33] = pu3, pODEProblemNames[34] = m22, pODEProblemNames[35] = p11, pODEProblemNames[36] = Ect_Expr_CERT_flag, pODEProblemNames[37] = a31, pODEProblemNames[38] = u3_bool1, pODEProblemNames[39] = a32
##parameterInfo.nominalValue[20] = pu2_C 
#parameterInfo.nominalValue[33] = std_yPKDt_C 
#parameterInfo.nominalValue[34] = std_yPI4K3Bt_C 
#parameterInfo.nominalValue[35] = std_yCERTt_C 


function compute_h(u::AbstractVector, t::Real, pODEProblem::AbstractVector, θ_observable::AbstractVector,
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

function compute_σ(u::AbstractVector, t::Real, θ_sd::AbstractVector, pODEProblem::AbstractVector, θ_nonDynamic::AbstractVector, 
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