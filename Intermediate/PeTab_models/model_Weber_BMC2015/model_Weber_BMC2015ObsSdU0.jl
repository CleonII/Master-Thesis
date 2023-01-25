function evalYmod(u, t, dynPar, obsPar, nonDynParam, paramData, observableId, mapObsParam) 

	CERTERa, PI4K3B, CERT, CERTTGNa, PKDDAGa, PKD, PI4K3Ba= u 
	a11, a12, a21, a22, a31, a32, a33, m11, m22, m31, m33, p11, p12, p13, p21, p22, p31, p32, p33, pu3, pu4, pu5, pu6, s12, s21, s31 = dynPar 
	pu2_C = paramData.nominalValue[20] 
	std_yPKDt_C = paramData.nominalValue[33] 
	std_yPI4K3Bt_C = paramData.nominalValue[34] 
	std_yCERTt_C = paramData.nominalValue[35] 

	if observableId == :yCERTpRN24 
		observableParameter1_yCERTpRN24 = getObsOrSdParam(obsPar, mapObsParam)
		return CERT * observableParameter1_yCERTpRN24 / ( CERT + CERTERa + CERTTGNa ) 
	end

	if observableId == :yCERTt 
		return CERT + CERTERa + CERTTGNa 
	end

	if observableId == :yPI4K3BpRN24 
		observableParameter1_yPI4K3BpRN24 = getObsOrSdParam(obsPar, mapObsParam)
		return PI4K3Ba * observableParameter1_yPI4K3BpRN24 / ( PI4K3B + PI4K3Ba ) 
	end

	if observableId == :yPI4K3Bt 
		return PI4K3B + PI4K3Ba 
	end

	if observableId == :yPKDpN0 
		observableParameter1_yPKDpN0 = getObsOrSdParam(obsPar, mapObsParam)
		return PKDDAGa * observableParameter1_yPKDpN0 
	end

	if observableId == :yPKDpN24 
		observableParameter1_yPKDpN24 = getObsOrSdParam(obsPar, mapObsParam)
		return PKDDAGa * observableParameter1_yPKDpN24 
	end

	if observableId == :yPKDpN25 
		observableParameter1_yPKDpN25 = getObsOrSdParam(obsPar, mapObsParam)
		return PKDDAGa * observableParameter1_yPKDpN25 
	end

	if observableId == :yPKDt 
		return PKD + PKDDAGa 
	end

end

function evalU0!(u0Vec, paramVec) 

	a21, m33, u4_bool1, u5_bool1, kb_NB142_70_dose, a22, kb_NB142_70_time, p12, p33, s31, PdBu_dose, a12, p22, a33, p13, cyt, pu5, u6_bool1, pu2, p31, pu4, PdBu_time, s12, m11, m31, u2, Ect_Expr_PI4K3beta_flag, a11, p32, p21, pu6, s21, pu3, m22, p11, Ect_Expr_CERT_flag, a31, u3_bool1, a32 = paramVec 

	CERTERa = 3.19483885902e7 
	PI4K3B = 1.5775405394e6 
	CERT = 160797.7364 
	CERTTGNa = 4.20828286681e7 
	PKDDAGa = 123.8608 
	PKD = 466534.7994 
	PI4K3Ba = 332054.5041 

	u0Vec .= CERTERa, PI4K3B, CERT, CERTTGNa, PKDDAGa, PKD, PI4K3Ba
end

function evalU0(paramVec) 

	a21, m33, u4_bool1, u5_bool1, kb_NB142_70_dose, a22, kb_NB142_70_time, p12, p33, s31, PdBu_dose, a12, p22, a33, p13, cyt, pu5, u6_bool1, pu2, p31, pu4, PdBu_time, s12, m11, m31, u2, Ect_Expr_PI4K3beta_flag, a11, p32, p21, pu6, s21, pu3, m22, p11, Ect_Expr_CERT_flag, a31, u3_bool1, a32 = paramVec 

	CERTERa = 3.19483885902e7 
	PI4K3B = 1.5775405394e6 
	CERT = 160797.7364 
	CERTTGNa = 4.20828286681e7 
	PKDDAGa = 123.8608 
	PKD = 466534.7994 
	PI4K3Ba = 332054.5041 

	 return [CERTERa, PI4K3B, CERT, CERTTGNa, PKDDAGa, PKD, PI4K3Ba]
end

function evalSd!(u, t, sdPar, dynPar, nonDynParam, paramData, observableId, mapSdParam) 

	CERTERa, PI4K3B, CERT, CERTTGNa, PKDDAGa, PKD, PI4K3Ba= u 
	a11, a12, a21, a22, a31, a32, a33, m11, m22, m31, m33, p11, p12, p13, p21, p22, p31, p32, p33, pu3, pu4, pu5, pu6, s12, s21, s31 = dynPar 
	pu2_C = paramData.nominalValue[20] 
	std_yPKDt_C = paramData.nominalValue[33] 
	std_yPI4K3Bt_C = paramData.nominalValue[34] 
	std_yCERTt_C = paramData.nominalValue[35] 

	if observableId == :yCERTpRN24 
		noiseParameter1_yCERTpRN24 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_yCERTpRN24 
	end

	if observableId == :yCERTt 
		noiseParameter1_yCERTt = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_yCERTt 
	end

	if observableId == :yPI4K3BpRN24 
		noiseParameter1_yPI4K3BpRN24 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_yPI4K3BpRN24 
	end

	if observableId == :yPI4K3Bt 
		noiseParameter1_yPI4K3Bt = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_yPI4K3Bt 
	end

	if observableId == :yPKDpN0 
		noiseParameter1_yPKDpN0 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_yPKDpN0 
	end

	if observableId == :yPKDpN24 
		noiseParameter1_yPKDpN24 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_yPKDpN24 
	end

	if observableId == :yPKDpN25 
		noiseParameter1_yPKDpN25 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_yPKDpN25 
	end

	if observableId == :yPKDt 
		noiseParameter1_yPKDt = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_yPKDt 
	end

end