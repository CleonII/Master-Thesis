function evalYmod(u, t, dynPar, obsPar, paramData, obsData, observableId, simulationId) 

	CERTERa, PI4K3B, CERT, CERTTGNa, PKDDAGa, PKD, PI4K3Ba, u4, u3, u5, dummyVariable= u 
	a11, a12, a21, a22, a31, a32, a33, m11, m22, m31, m33, p11, p12, p13, p21, p22, p31, p32, p33, pu3, pu4, pu5, pu6, s12, s21, s31 = dynPar 
	pu2_C = paramData.paramVal[20] 
	std_yPKDt_C = paramData.paramVal[33] 
	std_yPI4K3Bt_C = paramData.paramVal[34] 
	std_yCERTt_C = paramData.paramVal[35] 

	if observableId == "yCERTpRN24" 
		observableParameter1_yCERTpRN24 = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return CERT * observableParameter1_yCERTpRN24 / ( CERT + CERTERa + CERTTGNa ) 
	end

	if observableId == "yCERTt" 
		return CERT + CERTERa + CERTTGNa 
	end

	if observableId == "yPI4K3BpRN24" 
		observableParameter1_yPI4K3BpRN24 = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return PI4K3Ba * observableParameter1_yPI4K3BpRN24 / ( PI4K3B + PI4K3Ba ) 
	end

	if observableId == "yPI4K3Bt" 
		return PI4K3B + PI4K3Ba 
	end

	if observableId == "yPKDpN0" 
		observableParameter1_yPKDpN0 = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return PKDDAGa * observableParameter1_yPKDpN0 
	end

	if observableId == "yPKDpN24" 
		observableParameter1_yPKDpN24 = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return PKDDAGa * observableParameter1_yPKDpN24 
	end

	if observableId == "yPKDpN25" 
		observableParameter1_yPKDpN25 = getObsOrSdParam(obsPar, paramData, obsData, observableId, simulationId, t)
		return PKDDAGa * observableParameter1_yPKDpN25 
	end

	if observableId == "yPKDt" 
		return PKD + PKDDAGa 
	end

end

function evalU0!(u0Vec, paramVec) 

	pu4, s31, m31, p31, a31, p32, p22, pu3, s21, p21, a21, m22, p33, a32, m33, a33, kb_NB142_70_time, kb_NB142_70_dose, p11, p12, pu5, m11, p13, a12, pu6, u2, a11, pu2, s12, a22, PdBu_time, Ect_Expr_CERT_flag, Ect_Expr_PI4K3beta_flag, PdBu_dose = paramVec 

	CERTERa = 3.19483885902e7 
	PI4K3B = 1.5775405394e6 
	CERT = 160797.7364 
	CERTTGNa = 4.20828286681e7 
	PKDDAGa = 123.8608 
	PKD = 466534.7994 
	PI4K3Ba = 332054.5041 
	u4 = 0.0 
	u3 = 0.0 
	u5 = 0.0 
	dummyVariable = 0.0 

	u0Vec .= CERTERa, PI4K3B, CERT, CERTTGNa, PKDDAGa, PKD, PI4K3Ba, u4, u3, u5, dummyVariable
end

function evalSd!(u, t, sdPar, dynPar, paramData, obsData, observableId, simulationId) 

	CERTERa, PI4K3B, CERT, CERTTGNa, PKDDAGa, PKD, PI4K3Ba, u4, u3, u5, dummyVariable= u 
	a11, a12, a21, a22, a31, a32, a33, m11, m22, m31, m33, p11, p12, p13, p21, p22, p31, p32, p33, pu3, pu4, pu5, pu6, s12, s21, s31 = dynPar 
	pu2_C = paramData.paramVal[20] 
	std_yPKDt_C = paramData.paramVal[33] 
	std_yPI4K3Bt_C = paramData.paramVal[34] 
	std_yCERTt_C = paramData.paramVal[35] 

	if observableId == "yCERTpRN24" 
		noiseParameter1_yCERTpRN24 = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_yCERTpRN24 
	end

	if observableId == "yCERTt" 
		noiseParameter1_yCERTt = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_yCERTt 
	end

	if observableId == "yPI4K3BpRN24" 
		noiseParameter1_yPI4K3BpRN24 = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_yPI4K3BpRN24 
	end

	if observableId == "yPI4K3Bt" 
		noiseParameter1_yPI4K3Bt = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_yPI4K3Bt 
	end

	if observableId == "yPKDpN0" 
		noiseParameter1_yPKDpN0 = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_yPKDpN0 
	end

	if observableId == "yPKDpN24" 
		noiseParameter1_yPKDpN24 = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_yPKDpN24 
	end

	if observableId == "yPKDpN25" 
		noiseParameter1_yPKDpN25 = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_yPKDpN25 
	end

	if observableId == "yPKDt" 
		noiseParameter1_yPKDt = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_yPKDt 
	end

end