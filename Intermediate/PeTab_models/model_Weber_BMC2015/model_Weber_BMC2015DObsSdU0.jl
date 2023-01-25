function evalDYmodDu(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	CERTERa, PI4K3B, CERT, CERTTGNa, PKDDAGa, PKD, PI4K3Ba= u 
	a21, m33, u4_bool1, u5_bool1, kb_NB142_70_dose, a22, kb_NB142_70_time, p12, p33, s31, PdBu_dose, a12, p22, a33, p13, cyt, pu5, u6_bool1, pu2, p31, pu4, PdBu_time, s12, m11, m31, u2, Ect_Expr_PI4K3beta_flag, a11, p32, p21, pu6, s21, pu3, m22, p11, Ect_Expr_CERT_flag, a31, u3_bool1, a32 = p 
	if observableId == :yCERTpRN24 
		observableParameter1_yCERTpRN24 = getObsOrSdParam(obsPar, mapObsParam)
		out[1] = (-CERT*observableParameter1_yCERTpRN24) / ((CERT + CERTERa + CERTTGNa)^2)
		out[3] = (CERTERa*observableParameter1_yCERTpRN24 + CERTTGNa*observableParameter1_yCERTpRN24) / ((CERT + CERTERa + CERTTGNa)^2)
		out[4] = (-CERT*observableParameter1_yCERTpRN24) / ((CERT + CERTERa + CERTTGNa)^2)
		return nothing
	end

	if observableId == :yCERTt 
		out[1] = 1
		out[3] = 1
		out[4] = 1
		return nothing
	end

	if observableId == :yPI4K3BpRN24 
		observableParameter1_yPI4K3BpRN24 = getObsOrSdParam(obsPar, mapObsParam)
		out[2] = (-PI4K3Ba*observableParameter1_yPI4K3BpRN24) / ((PI4K3B + PI4K3Ba)^2)
		out[7] = (PI4K3B*observableParameter1_yPI4K3BpRN24) / ((PI4K3B + PI4K3Ba)^2)
		return nothing
	end

	if observableId == :yPI4K3Bt 
		out[2] = 1
		out[7] = 1
		return nothing
	end

	if observableId == :yPKDpN0 
		observableParameter1_yPKDpN0 = getObsOrSdParam(obsPar, mapObsParam)
		out[5] = observableParameter1_yPKDpN0
		return nothing
	end

	if observableId == :yPKDpN24 
		observableParameter1_yPKDpN24 = getObsOrSdParam(obsPar, mapObsParam)
		out[5] = observableParameter1_yPKDpN24
		return nothing
	end

	if observableId == :yPKDpN25 
		observableParameter1_yPKDpN25 = getObsOrSdParam(obsPar, mapObsParam)
		out[5] = observableParameter1_yPKDpN25
		return nothing
	end

	if observableId == :yPKDt 
		out[5] = 1
		out[6] = 1
		return nothing
	end

end

function evalDYmodDp(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	CERTERa, PI4K3B, CERT, CERTTGNa, PKDDAGa, PKD, PI4K3Ba= u 
	a21, m33, u4_bool1, u5_bool1, kb_NB142_70_dose, a22, kb_NB142_70_time, p12, p33, s31, PdBu_dose, a12, p22, a33, p13, cyt, pu5, u6_bool1, pu2, p31, pu4, PdBu_time, s12, m11, m31, u2, Ect_Expr_PI4K3beta_flag, a11, p32, p21, pu6, s21, pu3, m22, p11, Ect_Expr_CERT_flag, a31, u3_bool1, a32 = p 
	if observableId == :yCERTpRN24 
		return nothing
	end

	if observableId == :yCERTt 
		return nothing
	end

	if observableId == :yPI4K3BpRN24 
		return nothing
	end

	if observableId == :yPI4K3Bt 
		return nothing
	end

	if observableId == :yPKDpN0 
		return nothing
	end

	if observableId == :yPKDpN24 
		return nothing
	end

	if observableId == :yPKDpN25 
		return nothing
	end

	if observableId == :yPKDt 
		return nothing
	end

end

function evalDSdDu!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	CERTERa, PI4K3B, CERT, CERTTGNa, PKDDAGa, PKD, PI4K3Ba= u 
	a21, m33, u4_bool1, u5_bool1, kb_NB142_70_dose, a22, kb_NB142_70_time, p12, p33, s31, PdBu_dose, a12, p22, a33, p13, cyt, pu5, u6_bool1, pu2, p31, pu4, PdBu_time, s12, m11, m31, u2, Ect_Expr_PI4K3beta_flag, a11, p32, p21, pu6, s21, pu3, m22, p11, Ect_Expr_CERT_flag, a31, u3_bool1, a32 = p 
	if observableId == "yCERTpRN24" 
		return nothing
	end

	if observableId == "yCERTt" 
		return nothing
	end

	if observableId == "yPI4K3BpRN24" 
		return nothing
	end

	if observableId == "yPI4K3Bt" 
		return nothing
	end

	if observableId == "yPKDpN0" 
		return nothing
	end

	if observableId == "yPKDpN24" 
		return nothing
	end

	if observableId == "yPKDpN25" 
		return nothing
	end

	if observableId == "yPKDt" 
		return nothing
	end

end

function evalDSdDp!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	CERTERa, PI4K3B, CERT, CERTTGNa, PKDDAGa, PKD, PI4K3Ba= u 
	a21, m33, u4_bool1, u5_bool1, kb_NB142_70_dose, a22, kb_NB142_70_time, p12, p33, s31, PdBu_dose, a12, p22, a33, p13, cyt, pu5, u6_bool1, pu2, p31, pu4, PdBu_time, s12, m11, m31, u2, Ect_Expr_PI4K3beta_flag, a11, p32, p21, pu6, s21, pu3, m22, p11, Ect_Expr_CERT_flag, a31, u3_bool1, a32 = p 
	if observableId == "yCERTpRN24" 
		return nothing
	end

	if observableId == "yCERTt" 
		return nothing
	end

	if observableId == "yPI4K3BpRN24" 
		return nothing
	end

	if observableId == "yPI4K3Bt" 
		return nothing
	end

	if observableId == "yPKDpN0" 
		return nothing
	end

	if observableId == "yPKDpN24" 
		return nothing
	end

	if observableId == "yPKDpN25" 
		return nothing
	end

	if observableId == "yPKDt" 
		return nothing
	end

end

