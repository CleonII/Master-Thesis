#u[1] = CERTERa, u[2] = PI4K3B, u[3] = CERT, u[4] = CERTTGNa, u[5] = PKDDAGa, u[6] = PKD, u[7] = PI4K3Ba
#pODEProblem[1] = a21, pODEProblem[2] = m33, pODEProblem[3] = u4_bool1, pODEProblem[4] = u5_bool1, pODEProblem[5] = kb_NB142_70_dose, pODEProblem[6] = a22, pODEProblem[7] = kb_NB142_70_time, pODEProblem[8] = p12, pODEProblem[9] = p33, pODEProblem[10] = s31, pODEProblem[11] = PdBu_dose, pODEProblem[12] = a12, pODEProblem[13] = p22, pODEProblem[14] = a33, pODEProblem[15] = p13, pODEProblem[16] = cyt, pODEProblem[17] = pu5, pODEProblem[18] = u6_bool1, pODEProblem[19] = pu2, pODEProblem[20] = p31, pODEProblem[21] = pu4, pODEProblem[22] = PdBu_time, pODEProblem[23] = s12, pODEProblem[24] = m11, pODEProblem[25] = m31, pODEProblem[26] = u2, pODEProblem[27] = Ect_Expr_PI4K3beta_flag, pODEProblem[28] = a11, pODEProblem[29] = p32, pODEProblem[30] = p21, pODEProblem[31] = pu6, pODEProblem[32] = s21, pODEProblem[33] = pu3, pODEProblem[34] = m22, pODEProblem[35] = p11, pODEProblem[36] = Ect_Expr_CERT_flag, pODEProblem[37] = a31, pODEProblem[38] = u3_bool1, pODEProblem[39] = a32
#
function compute_∂h∂u!(u, t::Real, pODEProblem::AbstractVector, θ_observable::AbstractVector,                    
                       θ_nonDynamic::AbstractVector, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :yCERTpRN24 
		observableParameter1_yCERTpRN24 = getObsOrSdParam(θ_observable, parameterMap)
		out[1] = (-u[3]*observableParameter1_yCERTpRN24) / ((u[3] + u[1] + u[4])^2)
		out[3] = (u[1]*observableParameter1_yCERTpRN24 + u[4]*observableParameter1_yCERTpRN24) / ((u[3] + u[1] + u[4])^2)
		out[4] = (-u[3]*observableParameter1_yCERTpRN24) / ((u[3] + u[1] + u[4])^2)
		return nothing
	end

	if observableId == :yCERTt 
		out[1] = 1
		out[3] = 1
		out[4] = 1
		return nothing
	end

	if observableId == :yPI4K3BpRN24 
		observableParameter1_yPI4K3BpRN24 = getObsOrSdParam(θ_observable, parameterMap)
		out[2] = (-u[7]*observableParameter1_yPI4K3BpRN24) / ((u[2] + u[7])^2)
		out[7] = (u[2]*observableParameter1_yPI4K3BpRN24) / ((u[2] + u[7])^2)
		return nothing
	end

	if observableId == :yPI4K3Bt 
		out[2] = 1
		out[7] = 1
		return nothing
	end

	if observableId == :yPKDpN0 
		observableParameter1_yPKDpN0 = getObsOrSdParam(θ_observable, parameterMap)
		out[5] = observableParameter1_yPKDpN0
		return nothing
	end

	if observableId == :yPKDpN24 
		observableParameter1_yPKDpN24 = getObsOrSdParam(θ_observable, parameterMap)
		out[5] = observableParameter1_yPKDpN24
		return nothing
	end

	if observableId == :yPKDpN25 
		observableParameter1_yPKDpN25 = getObsOrSdParam(θ_observable, parameterMap)
		out[5] = observableParameter1_yPKDpN25
		return nothing
	end

	if observableId == :yPKDt 
		out[5] = 1
		out[6] = 1
		return nothing
	end

end

function compute_∂h∂p!(u, t::Real, pODEProblem::AbstractVector, θ_observable::AbstractVector,
                       θ_nonDynamic::AbstractVector, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
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

function compute_∂σ∂σu!(u, t::Real, θ_sd::AbstractVector, pODEProblem::AbstractVector, θ_nonDynamic::AbstractVector, 
                        parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
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

function compute_∂σ∂σp!(u, t::Real, θ_sd::AbstractVector, pODEProblem::AbstractVector, θ_nonDynamic::AbstractVector, 
                        parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
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

