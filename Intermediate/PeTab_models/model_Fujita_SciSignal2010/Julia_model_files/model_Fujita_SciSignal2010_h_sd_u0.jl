#u[1] = pAkt_S6, u[2] = pAkt, u[3] = pS6, u[4] = EGFR, u[5] = pEGFR_Akt, u[6] = pEGFR, u[7] = Akt, u[8] = S6, u[9] = EGF_EGFR
#θ_dynamicNames[1] = EGFR_turnover, θ_dynamicNames[2] = init_AKT, θ_dynamicNames[3] = init_EGFR, θ_dynamicNames[4] = init_S6, θ_dynamicNames[5] = reaction_1_k1, θ_dynamicNames[6] = reaction_1_k2, θ_dynamicNames[7] = reaction_2_k1, θ_dynamicNames[8] = reaction_2_k2, θ_dynamicNames[9] = reaction_3_k1, θ_dynamicNames[10] = reaction_4_k1, θ_dynamicNames[11] = reaction_5_k1, θ_dynamicNames[12] = reaction_5_k2, θ_dynamicNames[13] = reaction_6_k1, θ_dynamicNames[14] = reaction_7_k1, θ_dynamicNames[15] = reaction_8_k1, θ_dynamicNames[16] = reaction_9_k1
#

function compute_h(u::AbstractVector, t::Real, θ_dynamic::AbstractVector, θ_observable::AbstractVector,
                   θ_nonDynamic::AbstractVector, parameterInfo::ParametersInfo, observableId::Symbol, 
                      parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :pAkt_tot 
		observableParameter1_pAkt_tot = getObsOrSdParam(θ_observable, parameterMap)
		return observableParameter1_pAkt_tot * ( u[2] + u[1] ) 
	end

	if observableId == :pEGFR_tot 
		observableParameter1_pEGFR_tot = getObsOrSdParam(θ_observable, parameterMap)
		return observableParameter1_pEGFR_tot * ( u[6] + u[5] ) 
	end

	if observableId == :pS6_tot 
		observableParameter1_pS6_tot = getObsOrSdParam(θ_observable, parameterMap)
		return u[3] * observableParameter1_pS6_tot 
	end

end

function compute_u0!(u0::AbstractVector, pODEProblem::AbstractVector) 

	#pODEProblem[1] = EGF_end, pODEProblem[2] = reaction_5_k1, pODEProblem[3] = reaction_2_k2, pODEProblem[4] = init_AKT, pODEProblem[5] = init_EGFR, pODEProblem[6] = EGF_bool1, pODEProblem[7] = EGF_rate, pODEProblem[8] = EGFR_turnover, pODEProblem[9] = reaction_1_k1, pODEProblem[10] = reaction_1_k2, pODEProblem[11] = reaction_8_k1, pODEProblem[12] = reaction_4_k1, pODEProblem[13] = reaction_6_k1, pODEProblem[14] = reaction_2_k1, pODEProblem[15] = init_S6, pODEProblem[16] = reaction_7_k1, pODEProblem[17] = reaction_9_k1, pODEProblem[18] = reaction_3_k1, pODEProblem[19] = reaction_5_k2, pODEProblem[20] = Cell, pODEProblem[21] = EGF_0

	pAkt_S6 = 0.0 
	pAkt = 0.0 
	pS6 = 0.0 
	EGFR = pODEProblem[5] 
	pEGFR_Akt = 0.0 
	pEGFR = 0.0 
	Akt = pODEProblem[4] 
	S6 = pODEProblem[15] 
	EGF_EGFR = 0.0 

	u0 .= pAkt_S6, pAkt, pS6, EGFR, pEGFR_Akt, pEGFR, Akt, S6, EGF_EGFR
end

function compute_u0(pODEProblem::AbstractVector)::AbstractVector 

	#pODEProblem[1] = EGF_end, pODEProblem[2] = reaction_5_k1, pODEProblem[3] = reaction_2_k2, pODEProblem[4] = init_AKT, pODEProblem[5] = init_EGFR, pODEProblem[6] = EGF_bool1, pODEProblem[7] = EGF_rate, pODEProblem[8] = EGFR_turnover, pODEProblem[9] = reaction_1_k1, pODEProblem[10] = reaction_1_k2, pODEProblem[11] = reaction_8_k1, pODEProblem[12] = reaction_4_k1, pODEProblem[13] = reaction_6_k1, pODEProblem[14] = reaction_2_k1, pODEProblem[15] = init_S6, pODEProblem[16] = reaction_7_k1, pODEProblem[17] = reaction_9_k1, pODEProblem[18] = reaction_3_k1, pODEProblem[19] = reaction_5_k2, pODEProblem[20] = Cell, pODEProblem[21] = EGF_0

	pAkt_S6 = 0.0 
	pAkt = 0.0 
	pS6 = 0.0 
	EGFR = pODEProblem[5] 
	pEGFR_Akt = 0.0 
	pEGFR = 0.0 
	Akt = pODEProblem[4] 
	S6 = pODEProblem[15] 
	EGF_EGFR = 0.0 

	 return [pAkt_S6, pAkt, pS6, EGFR, pEGFR_Akt, pEGFR, Akt, S6, EGF_EGFR]
end

function compute_σ(u::AbstractVector, t::Real, θ_sd::AbstractVector, θ_dynamic::AbstractVector, θ_nonDynamic::AbstractVector, 
                   parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :pAkt_tot 
		noiseParameter1_pAkt_tot = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_pAkt_tot 
	end

	if observableId == :pEGFR_tot 
		noiseParameter1_pEGFR_tot = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_pEGFR_tot 
	end

	if observableId == :pS6_tot 
		noiseParameter1_pS6_tot = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_pS6_tot 
	end

end