#u[1] = STAT5A, u[2] = pApA, u[3] = nucpApB, u[4] = nucpBpB, u[5] = STAT5B, u[6] = pApB, u[7] = nucpApA, u[8] = pBpB
#θ_dynamicNames[1] = Epo_degradation_BaF3, θ_dynamicNames[2] = k_exp_hetero, θ_dynamicNames[3] = k_exp_homo, θ_dynamicNames[4] = k_imp_hetero, θ_dynamicNames[5] = k_imp_homo, θ_dynamicNames[6] = k_phos
##parameterInfo.nominalValue[7] = ratio_C 
#parameterInfo.nominalValue[11] = specC17_C 


function compute_h(u::AbstractVector, t::Real, θ_dynamic::AbstractVector, θ_observable::AbstractVector,
                   θ_nonDynamic::AbstractVector, parameterInfo::ParametersInfo, observableId::Symbol, 
                      parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :pSTAT5A_rel 
		return ( 100 * u[6] + 200 * u[2] * parameterInfo.nominalValue[11] ) / ( u[6] + u[1] * parameterInfo.nominalValue[11] + 2 * u[2] * parameterInfo.nominalValue[11] ) 
	end

	if observableId == :pSTAT5B_rel 
		return - ( 100 * u[6] - 200 * u[8] * ( parameterInfo.nominalValue[11] - 1 ) ) / ( ( u[5] * ( parameterInfo.nominalValue[11] - 1 ) - u[6] ) + 2 * u[8] * ( parameterInfo.nominalValue[11] - 1 ) ) 
	end

	if observableId == :rSTAT5A_rel 
		return ( 100 * u[6] + 100 * u[1] * parameterInfo.nominalValue[11] + 200 * u[2] * parameterInfo.nominalValue[11] ) / ( 2 * u[6] + u[1] * parameterInfo.nominalValue[11] + 2 * u[2] * parameterInfo.nominalValue[11] - u[5] * ( parameterInfo.nominalValue[11] - 1 ) - 2 * u[8] * ( parameterInfo.nominalValue[11] - 1 ) ) 
	end

end

function compute_u0!(u0::AbstractVector, pODEProblem::AbstractVector) 

	#pODEProblem[1] = ratio, pODEProblem[2] = k_imp_homo, pODEProblem[3] = k_exp_hetero, pODEProblem[4] = cyt, pODEProblem[5] = k_phos, pODEProblem[6] = specC17, pODEProblem[7] = Epo_degradation_BaF3, pODEProblem[8] = k_exp_homo, pODEProblem[9] = nuc, pODEProblem[10] = k_imp_hetero

	STAT5A = 207.6 * pODEProblem[1] 
	pApA = 0.0 
	nucpApB = 0.0 
	nucpBpB = 0.0 
	STAT5B = 207.6 - 207.6 * pODEProblem[1] 
	pApB = 0.0 
	nucpApA = 0.0 
	pBpB = 0.0 

	u0 .= STAT5A, pApA, nucpApB, nucpBpB, STAT5B, pApB, nucpApA, pBpB
end

function compute_u0(pODEProblem::AbstractVector)::AbstractVector 

	#pODEProblem[1] = ratio, pODEProblem[2] = k_imp_homo, pODEProblem[3] = k_exp_hetero, pODEProblem[4] = cyt, pODEProblem[5] = k_phos, pODEProblem[6] = specC17, pODEProblem[7] = Epo_degradation_BaF3, pODEProblem[8] = k_exp_homo, pODEProblem[9] = nuc, pODEProblem[10] = k_imp_hetero

	STAT5A = 207.6 * pODEProblem[1] 
	pApA = 0.0 
	nucpApB = 0.0 
	nucpBpB = 0.0 
	STAT5B = 207.6 - 207.6 * pODEProblem[1] 
	pApB = 0.0 
	nucpApA = 0.0 
	pBpB = 0.0 

	 return [STAT5A, pApA, nucpApB, nucpBpB, STAT5B, pApB, nucpApA, pBpB]
end

function compute_σ(u::AbstractVector, t::Real, θ_sd::AbstractVector, θ_dynamic::AbstractVector, θ_nonDynamic::AbstractVector, 
                   parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :pSTAT5A_rel 
		noiseParameter1_pSTAT5A_rel = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_pSTAT5A_rel 
	end

	if observableId == :pSTAT5B_rel 
		noiseParameter1_pSTAT5B_rel = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_pSTAT5B_rel 
	end

	if observableId == :rSTAT5A_rel 
		noiseParameter1_rSTAT5A_rel = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_rSTAT5A_rel 
	end

end