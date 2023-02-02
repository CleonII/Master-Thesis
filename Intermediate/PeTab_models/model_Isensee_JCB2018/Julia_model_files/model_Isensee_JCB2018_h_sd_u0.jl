#u[1] = pAC, u[2] = Rp8_Br_cAMPS, u[3] = Rp8_pCPT_cAMPS, u[4] = PDE, u[5] = Rp_cAMPS, u[6] = RII_2, u[7] = RIIp_Rp8_Br_cAMPS_C_2, u[8] = cAMP, u[9] = RIIp_Sp8_Br_cAMPS_C_2, u[10] = IBMX, u[11] = AC_Fsk, u[12] = RIIp_C_2, u[13] = Sp8_Br_cAMPS, u[14] = RII_C_2, u[15] = RIIp_Sp8_Br_cAMPS_2, u[16] = Csub, u[17] = Csub_H89, u[18] = AC, u[19] = RIIp_cAMP_2, u[20] = pAC_Fsk, u[21] = RIIp_Rp8_pCPT_cAMPS_C_2, u[22] = pPDE, u[23] = RIIp_Rp_cAMPS_C_2, u[24] = RIIp_2, u[25] = RIIp_cAMP_C_2
#θ_dynamicNames[1] = KD_Fsk, θ_dynamicNames[2] = KD_H89, θ_dynamicNames[3] = KD_IBMX, θ_dynamicNames[4] = KD_cAMP, θ_dynamicNames[5] = kdeg_cAMP, θ_dynamicNames[6] = kdeg_cAMP_free, θ_dynamicNames[7] = kf_Fsk, θ_dynamicNames[8] = kf_H89, θ_dynamicNames[9] = kf_RII_2__RII_C_2, θ_dynamicNames[10] = kf_RII_C_2__RII_2, θ_dynamicNames[11] = kf_RII_C_2__RIIp_C_2, θ_dynamicNames[12] = kf_RIIp_2__RII_2, θ_dynamicNames[13] = kf_RIIp_C_2__RII_C_2, θ_dynamicNames[14] = kf_RIIp_cAMP_C_2__RIIp_2, θ_dynamicNames[15] = kf_cAMP, θ_dynamicNames[16] = ki_IBMX, θ_dynamicNames[17] = ki_Rp8_Br_cAMPS_pAB, θ_dynamicNames[18] = ki_Rp8_pCPT_cAMPS_pAB, θ_dynamicNames[19] = ki_Rp_cAMPS_pAB, θ_dynamicNames[20] = ki_Sp8_Br_cAMPS_AM, θ_dynamicNames[21] = ks_AC_cAMP, θ_dynamicNames[22] = xi_AC_cAMP_Fsk, θ_dynamicNames[23] = xi_KD_Rp8_Br_cAMPS, θ_dynamicNames[24] = xi_KD_Rp8_pCPT_cAMPS, θ_dynamicNames[25] = xi_KD_Rp_cAMPS, θ_dynamicNames[26] = xi_KD_Sp8_Br_cAMPS, θ_dynamicNames[27] = xi_b_Rp8_Br_cAMPS, θ_dynamicNames[28] = xi_b_Rp8_pCPT_cAMPS, θ_dynamicNames[29] = xi_b_Rp_cAMPS, θ_dynamicNames[30] = xi_b_Sp8_Br_cAMPS, θ_dynamicNames[31] = xi_kf_RII_2__RII_C_2, θ_dynamicNames[32] = xi_kf_RII_C_2__RII_2
#θ_nonDynamic[1] = b_Calpha_global, θ_nonDynamic[2] = b_pRII_global, θ_nonDynamic[3] = rel_open, θ_nonDynamic[4] = s_Calpha_global, θ_nonDynamic[5] = s_pRII_Western, θ_nonDynamic[6] = s_pRII_global, θ_nonDynamic[7] = rho_Calpha_Microscopy, θ_nonDynamic[8] = rho_pRII_Microscopy, θ_nonDynamic[9] = rho_pRII_Western, θ_nonDynamic[10] = xi_rel_open
#parameterInfo.nominalValue[1] = AC_total_C 
#parameterInfo.nominalValue[5] = KD_PDE_Csub_C 
#parameterInfo.nominalValue[7] = PDE_total_C 
#parameterInfo.nominalValue[8] = RII2_total_C 
#parameterInfo.nominalValue[13] = kdp_AC_C 
#parameterInfo.nominalValue[16] = kf_PDE_Csub_C 
#parameterInfo.nominalValue[29] = kp_AC_C 
#parameterInfo.nominalValue[32] = s_Calpha_JI09_160201_Drg453_452_CycNuc_C 
#parameterInfo.nominalValue[52] = xi_i_Rp8_Br_cAMPS_pAB_C 
#parameterInfo.nominalValue[53] = xi_i_Rp8_pCPT_cAMPS_pAB_C 
#parameterInfo.nominalValue[54] = xi_i_Rp_cAMPS_pAB_C 
#parameterInfo.nominalValue[55] = xi_i_Sp8_Br_cAMPS_AM_C 
#parameterInfo.nominalValue[58] = xi_pAC_C 
#parameterInfo.nominalValue[59] = xi_pPDE_C 


function compute_h(u::AbstractVector, t::Real, θ_dynamic::AbstractVector, θ_observable::AbstractVector,
                   θ_nonDynamic::AbstractVector, parameterInfo::ParametersInfo, observableId::Symbol, 
                      parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :pRII_Microscopy 
		observableParameter1_pRII_Microscopy = getObsOrSdParam(θ_observable, parameterMap)
		return observableParameter1_pRII_Microscopy * ( θ_nonDynamic[2] + θ_nonDynamic[6] * ( 2 * u[24] + 2 * u[19] + 2 * u[15] + 2 * θ_nonDynamic[3] * ( u[12] + u[23] + u[7] + u[21] ) + 2 * ( u[25] + u[9] ) * ( θ_nonDynamic[3] - θ_nonDynamic[10] * ( θ_nonDynamic[3] - 1 ) ) ) ) 
	end

	if observableId == :pRII_Western 
		return θ_nonDynamic[5] * ( 2 * u[24] + 2 * u[12] + 2 * u[19] + 2 * u[25] + 2 * u[23] + 2 * u[15] + 2 * u[7] + 2 * u[9] + 2 * u[21] ) 
	end

	if observableId == :Calpha_Microscopy 
		observableParameter1_Calpha_Microscopy = getObsOrSdParam(θ_observable, parameterMap)
		return observableParameter1_Calpha_Microscopy * ( θ_nonDynamic[1] + θ_nonDynamic[4] * ( 2 * u[16] + 2 * u[17] + 2 * θ_nonDynamic[3] * ( u[12] + u[23] + u[7] + u[21] ) + 2 * ( u[25] + u[9] ) * ( θ_nonDynamic[3] - θ_nonDynamic[10] * ( θ_nonDynamic[3] - 1 ) ) ) ) 
	end

end

function compute_u0!(u0::AbstractVector, pODEProblem::AbstractVector) 

	#pODEProblem[1] = ki_Rp8_pCPT_cAMPS_pAB, pODEProblem[2] = xi_b_Rp_cAMPS, pODEProblem[3] = RII2_total, pODEProblem[4] = H89_level, pODEProblem[5] = fourABnOH_incubation_time, pODEProblem[6] = KD_Fsk, pODEProblem[7] = kdeg_cAMP_free, pODEProblem[8] = Rp8_Br_cAMPS_pAB_level, pODEProblem[9] = xi_KD_Rp8_Br_cAMPS, pODEProblem[10] = kf_PDE_Csub, pODEProblem[11] = Sp8_Br_cAMPS_AM_level, pODEProblem[12] = xi_kf_RII_C_2__RII_2, pODEProblem[13] = kf_RIIp_2__RII_2, pODEProblem[14] = kf_cAMP, pODEProblem[15] = IBMX_time, pODEProblem[16] = H89_bool1, pODEProblem[17] = Rp_cAMPS_pAB_incubation_time, pODEProblem[18] = kf_H89, pODEProblem[19] = kf_RII_C_2__RII_2, pODEProblem[20] = kf_RII_C_2__RIIp_C_2, pODEProblem[21] = xi_i_Rp8_pCPT_cAMPS_pAB, pODEProblem[22] = xi_pAC, pODEProblem[23] = fourABnOH_level, pODEProblem[24] = ki_Sp8_Br_cAMPS_AM, pODEProblem[25] = Rp8_pCPT_cAMPS_pAB_bool1, pODEProblem[26] = xi_b_Sp8_Br_cAMPS, pODEProblem[27] = xi_b_Rp8_Br_cAMPS, pODEProblem[28] = H89_time, pODEProblem[29] = Rp_cAMPS_pAB_bool1, pODEProblem[30] = kdeg_cAMP, pODEProblem[31] = xi_AC_cAMP_Fsk, pODEProblem[32] = xi_b_Rp8_pCPT_cAMPS, pODEProblem[33] = ki_Rp_cAMPS_pAB, pODEProblem[34] = xi_i_Rp_cAMPS_pAB, pODEProblem[35] = KD_PDE_Csub, pODEProblem[36] = ki_IBMX, pODEProblem[37] = IBMXex_bool1, pODEProblem[38] = Fsk_time, pODEProblem[39] = Rp8_pCPT_cAMPS_pAB_incubation_time, pODEProblem[40] = PDE_total, pODEProblem[41] = Rp8_Br_cAMPS_pAB_bool1, pODEProblem[42] = ki_Rp8_Br_cAMPS_pAB, pODEProblem[43] = xi_kf_RII_2__RII_C_2, pODEProblem[44] = default, pODEProblem[45] = xi_i_Rp8_Br_cAMPS_pAB, pODEProblem[46] = Sp8_Br_cAMPS_AM_time, pODEProblem[47] = xi_KD_Rp8_pCPT_cAMPS, pODEProblem[48] = kp_AC, pODEProblem[49] = xi_pPDE, pODEProblem[50] = xi_KD_Rp_cAMPS, pODEProblem[51] = Rp8_Br_cAMPS_pAB_incubation_time, pODEProblem[52] = AC_total, pODEProblem[53] = kf_RIIp_C_2__RII_C_2, pODEProblem[54] = xi_KD_Sp8_Br_cAMPS, pODEProblem[55] = nuc, pODEProblem[56] = Sp8_Br_cAMPS_AM_bool1, pODEProblem[57] = kf_RIIp_cAMP_C_2__RIIp_2, pODEProblem[58] = KD_cAMP, pODEProblem[59] = KD_IBMX, pODEProblem[60] = kf_Fsk, pODEProblem[61] = xi_i_Sp8_Br_cAMPS_AM, pODEProblem[62] = cyt, pODEProblem[63] = KD_H89, pODEProblem[64] = ks_AC_cAMP, pODEProblem[65] = fourABnOH_bool1, pODEProblem[66] = Rp8_pCPT_cAMPS_pAB_level, pODEProblem[67] = Rp_cAMPS_pAB_level, pODEProblem[68] = kf_RII_2__RII_C_2, pODEProblem[69] = Fsk_level, pODEProblem[70] = Fsk_bool1, pODEProblem[71] = kdp_AC, pODEProblem[72] = IBMX_level

	pAC = 0.0 
	Rp8_Br_cAMPS = 0.0 
	Rp8_pCPT_cAMPS = 0.0 
	PDE = 1.0 
	Rp_cAMPS = 0.0 
	RII_2 = 0.057671482854616 
	RIIp_Rp8_Br_cAMPS_C_2 = 0.0 
	cAMP = 0.0575218758977949 
	RIIp_Sp8_Br_cAMPS_C_2 = 0.0 
	IBMX = 0.0 
	AC_Fsk = 0.0 
	RIIp_C_2 = 0.287974352643203 
	Sp8_Br_cAMPS = 0.0 
	RII_C_2 = 0.514562213223039 
	RIIp_Sp8_Br_cAMPS_2 = 0.0 
	Csub = 0.192051405433538 
	Csub_H89 = 0.0 
	AC = 1.0 
	RIIp_cAMP_2 = 0.000400416634006059 
	pAC_Fsk = 0.0 
	RIIp_Rp8_pCPT_cAMPS_C_2 = 0.0 
	pPDE = 0.0 
	RIIp_Rp_cAMPS_C_2 = 0.0 
	RIIp_2 = 0.133979505944916 
	RIIp_cAMP_C_2 = 0.00541202870022029 

	u0 .= pAC, Rp8_Br_cAMPS, Rp8_pCPT_cAMPS, PDE, Rp_cAMPS, RII_2, RIIp_Rp8_Br_cAMPS_C_2, cAMP, RIIp_Sp8_Br_cAMPS_C_2, IBMX, AC_Fsk, RIIp_C_2, Sp8_Br_cAMPS, RII_C_2, RIIp_Sp8_Br_cAMPS_2, Csub, Csub_H89, AC, RIIp_cAMP_2, pAC_Fsk, RIIp_Rp8_pCPT_cAMPS_C_2, pPDE, RIIp_Rp_cAMPS_C_2, RIIp_2, RIIp_cAMP_C_2
end

function compute_u0(pODEProblem::AbstractVector)::AbstractVector 

	#pODEProblem[1] = ki_Rp8_pCPT_cAMPS_pAB, pODEProblem[2] = xi_b_Rp_cAMPS, pODEProblem[3] = RII2_total, pODEProblem[4] = H89_level, pODEProblem[5] = fourABnOH_incubation_time, pODEProblem[6] = KD_Fsk, pODEProblem[7] = kdeg_cAMP_free, pODEProblem[8] = Rp8_Br_cAMPS_pAB_level, pODEProblem[9] = xi_KD_Rp8_Br_cAMPS, pODEProblem[10] = kf_PDE_Csub, pODEProblem[11] = Sp8_Br_cAMPS_AM_level, pODEProblem[12] = xi_kf_RII_C_2__RII_2, pODEProblem[13] = kf_RIIp_2__RII_2, pODEProblem[14] = kf_cAMP, pODEProblem[15] = IBMX_time, pODEProblem[16] = H89_bool1, pODEProblem[17] = Rp_cAMPS_pAB_incubation_time, pODEProblem[18] = kf_H89, pODEProblem[19] = kf_RII_C_2__RII_2, pODEProblem[20] = kf_RII_C_2__RIIp_C_2, pODEProblem[21] = xi_i_Rp8_pCPT_cAMPS_pAB, pODEProblem[22] = xi_pAC, pODEProblem[23] = fourABnOH_level, pODEProblem[24] = ki_Sp8_Br_cAMPS_AM, pODEProblem[25] = Rp8_pCPT_cAMPS_pAB_bool1, pODEProblem[26] = xi_b_Sp8_Br_cAMPS, pODEProblem[27] = xi_b_Rp8_Br_cAMPS, pODEProblem[28] = H89_time, pODEProblem[29] = Rp_cAMPS_pAB_bool1, pODEProblem[30] = kdeg_cAMP, pODEProblem[31] = xi_AC_cAMP_Fsk, pODEProblem[32] = xi_b_Rp8_pCPT_cAMPS, pODEProblem[33] = ki_Rp_cAMPS_pAB, pODEProblem[34] = xi_i_Rp_cAMPS_pAB, pODEProblem[35] = KD_PDE_Csub, pODEProblem[36] = ki_IBMX, pODEProblem[37] = IBMXex_bool1, pODEProblem[38] = Fsk_time, pODEProblem[39] = Rp8_pCPT_cAMPS_pAB_incubation_time, pODEProblem[40] = PDE_total, pODEProblem[41] = Rp8_Br_cAMPS_pAB_bool1, pODEProblem[42] = ki_Rp8_Br_cAMPS_pAB, pODEProblem[43] = xi_kf_RII_2__RII_C_2, pODEProblem[44] = default, pODEProblem[45] = xi_i_Rp8_Br_cAMPS_pAB, pODEProblem[46] = Sp8_Br_cAMPS_AM_time, pODEProblem[47] = xi_KD_Rp8_pCPT_cAMPS, pODEProblem[48] = kp_AC, pODEProblem[49] = xi_pPDE, pODEProblem[50] = xi_KD_Rp_cAMPS, pODEProblem[51] = Rp8_Br_cAMPS_pAB_incubation_time, pODEProblem[52] = AC_total, pODEProblem[53] = kf_RIIp_C_2__RII_C_2, pODEProblem[54] = xi_KD_Sp8_Br_cAMPS, pODEProblem[55] = nuc, pODEProblem[56] = Sp8_Br_cAMPS_AM_bool1, pODEProblem[57] = kf_RIIp_cAMP_C_2__RIIp_2, pODEProblem[58] = KD_cAMP, pODEProblem[59] = KD_IBMX, pODEProblem[60] = kf_Fsk, pODEProblem[61] = xi_i_Sp8_Br_cAMPS_AM, pODEProblem[62] = cyt, pODEProblem[63] = KD_H89, pODEProblem[64] = ks_AC_cAMP, pODEProblem[65] = fourABnOH_bool1, pODEProblem[66] = Rp8_pCPT_cAMPS_pAB_level, pODEProblem[67] = Rp_cAMPS_pAB_level, pODEProblem[68] = kf_RII_2__RII_C_2, pODEProblem[69] = Fsk_level, pODEProblem[70] = Fsk_bool1, pODEProblem[71] = kdp_AC, pODEProblem[72] = IBMX_level

	pAC = 0.0 
	Rp8_Br_cAMPS = 0.0 
	Rp8_pCPT_cAMPS = 0.0 
	PDE = 1.0 
	Rp_cAMPS = 0.0 
	RII_2 = 0.057671482854616 
	RIIp_Rp8_Br_cAMPS_C_2 = 0.0 
	cAMP = 0.0575218758977949 
	RIIp_Sp8_Br_cAMPS_C_2 = 0.0 
	IBMX = 0.0 
	AC_Fsk = 0.0 
	RIIp_C_2 = 0.287974352643203 
	Sp8_Br_cAMPS = 0.0 
	RII_C_2 = 0.514562213223039 
	RIIp_Sp8_Br_cAMPS_2 = 0.0 
	Csub = 0.192051405433538 
	Csub_H89 = 0.0 
	AC = 1.0 
	RIIp_cAMP_2 = 0.000400416634006059 
	pAC_Fsk = 0.0 
	RIIp_Rp8_pCPT_cAMPS_C_2 = 0.0 
	pPDE = 0.0 
	RIIp_Rp_cAMPS_C_2 = 0.0 
	RIIp_2 = 0.133979505944916 
	RIIp_cAMP_C_2 = 0.00541202870022029 

	 return [pAC, Rp8_Br_cAMPS, Rp8_pCPT_cAMPS, PDE, Rp_cAMPS, RII_2, RIIp_Rp8_Br_cAMPS_C_2, cAMP, RIIp_Sp8_Br_cAMPS_C_2, IBMX, AC_Fsk, RIIp_C_2, Sp8_Br_cAMPS, RII_C_2, RIIp_Sp8_Br_cAMPS_2, Csub, Csub_H89, AC, RIIp_cAMP_2, pAC_Fsk, RIIp_Rp8_pCPT_cAMPS_C_2, pPDE, RIIp_Rp_cAMPS_C_2, RIIp_2, RIIp_cAMP_C_2]
end

function compute_σ(u::AbstractVector, t::Real, θ_sd::AbstractVector, θ_dynamic::AbstractVector, θ_nonDynamic::AbstractVector, 
                   parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :pRII_Microscopy 
		noiseParameter1_pRII_Microscopy = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_pRII_Microscopy * θ_nonDynamic[6] * θ_nonDynamic[8] 
	end

	if observableId == :pRII_Western 
		return θ_nonDynamic[5] * θ_nonDynamic[9] 
	end

	if observableId == :Calpha_Microscopy 
		return θ_nonDynamic[4] * θ_nonDynamic[7] 
	end

end