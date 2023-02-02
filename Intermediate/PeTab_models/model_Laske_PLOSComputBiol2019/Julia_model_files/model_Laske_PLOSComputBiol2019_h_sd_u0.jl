#u[1] = P_PA, u[2] = Vp_cyt, u[3] = Vp_nuc, u[4] = P_NA, u[5] = R_C_RdRp, u[6] = R_M6, u[7] = P_NEP, u[8] = V_end, u[9] = B_att_Hi, u[10] = R_M5, u[11] = Vp_cyt_M1, u[12] = V_rel, u[13] = R_C, u[14] = Vp_nuc_M1, u[15] = Cp, u[16] = V_att_Lo, u[17] = V_att_Hi, u[18] = B_att_Lo, u[19] = P_HA, u[20] = P_M1, u[21] = R_M1, u[22] = P_M2, u[23] = R_M7, u[24] = P_NP, u[25] = V_ex, u[26] = R_M4, u[27] = R_M8, u[28] = R_V_RdRp, u[29] = R_M2, u[30] = R_V, u[31] = P_B2, u[32] = P_RdRp, u[33] = R_M3, u[34] = P_B1
#θ_dynamicNames[1] = k_imp, θ_dynamicNames[2] = k_syn_R_M, θ_dynamicNames[3] = k_syn_R_C, θ_dynamicNames[4] = k_syn_R_V, θ_dynamicNames[5] = k_bind_M1, θ_dynamicNames[6] = k_rel
#

function compute_h(u::AbstractVector, t::Real, θ_dynamic::AbstractVector, θ_observable::AbstractVector,
                   θ_nonDynamic::AbstractVector, parameterInfo::ParametersInfo, observableId::Symbol, 
                      parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :RM5 
		return u[10] 
	end

	if observableId == :RVSegTot 
		observableParameter1_RVSegTot = getObsOrSdParam(θ_observable, parameterMap)
		return observableParameter1_RVSegTot + ( ( 8  *  ( u[17] + u[16] + u[8] )  + u[11] + u[14] + u[2] + u[3] + u[30] )  / 8 ) 
	end

	if observableId == :RCSegTot 
		return ( ( u[15] + u[5] + u[13] )  / 8 ) 
	end

	if observableId == :Vrel 
		return u[12] 
	end

	if observableId == :IntNucOffset 
		observableParameter1_IntNucOffset = getObsOrSdParam(θ_observable, parameterMap)
		return observableParameter1_IntNucOffset 
	end

	if observableId == :FracNucInt_1 
		observableParameter1_FracNucInt_1 = getObsOrSdParam(θ_observable, parameterMap)
		return ( ( u[14] + u[3] )  /  ( ( u[14] + u[3] )  +  ( 8  * u[8] + u[2] + u[11] ) )  * 100 ) + observableParameter1_FracNucInt_1 
	end

	if observableId == :FracNucInt_2 
		observableParameter1_FracNucInt_2 = getObsOrSdParam(θ_observable, parameterMap)
		return ( ( u[14] + u[3] )  /  ( ( u[14] + u[3] )  +  ( 8  * u[8] + u[2] + u[11] ) )  * 100 ) + observableParameter1_FracNucInt_2 
	end

	if observableId == :FracNucInt_3 
		observableParameter1_FracNucInt_3 = getObsOrSdParam(θ_observable, parameterMap)
		return ( ( u[14] + u[3] )  /  ( ( u[14] + u[3] )  +  ( 8  * u[8] + u[2] + u[11] ) )  * 100 ) + observableParameter1_FracNucInt_3 
	end

	if observableId == :FracNucInt_4 
		observableParameter1_FracNucInt_4 = getObsOrSdParam(θ_observable, parameterMap)
		return ( ( u[14] + u[3] )  /  ( ( u[14] + u[3] )  +  ( 8  * u[8] + u[2] + u[11] ) )  * 100 ) + observableParameter1_FracNucInt_4 
	end

	if observableId == :FracNucInt_5 
		observableParameter1_FracNucInt_5 = getObsOrSdParam(θ_observable, parameterMap)
		return ( ( u[14] + u[3] )  /  ( ( u[14] + u[3] )  +  ( 8  * u[8] + u[2] + u[11] ) )  * 100 ) + observableParameter1_FracNucInt_5 
	end

	if observableId == :FracNucInt_6 
		observableParameter1_FracNucInt_6 = getObsOrSdParam(θ_observable, parameterMap)
		return ( ( u[14] + u[3] )  /  ( ( u[14] + u[3] )  +  ( 8  * u[8] + u[2] + u[11] ) )  * 100 ) + observableParameter1_FracNucInt_6 
	end

	if observableId == :FracNucInt_7 
		observableParameter1_FracNucInt_7 = getObsOrSdParam(θ_observable, parameterMap)
		return ( ( u[14] + u[3] )  /  ( ( u[14] + u[3] )  +  ( 8  * u[8] + u[2] + u[11] ) )  * 100 ) + observableParameter1_FracNucInt_7 
	end

	if observableId == :FracNucInt_8 
		observableParameter1_FracNucInt_8 = getObsOrSdParam(θ_observable, parameterMap)
		return ( ( u[14] + u[3] )  /  ( ( u[14] + u[3] )  +  ( 8  * u[8] + u[2] + u[11] ) )  * 100 ) + observableParameter1_FracNucInt_8 
	end

end

function compute_u0!(u0::AbstractVector, pODEProblem::AbstractVector) 

	#pODEProblem[1] = k_deg_R_RdRp, pODEProblem[2] = ModelValue_104, pODEProblem[3] = ModelValue_63, pODEProblem[4] = ModelValue_114, pODEProblem[5] = ModelValue_108, pODEProblem[6] = ModelValue_64, pODEProblem[7] = L6, pODEProblem[8] = k_end, pODEProblem[9] = N_P_RdRp, pODEProblem[10] = D_rib, pODEProblem[11] = ModelValue_111, pODEProblem[12] = ModelValue_107, pODEProblem[13] = k_rel, pODEProblem[14] = ModelValue_105, pODEProblem[15] = compartment, pODEProblem[16] = L1, pODEProblem[17] = K_eq_Hi, pODEProblem[18] = N_P_NP, pODEProblem[19] = L2, pODEProblem[20] = ModelValue_113, pODEProblem[21] = k_bind_M1, pODEProblem[22] = ModelValue_116, pODEProblem[23] = F_Spl7, pODEProblem[24] = k_deg_Rnp, pODEProblem[25] = N_P_NA, pODEProblem[26] = K_eq_Lo, pODEProblem[27] = ModelValue_101, pODEProblem[28] = k_exp_Vp_nuc_M1, pODEProblem[29] = ModelValue_79, pODEProblem[30] = ModelValue_69, pODEProblem[31] = ModelValue_90, pODEProblem[32] = ModelValue_91, pODEProblem[33] = ModelValue_103, pODEProblem[34] = N_P_M1, pODEProblem[35] = k_syn_R_M, pODEProblem[36] = L8, pODEProblem[37] = ModelValue_84, pODEProblem[38] = ModelValue_85, pODEProblem[39] = ModelValue_115, pODEProblem[40] = k_att_Lo, pODEProblem[41] = k_fus, pODEProblem[42] = N_P_HA, pODEProblem[43] = k_deg_R, pODEProblem[44] = N_P_NEP, pODEProblem[45] = k_bind_RdRp, pODEProblem[46] = ModelValue_80, pODEProblem[47] = ModelValue_89, pODEProblem[48] = k_RdRp, pODEProblem[49] = L5, pODEProblem[50] = N_P_M2, pODEProblem[51] = k_syn_R_V, pODEProblem[52] = L4, pODEProblem[53] = ModelValue_87, pODEProblem[54] = k_bind_NP, pODEProblem[55] = L3, pODEProblem[56] = k_imp, pODEProblem[57] = ModelValue_106, pODEProblem[58] = F_Spl8, pODEProblem[59] = k_deg_R_M, pODEProblem[60] = k_syn_P, pODEProblem[61] = ModelValue_82, pODEProblem[62] = k_syn_R_C, pODEProblem[63] = L7, pODEProblem[64] = K_V_rel, pODEProblem[65] = ModelValue_86, pODEProblem[66] = ModelValue_102, pODEProblem[67] = F_fus, pODEProblem[68] = k_att_Hi, pODEProblem[69] = ModelValue_88

	P_PA = 0.0 
	Vp_cyt = 0.0 
	Vp_nuc = 0.0 
	P_NA = 0.0 
	R_C_RdRp = 0.0 
	R_M6 = 0.0 
	P_NEP = 0.0 
	V_end = 0.0 
	B_att_Hi = 150.0 
	R_M5 = 0.0 
	Vp_cyt_M1 = 0.0 
	V_rel = 0.0 
	R_C = 0.0 
	Vp_nuc_M1 = 0.0 
	Cp = 0.0 
	V_att_Lo = 0.0 
	V_att_Hi = 0.0 
	B_att_Lo = 1000.0 
	P_HA = 0.0 
	P_M1 = 0.0 
	R_M1 = 0.0 
	P_M2 = 0.0 
	R_M7 = 0.0 
	P_NP = 0.0 
	V_ex = 50.0 
	R_M4 = 0.0 
	R_M8 = 0.0 
	R_V_RdRp = 0.0 
	R_M2 = 0.0 
	R_V = 0.0 
	P_B2 = 0.0 
	P_RdRp = 0.0 
	R_M3 = 0.0 
	P_B1 = 0.0 

	u0 .= P_PA, Vp_cyt, Vp_nuc, P_NA, R_C_RdRp, R_M6, P_NEP, V_end, B_att_Hi, R_M5, Vp_cyt_M1, V_rel, R_C, Vp_nuc_M1, Cp, V_att_Lo, V_att_Hi, B_att_Lo, P_HA, P_M1, R_M1, P_M2, R_M7, P_NP, V_ex, R_M4, R_M8, R_V_RdRp, R_M2, R_V, P_B2, P_RdRp, R_M3, P_B1
end

function compute_u0(pODEProblem::AbstractVector)::AbstractVector 

	#pODEProblem[1] = k_deg_R_RdRp, pODEProblem[2] = ModelValue_104, pODEProblem[3] = ModelValue_63, pODEProblem[4] = ModelValue_114, pODEProblem[5] = ModelValue_108, pODEProblem[6] = ModelValue_64, pODEProblem[7] = L6, pODEProblem[8] = k_end, pODEProblem[9] = N_P_RdRp, pODEProblem[10] = D_rib, pODEProblem[11] = ModelValue_111, pODEProblem[12] = ModelValue_107, pODEProblem[13] = k_rel, pODEProblem[14] = ModelValue_105, pODEProblem[15] = compartment, pODEProblem[16] = L1, pODEProblem[17] = K_eq_Hi, pODEProblem[18] = N_P_NP, pODEProblem[19] = L2, pODEProblem[20] = ModelValue_113, pODEProblem[21] = k_bind_M1, pODEProblem[22] = ModelValue_116, pODEProblem[23] = F_Spl7, pODEProblem[24] = k_deg_Rnp, pODEProblem[25] = N_P_NA, pODEProblem[26] = K_eq_Lo, pODEProblem[27] = ModelValue_101, pODEProblem[28] = k_exp_Vp_nuc_M1, pODEProblem[29] = ModelValue_79, pODEProblem[30] = ModelValue_69, pODEProblem[31] = ModelValue_90, pODEProblem[32] = ModelValue_91, pODEProblem[33] = ModelValue_103, pODEProblem[34] = N_P_M1, pODEProblem[35] = k_syn_R_M, pODEProblem[36] = L8, pODEProblem[37] = ModelValue_84, pODEProblem[38] = ModelValue_85, pODEProblem[39] = ModelValue_115, pODEProblem[40] = k_att_Lo, pODEProblem[41] = k_fus, pODEProblem[42] = N_P_HA, pODEProblem[43] = k_deg_R, pODEProblem[44] = N_P_NEP, pODEProblem[45] = k_bind_RdRp, pODEProblem[46] = ModelValue_80, pODEProblem[47] = ModelValue_89, pODEProblem[48] = k_RdRp, pODEProblem[49] = L5, pODEProblem[50] = N_P_M2, pODEProblem[51] = k_syn_R_V, pODEProblem[52] = L4, pODEProblem[53] = ModelValue_87, pODEProblem[54] = k_bind_NP, pODEProblem[55] = L3, pODEProblem[56] = k_imp, pODEProblem[57] = ModelValue_106, pODEProblem[58] = F_Spl8, pODEProblem[59] = k_deg_R_M, pODEProblem[60] = k_syn_P, pODEProblem[61] = ModelValue_82, pODEProblem[62] = k_syn_R_C, pODEProblem[63] = L7, pODEProblem[64] = K_V_rel, pODEProblem[65] = ModelValue_86, pODEProblem[66] = ModelValue_102, pODEProblem[67] = F_fus, pODEProblem[68] = k_att_Hi, pODEProblem[69] = ModelValue_88

	P_PA = 0.0 
	Vp_cyt = 0.0 
	Vp_nuc = 0.0 
	P_NA = 0.0 
	R_C_RdRp = 0.0 
	R_M6 = 0.0 
	P_NEP = 0.0 
	V_end = 0.0 
	B_att_Hi = 150.0 
	R_M5 = 0.0 
	Vp_cyt_M1 = 0.0 
	V_rel = 0.0 
	R_C = 0.0 
	Vp_nuc_M1 = 0.0 
	Cp = 0.0 
	V_att_Lo = 0.0 
	V_att_Hi = 0.0 
	B_att_Lo = 1000.0 
	P_HA = 0.0 
	P_M1 = 0.0 
	R_M1 = 0.0 
	P_M2 = 0.0 
	R_M7 = 0.0 
	P_NP = 0.0 
	V_ex = 50.0 
	R_M4 = 0.0 
	R_M8 = 0.0 
	R_V_RdRp = 0.0 
	R_M2 = 0.0 
	R_V = 0.0 
	P_B2 = 0.0 
	P_RdRp = 0.0 
	R_M3 = 0.0 
	P_B1 = 0.0 

	 return [P_PA, Vp_cyt, Vp_nuc, P_NA, R_C_RdRp, R_M6, P_NEP, V_end, B_att_Hi, R_M5, Vp_cyt_M1, V_rel, R_C, Vp_nuc_M1, Cp, V_att_Lo, V_att_Hi, B_att_Lo, P_HA, P_M1, R_M1, P_M2, R_M7, P_NP, V_ex, R_M4, R_M8, R_V_RdRp, R_M2, R_V, P_B2, P_RdRp, R_M3, P_B1]
end

function compute_σ(u::AbstractVector, t::Real, θ_sd::AbstractVector, θ_dynamic::AbstractVector, θ_nonDynamic::AbstractVector, 
                   parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :RM5 
		noiseParameter1_RM5 = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_RM5 
	end

	if observableId == :RVSegTot 
		noiseParameter1_RVSegTot = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_RVSegTot 
	end

	if observableId == :RCSegTot 
		noiseParameter1_RCSegTot = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_RCSegTot 
	end

	if observableId == :Vrel 
		noiseParameter1_Vrel = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_Vrel 
	end

	if observableId == :IntNucOffset 
		noiseParameter1_IntNucOffset = getObsOrSdParam(θ_sd, parameterMap)
		return 13.0334748500416 * noiseParameter1_IntNucOffset 
	end

	if observableId == :FracNucInt_1 
		noiseParameter1_FracNucInt_1 = getObsOrSdParam(θ_sd, parameterMap)
		return 3.10956052629092 * noiseParameter1_FracNucInt_1 
	end

	if observableId == :FracNucInt_2 
		noiseParameter1_FracNucInt_2 = getObsOrSdParam(θ_sd, parameterMap)
		return 3.44650910342625 * noiseParameter1_FracNucInt_2 
	end

	if observableId == :FracNucInt_3 
		noiseParameter1_FracNucInt_3 = getObsOrSdParam(θ_sd, parameterMap)
		return 1.96540708251497 * noiseParameter1_FracNucInt_3 
	end

	if observableId == :FracNucInt_4 
		noiseParameter1_FracNucInt_4 = getObsOrSdParam(θ_sd, parameterMap)
		return 3.03676774438437 * noiseParameter1_FracNucInt_4 
	end

	if observableId == :FracNucInt_5 
		noiseParameter1_FracNucInt_5 = getObsOrSdParam(θ_sd, parameterMap)
		return 4.00311961683219 * noiseParameter1_FracNucInt_5 
	end

	if observableId == :FracNucInt_6 
		noiseParameter1_FracNucInt_6 = getObsOrSdParam(θ_sd, parameterMap)
		return 4.8783911282307 * noiseParameter1_FracNucInt_6 
	end

	if observableId == :FracNucInt_7 
		noiseParameter1_FracNucInt_7 = getObsOrSdParam(θ_sd, parameterMap)
		return 3.96353798182046 * noiseParameter1_FracNucInt_7 
	end

	if observableId == :FracNucInt_8 
		noiseParameter1_FracNucInt_8 = getObsOrSdParam(θ_sd, parameterMap)
		return 4.8820897165046 * noiseParameter1_FracNucInt_8 
	end

end