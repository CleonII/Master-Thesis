#u[1] = CumulativeDiagnosed, u[2] = Infected, u[3] = Extinct, u[4] = Diagnosed, u[5] = DiagnosedHealed, u[6] = Ailing, u[7] = Healed, u[8] = Recognized, u[9] = Susceptible, u[10] = Threatened
#θ_dynamicNames[1] = alpha_0, θ_dynamicNames[2] = beta_0, θ_dynamicNames[3] = delta_0, θ_dynamicNames[4] = epsilon_0, θ_dynamicNames[5] = eta_0, θ_dynamicNames[6] = gamma_0, θ_dynamicNames[7] = kappa_0, θ_dynamicNames[8] = lam_0, θ_dynamicNames[9] = mu_0, θ_dynamicNames[10] = nu_0, θ_dynamicNames[11] = rho_0, θ_dynamicNames[12] = sigma_0, θ_dynamicNames[13] = tau, θ_dynamicNames[14] = theta, θ_dynamicNames[15] = xi_0, θ_dynamicNames[16] = zeta_0, θ_dynamicNames[17] = alpha_4, θ_dynamicNames[18] = alpha_22, θ_dynamicNames[19] = alpha_28, θ_dynamicNames[20] = beta_4, θ_dynamicNames[21] = beta_22, θ_dynamicNames[22] = delta_4, θ_dynamicNames[23] = delta_22, θ_dynamicNames[24] = epsilon_12, θ_dynamicNames[25] = epsilon_38, θ_dynamicNames[26] = eta_22, θ_dynamicNames[27] = eta_38, θ_dynamicNames[28] = gamma_4, θ_dynamicNames[29] = gamma_22, θ_dynamicNames[30] = gamma_28, θ_dynamicNames[31] = kappa_22, θ_dynamicNames[32] = kappa_38, θ_dynamicNames[33] = lam_22, θ_dynamicNames[34] = mu_22, θ_dynamicNames[35] = nu_22, θ_dynamicNames[36] = rho_22, θ_dynamicNames[37] = rho_38, θ_dynamicNames[38] = sigma_22, θ_dynamicNames[39] = sigma_38, θ_dynamicNames[40] = xi_22, θ_dynamicNames[41] = xi_38, θ_dynamicNames[42] = zeta_22, θ_dynamicNames[43] = zeta_38
#

function compute_h(u::AbstractVector, t::Real, θ_dynamic::AbstractVector, θ_observable::AbstractVector,
                   θ_nonDynamic::AbstractVector, parameterInfo::ParametersInfo, observableId::Symbol, 
                      parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :observable_TotalCases 
		return u[1] 
	end

	if observableId == :observable_Deaths 
		return u[3] 
	end

	if observableId == :observable_Recovered 
		return u[5] 
	end

	if observableId == :observable_CurrentCases 
		return ( u[4] + u[8] + u[10] ) 
	end

	if observableId == :observable_DiagHome 
		return u[4] 
	end

	if observableId == :observable_Hospit 
		return u[8] 
	end

	if observableId == :observable_ICU 
		return u[10] 
	end

end

function compute_u0!(u0::AbstractVector, pODEProblem::AbstractVector) 

	#pODEProblem[1] = beta_bool3, pODEProblem[2] = alpha_bool12, pODEProblem[3] = rho_bool2, pODEProblem[4] = rho_bool3, pODEProblem[5] = kappa_bool7, pODEProblem[6] = rho_38, pODEProblem[7] = nu_bool4, pODEProblem[8] = kappa_0, pODEProblem[9] = alpha_bool6, pODEProblem[10] = zeta_bool7, pODEProblem[11] = rho_bool7, pODEProblem[12] = zeta_bool5, pODEProblem[13] = eta_bool2, pODEProblem[14] = rho_bool10, pODEProblem[15] = lam_bool6, pODEProblem[16] = xi_bool10, pODEProblem[17] = xi_38, pODEProblem[18] = delta_bool2, pODEProblem[19] = sigma_38, pODEProblem[20] = zeta_bool6, pODEProblem[21] = kappa_bool3, pODEProblem[22] = gamma_bool8, pODEProblem[23] = kappa_bool6, pODEProblem[24] = beta_bool8, pODEProblem[25] = zeta_bool9, pODEProblem[26] = kappa_bool8, pODEProblem[27] = delta_bool7, pODEProblem[28] = epsilon_bool8, pODEProblem[29] = nu_bool5, pODEProblem[30] = tau, pODEProblem[31] = alpha_bool8, pODEProblem[32] = zeta_bool4, pODEProblem[33] = beta_bool9, pODEProblem[34] = theta, pODEProblem[35] = xi_bool8, pODEProblem[36] = zeta_22, pODEProblem[37] = sigma_bool7, pODEProblem[38] = alpha_bool1, pODEProblem[39] = mu_bool4, pODEProblem[40] = nu_bool1, pODEProblem[41] = eta_0, pODEProblem[42] = zeta_38, pODEProblem[43] = xi_bool2, pODEProblem[44] = rho_bool5, pODEProblem[45] = alpha_50, pODEProblem[46] = beta_0, pODEProblem[47] = kappa_bool2, pODEProblem[48] = alpha_bool3, pODEProblem[49] = zeta_bool1, pODEProblem[50] = eta_bool10, pODEProblem[51] = beta_bool6, pODEProblem[52] = delta_0, pODEProblem[53] = beta_bool1, pODEProblem[54] = xi_bool6, pODEProblem[55] = gamma_22, pODEProblem[56] = gamma_bool9, pODEProblem[57] = lam_bool3, pODEProblem[58] = xi_bool1, pODEProblem[59] = sigma_0, pODEProblem[60] = beta_22, pODEProblem[61] = zeta_bool8, pODEProblem[62] = lam_bool4, pODEProblem[63] = gamma_28, pODEProblem[64] = alpha_0, pODEProblem[65] = mu_22, pODEProblem[66] = eta_bool6, pODEProblem[67] = gamma_bool4, pODEProblem[68] = gamma_bool11, pODEProblem[69] = epsilon_bool9, pODEProblem[70] = delta_bool5, pODEProblem[71] = gamma_bool1, pODEProblem[72] = kappa_bool5, pODEProblem[73] = zeta_bool10, pODEProblem[74] = nu_bool6, pODEProblem[75] = lam_bool5, pODEProblem[76] = eta_bool1, pODEProblem[77] = mu_50, pODEProblem[78] = lam_50, pODEProblem[79] = kappa_38, pODEProblem[80] = alpha_bool14, pODEProblem[81] = alpha_28, pODEProblem[82] = gamma_bool10, pODEProblem[83] = rho_0, pODEProblem[84] = sigma_bool8, pODEProblem[85] = gamma_0, pODEProblem[86] = epsilon_bool10, pODEProblem[87] = gamma_bool6, pODEProblem[88] = alpha_4, pODEProblem[89] = sigma_50, pODEProblem[90] = lam_0, pODEProblem[91] = beta_50, pODEProblem[92] = epsilon_bool2, pODEProblem[93] = rho_bool4, pODEProblem[94] = mu_bool5, pODEProblem[95] = sigma_bool9, pODEProblem[96] = delta_bool10, pODEProblem[97] = alpha_bool4, pODEProblem[98] = nu_0, pODEProblem[99] = rho_bool1, pODEProblem[100] = mu_bool2, pODEProblem[101] = delta_bool3, pODEProblem[102] = xi_22, pODEProblem[103] = xi_bool5, pODEProblem[104] = sigma_bool2, pODEProblem[105] = beta_bool2, pODEProblem[106] = nu_22, pODEProblem[107] = beta_bool4, pODEProblem[108] = eta_22, pODEProblem[109] = xi_50, pODEProblem[110] = lam_22, pODEProblem[111] = eta_bool7, pODEProblem[112] = beta_bool7, pODEProblem[113] = mu_bool1, pODEProblem[114] = sigma_bool3, pODEProblem[115] = epsilon_bool3, pODEProblem[116] = epsilon_38, pODEProblem[117] = lam_bool1, pODEProblem[118] = alpha_bool11, pODEProblem[119] = epsilon_0, pODEProblem[120] = sigma_bool5, pODEProblem[121] = Italy, pODEProblem[122] = rho_bool9, pODEProblem[123] = kappa_50, pODEProblem[124] = xi_bool3, pODEProblem[125] = mu_0, pODEProblem[126] = lam_bool2, pODEProblem[127] = eta_bool9, pODEProblem[128] = epsilon_bool5, pODEProblem[129] = epsilon_bool7, pODEProblem[130] = gamma_bool13, pODEProblem[131] = kappa_22, pODEProblem[132] = gamma_bool2, pODEProblem[133] = nu_bool3, pODEProblem[134] = epsilon_bool6, pODEProblem[135] = delta_bool8, pODEProblem[136] = nu_bool2, pODEProblem[137] = alpha_22, pODEProblem[138] = zeta_bool3, pODEProblem[139] = sigma_bool6, pODEProblem[140] = alpha_bool2, pODEProblem[141] = beta_bool5, pODEProblem[142] = delta_bool6, pODEProblem[143] = epsilon_bool1, pODEProblem[144] = eta_bool5, pODEProblem[145] = rho_22, pODEProblem[146] = eta_bool8, pODEProblem[147] = kappa_bool4, pODEProblem[148] = rho_bool6, pODEProblem[149] = delta_4, pODEProblem[150] = delta_bool9, pODEProblem[151] = xi_bool7, pODEProblem[152] = zeta_50, pODEProblem[153] = epsilon_50, pODEProblem[154] = xi_bool9, pODEProblem[155] = xi_0, pODEProblem[156] = delta_bool1, pODEProblem[157] = xi_bool4, pODEProblem[158] = gamma_bool5, pODEProblem[159] = beta_4, pODEProblem[160] = initalTimeManual, pODEProblem[161] = eta_bool3, pODEProblem[162] = zeta_bool2, pODEProblem[163] = kappa_bool10, pODEProblem[164] = epsilon_12, pODEProblem[165] = mu_bool6, pODEProblem[166] = gamma_bool14, pODEProblem[167] = gamma_bool3, pODEProblem[168] = delta_50, pODEProblem[169] = eta_38, pODEProblem[170] = eta_bool4, pODEProblem[171] = mu_bool3, pODEProblem[172] = epsilon_bool4, pODEProblem[173] = gamma_50, pODEProblem[174] = gamma_bool7, pODEProblem[175] = alpha_bool10, pODEProblem[176] = kappa_bool9, pODEProblem[177] = kappa_bool1, pODEProblem[178] = sigma_22, pODEProblem[179] = rho_bool8, pODEProblem[180] = eta_50, pODEProblem[181] = delta_22, pODEProblem[182] = alpha_bool9, pODEProblem[183] = zeta_0, pODEProblem[184] = sigma_bool1, pODEProblem[185] = beta_bool10, pODEProblem[186] = rho_50, pODEProblem[187] = alpha_bool7, pODEProblem[188] = nu_50, pODEProblem[189] = gamma_4, pODEProblem[190] = delta_bool4, pODEProblem[191] = gamma_bool12, pODEProblem[192] = sigma_bool10, pODEProblem[193] = sigma_bool4, pODEProblem[194] = alpha_bool5, pODEProblem[195] = alpha_bool13

	CumulativeDiagnosed = 3.6666666630000003e-7 
	Infected = 3.33333333e-6 
	Extinct = 0.0 
	Diagnosed = 3.33333333e-7 
	DiagnosedHealed = 0.0 
	Ailing = 1.66666666e-8 
	Healed = 0.0 
	Recognized = 3.33333333e-8 
	Susceptible = 0.9999963 
	Threatened = 0.0 

	u0 .= CumulativeDiagnosed, Infected, Extinct, Diagnosed, DiagnosedHealed, Ailing, Healed, Recognized, Susceptible, Threatened
end

function compute_u0(pODEProblem::AbstractVector)::AbstractVector 

	#pODEProblem[1] = beta_bool3, pODEProblem[2] = alpha_bool12, pODEProblem[3] = rho_bool2, pODEProblem[4] = rho_bool3, pODEProblem[5] = kappa_bool7, pODEProblem[6] = rho_38, pODEProblem[7] = nu_bool4, pODEProblem[8] = kappa_0, pODEProblem[9] = alpha_bool6, pODEProblem[10] = zeta_bool7, pODEProblem[11] = rho_bool7, pODEProblem[12] = zeta_bool5, pODEProblem[13] = eta_bool2, pODEProblem[14] = rho_bool10, pODEProblem[15] = lam_bool6, pODEProblem[16] = xi_bool10, pODEProblem[17] = xi_38, pODEProblem[18] = delta_bool2, pODEProblem[19] = sigma_38, pODEProblem[20] = zeta_bool6, pODEProblem[21] = kappa_bool3, pODEProblem[22] = gamma_bool8, pODEProblem[23] = kappa_bool6, pODEProblem[24] = beta_bool8, pODEProblem[25] = zeta_bool9, pODEProblem[26] = kappa_bool8, pODEProblem[27] = delta_bool7, pODEProblem[28] = epsilon_bool8, pODEProblem[29] = nu_bool5, pODEProblem[30] = tau, pODEProblem[31] = alpha_bool8, pODEProblem[32] = zeta_bool4, pODEProblem[33] = beta_bool9, pODEProblem[34] = theta, pODEProblem[35] = xi_bool8, pODEProblem[36] = zeta_22, pODEProblem[37] = sigma_bool7, pODEProblem[38] = alpha_bool1, pODEProblem[39] = mu_bool4, pODEProblem[40] = nu_bool1, pODEProblem[41] = eta_0, pODEProblem[42] = zeta_38, pODEProblem[43] = xi_bool2, pODEProblem[44] = rho_bool5, pODEProblem[45] = alpha_50, pODEProblem[46] = beta_0, pODEProblem[47] = kappa_bool2, pODEProblem[48] = alpha_bool3, pODEProblem[49] = zeta_bool1, pODEProblem[50] = eta_bool10, pODEProblem[51] = beta_bool6, pODEProblem[52] = delta_0, pODEProblem[53] = beta_bool1, pODEProblem[54] = xi_bool6, pODEProblem[55] = gamma_22, pODEProblem[56] = gamma_bool9, pODEProblem[57] = lam_bool3, pODEProblem[58] = xi_bool1, pODEProblem[59] = sigma_0, pODEProblem[60] = beta_22, pODEProblem[61] = zeta_bool8, pODEProblem[62] = lam_bool4, pODEProblem[63] = gamma_28, pODEProblem[64] = alpha_0, pODEProblem[65] = mu_22, pODEProblem[66] = eta_bool6, pODEProblem[67] = gamma_bool4, pODEProblem[68] = gamma_bool11, pODEProblem[69] = epsilon_bool9, pODEProblem[70] = delta_bool5, pODEProblem[71] = gamma_bool1, pODEProblem[72] = kappa_bool5, pODEProblem[73] = zeta_bool10, pODEProblem[74] = nu_bool6, pODEProblem[75] = lam_bool5, pODEProblem[76] = eta_bool1, pODEProblem[77] = mu_50, pODEProblem[78] = lam_50, pODEProblem[79] = kappa_38, pODEProblem[80] = alpha_bool14, pODEProblem[81] = alpha_28, pODEProblem[82] = gamma_bool10, pODEProblem[83] = rho_0, pODEProblem[84] = sigma_bool8, pODEProblem[85] = gamma_0, pODEProblem[86] = epsilon_bool10, pODEProblem[87] = gamma_bool6, pODEProblem[88] = alpha_4, pODEProblem[89] = sigma_50, pODEProblem[90] = lam_0, pODEProblem[91] = beta_50, pODEProblem[92] = epsilon_bool2, pODEProblem[93] = rho_bool4, pODEProblem[94] = mu_bool5, pODEProblem[95] = sigma_bool9, pODEProblem[96] = delta_bool10, pODEProblem[97] = alpha_bool4, pODEProblem[98] = nu_0, pODEProblem[99] = rho_bool1, pODEProblem[100] = mu_bool2, pODEProblem[101] = delta_bool3, pODEProblem[102] = xi_22, pODEProblem[103] = xi_bool5, pODEProblem[104] = sigma_bool2, pODEProblem[105] = beta_bool2, pODEProblem[106] = nu_22, pODEProblem[107] = beta_bool4, pODEProblem[108] = eta_22, pODEProblem[109] = xi_50, pODEProblem[110] = lam_22, pODEProblem[111] = eta_bool7, pODEProblem[112] = beta_bool7, pODEProblem[113] = mu_bool1, pODEProblem[114] = sigma_bool3, pODEProblem[115] = epsilon_bool3, pODEProblem[116] = epsilon_38, pODEProblem[117] = lam_bool1, pODEProblem[118] = alpha_bool11, pODEProblem[119] = epsilon_0, pODEProblem[120] = sigma_bool5, pODEProblem[121] = Italy, pODEProblem[122] = rho_bool9, pODEProblem[123] = kappa_50, pODEProblem[124] = xi_bool3, pODEProblem[125] = mu_0, pODEProblem[126] = lam_bool2, pODEProblem[127] = eta_bool9, pODEProblem[128] = epsilon_bool5, pODEProblem[129] = epsilon_bool7, pODEProblem[130] = gamma_bool13, pODEProblem[131] = kappa_22, pODEProblem[132] = gamma_bool2, pODEProblem[133] = nu_bool3, pODEProblem[134] = epsilon_bool6, pODEProblem[135] = delta_bool8, pODEProblem[136] = nu_bool2, pODEProblem[137] = alpha_22, pODEProblem[138] = zeta_bool3, pODEProblem[139] = sigma_bool6, pODEProblem[140] = alpha_bool2, pODEProblem[141] = beta_bool5, pODEProblem[142] = delta_bool6, pODEProblem[143] = epsilon_bool1, pODEProblem[144] = eta_bool5, pODEProblem[145] = rho_22, pODEProblem[146] = eta_bool8, pODEProblem[147] = kappa_bool4, pODEProblem[148] = rho_bool6, pODEProblem[149] = delta_4, pODEProblem[150] = delta_bool9, pODEProblem[151] = xi_bool7, pODEProblem[152] = zeta_50, pODEProblem[153] = epsilon_50, pODEProblem[154] = xi_bool9, pODEProblem[155] = xi_0, pODEProblem[156] = delta_bool1, pODEProblem[157] = xi_bool4, pODEProblem[158] = gamma_bool5, pODEProblem[159] = beta_4, pODEProblem[160] = initalTimeManual, pODEProblem[161] = eta_bool3, pODEProblem[162] = zeta_bool2, pODEProblem[163] = kappa_bool10, pODEProblem[164] = epsilon_12, pODEProblem[165] = mu_bool6, pODEProblem[166] = gamma_bool14, pODEProblem[167] = gamma_bool3, pODEProblem[168] = delta_50, pODEProblem[169] = eta_38, pODEProblem[170] = eta_bool4, pODEProblem[171] = mu_bool3, pODEProblem[172] = epsilon_bool4, pODEProblem[173] = gamma_50, pODEProblem[174] = gamma_bool7, pODEProblem[175] = alpha_bool10, pODEProblem[176] = kappa_bool9, pODEProblem[177] = kappa_bool1, pODEProblem[178] = sigma_22, pODEProblem[179] = rho_bool8, pODEProblem[180] = eta_50, pODEProblem[181] = delta_22, pODEProblem[182] = alpha_bool9, pODEProblem[183] = zeta_0, pODEProblem[184] = sigma_bool1, pODEProblem[185] = beta_bool10, pODEProblem[186] = rho_50, pODEProblem[187] = alpha_bool7, pODEProblem[188] = nu_50, pODEProblem[189] = gamma_4, pODEProblem[190] = delta_bool4, pODEProblem[191] = gamma_bool12, pODEProblem[192] = sigma_bool10, pODEProblem[193] = sigma_bool4, pODEProblem[194] = alpha_bool5, pODEProblem[195] = alpha_bool13

	CumulativeDiagnosed = 3.6666666630000003e-7 
	Infected = 3.33333333e-6 
	Extinct = 0.0 
	Diagnosed = 3.33333333e-7 
	DiagnosedHealed = 0.0 
	Ailing = 1.66666666e-8 
	Healed = 0.0 
	Recognized = 3.33333333e-8 
	Susceptible = 0.9999963 
	Threatened = 0.0 

	 return [CumulativeDiagnosed, Infected, Extinct, Diagnosed, DiagnosedHealed, Ailing, Healed, Recognized, Susceptible, Threatened]
end

function compute_σ(u::AbstractVector, t::Real, θ_sd::AbstractVector, θ_dynamic::AbstractVector, θ_nonDynamic::AbstractVector, 
                   parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :observable_TotalCases 
		noiseParameter1_observable_TotalCases = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_TotalCases 
	end

	if observableId == :observable_Deaths 
		noiseParameter1_observable_Deaths = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_Deaths 
	end

	if observableId == :observable_Recovered 
		noiseParameter1_observable_Recovered = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_Recovered 
	end

	if observableId == :observable_CurrentCases 
		noiseParameter1_observable_CurrentCases = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_CurrentCases 
	end

	if observableId == :observable_DiagHome 
		noiseParameter1_observable_DiagHome = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_DiagHome 
	end

	if observableId == :observable_Hospit 
		noiseParameter1_observable_Hospit = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_Hospit 
	end

	if observableId == :observable_ICU 
		noiseParameter1_observable_ICU = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_ICU 
	end

end