function evalYmod(u, t, dynPar, obsPar, nonDynParam, paramData, observableId, mapObsParam) 

	CumulativeDiagnosed, Infected, Extinct, Diagnosed, DiagnosedHealed, Ailing, Healed, Recognized, Susceptible, Threatened= u 
	alpha_0, beta_0, delta_0, epsilon_0, eta_0, gamma_0, kappa_0, lam_0, mu_0, nu_0, rho_0, sigma_0, tau, theta, xi_0, zeta_0, alpha_4, alpha_22, alpha_28, beta_4, beta_22, delta_4, delta_22, epsilon_12, epsilon_38, eta_22, eta_38, gamma_4, gamma_22, gamma_28, kappa_22, kappa_38, lam_22, mu_22, nu_22, rho_22, rho_38, sigma_22, sigma_38, xi_22, xi_38, zeta_22, zeta_38 = dynPar 

	CurrentDiagnosedInfected = Diagnosed + Recognized + Threatened

	if observableId == "observable_TotalCases" 
		return CumulativeDiagnosed 
	end

	if observableId == "observable_Deaths" 
		return Extinct 
	end

	if observableId == "observable_Recovered" 
		return DiagnosedHealed 
	end

	if observableId == "observable_CurrentCases" 
		return CurrentDiagnosedInfected 
	end

	if observableId == "observable_DiagHome" 
		return Diagnosed 
	end

	if observableId == "observable_Hospit" 
		return Recognized 
	end

	if observableId == "observable_ICU" 
		return Threatened 
	end

end

function evalU0!(u0Vec, paramVec) 

	beta_bool3, alpha_bool12, rho_bool2, rho_bool3, kappa_bool7, rho_38, nu_bool4, kappa_0, alpha_bool6, zeta_bool7, rho_bool7, zeta_bool5, eta_bool2, rho_bool10, lam_bool6, xi_bool10, xi_38, delta_bool2, sigma_38, zeta_bool6, kappa_bool3, gamma_bool8, kappa_bool6, beta_bool8, zeta_bool9, kappa_bool8, delta_bool7, epsilon_bool8, nu_bool5, tau, alpha_bool8, zeta_bool4, beta_bool9, theta, xi_bool8, zeta_22, sigma_bool7, alpha_bool1, mu_bool4, nu_bool1, eta_0, zeta_38, xi_bool2, rho_bool5, alpha_50, beta_0, kappa_bool2, alpha_bool3, zeta_bool1, eta_bool10, beta_bool6, delta_0, beta_bool1, xi_bool6, gamma_22, gamma_bool9, lam_bool3, xi_bool1, sigma_0, beta_22, zeta_bool8, lam_bool4, gamma_28, alpha_0, mu_22, eta_bool6, gamma_bool4, gamma_bool11, epsilon_bool9, delta_bool5, gamma_bool1, kappa_bool5, zeta_bool10, nu_bool6, lam_bool5, eta_bool1, mu_50, lam_50, kappa_38, alpha_bool14, alpha_28, gamma_bool10, rho_0, sigma_bool8, gamma_0, epsilon_bool10, gamma_bool6, alpha_4, sigma_50, lam_0, beta_50, epsilon_bool2, rho_bool4, mu_bool5, sigma_bool9, delta_bool10, alpha_bool4, nu_0, rho_bool1, mu_bool2, delta_bool3, xi_22, xi_bool5, sigma_bool2, beta_bool2, nu_22, beta_bool4, eta_22, xi_50, lam_22, eta_bool7, beta_bool7, mu_bool1, sigma_bool3, epsilon_bool3, epsilon_38, lam_bool1, alpha_bool11, epsilon_0, sigma_bool5, Italy, rho_bool9, kappa_50, xi_bool3, mu_0, lam_bool2, eta_bool9, epsilon_bool5, epsilon_bool7, gamma_bool13, kappa_22, gamma_bool2, nu_bool3, epsilon_bool6, delta_bool8, nu_bool2, alpha_22, zeta_bool3, sigma_bool6, alpha_bool2, beta_bool5, delta_bool6, epsilon_bool1, eta_bool5, rho_22, eta_bool8, kappa_bool4, rho_bool6, delta_4, delta_bool9, xi_bool7, zeta_50, epsilon_50, xi_bool9, xi_0, delta_bool1, xi_bool4, gamma_bool5, beta_4, initalTimeManual, eta_bool3, zeta_bool2, kappa_bool10, epsilon_12, mu_bool6, gamma_bool14, gamma_bool3, delta_50, eta_38, eta_bool4, mu_bool3, epsilon_bool4, gamma_50, gamma_bool7, alpha_bool10, kappa_bool9, kappa_bool1, sigma_22, rho_bool8, eta_50, delta_22, alpha_bool9, zeta_0, sigma_bool1, beta_bool10, rho_50, alpha_bool7, nu_50, gamma_4, delta_bool4, gamma_bool12, sigma_bool10, sigma_bool4, alpha_bool5, alpha_bool13 = paramVec 

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

	u0Vec .= CumulativeDiagnosed, Infected, Extinct, Diagnosed, DiagnosedHealed, Ailing, Healed, Recognized, Susceptible, Threatened
end

function evalU0(paramVec) 

	beta_bool3, alpha_bool12, rho_bool2, rho_bool3, kappa_bool7, rho_38, nu_bool4, kappa_0, alpha_bool6, zeta_bool7, rho_bool7, zeta_bool5, eta_bool2, rho_bool10, lam_bool6, xi_bool10, xi_38, delta_bool2, sigma_38, zeta_bool6, kappa_bool3, gamma_bool8, kappa_bool6, beta_bool8, zeta_bool9, kappa_bool8, delta_bool7, epsilon_bool8, nu_bool5, tau, alpha_bool8, zeta_bool4, beta_bool9, theta, xi_bool8, zeta_22, sigma_bool7, alpha_bool1, mu_bool4, nu_bool1, eta_0, zeta_38, xi_bool2, rho_bool5, alpha_50, beta_0, kappa_bool2, alpha_bool3, zeta_bool1, eta_bool10, beta_bool6, delta_0, beta_bool1, xi_bool6, gamma_22, gamma_bool9, lam_bool3, xi_bool1, sigma_0, beta_22, zeta_bool8, lam_bool4, gamma_28, alpha_0, mu_22, eta_bool6, gamma_bool4, gamma_bool11, epsilon_bool9, delta_bool5, gamma_bool1, kappa_bool5, zeta_bool10, nu_bool6, lam_bool5, eta_bool1, mu_50, lam_50, kappa_38, alpha_bool14, alpha_28, gamma_bool10, rho_0, sigma_bool8, gamma_0, epsilon_bool10, gamma_bool6, alpha_4, sigma_50, lam_0, beta_50, epsilon_bool2, rho_bool4, mu_bool5, sigma_bool9, delta_bool10, alpha_bool4, nu_0, rho_bool1, mu_bool2, delta_bool3, xi_22, xi_bool5, sigma_bool2, beta_bool2, nu_22, beta_bool4, eta_22, xi_50, lam_22, eta_bool7, beta_bool7, mu_bool1, sigma_bool3, epsilon_bool3, epsilon_38, lam_bool1, alpha_bool11, epsilon_0, sigma_bool5, Italy, rho_bool9, kappa_50, xi_bool3, mu_0, lam_bool2, eta_bool9, epsilon_bool5, epsilon_bool7, gamma_bool13, kappa_22, gamma_bool2, nu_bool3, epsilon_bool6, delta_bool8, nu_bool2, alpha_22, zeta_bool3, sigma_bool6, alpha_bool2, beta_bool5, delta_bool6, epsilon_bool1, eta_bool5, rho_22, eta_bool8, kappa_bool4, rho_bool6, delta_4, delta_bool9, xi_bool7, zeta_50, epsilon_50, xi_bool9, xi_0, delta_bool1, xi_bool4, gamma_bool5, beta_4, initalTimeManual, eta_bool3, zeta_bool2, kappa_bool10, epsilon_12, mu_bool6, gamma_bool14, gamma_bool3, delta_50, eta_38, eta_bool4, mu_bool3, epsilon_bool4, gamma_50, gamma_bool7, alpha_bool10, kappa_bool9, kappa_bool1, sigma_22, rho_bool8, eta_50, delta_22, alpha_bool9, zeta_0, sigma_bool1, beta_bool10, rho_50, alpha_bool7, nu_50, gamma_4, delta_bool4, gamma_bool12, sigma_bool10, sigma_bool4, alpha_bool5, alpha_bool13 = paramVec 

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

function evalSd!(u, t, sdPar, dynPar, nonDynParam, paramData, observableId, mapSdParam) 

	CumulativeDiagnosed, Infected, Extinct, Diagnosed, DiagnosedHealed, Ailing, Healed, Recognized, Susceptible, Threatened= u 
	alpha_0, beta_0, delta_0, epsilon_0, eta_0, gamma_0, kappa_0, lam_0, mu_0, nu_0, rho_0, sigma_0, tau, theta, xi_0, zeta_0, alpha_4, alpha_22, alpha_28, beta_4, beta_22, delta_4, delta_22, epsilon_12, epsilon_38, eta_22, eta_38, gamma_4, gamma_22, gamma_28, kappa_22, kappa_38, lam_22, mu_22, nu_22, rho_22, rho_38, sigma_22, sigma_38, xi_22, xi_38, zeta_22, zeta_38 = dynPar 

	if observableId == "observable_TotalCases" 
		noiseParameter1_observable_TotalCases = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_TotalCases 
	end

	if observableId == "observable_Deaths" 
		noiseParameter1_observable_Deaths = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_Deaths 
	end

	if observableId == "observable_Recovered" 
		noiseParameter1_observable_Recovered = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_Recovered 
	end

	if observableId == "observable_CurrentCases" 
		noiseParameter1_observable_CurrentCases = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_CurrentCases 
	end

	if observableId == "observable_DiagHome" 
		noiseParameter1_observable_DiagHome = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_DiagHome 
	end

	if observableId == "observable_Hospit" 
		noiseParameter1_observable_Hospit = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_Hospit 
	end

	if observableId == "observable_ICU" 
		noiseParameter1_observable_ICU = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_observable_ICU 
	end

end