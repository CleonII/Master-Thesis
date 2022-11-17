function evalYmod(u, t, dynPar, obsPar, nonDynParam, paramData, observableId, mapObsParam) 

	CumulativeDiagnosed, Infected, Extinct, Diagnosed, DiagnosedHealed, Ailing, Healed, Recognized, Susceptible, Threatened= u 
	alpha_0, beta_0, delta_0, epsilon_0, eta_0, gamma_0, kappa_0, lam_0, mu_0, nu_0, rho_0, sigma_0, tau, theta, xi_0, zeta_0, alpha_4, alpha_22, alpha_28, beta_4, beta_22, delta_4, delta_22, epsilon_12, epsilon_38, eta_22, eta_38, gamma_4, gamma_22, gamma_28, kappa_22, kappa_38, lam_22, mu_22, nu_22, rho_22, rho_38, sigma_22, sigma_38, xi_22, xi_38, zeta_22, zeta_38 = dynPar 

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
		return  
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

	rho_22, nu_0, delta_4, beta_22, zeta_50, rho_38, gamma_28, xi_22, kappa_0, epsilon_50, nu_22, alpha_0, xi_0, mu_22, eta_22, xi_50, beta_4, lam_22, initalTimeManual, xi_38, epsilon_12, epsilon_38, sigma_38, epsilon_0, delta_50, Italy, mu_50, eta_38, lam_50, kappa_50, kappa_38, gamma_50, tau, alpha_28, mu_0, theta, zeta_22, kappa_22, rho_0, gamma_0, sigma_22, eta_50, delta_22, eta_0, zeta_0, zeta_38, alpha_4, rho_50, alpha_50, sigma_50, nu_50, gamma_4, beta_0, lam_0, beta_50, alpha_22, delta_0, gamma_22, sigma_0 = paramVec 

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