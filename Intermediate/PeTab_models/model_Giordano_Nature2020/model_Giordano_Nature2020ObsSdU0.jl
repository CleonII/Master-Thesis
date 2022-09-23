function evalYmod(u, t, dynPar, obsPar, paramData, obsData, observableId, simulationId) 

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

	initalTimeManual, epsilon_50, epsilon_12, mu_22, epsilon_38, mu_50, mu_0, epsilon_0, theta, beta_50, lam_50, zeta_38, gamma_4, gamma_28, delta_0, gamma_0, alpha_22, lam_22, zeta_50, alpha_28, alpha_0, delta_50, alpha_4, beta_22, lam_0, zeta_22, alpha_50, beta_0, gamma_50, beta_4, zeta_0, delta_4, delta_22, gamma_22, tau, eta_22, rho_38, rho_0, rho_22, rho_50, eta_0, eta_50, eta_38, sigma_22, sigma_38, xi_38, xi_22, sigma_50, xi_50, sigma_0, xi_0, kappa_22, kappa_38, kappa_0, kappa_50, nu_22, nu_50, nu_0 = paramVec 

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

function evalSd!(u, t, sdPar, dynPar, paramData, obsData, observableId, simulationId) 

	CumulativeDiagnosed, Infected, Extinct, Diagnosed, DiagnosedHealed, Ailing, Healed, Recognized, Susceptible, Threatened= u 
	alpha_0, beta_0, delta_0, epsilon_0, eta_0, gamma_0, kappa_0, lam_0, mu_0, nu_0, rho_0, sigma_0, tau, theta, xi_0, zeta_0, alpha_4, alpha_22, alpha_28, beta_4, beta_22, delta_4, delta_22, epsilon_12, epsilon_38, eta_22, eta_38, gamma_4, gamma_22, gamma_28, kappa_22, kappa_38, lam_22, mu_22, nu_22, rho_22, rho_38, sigma_22, sigma_38, xi_22, xi_38, zeta_22, zeta_38 = dynPar 

	if observableId == "observable_TotalCases" 
		noiseParameter1_observable_TotalCases = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_TotalCases 
	end

	if observableId == "observable_Deaths" 
		noiseParameter1_observable_Deaths = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_Deaths 
	end

	if observableId == "observable_Recovered" 
		noiseParameter1_observable_Recovered = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_Recovered 
	end

	if observableId == "observable_CurrentCases" 
		noiseParameter1_observable_CurrentCases = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_CurrentCases 
	end

	if observableId == "observable_DiagHome" 
		noiseParameter1_observable_DiagHome = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_DiagHome 
	end

	if observableId == "observable_Hospit" 
		noiseParameter1_observable_Hospit = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_Hospit 
	end

	if observableId == "observable_ICU" 
		noiseParameter1_observable_ICU = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_ICU 
	end

end