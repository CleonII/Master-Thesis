function evalDYmodDu(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	CumulativeDiagnosed, Infected, Extinct, Diagnosed, DiagnosedHealed, Ailing, Healed, Recognized, Susceptible, Threatened= u 
	beta_bool3, alpha_bool12, rho_bool2, rho_bool3, kappa_bool7, rho_38, nu_bool4, kappa_0, alpha_bool6, zeta_bool7, rho_bool7, zeta_bool5, eta_bool2, rho_bool10, lam_bool6, xi_bool10, xi_38, delta_bool2, sigma_38, zeta_bool6, kappa_bool3, gamma_bool8, kappa_bool6, beta_bool8, zeta_bool9, kappa_bool8, delta_bool7, epsilon_bool8, nu_bool5, tau, alpha_bool8, zeta_bool4, beta_bool9, theta, xi_bool8, zeta_22, sigma_bool7, alpha_bool1, mu_bool4, nu_bool1, eta_0, zeta_38, xi_bool2, rho_bool5, alpha_50, beta_0, kappa_bool2, alpha_bool3, zeta_bool1, eta_bool10, beta_bool6, delta_0, beta_bool1, xi_bool6, gamma_22, gamma_bool9, lam_bool3, xi_bool1, sigma_0, beta_22, zeta_bool8, lam_bool4, gamma_28, alpha_0, mu_22, eta_bool6, gamma_bool4, gamma_bool11, epsilon_bool9, delta_bool5, gamma_bool1, kappa_bool5, zeta_bool10, nu_bool6, lam_bool5, eta_bool1, mu_50, lam_50, kappa_38, alpha_bool14, alpha_28, gamma_bool10, rho_0, sigma_bool8, gamma_0, epsilon_bool10, gamma_bool6, alpha_4, sigma_50, lam_0, beta_50, epsilon_bool2, rho_bool4, mu_bool5, sigma_bool9, delta_bool10, alpha_bool4, nu_0, rho_bool1, mu_bool2, delta_bool3, xi_22, xi_bool5, sigma_bool2, beta_bool2, nu_22, beta_bool4, eta_22, xi_50, lam_22, eta_bool7, beta_bool7, mu_bool1, sigma_bool3, epsilon_bool3, epsilon_38, lam_bool1, alpha_bool11, epsilon_0, sigma_bool5, Italy, rho_bool9, kappa_50, xi_bool3, mu_0, lam_bool2, eta_bool9, epsilon_bool5, epsilon_bool7, gamma_bool13, kappa_22, gamma_bool2, nu_bool3, epsilon_bool6, delta_bool8, nu_bool2, alpha_22, zeta_bool3, sigma_bool6, alpha_bool2, beta_bool5, delta_bool6, epsilon_bool1, eta_bool5, rho_22, eta_bool8, kappa_bool4, rho_bool6, delta_4, delta_bool9, xi_bool7, zeta_50, epsilon_50, xi_bool9, xi_0, delta_bool1, xi_bool4, gamma_bool5, beta_4, initalTimeManual, eta_bool3, zeta_bool2, kappa_bool10, epsilon_12, mu_bool6, gamma_bool14, gamma_bool3, delta_50, eta_38, eta_bool4, mu_bool3, epsilon_bool4, gamma_50, gamma_bool7, alpha_bool10, kappa_bool9, kappa_bool1, sigma_22, rho_bool8, eta_50, delta_22, alpha_bool9, zeta_0, sigma_bool1, beta_bool10, rho_50, alpha_bool7, nu_50, gamma_4, delta_bool4, gamma_bool12, sigma_bool10, sigma_bool4, alpha_bool5, alpha_bool13 = p 
	CurrentDiagnosedInfected = Diagnosed + Recognized + Threatened

	if observableId == "observable_TotalCases" 
		out[1] = 1
		return nothing
	end

	if observableId == "observable_Deaths" 
		out[3] = 1
		return nothing
	end

	if observableId == "observable_Recovered" 
		out[5] = 1
		return nothing
	end

	if observableId == "observable_CurrentCases" 
		return nothing
	end

	if observableId == "observable_DiagHome" 
		out[4] = 1
		return nothing
	end

	if observableId == "observable_Hospit" 
		out[8] = 1
		return nothing
	end

	if observableId == "observable_ICU" 
		out[10] = 1
		return nothing
	end

end

function evalDYmodDp(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	CumulativeDiagnosed, Infected, Extinct, Diagnosed, DiagnosedHealed, Ailing, Healed, Recognized, Susceptible, Threatened= u 
	beta_bool3, alpha_bool12, rho_bool2, rho_bool3, kappa_bool7, rho_38, nu_bool4, kappa_0, alpha_bool6, zeta_bool7, rho_bool7, zeta_bool5, eta_bool2, rho_bool10, lam_bool6, xi_bool10, xi_38, delta_bool2, sigma_38, zeta_bool6, kappa_bool3, gamma_bool8, kappa_bool6, beta_bool8, zeta_bool9, kappa_bool8, delta_bool7, epsilon_bool8, nu_bool5, tau, alpha_bool8, zeta_bool4, beta_bool9, theta, xi_bool8, zeta_22, sigma_bool7, alpha_bool1, mu_bool4, nu_bool1, eta_0, zeta_38, xi_bool2, rho_bool5, alpha_50, beta_0, kappa_bool2, alpha_bool3, zeta_bool1, eta_bool10, beta_bool6, delta_0, beta_bool1, xi_bool6, gamma_22, gamma_bool9, lam_bool3, xi_bool1, sigma_0, beta_22, zeta_bool8, lam_bool4, gamma_28, alpha_0, mu_22, eta_bool6, gamma_bool4, gamma_bool11, epsilon_bool9, delta_bool5, gamma_bool1, kappa_bool5, zeta_bool10, nu_bool6, lam_bool5, eta_bool1, mu_50, lam_50, kappa_38, alpha_bool14, alpha_28, gamma_bool10, rho_0, sigma_bool8, gamma_0, epsilon_bool10, gamma_bool6, alpha_4, sigma_50, lam_0, beta_50, epsilon_bool2, rho_bool4, mu_bool5, sigma_bool9, delta_bool10, alpha_bool4, nu_0, rho_bool1, mu_bool2, delta_bool3, xi_22, xi_bool5, sigma_bool2, beta_bool2, nu_22, beta_bool4, eta_22, xi_50, lam_22, eta_bool7, beta_bool7, mu_bool1, sigma_bool3, epsilon_bool3, epsilon_38, lam_bool1, alpha_bool11, epsilon_0, sigma_bool5, Italy, rho_bool9, kappa_50, xi_bool3, mu_0, lam_bool2, eta_bool9, epsilon_bool5, epsilon_bool7, gamma_bool13, kappa_22, gamma_bool2, nu_bool3, epsilon_bool6, delta_bool8, nu_bool2, alpha_22, zeta_bool3, sigma_bool6, alpha_bool2, beta_bool5, delta_bool6, epsilon_bool1, eta_bool5, rho_22, eta_bool8, kappa_bool4, rho_bool6, delta_4, delta_bool9, xi_bool7, zeta_50, epsilon_50, xi_bool9, xi_0, delta_bool1, xi_bool4, gamma_bool5, beta_4, initalTimeManual, eta_bool3, zeta_bool2, kappa_bool10, epsilon_12, mu_bool6, gamma_bool14, gamma_bool3, delta_50, eta_38, eta_bool4, mu_bool3, epsilon_bool4, gamma_50, gamma_bool7, alpha_bool10, kappa_bool9, kappa_bool1, sigma_22, rho_bool8, eta_50, delta_22, alpha_bool9, zeta_0, sigma_bool1, beta_bool10, rho_50, alpha_bool7, nu_50, gamma_4, delta_bool4, gamma_bool12, sigma_bool10, sigma_bool4, alpha_bool5, alpha_bool13 = p 
	CurrentDiagnosedInfected = Diagnosed + Recognized + Threatened

	if observableId == "observable_TotalCases" 
		return nothing
	end

	if observableId == "observable_Deaths" 
		return nothing
	end

	if observableId == "observable_Recovered" 
		return nothing
	end

	if observableId == "observable_CurrentCases" 
		return nothing
	end

	if observableId == "observable_DiagHome" 
		return nothing
	end

	if observableId == "observable_Hospit" 
		return nothing
	end

	if observableId == "observable_ICU" 
		return nothing
	end

end

function evalDSdDu!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	CumulativeDiagnosed, Infected, Extinct, Diagnosed, DiagnosedHealed, Ailing, Healed, Recognized, Susceptible, Threatened= u 
	beta_bool3, alpha_bool12, rho_bool2, rho_bool3, kappa_bool7, rho_38, nu_bool4, kappa_0, alpha_bool6, zeta_bool7, rho_bool7, zeta_bool5, eta_bool2, rho_bool10, lam_bool6, xi_bool10, xi_38, delta_bool2, sigma_38, zeta_bool6, kappa_bool3, gamma_bool8, kappa_bool6, beta_bool8, zeta_bool9, kappa_bool8, delta_bool7, epsilon_bool8, nu_bool5, tau, alpha_bool8, zeta_bool4, beta_bool9, theta, xi_bool8, zeta_22, sigma_bool7, alpha_bool1, mu_bool4, nu_bool1, eta_0, zeta_38, xi_bool2, rho_bool5, alpha_50, beta_0, kappa_bool2, alpha_bool3, zeta_bool1, eta_bool10, beta_bool6, delta_0, beta_bool1, xi_bool6, gamma_22, gamma_bool9, lam_bool3, xi_bool1, sigma_0, beta_22, zeta_bool8, lam_bool4, gamma_28, alpha_0, mu_22, eta_bool6, gamma_bool4, gamma_bool11, epsilon_bool9, delta_bool5, gamma_bool1, kappa_bool5, zeta_bool10, nu_bool6, lam_bool5, eta_bool1, mu_50, lam_50, kappa_38, alpha_bool14, alpha_28, gamma_bool10, rho_0, sigma_bool8, gamma_0, epsilon_bool10, gamma_bool6, alpha_4, sigma_50, lam_0, beta_50, epsilon_bool2, rho_bool4, mu_bool5, sigma_bool9, delta_bool10, alpha_bool4, nu_0, rho_bool1, mu_bool2, delta_bool3, xi_22, xi_bool5, sigma_bool2, beta_bool2, nu_22, beta_bool4, eta_22, xi_50, lam_22, eta_bool7, beta_bool7, mu_bool1, sigma_bool3, epsilon_bool3, epsilon_38, lam_bool1, alpha_bool11, epsilon_0, sigma_bool5, Italy, rho_bool9, kappa_50, xi_bool3, mu_0, lam_bool2, eta_bool9, epsilon_bool5, epsilon_bool7, gamma_bool13, kappa_22, gamma_bool2, nu_bool3, epsilon_bool6, delta_bool8, nu_bool2, alpha_22, zeta_bool3, sigma_bool6, alpha_bool2, beta_bool5, delta_bool6, epsilon_bool1, eta_bool5, rho_22, eta_bool8, kappa_bool4, rho_bool6, delta_4, delta_bool9, xi_bool7, zeta_50, epsilon_50, xi_bool9, xi_0, delta_bool1, xi_bool4, gamma_bool5, beta_4, initalTimeManual, eta_bool3, zeta_bool2, kappa_bool10, epsilon_12, mu_bool6, gamma_bool14, gamma_bool3, delta_50, eta_38, eta_bool4, mu_bool3, epsilon_bool4, gamma_50, gamma_bool7, alpha_bool10, kappa_bool9, kappa_bool1, sigma_22, rho_bool8, eta_50, delta_22, alpha_bool9, zeta_0, sigma_bool1, beta_bool10, rho_50, alpha_bool7, nu_50, gamma_4, delta_bool4, gamma_bool12, sigma_bool10, sigma_bool4, alpha_bool5, alpha_bool13 = p 
	if observableId == "observable_TotalCases" 
		return nothing
	end

	if observableId == "observable_Deaths" 
		return nothing
	end

	if observableId == "observable_Recovered" 
		return nothing
	end

	if observableId == "observable_CurrentCases" 
		return nothing
	end

	if observableId == "observable_DiagHome" 
		return nothing
	end

	if observableId == "observable_Hospit" 
		return nothing
	end

	if observableId == "observable_ICU" 
		return nothing
	end

end

function evalDSdDp!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	CumulativeDiagnosed, Infected, Extinct, Diagnosed, DiagnosedHealed, Ailing, Healed, Recognized, Susceptible, Threatened= u 
	beta_bool3, alpha_bool12, rho_bool2, rho_bool3, kappa_bool7, rho_38, nu_bool4, kappa_0, alpha_bool6, zeta_bool7, rho_bool7, zeta_bool5, eta_bool2, rho_bool10, lam_bool6, xi_bool10, xi_38, delta_bool2, sigma_38, zeta_bool6, kappa_bool3, gamma_bool8, kappa_bool6, beta_bool8, zeta_bool9, kappa_bool8, delta_bool7, epsilon_bool8, nu_bool5, tau, alpha_bool8, zeta_bool4, beta_bool9, theta, xi_bool8, zeta_22, sigma_bool7, alpha_bool1, mu_bool4, nu_bool1, eta_0, zeta_38, xi_bool2, rho_bool5, alpha_50, beta_0, kappa_bool2, alpha_bool3, zeta_bool1, eta_bool10, beta_bool6, delta_0, beta_bool1, xi_bool6, gamma_22, gamma_bool9, lam_bool3, xi_bool1, sigma_0, beta_22, zeta_bool8, lam_bool4, gamma_28, alpha_0, mu_22, eta_bool6, gamma_bool4, gamma_bool11, epsilon_bool9, delta_bool5, gamma_bool1, kappa_bool5, zeta_bool10, nu_bool6, lam_bool5, eta_bool1, mu_50, lam_50, kappa_38, alpha_bool14, alpha_28, gamma_bool10, rho_0, sigma_bool8, gamma_0, epsilon_bool10, gamma_bool6, alpha_4, sigma_50, lam_0, beta_50, epsilon_bool2, rho_bool4, mu_bool5, sigma_bool9, delta_bool10, alpha_bool4, nu_0, rho_bool1, mu_bool2, delta_bool3, xi_22, xi_bool5, sigma_bool2, beta_bool2, nu_22, beta_bool4, eta_22, xi_50, lam_22, eta_bool7, beta_bool7, mu_bool1, sigma_bool3, epsilon_bool3, epsilon_38, lam_bool1, alpha_bool11, epsilon_0, sigma_bool5, Italy, rho_bool9, kappa_50, xi_bool3, mu_0, lam_bool2, eta_bool9, epsilon_bool5, epsilon_bool7, gamma_bool13, kappa_22, gamma_bool2, nu_bool3, epsilon_bool6, delta_bool8, nu_bool2, alpha_22, zeta_bool3, sigma_bool6, alpha_bool2, beta_bool5, delta_bool6, epsilon_bool1, eta_bool5, rho_22, eta_bool8, kappa_bool4, rho_bool6, delta_4, delta_bool9, xi_bool7, zeta_50, epsilon_50, xi_bool9, xi_0, delta_bool1, xi_bool4, gamma_bool5, beta_4, initalTimeManual, eta_bool3, zeta_bool2, kappa_bool10, epsilon_12, mu_bool6, gamma_bool14, gamma_bool3, delta_50, eta_38, eta_bool4, mu_bool3, epsilon_bool4, gamma_50, gamma_bool7, alpha_bool10, kappa_bool9, kappa_bool1, sigma_22, rho_bool8, eta_50, delta_22, alpha_bool9, zeta_0, sigma_bool1, beta_bool10, rho_50, alpha_bool7, nu_50, gamma_4, delta_bool4, gamma_bool12, sigma_bool10, sigma_bool4, alpha_bool5, alpha_bool13 = p 
	if observableId == "observable_TotalCases" 
		return nothing
	end

	if observableId == "observable_Deaths" 
		return nothing
	end

	if observableId == "observable_Recovered" 
		return nothing
	end

	if observableId == "observable_CurrentCases" 
		return nothing
	end

	if observableId == "observable_DiagHome" 
		return nothing
	end

	if observableId == "observable_Hospit" 
		return nothing
	end

	if observableId == "observable_ICU" 
		return nothing
	end

end

