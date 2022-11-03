# Model name: model_Giordano_Nature2020
# Number of parameters: 72
# Number of species: 13
function getODEModel_model_Giordano_Nature2020()

    ### Define independent and dependent variables
    ModelingToolkit.@variables t CumulativeDiagnosed(t) Infected(t) Extinct(t) Diagnosed(t) DiagnosedHealed(t) Ailing(t) Healed(t) Recognized(t) Susceptible(t) Threatened(t)

    ### Define variable parameters

    ### Define potential algebraic variables
    ModelingToolkit.@variables nu(t) lam(t) mu(t) delta(t) sigma(t) kappa(t) xi(t) alpha(t) zeta(t) eta(t) rho(t) beta(t) epsilon(t) gamma(t)

    ### Define dummy variable
    ModelingToolkit.@variables dummyVariable(t)

    ### Define parameters
    ModelingToolkit.@parameters rho_22 nu_0 delta_4 beta_22 zeta_50 rho_38 gamma_28 xi_22 kappa_0 epsilon_50 nu_22 alpha_0 xi_0 mu_22 eta_22 xi_50 beta_4 lam_22 initalTimeManual xi_38 epsilon_12 epsilon_38 sigma_38 epsilon_0 delta_50 Italy mu_50 eta_38 lam_50 kappa_50 kappa_38 gamma_50 tau alpha_28 mu_0 theta zeta_22 kappa_22 rho_0 gamma_0 sigma_22 eta_50 delta_22 eta_0 zeta_0 zeta_38 alpha_4 rho_50 alpha_50 sigma_50 nu_50 gamma_4 beta_0 lam_0 beta_50 alpha_22 delta_0 gamma_22 sigma_0

    ### Define an operator for the differentiation w.r.t. time
    D = Differential(t)

    ### Continious events ###

    ### Discrete events ###

    ### Derivatives ###
    eqs = [
    D(CumulativeDiagnosed) ~ +1.0 * ( 1 /Italy ) * (Italy * epsilon * Infected)+1.0 * ( 1 /Italy ) * (Italy * theta * Ailing)+1.0 * ( 1 /Italy ) * (Italy * mu * Ailing),
    D(Infected) ~ +1.0 * ( 1 /Italy ) * (Italy * (Susceptible * (alpha * Infected + beta * Diagnosed + gamma * Ailing + delta * Recognized)))-1.0 * ( 1 /Italy ) * (Italy * epsilon * Infected)-1.0 * ( 1 /Italy ) * (Italy * zeta * Infected)-1.0 * ( 1 /Italy ) * (Italy * lam * Infected),
    D(Extinct) ~ +1.0 * ( 1 /Italy ) * (Italy * tau * Threatened),
    D(Diagnosed) ~ +1.0 * ( 1 /Italy ) * (Italy * epsilon * Infected)-1.0 * ( 1 /Italy ) * (Italy * eta * Diagnosed)-1.0 * ( 1 /Italy ) * (Italy * rho * Diagnosed),
    D(DiagnosedHealed) ~ +1.0 * ( 1 /Italy ) * (Italy * rho * Diagnosed)+1.0 * ( 1 /Italy ) * (Italy * xi * Recognized)+1.0 * ( 1 /Italy ) * (Italy * sigma * Threatened),
    D(Ailing) ~ +1.0 * ( 1 /Italy ) * (Italy * zeta * Infected)-1.0 * ( 1 /Italy ) * (Italy * theta * Ailing)-1.0 * ( 1 /Italy ) * (Italy * kappa * Ailing)-1.0 * ( 1 /Italy ) * (Italy * mu * Ailing),
    D(Healed) ~ +1.0 * ( 1 /Italy ) * (Italy * lam * Infected)+1.0 * ( 1 /Italy ) * (Italy * rho * Diagnosed)+1.0 * ( 1 /Italy ) * (Italy * kappa * Ailing)+1.0 * ( 1 /Italy ) * (Italy * xi * Recognized)+1.0 * ( 1 /Italy ) * (Italy * sigma * Threatened),
    D(Recognized) ~ +1.0 * ( 1 /Italy ) * (Italy * eta * Diagnosed)+1.0 * ( 1 /Italy ) * (Italy * theta * Ailing)-1.0 * ( 1 /Italy ) * (Italy * nu * Recognized)-1.0 * ( 1 /Italy ) * (Italy * xi * Recognized),
    D(Susceptible) ~ -1.0 * ( 1 /Italy ) * (Italy * (Susceptible * (alpha * Infected + beta * Diagnosed + gamma * Ailing + delta * Recognized))),
    D(Threatened) ~ +1.0 * ( 1 /Italy ) * (Italy * mu * Ailing)+1.0 * ( 1 /Italy ) * (Italy * nu * Recognized)-1.0 * ( 1 /Italy ) * (Italy * tau * Threatened)-1.0 * ( 1 /Italy ) * (Italy * sigma * Threatened),
    nu ~ ifelse(t <= 22 + initalTimeManual, nu_0, 0) + (ifelse(t > 22 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 50, 1.0, 0.0)) * (nu_22) + (1 - (ifelse(t > 22 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 50, 1.0, 0.0))) * (0) + ifelse(t > 50 + initalTimeManual, nu_50, 0),
    lam ~ ifelse(t <= 22 + initalTimeManual, lam_0, 0) + (ifelse(t > 22 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 50, 1.0, 0.0)) * (lam_22) + (1 - (ifelse(t > 22 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 50, 1.0, 0.0))) * (0) + ifelse(t > 50 + initalTimeManual, lam_50, 0),
    mu ~ ifelse(t <= 22 + initalTimeManual, mu_0, 0) + (ifelse(t > 22 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 50, 1.0, 0.0)) * (mu_22) + (1 - (ifelse(t > 22 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 50, 1.0, 0.0))) * (0) + ifelse(t > 50 + initalTimeManual, mu_50, 0),
    delta ~ ifelse(t <= 4 + initalTimeManual, delta_0, 0) + (ifelse(t > 4 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 22, 1.0, 0.0)) * (delta_4) + (1 - (ifelse(t > 4 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 22, 1.0, 0.0))) * (0) + (ifelse(t > 22 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 50, 1.0, 0.0)) * (delta_22) + (1 - (ifelse(t > 22 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 50, 1.0, 0.0))) * (0) + ifelse(t > 50 + initalTimeManual, delta_50, 0),
    sigma ~ ifelse(t <= 22 + initalTimeManual, sigma_0, 0) + (ifelse(t > 22 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 38, 1.0, 0.0)) * (sigma_22) + (1 - (ifelse(t > 22 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 38, 1.0, 0.0))) * (0) + (ifelse(t > 38 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 50, 1.0, 0.0)) * (sigma_38) + (1 - (ifelse(t > 38 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 50, 1.0, 0.0))) * (0) + ifelse(t > 50 + initalTimeManual, sigma_50, 0),
    kappa ~ ifelse(t <= 22 + initalTimeManual, kappa_0, 0) + (ifelse(t > 22 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 38, 1.0, 0.0)) * (kappa_22) + (1 - (ifelse(t > 22 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 38, 1.0, 0.0))) * (0) + (ifelse(t > 38 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 50, 1.0, 0.0)) * (kappa_38) + (1 - (ifelse(t > 38 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 50, 1.0, 0.0))) * (0) + ifelse(t > 50 + initalTimeManual, kappa_50, 0),
    xi ~ ifelse(t <= 22 + initalTimeManual, xi_0, 0) + (ifelse(t > 22 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 38, 1.0, 0.0)) * (xi_22) + (1 - (ifelse(t > 22 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 38, 1.0, 0.0))) * (0) + (ifelse(t > 38 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 50, 1.0, 0.0)) * (xi_38) + (1 - (ifelse(t > 38 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 50, 1.0, 0.0))) * (0) + ifelse(t > 50 + initalTimeManual, xi_50, 0),
    alpha ~ ifelse(t <= 4 + initalTimeManual, alpha_0, 0) + (ifelse(t > 4 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 22, 1.0, 0.0)) * (alpha_4) + (1 - (ifelse(t > 4 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 22, 1.0, 0.0))) * (0) + (ifelse(t > 22 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 28, 1.0, 0.0)) * (alpha_22) + (1 - (ifelse(t > 22 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 28, 1.0, 0.0))) * (0) + (ifelse(t > 28 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 50, 1.0, 0.0)) * (alpha_28) + (1 - (ifelse(t > 28 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 50, 1.0, 0.0))) * (0) + ifelse(t > 50 + initalTimeManual, alpha_50, 0),
    zeta ~ ifelse(t <= 22 + initalTimeManual, zeta_0, 0) + (ifelse(t > 22 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 38, 1.0, 0.0)) * (zeta_22) + (1 - (ifelse(t > 22 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 38, 1.0, 0.0))) * (0) + (ifelse(t > 38 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 50, 1.0, 0.0)) * (zeta_38) + (1 - (ifelse(t > 38 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 50, 1.0, 0.0))) * (0) + ifelse(t > 50 + initalTimeManual, zeta_50, 0),
    eta ~ ifelse(t <= 22 + initalTimeManual, eta_0, 0) + (ifelse(t > 22 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 38, 1.0, 0.0)) * (eta_22) + (1 - (ifelse(t > 22 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 38, 1.0, 0.0))) * (0) + (ifelse(t > 38 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 50, 1.0, 0.0)) * (eta_38) + (1 - (ifelse(t > 38 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 50, 1.0, 0.0))) * (0) + ifelse(t > 50 + initalTimeManual, eta_50, 0),
    rho ~ ifelse(t <= 22 + initalTimeManual, rho_0, 0) + (ifelse(t > 22 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 38, 1.0, 0.0)) * (rho_22) + (1 - (ifelse(t > 22 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 38, 1.0, 0.0))) * (0) + (ifelse(t > 38 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 50, 1.0, 0.0)) * (rho_38) + (1 - (ifelse(t > 38 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 50, 1.0, 0.0))) * (0) + ifelse(t > 50 + initalTimeManual, rho_50, 0),
    beta ~ ifelse(t <= 4 + initalTimeManual, beta_0, 0) + (ifelse(t > 4 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 22, 1.0, 0.0)) * (beta_4) + (1 - (ifelse(t > 4 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 22, 1.0, 0.0))) * (0) + (ifelse(t > 22 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 50, 1.0, 0.0)) * (beta_22) + (1 - (ifelse(t > 22 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 50, 1.0, 0.0))) * (0) + ifelse(t > 50 + initalTimeManual, beta_50, 0),
    epsilon ~ ifelse(t <= 12 + initalTimeManual, epsilon_0, 0) + (ifelse(t > 12 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 38, 1.0, 0.0)) * (epsilon_12) + (1 - (ifelse(t > 12 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 38, 1.0, 0.0))) * (0) + (ifelse(t > 38 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 50, 1.0, 0.0)) * (epsilon_38) + (1 - (ifelse(t > 38 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 50, 1.0, 0.0))) * (0) + ifelse(t > 50 + initalTimeManual, epsilon_50, 0),
    gamma ~ ifelse(t <= 4 + initalTimeManual, gamma_0, 0) + (ifelse(t > 4 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 22, 1.0, 0.0)) * (gamma_4) + (1 - (ifelse(t > 4 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 22, 1.0, 0.0))) * (0) + (ifelse(t > 22 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 28, 1.0, 0.0)) * (gamma_22) + (1 - (ifelse(t > 22 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 28, 1.0, 0.0))) * (0) + (ifelse(t > 28 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 50, 1.0, 0.0)) * (gamma_28) + (1 - (ifelse(t > 28 + initalTimeManual, 1.0, 0.0)) * (ifelse(t <= initalTimeManual + 50, 1.0, 0.0))) * (0) + ifelse(t > 50 + initalTimeManual, gamma_50, 0),
    D(dummyVariable) ~ 1e-60*( +rho_22+nu_0+delta_4+beta_22+zeta_50+rho_38+gamma_28+epsilon_50+xi_22+kappa_0+nu_22+alpha_0+xi_0+mu_22+eta_22+xi_50+beta_4+lam_22+initalTimeManual+xi_38+epsilon_12+epsilon_38+sigma_38+epsilon_0+delta_50+mu_50+eta_38+lam_50+kappa_50+kappa_38+gamma_50+alpha_28+mu_0+zeta_22+kappa_22+rho_0+gamma_0+sigma_22+eta_50+delta_22+eta_0+zeta_0+zeta_38+alpha_4+rho_50+alpha_50+sigma_50+nu_50+gamma_4+beta_0+lam_0+beta_50+alpha_22+delta_0+gamma_22+sigma_0)
    ]

    @named sys = ODESystem(eqs)

    ### Initial species concentrations ###
    initialSpeciesValues = [
    CumulativeDiagnosed => 3.33333333e-8 + 3.33333333e-7,
    Infected => 3.33333333e-6,
    Extinct => 0.0,
    Diagnosed => 3.33333333e-7,
    DiagnosedHealed => 0.0,
    Ailing => 1.66666666e-8,
    Healed => 0.0,
    Recognized => 3.33333333e-8,
    Susceptible => 0.9999963,
    Threatened => 0.0,
    dummyVariable => 0.0]

    ### SBML file parameter values ###
    trueParameterValues = [
    rho_22 => 0.017,
    nu_0 => 0.027,
    delta_4 => 0.0057,
    beta_22 => 0.005,
    zeta_50 => 0.025,
    rho_38 => 0.02,
    gamma_28 => 0.11,
    xi_22 => 0.017,
    kappa_0 => 0.017,
    epsilon_50 => 0.2,
    nu_22 => 0.015,
    alpha_0 => 0.57,
    xi_0 => 0.017,
    mu_22 => 0.008,
    eta_22 => 0.034,
    xi_50 => 0.02,
    beta_4 => 0.0057,
    lam_22 => 0.08,
    initalTimeManual => 0.0,
    xi_38 => 0.02,
    epsilon_12 => 0.143,
    epsilon_38 => 0.2,
    sigma_38 => 0.01,
    epsilon_0 => 0.171,
    delta_50 => 0.005,
    Italy => 1.0,
    mu_50 => 0.008,
    eta_38 => 0.025,
    lam_50 => 0.08,
    kappa_50 => 0.02,
    kappa_38 => 0.02,
    gamma_50 => 0.11,
    tau => 0.01,
    alpha_28 => 0.21,
    mu_0 => 0.017,
    theta => 0.371,
    zeta_22 => 0.034,
    kappa_22 => 0.017,
    rho_0 => 0.034,
    gamma_0 => 0.456,
    sigma_22 => 0.017,
    eta_50 => 0.025,
    delta_22 => 0.005,
    eta_0 => 0.125,
    zeta_0 => 0.125,
    zeta_38 => 0.025,
    alpha_4 => 0.422,
    rho_50 => 0.02,
    alpha_50 => 0.21,
    sigma_50 => 0.01,
    nu_50 => 0.015,
    gamma_4 => 0.285,
    beta_0 => 0.011,
    lam_0 => 0.034,
    beta_50 => 0.005,
    alpha_22 => 0.36,
    delta_0 => 0.011,
    gamma_22 => 0.2,
    sigma_0 => 0.017]

    return sys, initialSpeciesValues, trueParameterValues

end
