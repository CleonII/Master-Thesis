# Model name: model_Boehm_JProteomeRes2014
# Number of parameters: 9
# Number of species: 8
function getODEModel_model_Boehm_JProteomeRes2014()

    ### Define constant parameters
    cyt = 1.4
    nuc = 0.45

    ### Define independent and dependent variables
    ModelingToolkit.@variables t STAT5A(t) pApA(t) nucpApB(t) nucpBpB(t) STAT5B(t) pApB(t) nucpApA(t) pBpB(t)

    ### Define variable parameters

    ### Define dummy variable
    ModelingToolkit.@variables dummyVariable(t)

    ### Define parameters
    ModelingToolkit.@parameters ratio k_imp_homo k_exp_hetero k_phos specC17 Epo_degradation_BaF3 k_exp_homo k_imp_hetero

    ### Define an operator for the differentiation w.r.t. time
    D = Differential(t)

    ### Events ###

    ### Derivatives ###
    eqs = [
    D(STAT5A) ~ -2.0 * ( 1 /cyt ) * (cyt * (1.25e-7 * exp(-1 * Epo_degradation_BaF3 * t)) * (STAT5A)^(2) * k_phos)-1.0 * ( 1 /cyt ) * (cyt * (1.25e-7 * exp(-1 * Epo_degradation_BaF3 * t)) * STAT5A * STAT5B * k_phos)+2.0 * ( 1 /cyt ) * (nuc * k_exp_homo * nucpApA)+1.0 * ( 1 /cyt ) * (nuc * k_exp_hetero * nucpApB),
    D(pApA) ~ +1.0 * ( 1 /cyt ) * (cyt * (1.25e-7 * exp(-1 * Epo_degradation_BaF3 * t)) * (STAT5A)^(2) * k_phos)-1.0 * ( 1 /cyt ) * (cyt * k_imp_homo * pApA),
    D(nucpApB) ~ +1.0 * ( 1 /nuc ) * (cyt * k_imp_hetero * pApB)-1.0 * ( 1 /nuc ) * (nuc * k_exp_hetero * nucpApB),
    D(nucpBpB) ~ +1.0 * ( 1 /nuc ) * (cyt * k_imp_homo * pBpB)-1.0 * ( 1 /nuc ) * (nuc * k_exp_homo * nucpBpB),
    D(STAT5B) ~ -1.0 * ( 1 /cyt ) * (cyt * (1.25e-7 * exp(-1 * Epo_degradation_BaF3 * t)) * STAT5A * STAT5B * k_phos)-2.0 * ( 1 /cyt ) * (cyt * (1.25e-7 * exp(-1 * Epo_degradation_BaF3 * t)) * (STAT5B)^(2) * k_phos)+1.0 * ( 1 /cyt ) * (nuc * k_exp_hetero * nucpApB)+2.0 * ( 1 /cyt ) * (nuc * k_exp_homo * nucpBpB),
    D(pApB) ~ +1.0 * ( 1 /cyt ) * (cyt * (1.25e-7 * exp(-1 * Epo_degradation_BaF3 * t)) * STAT5A * STAT5B * k_phos)-1.0 * ( 1 /cyt ) * (cyt * k_imp_hetero * pApB),
    D(nucpApA) ~ +1.0 * ( 1 /nuc ) * (cyt * k_imp_homo * pApA)-1.0 * ( 1 /nuc ) * (nuc * k_exp_homo * nucpApA),
    D(pBpB) ~ +1.0 * ( 1 /cyt ) * (cyt * (1.25e-7 * exp(-1 * Epo_degradation_BaF3 * t)) * (STAT5B)^(2) * k_phos)-1.0 * ( 1 /cyt ) * (cyt * k_imp_homo * pBpB),
    D(dummyVariable) ~ +specC17+ratio
    ]

    @named sys = ODESystem(eqs)

    ### Initial species concentrations ###
    initialSpeciesValues = [
    STAT5A => 207.6 * ratio,
    pApA => 0.0,
    nucpApB => 0.0,
    nucpBpB => 0.0,
    STAT5B => 207.6 - 207.6 * ratio,
    pApB => 0.0,
    nucpApA => 0.0,
    pBpB => 0.0,
    dummyVariable => 0.0]

    ### True parameter values ###
    trueParameterValues = [
    ratio => 0.693,
    k_imp_homo => 96945.5391768823,
    k_exp_hetero => 1.00094251286741e-5,
    k_phos => 15766.8336642826,
    specC17 => 0.107,
    Epo_degradation_BaF3 => 0.0269738286367359,
    k_exp_homo => 0.00617193081581346,
    k_imp_hetero => 0.0163708512310568]

    return sys, initialSpeciesValues, trueParameterValues

end
