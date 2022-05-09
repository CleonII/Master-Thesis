# Model name: model_Elowitz_Nature2000
# Number of parameters: 18
# Number of species: 8
function getODEModel_model_Elowitz_Nature2000()

    ### Define independent and dependent variables
    ModelingToolkit.@variables t X_protein(t) GFP_mRNA(t) Y_mRNA(t) X_mRNA(t) Z_mRNA(t) Z_protein(t) Y_protein(t) GFP(t)

    ### Define variable parameters

    ### Define dummy variable
    ModelingToolkit.@variables dummyVariable(t)

    ### Define parameters
    ModelingToolkit.@parameters tau_mRNA tps_repr init_GFP n_Hill init_Y_mRNA init_Z_mRNA init_X_protein init_Y_protein tau_mRNA_GFP eff init_Z_protein tau_prot_GFP eff_GFP tps_active init_GFP_mRNA KM init_X_mRNA tau_prot cell

    ### Define an operator for the differentiation w.r.t. time
    D = Differential(t)

    ### Events ###

    ### Derivatives ###
    eqs = [
    D(X_protein) ~ +1.0 * (cell * (X_mRNA * eff / tau_mRNA))-1.0 * (cell * (X_protein * log(2) / tau_prot)),
    D(GFP_mRNA) ~ -1.0 * (cell * (GFP_mRNA * log(2) / tau_mRNA_GFP)),
    D(Y_mRNA) ~ -1.0 * (cell * (Y_mRNA * log(2) / tau_mRNA))+1.0 * (cell * (60 * tps_repr - (KM)^(n_Hill) * (60 * tps_repr - 60 * tps_active) / ((KM)^(n_Hill) + (X_protein)^(n_Hill)))),
    D(X_mRNA) ~ -1.0 * (cell * (X_mRNA * log(2) / tau_mRNA))+1.0 * (cell * (60 * tps_repr - (KM)^(n_Hill) * (60 * tps_repr - 60 * tps_active) / ((KM)^(n_Hill) + (Z_protein)^(n_Hill)))),
    D(Z_mRNA) ~ -1.0 * (cell * (Z_mRNA * log(2) / tau_mRNA))+1.0 * (cell * (60 * tps_repr - (KM)^(n_Hill) * (60 * tps_repr - 60 * tps_active) / ((KM)^(n_Hill) + (Y_protein)^(n_Hill)))),
    D(Z_protein) ~ +1.0 * (cell * (Z_mRNA * eff / tau_mRNA))-1.0 * (cell * (Z_protein * log(2) / tau_prot)),
    D(Y_protein) ~ +1.0 * (cell * (Y_mRNA * eff / tau_mRNA))-1.0 * (cell * (Y_protein * log(2) / tau_prot)),
    D(GFP) ~ +1.0 * (cell * (GFP_mRNA * eff_GFP / tau_mRNA_GFP))-1.0 * (cell * (GFP * log(2) / tau_prot_GFP))+1.0 * (cell * (60 * tps_repr - (KM)^(n_Hill) * (60 * tps_repr - 60 * tps_active) / ((KM)^(n_Hill) + (X_protein)^(n_Hill)))),
    D(dummyVariable) ~ +init_X_protein+init_GFP+init_GFP_mRNA+init_Y_mRNA+init_Y_protein+init_Z_protein+init_X_mRNA+init_Z_mRNA
    ]

    @named sys = ODESystem(eqs)

    ### Initial species concentrations ###
    initialSpeciesValues = [
    X_protein => init_X_protein,
    GFP_mRNA => init_GFP_mRNA,
    Y_mRNA => init_Y_mRNA,
    X_mRNA => init_X_mRNA,
    Z_mRNA => init_Z_mRNA,
    Z_protein => init_Z_protein,
    Y_protein => init_Y_protein,
    GFP => init_GFP,
    dummyVariable => 0.0]

    ### True parameter values ###
    trueParameterValues = [
    tau_mRNA => 12.624824570772,
    tps_repr => 1.00000000000883e-5,
    init_GFP => 3.38716998236181e-5,
    n_Hill => 1.52254567553984,
    init_Y_mRNA => 999.999945990876,
    init_Z_mRNA => 19.3670294497271,
    init_X_protein => 30.8087735629587,
    init_Y_protein => 0.000818268062900602,
    tau_mRNA_GFP => 0.000404064377463447,
    eff => 0.0301592803622868,
    init_Z_protein => 994.381959318238,
    tau_prot_GFP => 71.9312640934983,
    eff_GFP => 0.00195212197096568,
    tps_active => 0.612278991244992,
    init_GFP_mRNA => 131.435737789559,
    KM => 1.00013184764194e-5,
    init_X_mRNA => 2.55665758135759,
    tau_prot => 5.35926527470063,
    cell => 1.0]

    return sys, initialSpeciesValues, trueParameterValues

end
