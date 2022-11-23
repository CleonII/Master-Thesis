# Model name: model_Bachmann_MSB2011
# Number of parameters: 39
# Number of species: 25
function getODEModel_model_Bachmann_MSB2011()

    ### Define independent and dependent variables
    ModelingToolkit.@variables t p1EpoRpJAK2(t) pSTAT5(t) EpoRJAK2_CIS(t) SOCS3nRNA4(t) SOCS3RNA(t) SHP1(t) STAT5(t) EpoRJAK2(t) CISnRNA1(t) SOCS3nRNA1(t) SOCS3nRNA2(t) CISnRNA3(t) CISnRNA4(t) SOCS3(t) CISnRNA5(t) SOCS3nRNA5(t) SOCS3nRNA3(t) SHP1Act(t) npSTAT5(t) p12EpoRpJAK2(t) p2EpoRpJAK2(t) CIS(t) EpoRpJAK2(t) CISnRNA2(t) CISRNA(t)

    ### Store dependent variables in array for ODESystem command
    stateArray = [p1EpoRpJAK2, pSTAT5, EpoRJAK2_CIS, SOCS3nRNA4, SOCS3RNA, SHP1, STAT5, EpoRJAK2, CISnRNA1, SOCS3nRNA1, SOCS3nRNA2, CISnRNA3, CISnRNA4, SOCS3, CISnRNA5, SOCS3nRNA5, SOCS3nRNA3, SHP1Act, npSTAT5, p12EpoRpJAK2, p2EpoRpJAK2, CIS, EpoRpJAK2, CISnRNA2, CISRNA]

    ### Define variable parameters

    ### Define potential algebraic variables

    ### Define parameters
    ModelingToolkit.@parameters STAT5Exp STAT5Imp init_SOCS3_multiplier EpoRCISRemove STAT5ActEpoR SHP1ActEpoR JAK2EpoRDeaSHP1 CISTurn SOCS3Turn init_EpoRJAK2_CIS SOCS3Inh ActD init_CIS_multiplier CISRNAEqc JAK2ActEpo Epo SOCS3oe CISInh SHP1Dea SOCS3EqcOE CISRNADelay init_SHP1 CISEqcOE EpoRActJAK2 SOCS3RNAEqc CISEqc SHP1ProOE SOCS3RNADelay init_STAT5 CISoe CISRNATurn init_SHP1_multiplier init_EpoRJAK2 EpoRCISInh STAT5ActJAK2 SOCS3RNATurn SOCS3Eqc cyt nuc

    ### Store parameters in array for ODESystem command
    parameterArray = [STAT5Exp, STAT5Imp, init_SOCS3_multiplier, EpoRCISRemove, STAT5ActEpoR, SHP1ActEpoR, JAK2EpoRDeaSHP1, CISTurn, SOCS3Turn, init_EpoRJAK2_CIS, SOCS3Inh, ActD, init_CIS_multiplier, CISRNAEqc, JAK2ActEpo, Epo, SOCS3oe, CISInh, SHP1Dea, SOCS3EqcOE, CISRNADelay, init_SHP1, CISEqcOE, EpoRActJAK2, SOCS3RNAEqc, CISEqc, SHP1ProOE, SOCS3RNADelay, init_STAT5, CISoe, CISRNATurn, init_SHP1_multiplier, init_EpoRJAK2, EpoRCISInh, STAT5ActJAK2, SOCS3RNATurn, SOCS3Eqc, cyt, nuc]

    ### Define an operator for the differentiation w.r.t. time
    D = Differential(t)

    ### Continious events ###

    ### Discrete events ###

    ### Derivatives ###
    eqs = [
    D(p1EpoRpJAK2) ~ (EpoRActJAK2*EpoRpJAK2) / (1 + (SOCS3Inh*SOCS3) / (SOCS3Eqc*cyt)) + (-3EpoRActJAK2*p1EpoRpJAK2) / ((1 + (EpoRCISInh*EpoRJAK2_CIS) / cyt)*(1 + (SOCS3Inh*SOCS3) / (SOCS3Eqc*cyt))) + (-JAK2EpoRDeaSHP1*SHP1Act*p1EpoRpJAK2) / (cyt*init_SHP1),
    D(pSTAT5) ~ (STAT5ActEpoR*((p12EpoRpJAK2 / cyt + p1EpoRpJAK2 / cyt)^2.0)*STAT5) / ((init_EpoRJAK2^2.0)*(1 + (CISInh*CIS) / (CISEqc*cyt))*(1 + (SOCS3Inh*SOCS3) / (SOCS3Eqc*cyt))) + (STAT5ActJAK2*(EpoRpJAK2 / cyt + p12EpoRpJAK2 / cyt + p1EpoRpJAK2 / cyt + p2EpoRpJAK2 / cyt)*STAT5) / (init_EpoRJAK2*(1 + (SOCS3Inh*SOCS3) / (SOCS3Eqc*cyt))) - STAT5Imp*pSTAT5,
    D(EpoRJAK2_CIS) ~ -((EpoRCISRemove*(p12EpoRpJAK2 / cyt + p1EpoRpJAK2 / cyt)*EpoRJAK2_CIS) / init_EpoRJAK2),
    D(SOCS3nRNA4) ~ SOCS3RNADelay*SOCS3nRNA3 - SOCS3RNADelay*SOCS3nRNA4,
    D(SOCS3RNA) ~ SOCS3RNADelay*SOCS3nRNA5 - SOCS3RNATurn*SOCS3RNA,
    D(SHP1) ~ SHP1Dea*SHP1Act + (-SHP1ActEpoR*(EpoRpJAK2 / cyt + p12EpoRpJAK2 / cyt + p1EpoRpJAK2 / cyt + p2EpoRpJAK2 / cyt)*SHP1) / init_EpoRJAK2,
    D(STAT5) ~ STAT5Exp*npSTAT5 + (-STAT5ActJAK2*(EpoRpJAK2 / cyt + p12EpoRpJAK2 / cyt + p1EpoRpJAK2 / cyt + p2EpoRpJAK2 / cyt)*STAT5) / (init_EpoRJAK2*(1 + (SOCS3Inh*SOCS3) / (SOCS3Eqc*cyt))) + (-STAT5ActEpoR*((p12EpoRpJAK2 / cyt + p1EpoRpJAK2 / cyt)^2.0)*STAT5) / ((init_EpoRJAK2^2.0)*(1 + (CISInh*CIS) / (CISEqc*cyt))*(1 + (SOCS3Inh*SOCS3) / (SOCS3Eqc*cyt))),
    D(EpoRJAK2) ~ (JAK2EpoRDeaSHP1*EpoRpJAK2*SHP1Act) / (cyt*init_SHP1) + (JAK2EpoRDeaSHP1*SHP1Act*p12EpoRpJAK2) / (cyt*init_SHP1) + (JAK2EpoRDeaSHP1*SHP1Act*p1EpoRpJAK2) / (cyt*init_SHP1) + (JAK2EpoRDeaSHP1*SHP1Act*p2EpoRpJAK2) / (cyt*init_SHP1) + (-Epo*JAK2ActEpo*EpoRJAK2) / (1 + (SOCS3Inh*SOCS3) / (SOCS3Eqc*cyt)),
    D(CISnRNA1) ~ (ActD*CISRNAEqc*CISRNATurn*npSTAT5) / init_STAT5 - CISRNADelay*CISnRNA1,
    D(SOCS3nRNA1) ~ (ActD*SOCS3RNAEqc*SOCS3RNATurn*npSTAT5) / init_STAT5 - SOCS3RNADelay*SOCS3nRNA1,
    D(SOCS3nRNA2) ~ SOCS3RNADelay*SOCS3nRNA1 - SOCS3RNADelay*SOCS3nRNA2,
    D(CISnRNA3) ~ CISRNADelay*CISnRNA2 - CISRNADelay*CISnRNA3,
    D(CISnRNA4) ~ CISRNADelay*CISnRNA3 - CISRNADelay*CISnRNA4,
    D(SOCS3) ~ (SOCS3Eqc*SOCS3Turn*SOCS3RNA) / SOCS3RNAEqc + SOCS3Eqc*SOCS3EqcOE*SOCS3Turn*SOCS3oe*cyt - SOCS3Turn*SOCS3,
    D(CISnRNA5) ~ CISRNADelay*CISnRNA4 - CISRNADelay*CISnRNA5,
    D(SOCS3nRNA5) ~ SOCS3RNADelay*SOCS3nRNA4 - SOCS3RNADelay*SOCS3nRNA5,
    D(SOCS3nRNA3) ~ SOCS3RNADelay*SOCS3nRNA2 - SOCS3RNADelay*SOCS3nRNA3,
    D(SHP1Act) ~ (SHP1ActEpoR*(EpoRpJAK2 / cyt + p12EpoRpJAK2 / cyt + p1EpoRpJAK2 / cyt + p2EpoRpJAK2 / cyt)*SHP1) / init_EpoRJAK2 - SHP1Dea*SHP1Act,
    D(npSTAT5) ~ STAT5Imp*pSTAT5 - STAT5Exp*npSTAT5,
    D(p12EpoRpJAK2) ~ (EpoRActJAK2*p2EpoRpJAK2) / (1 + (SOCS3Inh*SOCS3) / (SOCS3Eqc*cyt)) + (3EpoRActJAK2*p1EpoRpJAK2) / ((1 + (EpoRCISInh*EpoRJAK2_CIS) / cyt)*(1 + (SOCS3Inh*SOCS3) / (SOCS3Eqc*cyt))) + (-JAK2EpoRDeaSHP1*SHP1Act*p12EpoRpJAK2) / (cyt*init_SHP1),
    D(p2EpoRpJAK2) ~ (3EpoRActJAK2*EpoRpJAK2) / ((1 + (EpoRCISInh*EpoRJAK2_CIS) / cyt)*(1 + (SOCS3Inh*SOCS3) / (SOCS3Eqc*cyt))) + (-EpoRActJAK2*p2EpoRpJAK2) / (1 + (SOCS3Inh*SOCS3) / (SOCS3Eqc*cyt)) + (-JAK2EpoRDeaSHP1*SHP1Act*p2EpoRpJAK2) / (cyt*init_SHP1),
    D(CIS) ~ (CISEqc*CISTurn*CISRNA) / CISRNAEqc + CISEqc*CISEqcOE*CISTurn*CISoe*cyt - CISTurn*CIS,
    D(EpoRpJAK2) ~ (Epo*JAK2ActEpo*EpoRJAK2) / (1 + (SOCS3Inh*SOCS3) / (SOCS3Eqc*cyt)) + (-EpoRActJAK2*EpoRpJAK2) / (1 + (SOCS3Inh*SOCS3) / (SOCS3Eqc*cyt)) + (-3EpoRActJAK2*EpoRpJAK2) / ((1 + (EpoRCISInh*EpoRJAK2_CIS) / cyt)*(1 + (SOCS3Inh*SOCS3) / (SOCS3Eqc*cyt))) + (-JAK2EpoRDeaSHP1*EpoRpJAK2*SHP1Act) / (cyt*init_SHP1),
    D(CISnRNA2) ~ CISRNADelay*CISnRNA1 - CISRNADelay*CISnRNA2,
    D(CISRNA) ~ CISRNADelay*CISnRNA5 - CISRNATurn*CISRNA
    ]

    @named sys = ODESystem(eqs, t, stateArray, parameterArray)

    ### Initial species concentrations ###
    initialSpeciesValues = [
    p1EpoRpJAK2 => 0.0,
    pSTAT5 => 0.0,
    EpoRJAK2_CIS => init_EpoRJAK2_CIS,
    SOCS3nRNA4 => 0.0,
    SOCS3RNA => 0.0,
    SHP1 => init_SHP1*(1 + SHP1ProOE*init_SHP1_multiplier),
    STAT5 => init_STAT5,
    EpoRJAK2 => init_EpoRJAK2,
    CISnRNA1 => 0.0,
    SOCS3nRNA1 => 0.0,
    SOCS3nRNA2 => 0.0,
    CISnRNA3 => 0.0,
    CISnRNA4 => 0.0,
    SOCS3 => SOCS3Eqc*SOCS3EqcOE*init_SOCS3_multiplier,
    CISnRNA5 => 0.0,
    SOCS3nRNA5 => 0.0,
    SOCS3nRNA3 => 0.0,
    SHP1Act => 0.0,
    npSTAT5 => 0.0,
    p12EpoRpJAK2 => 0.0,
    p2EpoRpJAK2 => 0.0,
    CIS => CISEqc*CISEqcOE*init_CIS_multiplier,
    EpoRpJAK2 => 0.0,
    CISnRNA2 => 0.0,
    CISRNA => 0.0
    ]

    ### SBML file parameter values ###
    trueParameterValues = [
    STAT5Exp => 0.0745150819016423,
    STAT5Imp => 0.0268865083829685,
    init_SOCS3_multiplier => 0.0,
    EpoRCISRemove => 5.42980693903448,
    STAT5ActEpoR => 38.9957991073948,
    SHP1ActEpoR => 0.00100000000000006,
    JAK2EpoRDeaSHP1 => 142.72332309738,
    CISTurn => 0.0083988695167017,
    SOCS3Turn => 9999.99999999912,
    init_EpoRJAK2_CIS => 0.0,
    SOCS3Inh => 10.4078649133666,
    ActD => 1.25e-7,
    init_CIS_multiplier => 0.0,
    CISRNAEqc => 1.0,
    JAK2ActEpo => 633167.430600806,
    Epo => 1.25e-7,
    SOCS3oe => 1.25e-7,
    CISInh => 7.85269991450496e8,
    SHP1Dea => 0.00816220490950374,
    SOCS3EqcOE => 0.679165515556864,
    CISRNADelay => 0.14477775532111,
    init_SHP1 => 26.7251164277109,
    CISEqcOE => 0.530264447119609,
    EpoRActJAK2 => 0.267304849333058,
    SOCS3RNAEqc => 1.0,
    CISEqc => 432.860413434913,
    SHP1ProOE => 2.82568153411555,
    SOCS3RNADelay => 1.06458446742251,
    init_STAT5 => 79.75363993771,
    CISoe => 1.25e-7,
    CISRNATurn => 999.999999999946,
    init_SHP1_multiplier => 1.0,
    init_EpoRJAK2 => 3.97622369384192,
    EpoRCISInh => 999999.999999912,
    STAT5ActJAK2 => 0.0781068855795467,
    SOCS3RNATurn => 0.00830917643120369,
    SOCS3Eqc => 173.64470023136,
    cyt => 0.4,
    nuc => 0.275
    ]

    return sys, initialSpeciesValues, trueParameterValues

end
