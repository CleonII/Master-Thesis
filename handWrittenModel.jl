using ModelingToolkit, OrdinaryDiffEq, BenchmarkTools
# Handwritten model of model_Boehm_JProteomeRes2014.xml

@parameters Epo_degradation_BaF3 k_exp_hetero k_exp_homo k_imp_hetero k_imp_homo k_phos ratio specC17 
@parameters cyt nuc
@variables t STAT5A(t) STAT5B(t) pApB(t) pApA(t) pBpB(t) nucpApA(t) nucpApB(t) nucpBpB(t)

D = Differential(t)

BaF3_Epo(t, Epo_degradation_BaF3) = 1.25 * 10^(-7) * exp(-1 * Epo_degradation_BaF3 * t)
@register BaF3_Epo(t, Epo_degradation_BaF3)

eqs = [
    D(STAT5A) ~ -2 * (cyt * BaF3_Epo(t, Epo_degradation_BaF3) * STAT5A^2 * k_phos)-1 * (cyt * BaF3_Epo(t, Epo_degradation_BaF3) * STAT5A * STAT5B * k_phos)+2 * (nuc * k_exp_homo * nucpApA)+1 * (nuc * k_exp_hetero * nucpApB),
    D(STAT5B) ~ -1 * (cyt * BaF3_Epo(t, Epo_degradation_BaF3) * STAT5A * STAT5B * k_phos)-2 * (cyt * BaF3_Epo(t, Epo_degradation_BaF3) * STAT5B^2 * k_phos)+1 * (nuc * k_exp_hetero * nucpApB)+2 * (nuc * k_exp_homo * nucpBpB),
    D(pApB) ~ +1 * (cyt * BaF3_Epo(t, Epo_degradation_BaF3) * STAT5A * STAT5B * k_phos)-1 * (cyt * k_imp_hetero * pApB),
    D(pApA) ~ +1 * (cyt * BaF3_Epo(t, Epo_degradation_BaF3) * STAT5A^2 * k_phos)-1 * (cyt * k_imp_homo * pApA),
    D(pBpB) ~ +1 * (cyt * BaF3_Epo(t, Epo_degradation_BaF3) * STAT5B^2 * k_phos)-1 * (cyt * k_imp_homo * pBpB),
    D(nucpApA) ~ +1 * (cyt * k_imp_homo * pApA)-1 * (nuc * k_exp_homo * nucpApA),
    D(nucpApB) ~ +1 * (cyt * k_imp_hetero * pApB)-1 * (nuc * k_exp_hetero * nucpApB),
    D(nucpBpB) ~ +1 * (cyt * k_imp_homo * pBpB)-1 * (nuc * k_exp_homo * nucpBpB)
]

@named sys = ODESystem(eqs)
sys = ode_order_lowering(sys)

u0 = [
    STAT5A => 207.6 * ratio,
    STAT5B => 207.6 - 207.6 * ratio,
    pApB => 0.0,
    pApA => 0.0,
    pBpB => 0.0,
    nucpApA => 0.0,
    nucpApB => 0.0,
    nucpBpB => 0.0
]

p = [
    Epo_degradation_BaF3 => 0.0269738286367359,
    k_exp_hetero => 1.00094251286741e-05,
    k_exp_homo => 0.00617193081581346,
    k_imp_hetero => 0.0163708512310568,
    k_imp_homo => 96945.5391768823,
    k_phos => 15766.8336642826,
    ratio => 0.693,
    specC17 => 0.107
]

c = [
    cyt => 1.4,
    nuc => 0.45
]

tspan = (0.0,10.0)
prob = @btime ODEProblem(sys,u0,tspan,[p;c],jac=true)
sol = @btime solve(prob,Rosenbrock23())
using Plots 
plot(sol)