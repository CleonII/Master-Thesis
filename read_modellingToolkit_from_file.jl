using ModelingToolkit, OrdinaryDiffEq

include("./testMTK.jl")

sys = ode_order_lowering(sys)

u0 = initialSpeciesValues

p = trueParameterValues

c = trueConstantsValues

tspan = (0.0,10.0)
prob = ODEProblem(sys,u0,tspan,[p;c],jac=true)
sol = solve(prob,TRBDF2())