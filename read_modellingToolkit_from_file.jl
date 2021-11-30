using ModelingToolkit, OrdinaryDiffEq

modelName = "model_Beer_MolBioSystems2014"
include("./RewrittenModels/" * modelName * ".jl")

sys = ode_order_lowering(sys)

u0 = initialSpeciesValues

p = trueParameterValues

c = trueConstantsValues

tspan = (0.0,10.0)
prob = ODEProblem(sys,u0,tspan,[p;c],jac=true)
sol = solve(prob,TRBDF2())
using Plots
plotly()
#plt1 = plot(sol,vars=(Cells))
plt2 = plot(sol)
display(plt2)