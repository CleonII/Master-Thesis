using ModelingToolkit, OrdinaryDiffEq, BenchmarkTools

modelName = "model_Okuonghae_ChaosSolitonsFractals2020"
include("./RewrittenModels/" * modelName * ".jl")

sys = ode_order_lowering(sys)

u0 = initialSpeciesValues

p = trueParameterValues

c = trueConstantsValues

tspan = (0.0,150.0)
prob = ODEProblem(sys,u0,tspan,[p;c],jac=true)

sol = solve(prob,Rosenbrock23())
#sol = @btime solve(prob,TRBDF2(), dtmax = 0.001)
#sol = solve(prob, Tsit5(), dtmax = 0.001)
using Plots
plotly()
#plt1 = plot(sol,vars=(Cells))
plt2 = plot(sol)
display(plt2)