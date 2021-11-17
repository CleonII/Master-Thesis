using BenchmarkTools
using OrdinaryDiffEq
#using Plots
include("./test.jl")

u0 = initialSpeciesValue
du = Vector{Float64}(undef, length(u0))
p = trueParameterValues
tspan = [0.0, 10.0]
prob = ODEProblem(ODE_model!, u0, tspan, p)
sol = @btime OrdinaryDiffEq.solve(prob, TRBDF2(), reltol=1e-18, abstol=1e-18, dtmax = 0.01)
print("done")
#plt1 = plot(sol)
#display(plt1)