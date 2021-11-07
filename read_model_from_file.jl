using PyCall
using BenchmarkTools
using OrdinaryDiffEq
using Plots
include("test.jl")

du = Vector{Float64}(undef, 3)
u0 = [1.0, 1.0, 1.0]
p = 1.0
tspan = [0.0, 10.0]
prob = ODEProblem(ODE_model!, u0, tspan, p)
sol = @btime OrdinaryDiffEq.solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
plt1 = plot(sol)
display(plt1)