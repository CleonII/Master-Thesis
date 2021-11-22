
using ModelingToolkit, OrdinaryDiffEq

@parameters t σ ρ β
@variables x(t) y(t) z(t) c(t) d(t)
D = Differential(t)

continuous_events = [ # This time we have a vector of pairs
    [t ~ 10] => [c ~ 10, d ~ 10]; [t ~ 11, t ~ 12] => [d ~ 3]
]
# write rules either as events or as functions 


fun_z(t, a, b) = t < 10 ? a : b
fun1(d) = begin; if d > 3; 1; else; 2; end; end
@register fun_z(t,a,b)
@register fun1(d)


eqs = [D(D(x)) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z + c*d * fun1(d),
       D(c) ~ 0,
       D(d) ~ 0]

@named sys = ODESystem(eqs, t, continuous_events = continuous_events)
sys = ode_order_lowering(sys)

u0 = [D(x) => 2.0,
      x => 1.0,
      y => 0.0,
      z => 0.0,
      c => 0.0,
      d => 5.0]

p  = [σ => 28.0,
      ρ => 10.0,
      β => 8/3]


tspan = (0.0,100.0)
prob = ODEProblem(sys,u0,tspan,p,jac=true)
sol = @time solve(prob,Tsit5())
using Plots; plot(sol,vars=(x,y))
tUse = (sol.t .> 8) .== (sol.t .< 12)
plot(sol.t[tUse], sol[4,tUse])

# piecewise function?
# begin if SN38_level > 0; 1; elseif SN38_level <= 0; 0; else; 0; end; end

#=
using ModelingToolkit, OrdinaryDiffEq, Plots
@variables t x(t)=1 y(t)=0 vx(t)=0 vy(t)=2
D = Differential(t)

continuous_events = [ # This time we have a vector of pairs
    [x ~ 0] => [vx ~ -vx]
    [y ~ -1.5, y ~ 1.5] => [vy ~ -vy]
]

@named ball = ODESystem([
    D(x)  ~ vx,
    D(y)  ~ vy,
    D(vx) ~ -9.8-0.1vx, # gravity + some small air resistance
    D(vy) ~ -0.1vy,
], t, continuous_events = continuous_events)


ball = structural_simplify(ball)

tspan = (0.0,10.0)
prob = ODEProblem(ball, Pair[], tspan)

sol = solve(prob,Tsit5())
@assert 0 <= minimum(sol[x]) <= 1e-10 # the ball never went through the floor but got very close
@assert minimum(sol[y]) > -1.5 # check wall conditions
@assert maximum(sol[y]) < 1.5  # check wall conditions

tv = sort([LinRange(0, 10, 200); sol.t])
plot(sol(tv)[y], sol(tv)[x], line_z=tv)
vline!([-1.5, 1.5], l=(:black, 5), primary=false)
hline!([0], l=(:black, 5), primary=false)
=#