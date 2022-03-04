using DifferentialEquations
using ForwardDiff
using Random
using BenchmarkTools
using Plots
using ProgressMeter


# Struct containing all necessary info
mutable struct Adam
    theta::AbstractArray{Float64} # Parameter array
    loss::Function                # Loss function
    grad::Function                # Gradient function
    m::AbstractArray{Float64}     # First moment
    v::AbstractArray{Float64}     # Second moment
    b1::Float64                   # Exp. decay first moment
    b2::Float64                   # Exp. decay second moment
    a::Float64                    # Step size
    eps::Float64                  # Epsilon for stability
    t::Int                        # Time step (iteration)
end
  
# Outer constructor
function Adam(theta::AbstractArray{Float64}, loss::Function, grad::Function; a::Float64=0.005)
    m   = zeros(size(theta))
    v   = zeros(size(theta))
    b1  = 0.9
    b2  = 0.999
    a   = a
    eps = 1e-8
    t   = 0
    Adam(theta, loss, grad, m, v, b1, b2, a, eps, t)
end
  

  # Step function with optional keyword arguments for the data passed to grad()
function step!(opt::Adam; data...)
    opt.t += 1
    gt    = opt.grad(opt.theta; data...)
    opt.m = opt.b1 .* opt.m + (1 - opt.b1) .* gt
    opt.v = opt.b2 .* opt.v + (1 - opt.b2) .* gt .^ 2
    mhat = opt.m ./ (1 - opt.b1^opt.t)
    vhat = opt.v ./ (1 - opt.b2^opt.t)
    opt.theta -= opt.a .* (mhat ./ (sqrt.(vhat) .+ opt.eps))
end


function simple_ode!(du, u, p, t)
    p1, p2 = p
    du[1] = -p1*u[1] + p2*u[2] 
    du[2] = p1*u[1] - p2*u[2] 
end


# F - function on p 
function unbiased_grad(p::T1, f::Function)::Array{Float64, 1} where T1<:Array{<:AbstractFloat, 1}

    v::Array{Float64, 1} = randn(length(p)) # v ~ MV(0, I)

    g = r -> f(p + r*v)
    return ForwardDiff.derivative(g, 0.0) * v
end


function get_data(time_span, p, u0, time_save)

    ode_prob = ODEProblem(simple_ode!, u0, time_span, p)
    sol = solve(ode_prob, saveat=time_save, solver=Rosenbrock23(), abstol=1e-6, reltol=1e-6)

    data_ret = sol[2, :] + randn(length(time_save)) * 0.2

    return data_ret
end


function calc_cost(p, ode_prob, data)

    _ode_prob = remake(ode_prob, u0=convert.(eltype(p),ode_prob.u0), p=p)
    sol = solve(_ode_prob, solver=Rosenbrock23(), saveat=1.0:0.1:15.0, abstol=1e-6, reltol=1e-6)

    cost = sum((data .- sol[2, :]).^2)

    return cost
end


function test_adam()

    Random.seed!(123)

    u0 = [50.0, 100.0]
    p = [0.1, 0.2]
    time_span = [1.0, 15.0]
    time_save = 1.0:0.1:15.0
    data_test = get_data(time_span, p, u0, time_save)

    ode_prob = ODEProblem(simple_ode!, u0, time_span, p)

    f = (x) -> calc_cost(x, ode_prob, data_test)
    # Unbiased gradient estimate 
    grad_g = (x) -> unbiased_grad(x, f)

    n_it = 15000
    theta_save = Array{Float64, 2}(undef, (2, n_it))
    loss_val = Array{Float64, 1}(undef, n_it)

    p::Array{Float64, 1} = [2.0, 4.0]
    adam_opt = Adam(p, f, grad_g, a=0.001)
    
    @showprogress 1 "Running ADAM... " for i in 1:n_it
        step!(adam_opt)
        theta_save[:, i] .= adam_opt.theta
        loss_val[i] = f(adam_opt.theta)
    end

    return theta_save, loss_val
end

theta_i, f_i = test_adam()

p1 = plot(theta_i[1, :], label = "c1")
p1 = plot!(theta_i[2, :], label = "c2")