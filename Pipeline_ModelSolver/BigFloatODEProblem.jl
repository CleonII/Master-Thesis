

function BigFloatODEProblem(sys::ODESystem, u0, tspan, par)::ODEProblem
    par = convert(Vector{Pair{Num, BigFloat}}, par)

    if typeof(u0) == Vector{Pair{Num, Num}}
        # This means that u0 has some initializations that depend on par, 
        # since the values of par is BigFloat, u0 will be BigFloat
    elseif typeof(u0) == Vector{Pair{Num, Float64}}
        u0 = convert(Vector{Pair{Num, BigFloat}}, u0)
    else
        println("Error: something went wrong!")
    end

    prob = ODEProblem(sys,u0,tspan,par,jac=true)
    return prob
end