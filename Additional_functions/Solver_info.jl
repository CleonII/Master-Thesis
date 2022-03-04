# TODO: fix sundials: KLU, needs a sparce jacobian
function getSolvers(;algSets = ["all"], algNum = nothing, algHintSet = true, algHintNum = nothing)
    algs_DiffEq = [Vern6(), Vern7(), Vern8(), Tsit5(), DP5(), DP8(), Feagin14(), VCABM(),
                   Rosenbrock23(), TRBDF2(), Rodas4(), Rodas4P(), Rodas4P2(), Rodas5(), QNDF(), FBDF(), 
                   Trapezoid(), KenCarp4(), Kvaerno5(), RadauIIA3(), RadauIIA5(), 
                   AutoTsit5(Rosenbrock23()), AutoVern7(Rodas5()), AutoVern9(Rodas4P()), AutoVern9(Rodas5())]
    algs_LSODA = [lsoda()]
    algs_Sundials = [CVODE_BDF(linear_solver=:Dense), CVODE_BDF(linear_solver=:LapackDense), CVODE_BDF(linear_solver=:GMRES), 
                     CVODE_Adams(linear_solver=:Dense), CVODE_Adams(linear_solver=:LapackDense), 
                     ARKODE(Sundials.Explicit(), order=4), ARKODE(Sundials.Explicit(), order=8), 
                     ARKODE(Sundials.Implicit(), order = 3), ARKODE(Sundials.Implicit(), order = 5)]
    algs_ODEInterface = [dopri5(), dop853(), radau(), radau5(), rodas()]
    algHints = [[:auto], [:nonstiff], [:stiff]]
    if algSets == ["all"]
        algs = [algs_DiffEq; algs_LSODA; algs_Sundials; algs_ODEInterface]
    else
        algs = Array{SciMLBase.AbstractODEAlgorithm, 1}(undef, 0)
        for algSet in algSets 
            if algSet == "algs_DiffEq"
                append!(algs, algs_DiffEq)
            elseif algSet == "algs_LSODA"
                append!(algs, algs_LSODA)
            elseif algSet == "algs_Sundials"
                append!(algs, algs_Sundials)
            elseif algSet == "algs_ODEInterface"
                append!(algs, algs_ODEInterface)
            else
                println("Error: No such set")
            end
        end
    end
    if algNum !== nothing
        algs = algs[algNum]
    end
    if algHintNum !== nothing
        algHints = algHints[algHintNum]
    end
    if algHintSet == true
        return algs, algHints
    else
        return algs
    end
end

function getHiAccSolver()
    algs = [AutoVern9(Rodas5()), Rodas4P()]
    return algs
end

function getSolverType()
    ["nonstiff", "nonstiff", "nonstiff", "nonstiff", "nonstiff", "nonstiff", "nonstiff", "nonstiff", 
    "stiff", "stiff", "stiff", "stiff", "stiff", "stiff", "stiff", "stiff",
    "stiff", "stiff", "stiff", "stiff", "stiff", 
    "composite", "composite", "composite", "composite",
    "composite",
    "stiff", "stiff", "stiff", 
    "nonstiff", "nonstiff", 
    "nonstiff", "nonstiff",
    "stiff", "stiff", 
    "nonstiff", "nonstiff", "stiff", "stiff", "stiff",
    "hint", "hint", "hint"]
end

function getSolverLibrary()
    ["OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", 
    "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", 
    "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", 
    "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", 
    "LSODA.jl", 
    "Sundials", "Sundials", "Sundials", 
    "Sundials", "Sundials", 
    "Sundials", "Sundials", 
    "Sundials", "Sundials", 
    "ODEInterface", "ODEInterface", "ODEInterface", "ODEInterface", "ODEInterface", 
    "DifferentialEquations", "DifferentialEquations", "DifferentialEquations"]
end