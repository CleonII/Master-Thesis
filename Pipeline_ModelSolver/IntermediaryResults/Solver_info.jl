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