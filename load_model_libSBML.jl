using PyCall
using BenchmarkTools
using OrdinaryDiffEq
using Plots
#ENV["PYTHON"] = "H:\\Program\\Python37\\python.exe"
#Pkg.build("PyCall")

libsbml = pyimport("libsbml")
reader = libsbml.SBMLReader()
document = reader[:readSBML]("Summer-project/reactionsystem_01.xml")
model = document[:getModel]() # Get the model


#1.665650 seconds (2.31 M allocations: 143.981 MiB, 1.32% gc time, 85.42% compilation time)
function ODE_model!(du, u, p, t)
    # reads the speciesNames from the model
    speciesName = [s[:getName]() for s in model[:getListOfSpecies]()]
    # creates a dictionary for the derivatives so that they can be referenced using the speciesName
    dus = Dict(speciesName .=> 0.0)
    # sets the concentrations for the species from u (creates a symbor from a string (e.g. "s1 = 1.0") and then evaluates the symbol)
    speciesConc = speciesName .* " = " .* string.(u)
    eval.(Meta.parse.(speciesConc))
    # sets the parameters according to p (same as above)
    parametersName = [p[:getName]() for p in model[:getListOfParameters]()]
    parameterValue = parametersName .* " = " .* string.(p)
    eval.(Meta.parse.(parameterValue))
    # sets compartmentsSize according to model
    compartmentsSize = [c[:getName]() * " = " * string(c[:getSize]()) for c in model[:getListOfCompartments]()]
    eval.(Meta.parse.(compartmentsSize))
    
    # iterates over the different reactions
    reactions = [(r, r[:getKineticLaw]()[:getFormula]()) for r in model[:getListOfReactions]()]
    for (reac, formula) in reactions
        f_eval = eval(Meta.parse(formula))
        products = [(p[:species], p[:getStoichiometry]()) for p in reac[:getListOfProducts]()]
        reactants = [(r[:species], r[:getStoichiometry]()) for r in reac[:getListOfReactants]()]
        for (rName, rStoich) in reactants
            dus[rName] += -rStoich * f_eval
        end
        for (pName, pStoich) in products
            dus[pName] += pStoich * f_eval
        end
    end
    for (sIndex, sName) in enumerate(speciesName)
        du[sIndex] = dus[sName]
    end
end


#=
function ODE_model!(du, u, p, t)
    # reads the speciesNames from the model
    speciesName = [s[:getName]() for s in model[:getListOfSpecies]()]
    # creates a dictionary for the derivatives so that they can be referenced using the speciesName
    dus = Dict(speciesName .=> 0.0)
    # sets the concentrations for the species from u (creates a symbor from a string (e.g. "s1 = 1.0") and then evaluates the symbol)
    speciesConc = Dict(speciesName .=> u)
    # sets the parameters according to p (same as above)
    parametersName = [p[:getName]() for p in model[:getListOfParameters]()]
    parameterValue = Dict(parametersName .=> p)
    # sets compartmentsSize according to model
    compartmentsSize = Dict(c[:getName]() => c[:getSize]() for c in model[:getListOfCompartments]())
    
    # iterates over the different reactions
    reactions = [(r, r[:getKineticLaw]()[:getFormula]()) for r in model[:getListOfReactions]()]
    for (reac, formula) in reactions
        parts = split(formula)
        products = [(p[:species], p[:getStoichiometry]()) for p in reac[:getListOfProducts]()]
        reactants = [(r[:species], r[:getStoichiometry]()) for r in reac[:getListOfReactants]()]
        for (rName, rStoich) in reactants
            dus[rName] += -rStoich * f_eval
        end
        for (pName, pStoich) in products
            dus[pName] += pStoich * f_eval
        end
    end
    for (sIndex, sName) in enumerate(speciesName)
        du[sIndex] = dus[sName]
    end
end
=#

du = Vector{Float64}(undef, 3)
u0 = [1.0, 1.0, 1.0]
p = [1.0]
tspan = [0.0, 10.0]
prob = ODEProblem(ODE_model!, u0, tspan, p)
sol = @btime OrdinaryDiffEq.solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
plt1 = plot(sol)
display(plt1)

# 0.913854 seconds (1.22 M allocations: 72.285 MiB, 2.35% gc time, 99.55% compilation time)
function ODE_model2!(du, u, p, t)
    du[1] = - 2 * p[1] * u[1] * u[3]
    du[2] = 2 * p[1] * u[1] * u[3]
    du[3] = - 2 * p[1] * u[1] * u[3]
end

prob2 = ODEProblem(ODE_model2!, u0, tspan, p)
sol2 = @btime OrdinaryDiffEq.solve(prob2, Tsit5(), reltol=1e-8, abstol=1e-8)
plt2 = plot(sol2)
display(plt2)

