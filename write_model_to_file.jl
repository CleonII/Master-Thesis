using PyCall
using BenchmarkTools
using OrdinaryDiffEq
using Plots

#ENV["PYTHON"] = "H:\\Program\\Python37\\python.exe"
#Pkg.build("PyCall")

libsbml = pyimport("libsbml")
reader = libsbml.SBMLReader()
document = reader[:readSBML]("Master-Thesis/reactionsystem_01.xml")
model = document[:getModel]() # Get the model

function writeODEModelToFile()

    f = open("Master-Thesis/test.jl", "w")
    println(f, "function ODE_model!(du, u, p, t)")

    parameters = ""
    for (pIndex, par) in enumerate(model[:getListOfParameters]())
        if pIndex == 1
            parameters = par[:getName]()
        else
            parameters = parameters * ", " * par[:getName]()
        end
    end
    parametersValue = parameters * " = p"
    println(f, parametersValue)

    compartmentsSize = [c[:getName]() * " = " * string(c[:getSize]()) for c in model[:getListOfCompartments]()]
    for compSize in compartmentsSize
        println(f, compSize)
    end

    species = ""
    for (sIndex, spec) in enumerate(model[:getListOfSpecies]())
        if sIndex == 1
            species = spec[:getName]()
        else
            species = species * ", " * spec[:getName]()
        end
    end
    speciesValue = species * " = u"
    println(f, speciesValue)

    dus = Dict()
    for (index, spec) in enumerate(model[:getListOfSpecies]())
        dus[spec[:getName]()] = "du[" * string(index) * "] = "
    end

    reactions = [(r, r[:getKineticLaw]()[:getFormula]()) for r in model[:getListOfReactions]()]
    for (reac, formula) in reactions
        products = [(p[:species], p[:getStoichiometry]()) for p in reac[:getListOfProducts]()]
        reactants = [(r[:species], r[:getStoichiometry]()) for r in reac[:getListOfReactants]()]
        for (rName, rStoich) in reactants
            dus[rName] = dus[rName] * "-" * string(rStoich) * " * (" * formula * ")"
        end
        for (pName, pStoich) in products
            dus[pName] = dus[pName] * "+" * string(pStoich) * " * (" * formula * ")"
        end
    end

    for spec in model[:getListOfSpecies]()
        println(f, dus[spec[:getName]()])
    end

    println(f, "nothing")
    println(f, "end")
    
    close(f)

end



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


writeODEModelToFile()



