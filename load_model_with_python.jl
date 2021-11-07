# in python consol:
# >>> import pip
# >>> import subprocess
# >>> import sys
# >>> def install(package):
# ...     subprocess.check_call([sys.executable, "-m", "pip", "install", package])
# ...
# >>> install(python-libsbml)

using PyCall

libsbml = pyimport("libsbml")
reader = libsbml.SBMLReader()
document = reader[:readSBML]("Master-Thesis/reactionsystem_02.xml")
model = document[:getModel]() # Get the model
species = [(s[:getName](), s[:getInitialConcentration](), s[:getCompartment]()) for s in model[:getListOfSpecies]()] # Use getName to extract the species names from getListOfSpecies, return an array.
reactions = [(r[:getName](), r[:getKineticLaw]()[:getFormula]()) for r in model[:getListOfReactions]()] # Get the names of all reactions in the SBML file. 
parameters = [(p[:getId](), p[:getValue]()) for p in model[:getListOfParameters]()] # Get a list of all parameters and their values.

println(model[:getReaction](0)) # Select reaction 3 (python indexing)
println(model[:getReaction](0)[:getKineticLaw]()[:getFormula]()) # Get a string for the mathematical formula (rate law). 
println(model[:getReaction](0)[:getListOfReactants]()[1][:getStoichiometry]())




fd = model[:getListOfFunctionDefinitions]()[1]
math = fd[:getMath]()
if math[:getNumChildren]() > 1
    print("(" )
    print(math[:getLeftChild]()[:getName]())
    for n in range(1, math[:getNumChildren]()-1, step = 1)
        arg = math[:getChild](n)[:getName]()
        if arg !== nothing
            print(", " * arg)
        end
    end
    println(") = ")
end
fd[:getArgument](0)[:getName]()

libsbml[:formulaToString](math) # Important


functionDefinitions = [(f[:getName](), f[:getMetaId](), f[:getSBOTerm](), f[:toXMLNode](), f[:getTypeCode](), f[:getMath]()) for f in model[:getListOfFunctionDefinitions]()]
tmp = [(f[:getAnnotation](), f[:getCVTerms]()) for f in model[:getListOfFunctionDefinitions]()]
funcMath = [(f[:getName](), f[:getMath]()) for f in model[:getListOfFunctionDefinitions]()]

reac = model[:getReaction](0)
reac[:getCompartment]()
reactants = reac[:getListOfReactants]()
spec = [s[:getName] for s in reac[:getListOfReactants]()]

reac[:getListOfAllElements]()
reac[:getListOfProducts]()
prod = [p[:species] for p in reac[:getListOfProducts]()]
reactants = [r[:species] for r in reac[:getListOfReactants]()]

model[:getListOfAllElements]() 

compartments = [comp[:getName]() for comp in model[:getListOfCompartments]()]
constraints = [con[:getName]() for con in model[:getListOfConstraints]()]
rules = [rule[:getName]() for rule in model[:getListOfRules]()]
units = [u[:getName]() for u in model[:getListOfUnitDefinitions]()]
