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
document = reader[:readSBML]("./b3.xml")
model = document[:getModel]() # Get the model
species = [(s[:getName](), s[:getInitialConcentration](), s[:getCompartment]()) for s in model[:getListOfSpecies]()] # Use getName to extract the species names from getListOfSpecies, return an array.
reactions = [(r[:getName](), r[:getKineticLaw]()[:getFormula]()) for r in model[:getListOfReactions]()] # Get the names of all reactions in the SBML file. 
parameters = [(p[:getId](), p[:getValue]()) for p in model[:getListOfParameters]()] # Get a list of all parameters and their values.

println(model[:getReaction](0)) # Select reaction 3 (python indexing)
println(model[:getReaction](0)[:getKineticLaw]()[:getFormula]()) # Get a string for the mathematical formula (rate law). 
println(model[:getReaction](0)[:getListOfReactants]()[1][:getStoichiometry]())


# Function definitions

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

tmp = libsbml[:formulaToString](math) # Important
tmp2 = split(split(tmp, '(')[2], ')')[1]
func = lstrip(split(tmp2, ',')[end])
args = tmp2[1:(end-length(func)-2)]
funcName = fd[:getId]()
funcName * "(" * args * ") = " * func


# rules
rule = model[:getListOfRules]()[7]
rule[:getElementName]() # type
rule[:getVariable]() # variable
rule[:getFormula]() # need to add piecewise function (and others)


# events
event = model[:getListOfEvents]()[1]
event[:getElementName]()
event[:getName]()
trigger = event[:getTrigger]()
triggerMath = trigger[:getMath]()
triggerFormula = libsbml[:formulaToString](triggerMath)
elements = event[:getListOfAllElements]()
elements = event[:getListOfEventAssignments]()
elements[1]


reac = model[:getReaction](0)
reac[:getCompartment]()
reactants = reac[:getListOfReactants]()
spec = [s[:getName] for s in reac[:getListOfReactants]()]

reac[:getListOfAllElements]()
reac[:getListOfProducts]()
prod = [p[:species] for p in reac[:getListOfProducts]()]
reactants = [r[:species] for r in reac[:getListOfReactants]()]

model[:getListOfAllElements]() 

constraints = [con[:getId]() for con in model[:getListOfConstraints]()]
units = [u[:getId]() for u in model[:getListOfUnitDefinitions]()]
