using SBMLToolkit, ModelingToolkit, OrdinaryDiffEq

tspan = (0., 140.0)
solver = Tsit5()
include("/home/damiano/Master-Thesis/Intermediate/PeTab_models/model_Bachmann_MSB2011/model_Bachmann_MSB2011.jl")
#include("/home/damiano/Master-Thesis/tests/model_Bachmann_test.jl")
sys, initialSpeciesValues, trueParameterValues = getODEModel_model_Bachmann_MSB2011()
prob = ODEProblem(sys, initialSpeciesValues, tspan, trueParameterValues)
trueSol = solve(prob, solver)


modelName = "model_Bachmann_MSB2011"
modelPath = joinpath(pwd(), "Intermediate", "PeTab_models", modelName)
xmlPath = joinpath(modelPath, modelName * ".xml")

mdl = readSBML(xmlPath)
#odesys = ODESystem(mdl)

rs = ReactionSystem(mdl)  # If you want to create a reaction system
odesys = convert(ODESystem, rs)  # Alternatively: ODESystem(mdl)

prob = ODEProblem(odesys, initialSpeciesValues, tspan, trueParameterValues)
newSol = solve(prob, solver)



sum(newSol-trueSol)
