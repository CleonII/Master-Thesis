
using ModelingToolkit, DifferentialEquations, BenchmarkTools, DataFrames, CSV, LSODA, Sundials, ODEInterface, ODEInterfaceDiffEq
using Plots
plotly()


include(joinpath(pwd(), "Pipeline_ModelSolver", "BigFloatODEProblem.jl"))
include(joinpath(pwd(), "Additional_functions", "additional_tools.jl"))
include(joinpath(pwd(), "Additional_functions", "Solver_info.jl"))

allModelFunctionVector = includeAllModels(getModelFiles(joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels")), 
        joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels"))


function modelSolver(modelFunction, timeEnd, solver, tol)
    sys, initialSpeciesValues, trueParameterValues = modelFunction()
    new_sys = ode_order_lowering(sys)
    u0 = initialSpeciesValues
    pars = trueParameterValues 
    tspan = (0.0, timeEnd)

    prob = ODEProblem(new_sys,u0,tspan,pars,jac=true)

    sol = solve(prob, solver, reltol = tol, abstol = tol)

    #plt = plot(sol)
    #display(plt)
end



function viewModelDynamics(modelFile, solver, tol)
    readPath = joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels")
    writePath = joinpath(pwd(), "Pipeline_ModelSolver", "IntermediaryResults")
    fixDirectories(writePath)
        
    allModelFiles = getModelFiles(readPath)
        
    usedModelFunctionVector = allModelFunctionVector[[allModelFile .== modelFile for allModelFile in allModelFiles]][1]
            
    timeEnds = CSV.read(joinpath(writePath, "timeScales.csv"), DataFrame)
    timeEnd = timeEnds[timeEnds[:,1] .== modelFile, 2][1]
            
    modelSolver(usedModelFunctionVector, timeEnd, solver, tol)
end

viewModelDynamics("model_Chen_MSB2009.jl", Rodas4P(), 1e-9)