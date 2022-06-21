using ModelingToolkit, DifferentialEquations, BenchmarkTools, DataFrames, CSV, LSODA, Sundials, ODEInterface, ODEInterfaceDiffEq
using Plots
using StatsBase
using Random
plotly()

include(joinpath(pwd(), "Pipeline_ModelSolver", "BigFloatODEProblem.jl"))
include(joinpath(pwd(), "Additional_functions", "additional_tools.jl"))
include(joinpath(pwd(), "Additional_functions", "Solver_info.jl"))
include(joinpath(pwd(), "Additional_functions", "benchmarkSolvers.jl"))

allModelFunctionVector = includeAllModels(getModelFiles(joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels")), 
        joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels"))


function viewModelDynamics(modelFile, solver, tol)
    
    modelName = replace.(modelFile, ".jl" => "")
    readPath = joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels")
    writePath = joinpath(pwd(), "Pipeline_ModelSolver", "IntermediaryResults")
    fixDirectories(writePath)
        
    allModelFiles = getModelFiles(readPath)    
    usedModelFunctionVector = allModelFunctionVector[[allModelFile .== modelFile for allModelFile in allModelFiles]][1]
                        
    ode_prob = modelSolver(usedModelFunctionVector, modelName, solver, tol)

    return ode_prob
end


function modelSolver(modelFunction, modelName, solver, tol)
    
    sys, initialSpeciesValues, trueParameterValues = modelFunction()
    new_sys = ode_order_lowering(sys)

    # Get informaiton on conditions for model simulations 
    modelPath = joinpath(pwd(), "Pipeline_SBMLImporter", "JuliaModels")
    modelFile = modelName * ".jl" 
    experimentalConditions, measurementData, parameterBounds = readDataFiles(modelName)
    paramData = processParameterData(parameterBounds) # Model data in convient structure 
    # Set stateMap and paramMap non condition parameter to the nominal reporeted values n parameter-files
    stateMap = initialSpeciesValues
    paramMap = trueParameterValues     
    setParamToParamFileVal!(paramMap, stateMap, paramData)

    # Get data on experimental conditions 
    firstExpIds, shiftExpIds, simulateSS, parameterNames, stateNames = getSimulationInfo(measurementData, sys)
    
    # Set up for first experimental condtition 
    changeToCondUse! = (pVec, u0Vec, expID) -> changeToCond!(pVec, u0Vec, expID, paramData, experimentalConditions, parameterNames, stateNames, paramMap, stateMap)
    prob = ODEProblem(new_sys, stateMap, (0.0, 5e3), paramMap, jac=true)
    prob = remake(prob, p = convert.(Float64, prob.p), u0 = convert.(Float64, prob.u0))
    
    # Base soluation on wheter or not steady state pre simulatioln occurs 
    solArray, success = solveOdeModelAllCond(prob, changeToCondUse!, simulateSS, measurementData, firstExpIds, shiftExpIds, tol, solver, nTSave=100)
    println("Success = ", success)

    sqErr = calcSqErr(prob, changeToCondUse!, solArray, simulateSS, measurementData, firstExpIds, shiftExpIds, tol, Rosenbrock23())

    GC.gc()    
        
    return prob, sys, solArray, firstExpIds, shiftExpIds, sqErr
end

# TODO: Succesfull integration on terminated

Random.seed!(123)
modelList = ["model_Beer_MolBioSystems2014.jl", "model_Weber_BMC2015.jl", "model_Schwen_PONE2014.jl", "model_Alkan_SciSignal2018.jl", 
    "model_Bachmann_MSB2011.jl", "model_Bertozzi_PNAS2020.jl", "model_Blasi_CellSystems2016.jl", "model_Boehm_JProteomeRes2014.jl", 
    "model_Borghans_BiophysChem1997.jl", "model_Brannmark_JBC2010.jl", "model_Bruno_JExpBot2016.jl", "model_Crauste_CellSystems2017.jl", 
    "model_Elowitz_Nature2000.jl", "model_Fiedler_BMC2016.jl", "model_Fujita_SciSignal2010.jl", "model_Giordano_Nature2020.jl", 
    "model_Isensee_JCB2018.jl", "model_Laske_PLOSComputBiol2019.jl", "model_Lucarelli_CellSystems2018.jl", "model_Okuonghae_ChaosSolitonsFractals2020.jl", 
    "model_Oliveira_NatCommun2021.jl", "model_Perelson_Science1996.jl", "model_Rahman_MBS2016.jl", "model_Raimundez_PCB2020.jl", 
    "model_SalazarCavazos_MBoC2020.jl", "model_Sneyd_PNAS2002.jl", "model_Zhao_QuantBiol2020.jl", "model_Zheng_PNAS2012.jl"]

prob, sys, solArray, firstExpIds, shiftExpIds, sqErr = viewModelDynamics("model_Blasi_CellSystems2016.jl", Rodas4P(), 1e-9)

a = 1

