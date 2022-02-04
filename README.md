# Master-Thesis
To benchmark the ODE-solvers for the models, first run the "main.jl" file in "/Pipeline_SBMLImporter/". 
The models in SBML format should be located in "/Pipeline_SBMLImporter/SBML". 
The main file calls "write_modellingToolkit_to_file.jl" for each model, 
which converts the model into a julia script file with the model written in the ModellingToolkit framework.
The rewritten models will be placed in "/Pipeline_SBMLImporter/JuliaModels", which is created during the run.
For each model, the time span for the model have to be supplied in "timeScales.csv" in "/Pipeline_ModelSolver/IntermediaryResults/"

When the models have been imported, run the "main.jl" file in "/Pipeline_ModelSolver/"
The solving of the models will then be benchmarked with the different solvers specified in the "getSolvers()" functions in the main file.
The solvers can be taken from "DifferentialEquations.jl", "LSODA.jl", "Sundials.jl" and "ODEInterface" (also using "ODEInterfaceDiffEq" to ease use.)
The models will be solved using different relative and absolute tolerances and they will also be solved multiple times for the same setting to get reliable results. 
Each model will be solved with a high accuracy solver, the result from thid will be used to calculate the "error" of each solver. 
For this high accuracy solver, BigFloat will be used.
The solvers will be benchmarked for speed, memory use and number of allocations. 
The result will be saved in a tidy format in a CSV file. This file is stored in "/Pipeline_ModelSolver/IntermediaryResults/", previouse files will not be overwritten.

The enviroment used for the code can be recreated with the "Project.toml" file.



