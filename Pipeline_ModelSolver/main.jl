using ModelingToolkit, DifferentialEquations, BenchmarkTools, DataFrames, CSV

# TODO:
# get time by finding time to steady state
# use a high precission solver to estimate the ground truth for each model

function getModelFiles(path)
    if isdir(path)
        modelFiles = readdir(path)
        return modelFiles[1:2]
    else
        println("No such directory")
        return nothing
    end
end

function getSolvers()
    solvers = [Vern6(), Tsit5()]
    return solvers
end

function getTolerances()
    relTols = [1e-5, 1e-8, 1e-14]
    absTols = [1e-5, 1e-8, 1e-14]
    return relTols, absTols
end

function fixDirectories(path)
    if ~isdir(path)
        println("Create directory: " * path)
        mkdir(path)
    end
end

function getNumberOfFiles(path)
    if isdir(path)
        files = readdir(path)
        return length(files)
    else
        mkdir(path)
        return 0
    end
end

function modelSolver(modelFile, solver, relTol, absTol, iterations, readPath, writefile)
    modelPath = readPath * "\\" * modelFile
    include(modelPath)

    new_sys = ode_order_lowering(sys)
    u0 = initialSpeciesValues

    p = trueParameterValues

    c = trueConstantsValues

    steadyStateProblem = SteadyStateProblem{true}(new_sys,u0,[p;c])

    #steadyStateProblem = ODEProblem(new_sys, u0, [p;c], jac=true)

    simulationTime = solve(steadyStateProblem, DynamicSS(Tsit5()))

    tspan = (0.0,500.0)
    prob = ODEProblem(new_sys,u0,tspan,[p;c],jac=true)

    runTime = Vector{Float64}(undef, iterations)

    for i in 1:iterations
        println("iteration: " * string(i))
        b = @benchmark solve($prob, $solver, reltol = $relTol, abstol = $absTol) 
        bm = mean(b)
        runTime[i] = bm.time 
    end
    data = DataFrame(model=modelFile, solver=solver, relTol=relTol, absTol=absTol, runTime=runTime, iteration=1:iterations)
    if isfile(writefile)
        CSV.write(writefile, data, append = true)
    else
        CSV.write(writefile, data)
    end
end

function modelSolverIterator(modelFiles, solvers, relTols, absTols, iterations, readPath, writefile)
    for modelFile in modelFiles
        
        println(modelFile)
        for solver in solvers
            println(solver)
            for relTol in relTols
                println(relTol)
                for absTol in absTols
                    println(absTol)
                    modelSolver(modelFile, solver, relTol, absTol, iterations, readPath, writefile)
                end
            end
        end
    end
end


function main()
    readPath = pwd() * "\\Pipeline_SBMLImporter\\JuliaModels"
    writePath = pwd() * "\\Pipeline_ModelSolver\\IntermediaryResults"
    fixDirectories(writePath)
    writefile = writePath * "\\benchmark_" * string(getNumberOfFiles(writePath) + 1) * ".csv"

    modelFiles = getModelFiles(readPath)
    solvers = getSolvers()
    relTols, absTols = getTolerances()
    iterations = 5
    
    modelSolverIterator(modelFiles, solvers, relTols, absTols, iterations, readPath, writefile)

    #modelSolver(modelFiles[1], solvers[1], relTols[1], absTols[1], iterations, readPath, writefile)
    

end

main()