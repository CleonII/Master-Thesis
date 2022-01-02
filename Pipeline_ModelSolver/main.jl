using ModelingToolkit, OrdinaryDiffEq, BenchmarkTools, DataFrames, CSV

function getModelFiles(path)
    modelFiles = readdir(path)
    return modelFiles[1:2]
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
        mkdir(path)
    end
end

function modelSolver(modelFile, solver, relTol, absTol, iterations, readPath, writePath)
    modelPath = readPath * "\\" * modelFile
    include(modelPath)

    println("Here")
    new_sys = ode_order_lowering(sys)
    u0 = initialSpeciesValues

    p = trueParameterValues

    c = trueConstantsValues

    tspan = (0.0,150.0)
    prob = ODEProblem(new_sys,u0,tspan,[p;c],jac=true)

    runTime = Vector{Float64}(undef, iterations)

    for i in 1:iterations
        println("iteration: " * string(i))
        b = @benchmark solve($prob, $solver, reltol = $relTol, abstol = $absTol) 
        bm = mean(b)
        runTime[i] = bm.time 
    end
    println(runTime)
    
end

function modelSolverIterator(modelFiles, solvers, relTols, absTols, iterations, readPath, writePath)
    minRelTol = minimum(relTols)
    minAbsTol = minimum(absTols)

    for modelFile in modelFiles
        println(modelFile)
        for solver in solvers
            println(solver)
            modelSolver(modelFile, solver, minRelTol, minAbsTol, iterations, readPath, writePath)
            for relTol in relTols
                println(relTol)
                for absTol in absTols
                    println(absTol)
                    modelSolver(modelFile, solver, relTol, absTol, iterations, readPath, writePath)
                end
            end
        end
    end
end


function main()
    readPath = pwd() * "\\Pipeline_SBMLImporter\\JuliaModels"
    writePath = pwd() * "\\Pipeline_ModelSolver\\IntermediaryResults"
    fixDirectories(writePath)

    modelFiles = getModelFiles(readPath)
    println(modelFiles)
    solvers = getSolvers()
    println(solvers)
    relTols, absTols = getTolerances()
    println(relTols)
    println(absTols)
    iterations = 40
    
    #modelSolverIterator(modelFiles, solvers, relTols, absTols, iterations, readPath, writePath)

    modelSolver(modelFiles[1], solvers[1], relTols[1], absTols[1], iterations, readPath, writePath)
    

end

main()