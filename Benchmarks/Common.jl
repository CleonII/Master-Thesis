
function getPathYmlFile(pathDir::String)::String
    iFile = findfirst(x -> splitext(x)[end] == ".yaml", readdir(pathDir))
    return joinpath(pathDir, readdir(pathDir)[iFile])
end


function getRandomModelParameters(petabModel::PEtabModel, 
                                  odeSolver::SciMLAlgorithm, 
                                  iCube;
                                  nParamCube=100,
                                  tol::Float64=1e-8,
                                  odeSolvers::Bool=true)
  
    # Generate Cube if lacking and choose one parameter vector 
    petabProblem = setUpPEtabODEProblem(petabModel, odeSolver, solverAbsTol=tol, solverRelTol=tol)
    if odeSolvers == true
        pathCube = joinpath(petabModel.dirJulia, "Cube_ode_solvers.csv")
    else
        pathCube = joinpath(petabModel.dirJulia, "Cost_grad_hess.csv")
    end
    createCube(pathCube, petabProblem, nParamCube, seed=123, verbose=false)
    cube = Matrix(CSV.read(pathCube, DataFrame))
    θ_est = cube[iCube, :]

    # Transform model parameters to correct scale for ODE-problem
    experimentalConditions, measurementsData, parametersData, observablesData = readPEtabFiles(petabModel)
    parameterInfo = processParameters(parametersData) 
    measurementInfo = processMeasurements(measurementsData, observablesData) 
    θ_indices = computeIndicesθ(parameterInfo, measurementInfo, petabModel)

    if odeSolvers == true
        transformθ!(θ_est, θ_indices.θ_estNames, θ_indices)
        return θ_est[θ_indices.iθ_dynamic]
    else
        return θ_est
    end
end


function getNominalθ(petabModel::PEtabModel)
  
    # Transform model parameters to correct scale for ODE-problem
    experimentalConditions, measurementsData, parametersData, observablesData = readPEtabFiles(petabModel)
    parameterInfo = processParameters(parametersData) 
    measurementInfo = processMeasurements(measurementsData, observablesData) 
    θ_indices = computeIndicesθ(parameterInfo, measurementInfo, petabModel)

    θ_estNames = θ_indices.θ_estNames               
    θ_nominal = [parameterInfo.nominalValue[findfirst(x -> x == θ_estNames[i], parameterInfo.parameterId)] for i in eachindex(θ_estNames)]
    transformθ!(θ_nominal, θ_indices.θ_estNames, θ_indices, reverseTransform=true)

    return θ_nominal
end


function getNominalODEValues(petabModel::PEtabModel)

    # Transform model parameters to correct scale for ODE-problem
    experimentalConditions, measurementsData, parametersData, observablesData = readPEtabFiles(petabModel)
    parameterInfo = processParameters(parametersData) 
    measurementInfo = processMeasurements(measurementsData, observablesData) 
    θ_indices = computeIndicesθ(parameterInfo, measurementInfo, petabModel)

    θ_estNames = θ_indices.θ_estNames
    paramVecNominal = [parameterInfo.nominalValue[findfirst(x -> x == θ_estNames[i], parameterInfo.parameterId)] for i in eachindex(θ_estNames)]

    return paramVecNominal[θ_indices.iθ_dynamic]
end



# For a PEtab model create a tempory directory where the PEtab files is changed such that nParamFixate are 
# fixed (not considered to be gradient related parameters )
function getPEtabModelNparamFixed(petabModel::PEtabModel, nParamFixate::Integer)::PEtabModel
    
    dirNew = joinpath(petabModel.dirModel, "Fewer_param" * petabModel.modelName)
    if !isdir(dirNew)
        mkdir(dirNew)
    end

    # Copy PEtab files to new directory 
    pathParameters = copyFileToDest(petabModel.pathParameters, dirNew)
    pathConditions = copyFileToDest(petabModel.pathConditions, dirNew)
    pathObservables = copyFileToDest(petabModel.pathObservables, dirNew)
    pathMeasurements = copyFileToDest(petabModel.pathMeasurements, dirNew)
    pathSBML = copyFileToDest(petabModel.pathSBML, dirNew)
    pathYAML = copyFileToDest(petabModel.pathYAML, dirNew)

    experimentalConditions, measurementsData, parametersData, observablesData = readPEtabFiles(pathYAML)

    parameterInfo = processParameters(parametersData) 
    measurementInfo = processMeasurements(measurementsData, observablesData) 
    θ_indices = computeIndicesθ(parameterInfo, measurementInfo, petabModel)

    parameterToFixate = string.(sample(θ_indices.θ_dynamicNames, nParamFixate, replace=false))
    for i in 1:nrow(parametersData)
        if parametersData[i, :parameterId] ∈ parameterToFixate
            parametersData[i, :estimate] = 0
        end
    end
    rm(pathParameters)
    CSV.write(pathParameters, parametersData, delim = '\t')

    return readPEtabModel(pathYAML)
end


# Get a PEtab model with permuted parmeters 
function getPEtabModelParamPermuted(petabModel::PEtabModel; seed=123)::PEtabModel
    
    Random.seed!(seed)
    dirNew = joinpath(petabModel.dirModel, "Permuted_param" * petabModel.modelName)
    if !isdir(dirNew)
        mkdir(dirNew)
    end

    # Copy PEtab files to new directory 
    pathParameters = copyFileToDest(petabModel.pathParameters, dirNew)
    pathConditions = copyFileToDest(petabModel.pathConditions, dirNew)
    pathObservables = copyFileToDest(petabModel.pathObservables, dirNew)
    pathMeasurements = copyFileToDest(petabModel.pathMeasurements, dirNew)
    pathSBML = copyFileToDest(petabModel.pathSBML, dirNew)
    pathYAML = copyFileToDest(petabModel.pathYAML, dirNew)

    experimentalConditions, measurementsData, parametersData, observablesData = readPEtabFiles(pathYAML)

    parameterInfo = processParameters(parametersData) 
    measurementInfo = processMeasurements(measurementsData, observablesData) 
    θ_indices = computeIndicesθ(parameterInfo, measurementInfo, petabModel)

    iEstimate = findall(x ->  x == 1, parametersData[!, :estimate])
    parametersData[iEstimate, :] = parametersData[shuffle(iEstimate), :]
    rm(pathParameters)
    CSV.write(pathParameters, parametersData, delim = '\t')

    return readPEtabModel(pathYAML)
end


# Will given a filePath copy the file (with the same name) to dirDest
function copyFileToDest(filePath::String, dirDest::String)
    newPath = joinpath(dirDest, splitdir(filePath)[end])
    cp(filePath, newPath, force=true)
    return newPath
end
