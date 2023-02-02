
function getPathYmlFile(pathDir::String)::String
    iFile = findfirst(x -> splitext(x)[end] == ".yaml", readdir(pathDir))
    return joinpath(pathDir, readdir(pathDir)[iFile])
end


function getRandomModelParameters(petabModel::PEtabModel, 
                                  odeSolver::SciMLAlgorithm, 
                                  iCube;
                                  nParamCube=100,
                                  tol::Float64=1e-8)
  
    # Generate Cube if lacking and choose one parameter vector 
    petabProblem = setUpPEtabODEProblem(petabModel, odeSolver, solverAbsTol=tol, solverRelTol=tol)
    pathCube = petabModel.dirJulia * "Cube_ode_solvers.csv"
    createCube(pathCube, petabProblem, nParamCube, seed=123, verbose=false)
    cube = Matrix(CSV.read(pathCube, DataFrame))
    θ_est = cube[iCube, :]

    # Transform model parameters to correct scale for ODE-problem
    experimentalConditions, measurementsData, parametersData, observablesData = readPEtabFiles(petabModel)
    parameterInfo = processParameters(parametersData) 
    measurementInfo = processMeasurements(measurementsData, observablesData) 
    θ_indices = computeIndicesθ(parameterInfo, measurementInfo, petabModel.odeSystem, experimentalConditions)
    transformθ!(θ_est, θ_indices.θ_estNames)

    return θ_est[θ_indices.iθ_dynamic]
end


function getNominalODEValues(petabModel::PEtabModel)

    # Transform model parameters to correct scale for ODE-problem
    experimentalConditions, measurementsData, parametersData, observablesData = readPEtabFiles(petabModel)
    parameterInfo = processParameters(parametersData) 
    measurementInfo = processMeasurements(measurementsData, observablesData) 
    θ_indices = computeIndicesθ(parameterInfo, measurementInfo, petabModel.odeSystem, experimentalConditions)

    θ_estNames = θ_indices.θ_estNames
    paramVecNominal = [parameterInfo.nominalValue[findfirst(x -> x == θ_estNames[i], parameterInfo.parameterId)] for i in eachindex(θ_estNames)]

    return paramVecNominal[θ_indices.θ_dynamicNames]
end