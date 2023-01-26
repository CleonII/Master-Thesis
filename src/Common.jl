# Functions used by both the ODE-solvers and PeTab importer. 


"""
    setParamToFileValues!(paramMap, stateMap, paramData::ParamData)

    Function that sets the parameter and state values in paramMap and stateMap 
    to those in the PeTab parameters file. 

    Used when setting up the PeTab cost function, and when solving the ODE-system 
    for the values in the parameters-file. 
"""
function setParamToFileValues!(paramMap, stateMap, paramData::ParamData)

    parameterNames = paramData.parameterID
    parameterNamesStr = string.([paramMap[i].first for i in eachindex(paramMap)])
    stateNamesStr = replace.(string.([stateMap[i].first for i in eachindex(stateMap)]), "(t)" => "")
    for i in eachindex(parameterNames)
        
        parameterName = parameterNames[i]
        valChangeTo = paramData.paramVal[i]
        
        # Check for value to change to in parameter file 
        i_param = findfirst(x -> x == parameterName, parameterNamesStr)
        i_state = findfirst(x -> x == parameterName, stateNamesStr)

        if !isnothing(i_param)
            paramMap[i_param] = Pair(paramMap[i_param].first, valChangeTo) 
        elseif !isnothing(i_state)
            stateMap[i_state] = Pair(stateMap[i_state].first, valChangeTo)
        end
    end

end


function getRandomODEParameters(peTabModel::PeTabModel, 
                                solver::SciMLAlgorithm, 
                                iCube;
                                nParamCube=20,
                                costGradHess::Bool=false,
                                tol::Float64=1e-8)
  
    # Generate Cube if lacking and choose one parameter vector 
    peTabOpt = setUpCostGradHess(peTabModel, solver, tol)
    if costGradHess == true
        pathCube = peTabModel.dirModel * "Cube_cost_grad_hess.csv"
    else
        pathCube = peTabModel.dirModel * "Cube_ode_solve.csv"
    end
    createCube(pathCube, peTabOpt, nParamCube, seed=123, verbose=true)
    cube = Matrix(CSV.read(pathCube, DataFrame))
    paramEst = cube[iCube, :]

    # Change model parameters 
    experimentalConditionsFile, measurementDataFile, parameterDataFile, observablesDataFile = readDataFiles(peTabModel.dirModel, readObs=true)
    parameterData = processParameterData(parameterDataFile)
    measurementData = processMeasurementData(measurementDataFile, observablesDataFile) 
    paramEstIndices = getIndicesParam(parameterData, measurementData, peTabModel.odeSystem, experimentalConditionsFile)

    if costGradHess == false
        transformParamVec!(paramEst, paramEstIndices.namesParamEst, parameterData)
        return paramEst[paramEstIndices.iDynParam]
    else
        return paramEst
    end
end