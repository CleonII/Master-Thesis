# Functions for solving a PeTab ODE-system for all experimental conditions. Can handle models with and without 
# a preequilibration-criteria (where the model must be simulated to a steady state).
# TODO: Pre-compute better with indices and when extracting the t-max value.


"""
    solveOdeModelAtFileValues(peTabModel::PeTabModel)

    For a peTab model solve the ODE:s at the values in the parameters-PeTab-file 
    using a specific solver with absTol=relTol=tol. Returns an array with the 
    solution for each experimenta condition along with a vector with the condition 
    name for each solution 
"""
function solveOdeModelAtFileValues(peTabModel::PeTabModel, solver, tol::Float64; 
                                   nTSave::Int64=0, denseSol::Bool=true, absTolSS=1e-8, relTolSS=1e-6)

    # Process PeTab files into type-stable Julia structs 
    experimentalConditionsFile, measurementDataFile, parameterDataFile, observablesDataFile = readDataFiles(peTabModel.dirModel, readObs=true)
    parameterData = processParameterData(parameterDataFile)
    simulationInfo = getSimulationInfo(measurementDataFile, absTolSS=absTolSS, relTolSS=relTolSS)
    
    # Set model parameter values to those in the PeTab parameter data ensuring correct value of constant parameters 
    setParamToFileValues!(peTabModel.paramMap, peTabModel.stateMap, parameterData)

    # The time-span 5e3 is overwritten when performing actual forward simulations 
    odeProb = ODEProblem(peTabModel.odeSystem, peTabModel.stateMap, (0.0, 5e3), peTabModel.paramMap, jac=true)
    odeProb = remake(odeProb, p = convert.(Float64, odeProb.p), u0 = convert.(Float64, odeProb.u0))

    # Functions to map experimental conditions and parameters correctly to the ODE model 
    changeToExperimentalCondUse! = (pVec, u0Vec, expID) -> changeExperimentalCond!(pVec, u0Vec, expID, parameterData, experimentalConditionsFile, peTabModel)

    # Set up function which solves the ODE model for all conditions and stores result 
    solArray, status = solveOdeModelAllExperimentalCond(odeProb, changeToExperimentalCondUse!, measurementDataFile, simulationInfo, solver, tol, nTSave=nTSave, denseSol=denseSol)

    return solArray, simulationInfo
end


"""
    solveOdeModelAllExperimentalCond!(solArray::Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1},
                                      prob::ODEProblem, 
                                      changeToExperimentalCondUse!::Function, 
                                      measurementData::DataFrame,
                                      simulationInfo::SimulationInfo,
                                      solver,
                                      tol::Float64;
                                      nTSave::Int64=0, 
                                      denseSol::Bool=true)::Bool

    Solve a PeTab ODE model for all experimental conditions specified in the PeTab experimentaCondition-file 
    and store each ODE-soluation into solArray. Returns true if model could be solved succesfully for all 
    conditions, else returns false.
    
    All inputs are automatically computed by the PeTab importer. Each experimental condition is simulated to the 
    maximum time for that condition found in the the measurementData PeTab-file. In case a pre-equlibration condition 
    exists the model is first simulated to a steady-state. Then starting from the steady state the solution stored in 
    solArray is calculcated.

    # Args
    `solArray`: array storing ODE-solution for each experimental condition. Is pre-allocated by PeTab importer.
    `prob`: ODEProblem struct for the model to be simulated.
    `changeToExperimentalCondUse!`: function that changes the parameters for the ODE-problem to those for a specific experimental condition 
    `measurementData`: the PeTab measurementData file 
    `simulationInfo`: struct storing simulation info like experimental conditions, whether or not to simulate to a steady state.
    `solver`: ode-solver for solving the ode-problem. Any Julia solver or alg-hint works.
    `tol`: rel- and abs-tol for the ODE solver 
    `nTSave`: number of equidistant data-points to save for each ODE-solution. If 0 the ODE-solver outputs a dense solution.
    `denseSol`: whether to have a dense (true) or none-dense (false) ODE-solution for each experimental condition.
    
    See also: [`setUpCostFunc`, `SimulationInfo`, `changeToExperimentalCond!`]
"""
function solveOdeModelAllExperimentalCond!(solArray::Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1},
                                           prob::ODEProblem, 
                                           changeToExperimentalCondUse!::Function, 
                                           measurementData::DataFrame,
                                           simulationInfo::SimulationInfo,
                                           solver,
                                           tol::Float64;
                                           nTSave::Int64=0, 
                                           denseSol::Bool=true)::Bool

    local sucess = true
    # In case the model is first simulated to a steady state 
    if simulationInfo.simulateSS == true
        k = 1
        @inbounds for i in eachindex(simulationInfo.firstExpIds)
            for j in eachindex(simulationInfo.shiftExpIds[i])

                firstExpId = simulationInfo.firstExpIds[i]
                shiftExpId = simulationInfo.shiftExpIds[i][j]
                t_max_ss = getTimeMax(measurementData, shiftExpId)
                solArray[k] = solveOdeSS(prob, 
                                         changeToExperimentalCondUse!, 
                                         firstExpId, 
                                         shiftExpId, 
                                         tol, 
                                         t_max_ss, 
                                         solver, 
                                         nTSave=nTSave, 
                                         denseSol=denseSol, 
                                         absTolSS=simulationInfo.absTolSS, 
                                         relTolSS=simulationInfo.relTolSS)
                
                # Keep index of which forward solution index k corresponds for calculating cost 
                simulationInfo.conditionIdSol[k] = firstExpId * shiftExpId

                if solArray[k].retcode != :Success
                    sucess = false
                    break 
                end
                k += 1
            end
        end

    # In case the model is not first simulated to a steady state 
    elseif simulationInfo.simulateSS == false

        @inbounds for i in eachindex(simulationInfo.firstExpIds)
            
            firstExpId = simulationInfo.firstExpIds[i]
            t_max = getTimeMax(measurementData, firstExpId)
            solArray[i] = solveOdeNoSS(prob, 
                                       changeToExperimentalCondUse!, 
                                       firstExpId, 
                                       tol, 
                                       solver, 
                                       t_max, 
                                       nTSave=nTSave, 
                                       denseSol=denseSol, 
                                       absTolSS=simulationInfo.absTolSS, 
                                       relTolSS=simulationInfo.relTolSS)

            # Keep index of which forward solution index i corresponds for calculating cost 
            simulationInfo.conditionIdSol[i] = firstExpId

            if !(solArray[i].retcode == :Success || solArray[i].retcode == :Terminated)
                sucess = false
                break 
            end
        end
    end

    return sucess
end


"""
    solveOdeModelAllExperimentalCond(prob::ODEProblem, 
                                     changeToExperimentalCondUse!::Function, 
                                     measurementData::DataFrame,
                                     simulationInfo::SimulationInfo,
                                     solver,
                                     tol::Float64;
                                     nTSave::Int64=0, 
                                     denseSol::Bool=true)

    Solve a PeTab ODE model for all experimental conditions specified in the PeTab experimentaCondition-file 
    and return each ODE-soluation into solArray. Arguments the same as for `solveOdeModelAllExperimentalCond!`
     
    See also: [`solveOdeModelAllExperimentalCond!`]
"""
function solveOdeModelAllExperimentalCond(prob::ODEProblem, 
                                          changeToExperimentalCondUse!::Function, 
                                          measurementData::DataFrame,
                                          simulationInfo::SimulationInfo,
                                          solver, 
                                          tol::Float64;
                                          nTSave::Int64=0, 
                                          denseSol::Bool=true)

    local solArray
    # Compute the number of ODE-solutions to cover all experimental conditions 
    if simulationInfo.simulateSS == true
        nShiftId = Int(sum([length(simulationInfo.shiftExpIds[i]) for i in eachindex(simulationInfo.shiftExpIds)]))
        solArray = Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1}(undef, nShiftId)
    else
        nExperimentalCond = Int64(length(simulationInfo.firstExpIds))
        solArray = Array{Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}, 1}(undef, nExperimentalCond)
    end
        
    success = solveOdeModelAllExperimentalCond!(solArray, 
                                                prob, 
                                                changeToExperimentalCondUse!, 
                                                measurementData, 
                                                simulationInfo, 
                                                solver, 
                                                tol, 
                                                nTSave=nTSave, 
                                                denseSol=denseSol)

    return solArray, success
end


"""
    solveOdeSS(prob::ODEProblem, 
               changeToExperimentalCondUse!::Function, 
               firstExpId::String, 
               shiftExpId::String,
               tol::Float64, 
               t_max_ss::Float64,
               solver;
               tSave=Float64[], 
               nTSave=0, 
               denseSol=true)::Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}

    For an experimentaCondition specifed by firstExpId (preequilibration ID) and shiftExpId
    (postequilibration ID) solve a PeTab ODE-problem using any Julia ODE-solver (solver can 
    also be an alg-hint) using absTol and relTol value given by tol. Returns an ODE-solution.

    Here the model is first simualted to a steady state using the experimental-condition parameters 
    corresponding to firstExpId in the experimentaCondition peTab-file. Then using the experimental-condition 
    parameters specifed by shiftExpId the actual simulation that is saved and returned is produced.
    Optional args are the save as for `solveOdeModelAllExperimentalCond!`, and this function is called 
    by solveOdeModelAllExperimentalCond!. 
     
    See also: [`solveOdeModelAllExperimentalCond!`]
"""
function solveOdeSS(prob::ODEProblem, 
                    changeToExperimentalCondUse!::Function, 
                    firstExpId::String, 
                    shiftExpId::String,
                    tol::Float64, 
                    t_max_ss::Float64,
                    solver;
                    tSave=Float64[], 
                    nTSave=0, 
                    denseSol::Bool=true,
                    absTolSS::Float64=1e-8, 
                    relTolSS::Float64=1e-6)::Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}

    # Sanity check input. Can only provide either nTsave (points to save solution at) or tSave (number of points to save)
    if length(tSave) != 0 && nTSave != 0
        println("Error : Can only provide tSave (vector to save at) or nTSave as saveat argument to solvers")
    elseif nTSave != 0
        saveAtVec = collect(LinRange(0.0, t_max_ss, nTSave))
    else
        saveAtVec = tSave
    end

    # Sanity check input. Both tSave and nTsave must be empty or zero in order to be able to output a dense solution.
    if (isempty(tSave) && nTSave == 0) && denseSol == true
        dense = true
    else
        dense = false
    end
    
    # Different funcion calls to solve are required if a solver or a Alg-hint are provided. 
    # The preequilibration simulations are terminated upon a steady state using the TerminateSteadyState callback.
    if typeof(solver) <: Vector{Symbol} # Alg-hint case
        solveCallPre = (prob) -> solve(prob, alg_hints=solver, abstol=tol, reltol=tol, dense=false, 
                                       callback=TerminateSteadyState(absTolSS, relTolSS))
        solveCallPost = (prob) -> solve(prob, alg_hints=solver, abstol=tol, reltol=tol, saveat=saveAtVec, dense=dense, isoutofdomain = (u,p,t)->any(x->x<0,u))

    else # Julia solver case
        solveCallPre = (prob) -> solve(prob, solver, abstol=tol, reltol=tol, dense=false, 
            callback=TerminateSteadyState(absTolSS, relTolSS))
        solveCallPost = (prob) -> solve(prob, solver, abstol=tol, reltol=tol, saveat=saveAtVec, dense=dense)
    end

    # Change to parameters for the preequilibration simulations 
    changeToExperimentalCondUse!(prob.p, prob.u0, firstExpId)
    u0_pre = deepcopy(prob.u0)
    prob = remake(prob, tspan = (0.0, 1e8), p = prob.p[:], u0 = prob.u0[:])

    # Terminate if a steady state was not reached in preequilibration simulations 
    sol_pre = solveCallPre(prob)
    if sol_pre.retcode != :Terminated
        return sol_pre
    end

    # Change to parameters for the post steady state parameters 
    changeToExperimentalCondUse!(prob.p, prob.u0, shiftExpId)
    # Sometimes the experimentaCondition-file changes the initial values for a state 
    # whose value was changed in the preequilibration-simulation. The experimentaCondition
    # value is prioritized by only changing u0 to the steady state value for those states  
    # that were not affected by change to shiftExpId.
    has_not_changed = (prob.u0 .== u0_pre)
    prob.u0[has_not_changed] .= sol_pre.u[end][has_not_changed]
    prob = remake(prob, tspan = (0.0, t_max_ss))
    
    sol = solveCallPost(prob)

    return sol
end


"""
    solveOdeNoSS(prob::ODEProblem, 
                 changeToExperimentalCondUse!::Function, 
                 firstExpId::String, 
                 tol::Float64, 
                 t_max::Float64,
                 solver;
                 tSave=Float64[], 
                 nTSave=0, 
                 denseSol=true)::Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}

    For an experimentaCondition specifed by firstExpId solve a PeTab ODE-problem using any Julia ODE-solver 
    (solver can also be an alg-hint) using absTol and relTol value given by tol. Returns an ODE-solution.

    Here the model is not simulated to steady state first first simualted to a steady state
    Optional args are the save as for `solveOdeModelAllExperimentalCond!`, and this function 
    is called by `solveOdeModelAllExperimentalCond!`. 
     
    See also: [`solveOdeModelAllExperimentalCond!`]
"""
function solveOdeNoSS(prob::ODEProblem, 
                      changeToExperimentalCondUse!::Function, 
                      firstExpId::String, 
                      tol::Float64, 
                      solver, 
                      t_max::Float64; 
                      tSave=Float64[], 
                      nTSave::Int64=0, 
                      denseSol::Bool=true, 
                      absTolSS::Float64=1e-8, 
                      relTolSS::Float64=1e-6)::Union{OrdinaryDiffEq.ODECompositeSolution, ODESolution}

    # Sanity check input. Can only provide either nTsave (points to save solution at) or tSave (number of points to save)
    if length(tSave) != 0 && nTSave != 0
        println("Error : Can only provide tSave (vector to save at) or nTSave as saveat argument to solvers")
    elseif nTSave != 0
        saveAtVec = collect(LinRange(0.0, t_max, nTSave))
    else
        saveAtVec = tSave
    end

    # Sanity check input. Both tSave and nTsave must be empty or zero in order to be able to output a dense solution.
    if (isempty(tSave) && nTSave == 0) && denseSol == true
        dense = true
    else
        dense = false
    end

    # Change to parameters to the relevant experimental condition 
    t_max_use = isinf(t_max) ? 1e6 : t_max
    changeToExperimentalCondUse!(prob.p, prob.u0, firstExpId)
    probUse = remake(prob, tspan=(0.0, t_max_use), u0 = prob.u0[:], p = prob.p[:])

    # Different funcion calls to solve are required if a solver or a Alg-hint are provided. 
    # If t_max = inf the model is simulated to steady state using the TerminateSteadyState callback.
    if typeof(solver) <: Vector{Symbol} && isinf(t_max)
        solveCall = (probArg) -> solve(probArg, alg_hints=solver, abstol=tol, reltol=tol, save_on=false, save_end=true, dense=dense, 
            callback=TerminateSteadyState(absTolSS, relTolSS))
    elseif typeof(solver) <: Vector{Symbol} && !isinf(t_max)
        solveCall = (probArg) -> solve(probArg, alg_hints=solver, abstol=tol, reltol=tol, saveat=saveAtVec, dense=dense)
    elseif !(typeof(solver) <: Vector{Symbol}) && isinf(t_max)
        solveCall = (probArg) -> solve(probArg, solver, abstol=tol, reltol=tol, save_on=false, save_end=true, dense=dense, 
            callback=TerminateSteadyState(absTolSS, relTolSS))
    elseif !(typeof(solver) <: Vector{Symbol}) && !isinf(t_max)
        solveCall = (probArg) -> solve(probArg, solver, abstol=tol, reltol=tol, saveat=saveAtVec, dense=dense)
    else
        println("Error : Solver option does not exist")        
    end

    sol = solveCall(probUse)

    return sol
end


"""
    changeExperimentalCond!(paramVec, 
                            stateVec, 
                            expID::String, 
                            parameterData::ParamData,
                            experimentalConditions::DataFrame,
                            peTabModel::PeTabModel)

    Change the ODE parameter vector (paramVec) and initial value vector (stateVec)
    values to the values specified for the experimental ID expID given by the 
    experimentalConditions peTab-file for a specific peTabModel. 

    parameterData is needed to correctly map the parameters. The function can 
    handle that paramVec is a Float64 vector or a vector of Duals for the 
    gradient calculations. This function is used by the `solveOdeModelAllExperimentalCond!`. 
     
    See also: [`solveOdeModelAllExperimentalCond!`]
"""
function changeExperimentalCond!(paramVec, 
                                 stateVec, 
                                 expID::String, 
                                 parameterData::ParamData,
                                 experimentalConditions::DataFrame,
                                 peTabModel::PeTabModel)

    # TODO : Several things can be precomputed for this function

    # When computing the gradient the paramMap must be able to handle dual 
    # numbers, hence creating a paramMapUse
    paramMapUse = convert.(Pair{Num, eltype(paramVec)}, peTabModel.paramMap)

    # Extract names of parameters to change for specific experimental condition 
    colNames = names(experimentalConditions)
    i_start = "conditionName" in colNames ? 3 : 2
    paramStateChange = colNames[i_start:end]
    if isempty(paramStateChange)
        return 
    end

    # As values to change to can be a parameter or value they storing them as string initally is required 
    valsChangeTo = string.(Vector(experimentalConditions[getRowExpId(expID, experimentalConditions), i_start:end]))
    
    # To help with mapping extract parameter names as string
    parameterNamesStr = string.([paramMapUse[i].first for i in eachindex(paramMapUse)])
    stateNamesStr = replace.(string.(peTabModel.stateNames), "(t)" => "")
    
    # Get number of states and parameters to change 
    nParamChange = length(intersect(paramStateChange, parameterNamesStr))
    nStateChange = length(intersect(paramStateChange, stateNamesStr))
        
    # Keep tab of which parameters are changed.
    iParamChange = Array{Int64, 1}(undef, nParamChange)
    iStateChange = Array{Int64, 1}(undef, nStateChange)
    valChangeU0 = Array{Float64, 1}(undef, nStateChange)
    iP, iS = 1, 1
    changeParam::Bool = true
    for i in eachindex(paramStateChange)
        
        variable = paramStateChange[i]
        # If param is a model parameter change said parameter. If param is one state according to PeTab 
        # standard the initial value for said state should be changed. 
        iChangeP = findfirst(x -> x == variable, string.(parameters(peTabModel.odeSystem))) # Do not change to map correctly to ODE-sys
        iChangeS = findfirst(x -> x == variable, stateNamesStr) # Can be precomputed but is not expansive
        if !isnothing(iChangeP)
            iParamChange[iP] = iChangeP
            changeParam = true
            iP += 1
        elseif !isnothing(iChangeS)
            iStateChange[iS] = iChangeS
            changeParam = false
        else
            println("Error : $variable cannot be mapped to experimental condition")
        end

        # Extract value param should be changed to 
        valChangeTo::Float64 = 0.0
        if isNumber(valsChangeTo[i])
            valChangeTo = parse(Float64, valsChangeTo[i])

        # In case the value to change to is given as parameter look for said value in parameterData struct 
        # (where all parameters are stored)
        elseif findfirst(x -> x == valsChangeTo[i], parameterData.parameterID) != nothing
            iVal = findfirst(x -> x == valsChangeTo[i], parameterData.parameterID)
            valChangeTo = parameterData.paramVal[iVal]

        else
            println("Error : Simulation parameter not found for experimental condition $expID")
            println("valsChangeTo[i] = ", valsChangeTo[i])
        end

        # Identify which index param corresponds to the in paramMap 
        if changeParam == true        
            iParam = findfirst(x -> x == variable, parameterNamesStr)
            if !isnothing(iParam)
                paramMapUse[iParam] = Pair(paramMapUse[iParam].first, valChangeTo) 
            else
                println("Error : Simulation parameter to change not found for experimental condition $expID")
            end

        # In case a state is changed 
        else
            valChangeU0[iS] = valChangeTo
            iS += 1
        end
    end

    # To prevent that all parameter values are reset to their defualt values  
    # only change the parameter values for the parameters change with the 
    # new experimental condition. 
    newVal = ModelingToolkit.varmap_to_vars(paramMapUse, peTabModel.paramNames)
    paramVec[iParamChange] .= newVal[iParamChange]
    peTabModel.evalU0!(stateVec, paramVec) 

    # In case an experimental condition maps directly to the initial value of a state. 
    if !isempty(iStateChange)
        stateVec[iStateChange] .= valChangeU0
    end

    return nothing
end


"""
    getRowExpId(expId::String, data::DataFrame; colSearch="conditionId") 

    Small helper function to get which row in a DataFrame corresponds to specific ExpId 
"""
function getRowExpId(expId::String, data::DataFrame; colSearch="conditionId") 
    return findfirst(x -> x == expId, data[!, colSearch])
end


"""
    getTimeMax(measurementData::DataFrame, expId::String)::Float64

    Small helper function to get the time-max value for a specific simulationConditionId when simulating 
    the PeTab ODE-model 
"""
function getTimeMax(measurementData::DataFrame, expId::String)::Float64
    return Float64(maximum(measurementData[findall(x -> x == expId, measurementData[!, "simulationConditionId"]), "time"]))
end