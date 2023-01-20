# PEtab structs 
include(joinpath(pwd(), "src", "PeTab_structs.jl"))

# Files related to computing the cost (likelihood)
include(joinpath(pwd(), "src", "Compute_cost", "Compute_cost.jl"))
include(joinpath(pwd(), "src", "Compute_cost", "Compute_cost_zygote.jl"))

# Files related to computing derivatives 
include(joinpath(pwd(), "src", "Derivatives", "Hessian.jl"))
include(joinpath(pwd(), "src", "Derivatives", "Gradient.jl"))
include(joinpath(pwd(), "src", "Derivatives", "Adjoint_sensitivity_analysis.jl"))
include(joinpath(pwd(), "src", "Derivatives", "Forward_sensitivity_equations.jl"))
include(joinpath(pwd(), "src", "Derivatives", "Gauss_newton.jl"))
include(joinpath(pwd(), "src", "Derivatives", "Common.jl"))

# Files related to solving the ODE-system 
include(joinpath(pwd(), "src", "Solve_ODE", "Solve_ode_model.jl"))

# Files related to distributed computing 
include(joinpath(pwd(), "src", "Distributed", "Distributed.jl"))

# Files related to processing PEtab files 
include(joinpath(pwd(), "src", "Process_PEtab_files", "Common.jl"))
include(joinpath(pwd(), "src", "Process_PEtab_files", "Map_parameters.jl"))
include(joinpath(pwd(), "src", "Process_PEtab_files", "Create_obs_u0_sd_common.jl"))
include(joinpath(pwd(), "src", "Process_PEtab_files", "Create_obs_u0_sd_functions.jl"))
include(joinpath(pwd(), "src", "Process_PEtab_files", "Create_obs_u0_sd_derivatives.jl"))
include(joinpath(pwd(), "src", "Process_PEtab_files", "Process_PeTab_files.jl"))
include(joinpath(pwd(), "src", "Common.jl"))


"""
    setUpCostGradHess(peTabModel::PeTabModel, solver, tol::Float64)

    For a PeTab-model set up functions for computing i) the likelihood, ii) likelhood gradient, 
    and iii) likelhood Hessian block approximation. The functions are stored in PeTabOpt-struct 
    that can be used as input to the optimizers. 

    Currently the gradient for dynamic parameters (part of ODE-system) is computed via ForwardDiff, 
    and ReverseDiff is used for observable and sd parameters. The hessian approximation assumes the 
    interaction betweeen dynamic and (observable, sd) parameters is zero.
"""
function setUpCostGradHess(peTabModel::PeTabModel, 
                           solver::SciMLAlgorithm, 
                           tol::Float64; 
                           sensealg=ForwardDiffSensitivity(),
                           sensealgForward::Union{Symbol, SciMLSensitivity.AbstractForwardSensitivityAlgorithm}=ForwardSensitivity(),
                           sparseJac::Bool=false, 
                           absTolSS::Float64=1e-8, 
                           relTolSS::Float64=1e-6, 
                           nProcs::Signed=1,
                           adjSolver=Rodas5P(), 
                           solverForward=Rodas5P(autodiff=false),
                           adjTol::Float64=1e-6, 
                           adjSensealg::SciMLSensitivity.AbstractAdjointSensitivityAlgorithm=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(false)), 
                           adjSensealgSS::SciMLSensitivity.AbstractAdjointSensitivityAlgorithm=SteadyStateAdjoint())::PeTabOpt

    if !(typeof(adjSensealgSS) <: SteadyStateAdjoint)
        println("If you are using adjoint sensitivity analysis for a model with PreEq-criteria the most the most efficient adjSensealgSS is usually SteadyStateAdjoint. The algorithm you have provided, ", adjSensealgSS, "might not work (as there are some bugs here). In case it does not work, and SteadyStateAdjoint fails (because a dependancy on time) a good choice might be QuadratureAdjoint(autodiff=false, autojacvec=false)")
    end

    # Process PeTab files into type-stable Julia structs 
    experimentalConditionsFile, measurementDataFile, parameterDataFile, observablesDataFile = readDataFiles(peTabModel.dirModel, readObs=true)
    parameterData = processParameterData(parameterDataFile)
    measurementData = processMeasurementData(measurementDataFile, observablesDataFile) 
    simulationInfo = getSimulationInfo(peTabModel, measurementDataFile, measurementData, sensealg=sensealg, absTolSS=absTolSS, relTolSS=relTolSS)

    # Indices for mapping parameter-estimation vector to dynamic, observable and sd parameters correctly when calculating cost
    paramEstIndices = getIndicesParam(parameterData, measurementData, peTabModel.odeSystem, experimentalConditionsFile)
    
    # Set up potential prior for the parameters to estimate 
    priorInfo::PriorInfo = getPriorInfo(paramEstIndices, parameterDataFile)

    # Set model parameter values to those in the PeTab parameter data ensuring correct value of constant parameters 
    setParamToFileValues!(peTabModel.paramMap, peTabModel.stateMap, parameterData)

    # The time-span 5e3 is overwritten when performing actual forward simulations 
    odeProb = ODEProblem(peTabModel.odeSystem, peTabModel.stateMap, (0.0, 5e3), peTabModel.paramMap, jac=true, sparse=sparseJac)
    odeProb = remake(odeProb, p = convert.(Float64, odeProb.p), u0 = convert.(Float64, odeProb.u0))
    if sensealgForward == :AutoDiffForward
        odeProbSenseEq = deepcopy(odeProb)
    else
        odeProbSenseEq = ODEForwardSensitivityProblem(odeProb.f, odeProb.u0, odeProb.tspan, odeProb.p, 
                                                      sensealg=sensealgForward)
    end

    # Functions to map experimental conditions and parameters correctly to the ODE model 
    changeToExperimentalCondUse! = (pVec, u0Vec, expID, dynParamEst) -> changeExperimentalCondEst!(pVec, u0Vec, expID, dynParamEst, peTabModel, paramEstIndices)
    changeToExperimentalCondSenseEqUse! = (pVec, u0Vec, expID, dynParamEst) -> changeExperimentalCondEstSenseEq!(pVec, u0Vec, expID, dynParamEst, peTabModel, paramEstIndices)
    changeToExperimentalCondUse = (pVec, u0Vec, expID, dynParamEst) -> changeExperimentalCondEst(pVec, u0Vec, expID, dynParamEst, peTabModel, paramEstIndices)
    changeModelParamUse! = (pVec, u0Vec, paramEst) -> changeODEProblemParameters!(pVec, u0Vec, paramEst, paramEstIndices, peTabModel)
    changeModelParamUse = (pVec, paramEst) -> changeODEProblemParameters(pVec, paramEst, paramEstIndices, peTabModel)

    # Set up function which solves the ODE model for all conditions and stores result 
    solveOdeModelAllCondUse! = (solArrayArg, odeProbArg, dynParamEst, expIDSolveArg) -> solveOdeModelAllExperimentalCond!(solArrayArg, odeProbArg, dynParamEst, changeToExperimentalCondUse!, simulationInfo, solver, tol, tol, peTabModel.getTStops, onlySaveAtTobs=true, expIDSolve=expIDSolveArg)
    solveOdeModelAllCondAdjUse! = (solArrayArg, odeProbArg, dynParamEst, expIDSolveArg) -> solveOdeModelAllExperimentalCond!(solArrayArg, odeProbArg, dynParamEst, changeToExperimentalCondUse!, simulationInfo, solver, tol, tol, peTabModel.getTStops, denseSol=true, expIDSolve=expIDSolveArg, trackCallback=true)
    solveOdeModelAtCondZygoteUse = (odeProbArg, conditionId, dynParamEst, t_max) -> solveOdeModelAtExperimentalCondZygote(odeProbArg, conditionId, dynParamEst, t_max, changeToExperimentalCondUse, measurementData, simulationInfo, solver, tol, tol, sensealg, peTabModel.getTStops)
    if sensealgForward == :AutoDiffForward
        solveOdeModelAllCondForwardEq! = (solArrayArg, SMat, odeProbArg, dynParamEst, expIDSolveArg) -> solveOdeModelAllExperimentalCond!(solArrayArg, SMat, odeProbArg, dynParamEst, changeToExperimentalCondUse!, changeModelParamUse!, simulationInfo, solverForward, tol, tol, peTabModel.getTStops, onlySaveAtTobs=true, expIDSolve=expIDSolveArg)                                           
    else
        solveOdeModelAllCondForwardEq! = (solArrayArg, odeProbArg, dynParamEst, expIDSolveArg) -> solveOdeModelAllExperimentalCond!(solArrayArg, odeProbArg, dynParamEst, changeToExperimentalCondSenseEqUse!, simulationInfo, solverForward, tol, tol, peTabModel.getTStops, onlySaveAtTobs=true, expIDSolve=expIDSolveArg)
    end
    solveOdeModelAllCondGuassNewtonForwardEq! = (solArrayArg, SMat, odeProbArg, dynParamEst, expIDSolveArg) -> solveOdeModelAllExperimentalCond!(solArrayArg, SMat, odeProbArg, dynParamEst, changeToExperimentalCondUse!, changeModelParamUse!, simulationInfo, solverForward, tol, tol, peTabModel.getTStops, onlySaveAtTobs=true, expIDSolve=expIDSolveArg)                                           

    if nProcs > 1 && nprocs() != nProcs
        println("Error : PEtab importer was set to build the cost, grad and hessian with $nProcs processes, 
                 however, Julia is currently running with ", nprocs(), " processes which does not match input 
                 value. Input argument nProcs must match nprocs()")
    elseif nProcs == 1
        evalF = (paramVecEst) -> computeCost(paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, priorInfo)
        evalGradF = (grad, paramVecEst) -> computeGradientAutoDiff!(grad, paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, priorInfo)    
        evalGradFAdjoint = (grad, paramVecEst) -> computeGradientAdjointEquations!(grad, paramVecEst, adjSolver, adjSensealg, adjSensealgSS, adjTol, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondAdjUse!, priorInfo) 
        evalGradFForwardEq = (grad, paramVecEst) -> computeGradientForwardEquations!(grad, paramVecEst, peTabModel, odeProbSenseEq, sensealgForward, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondForwardEq!, priorInfo) 
        evalHessApprox = (hessianMat, paramVecEst) -> computeHessianBlockApproximation!(hessianMat, paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, priorInfo)
        evalHess = (hessianMat, paramVecEst) -> computeHessian(hessianMat, paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondUse!, priorInfo)
    elseif nProcs > 1 && nprocs() == nProcs
        evalF, evalGradF, evalGradFForwardEq, evalGradFAdjoint, evalHess, evalHessApprox = setUpPEtabOptDistributed(peTabModel, solver, tol, 
                                                                                                                    adjSolver, adjSensealg, adjSensealgSS, adjTol,
                                                                                                                    solverForward, sensealgForward, 
                                                                                                                    parameterData, measurementData, 
                                                                                                                    simulationInfo, paramEstIndices, priorInfo, odeProb)
    end
    evalFZygote = (paramVecEst) -> computeCostZygote(paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse, solveOdeModelAtCondZygoteUse, priorInfo)
    evalGradFZygote = (grad, paramVecEst) -> computeGradientZygote(grad, paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse, solveOdeModelAtCondZygoteUse, priorInfo)
    evalHessGaussNewton = (hessian, paramVecEst) -> computeGaussNewtonHessianApproximation!(hessian, paramVecEst, odeProb, peTabModel, simulationInfo, paramEstIndices, measurementData, parameterData, changeModelParamUse!, solveOdeModelAllCondGuassNewtonForwardEq!, priorInfo)      

    # Lower and upper bounds for parameters to estimate 
    namesParamEst = paramEstIndices.namesParamEst
    lowerBounds = [parameterData.lowerBounds[findfirst(x -> x == namesParamEst[i], parameterData.parameterID)] for i in eachindex(namesParamEst)] 
    upperBounds = [parameterData.upperBounds[findfirst(x -> x == namesParamEst[i], parameterData.parameterID)] for i in eachindex(namesParamEst)] 
    # Parameter with nominal values in PeTab file 
    paramVecNominal = [parameterData.paramVal[findfirst(x -> x == namesParamEst[i], parameterData.parameterID)] for i in eachindex(namesParamEst)]

    # Transform upper and lower bounds if the case 
    transformθ!(lowerBounds, namesParamEst, parameterData, reverseTransform=true)
    transformθ!(upperBounds, namesParamEst, parameterData, reverseTransform=true)
    paramVecNominalTransformed = transformθ(paramVecNominal, namesParamEst, parameterData, reverseTransform=true)

    peTabOpt = PeTabOpt(evalF, 
                        evalFZygote,
                        evalGradF, 
                        evalGradFZygote,
                        evalGradFAdjoint,
                        evalGradFForwardEq,
                        evalHess,
                        evalHessApprox, 
                        evalHessGaussNewton,
                        Int64(length(namesParamEst)), 
                        namesParamEst, 
                        paramVecNominal, 
                        paramVecNominalTransformed, 
                        lowerBounds, 
                        upperBounds, 
                        peTabModel.dirModel * "Cube" * peTabModel.modelName * ".csv",
                        peTabModel)
    return peTabOpt
end

# Evaluate contribution of potential prior to the final cost function value. Not mutating so works with both Zygote 
# and forwardiff.
function evalPriors(paramVecTransformed, 
                    paramVecNotTransformed,
                    namesParamVec::Array{String, 1}, 
                    paramEstIndices::ParameterIndices, 
                    priorInfo::PriorInfo)::Real

    if priorInfo.hasPriors == false
        return 0.0
    end

    k = 0
    priorContribution = 0.0
    for i in eachindex(paramVecNotTransformed)
        iParam = findfirst(x -> x == namesParamVec[i], paramEstIndices.namesParamEst)
        if priorInfo.priorOnParamScale[iParam] == true
            pInput = paramVecNotTransformed[i]
        else
            pInput = paramVecTransformed[i]
        end
        priorContribution += priorInfo.logpdf[iParam](pInput)

        if priorInfo.logpdf[iParam](pInput) != 0
            k += 1
        end
    end

    return priorContribution
end



