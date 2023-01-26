using ModelingToolkit 
using DifferentialEquations
using DataFrames
using CSV 
using ForwardDiff
using ReverseDiff
using StatsBase
using Random
using LinearAlgebra
using Distributions
using Printf
using Zygote
using SciMLSensitivity
using Sundials


# Relevant PeTab structs for compuations 
include(joinpath(pwd(), "src", "PeTab_structs.jl"))

# PeTab importer to get cost, grad etc 
include(joinpath(pwd(), "src", "Process_PEtab_files", "Create_cost_grad_hessian.jl"))

# HyperCube sampling 
include(joinpath(pwd(), "src", "Optimizers", "Lathin_hypercube.jl"))

# For converting to SBML 
include(joinpath(pwd(), "src", "SBML", "SBML_to_ModellingToolkit.jl"))


"""
    testCostGradHess(solver, tol; printRes::Bool=false)

    Compare cost and gradient for Julia PeTab importer against the cost and 
    gradient computed in PyPesto for Boehm model. 
"""
function compareAgainstPyPesto(peTabModel::PeTabModel, solver, tol; printRes::Bool=false)

    peTabOpt = setUpCostGradHess(peTabModel, solver, tol, sensealg = ForwardDiffSensitivity(), 
                                 sensealgForward = ForwardSensitivity(), solverForward=CVODE_BDF(),
                                 adjSolver=solver, adjTol=tol, 
                                 adjSensealg=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(true)))
    peTabOptAlt = setUpCostGradHess(peTabModel, solver, tol, 
                                    sensealgForward = :AutoDiffForward, solverForward=solver)                                 

    paramVals = CSV.read(pwd() * "/tests/Boehm/Params.csv", DataFrame)
    paramMat = paramVals[!, Not([:Id, :ratio, :specC17])]
    nParam = ncol(paramMat)

    costPython = (CSV.read(pwd() * "/tests/Boehm/Cost.csv", DataFrame))[!, :Cost]
    gradPythonMat = CSV.read(pwd() * "/tests/Boehm/Grad.csv", DataFrame)
    gradPythonMat = gradPythonMat[!, Not([:Id, :ratio, :specC17])]
    hessPythonMat = CSV.read(pwd() * "/tests/Boehm/Hess.csv", DataFrame)
    hessPythonMat = hessPythonMat[!, Not([:Id, :Epo_degradation_BaF3ratio, :k_exp_heteroratio, :k_exp_homoratio, :k_imp_heteroratio, :k_imp_homoratio, :k_phosratio, :ratioEpo_degradation_BaF3, :ratiok_exp_hetero, :ratiok_exp_homo, :ratiok_imp_hetero, :ratiok_imp_homo, :ratiok_phos, :ratioratio, :ratiosd_pSTAT5A_rel, :ratiosd_pSTAT5B_rel, :ratiosd_rSTAT5A_rel, :ratiospecC17, :sd_pSTAT5A_relratio, :sd_pSTAT5B_relratio, :sd_rSTAT5A_relratio, :specC17ratio, :Epo_degradation_BaF3specC17, :k_exp_heterospecC17, :k_exp_homospecC17, :k_imp_heterospecC17, :k_imp_homospecC17, :k_phosspecC17, :sd_pSTAT5A_relspecC17, :sd_pSTAT5B_relspecC17, :sd_rSTAT5A_relspecC17, :specC17Epo_degradation_BaF3, :specC17k_exp_hetero, :specC17k_exp_homo, :specC17k_imp_hetero, :specC17k_imp_homo, :specC17k_phos, :specC17sd_pSTAT5A_rel, :specC17sd_pSTAT5B_rel, :specC17sd_rSTAT5A_rel, :specC17specC17])]

    for i in 1:nrow(paramMat)
        paramVec = collect(paramMat[i, :])

        costJulia = peTabOpt.evalF(paramVec)
        sqDiffCost = (costJulia - costPython[i])^2
        if sqDiffCost > 1e-5
            @printf("sqDiffCost = %.3e\n", sqDiffCost)
            @printf("Does not pass test on cost\n")
            return false
        end

        gradJulia = zeros(nParam); peTabOpt.evalGradF(gradJulia, paramVec)
        gradPython = collect(gradPythonMat[i, :])
        sqDiffGrad = sum((gradJulia - gradPython).^2)
        if sqDiffGrad > 1e-5
            @printf("sqDiffGrad = %.3e\n", sqDiffGrad)
            @printf("Does not pass test on gradient\n")
            return false
        end

        costZygote = peTabOpt.evalFZygote(paramVec)
        sqDiffZygote = (costZygote - costPython[i])^2
        if sqDiffZygote > 1e-5
            @printf("sqDiffCostZygote = %.3e\n", sqDiffZygote)
            @printf("Does not pass test on cost for Zygote\n")
            return false
        end

        gradZygote = zeros(nParam)
        peTabOpt.evalGradFZygote(gradZygote, paramVec)
        gradPython = collect(gradPythonMat[i, :])
        sqDiffGradZygote = sum((gradZygote - gradPython).^2)
        if sqDiffGradZygote > 1e-5
            @printf("sqDiffGradZygote = %.3e\n", sqDiffGradZygote)
            @printf("Does not pass test on gradient from Zygote\n")
            return false
        end

        # Evaluate lower level adjoint sensitivity interfance gradient 
        gradAdj = zeros(nParam)
        peTabOpt.evalGradFAdjoint(gradAdj, paramVec)
        sqDiffGradAdjoint = sum((gradAdj - gradPython).^2)
        if sqDiffGradAdjoint > 1e-4
            @printf("sqDiffGradAdjoint = %.3e\n", sqDiffGradAdjoint)
            @printf("Does not pass test on adjoint gradient gradient\n")
            return false
        end

        # Forward sensitivity equations using autodiff 
        gradForwardEqAuto = zeros(nParam)
        peTabOptAlt.evalGradFForwardEq(gradForwardEqAuto, paramVec)
        sqDiffGradForwardEqAuto = sum((gradForwardEqAuto - gradPython).^2)
        if sqDiffGradForwardEqAuto > 1e-5
            @printf("sqDiffGradForwardEqAuto = %.3e\n", sqDiffGradForwardEqAuto)
            @printf("Does not pass test on gradient from Forward sensitivity equations\n")
            return false
        end

        # Forward sensitivity equations gradient 
        gradForwardEq = zeros(nParam)
        peTabOpt.evalGradFForwardEq(gradForwardEq, paramVec)
        gradPython = collect(gradPythonMat[i, :])
        sqDiffGradForwardEq = sum((gradForwardEq - gradPython).^2)
        if sqDiffGradForwardEq > 1e-5
            @printf("sqDiffGradForwardEq = %.3e\n", sqDiffGradForwardEq)
            @printf("Does not pass test on gradient from Forward sensitivity equations\n")
            return false
        end

        hessJulia = zeros(nParam,nParam); peTabOptAlt.evalHessGaussNewton(hessJulia, paramVec)
        hessJulia = hessJulia[:]
        hessPython = collect(hessPythonMat[i, :])
        sqDiffHess = sum((hessJulia - hessPython).^2)
        if sqDiffHess > 1e-5
            @printf("sqDiffHess = %.3e\n", sqDiffHess)
            @printf("Does not pass test on hessian\n")
            return false
        end

        if printRes == true
            @printf("sqDiffCost = %.3e\n", sqDiffCost)
            @printf("sqDiffGrad = %.3e\n", sqDiffGrad)
            @printf("sqDiffHess = %.3e\n", sqDiffHess)
            @printf("sqDiffCostZygote = %.3e\n", sqDiffZygote)
            @printf("sqDiffGradZygote = %.3e\n", sqDiffGradZygote)
            @printf("sqDiffGradAdjoint = %.3e\n", sqDiffGradAdjoint)
            @printf("sqDiffGradForwardEq = %.3e\n", sqDiffGradForwardEq)
            @printf("sqDiffGradForwardEqAuto = %.3e\n", sqDiffGradForwardEqAuto)
        end
    end

    return true 
end


peTabModel = setUpPeTabModel("Boehm_JProteomeRes2014", pwd() * "/tests/Boehm/", forceBuildJlFile=true)
passTest = compareAgainstPyPesto(peTabModel, Rodas4P(), 1e-12, printRes=true)
if passTest == true
    @printf("Passed test Boehm against PyPesto\n")
else
    @printf("Did not pass test Boehm against PyPesto\n")
end