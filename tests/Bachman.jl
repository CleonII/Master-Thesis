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
    gradient computed in PyPesto for Bachman model. 
"""
function compareAgainstPyPesto(peTabModel::PeTabModel, solver, tol; printRes::Bool=false)

    peTabOpt = setUpCostGradHess(peTabModel, solver, tol, sensealg = ForwardDiffSensitivity(), 
                                 sensealgForward = ForwardSensitivity(), solverForward=CVODE_BDF(),
                                 adjSolver=solver, adjTol=1e-10, adjSensealg=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(true)))
    peTabOptAlt = setUpCostGradHess(peTabModel, solver, tol, 
                                    sensealgForward = :AutoDiffForward, solverForward=solver)                                                                  

    # Parameter values to test gradient at 
    paramVals = CSV.read(pwd() * "/tests/Bachman/Params.csv", DataFrame)
    paramMat = paramVals[!, Not([:Id, :SOCS3RNAEqc, :CISRNAEqc])]

    # PyPesto gradient and cost values 
    costPython = (CSV.read(pwd() * "/tests/Bachman/Cost.csv", DataFrame))[!, :Cost]
    gradPythonMat = CSV.read(pwd() * "/tests/Bachman/Grad.csv", DataFrame)
    gradPythonMat = gradPythonMat[!, Not([:Id, :SOCS3RNAEqc, :CISRNAEqc])]
    hessPythonMat = CSV.read(pwd() * "/tests/Bachman/Hess.csv", DataFrame)
    hessPythonMat = hessPythonMat[!, Not([:Id, :CISEqcSOCS3RNAEqc, :CISEqcOESOCS3RNAEqc, :CISInhSOCS3RNAEqc, :CISRNADelaySOCS3RNAEqc, :CISRNAEqcSOCS3RNAEqc, :CISRNATurnSOCS3RNAEqc, :CISTurnSOCS3RNAEqc, :EpoRActJAK2SOCS3RNAEqc, :EpoRCISInhSOCS3RNAEqc, :EpoRCISRemoveSOCS3RNAEqc, :JAK2ActEpoSOCS3RNAEqc, :JAK2EpoRDeaSHP1SOCS3RNAEqc, :SHP1ActEpoRSOCS3RNAEqc, :SHP1DeaSOCS3RNAEqc, :SHP1ProOESOCS3RNAEqc, :SOCS3EqcSOCS3RNAEqc, :SOCS3EqcOESOCS3RNAEqc, :SOCS3InhSOCS3RNAEqc, :SOCS3RNADelaySOCS3RNAEqc, :SOCS3RNAEqcCISEqc, :SOCS3RNAEqcCISEqcOE, :SOCS3RNAEqcCISInh, :SOCS3RNAEqcCISRNADelay, :SOCS3RNAEqcCISRNAEqc, :SOCS3RNAEqcCISRNATurn, :SOCS3RNAEqcCISTurn, :SOCS3RNAEqcEpoRActJAK2, :SOCS3RNAEqcEpoRCISInh, :SOCS3RNAEqcEpoRCISRemove, :SOCS3RNAEqcJAK2ActEpo, :SOCS3RNAEqcJAK2EpoRDeaSHP1, :SOCS3RNAEqcSHP1ActEpoR, :SOCS3RNAEqcSHP1Dea, :SOCS3RNAEqcSHP1ProOE, :SOCS3RNAEqcSOCS3Eqc, :SOCS3RNAEqcSOCS3EqcOE, :SOCS3RNAEqcSOCS3Inh, :SOCS3RNAEqcSOCS3RNADelay, :SOCS3RNAEqcSOCS3RNAEqc, :SOCS3RNAEqcSOCS3RNATurn, :SOCS3RNAEqcSOCS3Turn, :SOCS3RNAEqcSTAT5ActEpoR, :SOCS3RNAEqcSTAT5ActJAK2, :SOCS3RNAEqcSTAT5Exp, :SOCS3RNAEqcSTAT5Imp, :SOCS3RNAEqcinit_EpoRJAK2, :SOCS3RNAEqcinit_SHP1, :SOCS3RNAEqcinit_STAT5, :SOCS3RNAEqcoffset_CIS_actd, :SOCS3RNAEqcoffset_CIS_cisoe, :SOCS3RNAEqcoffset_CIS_long, :SOCS3RNAEqcoffset_CIS_shp1oe, :SOCS3RNAEqcoffset_CIS_socs3oe, :SOCS3RNAEqcoffset_SOCS3_cisoe, :SOCS3RNAEqcoffset_SOCS3_long, :SOCS3RNAEqcoffset_SOCS3_socs3oe, :SOCS3RNAEqcoffset_pEpoR_actd, :SOCS3RNAEqcoffset_pEpoR_cisoe, :SOCS3RNAEqcoffset_pEpoR_cisoe_pepor, :SOCS3RNAEqcoffset_pEpoR_dr30, :SOCS3RNAEqcoffset_pEpoR_dr7, :SOCS3RNAEqcoffset_pEpoR_fine, :SOCS3RNAEqcoffset_pEpoR_long, :SOCS3RNAEqcoffset_pEpoR_shp1oe, :SOCS3RNAEqcoffset_pEpoR_socs3oe, :SOCS3RNAEqcoffset_pJAK2_actd, :SOCS3RNAEqcoffset_pJAK2_cisoe, :SOCS3RNAEqcoffset_pJAK2_dr30, :SOCS3RNAEqcoffset_pJAK2_dr7, :SOCS3RNAEqcoffset_pJAK2_fine, :SOCS3RNAEqcoffset_pJAK2_long, :SOCS3RNAEqcoffset_pJAK2_shp1oe, :SOCS3RNAEqcoffset_pJAK2_socs3oe, :SOCS3RNAEqcoffset_pSTAT5_actd, :SOCS3RNAEqcoffset_pSTAT5_cisoe, :SOCS3RNAEqcoffset_pSTAT5_conc, :SOCS3RNAEqcoffset_pSTAT5_long, :SOCS3RNAEqcoffset_pSTAT5_shp1oe, :SOCS3RNAEqcoffset_pSTAT5_socs3oe, :SOCS3RNAEqcscale1_CIS_dr90, :SOCS3RNAEqcscale2_CIS_dr90, :SOCS3RNAEqcscale_CISRNA_foldA, :SOCS3RNAEqcscale_CISRNA_foldB, :SOCS3RNAEqcscale_CISRNA_foldC, :SOCS3RNAEqcscale_CIS_actd, :SOCS3RNAEqcscale_CIS_cisoe, :SOCS3RNAEqcscale_CIS_long, :SOCS3RNAEqcscale_CIS_shp1oe, :SOCS3RNAEqcscale_CIS_socs3oe, :SOCS3RNAEqcscale_SHP1_shp1oe, :SOCS3RNAEqcscale_SOCS3RNA_foldA, :SOCS3RNAEqcscale_SOCS3RNA_foldB, :SOCS3RNAEqcscale_SOCS3RNA_foldC, :SOCS3RNAEqcscale_SOCS3_cisoe, :SOCS3RNAEqcscale_SOCS3_long, :SOCS3RNAEqcscale_SOCS3_socs3oe, :SOCS3RNAEqcscale_pEpoR_actd, :SOCS3RNAEqcscale_pEpoR_cisoe, :SOCS3RNAEqcscale_pEpoR_cisoe_pepor, :SOCS3RNAEqcscale_pEpoR_dr30, :SOCS3RNAEqcscale_pEpoR_dr7, :SOCS3RNAEqcscale_pEpoR_fine, :SOCS3RNAEqcscale_pEpoR_long, :SOCS3RNAEqcscale_pEpoR_shp1oe, :SOCS3RNAEqcscale_pEpoR_socs3oe, :SOCS3RNAEqcscale_pJAK2_actd, :SOCS3RNAEqcscale_pJAK2_cisoe, :SOCS3RNAEqcscale_pJAK2_dr30, :SOCS3RNAEqcscale_pJAK2_dr7, :SOCS3RNAEqcscale_pJAK2_fine, :SOCS3RNAEqcscale_pJAK2_long, :SOCS3RNAEqcscale_pJAK2_shp1oe, :SOCS3RNAEqcscale_pJAK2_socs3oe, :SOCS3RNAEqcscale_pSTAT5_actd, :SOCS3RNAEqcscale_pSTAT5_cisoe, :SOCS3RNAEqcscale_pSTAT5_dr10, :SOCS3RNAEqcscale_pSTAT5_long, :SOCS3RNAEqcscale_pSTAT5_shp1oe, :SOCS3RNAEqcscale_pSTAT5_socs3oe, :SOCS3RNAEqcscale_tSTAT5_actd, :SOCS3RNAEqcscale_tSTAT5_long, :SOCS3RNAEqcscale_tSTAT5_shp1oe, :SOCS3RNAEqcsd_CIS_abs, :SOCS3RNAEqcsd_CIS_au, :SOCS3RNAEqcsd_JAK2EpoR_au, :SOCS3RNAEqcsd_RNA_fold, :SOCS3RNAEqcsd_SHP1_abs, :SOCS3RNAEqcsd_SHP1_au, :SOCS3RNAEqcsd_SOCS3_abs, :SOCS3RNAEqcsd_SOCS3_au, :SOCS3RNAEqcsd_STAT5_abs, :SOCS3RNAEqcsd_STAT5_au, :SOCS3RNAEqcsd_pSTAT5_rel, :SOCS3RNAEqcsd_pSTAT5_socs3oe, :SOCS3RNATurnSOCS3RNAEqc, :SOCS3TurnSOCS3RNAEqc, :STAT5ActEpoRSOCS3RNAEqc, :STAT5ActJAK2SOCS3RNAEqc, :STAT5ExpSOCS3RNAEqc, :STAT5ImpSOCS3RNAEqc, :init_EpoRJAK2SOCS3RNAEqc, :init_SHP1SOCS3RNAEqc, :init_STAT5SOCS3RNAEqc, :offset_CIS_actdSOCS3RNAEqc, :offset_CIS_cisoeSOCS3RNAEqc, :offset_CIS_longSOCS3RNAEqc, :offset_CIS_shp1oeSOCS3RNAEqc, :offset_CIS_socs3oeSOCS3RNAEqc, :offset_SOCS3_cisoeSOCS3RNAEqc, :offset_SOCS3_longSOCS3RNAEqc, :offset_SOCS3_socs3oeSOCS3RNAEqc, :offset_pEpoR_actdSOCS3RNAEqc, :offset_pEpoR_cisoeSOCS3RNAEqc, :offset_pEpoR_cisoe_peporSOCS3RNAEqc, :offset_pEpoR_dr30SOCS3RNAEqc, :offset_pEpoR_dr7SOCS3RNAEqc, :offset_pEpoR_fineSOCS3RNAEqc, :offset_pEpoR_longSOCS3RNAEqc, :offset_pEpoR_shp1oeSOCS3RNAEqc, :offset_pEpoR_socs3oeSOCS3RNAEqc, :offset_pJAK2_actdSOCS3RNAEqc, :offset_pJAK2_cisoeSOCS3RNAEqc, :offset_pJAK2_dr30SOCS3RNAEqc, :offset_pJAK2_dr7SOCS3RNAEqc, :offset_pJAK2_fineSOCS3RNAEqc, :offset_pJAK2_longSOCS3RNAEqc, :offset_pJAK2_shp1oeSOCS3RNAEqc, :offset_pJAK2_socs3oeSOCS3RNAEqc, :offset_pSTAT5_actdSOCS3RNAEqc, :offset_pSTAT5_cisoeSOCS3RNAEqc, :offset_pSTAT5_concSOCS3RNAEqc, :offset_pSTAT5_longSOCS3RNAEqc, :offset_pSTAT5_shp1oeSOCS3RNAEqc, :offset_pSTAT5_socs3oeSOCS3RNAEqc, :scale1_CIS_dr90SOCS3RNAEqc, :scale2_CIS_dr90SOCS3RNAEqc, :scale_CISRNA_foldASOCS3RNAEqc, :scale_CISRNA_foldBSOCS3RNAEqc, :scale_CISRNA_foldCSOCS3RNAEqc, :scale_CIS_actdSOCS3RNAEqc, :scale_CIS_cisoeSOCS3RNAEqc, :scale_CIS_longSOCS3RNAEqc, :scale_CIS_shp1oeSOCS3RNAEqc, :scale_CIS_socs3oeSOCS3RNAEqc, :scale_SHP1_shp1oeSOCS3RNAEqc, :scale_SOCS3RNA_foldASOCS3RNAEqc, :scale_SOCS3RNA_foldBSOCS3RNAEqc, :scale_SOCS3RNA_foldCSOCS3RNAEqc, :scale_SOCS3_cisoeSOCS3RNAEqc, :scale_SOCS3_longSOCS3RNAEqc, :scale_SOCS3_socs3oeSOCS3RNAEqc, :scale_pEpoR_actdSOCS3RNAEqc, :scale_pEpoR_cisoeSOCS3RNAEqc, :scale_pEpoR_cisoe_peporSOCS3RNAEqc, :scale_pEpoR_dr30SOCS3RNAEqc, :scale_pEpoR_dr7SOCS3RNAEqc, :scale_pEpoR_fineSOCS3RNAEqc, :scale_pEpoR_longSOCS3RNAEqc, :scale_pEpoR_shp1oeSOCS3RNAEqc, :scale_pEpoR_socs3oeSOCS3RNAEqc, :scale_pJAK2_actdSOCS3RNAEqc, :scale_pJAK2_cisoeSOCS3RNAEqc, :scale_pJAK2_dr30SOCS3RNAEqc, :scale_pJAK2_dr7SOCS3RNAEqc, :scale_pJAK2_fineSOCS3RNAEqc, :scale_pJAK2_longSOCS3RNAEqc, :scale_pJAK2_shp1oeSOCS3RNAEqc, :scale_pJAK2_socs3oeSOCS3RNAEqc, :scale_pSTAT5_actdSOCS3RNAEqc, :scale_pSTAT5_cisoeSOCS3RNAEqc, :scale_pSTAT5_dr10SOCS3RNAEqc, :scale_pSTAT5_longSOCS3RNAEqc, :scale_pSTAT5_shp1oeSOCS3RNAEqc, :scale_pSTAT5_socs3oeSOCS3RNAEqc, :scale_tSTAT5_actdSOCS3RNAEqc, :scale_tSTAT5_longSOCS3RNAEqc, :scale_tSTAT5_shp1oeSOCS3RNAEqc, :sd_CIS_absSOCS3RNAEqc, :sd_CIS_auSOCS3RNAEqc, :sd_JAK2EpoR_auSOCS3RNAEqc, :sd_RNA_foldSOCS3RNAEqc, :sd_SHP1_absSOCS3RNAEqc, :sd_SHP1_auSOCS3RNAEqc, :sd_SOCS3_absSOCS3RNAEqc, :sd_SOCS3_auSOCS3RNAEqc, :sd_STAT5_absSOCS3RNAEqc, :sd_STAT5_auSOCS3RNAEqc, :sd_pSTAT5_relSOCS3RNAEqc, :sd_pSTAT5_socs3oeSOCS3RNAEqc, :CISEqcCISRNAEqc, :CISEqcOECISRNAEqc, :CISInhCISRNAEqc, :CISRNADelayCISRNAEqc, :CISRNAEqcCISEqc, :CISRNAEqcCISEqcOE, :CISRNAEqcCISInh, :CISRNAEqcCISRNADelay, :CISRNAEqcCISRNAEqc, :CISRNAEqcCISRNATurn, :CISRNAEqcCISTurn, :CISRNAEqcEpoRActJAK2, :CISRNAEqcEpoRCISInh, :CISRNAEqcEpoRCISRemove, :CISRNAEqcJAK2ActEpo, :CISRNAEqcJAK2EpoRDeaSHP1, :CISRNAEqcSHP1ActEpoR, :CISRNAEqcSHP1Dea, :CISRNAEqcSHP1ProOE, :CISRNAEqcSOCS3Eqc, :CISRNAEqcSOCS3EqcOE, :CISRNAEqcSOCS3Inh, :CISRNAEqcSOCS3RNADelay, :CISRNAEqcSOCS3RNATurn, :CISRNAEqcSOCS3Turn, :CISRNAEqcSTAT5ActEpoR, :CISRNAEqcSTAT5ActJAK2, :CISRNAEqcSTAT5Exp, :CISRNAEqcSTAT5Imp, :CISRNAEqcinit_EpoRJAK2, :CISRNAEqcinit_SHP1, :CISRNAEqcinit_STAT5, :CISRNAEqcoffset_CIS_actd, :CISRNAEqcoffset_CIS_cisoe, :CISRNAEqcoffset_CIS_long, :CISRNAEqcoffset_CIS_shp1oe, :CISRNAEqcoffset_CIS_socs3oe, :CISRNAEqcoffset_SOCS3_cisoe, :CISRNAEqcoffset_SOCS3_long, :CISRNAEqcoffset_SOCS3_socs3oe, :CISRNAEqcoffset_pEpoR_actd, :CISRNAEqcoffset_pEpoR_cisoe, :CISRNAEqcoffset_pEpoR_cisoe_pepor, :CISRNAEqcoffset_pEpoR_dr30, :CISRNAEqcoffset_pEpoR_dr7, :CISRNAEqcoffset_pEpoR_fine, :CISRNAEqcoffset_pEpoR_long, :CISRNAEqcoffset_pEpoR_shp1oe, :CISRNAEqcoffset_pEpoR_socs3oe, :CISRNAEqcoffset_pJAK2_actd, :CISRNAEqcoffset_pJAK2_cisoe, :CISRNAEqcoffset_pJAK2_dr30, :CISRNAEqcoffset_pJAK2_dr7, :CISRNAEqcoffset_pJAK2_fine, :CISRNAEqcoffset_pJAK2_long, :CISRNAEqcoffset_pJAK2_shp1oe, :CISRNAEqcoffset_pJAK2_socs3oe, :CISRNAEqcoffset_pSTAT5_actd, :CISRNAEqcoffset_pSTAT5_cisoe, :CISRNAEqcoffset_pSTAT5_conc, :CISRNAEqcoffset_pSTAT5_long, :CISRNAEqcoffset_pSTAT5_shp1oe, :CISRNAEqcoffset_pSTAT5_socs3oe, :CISRNAEqcscale1_CIS_dr90, :CISRNAEqcscale2_CIS_dr90, :CISRNAEqcscale_CISRNA_foldA, :CISRNAEqcscale_CISRNA_foldB, :CISRNAEqcscale_CISRNA_foldC, :CISRNAEqcscale_CIS_actd, :CISRNAEqcscale_CIS_cisoe, :CISRNAEqcscale_CIS_long, :CISRNAEqcscale_CIS_shp1oe, :CISRNAEqcscale_CIS_socs3oe, :CISRNAEqcscale_SHP1_shp1oe, :CISRNAEqcscale_SOCS3RNA_foldA, :CISRNAEqcscale_SOCS3RNA_foldB, :CISRNAEqcscale_SOCS3RNA_foldC, :CISRNAEqcscale_SOCS3_cisoe, :CISRNAEqcscale_SOCS3_long, :CISRNAEqcscale_SOCS3_socs3oe, :CISRNAEqcscale_pEpoR_actd, :CISRNAEqcscale_pEpoR_cisoe, :CISRNAEqcscale_pEpoR_cisoe_pepor, :CISRNAEqcscale_pEpoR_dr30, :CISRNAEqcscale_pEpoR_dr7, :CISRNAEqcscale_pEpoR_fine, :CISRNAEqcscale_pEpoR_long, :CISRNAEqcscale_pEpoR_shp1oe, :CISRNAEqcscale_pEpoR_socs3oe, :CISRNAEqcscale_pJAK2_actd, :CISRNAEqcscale_pJAK2_cisoe, :CISRNAEqcscale_pJAK2_dr30, :CISRNAEqcscale_pJAK2_dr7, :CISRNAEqcscale_pJAK2_fine, :CISRNAEqcscale_pJAK2_long, :CISRNAEqcscale_pJAK2_shp1oe, :CISRNAEqcscale_pJAK2_socs3oe, :CISRNAEqcscale_pSTAT5_actd, :CISRNAEqcscale_pSTAT5_cisoe, :CISRNAEqcscale_pSTAT5_dr10, :CISRNAEqcscale_pSTAT5_long, :CISRNAEqcscale_pSTAT5_shp1oe, :CISRNAEqcscale_pSTAT5_socs3oe, :CISRNAEqcscale_tSTAT5_actd, :CISRNAEqcscale_tSTAT5_long, :CISRNAEqcscale_tSTAT5_shp1oe, :CISRNAEqcsd_CIS_abs, :CISRNAEqcsd_CIS_au, :CISRNAEqcsd_JAK2EpoR_au, :CISRNAEqcsd_RNA_fold, :CISRNAEqcsd_SHP1_abs, :CISRNAEqcsd_SHP1_au, :CISRNAEqcsd_SOCS3_abs, :CISRNAEqcsd_SOCS3_au, :CISRNAEqcsd_STAT5_abs, :CISRNAEqcsd_STAT5_au, :CISRNAEqcsd_pSTAT5_rel, :CISRNAEqcsd_pSTAT5_socs3oe, :CISRNATurnCISRNAEqc, :CISTurnCISRNAEqc, :EpoRActJAK2CISRNAEqc, :EpoRCISInhCISRNAEqc, :EpoRCISRemoveCISRNAEqc, :JAK2ActEpoCISRNAEqc, :JAK2EpoRDeaSHP1CISRNAEqc, :SHP1ActEpoRCISRNAEqc, :SHP1DeaCISRNAEqc, :SHP1ProOECISRNAEqc, :SOCS3EqcCISRNAEqc, :SOCS3EqcOECISRNAEqc, :SOCS3InhCISRNAEqc, :SOCS3RNADelayCISRNAEqc, :SOCS3RNATurnCISRNAEqc, :SOCS3TurnCISRNAEqc, :STAT5ActEpoRCISRNAEqc, :STAT5ActJAK2CISRNAEqc, :STAT5ExpCISRNAEqc, :STAT5ImpCISRNAEqc, :init_EpoRJAK2CISRNAEqc, :init_SHP1CISRNAEqc, :init_STAT5CISRNAEqc, :offset_CIS_actdCISRNAEqc, :offset_CIS_cisoeCISRNAEqc, :offset_CIS_longCISRNAEqc, :offset_CIS_shp1oeCISRNAEqc, :offset_CIS_socs3oeCISRNAEqc, :offset_SOCS3_cisoeCISRNAEqc, :offset_SOCS3_longCISRNAEqc, :offset_SOCS3_socs3oeCISRNAEqc, :offset_pEpoR_actdCISRNAEqc, :offset_pEpoR_cisoeCISRNAEqc, :offset_pEpoR_cisoe_peporCISRNAEqc, :offset_pEpoR_dr30CISRNAEqc, :offset_pEpoR_dr7CISRNAEqc, :offset_pEpoR_fineCISRNAEqc, :offset_pEpoR_longCISRNAEqc, :offset_pEpoR_shp1oeCISRNAEqc, :offset_pEpoR_socs3oeCISRNAEqc, :offset_pJAK2_actdCISRNAEqc, :offset_pJAK2_cisoeCISRNAEqc, :offset_pJAK2_dr30CISRNAEqc, :offset_pJAK2_dr7CISRNAEqc, :offset_pJAK2_fineCISRNAEqc, :offset_pJAK2_longCISRNAEqc, :offset_pJAK2_shp1oeCISRNAEqc, :offset_pJAK2_socs3oeCISRNAEqc, :offset_pSTAT5_actdCISRNAEqc, :offset_pSTAT5_cisoeCISRNAEqc, :offset_pSTAT5_concCISRNAEqc, :offset_pSTAT5_longCISRNAEqc, :offset_pSTAT5_shp1oeCISRNAEqc, :offset_pSTAT5_socs3oeCISRNAEqc, :scale1_CIS_dr90CISRNAEqc, :scale2_CIS_dr90CISRNAEqc, :scale_CISRNA_foldACISRNAEqc, :scale_CISRNA_foldBCISRNAEqc, :scale_CISRNA_foldCCISRNAEqc, :scale_CIS_actdCISRNAEqc, :scale_CIS_cisoeCISRNAEqc, :scale_CIS_longCISRNAEqc, :scale_CIS_shp1oeCfunction peTabFormulaToJulia(formula::String, stateNames, paramData::ParamData, namesParamDyn::Array{String, 1}ISRNAEqc, :scale_CIS_socs3oeCISRNAEqc, :scale_SHP1_shp1oeCISRNAEqc, :scale_SOCS3RNA_foldACISRNAEqc, :scale_SOCS3RNA_foldBCISRNAEqc, :scale_SOCS3RNA_foldCCISRNAEqc, :scale_SOCS3_cisoeCISRNAEqc, :scale_SOCS3_longCISRNAEqc, :scale_SOCS3_socs3oeCISRNAEqc, :scale_pEpoR_actdCISRNAEqc, :scale_pEpoR_cisoeCISRNAEqc, :scale_pEpoR_cisoe_peporCISRNAEqc, :scale_pEpoR_dr30CISRNAEqc, :scale_pEpoR_dr7CISRNAEqc, :scale_pEpoR_fineCISRNAEqc, :scale_pEpoR_longCISRNAEqc, :scale_pEpoR_shp1oeCISRNAEqc, :scale_pEpoR_socs3oeCISRNAEqc, :scale_pJAK2_actdCISRNAEqc, :scale_pJAK2_cisoeCISRNAEqc, :scale_pJAK2_dr30CISRNAEqc, :scale_pJAK2_dr7CISRNAEqc, :scale_pJAK2_fineCISRNAEqc, :scale_pJAK2_longCISRNAEqc, :scale_pJAK2_shp1oeCISRNAEqc, :scale_pJAK2_socs3oeCISRNAEqc, :scale_pSTAT5_actdCISRNAEqc, :scale_pSTAT5_cisoeCISRNAEqc, :scale_pSTAT5_dr10CISRNAEqc, :scale_pSTAT5_longCISRNAEqc, :scale_pSTAT5_shp1oeCISRNAEqc, :scale_pSTAT5_socs3oeCISRNAEqc, :scale_tSTAT5_actdCISRNAEqc, :scale_tSTAT5_longCISRNAEqc, :scale_tSTAT5_shp1oeCISRNAEqc, :sd_CIS_absCISRNAEqc, :sd_CIS_auCISRNAEqc, :sd_JAK2EpoR_auCISRNAEqc, :sd_RNA_foldCISRNAEqc, :sd_SHP1_absCISRNAEqc, :sd_SHP1_auCISRNAEqc, :sd_SOCS3_absCISRNAEqc, :sd_SOCS3_auCISRNAEqc, :sd_STAT5_absCISRNAEqc, :sd_STAT5_auCISRNAEqc, :sd_pSTAT5_relCISRNAEqc, :sd_pSTAT5_socs3oeCISRNAEqc])]

    # For correct indexing when comparing gradient or when inputing PyPesto vector to Julia 
    iUse = [findfirst(x -> x == string(peTabOpt.namesParam[i]), names(paramMat)) for i in eachindex(peTabOpt.namesParam)]
    nParam = ncol(paramMat)

    # For correct indexin when comparing hessian
    namesParamHess = Array{String, 1}(undef,nParam*nParam)
    for k in eachindex(namesParamHess)
        row = 1 + (k-1) % nParam
        col = 1 + (k-1) ÷ nParam
        namesParamHess[k] = peTabOpt.namesParam[row] * peTabOpt.namesParam[col]
    end
    iUseHess = [findfirst(x -> x == namesParamHess[i], names(hessPythonMat)) for i in eachindex(namesParamHess)]

    for i in 1:nrow(paramMat)
        
        paramVec = collect(paramMat[i, :])[iUse]

        # Test cost 
        costJulia = peTabOpt.evalF(paramVec)
        sqDiffCost = (costJulia - costPython[i])^2
        if sqDiffCost > 1e-5
            @printf("sqDiffCost = %.3e\n", sqDiffCost)
            @printf("Does not pass test on cost\n")
            return false
        end

        grad = zeros(nParam)
        gradPython = collect(gradPythonMat[i, :])
        peTabOpt.evalGradF(grad, paramVec)
        sqDiffGrad = sum((grad - gradPython[iUse]).^2)
        if sqDiffGrad > 1e-5
            @printf("sqDiffGrad = %.3e\n", sqDiffGrad)
            @printf("Does not pass test on grad\n")
            return false
        end

        hessJulia = zeros(nParam,nParam); peTabOptAlt.evalHessGaussNewton(hessJulia, paramVec)
        hessJulia = hessJulia[:]
        hessPython = collect(hessPythonMat[i, :])
        sqDiffHess = sum((hessJulia - hessPython[iUseHess]).^2)
        if sqDiffHess > 1e-5
            @printf("sqDiffHess = %.3e\n", sqDiffHess)
            @printf("Does not pass test on hessian\n")
            return false
        end

        costZygote = peTabOpt.evalFZygote(paramVec)
        sqDiffZygote = (costZygote - costPython[i])^2
        if sqDiffZygote > 1e-5
            @printf("sqDiffCostZygote = %.3e\n", sqDiffZygote)
            @printf("Does not pass test on cost for Zygote\n")
            return false
        end

        # Evaluate lower level adjoint sensitivity interface gradient 
        gradAdj = zeros(length(paramVec))
        peTabOpt.evalGradFAdjoint(gradAdj, paramVec)
        sqDiffGradAdjoint = sum((gradAdj - gradPython[iUse]).^2)
        if sqDiffGradAdjoint > 1e-4
            @printf("sqDiffGradAdjoint = %.3e\n", sqDiffGradAdjoint)
            @printf("Does not pass test on adjoint gradient gradient\n")
            return false
        end

        # Evaluate forward sensitivity equations 
        gradForwardEq = zeros(length(paramVec))
        peTabOpt.evalGradFForwardEq(gradForwardEq, paramVec)
        gradPython = collect(gradPythonMat[i, :])
        sqDiffGradForwardEq = sum((gradForwardEq - gradPython[iUse]).^2)
        if sqDiffGradForwardEq > 1e-5
            @printf("sqDiffGradForwardEq = %.3e\n", sqDiffGradForwardEq)
            @printf("Does not pass test on gradient from Forward sensitivity equations\n")
            return false
        end

        # Forward sensitivity equations using autodiff 
        gradForwardEqAuto = zeros(length(paramVec))
        peTabOptAlt.evalGradFForwardEq(gradForwardEqAuto, paramVec)
        sqDiffGradForwardEqAuto = sum((gradForwardEq - gradPython[iUse]).^2)
        if sqDiffGradForwardEqAuto > 1e-5
            @printf("sqDiffGradForwardEqAuto = %.3e\n", sqDiffGradForwardEqAuto)
            @printf("Does not pass test on gradient from Forward sensitivity equations\n")
            return false
        end

        # This currently takes considerble time (hence it does not run for all passes)
        if i == 1
            gradZygote = zeros(length(paramVec))
            peTabOpt.evalGradFZygote(gradZygote, paramVec)
            gradPython = collect(gradPythonMat[i, :])
            sqDiffGradZygote = sum((gradZygote - gradPython[iUse]).^2)
            if sqDiffGradZygote > 1e-5
                @printf("sqDiffGradZygote = %.3e\n", sqDiffGradZygote)
                @printf("Does not pass test on gradient from Zygote\n")
                return false
            end
            @printf("sqDiffGradZygote = %.3e\n", sqDiffGradZygote)
        end

        if printRes == true
            @printf("sqDiffCost = %.3e\n", sqDiffCost)
            @printf("sqDiffGrad = %.3e\n", sqDiffGrad)
            @printf("sqDiffHess = %.3e\n", sqDiffHess)
            @printf("sqDiffCostZygote = %.3e\n", sqDiffZygote)
            @printf("sqDiffGradAdjoint = %.3e\n", sqDiffGradAdjoint)
            @printf("sqDiffGradForwardEq = %.3e\n", sqDiffGradForwardEq)
            @printf("sqDiffGradForwardEqAuto = %.3e\n", sqDiffGradForwardEqAuto)
        end
    end

    return true
end


peTabModel = setUpPeTabModel("Bachmann_MSB2011", pwd() * "/tests/Bachman/", forceBuildJlFile=true)
passTest = compareAgainstPyPesto(peTabModel, Rodas5(), 1e-12, printRes=true)
if passTest == true
    @printf("Passed test Bachaman against PyPesto\n")
else
    @printf("Did not pass test Bachman against PyPesto\n")
end
