#u[1] = P_PA, u[2] = Vp_cyt, u[3] = Vp_nuc, u[4] = P_NA, u[5] = R_C_RdRp, u[6] = R_M6, u[7] = P_NEP, u[8] = V_end, u[9] = B_att_Hi, u[10] = R_M5, u[11] = Vp_cyt_M1, u[12] = V_rel, u[13] = R_C, u[14] = Vp_nuc_M1, u[15] = Cp, u[16] = V_att_Lo, u[17] = V_att_Hi, u[18] = B_att_Lo, u[19] = P_HA, u[20] = P_M1, u[21] = R_M1, u[22] = P_M2, u[23] = R_M7, u[24] = P_NP, u[25] = V_ex, u[26] = R_M4, u[27] = R_M8, u[28] = R_V_RdRp, u[29] = R_M2, u[30] = R_V, u[31] = P_B2, u[32] = P_RdRp, u[33] = R_M3, u[34] = P_B1
#pODEProblem[1] = k_deg_R_RdRp, pODEProblem[2] = ModelValue_104, pODEProblem[3] = ModelValue_63, pODEProblem[4] = ModelValue_114, pODEProblem[5] = ModelValue_108, pODEProblem[6] = ModelValue_64, pODEProblem[7] = L6, pODEProblem[8] = k_end, pODEProblem[9] = N_P_RdRp, pODEProblem[10] = D_rib, pODEProblem[11] = ModelValue_111, pODEProblem[12] = ModelValue_107, pODEProblem[13] = k_rel, pODEProblem[14] = ModelValue_105, pODEProblem[15] = compartment, pODEProblem[16] = L1, pODEProblem[17] = K_eq_Hi, pODEProblem[18] = N_P_NP, pODEProblem[19] = L2, pODEProblem[20] = ModelValue_113, pODEProblem[21] = k_bind_M1, pODEProblem[22] = ModelValue_116, pODEProblem[23] = F_Spl7, pODEProblem[24] = k_deg_Rnp, pODEProblem[25] = N_P_NA, pODEProblem[26] = K_eq_Lo, pODEProblem[27] = ModelValue_101, pODEProblem[28] = k_exp_Vp_nuc_M1, pODEProblem[29] = ModelValue_79, pODEProblem[30] = ModelValue_69, pODEProblem[31] = ModelValue_90, pODEProblem[32] = ModelValue_91, pODEProblem[33] = ModelValue_103, pODEProblem[34] = N_P_M1, pODEProblem[35] = k_syn_R_M, pODEProblem[36] = L8, pODEProblem[37] = ModelValue_84, pODEProblem[38] = ModelValue_85, pODEProblem[39] = ModelValue_115, pODEProblem[40] = k_att_Lo, pODEProblem[41] = k_fus, pODEProblem[42] = N_P_HA, pODEProblem[43] = k_deg_R, pODEProblem[44] = N_P_NEP, pODEProblem[45] = k_bind_RdRp, pODEProblem[46] = ModelValue_80, pODEProblem[47] = ModelValue_89, pODEProblem[48] = k_RdRp, pODEProblem[49] = L5, pODEProblem[50] = N_P_M2, pODEProblem[51] = k_syn_R_V, pODEProblem[52] = L4, pODEProblem[53] = ModelValue_87, pODEProblem[54] = k_bind_NP, pODEProblem[55] = L3, pODEProblem[56] = k_imp, pODEProblem[57] = ModelValue_106, pODEProblem[58] = F_Spl8, pODEProblem[59] = k_deg_R_M, pODEProblem[60] = k_syn_P, pODEProblem[61] = ModelValue_82, pODEProblem[62] = k_syn_R_C, pODEProblem[63] = L7, pODEProblem[64] = K_V_rel, pODEProblem[65] = ModelValue_86, pODEProblem[66] = ModelValue_102, pODEProblem[67] = F_fus, pODEProblem[68] = k_att_Hi, pODEProblem[69] = ModelValue_88
#
function compute_∂h∂u!(u, t::Real, pODEProblem::AbstractVector, θ_observable::AbstractVector,                    
                       θ_nonDynamic::AbstractVector, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :RM5 
		out[10] = 1
		return nothing
	end

	if observableId == :RVSegTot 
		observableParameter1_RVSegTot = getObsOrSdParam(θ_observable, parameterMap)
		out[2] = 1//8
		out[3] = 1//8
		out[8] = 1
		out[11] = 1//8
		out[14] = 1//8
		out[16] = 1
		out[17] = 1
		out[30] = 1//8
		return nothing
	end

	if observableId == :RCSegTot 
		out[5] = 1//8
		out[13] = 1//8
		out[15] = 1//8
		return nothing
	end

	if observableId == :Vrel 
		out[12] = 1
		return nothing
	end

	if observableId == :IntNucOffset 
		return nothing
	end

	if observableId == :FracNucInt_1 
		observableParameter1_FracNucInt_1 = getObsOrSdParam(θ_observable, parameterMap)
		out[2] = (-100u[3] - 100u[14]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[3] = (800.0u[8] + 100.0u[2] + 100.0u[11]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[8] = (-800u[3] - 800u[14]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[11] = (-100u[3] - 100u[14]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[14] = (800.0u[8] + 100.0u[2] + 100.0u[11]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		return nothing
	end

	if observableId == :FracNucInt_2 
		observableParameter1_FracNucInt_2 = getObsOrSdParam(θ_observable, parameterMap)
		out[2] = (-100u[3] - 100u[14]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[3] = (800.0u[8] + 100.0u[2] + 100.0u[11]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[8] = (-800u[3] - 800u[14]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[11] = (-100u[3] - 100u[14]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[14] = (800.0u[8] + 100.0u[2] + 100.0u[11]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		return nothing
	end

	if observableId == :FracNucInt_3 
		observableParameter1_FracNucInt_3 = getObsOrSdParam(θ_observable, parameterMap)
		out[2] = (-100u[3] - 100u[14]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[3] = (800.0u[8] + 100.0u[2] + 100.0u[11]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[8] = (-800u[3] - 800u[14]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[11] = (-100u[3] - 100u[14]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[14] = (800.0u[8] + 100.0u[2] + 100.0u[11]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		return nothing
	end

	if observableId == :FracNucInt_4 
		observableParameter1_FracNucInt_4 = getObsOrSdParam(θ_observable, parameterMap)
		out[2] = (-100u[3] - 100u[14]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[3] = (800.0u[8] + 100.0u[2] + 100.0u[11]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[8] = (-800u[3] - 800u[14]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[11] = (-100u[3] - 100u[14]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[14] = (800.0u[8] + 100.0u[2] + 100.0u[11]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		return nothing
	end

	if observableId == :FracNucInt_5 
		observableParameter1_FracNucInt_5 = getObsOrSdParam(θ_observable, parameterMap)
		out[2] = (-100u[3] - 100u[14]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[3] = (800.0u[8] + 100.0u[2] + 100.0u[11]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[8] = (-800u[3] - 800u[14]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[11] = (-100u[3] - 100u[14]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[14] = (800.0u[8] + 100.0u[2] + 100.0u[11]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		return nothing
	end

	if observableId == :FracNucInt_6 
		observableParameter1_FracNucInt_6 = getObsOrSdParam(θ_observable, parameterMap)
		out[2] = (-100u[3] - 100u[14]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[3] = (800.0u[8] + 100.0u[2] + 100.0u[11]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[8] = (-800u[3] - 800u[14]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[11] = (-100u[3] - 100u[14]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[14] = (800.0u[8] + 100.0u[2] + 100.0u[11]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		return nothing
	end

	if observableId == :FracNucInt_7 
		observableParameter1_FracNucInt_7 = getObsOrSdParam(θ_observable, parameterMap)
		out[2] = (-100u[3] - 100u[14]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[3] = (800.0u[8] + 100.0u[2] + 100.0u[11]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[8] = (-800u[3] - 800u[14]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[11] = (-100u[3] - 100u[14]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[14] = (800.0u[8] + 100.0u[2] + 100.0u[11]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		return nothing
	end

	if observableId == :FracNucInt_8 
		observableParameter1_FracNucInt_8 = getObsOrSdParam(θ_observable, parameterMap)
		out[2] = (-100u[3] - 100u[14]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[3] = (800.0u[8] + 100.0u[2] + 100.0u[11]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[8] = (-800u[3] - 800u[14]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[11] = (-100u[3] - 100u[14]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		out[14] = (800.0u[8] + 100.0u[2] + 100.0u[11]) / ((u[2] + u[11] + u[3] + u[14] + 8u[8])^2)
		return nothing
	end

end

function compute_∂h∂p!(u, t::Real, pODEProblem::AbstractVector, θ_observable::AbstractVector,
                       θ_nonDynamic::AbstractVector, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :RM5 
		return nothing
	end

	if observableId == :RVSegTot 
		return nothing
	end

	if observableId == :RCSegTot 
		return nothing
	end

	if observableId == :Vrel 
		return nothing
	end

	if observableId == :IntNucOffset 
		return nothing
	end

	if observableId == :FracNucInt_1 
		return nothing
	end

	if observableId == :FracNucInt_2 
		return nothing
	end

	if observableId == :FracNucInt_3 
		return nothing
	end

	if observableId == :FracNucInt_4 
		return nothing
	end

	if observableId == :FracNucInt_5 
		return nothing
	end

	if observableId == :FracNucInt_6 
		return nothing
	end

	if observableId == :FracNucInt_7 
		return nothing
	end

	if observableId == :FracNucInt_8 
		return nothing
	end

end

function compute_∂σ∂σu!(u, t::Real, θ_sd::AbstractVector, pODEProblem::AbstractVector, θ_nonDynamic::AbstractVector, 
                        parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :RM5 
		return nothing
	end

	if observableId == :RVSegTot 
		return nothing
	end

	if observableId == :RCSegTot 
		return nothing
	end

	if observableId == :Vrel 
		return nothing
	end

	if observableId == :IntNucOffset 
		return nothing
	end

	if observableId == :FracNucInt_1 
		return nothing
	end

	if observableId == :FracNucInt_2 
		return nothing
	end

	if observableId == :FracNucInt_3 
		return nothing
	end

	if observableId == :FracNucInt_4 
		return nothing
	end

	if observableId == :FracNucInt_5 
		return nothing
	end

	if observableId == :FracNucInt_6 
		return nothing
	end

	if observableId == :FracNucInt_7 
		return nothing
	end

	if observableId == :FracNucInt_8 
		return nothing
	end

end

function compute_∂σ∂σp!(u, t::Real, θ_sd::AbstractVector, pODEProblem::AbstractVector, θ_nonDynamic::AbstractVector, 
                        parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :RM5 
		return nothing
	end

	if observableId == :RVSegTot 
		return nothing
	end

	if observableId == :RCSegTot 
		return nothing
	end

	if observableId == :Vrel 
		return nothing
	end

	if observableId == :IntNucOffset 
		return nothing
	end

	if observableId == :FracNucInt_1 
		return nothing
	end

	if observableId == :FracNucInt_2 
		return nothing
	end

	if observableId == :FracNucInt_3 
		return nothing
	end

	if observableId == :FracNucInt_4 
		return nothing
	end

	if observableId == :FracNucInt_5 
		return nothing
	end

	if observableId == :FracNucInt_6 
		return nothing
	end

	if observableId == :FracNucInt_7 
		return nothing
	end

	if observableId == :FracNucInt_8 
		return nothing
	end

end

