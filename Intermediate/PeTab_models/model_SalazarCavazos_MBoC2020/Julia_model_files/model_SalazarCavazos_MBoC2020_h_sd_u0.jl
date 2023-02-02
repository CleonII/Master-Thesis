#u[1] = species_61, u[2] = species_10, u[3] = species_21, u[4] = species_47, u[5] = species_17, u[6] = species_1, u[7] = species_43, u[8] = species_52, u[9] = species_45, u[10] = species_20, u[11] = species_49, u[12] = species_12, u[13] = species_35, u[14] = species_16, u[15] = species_26, u[16] = species_73, u[17] = species_59, u[18] = species_13, u[19] = species_27, u[20] = species_71, u[21] = species_74, u[22] = species_58, u[23] = species_4, u[24] = species_28, u[25] = species_44, u[26] = species_3, u[27] = species_68, u[28] = species_22, u[29] = species_72, u[30] = species_7, u[31] = species_14, u[32] = species_66, u[33] = species_25, u[34] = species_60, u[35] = species_6, u[36] = species_30, u[37] = species_75, u[38] = species_36, u[39] = species_50, u[40] = species_19, u[41] = species_5, u[42] = species_64, u[43] = species_33, u[44] = species_69, u[45] = species_39, u[46] = species_55, u[47] = species_24, u[48] = species_48, u[49] = species_62, u[50] = species_53, u[51] = species_41, u[52] = species_51, u[53] = species_54, u[54] = species_18, u[55] = species_67, u[56] = species_57, u[57] = species_63, u[58] = species_23, u[59] = species_56, u[60] = species_37, u[61] = species_42, u[62] = species_32, u[63] = species_2, u[64] = species_46, u[65] = species_65, u[66] = species_31, u[67] = species_29, u[68] = species_40, u[69] = species_38, u[70] = species_8, u[71] = species_34, u[72] = species_9, u[73] = species_15, u[74] = species_11, u[75] = species_70
#θ_dynamicNames[1] = GRB2_total__FREE, θ_dynamicNames[2] = SHC1_total__FREE, θ_dynamicNames[3] = kdephosY1068__FREE, θ_dynamicNames[4] = kdephosYN__FREE, θ_dynamicNames[5] = ratio_kpkd_Y1068__FREE, θ_dynamicNames[6] = ratio_kpkd_YN__FREE
#

function compute_h(u::AbstractVector, t::Real, θ_dynamic::AbstractVector, θ_observable::AbstractVector,
                   θ_nonDynamic::AbstractVector, parameterInfo::ParametersInfo, observableId::Symbol, 
                      parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :observable_pY1068_percent 
		return ( 100 * ( u[70] + u[74] + u[31] + u[5] + u[54] + u[10] + u[3] + u[58] + u[47] + u[19] + u[36] + u[66] + u[43] + u[13] + u[38] + u[60] + u[69] + u[68] + u[51] + u[61] + u[7] + u[9] + u[64] + u[48] + u[39] + u[52] + u[8] + u[50] + u[46] + u[59] + u[56] + u[22] + u[17] + u[34] + u[1] + u[49] + u[57] + u[42] + u[32] + u[55] + u[27] + u[44] + u[75] + u[20] + u[29] + u[16] + u[21] + u[37] ) ) / ( u[63] + u[41] + u[35] + u[30] + u[70] + u[72] + u[2] + u[74] + u[12] + u[18] + u[31] + u[73] + u[14] + u[5] + u[54] + u[40] + u[10] + u[3] + u[28] + u[58] + u[47] + u[33] + u[15] + u[19] + u[24] + u[67] + u[36] + u[66] + u[62] + u[43] + u[71] + u[13] + u[38] + u[60] + u[69] + u[45] + u[68] + u[51] + u[61] + u[7] + u[25] + u[9] + u[64] + u[4] + u[48] + u[11] + u[39] + u[52] + u[8] + u[50] + u[53] + u[46] + u[59] + u[56] + u[22] + u[17] + u[34] + u[1] + u[49] + u[57] + u[42] + u[65] + u[32] + u[55] + u[27] + u[44] + u[75] + u[20] + u[29] + u[16] + u[21] + u[37] ) 
	end

	if observableId == :observable_pY1173_percent 
		return ( 100 * ( u[72] + u[12] + u[73] + u[5] + u[40] + u[10] + u[28] + u[33] + u[15] + u[24] + u[36] + u[62] + u[71] + u[13] + u[38] + u[60] + u[45] + u[68] + u[51] + u[61] + u[25] + u[9] + u[4] + u[11] + u[39] + u[52] + u[8] + u[53] + u[46] + u[59] + u[56] + u[22] + u[17] + u[34] + u[1] + u[49] + u[57] + u[65] + u[32] + u[55] + u[27] + u[44] + u[75] + u[20] + u[29] + u[16] + u[21] + u[37] ) ) / ( u[63] + u[41] + u[35] + u[30] + u[70] + u[72] + u[2] + u[74] + u[12] + u[18] + u[31] + u[73] + u[14] + u[5] + u[54] + u[40] + u[10] + u[3] + u[28] + u[58] + u[47] + u[33] + u[15] + u[19] + u[24] + u[67] + u[36] + u[66] + u[62] + u[43] + u[71] + u[13] + u[38] + u[60] + u[69] + u[45] + u[68] + u[51] + u[61] + u[7] + u[25] + u[9] + u[64] + u[4] + u[48] + u[11] + u[39] + u[52] + u[8] + u[50] + u[53] + u[46] + u[59] + u[56] + u[22] + u[17] + u[34] + u[1] + u[49] + u[57] + u[42] + u[65] + u[32] + u[55] + u[27] + u[44] + u[75] + u[20] + u[29] + u[16] + u[21] + u[37] ) 
	end

	if observableId == :observable_unphosR_per 
		return ( 100 * ( u[63] + u[41] + u[35] + u[30] ) ) / ( u[63] + u[41] + u[35] + u[30] + u[70] + u[72] + u[2] + u[74] + u[12] + u[18] + u[31] + u[73] + u[14] + u[5] + u[54] + u[40] + u[10] + u[3] + u[28] + u[58] + u[47] + u[33] + u[15] + u[19] + u[24] + u[67] + u[36] + u[66] + u[62] + u[43] + u[71] + u[13] + u[38] + u[60] + u[69] + u[45] + u[68] + u[51] + u[61] + u[7] + u[25] + u[9] + u[64] + u[4] + u[48] + u[11] + u[39] + u[52] + u[8] + u[50] + u[53] + u[46] + u[59] + u[56] + u[22] + u[17] + u[34] + u[1] + u[49] + u[57] + u[42] + u[65] + u[32] + u[55] + u[27] + u[44] + u[75] + u[20] + u[29] + u[16] + u[21] + u[37] ) 
	end

	if observableId == :observable_phosR_per 
		return ( 100 * ( ( u[63] + u[41] + u[35] + u[30] + u[70] + u[72] + u[2] + u[74] + u[12] + u[18] + u[31] + u[73] + u[14] + u[5] + u[54] + u[40] + u[10] + u[3] + u[28] + u[58] + u[47] + u[33] + u[15] + u[19] + u[24] + u[67] + u[36] + u[66] + u[62] + u[43] + u[71] + u[13] + u[38] + u[60] + u[69] + u[45] + u[68] + u[51] + u[61] + u[7] + u[25] + u[9] + u[64] + u[4] + u[48] + u[11] + u[39] + u[52] + u[8] + u[50] + u[53] + u[46] + u[59] + u[56] + u[22] + u[17] + u[34] + u[1] + u[49] + u[57] + u[42] + u[65] + u[32] + u[55] + u[27] + u[44] + u[75] + u[20] + u[29] + u[16] + u[21] + u[37] ) - ( u[63] + u[41] + u[35] + u[30] ) ) ) / ( u[63] + u[41] + u[35] + u[30] + u[70] + u[72] + u[2] + u[74] + u[12] + u[18] + u[31] + u[73] + u[14] + u[5] + u[54] + u[40] + u[10] + u[3] + u[28] + u[58] + u[47] + u[33] + u[15] + u[19] + u[24] + u[67] + u[36] + u[66] + u[62] + u[43] + u[71] + u[13] + u[38] + u[60] + u[69] + u[45] + u[68] + u[51] + u[61] + u[7] + u[25] + u[9] + u[64] + u[4] + u[48] + u[11] + u[39] + u[52] + u[8] + u[50] + u[53] + u[46] + u[59] + u[56] + u[22] + u[17] + u[34] + u[1] + u[49] + u[57] + u[42] + u[65] + u[32] + u[55] + u[27] + u[44] + u[75] + u[20] + u[29] + u[16] + u[21] + u[37] ) 
	end

end

function compute_u0!(u0::AbstractVector, pODEProblem::AbstractVector) 

	#pODEProblem[1] = volMedia, pODEProblem[2] = SHC1_total__FREE, pODEProblem[3] = kp_SE, pODEProblem[4] = kon_EGF__, pODEProblem[5] = km_GE, pODEProblem[6] = km_dim_L_L, pODEProblem[7] = kon_EGF, pODEProblem[8] = ratio_kphos_receiver, pODEProblem[9] = kp_EGF, pODEProblem[10] = kphos_f, pODEProblem[11] = kp_GE, pODEProblem[12] = kon, pODEProblem[13] = c1, pODEProblem[14] = ratio_kpkd_YN__FREE, pODEProblem[15] = onrate_f, pODEProblem[16] = kon__, pODEProblem[17] = ratio_kpkd_Y1068__FREE, pODEProblem[18] = Kd_EGF, pODEProblem[19] = km_SE, pODEProblem[20] = EGFR_total, pODEProblem[21] = ratio_kphosY1173, pODEProblem[22] = offrate_f, pODEProblem[23] = kp_dim_L_L_pre, pODEProblem[24] = numCells, pODEProblem[25] = Vextra, pODEProblem[26] = EGFconc, pODEProblem[27] = GRB2_total__FREE, pODEProblem[28] = kdephosY1068_f, pODEProblem[29] = f, pODEProblem[30] = Kd_EGF__, pODEProblem[31] = kdephosYN__FREE, pODEProblem[32] = NA, pODEProblem[33] = GRB2_f, pODEProblem[34] = km_EGF, pODEProblem[35] = Kd_GE, pODEProblem[36] = km_dim_L_L_pre, pODEProblem[37] = kdephosY1068__FREE, pODEProblem[38] = kdephosY1173_f, pODEProblem[39] = ratio_kdephosY1173, pODEProblem[40] = KD_dim, pODEProblem[41] = Vc, pODEProblem[42] = kp_dim_L_L, pODEProblem[43] = _rateLaw1, pODEProblem[44] = Kd_SE

	species_61 = 0.0 
	species_10 = 0.0 
	species_21 = 0.0 
	species_47 = 0.0 
	species_17 = 0.0 
	species_1 = pODEProblem[26] * pODEProblem[32] * pODEProblem[25] 
	species_43 = 0.0 
	species_52 = 0.0 
	species_45 = 0.0 
	species_20 = 0.0 
	species_49 = 0.0 
	species_12 = 0.0 
	species_35 = 0.0 
	species_16 = 0.0 
	species_26 = 0.0 
	species_73 = 0.0 
	species_59 = 0.0 
	species_13 = 0.0 
	species_27 = 0.0 
	species_71 = 0.0 
	species_74 = 0.0 
	species_58 = 0.0 
	species_4 = pODEProblem[2] * pODEProblem[29] 
	species_28 = 0.0 
	species_44 = 0.0 
	species_3 = pODEProblem[33] * pODEProblem[27] * pODEProblem[29] 
	species_68 = 0.0 
	species_22 = 0.0 
	species_72 = 0.0 
	species_7 = 0.0 
	species_14 = 0.0 
	species_66 = 0.0 
	species_25 = 0.0 
	species_60 = 0.0 
	species_6 = 0.0 
	species_30 = 0.0 
	species_75 = 0.0 
	species_36 = 0.0 
	species_50 = 0.0 
	species_19 = 0.0 
	species_5 = 0.0 
	species_64 = 0.0 
	species_33 = 0.0 
	species_69 = 0.0 
	species_39 = 0.0 
	species_55 = 0.0 
	species_24 = 0.0 
	species_48 = 0.0 
	species_62 = 0.0 
	species_53 = 0.0 
	species_41 = 0.0 
	species_51 = 0.0 
	species_54 = 0.0 
	species_18 = 0.0 
	species_67 = 0.0 
	species_57 = 0.0 
	species_63 = 0.0 
	species_23 = 0.0 
	species_56 = 0.0 
	species_37 = 0.0 
	species_42 = 0.0 
	species_32 = 0.0 
	species_2 = pODEProblem[20] 
	species_46 = 0.0 
	species_65 = 0.0 
	species_31 = 0.0 
	species_29 = 0.0 
	species_40 = 0.0 
	species_38 = 0.0 
	species_8 = 0.0 
	species_34 = 0.0 
	species_9 = 0.0 
	species_15 = 0.0 
	species_11 = 0.0 
	species_70 = 0.0 

	u0 .= species_61, species_10, species_21, species_47, species_17, species_1, species_43, species_52, species_45, species_20, species_49, species_12, species_35, species_16, species_26, species_73, species_59, species_13, species_27, species_71, species_74, species_58, species_4, species_28, species_44, species_3, species_68, species_22, species_72, species_7, species_14, species_66, species_25, species_60, species_6, species_30, species_75, species_36, species_50, species_19, species_5, species_64, species_33, species_69, species_39, species_55, species_24, species_48, species_62, species_53, species_41, species_51, species_54, species_18, species_67, species_57, species_63, species_23, species_56, species_37, species_42, species_32, species_2, species_46, species_65, species_31, species_29, species_40, species_38, species_8, species_34, species_9, species_15, species_11, species_70
end

function compute_u0(pODEProblem::AbstractVector)::AbstractVector 

	#pODEProblem[1] = volMedia, pODEProblem[2] = SHC1_total__FREE, pODEProblem[3] = kp_SE, pODEProblem[4] = kon_EGF__, pODEProblem[5] = km_GE, pODEProblem[6] = km_dim_L_L, pODEProblem[7] = kon_EGF, pODEProblem[8] = ratio_kphos_receiver, pODEProblem[9] = kp_EGF, pODEProblem[10] = kphos_f, pODEProblem[11] = kp_GE, pODEProblem[12] = kon, pODEProblem[13] = c1, pODEProblem[14] = ratio_kpkd_YN__FREE, pODEProblem[15] = onrate_f, pODEProblem[16] = kon__, pODEProblem[17] = ratio_kpkd_Y1068__FREE, pODEProblem[18] = Kd_EGF, pODEProblem[19] = km_SE, pODEProblem[20] = EGFR_total, pODEProblem[21] = ratio_kphosY1173, pODEProblem[22] = offrate_f, pODEProblem[23] = kp_dim_L_L_pre, pODEProblem[24] = numCells, pODEProblem[25] = Vextra, pODEProblem[26] = EGFconc, pODEProblem[27] = GRB2_total__FREE, pODEProblem[28] = kdephosY1068_f, pODEProblem[29] = f, pODEProblem[30] = Kd_EGF__, pODEProblem[31] = kdephosYN__FREE, pODEProblem[32] = NA, pODEProblem[33] = GRB2_f, pODEProblem[34] = km_EGF, pODEProblem[35] = Kd_GE, pODEProblem[36] = km_dim_L_L_pre, pODEProblem[37] = kdephosY1068__FREE, pODEProblem[38] = kdephosY1173_f, pODEProblem[39] = ratio_kdephosY1173, pODEProblem[40] = KD_dim, pODEProblem[41] = Vc, pODEProblem[42] = kp_dim_L_L, pODEProblem[43] = _rateLaw1, pODEProblem[44] = Kd_SE

	species_61 = 0.0 
	species_10 = 0.0 
	species_21 = 0.0 
	species_47 = 0.0 
	species_17 = 0.0 
	species_1 = pODEProblem[26] * pODEProblem[32] * pODEProblem[25] 
	species_43 = 0.0 
	species_52 = 0.0 
	species_45 = 0.0 
	species_20 = 0.0 
	species_49 = 0.0 
	species_12 = 0.0 
	species_35 = 0.0 
	species_16 = 0.0 
	species_26 = 0.0 
	species_73 = 0.0 
	species_59 = 0.0 
	species_13 = 0.0 
	species_27 = 0.0 
	species_71 = 0.0 
	species_74 = 0.0 
	species_58 = 0.0 
	species_4 = pODEProblem[2] * pODEProblem[29] 
	species_28 = 0.0 
	species_44 = 0.0 
	species_3 = pODEProblem[33] * pODEProblem[27] * pODEProblem[29] 
	species_68 = 0.0 
	species_22 = 0.0 
	species_72 = 0.0 
	species_7 = 0.0 
	species_14 = 0.0 
	species_66 = 0.0 
	species_25 = 0.0 
	species_60 = 0.0 
	species_6 = 0.0 
	species_30 = 0.0 
	species_75 = 0.0 
	species_36 = 0.0 
	species_50 = 0.0 
	species_19 = 0.0 
	species_5 = 0.0 
	species_64 = 0.0 
	species_33 = 0.0 
	species_69 = 0.0 
	species_39 = 0.0 
	species_55 = 0.0 
	species_24 = 0.0 
	species_48 = 0.0 
	species_62 = 0.0 
	species_53 = 0.0 
	species_41 = 0.0 
	species_51 = 0.0 
	species_54 = 0.0 
	species_18 = 0.0 
	species_67 = 0.0 
	species_57 = 0.0 
	species_63 = 0.0 
	species_23 = 0.0 
	species_56 = 0.0 
	species_37 = 0.0 
	species_42 = 0.0 
	species_32 = 0.0 
	species_2 = pODEProblem[20] 
	species_46 = 0.0 
	species_65 = 0.0 
	species_31 = 0.0 
	species_29 = 0.0 
	species_40 = 0.0 
	species_38 = 0.0 
	species_8 = 0.0 
	species_34 = 0.0 
	species_9 = 0.0 
	species_15 = 0.0 
	species_11 = 0.0 
	species_70 = 0.0 

	 return [species_61, species_10, species_21, species_47, species_17, species_1, species_43, species_52, species_45, species_20, species_49, species_12, species_35, species_16, species_26, species_73, species_59, species_13, species_27, species_71, species_74, species_58, species_4, species_28, species_44, species_3, species_68, species_22, species_72, species_7, species_14, species_66, species_25, species_60, species_6, species_30, species_75, species_36, species_50, species_19, species_5, species_64, species_33, species_69, species_39, species_55, species_24, species_48, species_62, species_53, species_41, species_51, species_54, species_18, species_67, species_57, species_63, species_23, species_56, species_37, species_42, species_32, species_2, species_46, species_65, species_31, species_29, species_40, species_38, species_8, species_34, species_9, species_15, species_11, species_70]
end

function compute_σ(u::AbstractVector, t::Real, θ_sd::AbstractVector, θ_dynamic::AbstractVector, θ_nonDynamic::AbstractVector, 
                   parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :observable_pY1068_percent 
		noiseParameter1_observable_pY1068_percent = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_pY1068_percent 
	end

	if observableId == :observable_pY1173_percent 
		noiseParameter1_observable_pY1173_percent = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_pY1173_percent 
	end

	if observableId == :observable_unphosR_per 
		noiseParameter1_observable_unphosR_per = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_unphosR_per 
	end

	if observableId == :observable_phosR_per 
		noiseParameter1_observable_phosR_per = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_observable_phosR_per 
	end

end