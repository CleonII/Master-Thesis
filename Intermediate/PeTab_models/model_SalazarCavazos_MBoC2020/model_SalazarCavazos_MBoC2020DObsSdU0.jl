function evalDYmodDu(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	species_61, species_10, species_21, species_47, species_17, species_1, species_43, species_52, species_45, species_20, species_49, species_12, species_35, species_16, species_26, species_73, species_59, species_13, species_27, species_71, species_74, species_58, species_4, species_28, species_44, species_3, species_68, species_22, species_72, species_7, species_14, species_66, species_25, species_60, species_6, species_30, species_75, species_36, species_50, species_19, species_5, species_64, species_33, species_69, species_39, species_55, species_24, species_48, species_62, species_53, species_41, species_51, species_54, species_18, species_67, species_57, species_63, species_23, species_56, species_37, species_42, species_32, species_2, species_46, species_65, species_31, species_29, species_40, species_38, species_8, species_34, species_9, species_15, species_11, species_70= u 
	volMedia, SHC1_total__FREE, kp_SE, kon_EGF__, km_GE, km_dim_L_L, kon_EGF, ratio_kphos_receiver, kp_EGF, kphos_f, kp_GE, kon, c1, ratio_kpkd_YN__FREE, onrate_f, kon__, ratio_kpkd_Y1068__FREE, Kd_EGF, km_SE, EGFR_total, ratio_kphosY1173, offrate_f, kp_dim_L_L_pre, numCells, Vextra, EGFconc, GRB2_total__FREE, kdephosY1068_f, f, Kd_EGF__, kdephosYN__FREE, NA, GRB2_f, km_EGF, Kd_GE, km_dim_L_L_pre, kdephosY1068__FREE, kdephosY1173_f, ratio_kdephosY1173, KD_dim, Vc, kp_dim_L_L, _rateLaw1, Kd_SE = p 
	unphosR = species_2 + species_5 + species_6 + species_7
	pY1068 = species_8 + species_11 + species_14 + species_17 + species_18 + species_20 + species_21 + species_23 + species_24 + species_27 + species_30 + species_31 + species_33 + species_35 + species_36 + species_37 + species_38 + species_40 + species_41 + species_42 + species_43 + species_45 + species_46 + species_48 + species_50 + species_51 + species_52 + species_53 + species_55 + species_56 + species_57 + species_58 + species_59 + species_60 + species_61 + species_62 + species_63 + species_64 + species_66 + species_67 + species_68 + species_69 + species_70 + species_71 + species_72 + species_73 + species_74 + species_75
	EGFRtot = species_2 + species_5 + species_6 + species_7 + species_8 + species_9 + species_10 + species_11 + species_12 + species_13 + species_14 + species_15 + species_16 + species_17 + species_18 + species_19 + species_20 + species_21 + species_22 + species_23 + species_24 + species_25 + species_26 + species_27 + species_28 + species_29 + species_30 + species_31 + species_32 + species_33 + species_34 + species_35 + species_36 + species_37 + species_38 + species_39 + species_40 + species_41 + species_42 + species_43 + species_44 + species_45 + species_46 + species_47 + species_48 + species_49 + species_50 + species_51 + species_52 + species_53 + species_54 + species_55 + species_56 + species_57 + species_58 + species_59 + species_60 + species_61 + species_62 + species_63 + species_64 + species_65 + species_66 + species_67 + species_68 + species_69 + species_70 + species_71 + species_72 + species_73 + species_74 + species_75
	pY1173 = species_9 + species_12 + species_15 + species_17 + species_19 + species_20 + species_22 + species_25 + species_26 + species_28 + species_30 + species_32 + species_34 + species_35 + species_36 + species_37 + species_39 + species_40 + species_41 + species_42 + species_44 + species_45 + species_47 + species_49 + species_50 + species_51 + species_52 + species_54 + species_55 + species_56 + species_57 + species_58 + species_59 + species_60 + species_61 + species_62 + species_63 + species_65 + species_66 + species_67 + species_68 + species_69 + species_70 + species_71 + species_72 + species_73 + species_74 + species_75

	if observableId == "observable_pY1068_percent" 
		return nothing
	end

	if observableId == "observable_pY1173_percent" 
		return nothing
	end

	if observableId == "observable_unphosR_per" 
		return nothing
	end

	if observableId == "observable_phosR_per" 
		return nothing
	end

end

function evalDYmodDp(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	species_61, species_10, species_21, species_47, species_17, species_1, species_43, species_52, species_45, species_20, species_49, species_12, species_35, species_16, species_26, species_73, species_59, species_13, species_27, species_71, species_74, species_58, species_4, species_28, species_44, species_3, species_68, species_22, species_72, species_7, species_14, species_66, species_25, species_60, species_6, species_30, species_75, species_36, species_50, species_19, species_5, species_64, species_33, species_69, species_39, species_55, species_24, species_48, species_62, species_53, species_41, species_51, species_54, species_18, species_67, species_57, species_63, species_23, species_56, species_37, species_42, species_32, species_2, species_46, species_65, species_31, species_29, species_40, species_38, species_8, species_34, species_9, species_15, species_11, species_70= u 
	volMedia, SHC1_total__FREE, kp_SE, kon_EGF__, km_GE, km_dim_L_L, kon_EGF, ratio_kphos_receiver, kp_EGF, kphos_f, kp_GE, kon, c1, ratio_kpkd_YN__FREE, onrate_f, kon__, ratio_kpkd_Y1068__FREE, Kd_EGF, km_SE, EGFR_total, ratio_kphosY1173, offrate_f, kp_dim_L_L_pre, numCells, Vextra, EGFconc, GRB2_total__FREE, kdephosY1068_f, f, Kd_EGF__, kdephosYN__FREE, NA, GRB2_f, km_EGF, Kd_GE, km_dim_L_L_pre, kdephosY1068__FREE, kdephosY1173_f, ratio_kdephosY1173, KD_dim, Vc, kp_dim_L_L, _rateLaw1, Kd_SE = p 
	unphosR = species_2 + species_5 + species_6 + species_7
	pY1068 = species_8 + species_11 + species_14 + species_17 + species_18 + species_20 + species_21 + species_23 + species_24 + species_27 + species_30 + species_31 + species_33 + species_35 + species_36 + species_37 + species_38 + species_40 + species_41 + species_42 + species_43 + species_45 + species_46 + species_48 + species_50 + species_51 + species_52 + species_53 + species_55 + species_56 + species_57 + species_58 + species_59 + species_60 + species_61 + species_62 + species_63 + species_64 + species_66 + species_67 + species_68 + species_69 + species_70 + species_71 + species_72 + species_73 + species_74 + species_75
	EGFRtot = species_2 + species_5 + species_6 + species_7 + species_8 + species_9 + species_10 + species_11 + species_12 + species_13 + species_14 + species_15 + species_16 + species_17 + species_18 + species_19 + species_20 + species_21 + species_22 + species_23 + species_24 + species_25 + species_26 + species_27 + species_28 + species_29 + species_30 + species_31 + species_32 + species_33 + species_34 + species_35 + species_36 + species_37 + species_38 + species_39 + species_40 + species_41 + species_42 + species_43 + species_44 + species_45 + species_46 + species_47 + species_48 + species_49 + species_50 + species_51 + species_52 + species_53 + species_54 + species_55 + species_56 + species_57 + species_58 + species_59 + species_60 + species_61 + species_62 + species_63 + species_64 + species_65 + species_66 + species_67 + species_68 + species_69 + species_70 + species_71 + species_72 + species_73 + species_74 + species_75
	pY1173 = species_9 + species_12 + species_15 + species_17 + species_19 + species_20 + species_22 + species_25 + species_26 + species_28 + species_30 + species_32 + species_34 + species_35 + species_36 + species_37 + species_39 + species_40 + species_41 + species_42 + species_44 + species_45 + species_47 + species_49 + species_50 + species_51 + species_52 + species_54 + species_55 + species_56 + species_57 + species_58 + species_59 + species_60 + species_61 + species_62 + species_63 + species_65 + species_66 + species_67 + species_68 + species_69 + species_70 + species_71 + species_72 + species_73 + species_74 + species_75

	if observableId == "observable_pY1068_percent" 
		return nothing
	end

	if observableId == "observable_pY1173_percent" 
		return nothing
	end

	if observableId == "observable_unphosR_per" 
		return nothing
	end

	if observableId == "observable_phosR_per" 
		return nothing
	end

end

function evalDSdDu!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	species_61, species_10, species_21, species_47, species_17, species_1, species_43, species_52, species_45, species_20, species_49, species_12, species_35, species_16, species_26, species_73, species_59, species_13, species_27, species_71, species_74, species_58, species_4, species_28, species_44, species_3, species_68, species_22, species_72, species_7, species_14, species_66, species_25, species_60, species_6, species_30, species_75, species_36, species_50, species_19, species_5, species_64, species_33, species_69, species_39, species_55, species_24, species_48, species_62, species_53, species_41, species_51, species_54, species_18, species_67, species_57, species_63, species_23, species_56, species_37, species_42, species_32, species_2, species_46, species_65, species_31, species_29, species_40, species_38, species_8, species_34, species_9, species_15, species_11, species_70= u 
	volMedia, SHC1_total__FREE, kp_SE, kon_EGF__, km_GE, km_dim_L_L, kon_EGF, ratio_kphos_receiver, kp_EGF, kphos_f, kp_GE, kon, c1, ratio_kpkd_YN__FREE, onrate_f, kon__, ratio_kpkd_Y1068__FREE, Kd_EGF, km_SE, EGFR_total, ratio_kphosY1173, offrate_f, kp_dim_L_L_pre, numCells, Vextra, EGFconc, GRB2_total__FREE, kdephosY1068_f, f, Kd_EGF__, kdephosYN__FREE, NA, GRB2_f, km_EGF, Kd_GE, km_dim_L_L_pre, kdephosY1068__FREE, kdephosY1173_f, ratio_kdephosY1173, KD_dim, Vc, kp_dim_L_L, _rateLaw1, Kd_SE = p 
	if observableId == "observable_pY1068_percent" 
		return nothing
	end

	if observableId == "observable_pY1173_percent" 
		return nothing
	end

	if observableId == "observable_unphosR_per" 
		return nothing
	end

	if observableId == "observable_phosR_per" 
		return nothing
	end

end

function evalDSdDp!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	species_61, species_10, species_21, species_47, species_17, species_1, species_43, species_52, species_45, species_20, species_49, species_12, species_35, species_16, species_26, species_73, species_59, species_13, species_27, species_71, species_74, species_58, species_4, species_28, species_44, species_3, species_68, species_22, species_72, species_7, species_14, species_66, species_25, species_60, species_6, species_30, species_75, species_36, species_50, species_19, species_5, species_64, species_33, species_69, species_39, species_55, species_24, species_48, species_62, species_53, species_41, species_51, species_54, species_18, species_67, species_57, species_63, species_23, species_56, species_37, species_42, species_32, species_2, species_46, species_65, species_31, species_29, species_40, species_38, species_8, species_34, species_9, species_15, species_11, species_70= u 
	volMedia, SHC1_total__FREE, kp_SE, kon_EGF__, km_GE, km_dim_L_L, kon_EGF, ratio_kphos_receiver, kp_EGF, kphos_f, kp_GE, kon, c1, ratio_kpkd_YN__FREE, onrate_f, kon__, ratio_kpkd_Y1068__FREE, Kd_EGF, km_SE, EGFR_total, ratio_kphosY1173, offrate_f, kp_dim_L_L_pre, numCells, Vextra, EGFconc, GRB2_total__FREE, kdephosY1068_f, f, Kd_EGF__, kdephosYN__FREE, NA, GRB2_f, km_EGF, Kd_GE, km_dim_L_L_pre, kdephosY1068__FREE, kdephosY1173_f, ratio_kdephosY1173, KD_dim, Vc, kp_dim_L_L, _rateLaw1, Kd_SE = p 
	if observableId == "observable_pY1068_percent" 
		return nothing
	end

	if observableId == "observable_pY1173_percent" 
		return nothing
	end

	if observableId == "observable_unphosR_per" 
		return nothing
	end

	if observableId == "observable_phosR_per" 
		return nothing
	end

end

