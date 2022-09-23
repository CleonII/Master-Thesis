function evalYmod(u, t, dynPar, obsPar, paramData, obsData, observableId, simulationId) 

	species_61, species_10, species_21, species_47, species_17, species_1, species_43, species_52, species_45, species_20, species_49, species_12, species_35, species_16, species_26, species_73, species_59, species_13, species_27, species_71, species_74, species_58, species_4, species_28, species_44, species_3, species_68, species_22, species_72, species_7, species_14, species_66, species_25, species_60, species_6, species_30, species_75, species_36, species_50, species_19, species_5, species_64, species_33, species_69, species_39, species_55, species_24, species_48, species_62, species_53, species_41, species_51, species_54, species_18, species_67, species_57, species_63, species_23, species_56, species_37, species_42, species_32, species_2, species_46, species_65, species_31, species_29, species_40, species_38, species_8, species_34, species_9, species_15, species_11, species_70, dummyVariable= u 
	GRB2_total__FREE, SHC1_total__FREE, kdephosY1068__FREE, kdephosYN__FREE, ratio_kpkd_Y1068__FREE, ratio_kpkd_YN__FREE = dynPar 

	EGFRtot = species_2+species_5+species_6+species_7+species_8+species_9+species_10+species_11+species_12+species_13+species_14+species_15+species_16+species_17+species_18+species_19+species_20+species_21+species_22+species_23+species_24+species_25+species_26+species_27+species_28+species_29+species_30+species_31+species_32+species_33+species_34+species_35+species_36+species_37+species_38+species_39+species_40+species_41+species_42+species_43+species_44+species_45+species_46+species_47+species_48+species_49+species_50+species_51+species_52+species_53+species_54+species_55+species_56+species_57+species_58+species_59+species_60+species_61+species_62+species_63+species_64+species_65+species_66+species_67+species_68+species_69+species_70+species_71+species_72+species_73+species_74+species_75
	pY1068 = species_8+species_11+species_14+species_17+species_18+species_20+species_21+species_23+species_24+species_27+species_30+species_31+species_33+species_35+species_36+species_37+species_38+species_40+species_41+species_42+species_43+species_45+species_46+species_48+species_50+species_51+species_52+species_53+species_55+species_56+species_57+species_58+species_59+species_60+species_61+species_62+species_63+species_64+species_66+species_67+species_68+species_69+species_70+species_71+species_72+species_73+species_74+species_75
	pY1173 = species_9+species_12+species_15+species_17+species_19+species_20+species_22+species_25+species_26+species_28+species_30+species_32+species_34+species_35+species_36+species_37+species_39+species_40+species_41+species_42+species_44+species_45+species_47+species_49+species_50+species_51+species_52+species_54+species_55+species_56+species_57+species_58+species_59+species_60+species_61+species_62+species_63+species_65+species_66+species_67+species_68+species_69+species_70+species_71+species_72+species_73+species_74+species_75
	unphosR = species_2+species_5+species_6+species_7

	if observableId == "observable_pY1068_percent" 
		return ( 100 *  pY1068) /  EGFRtot
	end

	if observableId == "observable_pY1173_percent" 
		return ( 100 * pY1173 ) /  EGFRtot
	end

	if observableId == "observable_unphosR_per" 
		return ( 100 *  unphosR) /  EGFRtot
	end

	if observableId == "observable_phosR_per" 
		return (100 * (EGFRtot - unphosR)) / EGFRtot
	end

end

function evalU0!(u0Vec, paramVec) 

	km_GE, kdephosYN__FREE, ratio_kdephosY1173, km_SE, km_EGF, c1, kp_SE, kp_EGF, kdephosY1173_f, kdephosY1068_f, kp_GE, kdephosY1068__FREE, ratio_kpkd_Y1068__FREE, ratio_kpkd_YN__FREE, ratio_kphosY1173, km_dim_L_L, _rateLaw1, kphos_f, ratio_kphos_receiver, Kd_EGF, Kd_SE, numCells, EGFconc, Vextra, kp_dim_L_L, SHC1_total__FREE, Vc, kon, kp_dim_L_L_pre, kon__, onrate_f, GRB2_f, kon_EGF__, EGFR_total, offrate_f, Kd_EGF__, NA, KD_dim, Kd_GE, f, km_dim_L_L_pre, volMedia, GRB2_total__FREE, kon_EGF = paramVec 

	species_61 = 0.0 
	species_10 = 0.0 
	species_21 = 0.0 
	species_47 = 0.0 
	species_17 = 0.0 
	species_1 = EGFconc * NA * Vextra 
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
	species_4 = SHC1_total__FREE * f 
	species_28 = 0.0 
	species_44 = 0.0 
	species_3 = GRB2_f * GRB2_total__FREE * f 
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
	species_2 = EGFR_total 
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
	dummyVariable = 0.0 

	u0Vec .= species_61, species_10, species_21, species_47, species_17, species_1, species_43, species_52, species_45, species_20, species_49, species_12, species_35, species_16, species_26, species_73, species_59, species_13, species_27, species_71, species_74, species_58, species_4, species_28, species_44, species_3, species_68, species_22, species_72, species_7, species_14, species_66, species_25, species_60, species_6, species_30, species_75, species_36, species_50, species_19, species_5, species_64, species_33, species_69, species_39, species_55, species_24, species_48, species_62, species_53, species_41, species_51, species_54, species_18, species_67, species_57, species_63, species_23, species_56, species_37, species_42, species_32, species_2, species_46, species_65, species_31, species_29, species_40, species_38, species_8, species_34, species_9, species_15, species_11, species_70, dummyVariable
end

function evalSd!(u, t, sdPar, dynPar, paramData, obsData, observableId, simulationId) 

	species_61, species_10, species_21, species_47, species_17, species_1, species_43, species_52, species_45, species_20, species_49, species_12, species_35, species_16, species_26, species_73, species_59, species_13, species_27, species_71, species_74, species_58, species_4, species_28, species_44, species_3, species_68, species_22, species_72, species_7, species_14, species_66, species_25, species_60, species_6, species_30, species_75, species_36, species_50, species_19, species_5, species_64, species_33, species_69, species_39, species_55, species_24, species_48, species_62, species_53, species_41, species_51, species_54, species_18, species_67, species_57, species_63, species_23, species_56, species_37, species_42, species_32, species_2, species_46, species_65, species_31, species_29, species_40, species_38, species_8, species_34, species_9, species_15, species_11, species_70, dummyVariable= u 
	GRB2_total__FREE, SHC1_total__FREE, kdephosY1068__FREE, kdephosYN__FREE, ratio_kpkd_Y1068__FREE, ratio_kpkd_YN__FREE = dynPar 

	if observableId == "observable_pY1068_percent" 
		noiseParameter1_observable_pY1068_percent = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_pY1068_percent 
	end

	if observableId == "observable_pY1173_percent" 
		noiseParameter1_observable_pY1173_percent = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_pY1173_percent 
	end

	if observableId == "observable_unphosR_per" 
		noiseParameter1_observable_unphosR_per = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_unphosR_per 
	end

	if observableId == "observable_phosR_per" 
		noiseParameter1_observable_phosR_per = getObsOrSdParam(sdPar, paramData, obsData, observableId, simulationId, t, getObsPar=false)
		return noiseParameter1_observable_phosR_per 
	end

end