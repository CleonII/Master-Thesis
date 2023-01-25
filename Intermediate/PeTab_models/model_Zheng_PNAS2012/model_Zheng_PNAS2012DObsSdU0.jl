function evalDYmodDu(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	K27me0K36me0, K27me2K36me3, K27me2K36me0, K27me0K36me1, K27me2K36me1, K27me0K36me2, K27me1K36me2, K27me3K36me2, K27me1K36me1, K27me1K36me3, K27me2K36me2, K27me0K36me3, K27me3K36me1, K27me3K36me0, K27me1K36me0= u 
	k20_10, k11_12, k22_32, k21_11, k13_12, k01_00, k11_01, k22_21, k22_12, k11_10, k12_11, k21_31, k02_12, k30_20, dilution, k31_21, k23_22, k02_01, k03_13, k10_11, k21_20, k10_00, k30_31, k20_30, k20_21, k00_10, k12_13, k01_11, k02_03, k00_01, k03_02, k23_13, k32_31, default, k13_03, k31_30, k12_02, k31_32, k01_02, k13_23, k21_22, k12_22, k11_21, k22_23, k10_20, inflowp, k32_22 = p 
	if observableId == :observable_K27me0K36me0 
		out[1] = (K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[2] = (-K27me0K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[3] = (-K27me0K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[4] = (-K27me0K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[5] = (-K27me0K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[6] = (-K27me0K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[7] = (-K27me0K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[8] = (-K27me0K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[9] = (-K27me0K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[10] = (-K27me0K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[11] = (-K27me0K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[12] = (-K27me0K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[13] = (-K27me0K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[14] = (-K27me0K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[15] = (-K27me0K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		return nothing
	end

	if observableId == :observable_K27me0K36me1 
		out[1] = (-K27me0K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[2] = (-K27me0K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[3] = (-K27me0K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[4] = (K27me0K36me0 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[5] = (-K27me0K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[6] = (-K27me0K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[7] = (-K27me0K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[8] = (-K27me0K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[9] = (-K27me0K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[10] = (-K27me0K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[11] = (-K27me0K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[12] = (-K27me0K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[13] = (-K27me0K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[14] = (-K27me0K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[15] = (-K27me0K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		return nothing
	end

	if observableId == :observable_K27me0K36me2 
		out[1] = (-K27me0K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[2] = (-K27me0K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[3] = (-K27me0K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[4] = (-K27me0K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[5] = (-K27me0K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[6] = (K27me0K36me0 + K27me0K36me1 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[7] = (-K27me0K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[8] = (-K27me0K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[9] = (-K27me0K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[10] = (-K27me0K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[11] = (-K27me0K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[12] = (-K27me0K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[13] = (-K27me0K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[14] = (-K27me0K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[15] = (-K27me0K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		return nothing
	end

	if observableId == :observable_K27me0K36me3 
		out[1] = (-K27me0K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[2] = (-K27me0K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[3] = (-K27me0K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[4] = (-K27me0K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[5] = (-K27me0K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[6] = (-K27me0K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[7] = (-K27me0K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[8] = (-K27me0K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[9] = (-K27me0K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[10] = (-K27me0K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[11] = (-K27me0K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[12] = (K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[13] = (-K27me0K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[14] = (-K27me0K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[15] = (-K27me0K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		return nothing
	end

	if observableId == :observable_K27me1K36me0 
		out[1] = (-K27me1K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[2] = (-K27me1K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[3] = (-K27me1K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[4] = (-K27me1K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[5] = (-K27me1K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[6] = (-K27me1K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[7] = (-K27me1K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[8] = (-K27me1K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[9] = (-K27me1K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[10] = (-K27me1K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[11] = (-K27me1K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[12] = (-K27me1K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[13] = (-K27me1K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[14] = (-K27me1K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[15] = (K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		return nothing
	end

	if observableId == :observable_K27me1K36me1 
		out[1] = (-K27me1K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[2] = (-K27me1K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[3] = (-K27me1K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[4] = (-K27me1K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[5] = (-K27me1K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[6] = (-K27me1K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[7] = (-K27me1K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[8] = (-K27me1K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[9] = (K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[10] = (-K27me1K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[11] = (-K27me1K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[12] = (-K27me1K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[13] = (-K27me1K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[14] = (-K27me1K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[15] = (-K27me1K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		return nothing
	end

	if observableId == :observable_K27me1K36me2 
		out[1] = (-K27me1K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[2] = (-K27me1K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[3] = (-K27me1K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[4] = (-K27me1K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[5] = (-K27me1K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[6] = (-K27me1K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[7] = (K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[8] = (-K27me1K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[9] = (-K27me1K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[10] = (-K27me1K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[11] = (-K27me1K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[12] = (-K27me1K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[13] = (-K27me1K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[14] = (-K27me1K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[15] = (-K27me1K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		return nothing
	end

	if observableId == :observable_K27me1K36me3 
		out[1] = (-K27me1K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[2] = (-K27me1K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[3] = (-K27me1K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[4] = (-K27me1K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[5] = (-K27me1K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[6] = (-K27me1K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[7] = (-K27me1K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[8] = (-K27me1K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[9] = (-K27me1K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[10] = (K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[11] = (-K27me1K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[12] = (-K27me1K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[13] = (-K27me1K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[14] = (-K27me1K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[15] = (-K27me1K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		return nothing
	end

	if observableId == :observable_K27me2K36me0 
		out[1] = (-K27me2K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[2] = (-K27me2K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[3] = (K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[4] = (-K27me2K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[5] = (-K27me2K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[6] = (-K27me2K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[7] = (-K27me2K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[8] = (-K27me2K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[9] = (-K27me2K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[10] = (-K27me2K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[11] = (-K27me2K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[12] = (-K27me2K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[13] = (-K27me2K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[14] = (-K27me2K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[15] = (-K27me2K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		return nothing
	end

	if observableId == :observable_K27me2K36me1 
		out[1] = (-K27me2K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[2] = (-K27me2K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[3] = (-K27me2K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[4] = (-K27me2K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[5] = (K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[6] = (-K27me2K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[7] = (-K27me2K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[8] = (-K27me2K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[9] = (-K27me2K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[10] = (-K27me2K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[11] = (-K27me2K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[12] = (-K27me2K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[13] = (-K27me2K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[14] = (-K27me2K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[15] = (-K27me2K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		return nothing
	end

	if observableId == :observable_K27me2K36me2 
		out[1] = (-K27me2K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[2] = (-K27me2K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[3] = (-K27me2K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[4] = (-K27me2K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[5] = (-K27me2K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[6] = (-K27me2K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[7] = (-K27me2K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[8] = (-K27me2K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[9] = (-K27me2K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[10] = (-K27me2K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[11] = (K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[12] = (-K27me2K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[13] = (-K27me2K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[14] = (-K27me2K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[15] = (-K27me2K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		return nothing
	end

	if observableId == :observable_K27me2K36me3 
		out[1] = (-K27me2K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[2] = (K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[3] = (-K27me2K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[4] = (-K27me2K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[5] = (-K27me2K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[6] = (-K27me2K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[7] = (-K27me2K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[8] = (-K27me2K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[9] = (-K27me2K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[10] = (-K27me2K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[11] = (-K27me2K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[12] = (-K27me2K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[13] = (-K27me2K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[14] = (-K27me2K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[15] = (-K27me2K36me3) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		return nothing
	end

	if observableId == :observable_K27me3K36me0 
		out[1] = (-K27me3K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[2] = (-K27me3K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[3] = (-K27me3K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[4] = (-K27me3K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[5] = (-K27me3K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[6] = (-K27me3K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[7] = (-K27me3K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[8] = (-K27me3K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[9] = (-K27me3K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[10] = (-K27me3K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[11] = (-K27me3K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[12] = (-K27me3K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[13] = (-K27me3K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[14] = (K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me1 + K27me3K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[15] = (-K27me3K36me0) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		return nothing
	end

	if observableId == :observable_K27me3K36me1 
		out[1] = (-K27me3K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[2] = (-K27me3K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[3] = (-K27me3K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[4] = (-K27me3K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[5] = (-K27me3K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[6] = (-K27me3K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[7] = (-K27me3K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[8] = (-K27me3K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[9] = (-K27me3K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[10] = (-K27me3K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[11] = (-K27me3K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[12] = (-K27me3K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[13] = (K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[14] = (-K27me3K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[15] = (-K27me3K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		return nothing
	end

	if observableId == :observable_K27me3K36me2 
		out[1] = (-K27me3K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[2] = (-K27me3K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[3] = (-K27me3K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[4] = (-K27me3K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[5] = (-K27me3K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[6] = (-K27me3K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[7] = (-K27me3K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[8] = (K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[9] = (-K27me3K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[10] = (-K27me3K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[11] = (-K27me3K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[12] = (-K27me3K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[13] = (-K27me3K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[14] = (-K27me3K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		out[15] = (-K27me3K36me2) / ((K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2)^2)
		return nothing
	end

end

function evalDYmodDp(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	K27me0K36me0, K27me2K36me3, K27me2K36me0, K27me0K36me1, K27me2K36me1, K27me0K36me2, K27me1K36me2, K27me3K36me2, K27me1K36me1, K27me1K36me3, K27me2K36me2, K27me0K36me3, K27me3K36me1, K27me3K36me0, K27me1K36me0= u 
	k20_10, k11_12, k22_32, k21_11, k13_12, k01_00, k11_01, k22_21, k22_12, k11_10, k12_11, k21_31, k02_12, k30_20, dilution, k31_21, k23_22, k02_01, k03_13, k10_11, k21_20, k10_00, k30_31, k20_30, k20_21, k00_10, k12_13, k01_11, k02_03, k00_01, k03_02, k23_13, k32_31, default, k13_03, k31_30, k12_02, k31_32, k01_02, k13_23, k21_22, k12_22, k11_21, k22_23, k10_20, inflowp, k32_22 = p 
	if observableId == :observable_K27me0K36me0 
		return nothing
	end

	if observableId == :observable_K27me0K36me1 
		return nothing
	end

	if observableId == :observable_K27me0K36me2 
		return nothing
	end

	if observableId == :observable_K27me0K36me3 
		return nothing
	end

	if observableId == :observable_K27me1K36me0 
		return nothing
	end

	if observableId == :observable_K27me1K36me1 
		return nothing
	end

	if observableId == :observable_K27me1K36me2 
		return nothing
	end

	if observableId == :observable_K27me1K36me3 
		return nothing
	end

	if observableId == :observable_K27me2K36me0 
		return nothing
	end

	if observableId == :observable_K27me2K36me1 
		return nothing
	end

	if observableId == :observable_K27me2K36me2 
		return nothing
	end

	if observableId == :observable_K27me2K36me3 
		return nothing
	end

	if observableId == :observable_K27me3K36me0 
		return nothing
	end

	if observableId == :observable_K27me3K36me1 
		return nothing
	end

	if observableId == :observable_K27me3K36me2 
		return nothing
	end

end

function evalDSdDu!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	K27me0K36me0, K27me2K36me3, K27me2K36me0, K27me0K36me1, K27me2K36me1, K27me0K36me2, K27me1K36me2, K27me3K36me2, K27me1K36me1, K27me1K36me3, K27me2K36me2, K27me0K36me3, K27me3K36me1, K27me3K36me0, K27me1K36me0= u 
	k20_10, k11_12, k22_32, k21_11, k13_12, k01_00, k11_01, k22_21, k22_12, k11_10, k12_11, k21_31, k02_12, k30_20, dilution, k31_21, k23_22, k02_01, k03_13, k10_11, k21_20, k10_00, k30_31, k20_30, k20_21, k00_10, k12_13, k01_11, k02_03, k00_01, k03_02, k23_13, k32_31, default, k13_03, k31_30, k12_02, k31_32, k01_02, k13_23, k21_22, k12_22, k11_21, k22_23, k10_20, inflowp, k32_22 = p 
	if observableId == "observable_K27me0K36me0" 
		return nothing
	end

	if observableId == "observable_K27me0K36me1" 
		return nothing
	end

	if observableId == "observable_K27me0K36me2" 
		return nothing
	end

	if observableId == "observable_K27me0K36me3" 
		return nothing
	end

	if observableId == "observable_K27me1K36me0" 
		return nothing
	end

	if observableId == "observable_K27me1K36me1" 
		return nothing
	end

	if observableId == "observable_K27me1K36me2" 
		return nothing
	end

	if observableId == "observable_K27me1K36me3" 
		return nothing
	end

	if observableId == "observable_K27me2K36me0" 
		return nothing
	end

	if observableId == "observable_K27me2K36me1" 
		return nothing
	end

	if observableId == "observable_K27me2K36me2" 
		return nothing
	end

	if observableId == "observable_K27me2K36me3" 
		return nothing
	end

	if observableId == "observable_K27me3K36me0" 
		return nothing
	end

	if observableId == "observable_K27me3K36me1" 
		return nothing
	end

	if observableId == "observable_K27me3K36me2" 
		return nothing
	end

end

function evalDSdDp!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	K27me0K36me0, K27me2K36me3, K27me2K36me0, K27me0K36me1, K27me2K36me1, K27me0K36me2, K27me1K36me2, K27me3K36me2, K27me1K36me1, K27me1K36me3, K27me2K36me2, K27me0K36me3, K27me3K36me1, K27me3K36me0, K27me1K36me0= u 
	k20_10, k11_12, k22_32, k21_11, k13_12, k01_00, k11_01, k22_21, k22_12, k11_10, k12_11, k21_31, k02_12, k30_20, dilution, k31_21, k23_22, k02_01, k03_13, k10_11, k21_20, k10_00, k30_31, k20_30, k20_21, k00_10, k12_13, k01_11, k02_03, k00_01, k03_02, k23_13, k32_31, default, k13_03, k31_30, k12_02, k31_32, k01_02, k13_23, k21_22, k12_22, k11_21, k22_23, k10_20, inflowp, k32_22 = p 
	if observableId == "observable_K27me0K36me0" 
		return nothing
	end

	if observableId == "observable_K27me0K36me1" 
		return nothing
	end

	if observableId == "observable_K27me0K36me2" 
		return nothing
	end

	if observableId == "observable_K27me0K36me3" 
		return nothing
	end

	if observableId == "observable_K27me1K36me0" 
		return nothing
	end

	if observableId == "observable_K27me1K36me1" 
		return nothing
	end

	if observableId == "observable_K27me1K36me2" 
		return nothing
	end

	if observableId == "observable_K27me1K36me3" 
		return nothing
	end

	if observableId == "observable_K27me2K36me0" 
		return nothing
	end

	if observableId == "observable_K27me2K36me1" 
		return nothing
	end

	if observableId == "observable_K27me2K36me2" 
		return nothing
	end

	if observableId == "observable_K27me2K36me3" 
		return nothing
	end

	if observableId == "observable_K27me3K36me0" 
		return nothing
	end

	if observableId == "observable_K27me3K36me1" 
		return nothing
	end

	if observableId == "observable_K27me3K36me2" 
		return nothing
	end

end

