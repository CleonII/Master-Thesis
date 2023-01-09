function evalDYmodDu(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	x_k5k12k16, x_k8, x_k16, x_0ac, x_k12, x_k5k8, x_k5k12, x_k12k16, x_k8k12k16, x_k5, x_k5k16, x_k5k8k12, x_k8k12, x_4ac, x_k8k16, x_k5k8k16= u 
	a_k5_k5k12, a_k8, d, default, a_k12k16_k8k12k16, a_basal, a_k5k12_k5k8k12, a_k8k12k16_4ac, a_k12_k5k12, a_k16_k12k16 = p 
	if observableId == "observable_0ac" 
		out[4] = 1
		return nothing
	end

	if observableId == "observable_4ac" 
		out[14] = 1
		return nothing
	end

	if observableId == "observable_k12" 
		out[5] = 1
		return nothing
	end

	if observableId == "observable_k12k16" 
		out[8] = 1
		return nothing
	end

	if observableId == "observable_k16" 
		out[3] = 1
		return nothing
	end

	if observableId == "observable_k5" 
		out[10] = 1
		return nothing
	end

	if observableId == "observable_k5k12" 
		out[7] = 1
		return nothing
	end

	if observableId == "observable_k5k12k16" 
		out[1] = 1
		return nothing
	end

	if observableId == "observable_k5k16" 
		out[11] = 1
		return nothing
	end

	if observableId == "observable_k5k8" 
		out[6] = 1
		return nothing
	end

	if observableId == "observable_k5k8k12" 
		out[12] = 1
		return nothing
	end

	if observableId == "observable_k5k8k16" 
		out[16] = 1
		return nothing
	end

	if observableId == "observable_k8" 
		out[2] = 1
		return nothing
	end

	if observableId == "observable_k8k12" 
		out[13] = 1
		return nothing
	end

	if observableId == "observable_k8k12k16" 
		out[9] = 1
		return nothing
	end

	if observableId == "observable_k8k16" 
		out[15] = 1
		return nothing
	end

end

function evalDYmodDp(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	x_k5k12k16, x_k8, x_k16, x_0ac, x_k12, x_k5k8, x_k5k12, x_k12k16, x_k8k12k16, x_k5, x_k5k16, x_k5k8k12, x_k8k12, x_4ac, x_k8k16, x_k5k8k16= u 
	a_k5_k5k12, a_k8, d, default, a_k12k16_k8k12k16, a_basal, a_k5k12_k5k8k12, a_k8k12k16_4ac, a_k12_k5k12, a_k16_k12k16 = p 
	if observableId == "observable_0ac" 
		return nothing
	end

	if observableId == "observable_4ac" 
		return nothing
	end

	if observableId == "observable_k12" 
		return nothing
	end

	if observableId == "observable_k12k16" 
		return nothing
	end

	if observableId == "observable_k16" 
		return nothing
	end

	if observableId == "observable_k5" 
		return nothing
	end

	if observableId == "observable_k5k12" 
		return nothing
	end

	if observableId == "observable_k5k12k16" 
		return nothing
	end

	if observableId == "observable_k5k16" 
		return nothing
	end

	if observableId == "observable_k5k8" 
		return nothing
	end

	if observableId == "observable_k5k8k12" 
		return nothing
	end

	if observableId == "observable_k5k8k16" 
		return nothing
	end

	if observableId == "observable_k8" 
		return nothing
	end

	if observableId == "observable_k8k12" 
		return nothing
	end

	if observableId == "observable_k8k12k16" 
		return nothing
	end

	if observableId == "observable_k8k16" 
		return nothing
	end

end

function evalDSdDu!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	x_k5k12k16, x_k8, x_k16, x_0ac, x_k12, x_k5k8, x_k5k12, x_k12k16, x_k8k12k16, x_k5, x_k5k16, x_k5k8k12, x_k8k12, x_4ac, x_k8k16, x_k5k8k16= u 
	a_k5_k5k12, a_k8, d, default, a_k12k16_k8k12k16, a_basal, a_k5k12_k5k8k12, a_k8k12k16_4ac, a_k12_k5k12, a_k16_k12k16 = p 
	if observableId == "observable_0ac" 
		return nothing
	end

	if observableId == "observable_4ac" 
		return nothing
	end

	if observableId == "observable_k12" 
		return nothing
	end

	if observableId == "observable_k12k16" 
		return nothing
	end

	if observableId == "observable_k16" 
		return nothing
	end

	if observableId == "observable_k5" 
		return nothing
	end

	if observableId == "observable_k5k12" 
		return nothing
	end

	if observableId == "observable_k5k12k16" 
		return nothing
	end

	if observableId == "observable_k5k16" 
		return nothing
	end

	if observableId == "observable_k5k8" 
		return nothing
	end

	if observableId == "observable_k5k8k12" 
		return nothing
	end

	if observableId == "observable_k5k8k16" 
		return nothing
	end

	if observableId == "observable_k8" 
		return nothing
	end

	if observableId == "observable_k8k12" 
		return nothing
	end

	if observableId == "observable_k8k12k16" 
		return nothing
	end

	if observableId == "observable_k8k16" 
		return nothing
	end

end

function evalDSdDp!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	x_k5k12k16, x_k8, x_k16, x_0ac, x_k12, x_k5k8, x_k5k12, x_k12k16, x_k8k12k16, x_k5, x_k5k16, x_k5k8k12, x_k8k12, x_4ac, x_k8k16, x_k5k8k16= u 
	a_k5_k5k12, a_k8, d, default, a_k12k16_k8k12k16, a_basal, a_k5k12_k5k8k12, a_k8k12k16_4ac, a_k12_k5k12, a_k16_k12k16 = p 
	if observableId == "observable_0ac" 
		return nothing
	end

	if observableId == "observable_4ac" 
		return nothing
	end

	if observableId == "observable_k12" 
		return nothing
	end

	if observableId == "observable_k12k16" 
		return nothing
	end

	if observableId == "observable_k16" 
		return nothing
	end

	if observableId == "observable_k5" 
		return nothing
	end

	if observableId == "observable_k5k12" 
		return nothing
	end

	if observableId == "observable_k5k12k16" 
		return nothing
	end

	if observableId == "observable_k5k16" 
		return nothing
	end

	if observableId == "observable_k5k8" 
		return nothing
	end

	if observableId == "observable_k5k8k12" 
		return nothing
	end

	if observableId == "observable_k5k8k16" 
		return nothing
	end

	if observableId == "observable_k8" 
		return nothing
	end

	if observableId == "observable_k8k12" 
		return nothing
	end

	if observableId == "observable_k8k12k16" 
		return nothing
	end

	if observableId == "observable_k8k16" 
		return nothing
	end

end

