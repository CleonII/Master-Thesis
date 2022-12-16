function evalDYmodDu(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	b10, bio, ohbio, zea, bcry, ohb10, bcar= u 
	kc2_multiplier, init_zea, kc4_multiplier, cyt, k5_multiplier, kc1_multiplier, init_b10, init_bcry, kb1_multiplier, kb2_multiplier, kc1, kc4, init_ohb10, init_bcar, kc2, kb2, k5, kb1 = p 
	if observableId == "ob10" 
		out[1] = 1
		return nothing
	end

	if observableId == "obcar" 
		out[7] = 1
		return nothing
	end

	if observableId == "obcry" 
		out[5] = 1
		return nothing
	end

	if observableId == "obio" 
		out[2] = 1
		return nothing
	end

	if observableId == "oohb10" 
		out[6] = 1
		return nothing
	end

	if observableId == "ozea" 
		out[4] = 1
		return nothing
	end

end

function evalDYmodDp(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	b10, bio, ohbio, zea, bcry, ohb10, bcar= u 
	kc2_multiplier, init_zea, kc4_multiplier, cyt, k5_multiplier, kc1_multiplier, init_b10, init_bcry, kb1_multiplier, kb2_multiplier, kc1, kc4, init_ohb10, init_bcar, kc2, kb2, k5, kb1 = p 
	if observableId == "ob10" 
		return nothing
	end

	if observableId == "obcar" 
		return nothing
	end

	if observableId == "obcry" 
		return nothing
	end

	if observableId == "obio" 
		return nothing
	end

	if observableId == "oohb10" 
		return nothing
	end

	if observableId == "ozea" 
		return nothing
	end

end

function evalDSdDu!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	b10, bio, ohbio, zea, bcry, ohb10, bcar= u 
	kc2_multiplier, init_zea, kc4_multiplier, cyt, k5_multiplier, kc1_multiplier, init_b10, init_bcry, kb1_multiplier, kb2_multiplier, kc1, kc4, init_ohb10, init_bcar, kc2, kb2, k5, kb1 = p 
	if observableId == "ob10" 
		return nothing
	end

	if observableId == "obcar" 
		return nothing
	end

	if observableId == "obcry" 
		return nothing
	end

	if observableId == "obio" 
		return nothing
	end

	if observableId == "oohb10" 
		return nothing
	end

	if observableId == "ozea" 
		return nothing
	end

end

function evalDSdDp!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	b10, bio, ohbio, zea, bcry, ohb10, bcar= u 
	kc2_multiplier, init_zea, kc4_multiplier, cyt, k5_multiplier, kc1_multiplier, init_b10, init_bcry, kb1_multiplier, kb2_multiplier, kc1, kc4, init_ohb10, init_bcar, kc2, kb2, k5, kb1 = p 
	if observableId == "ob10" 
		return nothing
	end

	if observableId == "obcar" 
		return nothing
	end

	if observableId == "obcry" 
		return nothing
	end

	if observableId == "obio" 
		return nothing
	end

	if observableId == "oohb10" 
		return nothing
	end

	if observableId == "ozea" 
		return nothing
	end

end

