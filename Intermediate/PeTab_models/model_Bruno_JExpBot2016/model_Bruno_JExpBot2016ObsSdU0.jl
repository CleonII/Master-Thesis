function evalYmod(u, t, dynPar, obsPar, nonDynParam, paramData, observableId, mapObsParam) 

	b10, bio, ohbio, zea, bcry, ohb10, bcar, dummyVariable= u 
	k5, kb1, kb2, kc1, kc2, kc4 = dynPar 
	init_b10_1, init_bcar1, init_bcar2, init_bcry_1, init_ohb10_1, init_zea_1, szea = nonDynParam 

	if observableId == "ob10" 
		return b10 
	end

	if observableId == "obcar" 
		return bcar 
	end

	if observableId == "obcry" 
		return bcry 
	end

	if observableId == "obio" 
		return bio 
	end

	if observableId == "oohb10" 
		return ohb10 
	end

	if observableId == "ozea" 
		return zea 
	end

end

function evalU0!(u0Vec, paramVec) 

	kb1, kb1_multiplier, kc1, kc1_multiplier, kb2_multiplier, kb2, kc2, kc2_multiplier, k5_multiplier, kc4_multiplier, k5, kc4, init_ohb10, init_bcry, init_zea, init_bcar, init_b10 = paramVec 

	b10 = init_b10 
	bio = 0.0 
	ohbio = 0.0 
	zea = init_zea 
	bcry = init_bcry 
	ohb10 = init_ohb10 
	bcar = init_bcar 
	dummyVariable = 0.0 

	u0Vec .= b10, bio, ohbio, zea, bcry, ohb10, bcar, dummyVariable
end

function evalSd!(u, t, sdPar, dynPar, nonDynParam, paramData, observableId, mapSdParam) 

	b10, bio, ohbio, zea, bcry, ohb10, bcar, dummyVariable= u 
	k5, kb1, kb2, kc1, kc2, kc4 = dynPar 
	init_b10_1, init_bcar1, init_bcar2, init_bcry_1, init_ohb10_1, init_zea_1, szea = nonDynParam 

	if observableId == "ob10" 
		noiseParameter1_ob10 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_ob10 
	end

	if observableId == "obcar" 
		noiseParameter1_obcar = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_obcar 
	end

	if observableId == "obcry" 
		noiseParameter1_obcry = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_obcry 
	end

	if observableId == "obio" 
		noiseParameter1_obio = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_obio 
	end

	if observableId == "oohb10" 
		noiseParameter1_oohb10 = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_oohb10 
	end

	if observableId == "ozea" 
		noiseParameter1_ozea = getObsOrSdParam(sdPar, mapSdParam)
		return noiseParameter1_ozea 
	end

end