function evalDYmodDu(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	IR2, IR2in, Rec2, IR1in, Uptake1, Uptake2, InsulinFragments, IR1, Rec1, Ins, BoundUnspec= u 
	ka1, ini_R2fold, kout, ini_R1, kout_frag, koff_unspec, kin, ka2fold, kin2, kd1, kon_unspec, init_Ins, kd2fold, ExtracellularMedium, kout2 = p 
	if observableId == "observable_IR1" 
		observableParameter1_observable_IR1, observableParameter2_observable_IR1 = getObsOrSdParam(obsPar, mapObsParam)
		out[4] = observableParameter2_observable_IR1
		out[8] = observableParameter2_observable_IR1
		return nothing
	end

	if observableId == "observable_IR2" 
		observableParameter1_observable_IR2, observableParameter2_observable_IR2 = getObsOrSdParam(obsPar, mapObsParam)
		out[1] = observableParameter2_observable_IR2
		out[2] = observableParameter2_observable_IR2
		return nothing
	end

	if observableId == "observable_IRsum" 
		observableParameter1_observable_IRsum, observableParameter2_observable_IRsum = getObsOrSdParam(obsPar, mapObsParam)
		out[1] = 0.395observableParameter2_observable_IRsum
		out[2] = 0.395observableParameter2_observable_IRsum
		out[4] = 0.605observableParameter2_observable_IRsum
		out[8] = 0.605observableParameter2_observable_IRsum
		return nothing
	end

	if observableId == "observable_Insulin" 
		observableParameter1_observable_Insulin, observableParameter2_observable_Insulin, observableParameter3_observable_Insulin, observableParameter4_observable_Insulin = getObsOrSdParam(obsPar, mapObsParam)
		out[7] = (observableParameter2_observable_Insulin*observableParameter3_observable_Insulin*(observableParameter4_observable_Insulin^2)*(((Ins + observableParameter4_observable_Insulin + InsulinFragments*observableParameter3_observable_Insulin) / observableParameter4_observable_Insulin)^2) - observableParameter2_observable_Insulin*observableParameter3_observable_Insulin*(Ins^2) - observableParameter2_observable_Insulin*(InsulinFragments^2)*(observableParameter3_observable_Insulin^3) - 2Ins*InsulinFragments*observableParameter2_observable_Insulin*(observableParameter3_observable_Insulin^2) - Ins*observableParameter2_observable_Insulin*observableParameter3_observable_Insulin*observableParameter4_observable_Insulin - InsulinFragments*observableParameter2_observable_Insulin*observableParameter4_observable_Insulin*(observableParameter3_observable_Insulin^2)) / (observableParameter4_observable_Insulin*(Ins + observableParameter4_observable_Insulin + InsulinFragments*observableParameter3_observable_Insulin)*(((Ins + observableParameter4_observable_Insulin + InsulinFragments*observableParameter3_observable_Insulin) / observableParameter4_observable_Insulin)^2))
		out[10] = (observableParameter2_observable_Insulin*(observableParameter4_observable_Insulin^2)*(((Ins + observableParameter4_observable_Insulin + InsulinFragments*observableParameter3_observable_Insulin) / observableParameter4_observable_Insulin)^2) - observableParameter2_observable_Insulin*(Ins^2) - Ins*observableParameter2_observable_Insulin*observableParameter4_observable_Insulin - observableParameter2_observable_Insulin*(InsulinFragments^2)*(observableParameter3_observable_Insulin^2) - 2Ins*InsulinFragments*observableParameter2_observable_Insulin*observableParameter3_observable_Insulin - InsulinFragments*observableParameter2_observable_Insulin*observableParameter3_observable_Insulin*observableParameter4_observable_Insulin) / (observableParameter4_observable_Insulin*(Ins + observableParameter4_observable_Insulin + InsulinFragments*observableParameter3_observable_Insulin)*(((Ins + observableParameter4_observable_Insulin + InsulinFragments*observableParameter3_observable_Insulin) / observableParameter4_observable_Insulin)^2))
		return nothing
	end

end

function evalDYmodDp(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	IR2, IR2in, Rec2, IR1in, Uptake1, Uptake2, InsulinFragments, IR1, Rec1, Ins, BoundUnspec= u 
	ka1, ini_R2fold, kout, ini_R1, kout_frag, koff_unspec, kin, ka2fold, kin2, kd1, kon_unspec, init_Ins, kd2fold, ExtracellularMedium, kout2 = p 
	if observableId == "observable_IR1" 
		return nothing
	end

	if observableId == "observable_IR2" 
		return nothing
	end

	if observableId == "observable_IRsum" 
		return nothing
	end

	if observableId == "observable_Insulin" 
		return nothing
	end

end

function evalDSdDu!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	IR2, IR2in, Rec2, IR1in, Uptake1, Uptake2, InsulinFragments, IR1, Rec1, Ins, BoundUnspec= u 
	ka1, ini_R2fold, kout, ini_R1, kout_frag, koff_unspec, kin, ka2fold, kin2, kd1, kon_unspec, init_Ins, kd2fold, ExtracellularMedium, kout2 = p 
	if observableId == "observable_IR1" 
		return nothing
	end

	if observableId == "observable_IR2" 
		return nothing
	end

	if observableId == "observable_IRsum" 
		return nothing
	end

	if observableId == "observable_Insulin" 
		return nothing
	end

end

function evalDSdDp!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	IR2, IR2in, Rec2, IR1in, Uptake1, Uptake2, InsulinFragments, IR1, Rec1, Ins, BoundUnspec= u 
	ka1, ini_R2fold, kout, ini_R1, kout_frag, koff_unspec, kin, ka2fold, kin2, kd1, kon_unspec, init_Ins, kd2fold, ExtracellularMedium, kout2 = p 
	if observableId == "observable_IR1" 
		return nothing
	end

	if observableId == "observable_IR2" 
		return nothing
	end

	if observableId == "observable_IRsum" 
		return nothing
	end

	if observableId == "observable_Insulin" 
		return nothing
	end

end

