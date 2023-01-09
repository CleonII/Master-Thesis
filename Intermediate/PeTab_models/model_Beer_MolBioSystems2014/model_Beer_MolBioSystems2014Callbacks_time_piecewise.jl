function getCallbacks_model_Beer_MolBioSystems2014()
	cb_lag_bool1 = ContinuousCallback(condition_lag_bool1, affect_lag_bool1!, save_positions=(false, false))

	return CallbackSet(cb_lag_bool1), [activeAtTime0_lag_bool1!]
end


function condition_lag_bool1(u, t, integrator)
	t - integrator.p[7] - 0
end
function affect_lag_bool1!(integrator)
	integrator.p[1] = 1.0
end
function activeAtTime0_lag_bool1!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[1] = 0.0 # Default to being off
	if !(t - p[7] < 0)
		p[1] = 1.0
	end
end

function getTstops(u, p)
	 return [dualToFloat(p[7])]
end