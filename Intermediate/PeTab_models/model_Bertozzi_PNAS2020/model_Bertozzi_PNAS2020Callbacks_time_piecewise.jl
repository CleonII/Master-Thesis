function getCallbacks_model_Bertozzi_PNAS2020()
	cb_Ro_bool2 = DiscreteCallback(condition_Ro_bool2, affect_Ro_bool2!, save_positions=(false, false))

	cb_Ro_bool5 = DiscreteCallback(condition_Ro_bool5, affect_Ro_bool5!, save_positions=(false, false))

	cb_Ro_bool4 = DiscreteCallback(condition_Ro_bool4, affect_Ro_bool4!, save_positions=(false, false))

	cb_Ro_bool1 = DiscreteCallback(condition_Ro_bool1, affect_Ro_bool1!, save_positions=(false, false))

	return CallbackSet(cb_Ro_bool2, cb_Ro_bool5, cb_Ro_bool4, cb_Ro_bool1), [activeAtTime0_Ro_bool2!, activeAtTime0_Ro_bool5!, activeAtTime0_Ro_bool4!, activeAtTime0_Ro_bool1!]
end


function condition_Ro_bool2(u, t, integrator)
	t == integrator.p[9] * integrator.p[1] + (1 - integrator.p[9]) * integrator.p[15]
end
function affect_Ro_bool2!(integrator)
	integrator.p[4] = 1.0
end
function activeAtTime0_Ro_bool2!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[4] = 0.0 # Default to being off
	if !(t ≤ p[9] * p[1] + (1 - p[9]) * p[15])
		p[4] = 1.0
	end
end


function condition_Ro_bool5(u, t, integrator)
	t == integrator.p[9] * integrator.p[1] + (1 - integrator.p[9]) * integrator.p[15]
end
function affect_Ro_bool5!(integrator)
	integrator.p[16] = 1.0
end
function activeAtTime0_Ro_bool5!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[16] = 0.0 # Default to being off
	if !(t ≤ p[9] * p[1] + (1 - p[9]) * p[15])
		p[16] = 1.0
	end
end


function condition_Ro_bool4(u, t, integrator)
	t == integrator.p[9] * integrator.p[19] + (1 - integrator.p[9]) * integrator.p[11]
end
function affect_Ro_bool4!(integrator)
	integrator.p[6] = 1.0
end
function activeAtTime0_Ro_bool4!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[6] = 0.0 # Default to being off
	if (t > p[9] * p[19] + (1 - p[9]) * p[11])
		p[6] = 1.0
	end
end


function condition_Ro_bool1(u, t, integrator)
	t == integrator.p[9] * integrator.p[19] + (1 - integrator.p[9]) * integrator.p[11]
end
function affect_Ro_bool1!(integrator)
	integrator.p[3] = 1.0
end
function activeAtTime0_Ro_bool1!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[3] = 0.0 # Default to being off
	if (t > p[9] * p[19] + (1 - p[9]) * p[11])
		p[3] = 1.0
	end
end

function getTstops(u, p)
	 return [dualToFloat(p[1]*p[9] + p[15]*(1 - p[9])), dualToFloat(p[1]*p[9] + p[15]*(1 - p[9])), dualToFloat(p[19]*p[9] + p[11]*(1 - p[9])), dualToFloat(p[19]*p[9] + p[11]*(1 - p[9]))]
end