function getCallbacks_model_Bertozzi_PNAS2020()
	cb_Ro_bool2 = DiscreteCallback(condition_Ro_bool2, affect_Ro_bool2!, save_positions=(false, false))

	cb_Ro_bool3 = DiscreteCallback(condition_Ro_bool3, affect_Ro_bool3!, save_positions=(false, false))

	cb_Ro_bool4 = DiscreteCallback(condition_Ro_bool4, affect_Ro_bool4!, save_positions=(false, false))

	cb_Ro_bool1 = DiscreteCallback(condition_Ro_bool1, affect_Ro_bool1!, save_positions=(false, false))

	return CallbackSet(cb_Ro_bool2, cb_Ro_bool3, cb_Ro_bool4, cb_Ro_bool1), [isActiveAtTime0_Ro_bool2!, isActiveAtTime0_Ro_bool3!, isActiveAtTime0_Ro_bool4!, isActiveAtTime0_Ro_bool1!]
end


function condition_Ro_bool2(u, t, integrator)
	t == integrator.p[10] * integrator.p[1] + (1 - integrator.p[10]) * integrator.p[16]
end

function affect_Ro_bool2!(integrator)
	integrator.p[4] = 1.0
end

function isActiveAtTime0_Ro_bool2!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[4] = 0.0 # Default to being off
	if !(t ≤ p[10] * p[1] + (1 - p[10]) * p[16])
		p[4] = 1.0
	end
end



function condition_Ro_bool3(u, t, integrator)
	t == integrator.p[10] * integrator.p[19] + (1 - integrator.p[10]) * integrator.p[12]
end

function affect_Ro_bool3!(integrator)
	integrator.p[6] = 1.0
end

function isActiveAtTime0_Ro_bool3!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[6] = 0.0 # Default to being off
	if (t > p[10] * p[19] + (1 - p[10]) * p[12])
		p[6] = 1.0
	end
end



function condition_Ro_bool4(u, t, integrator)
	t == integrator.p[10] * integrator.p[1] + (1 - integrator.p[10]) * integrator.p[16]
end

function affect_Ro_bool4!(integrator)
	integrator.p[7] = 1.0
end

function isActiveAtTime0_Ro_bool4!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[7] = 0.0 # Default to being off
	if !(t ≤ p[10] * p[1] + (1 - p[10]) * p[16])
		p[7] = 1.0
	end
end



function condition_Ro_bool1(u, t, integrator)
	t == integrator.p[10] * integrator.p[19] + (1 - integrator.p[10]) * integrator.p[12]
end

function affect_Ro_bool1!(integrator)
	integrator.p[3] = 1.0
end

function isActiveAtTime0_Ro_bool1!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[3] = 0.0 # Default to being off
	if (t > p[10] * p[19] + (1 - p[10]) * p[12])
		p[3] = 1.0
	end
end


function computeTstops(u::AbstractVector, p::AbstractVector)
	 return Float64[dualToFloat(p[1]*p[10] + p[16]*(1 - p[10])), dualToFloat(p[19]*p[10] + p[12]*(1 - p[10])), dualToFloat(p[1]*p[10] + p[16]*(1 - p[10])), dualToFloat(p[19]*p[10] + p[12]*(1 - p[10]))]
end