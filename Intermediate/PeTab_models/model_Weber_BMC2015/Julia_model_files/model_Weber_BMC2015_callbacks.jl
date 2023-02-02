function getCallbacks_model_Weber_BMC2015()
	cb_u4_bool1 = DiscreteCallback(condition_u4_bool1, affect_u4_bool1!, save_positions=(false, false))

	cb_u6_bool1 = DiscreteCallback(condition_u6_bool1, affect_u6_bool1!, save_positions=(false, false))

	cb_u5_bool1 = DiscreteCallback(condition_u5_bool1, affect_u5_bool1!, save_positions=(false, false))

	cb_u3_bool1 = DiscreteCallback(condition_u3_bool1, affect_u3_bool1!, save_positions=(false, false))

	return CallbackSet(cb_u4_bool1, cb_u6_bool1, cb_u5_bool1, cb_u3_bool1), [isActiveAtTime0_u4_bool1!, isActiveAtTime0_u6_bool1!, isActiveAtTime0_u5_bool1!, isActiveAtTime0_u3_bool1!]
end


function condition_u4_bool1(u, t, integrator)
	t == 0
end

function affect_u4_bool1!(integrator)
	integrator.p[3] = 1.0
end

function isActiveAtTime0_u4_bool1!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[3] = 0.0 # Default to being off
	if !(t < 0)
		p[3] = 1.0
	end
end



function condition_u6_bool1(u, t, integrator)
	t - integrator.p[7] == 0
end

function affect_u6_bool1!(integrator)
	integrator.p[18] = 1.0
end

function isActiveAtTime0_u6_bool1!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[18] = 0.0 # Default to being off
	if !(t - p[7] < 0)
		p[18] = 1.0
	end
end



function condition_u5_bool1(u, t, integrator)
	t - integrator.p[22] == 0
end

function affect_u5_bool1!(integrator)
	integrator.p[4] = 1.0
end

function isActiveAtTime0_u5_bool1!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[4] = 0.0 # Default to being off
	if !(t - p[22] < 0)
		p[4] = 1.0
	end
end



function condition_u3_bool1(u, t, integrator)
	t == 0
end

function affect_u3_bool1!(integrator)
	integrator.p[38] = 1.0
end

function isActiveAtTime0_u3_bool1!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[38] = 0.0 # Default to being off
	if !(t < 0)
		p[38] = 1.0
	end
end


function computeTstops(u::AbstractVector, p::AbstractVector)
	 return Float64[dualToFloat(-0.0), dualToFloat(p[7]), dualToFloat(p[22]), dualToFloat(-0.0)]
end