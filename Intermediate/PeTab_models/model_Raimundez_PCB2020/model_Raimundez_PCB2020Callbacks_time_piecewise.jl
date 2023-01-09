function getCallbacks_model_Raimundez_PCB2020()
	cb_CET_bool1 = DiscreteCallback(condition_CET_bool1, affect_CET_bool1!, save_positions=(false, false))

	cb_EGF_bool1 = DiscreteCallback(condition_EGF_bool1, affect_EGF_bool1!, save_positions=(false, false))

	cb_METinh_bool1 = DiscreteCallback(condition_METinh_bool1, affect_METinh_bool1!, save_positions=(false, false))

	return CallbackSet(cb_CET_bool1, cb_EGF_bool1, cb_METinh_bool1), [activeAtTime0_CET_bool1!, activeAtTime0_EGF_bool1!, activeAtTime0_METinh_bool1!]
end


function condition_CET_bool1(u, t, integrator)
	t == 0
end
function affect_CET_bool1!(integrator)
	integrator.p[4] = 1.0
end
function activeAtTime0_CET_bool1!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[4] = 0.0 # Default to being off
	if !(t < 0)
		p[4] = 1.0
	end
end


function condition_EGF_bool1(u, t, integrator)
	t == 0
end
function affect_EGF_bool1!(integrator)
	integrator.p[5] = 1.0
end
function activeAtTime0_EGF_bool1!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[5] = 0.0 # Default to being off
	if !(t < 0)
		p[5] = 1.0
	end
end


function condition_METinh_bool1(u, t, integrator)
	t == 0
end
function affect_METinh_bool1!(integrator)
	integrator.p[24] = 1.0
end
function activeAtTime0_METinh_bool1!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[24] = 0.0 # Default to being off
	if !(t < 0)
		p[24] = 1.0
	end
end

function getTstops(u, p)
	 return [dualToFloat(-0.0), dualToFloat(-0.0), dualToFloat(-0.0)]
end