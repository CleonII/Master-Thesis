function getCallbacks_model_Oliveira_NatCommun2021()
	cb_beta_bool2 = ContinuousCallback(condition_beta_bool2, affect_beta_bool2!, save_positions=(false, false))

	cb_beta_bool1 = ContinuousCallback(condition_beta_bool1, affect_beta_bool1!, save_positions=(false, false))

	return CallbackSet(cb_beta_bool2, cb_beta_bool1), [isActiveAtTime0_beta_bool2!, isActiveAtTime0_beta_bool1!]
end


function condition_beta_bool2(u, t, integrator)
	t - integrator.p[11]
end

function affect_beta_bool2!(integrator)
	integrator.p[21] = 1.0
end

function isActiveAtTime0_beta_bool2!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[21] = 0.0 # Default to being off
	if !(t < p[11])
		p[21] = 1.0
	end
end



function condition_beta_bool1(u, t, integrator)
	t - integrator.p[3]
end

function affect_beta_bool1!(integrator)
	integrator.p[17] = 1.0
end

function isActiveAtTime0_beta_bool1!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[17] = 0.0 # Default to being off
	if !(t < p[3])
		p[17] = 1.0
	end
end


function computeTstops(u::AbstractVector, p::AbstractVector)
	 return Float64[]
end