function getCallbacks_Smith_BMCSystBiol2013()
	cb_event1 = DiscreteCallback(condition_event1, affect_event1!, save_positions=(false, false))

	cb_event3 = DiscreteCallback(condition_event3, affect_event3!, save_positions=(false, false))

	cb_event2 = DiscreteCallback(condition_event2, affect_event2!, save_positions=(false, false))

	return CallbackSet(cb_event1, cb_event3, cb_event2), []
end


function condition_event1(u, t, integrator)
	t == integrator.p[68]
end

function affect_event1!(integrator)
	integrator.u[103] = 0
end


function condition_event3(u, t, integrator)
	t == 2895
end

function affect_event3!(integrator)
	integrator.u[103] = 0
end


function condition_event2(u, t, integrator)
	t == 2880
end

function affect_event2!(integrator)
	integrator.u[103] = 499999
end

function computeTstops(u::AbstractVector, p::AbstractVector)
	 return Float64[dualToFloat(p[68]), dualToFloat(2895.0), dualToFloat(2880.0)]
end