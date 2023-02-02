function getCallbacks_model_Isensee_JCB2018()
	cb_H89_bool1 = DiscreteCallback(condition_H89_bool1, affect_H89_bool1!, save_positions=(false, false))

	cb_Sp8_Br_cAMPS_AM_bool1 = DiscreteCallback(condition_Sp8_Br_cAMPS_AM_bool1, affect_Sp8_Br_cAMPS_AM_bool1!, save_positions=(false, false))

	cb_Rp_cAMPS_pAB_bool1 = DiscreteCallback(condition_Rp_cAMPS_pAB_bool1, affect_Rp_cAMPS_pAB_bool1!, save_positions=(false, false))

	cb_Fsk_bool1 = DiscreteCallback(condition_Fsk_bool1, affect_Fsk_bool1!, save_positions=(false, false))

	cb_IBMXex_bool1 = DiscreteCallback(condition_IBMXex_bool1, affect_IBMXex_bool1!, save_positions=(false, false))

	cb_fourABnOH_bool1 = DiscreteCallback(condition_fourABnOH_bool1, affect_fourABnOH_bool1!, save_positions=(false, false))

	cb_Rp8_Br_cAMPS_pAB_bool1 = DiscreteCallback(condition_Rp8_Br_cAMPS_pAB_bool1, affect_Rp8_Br_cAMPS_pAB_bool1!, save_positions=(false, false))

	cb_Rp8_pCPT_cAMPS_pAB_bool1 = DiscreteCallback(condition_Rp8_pCPT_cAMPS_pAB_bool1, affect_Rp8_pCPT_cAMPS_pAB_bool1!, save_positions=(false, false))

	return CallbackSet(cb_H89_bool1, cb_Sp8_Br_cAMPS_AM_bool1, cb_Rp_cAMPS_pAB_bool1, cb_Fsk_bool1, cb_IBMXex_bool1, cb_fourABnOH_bool1, cb_Rp8_Br_cAMPS_pAB_bool1, cb_Rp8_pCPT_cAMPS_pAB_bool1), [isActiveAtTime0_H89_bool1!, isActiveAtTime0_Sp8_Br_cAMPS_AM_bool1!, isActiveAtTime0_Rp_cAMPS_pAB_bool1!, isActiveAtTime0_Fsk_bool1!, isActiveAtTime0_IBMXex_bool1!, isActiveAtTime0_fourABnOH_bool1!, isActiveAtTime0_Rp8_Br_cAMPS_pAB_bool1!, isActiveAtTime0_Rp8_pCPT_cAMPS_pAB_bool1!]
end


function condition_H89_bool1(u, t, integrator)
	t - integrator.p[28] == 0
end

function affect_H89_bool1!(integrator)
	integrator.p[16] = 1.0
end

function isActiveAtTime0_H89_bool1!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[16] = 0.0 # Default to being off
	if !(t - p[28] < 0)
		p[16] = 1.0
	end
end



function condition_Sp8_Br_cAMPS_AM_bool1(u, t, integrator)
	t - integrator.p[46] == 0
end

function affect_Sp8_Br_cAMPS_AM_bool1!(integrator)
	integrator.p[56] = 1.0
end

function isActiveAtTime0_Sp8_Br_cAMPS_AM_bool1!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[56] = 0.0 # Default to being off
	if !(t - p[46] < 0)
		p[56] = 1.0
	end
end



function condition_Rp_cAMPS_pAB_bool1(u, t, integrator)
	t - integrator.p[17] == 0
end

function affect_Rp_cAMPS_pAB_bool1!(integrator)
	integrator.p[29] = 1.0
end

function isActiveAtTime0_Rp_cAMPS_pAB_bool1!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[29] = 0.0 # Default to being off
	if !(t - p[17] < 0)
		p[29] = 1.0
	end
end



function condition_Fsk_bool1(u, t, integrator)
	t - integrator.p[38] == 0
end

function affect_Fsk_bool1!(integrator)
	integrator.p[70] = 1.0
end

function isActiveAtTime0_Fsk_bool1!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[70] = 0.0 # Default to being off
	if !(t - p[38] < 0)
		p[70] = 1.0
	end
end



function condition_IBMXex_bool1(u, t, integrator)
	t - integrator.p[15] == 0
end

function affect_IBMXex_bool1!(integrator)
	integrator.p[37] = 1.0
end

function isActiveAtTime0_IBMXex_bool1!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[37] = 0.0 # Default to being off
	if !(t - p[15] < 0)
		p[37] = 1.0
	end
end



function condition_fourABnOH_bool1(u, t, integrator)
	t - integrator.p[5] == 0
end

function affect_fourABnOH_bool1!(integrator)
	integrator.p[65] = 1.0
end

function isActiveAtTime0_fourABnOH_bool1!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[65] = 0.0 # Default to being off
	if !(t - p[5] < 0)
		p[65] = 1.0
	end
end



function condition_Rp8_Br_cAMPS_pAB_bool1(u, t, integrator)
	t - integrator.p[51] == 0
end

function affect_Rp8_Br_cAMPS_pAB_bool1!(integrator)
	integrator.p[41] = 1.0
end

function isActiveAtTime0_Rp8_Br_cAMPS_pAB_bool1!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[41] = 0.0 # Default to being off
	if !(t - p[51] < 0)
		p[41] = 1.0
	end
end



function condition_Rp8_pCPT_cAMPS_pAB_bool1(u, t, integrator)
	t - integrator.p[39] == 0
end

function affect_Rp8_pCPT_cAMPS_pAB_bool1!(integrator)
	integrator.p[25] = 1.0
end

function isActiveAtTime0_Rp8_pCPT_cAMPS_pAB_bool1!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[25] = 0.0 # Default to being off
	if !(t - p[39] < 0)
		p[25] = 1.0
	end
end


function computeTstops(u::AbstractVector, p::AbstractVector)
	 return Float64[dualToFloat(p[28]), dualToFloat(p[46]), dualToFloat(p[17]), dualToFloat(p[38]), dualToFloat(p[15]), dualToFloat(p[5]), dualToFloat(p[51]), dualToFloat(p[39])]
end