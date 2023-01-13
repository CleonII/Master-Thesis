function getCallbacks_model_Chen_MSB2009()
	cb_Step1_bool3 = DiscreteCallback(condition_Step1_bool3, affect_Step1_bool3!, save_positions=(false, false))

	cb_StepMini_bool2 = DiscreteCallback(condition_StepMini_bool2, affect_StepMini_bool2!, save_positions=(false, false))

	cb_Step33_bool1 = DiscreteCallback(condition_Step33_bool1, affect_Step33_bool1!, save_positions=(false, false))

	cb_Step1_bool2 = DiscreteCallback(condition_Step1_bool2, affect_Step1_bool2!, save_positions=(false, false))

	cb_Step33_bool3 = DiscreteCallback(condition_Step33_bool3, affect_Step33_bool3!, save_positions=(false, false))

	cb_StepLate_bool1 = DiscreteCallback(condition_StepLate_bool1, affect_StepLate_bool1!, save_positions=(false, false))

	cb_Step1_bool1 = DiscreteCallback(condition_Step1_bool1, affect_Step1_bool1!, save_positions=(false, false))

	cb_Step10Mio_bool2 = DiscreteCallback(condition_Step10Mio_bool2, affect_Step10Mio_bool2!, save_positions=(false, false))

	cb_Step10Mio_bool1 = DiscreteCallback(condition_Step10Mio_bool1, affect_Step10Mio_bool1!, save_positions=(false, false))

	cb_Step33_bool2 = DiscreteCallback(condition_Step33_bool2, affect_Step33_bool2!, save_positions=(false, false))

	cb_StepLate_bool3 = DiscreteCallback(condition_StepLate_bool3, affect_StepLate_bool3!, save_positions=(false, false))

	cb_Step10Mio_bool3 = DiscreteCallback(condition_Step10Mio_bool3, affect_Step10Mio_bool3!, save_positions=(false, false))

	cb_Step5_bool1 = DiscreteCallback(condition_Step5_bool1, affect_Step5_bool1!, save_positions=(false, false))

	cb_Step5_bool3 = DiscreteCallback(condition_Step5_bool3, affect_Step5_bool3!, save_positions=(false, false))

	cb_StepMini_bool1 = DiscreteCallback(condition_StepMini_bool1, affect_StepMini_bool1!, save_positions=(false, false))

	cb_Step5_bool2 = DiscreteCallback(condition_Step5_bool2, affect_Step5_bool2!, save_positions=(false, false))

	cb_StepMini_bool3 = DiscreteCallback(condition_StepMini_bool3, affect_StepMini_bool3!, save_positions=(false, false))

	cb_StepLate_bool2 = DiscreteCallback(condition_StepLate_bool2, affect_StepLate_bool2!, save_positions=(false, false))

	return CallbackSet(cb_Step1_bool3, cb_StepMini_bool2, cb_Step33_bool1, cb_Step1_bool2, cb_Step33_bool3, cb_StepLate_bool1, cb_Step1_bool1, cb_Step10Mio_bool2, cb_Step10Mio_bool1, cb_Step33_bool2, cb_StepLate_bool3, cb_Step10Mio_bool3, cb_Step5_bool1, cb_Step5_bool3, cb_StepMini_bool1, cb_Step5_bool2, cb_StepMini_bool3, cb_StepLate_bool2), [activeAtTime0_Step1_bool3!, activeAtTime0_StepMini_bool2!, activeAtTime0_Step33_bool1!, activeAtTime0_Step1_bool2!, activeAtTime0_Step33_bool3!, activeAtTime0_StepLate_bool1!, activeAtTime0_Step1_bool1!, activeAtTime0_Step10Mio_bool2!, activeAtTime0_Step10Mio_bool1!, activeAtTime0_Step33_bool2!, activeAtTime0_StepLate_bool3!, activeAtTime0_Step10Mio_bool3!, activeAtTime0_Step5_bool1!, activeAtTime0_Step5_bool3!, activeAtTime0_StepMini_bool1!, activeAtTime0_Step5_bool2!, activeAtTime0_StepMini_bool3!, activeAtTime0_StepLate_bool2!]
end


function condition_Step1_bool3(u, t, integrator)
	t - 1800 == 0
end
function affect_Step1_bool3!(integrator)
	integrator.p[102] = 1.0
end
function activeAtTime0_Step1_bool3!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[102] = 0.0 # Default to being off
	if !(t - 1800 < 0)
		p[102] = 1.0
	end
end


function condition_StepMini_bool2(u, t, integrator)
	t - 1799.99 == 0
end
function affect_StepMini_bool2!(integrator)
	integrator.p[135] = 1.0
end
function activeAtTime0_StepMini_bool2!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[135] = 0.0 # Default to being off
	if !(t - 1799.99 < 0)
		p[135] = 1.0
	end
end


function condition_Step33_bool1(u, t, integrator)
	t - 1800 == 0
end
function affect_Step33_bool1!(integrator)
	integrator.p[37] = 1.0
end
function activeAtTime0_Step33_bool1!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[37] = 0.0 # Default to being off
	if !(t - 1800 < 0)
		p[37] = 1.0
	end
end


function condition_Step1_bool2(u, t, integrator)
	t - 1799.99 == 0
end
function affect_Step1_bool2!(integrator)
	integrator.p[194] = 1.0
end
function activeAtTime0_Step1_bool2!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[194] = 0.0 # Default to being off
	if !(t - 1799.99 < 0)
		p[194] = 1.0
	end
end


function condition_Step33_bool3(u, t, integrator)
	t - 1800 == 0
end
function affect_Step33_bool3!(integrator)
	integrator.p[29] = 1.0
end
function activeAtTime0_Step33_bool3!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[29] = 0.0 # Default to being off
	if !(t - 1800 < 0)
		p[29] = 1.0
	end
end


function condition_StepLate_bool1(u, t, integrator)
	t - 2660 == 0
end
function affect_StepLate_bool1!(integrator)
	integrator.p[115] = 1.0
end
function activeAtTime0_StepLate_bool1!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[115] = 0.0 # Default to being off
	if !(t - 2660 < 0)
		p[115] = 1.0
	end
end


function condition_Step1_bool1(u, t, integrator)
	t - 1800 == 0
end
function affect_Step1_bool1!(integrator)
	integrator.p[60] = 1.0
end
function activeAtTime0_Step1_bool1!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[60] = 0.0 # Default to being off
	if !(t - 1800 < 0)
		p[60] = 1.0
	end
end


function condition_Step10Mio_bool2(u, t, integrator)
	t - 1799.99 == 0
end
function affect_Step10Mio_bool2!(integrator)
	integrator.p[31] = 1.0
end
function activeAtTime0_Step10Mio_bool2!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[31] = 0.0 # Default to being off
	if !(t - 1799.99 < 0)
		p[31] = 1.0
	end
end


function condition_Step10Mio_bool1(u, t, integrator)
	t - 1800 == 0
end
function affect_Step10Mio_bool1!(integrator)
	integrator.p[129] = 1.0
end
function activeAtTime0_Step10Mio_bool1!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[129] = 0.0 # Default to being off
	if !(t - 1800 < 0)
		p[129] = 1.0
	end
end


function condition_Step33_bool2(u, t, integrator)
	t - 1799.99 == 0
end
function affect_Step33_bool2!(integrator)
	integrator.p[170] = 1.0
end
function activeAtTime0_Step33_bool2!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[170] = 0.0 # Default to being off
	if !(t - 1799.99 < 0)
		p[170] = 1.0
	end
end


function condition_StepLate_bool3(u, t, integrator)
	t - 2660 == 0
end
function affect_StepLate_bool3!(integrator)
	integrator.p[185] = 1.0
end
function activeAtTime0_StepLate_bool3!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[185] = 0.0 # Default to being off
	if !(t - 2660 < 0)
		p[185] = 1.0
	end
end


function condition_Step10Mio_bool3(u, t, integrator)
	t - 1800 == 0
end
function affect_Step10Mio_bool3!(integrator)
	integrator.p[49] = 1.0
end
function activeAtTime0_Step10Mio_bool3!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[49] = 0.0 # Default to being off
	if !(t - 1800 < 0)
		p[49] = 1.0
	end
end


function condition_Step5_bool1(u, t, integrator)
	t - 1800 == 0
end
function affect_Step5_bool1!(integrator)
	integrator.p[63] = 1.0
end
function activeAtTime0_Step5_bool1!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[63] = 0.0 # Default to being off
	if !(t - 1800 < 0)
		p[63] = 1.0
	end
end


function condition_Step5_bool3(u, t, integrator)
	t - 1800 == 0
end
function affect_Step5_bool3!(integrator)
	integrator.p[188] = 1.0
end
function activeAtTime0_Step5_bool3!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[188] = 0.0 # Default to being off
	if !(t - 1800 < 0)
		p[188] = 1.0
	end
end


function condition_StepMini_bool1(u, t, integrator)
	t - 1800 == 0
end
function affect_StepMini_bool1!(integrator)
	integrator.p[100] = 1.0
end
function activeAtTime0_StepMini_bool1!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[100] = 0.0 # Default to being off
	if !(t - 1800 < 0)
		p[100] = 1.0
	end
end


function condition_Step5_bool2(u, t, integrator)
	t - 1799.99 == 0
end
function affect_Step5_bool2!(integrator)
	integrator.p[84] = 1.0
end
function activeAtTime0_Step5_bool2!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[84] = 0.0 # Default to being off
	if !(t - 1799.99 < 0)
		p[84] = 1.0
	end
end


function condition_StepMini_bool3(u, t, integrator)
	t - 1800 == 0
end
function affect_StepMini_bool3!(integrator)
	integrator.p[23] = 1.0
end
function activeAtTime0_StepMini_bool3!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[23] = 0.0 # Default to being off
	if !(t - 1800 < 0)
		p[23] = 1.0
	end
end


function condition_StepLate_bool2(u, t, integrator)
	t - 2659.99 == 0
end
function affect_StepLate_bool2!(integrator)
	integrator.p[175] = 1.0
end
function activeAtTime0_StepLate_bool2!(u, p)
	t = 0.0 # Used to check conditions activated at t0=0
	p[175] = 0.0 # Default to being off
	if !(t - 2659.99 < 0)
		p[175] = 1.0
	end
end

function getTstops(u, p)
	 return [dualToFloat(1800.0), dualToFloat(1799.99), dualToFloat(1800.0), dualToFloat(1799.99), dualToFloat(1800.0), dualToFloat(2660.0), dualToFloat(1800.0), dualToFloat(1799.99), dualToFloat(1800.0), dualToFloat(1799.99), dualToFloat(2660.0), dualToFloat(1800.0), dualToFloat(1800.0), dualToFloat(1800.0), dualToFloat(1800.0), dualToFloat(1799.99), dualToFloat(1800.0), dualToFloat(2659.99)]
end