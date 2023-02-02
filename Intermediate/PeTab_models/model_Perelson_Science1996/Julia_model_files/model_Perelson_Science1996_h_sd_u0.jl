#u[1] = Vni, u[2] = V, u[3] = Vin, u[4] = Tstar
#θ_dynamicNames[1] = c, θ_dynamicNames[2] = delta
##parameterInfo.nominalValue[1] = NN_C 
#parameterInfo.nominalValue[2] = T0_C 
#parameterInfo.nominalValue[5] = K0_C 


function compute_h(u::AbstractVector, t::Real, θ_dynamic::AbstractVector, θ_observable::AbstractVector,
                   θ_nonDynamic::AbstractVector, parameterInfo::ParametersInfo, observableId::Symbol, 
                      parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :task0_model0_perelson1_V 
		return u[2] 
	end

end

function compute_u0!(u0::AbstractVector, pODEProblem::AbstractVector) 

	#pODEProblem[1] = c, pODEProblem[2] = T0, pODEProblem[3] = default, pODEProblem[4] = K0, pODEProblem[5] = NN, pODEProblem[6] = delta

	Vni = 0.0 
	V = 1.86e6 
	Vin = 1.86e6 
	Tstar = 15061.32075 

	u0 .= Vni, V, Vin, Tstar
end

function compute_u0(pODEProblem::AbstractVector)::AbstractVector 

	#pODEProblem[1] = c, pODEProblem[2] = T0, pODEProblem[3] = default, pODEProblem[4] = K0, pODEProblem[5] = NN, pODEProblem[6] = delta

	Vni = 0.0 
	V = 1.86e6 
	Vin = 1.86e6 
	Tstar = 15061.32075 

	 return [Vni, V, Vin, Tstar]
end

function compute_σ(u::AbstractVector, t::Real, θ_sd::AbstractVector, θ_dynamic::AbstractVector, θ_nonDynamic::AbstractVector, 
                   parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap)::Real 
	if observableId == :task0_model0_perelson1_V 
		noiseParameter1_task0_model0_perelson1_V = getObsOrSdParam(θ_sd, parameterMap)
		return noiseParameter1_task0_model0_perelson1_V 
	end

end