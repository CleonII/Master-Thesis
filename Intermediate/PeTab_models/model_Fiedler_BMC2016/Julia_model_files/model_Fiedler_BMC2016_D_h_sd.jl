#u[1] = RAF, u[2] = MEK, u[3] = pMEK, u[4] = pERK, u[5] = pRAF, u[6] = ERK
#pODEProblem[1] = ERK_total, pODEProblem[2] = UO126, pODEProblem[3] = k10, pODEProblem[4] = RAF_total, pODEProblem[5] = K_3, pODEProblem[6] = cyt, pODEProblem[7] = k4, pODEProblem[8] = Sorafenib, pODEProblem[9] = K_2, pODEProblem[10] = k6, pODEProblem[11] = k11, pODEProblem[12] = tau1, pODEProblem[13] = MEK_total, pODEProblem[14] = K_1, pODEProblem[15] = k3, pODEProblem[16] = tau2, pODEProblem[17] = k5, pODEProblem[18] = k2
#
function compute_∂h∂u!(u, t::Real, pODEProblem::AbstractVector, θ_observable::AbstractVector,                    
                       θ_nonDynamic::AbstractVector, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :pErk 
		observableParameter1_pErk = getObsOrSdParam(θ_observable, parameterMap)
		out[4] = observableParameter1_pErk
		return nothing
	end

	if observableId == :pMek 
		observableParameter1_pMek = getObsOrSdParam(θ_observable, parameterMap)
		out[3] = observableParameter1_pMek
		return nothing
	end

end

function compute_∂h∂p!(u, t::Real, pODEProblem::AbstractVector, θ_observable::AbstractVector,
                       θ_nonDynamic::AbstractVector, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :pErk 
		return nothing
	end

	if observableId == :pMek 
		return nothing
	end

end

function compute_∂σ∂σu!(u, t::Real, θ_sd::AbstractVector, pODEProblem::AbstractVector, θ_nonDynamic::AbstractVector, 
                        parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :pErk 
		return nothing
	end

	if observableId == :pMek 
		return nothing
	end

end

function compute_∂σ∂σp!(u, t::Real, θ_sd::AbstractVector, pODEProblem::AbstractVector, θ_nonDynamic::AbstractVector, 
                        parameterInfo::ParametersInfo, observableId::Symbol, parameterMap::θObsOrSdParameterMap, out) 
	if observableId == :pErk 
		return nothing
	end

	if observableId == :pMek 
		return nothing
	end

end

