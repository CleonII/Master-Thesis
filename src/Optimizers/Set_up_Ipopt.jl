# Evaluate hessian for Ipopt. Notice, the calc_hess function requires input on the 
# format (hessian, inputVec). 
function eval_h(x_arg::Vector{Float64}, 
                rows::Vector{Int32}, 
                cols::Vector{Int32}, 
                obj_factor::Float64, 
                lambda::Vector{Float64}, 
                values::Union{Nothing,Vector{Float64}}, 
                n_param, 
                calc_hess::Function)

    if values === nothing
        # Symmetric matrix, fill the lower left triangle only
        idx = 1
        for row in 1:n_param
            for col in 1:row
                rows[idx] = row
                cols[idx] = col
                idx += 1
            end
        end
    else
        # Again, only lower left triangle
        # Objective
        hessian_mat = zeros(n_param, n_param)
        calc_hess(hessian_mat, x_arg)
        idx = 1
        for row in 1:n_param
            for col in 1:row
                values[idx] = hessian_mat[row, col] * obj_factor
                idx += 1
            end
        end
    end
    return
end
# In case user does not want to employ a hessian 
function eval_h_empty(x_arg::Vector{Float64}, 
                      rows::Vector{Int32}, 
                      cols::Vector{Int32}, 
                      obj_factor::Float64, 
                      lambda::Vector{Float64}, 
                      values::Union{Nothing,Vector{Float64}})    
    return nothing
end
function eval_jac_g(x::Vector{Float64}, rows::Vector{Int32}, cols::Vector{Int32}, values::Union{Nothing,Vector{Float64}})
    return 
end
function eval_g(x::Vector{Float64}, g::Vector{Float64})
    return 
end
# Callback for Ipopt allowing the number of iterations to be extracted from the solver
function intermediate(alg_mod::Cint,
                      iter_count::Cint,
                      obj_value::Float64,
                      inf_pr::Float64,
                      inf_du::Float64,
                      mu::Float64,
                      d_norm::Float64,
                      regularization_size::Float64,
                      alpha_du::Float64,
                      alpha_pr::Float64,
                      ls_trials::Cint, 
                      iterArr)
    iterArr[1] = Int(iter_count)

    return true 
end


# Function will change when a PeTabOpt struct is created.
function createIpoptProbNew(fPre::Function, 
                            fGradPre::Function, 
                            fHessPre::Function, 
                            lowerBounds::Array{<:AbstractFloat, 1}, 
                            upperBounds::Array{<:AbstractFloat, 1}; 
                            emptyH::Bool=false)

    nParam = length(lowerBounds)
    if emptyH == false
        evalHUse = (x_arg, rows, cols, obj_factor, lambda, values) -> eval_h(x_arg, rows, cols, obj_factor, lambda, values, nParam, fHessPre)
    else
        evalHUse = eval_h_empty
    end

    iterArr = ones(Int64, 1) .* 20
    intermediateUse = (alg_mod, iter_count, obj_value, inf_pr, inf_du, mu, d_norm, regularization_size, alpha_du, alpha_pr, ls_trials) -> intermediate(alg_mod, iter_count, obj_value, inf_pr, inf_du, mu, d_norm, regularization_size, alpha_du, alpha_pr, ls_trials, iterArr)

    m = 0
    nParamHess = Int(nParam*(nParam + 1) / 2)
    # No inequality constraints assumed 
    g_L = Float64[]
    g_U = Float64[] 
    prob = Ipopt.CreateIpoptProblem(nParam, lowerBounds, upperBounds, m, g_L, g_U, 0, nParamHess, fPre, eval_g, fGradPre, eval_jac_g, evalHUse)
    Ipopt.SetIntermediateCallback(prob, intermediateUse) # Allow iterations to be retrevied 

    if emptyH == true
        Ipopt.AddIpoptStrOption(prob, "hessian_approximation", "limited-memory")
    end 
    Ipopt.AddIpoptIntOption(prob, "print_level", 0)
    Ipopt.AddIpoptIntOption(prob, "max_iter", 1000)

    return prob, iterArr
end

