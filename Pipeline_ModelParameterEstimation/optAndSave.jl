
function optModelSaveResults(model, p, doLogSearch, option1, option2, option3, iStartPar, iterations, methodFile, write)

    benchRunTime = Vector{Float64}(undef, iterations)
    success = true
    p_opt = Vector{Float64}(undef, length(p))
    cost_opt = 0.0 

    p_save = copy(p)

    try
        JuMP.optimize!(model)
        cost_opt = objective_value(model)
        if cost_opt == 0.0 && occursin("error", lowercase(termination_status(model)))
            p_opt .= NaN
            success = false
        else
            p_opt[:] = [value(p[i]) for i=1:length(p)]
            view(p_opt, doLogSearch) .= exp10.(view(p_opt, doLogSearch))
            success = true 
        end
        
    catch err
        println(err)
        success = false
        p_opt .= NaN
        cost_opt = Inf
    end

    terminationStatus = termination_status(model)
    primalStatus = primal_status(model)

    if success
        for iter in 1:iterations
            p = copy(p_save)
            benchRunTime[iter] = @elapsed JuMP.optimize!(model) # in seconds
        
        end
    else
        benchRunTime .= NaN
    end

    data = DataFrame(method = methodFile[1:end-3], option1 = option1, option2 = option2, option3 = option3, startParameterIndex = iStartPar, 
            success = success, terminationStatus = terminationStatus, primalStatus = primalStatus, cost = cost_opt, 
            runTime = benchRunTime, iteration = 1:iterations)
    if isfile(write)
        CSV.write(write, data, append = true)
    else
        CSV.write(write, data)
    end

    data_pars = DataFrame(p_opt', :auto)
    insertcols!(data_pars, 1, :method => methodFile[1:end-3])
    insertcols!(data_pars, 2, :option1 => option1)
    insertcols!(data_pars, 3, :option2 => option2)
    insertcols!(data_pars, 4, :option3 => option3)
    insertcols!(data_pars, 5, :startParameterIndex => iStartPar)
    insertcols!(data_pars, 6, :iteration => 1)
    write_pars = write[1:end-4] * "_pars.csv"
    if isfile(write_pars)
        CSV.write(write_pars, data_pars, append = true)
    else
        CSV.write(write_pars, data_pars)
    end
            
end


function optAdamSaveResults(step, adam_opt, doLogSearch, n_it, b2, stepRange, iStartPar, iterations, methodFile, write)
    loss_val = zeros(Float64, n_it)
    theta_save = copy(adam_opt.theta)

    benchRunTime = Vector{Float64}(undef, iterations)
    benchMemory = Vector{Float64}(undef, iterations)
    benchAllocs = Vector{Float64}(undef, iterations)
    success = Vector{Bool}(undef, iterations)
    success .= true
    p_opt = Vector{Float64}(undef, length(theta_save))
    cost_opt = Vector{Float64}(undef, iterations)
    
    for iter in 1:iterations 
        adam_opt.theta = copy(theta_save)
        adam_opt.loss = 0.0
        adam_opt.m = zeros(length(theta_save))
        adam_opt.v = zeros(length(theta_save))
        adam_opt.t = 0
        adam_opt.β = 1.0
        adam_opt.fail = 0
        benchRunTime[iter] = @elapsed begin # in seconds
            for i in 1:n_it
                step()
                loss_val[i] = adam_opt.loss
    
                if adam_opt.fail >= 10
                    success[iter] = false
                    break
                end
            end
        end

        cost_opt[iter] = adam_opt.loss
        p_opt[:] = adam_opt.theta
        view(p_opt, doLogSearch) .= exp10.(view(p_opt, doLogSearch))

        data_pars = DataFrame(p_opt', :auto)
        insertcols!(data_pars, 1, :method => methodFile[1:end-3])
        insertcols!(data_pars, 2, :option1 => n_it)
        insertcols!(data_pars, 3, :option2 => b2)
        insertcols!(data_pars, 4, :option3 => string(stepRange))
        insertcols!(data_pars, 5, :startParameterIndex => iStartPar)
        insertcols!(data_pars, 6, :iteration => iter)
        write_pars = write[1:end-4] * "_pars.csv"
        if isfile(write_pars)
            CSV.write(write_pars, data_pars, append = true)
        else
            CSV.write(write_pars, data_pars)
        end

        data_cost = DataFrame(loss_val', :auto)
        insertcols!(data_cost, 1, :method => methodFile[1:end-3])
        insertcols!(data_cost, 2, :option1 => n_it)
        insertcols!(data_cost, 3, :option2 => b2)
        insertcols!(data_cost, 4, :option3 => string(stepRange))
        insertcols!(data_cost, 5, :startParameterIndex => iStartPar)
        insertcols!(data_cost, 6, :iteration => iter)
        write_cost = write[1:end-4] * "_forwGrad_cost.csv"
        if isfile(write_cost)
            CSV.write(write_cost, data_cost, append = true)
        else
            CSV.write(write_cost, data_cost)
        end

    end

    data = DataFrame(method = methodFile[1:end-3], option1 = n_it, option2 = b2, option3 = string(stepRange), startParameterIndex = iStartPar, 
            success = success, terminationStatus = "-", primalStatus = "-", cost = cost_opt,
            runTime = benchRunTime, iteration = 1:iterations)
    if isfile(write)
        CSV.write(write, data, append = true)
    else
        CSV.write(write, data)
    end

end