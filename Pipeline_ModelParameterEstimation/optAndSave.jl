
function optModelSaveResults(model, p, doLogSearch, option1, option2, option3, option4, iStartPar, methodFile, write)

    benchRunTime = 0.0
    success = true
    p_opt = Vector{Float64}(undef, length(p))
    cost_opt = 0.0 

    println("Starting optimizing...")
    try
        benchRunTime = @elapsed JuMP.optimize!(model)
        cost_opt = objective_value(model)
        if cost_opt == 0.0 && occursin("error", lowercase(termination_status(model)))
            p_opt .= NaN
            benchRunTime = NaN
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
        benchRunTime = NaN
    end
    println("Done optimizing!")

    terminationStatus = termination_status(model)
    primalStatus = primal_status(model)

    data = DataFrame(method = methodFile[1:end-3], option1 = option1, option2 = option2, option3 = option3, option4 = option4, startParameterIndex = iStartPar, 
            success = success, terminationStatus = terminationStatus, primalStatus = primalStatus, cost = cost_opt, 
            runTime = benchRunTime, iteration = 1)
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
    insertcols!(data_pars, 5, :option4 => option4)
    insertcols!(data_pars, 6, :startParameterIndex => iStartPar)
    insertcols!(data_pars, 7, :iteration => 1)
    write_pars = write[1:end-4] * "_pars.csv"
    if isfile(write_pars)
        CSV.write(write_pars, data_pars, append = true)
    else
        CSV.write(write_pars, data_pars)
    end
            
end


function optAdamSaveResults(step, adam_opt, doLogSearch, option1, n_it, b2, stepRange, iStartPar, iterations, methodFile, write)
    loss_val = zeros(Float64, n_it)
    theta_save = copy(adam_opt.theta)

    benchRunTime = Vector{Float64}(undef, iterations)
    success = Vector{Bool}(undef, iterations)
    success .= true
    p_opt = Vector{Float64}(undef, length(theta_save))
    cost_opt = fill(Inf64, iterations)
    
    for iter in 1:iterations 
        adam_opt.theta = copy(theta_save)
        adam_opt.loss = 0.0
        adam_opt.m = zeros(length(theta_save))
        adam_opt.v = zeros(length(theta_save))
        adam_opt.t = 0
        adam_opt.β = 1.0
        adam_opt.fail = 0

        println("Starting optimizing...")
        benchRunTime[iter] = @elapsed begin # in seconds
            for i in 1:n_it
                step()
                loss_val[i] = adam_opt.loss
    
                if adam_opt.fail >= 10
                    success[iter] = false
                    loss_val[i+1:end] .= NaN
                    break
                end

                if adam_opt.loss < cost_opt[iter]
                    cost_opt[iter] = adam_opt.loss
                    p_opt .= adam_opt.theta
                end
            end
        end
        println("Done optimizing!")
        
        view(p_opt, doLogSearch) .= exp10.(view(p_opt, doLogSearch))

        data_pars = DataFrame(p_opt', :auto)
        insertcols!(data_pars, 1, :method => methodFile[1:end-3])
        insertcols!(data_pars, 2, :option1 => option1)
        insertcols!(data_pars, 3, :option2 => n_it)
        insertcols!(data_pars, 4, :option3 => b2)
        insertcols!(data_pars, 5, :option4 => string(stepRange))
        insertcols!(data_pars, 6, :startParameterIndex => iStartPar)
        insertcols!(data_pars, 7, :iteration => iter)
        write_pars = write[1:end-4] * "_pars.csv"
        if isfile(write_pars)
            CSV.write(write_pars, data_pars, append = true)
        else
            CSV.write(write_pars, data_pars)
        end

        data_cost = DataFrame(loss_val', :auto)
        insertcols!(data_cost, 1, :method => methodFile[1:end-3])
        insertcols!(data_cost, 2, :option1 => option1)
        insertcols!(data_cost, 3, :option2 => n_it)
        insertcols!(data_cost, 4, :option3 => b2)
        insertcols!(data_cost, 5, :option4 => string(stepRange))
        insertcols!(data_cost, 6, :startParameterIndex => iStartPar)
        insertcols!(data_cost, 7, :iteration => iter)
        write_cost = write[1:end-4] * "_forwGrad_cost.csv"
        if isfile(write_cost)
            CSV.write(write_cost, data_cost, append = true)
        else
            CSV.write(write_cost, data_cost)
        end

    end

    data = DataFrame(method = methodFile[1:end-3], option1 = option1, option2 = n_it, option3 = b2, option4 = string(stepRange), startParameterIndex = iStartPar, 
            success = success, terminationStatus = "-", primalStatus = "-", cost = cost_opt,
            runTime = benchRunTime, iteration = 1:iterations)
    if isfile(write)
        CSV.write(write, data, append = true)
    else
        CSV.write(write, data)
    end

end