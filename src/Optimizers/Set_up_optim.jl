"""
    createOptimProb(peTabOpt::PeTabOpt,
                    optimAlg;
                    hessianUse::Symbol=:blockAutoDiff)

    For a PeTab model optimization struct (peTabOpt) create an Optim optmization (evalOptim) 
    function using as optimAlg IPNewton (interior point Newton) or LBFGS, BFGS, ConjugateGradient. 

    For IPNewton the hessian is computed via eiter autoDiff (:autoDiff), or approximated 
    with blockAutoDiff (:blockAutoDiff). All optimizer struct can take their default 
    arguments, for example, LBFGS(linesearch = LineSearches.HagerZhang()) is a valid 
    argument for LBFGS.

    To run the optimization just call evalOptim(paramVec, showTrace=true/false)
"""
function createOptimProb(peTabOpt::PeTabOpt,
                         optimAlg;
                         hessianUse::Symbol=:blockAutoDiff, 
                         options=Optim.Options(iterations = 1000, 
                                               show_trace = false, 
                                               allow_f_increases=true, 
                                               successive_f_tol = 3, 
                                               f_tol=1e-8, 
                                               g_tol=1e-6, 
                                               x_tol=0.0))
    
    if typeof(optimAlg) <: IPNewton
        return createOptimInteriorNewton(peTabOpt, hessianUse=hessianUse, options=options)
    elseif typeof(optimAlg) <: LBFGS || typeof(optimAlg) <: BFGS || typeof(optimAlg) <: ConjugateGradient
        return createOptimFminbox(peTabOpt, optimAlg, options)
    else
        println("Error : optimAlg $optimAlg is not supported")
        println("Supported methods are IPNewton, ConjugateGradient, LBFGS and BFGS")
    end

end


"""
    createOptimInteriorNewton(peTabOpt::PeTabOpt;
                              hessianUse::Symbol=:blockAutoDiff)

    For a PeTab model optimization struct (peTabOpt) create an Optim interior point Newton 
    function struct where the hessian is computed via eiter autoDiff (:autoDiff), or approximated 
    with blockAutoDiff (:blockAutoDiff). 
"""
function createOptimInteriorNewton(peTabOpt::PeTabOpt;
                                   hessianUse::Symbol=:blockAutoDiff, 
                                   options)
                                   
    lowerBounds = peTabOpt.lowerBounds
    upperBounds = peTabOpt.upperBounds

    nParam = length(lowerBounds)
    if hessianUse == :autoDiff
        evalHessian = peTabOpt.evalHess
    elseif hessianUse == :blockAutoDiff
        evalHessian = peTabOpt.evalHessApprox
    else
        println("Error : For optim interior point Newton availble hessianUse options are :autoDiff, :blockAutoDiff not $hessianUse")
    end

    x0 = zeros(Float64, nParam)
    df = TwiceDifferentiable(peTabOpt.evalF, peTabOpt.evalGradF, evalHessian, x0)
    dfc = TwiceDifferentiableConstraints(lowerBounds, upperBounds)

    evalOptim = (p0; showTrace=false) -> begin 
                                               # Move points within bounds 
                                               iBelow = p0 .<= peTabOpt.lowerBounds
                                               iAbove = p0 .>= peTabOpt.upperBounds 
                                               p0[iBelow] .= peTabOpt.lowerBounds[iBelow] .+ 0.001
                                               p0[iAbove] .= peTabOpt.upperBounds[iAbove] .- 0.001
                                               df.f(p0)
                                               println("p0 = $p0")
                                               return Optim.optimize(df, 
                                                                     dfc, 
                                                                     p0, 
                                                                     IPNewton(), 
                                                                     options)
                                         end

    return evalOptim
end


"""
    createOptimFminbox(peTabOpt::PeTabOpt;
                       lineSearch=LineSearches.HagerZhang())

    For a PeTab model optimization struct (peTabOpt) create an Optim Fminbox optimizer where the 
    inner optimizer is either LBFGS or BFGS using lineSearch. 
"""
function createOptimFminbox(peTabOpt::PeTabOpt, 
                            optimAlg, 
                            options)

    lowerBounds = peTabOpt.lowerBounds
    upperBounds = peTabOpt.upperBounds

    evalOptim = (p0; showTrace=false) -> Optim.optimize(peTabOpt.evalF, 
                                                        peTabOpt.evalGradF, 
                                                        lowerBounds, 
                                                        upperBounds, 
                                                        p0, 
                                                        Fminbox(optimAlg), 
                                                        options)
    return evalOptim
end