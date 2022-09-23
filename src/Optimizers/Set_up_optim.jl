"""
    createOptimInteriorNewton(peTabOpt::PeTabOpt;
                              hessianUse::Symbol=:blockAutoDiff)

    For a PeTab model optimization struct (peTabOpt) create an Optim interior point Newton 
    function struct where the hessian is computed via eiter autoDiff (:autoDiff), or approximated 
    with blockAutoDiff (:blockAutoDiff). 

    To run the optimization just call evalOptim(paramVec, showTrace=true/false)
"""
function createOptimInteriorNewton(peTabOpt::PeTabOpt;
                                   hessianUse::Symbol=:blockAutoDiff)
    
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
    dfc = TwiceDifferentiableConstraints(lowerBounds .- 0.01, upperBounds .+ 0.01)

    evalOptim = (p0; showTrace=false) -> Optim.optimize(df, 
                                                        dfc, 
                                                        p0, 
                                                        IPNewton(), 
                                                        Optim.Options(iterations = 1000, show_trace = showTrace))

    return evalOptim
end