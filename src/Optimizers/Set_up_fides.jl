using PyCall
import Pkg


function setUpFides(peTabOpt::PeTabOpt,
                    autoDiffHess::Symbol; 
                    fidesHessApprox=py"None"o, 
                    verbose::Integer=1,
                    options=py"{'maxiter' : 1000}"o,
                    funargs=py"None"o,
                    resfun::Bool=false)

    nParam = length(peTabOpt.lowerBounds)
    if autoDiffHess == :autoDiff
        useHessApprox = false
        hessMat = zeros(Float64, (nParam, nParam))
        evalHessian = (pVec) -> evalAutoDiffHess(pVec, peTabOpt.evalHess, hessMat)
    elseif autoDiffHess == :blockAutoDiff
        useHessApprox = false
        hessMat = zeros(Float64, (nParam, nParam))
        evalHessian = (pVec) -> evalAutoDiffHess(pVec, peTabOpt.evalHessApprox, hessMat)        
    elseif autoDiffHess == :None
        useHessApprox = true
    else
        println("Error : For Fides availble autoDiffHess options are :autoDiff, :blockAutoDiff or :None" )
        println("User provided : $autoDiffHess")
    end

    if fidesHessApprox == py"None"o && autoDiffHess == :None
        println("Error : User must provide either a fides hessian approximation or specify if the hessian should 
                 be computed via autodiff.")
    end 

    gradient = zeros(Float64, nParam)
    evalGradF = (pVec) -> evalAutoDiffGrad(pVec, peTabOpt.evalGradF, gradient)

    # Runnable fides function 
    if useHessApprox == false
        fidesFunc = (pVec) -> fidesObjHess(pVec, peTabOpt.evalF, evalGradF, evalHessian)
    else
        fidesFunc = (pVec) -> fidesObjApprox(pVec, peTabOpt.evalF, evalGradF)
    end

    fidesObj = setUpFidesClass(fidesFunc, peTabOpt.upperBounds, peTabOpt.lowerBounds, 
                               verbose, 
                               options, 
                               funargs, 
                               fidesHessApprox, 
                               resfun)

    return fidesObj
end


function setUpFidesClass(fun, 
                         ub, 
                         lb,   
                         verbose,
                         options,
                         funargs,
                         hessian_update,
                         resfun::Bool)

    py"""
    import numpy as np
    import fides
    import logging

    def run_fides_python(x0, fun, ub, lb, verbose, options, funargs, hessian_update, resfun):
        fides_opt = fides.Optimizer(fun, ub, lb, 
                                    verbose=verbose,
                                    options=options,
                                    funargs=funargs,
                                    hessian_update=hessian_update, 
                                    resfun=resfun)

        opt_res = fides_opt.minimize(x0)
        n_iter = fides_opt.iteration
        converged = fides_opt.converged
        return opt_res, n_iter, converged

    """

    runFidesJulia = (x0; verbose=verbose, options=options, funargs=funargs, hessian_update=hessian_update, resfun=resfun) -> py"run_fides_python"(x0, fun, ub, lb, verbose, options, funargs, hessian_update, resfun)
    return runFidesJulia
end


# Helper function ensuring the hessian is returned as required by fides 
function evalAutoDiffHess(p,
                          evalHess!::Function, 
                          hessMat::Matrix{Float64})
    hessMat .= 0
    evalHess!(hessMat, p)
    return hessMat[:, :]
end


# Helper function ensuring the gradient is returned as required by fides 
function evalAutoDiffGrad(p, 
                          evelGrad!::Function, 
                          gradient::Vector{Float64})
    gradient .= 0.0
    evelGrad!(gradient, p)
    return gradient[:]
end


# Helper functions wrapping the output into a Fides acceptable format when hessian is computed by autodiff 
function fidesObjHess(p, evalF::Function, evalGradF::Function, evalHessF::Function)
    fVal = evalF(p)
    fGradVal = evalGradF(p)
    fHessVal = evalHessF(p)
    return (fVal, fGradVal, fHessVal)
end


# Helper functions wrapping the output into a Fides acceptable format when hessian approximation is used
function fidesObjApprox(p, evalF::Function, evalGradF::Function)
    fVal = evalF(p)
    fGradVal = evalGradF(p)
    return (fVal, fGradVal)
end


# Helper function to load the correct python environment 
function loadFidesFromPython(pathPythonExe::String)

    println("Loading Fides into Julia using the python environment at:")
    println(pathPythonExe)

    ENV["PYTHON"] = pathPythonExe
    Pkg.build("PyCall")
end
