# Function will change when a PeTabOpt struct is created.
function createOptimProb(fPre::Function, 
                         fGradPre::Function, 
                         fHessPre::Function, 
                         lowerBounds::Array{<:AbstractFloat, 1}, 
                         upperBounds::Array{<:AbstractFloat}; 
                         showTrace::Bool=false)
    
    x0 = zeros(Float64, length(lowerBounds))
    fGradUse = (grad, x) -> fGradPre(x, grad)
    fHessUse = (hess, x) -> fHessPre(hess, x)
    df = TwiceDifferentiable(fPre, fGradUse, fHessUse, x0)
    dfc = TwiceDifferentiableConstraints(lowerBounds .- 0.01, upperBounds .+ 0.01)

    evalOptim = (p0) -> Optim.optimize(df, dfc, p0, IPNewton(), Optim.Options(iterations = 1000, show_trace = showTrace))

    return evalOptim
end