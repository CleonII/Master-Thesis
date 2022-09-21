"""
    isNumber(x::String)::Bool

    Check if a string x is a number (Float).
"""
function isNumber(x::String)::Bool
    return tryparse(Float64, x) !== nothing
end
"""
    isNumber(x::SubString{String})::Bool
"""
function isNumber(x::SubString{String})::Bool
    return tryparse(Float64, x) !== nothing
end


"""
    transformYobsOrYmodArr!(vals, transformationArr::Array{Symbol, 1})

    Transform the Yobs or Ymod arrays (vals) in place using for each value 
    in vals the transformation specifed in transformationArr.

    Currently :lin, :log10 and :log transforamtions are supported, see
    `setUpCostFunc`.
"""
function transformYobsOrYmodArr!(vals, transformationArr::Array{Symbol, 1})
    for i in eachindex(vals)
        vals[i] = transformObsOrData(vals[i], transformationArr[i])
    end
end


"""
    transformYobsOrYmodArr!(vals, transformationArr::Array{Symbol, 1})

    Transform val using either :lin (identify), :log10 and :log transforamtions.
"""
function transformObsOrData(val, transform::Symbol)
    if transform == :lin
        return val
    elseif transform == :log10
        return log10(val)
    elseif transform == :log
        return log(val)
    else
        println("Error : $transform is not an allowed transformation")
        println("Only :lin, :log10 and :log are supported.")
    end
end
