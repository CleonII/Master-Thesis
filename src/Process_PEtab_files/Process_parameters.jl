"""
    processParameters(parametersFile::DataFrame)::ParameterInfo

    Process the PeTab parametersFile file into a type-stable Julia struct.
"""
function processParameters(parametersFile::DataFrame)::ParametersInfo

    nParameters = length(parametersFile[!, "estimate"])

    # Pre-allocate arrays to hold data 
    lowerBound::Vector{Float64} = zeros(Float64, nParameters) 
    upperBound::Vector{Float64} = zeros(Float64, nParameters)
    nominalValue::Vector{Float64} = zeros(Float64, nParameters) # Vector with Nominal value in PeTab-file
    parameterScale::Vector{Symbol} = Vector{Symbol}(undef, nParameters)
    parameterId::Vector{Symbol} = Vector{Symbol}(undef, nParameters) # Faster to do comparisons with Symbols than Strings
    estimate::Vector{Bool} = Vector{Bool}(undef, nParameters) # 

    for i in eachindex(estimate)

        # If upper or lower bounds are missing assume +Inf and -Inf respectively.
        if ismissing(parametersFile[i, "lowerBound"])
            lowerBound[i] = -Inf
        else
            lowerBound[i] = parametersFile[i, "lowerBound"]
        end
        if ismissing(parametersFile[i, "upperBound"])
            upperBound[i] = Inf
        else
            upperBound[i] = parametersFile[i, "upperBound"]
        end

        nominalValue[i] = parametersFile[i, "nominalValue"]
        parameterId[i] = Symbol(string(parametersFile[i, "parameterId"]))
        # Currently only supports parameters on log10-scale -> TODO: Change this 
        if parametersFile[i, "parameterScale"] == "log10"
            parameterScale[i] = :log10
        elseif parametersFile[i, "parameterScale"] == "log"
            parameterScale[i] = :log
        elseif parametersFile[i, "parameterScale"] == "lin"            
            parameterScale[i] = :lin
        else
            errorStr = "Parameter scale " * parametersFile[i, "parameterScale"] * "not supported. Only log10, log and lin are supported in the Parameters PEtab file under the parameterScale column."
            throw(PEtabFileError(errorStr)) 
        end
         
        estimate[i] = parametersFile[i, "estimate"] == 1 ? true : false
    end
    nParametersToEstimate::Int64 = Int64(sum(estimate))

    return ParametersInfo(nominalValue, lowerBound, upperBound, parameterId, parameterScale, estimate, nParametersToEstimate)
end