function init_time_solver!(ts::AbstractTimeSolver, dh::AbstractDataHandler)
    throw(MethodError(init_time_solver!, ts, dh))
end

function solve!(dh::AbstractDataHandler, ts::AbstractTimeSolver, options)
    msg = "interface function `solve!` not correctly implemented!\n"
    msg *= "No method for arg types:\n"
    msg *= "   dh::$(typeof(dh))\n"
    msg *= "   ts::$(typeof(ts))\n"
    error(msg)
end

function required_fields_timesolver(ts::Type{TS}) where {TS}
    throw(MethodError(required_fields_timesolver, ts))
    return nothing
end

function req_storage_fields_timesolver(::Type{S}, ::Type{TS}) where {S,TS}
    parameters = fieldnames(S)
    for req_field in required_fields_timesolver(TS)
        if !in(req_field, parameters)
            msg = "required field $req_field not found in $(S)!\n"
            msg *= "The field is required for the `$(TS)` time solver!\n"
            error(msg)
        end
    end
    return nothing
end
