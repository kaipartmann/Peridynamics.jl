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
