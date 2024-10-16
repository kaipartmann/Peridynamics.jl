const REGISTERED_SOLVERS = DataType[]

@inline registered_solvers() = REGISTERED_SOLVERS

function register_solver!(::Type{Solver}) where {Solver<:AbstractTimeSolver}
    push!(REGISTERED_SOLVERS, Solver)
    return nothing
end

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

function required_fields_timesolvers()
    all_fields = Vector{Symbol}()
    for solver in registered_solvers()
        push!(all_fields, req_point_data_fields_timesolver(solver)...)
        push!(all_fields, req_bond_data_fields_timesolver(solver)...)
        push!(all_fields, req_data_fields_timesolver(solver)...)
    end
    return Tuple(unique(all_fields))
end

function required_point_data_fields_timesolvers()
    all_fields = Vector{Symbol}()
    for solver in registered_solvers()
        push!(all_fields, req_point_data_fields_timesolver(solver)...)
    end
    return Tuple(unique(all_fields))
end

function required_bond_data_fields_timesolvers()
    all_fields = Vector{Symbol}()
    for solver in registered_solvers()
        push!(all_fields, req_bond_data_fields_timesolver(solver)...)
    end
    return Tuple(unique(all_fields))
end

function required_data_fields_timesolvers()
    all_fields = Vector{Symbol}()
    for solver in registered_solvers()
        push!(all_fields, req_data_fields_timesolver(solver)...)
    end
    return Tuple(unique(all_fields))
end

function req_point_data_fields_timesolver(::Type{TS}) where {TS}
    return throw(InterfaceError(TS, "req_point_data_fields_timesolver"))
end

function req_bond_data_fields_timesolver(::Type{TS}) where {TS}
    return throw(InterfaceError(TS, "req_bond_data_fields_timesolver"))
end

function req_data_fields_timesolver(::Type{TS}) where {TS}
    return throw(InterfaceError(TS, "req_data_fields_timesolver"))
end

function log_timesolver(::AbstractJobOptions, ::AbstractTimeSolver)
    return nothing
end

@inline function init_field_solver(solver, system, field)
    return nothing
end
