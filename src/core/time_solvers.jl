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

"""
    required_fields_timesolver(::Type{TS})

$(internal_api_warning())

Return all storage fields required by the time solver type `TS`.
"""
function required_fields_timesolver(::Type{TS}) where {TS<:AbstractTimeSolver}
    return (req_point_data_fields_timesolver(TS)...,
            req_bond_data_fields_timesolver(TS)...,
            req_data_fields_timesolver(TS)...)
end

function req_storage_fields(solver::AbstractTimeSolver)
    # A custom time solver is not required to declare its fields. If the interface methods
    # are not defined, then nothing can be checked for this solver.
    return try
        required_fields_timesolver(typeof(solver))
    catch err
        err isa InterfaceError || rethrow()
        ()
    end
end

# Informational only: the union of the fields of all registered solvers. Storages are
# checked against the solver that is actually used, see `check_storage_contract`, so
# registering a solver cannot invalidate storages of other solvers anymore.
function required_fields_timesolvers()
    all_fields = Vector{Symbol}()
    for solver in registered_solvers()
        push!(all_fields, required_fields_timesolver(solver)...)
    end
    return Tuple(unique(all_fields))
end

function req_point_data_fields_timesolver(::Type{TS}) where {TS}
    return throw(InterfaceError(TS, "req_point_data_fields_timesolver",
                                timesolver_fields_hint(TS, "req_point_data_fields_timesolver",
                                                       "(:velocity, :b_int)")))
end

function req_bond_data_fields_timesolver(::Type{TS}) where {TS}
    return throw(InterfaceError(TS, "req_bond_data_fields_timesolver",
                                timesolver_fields_hint(TS, "req_bond_data_fields_timesolver",
                                                       "(:bond_active,)")))
end

function req_data_fields_timesolver(::Type{TS}) where {TS}
    return throw(InterfaceError(TS, "req_data_fields_timesolver",
                                timesolver_fields_hint(TS, "req_data_fields_timesolver",
                                                       "(:residual,)")))
end

function timesolver_fields_hint(TS, method, example)
    msg = "A time solver declares the storage fields it works with, so that the storage "
    msg *= "contract\n    can be checked when a `Job` is created:\n"
    msg *= "        Peridynamics.$(method)(::Type{<:$(nameof(TS))}) = $(example)\n"
    msg *= "    Return an empty tuple `()` if the solver works with no field of this kind."
    return msg
end

function log_timesolver(::AbstractJobOptions, ::AbstractTimeSolver)
    return nothing
end

@inline function init_field_solver(solver, system, field)
    return nothing
end
