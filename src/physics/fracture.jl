"""
    failure_permit!(body, fail_permit)
    failure_permit!(body, set_name, fail_permit)

Set the failure permission for points of a `body`. By default, failure is allowed for
all points. If no `set_name` is specified, then the permission `fail_permit` is set for all
points of the body.

# Arguments

- `body::AbstractBody`: [`Body`](@ref) where the failure permission will be set.
- `set_name::Symbol`: The name of a point set of this body.
- `fail_permit::Bool`: If `true`, failure is allowed, and if `false` then no bonds of this
    point are allowed to break during the simulation.

# Throws

- Errors if the body does not contain a set with `set_name`.

# Examples

```julia-repl
julia> failure_permit!(body, false)

julia> body
1000-point Body{BBMaterial{NoCorrection}}:
  1 point set(s):
    1000-point set `all_points`
  1000 points with no failure permission
```
"""
function failure_permit! end

function failure_permit!(body::AbstractBody, fail_permit::Bool)
    body.fail_permit .= fail_permit
    return nothing
end

function failure_permit!(body::AbstractBody, set_name::Symbol, fail_permit::Bool)
    check_if_set_is_defined(body.point_sets, set_name)
    body.fail_permit[body.point_sets[set_name]] .= fail_permit
    return nothing
end

function get_frac_params(p::Dict{Symbol,Any}, params::NamedTuple)
    (; δ, K) = params
    local Gc::Float64
    local εc::Float64
    if haskey(p, :Gc) && !haskey(p, :epsilon_c)
        Gc = float(p[:Gc])
        εc = sqrt(5.0 * Gc / (9.0 * K * δ))
    elseif !haskey(p, :Gc) && haskey(p, :epsilon_c)
        εc = float(p[:epsilon_c])
        Gc = 9.0 / 5.0 * K * δ * εc^2
    elseif haskey(p, :Gc) && haskey(p, :epsilon_c)
        msg = "insufficient keywords for calculation of fracture parameters!\n"
        msg *= "Define either Gc or epsilon_c, not both!\n"
        throw(ArgumentError(msg))
    else
        msg = "insufficient keywords for calculation of fracture parameters!\n"
        msg *= "Define either Gc or epsilon_c!\n"
        throw(ArgumentError(msg))
    end
    return (; Gc, εc)
end


function required_fields_fracture()
    return (req_point_data_fields_fracture()..., req_data_fields_fracture()...)
end

function req_point_data_fields_fracture()
    return (:damage, :n_active_bonds)
end

function req_data_fields_fracture()
    return ()
end

function req_storage_fields_fracture(::Type{Storage}) where {Storage}
    parameters = fieldnames(Storage)
    for req_field in required_fields_fracture()
        if !in(req_field, parameters)
            msg = "required field $req_field not found in $(Storage)!\n"
            msg *= "The field is required for the fracture calculation!\n"
            error(msg)
        end
    end
    return nothing
end
