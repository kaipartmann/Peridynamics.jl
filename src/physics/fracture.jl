

"""
    failure_permit!(body, fail_permit)
    failure_permit!(body, set, fail_permit)

Determines whether failure is permitted `fail_permit = true` or prohibited
`fail_permit = false` in the `body` or the point set `set` of the `body`.

# Arguments

- `body::AbstractBody`: Peridynamic body
- `set::Symbol`: Point set on `body`
- `fail_permit::Bool`: Decides if failure is allowed on considered body or point set.

# Throws

- Error if no point set called `set` exists

# Example

```julia-repl
julia> failure_permit!(b, :set_bottom, false)
julia> b.fail_permit
12500-element Vector{Bool}:
 0
 0
 0
 ⋮
 1
 1
```
"""
function failure_permit! end

function failure_permit!(b::AbstractBody, fail_permit::Bool)
    b.fail_permit .= fail_permit
    return nothing
end

function failure_permit!(b::AbstractBody, name::Symbol, fail_permit::Bool)
    check_if_set_is_defined(b.point_sets, name)
    b.fail_permit[b.point_sets[name]] .= fail_permit
    return nothing
end


function get_frac_params(p::Dict{Symbol,Any}, δ::Float64, K::Float64)
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

    return Gc, εc
end

@inline function calc_damage!(b::AbstractBodyChunk)
    for point_id in each_point_idx(b)
        dmg = 1 - b.storage.n_active_bonds[point_id] / b.system.n_neighbors[point_id]
        b.storage.damage[point_id] = dmg
    end
    return nothing
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
