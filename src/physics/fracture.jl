@inline fracture_kwargs() = (:Gc, :epsilon_c)

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

!!! danger "Overwriting failure permission with `material!` and `failure_permit!`"
    The function `material!` calls `failure_permit!`, so if it is used afterwards,
    previously set failure permissions might be overwritten!

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
        Gc = 0.0;
        εc = 0.0;
    end

    return (; Gc, εc)
end

function set_failure_permissions!(body::AbstractBody, set_name::Symbol,
                                  params::AbstractPointParameters)
    if has_fracture(body.mat, params)
        failure_permit!(body, set_name, true)
    else
        failure_permit!(body, set_name, false)
    end
    return nothing
end

function has_fracture(mat::AbstractMaterial, params::AbstractPointParameters)
    if isapprox(params.Gc, 0; atol=eps()) || isapprox(params.εc, 0; atol=eps())
        return false
    else
        return true
    end
end

function required_fields_fracture(::Type{Material}) where {Material<:AbstractMaterial}
    fields = (req_point_data_fields_fracture(Material)...,
              req_bond_data_fields_fracture(Material)...,
              req_data_fields_fracture(Material)...)
    return fields
end

# function required_fields_fracture(::Any)
#     return ()
# end

function req_point_data_fields_fracture(::Type{Material}) where {Material<:AbstractMaterial}
    return ()
end

function req_bond_data_fields_fracture(::Type{Material}) where {Material<:AbstractMaterial}
    return ()
end

function req_data_fields_fracture(::Type{Material}) where {Material<:AbstractMaterial}
    return ()
end
