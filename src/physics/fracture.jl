"""
    CriticalStretch

A damage model based on the stretch of the bond. The bond is considered to be broken
if the stretch exceeds a critical value.
The critical value can be defined via the fracture energy `Gc` or the critical stretch `εc`
using the [`material!`](@ref) function. The damage model is defined globally for the whole
body as part of the material.
"""
struct CriticalStretch <: AbstractDamageModel end

@inline fracture_kwargs() = (:Gc, :epsilon_c)

"""
    failure_permit!(body, set_name, fail_permit)

$(internal_api_warning())

Set the failure permission for points of the set `set_name` of a `body`.

# Arguments

- `body::AbstractBody`: [`Body`](@ref) where the failure permission will be set.
- `set_name::Symbol`: The name of a point set of this body.
- `fail_permit::Bool`: If `true`, failure is allowed, and if `false` then no bonds of this
    point are allowed to break during the simulation.

!!! danger "Overwriting failure permission with `material!` and `failure_permit!`"
    The function `material!` calls `failure_permit!`, so if it is used afterwards,
    previously set failure permissions might be overwritten!

# Throws

- Error if the body does not contain a set with `set_name`.
"""
function failure_permit! end

function failure_permit!(body::AbstractBody, set_name::Symbol, fail_permit::Bool)
    check_if_set_is_defined(body.point_sets, set_name)
    body.fail_permit[body.point_sets[set_name]] .= fail_permit
    return nothing
end

"""
    no_failure!(body::AbstractBody, set_name::Symbol)
    no_failure!(body::AbstractBody)

Disallow failure for all points of the point set `set_name` of the `body`.
If no `set_name` is specified, failure is prohibited for the whole `body`.

# Arguments

- `body::AbstractBody`: [`Body`](@ref) for which failure is prohibited.
- `set_name::Symbol`: The name of a point set of this body.

!!! danger "Overwriting failure permission with `material!` and `no_failure!`"
    The function `material!` sets failure permissions due to the provided input parameters,
    so if it is used afterwards, previously set failure prohibitions might be overwritten!

# Throws

- Error if the body does not contain a set with `set_name`.

# Examples

```julia-repl
julia> no_failure!(body)

julia> body
1000-point Body{BBMaterial{NoCorrection}}:
  1 point set(s):
    1000-point set `all_points`
  1000 points with failure prohibited
```
"""
function no_failure! end

function no_failure!(body::AbstractBody, set_name::Symbol)
    check_if_set_is_defined(body.point_sets, set_name)
    points_without_material = findfirst(x -> x == 0,
                                        body.params_map[body.point_sets[set_name]])
    if points_without_material !== nothing
        msg = string("Not all points of point set \":", set_name,
                     "\" have material parameters!\n")
        msg *= "Please use the `material!` function to define them before prohibiting "
        msg *= "failure, otherwise failure permissions might be overwritten!\n"
        throw(ArgumentError(msg))
    end
    failure_permit!(body, set_name, false)
    return nothing
end

function no_failure!(body::AbstractBody)
    no_failure!(body, :all_points)
    return nothing
end

function get_frac_params(::CriticalStretch, p::Dict{Symbol,Any}, δ::Float64, K::Float64)
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

function get_frac_params(::AbstractDamageModel, p, δ, K)
    return (; )
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
    return has_fracture(mat.dmgmodel, params)
end

function has_fracture(::CriticalStretch, params::AbstractPointParameters)
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

function req_point_data_fields_fracture(::Type{Material}) where {Material<:AbstractMaterial}
    return ()
end

function req_bond_data_fields_fracture(::Type{Material}) where {Material<:AbstractMaterial}
    return ()
end

function req_data_fields_fracture(::Type{Material}) where {Material<:AbstractMaterial}
    return ()
end

function log_param_property(::Val{:Gc}, param; indentation)
    return msg_qty("critical energy release rate", param.Gc; indentation)
end

function log_param_property(::Val{:εc}, param; indentation)
    return msg_qty("critical stretch", param.εc; indentation)
end
