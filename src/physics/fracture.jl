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


"""
    get_frac_params(dmgmodel, δ, K; Gc=nothing, epsilon_c=nothing, kwargs...)

$(extension_api_note())

Read or calculate the fracture parameters of a damage model from the fracture keywords of
[`material!`](@ref). This function has to be defined when creating a new damage model.
Otherwise, a default method returns an empty named tuple `(; )`.

A keyword the user did not specify is `nothing`, which is how a damage model decides which
of the fracture keywords it accepts and how it converts them into each other. Every fracture
keyword is passed, so a method should end in `kwargs...` to stay valid when a keyword is
added.

# Arguments
- `dmgmodel::AbstractDamageModel`: The damage model
- `δ::Float64`: Horizon
- `K::Float64`: Bulk modulus

# Keywords
- `Gc`: Critical energy release rate
- `epsilon_c`: Critical strain

# Example
```julia
function Peridynamics.get_frac_params(::MyDamage, δ, K; Gc=nothing, epsilon_c=nothing,
                                      kwargs...)
    isnothing(Gc) && return (; Gc=0.0, εc=0.0)
    return (; Gc, εc=sqrt(5.0 * Gc / (9.0 * K * δ)))
end
```
"""
function get_frac_params end

function get_frac_params(::CriticalStretch, δ::Float64, K::Float64; Gc=nothing,
                         epsilon_c=nothing, kwargs...)
    local _Gc::Float64
    local εc::Float64

    if !isnothing(Gc) && isnothing(epsilon_c)
        _Gc = float(Gc)
        εc = sqrt(5.0 * _Gc / (9.0 * K * δ))
    elseif isnothing(Gc) && !isnothing(epsilon_c)
        εc = float(epsilon_c)
        _Gc = 9.0 / 5.0 * K * δ * εc^2
    elseif !isnothing(Gc) && !isnothing(epsilon_c)
        msg = "insufficient keywords for calculation of fracture parameters!\n"
        msg *= "Define either Gc or epsilon_c, not both!\n"
        throw(ArgumentError(msg))
    else
        _Gc = 0.0;
        εc = 0.0;
    end

    return (; Gc=_Gc, εc)
end

function get_frac_params(::AbstractDamageModel, δ, K; kwargs...)
    return (; )
end

# the point parameter constructors of `@params` still read the fracture keywords from the
# `Dict` that `material!` collects; this bridge goes when the parameters are keyword-based
function get_frac_params(dmgmodel::AbstractDamageModel, p::Dict{Symbol,Any}, δ, K)
    return get_frac_params(dmgmodel, δ, K;
                           (kw => get(p, kw, nothing) for kw in fracture_kwargs())...)
end

"""
    set_failure_permissions!(body, set_name, params)

$(internal_api_warning())

Grant or prohibit failure permission depending on the submitted fracture parameters by
calling [`failure_permit!`](@ref).

If fracture parameters are found, failure is allowed. If no fracture parameters are found,
failure is not allowed.
"""
function set_failure_permissions!(body::AbstractBody, set_name::Symbol,
                                  params::AbstractPointParameters)
    if has_fracture(body.mat, params)
        failure_permit!(body, set_name, true)
    else
        failure_permit!(body, set_name, false)
    end
    return nothing
end

"""
    calc_failure!(storage, system, mat, dmgmodel, paramsetup, i)

$(extension_api_note())

Decide which bonds of point `i` have failed and update the fracture bookkeeping of the
storage accordingly. This is the one method a damage model has to define; it is called once
per local point and per time step, right before [`force_density_point!`](@ref).

A method sets `storage.bond_active[bond_id] = false` for every bond that fails and keeps
`storage.n_active_bonds[i]` in sync, because that count is what `calc_damage!` turns into the
damage of the point. A bond whose `bond.fail_permit` is `false` must never fail, which is how
`no_failure!` and the pre-cracks are honored.

# Arguments

- `storage`: The storage of the body chunk.
- `system`: The system of the body chunk.
- `mat`: The material.
- `dmgmodel`: The damage model, i.e. what a new model dispatches on.
- `paramsetup`: The parameters of the body chunk; resolve them with [`get_params`](@ref).
- `i::Int`: The index of the local point that is evaluated.

# Example

```julia
struct MyDamage <: Peridynamics.AbstractDamageModel end

function Peridynamics.calc_failure!(storage, system, mat, ::MyDamage, paramsetup, i)
    params = Peridynamics.get_params(paramsetup, i)
    for bond_id in Peridynamics.each_bond_idx(system, i)
        bond = system.bonds[bond_id]
        ε = ... # the stretch of the bond
        if ε > params.εc && bond.fail_permit
            storage.bond_active[bond_id] = false
        end
        storage.n_active_bonds[i] += storage.bond_active[bond_id]
    end
    return nothing
end
```

See also [`get_frac_params`](@ref), [`has_fracture`](@ref), `AbstractDamageModel`.
"""
function calc_failure! end

"""
    has_fracture(mat, params)

$(internal_api_warning())

Return `true` if at least one fracture parameter is set `!=0` in `params` and the system
therefore is supposed to have failure allowed or return `false` if not.
"""
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

"""
    get_dmgmodel(mat::AbstractMaterial)

$(extension_api_note())

Return the damage model of `mat`, or `nothing` if the material has no damage model.
"""
@inline function get_dmgmodel(mat::AbstractMaterial)
    hasproperty(mat, :dmgmodel) || return nothing
    return mat.dmgmodel
end

# a damage model that declares a state with `@dmg_storage` needs a storage that carries it
function req_storage_fields(::AbstractMaterial, dmgmodel::AbstractDamageModel)
    return damage_storage_type(dmgmodel) === Nothing ? () : (:dmg_state,)
end
req_storage_fields(::AbstractMaterial, ::Nothing) = ()

# --------------------------------------------------------------------------------------
# state of a damage model
#
# The damage-model twin of the constitutive-model state (`cm_state::ConstitutiveState`):
# a damage model brings the per-bond variables it needs instead of every material having
# to allocate them for it. Note what is deliberately *absent*: nothing here marks a model
# as history dependent. A damage model integrates its state in `calc_failure!`, which runs
# exactly once per step under every solver that supports fracture, so a stateful damage
# model stays compatible with solvers that evaluate the force density several times per
# step.
# --------------------------------------------------------------------------------------

"""
    damage_storage_type(dmgmodel)
    damage_storage_type(dmgmodel, ::Type{FT})

$(extension_api_note())

Return the type of the state of a damage model, instantiated for the float type `FT` of the
simulation, or `Nothing` for a model without state. This method is generated by
[`@dmg_storage`](@ref) and is the damage-model analogue of [`storage_type`](@ref).
"""
function damage_storage_type end

damage_storage_type(dmgmodel) = Nothing
damage_storage_type(dmgmodel, ::Type) = damage_storage_type(dmgmodel)

"""
    get_dmg_storage(dmgmodel, solver, system)

$(internal_api_warning())

Allocate the state of a damage model for one body chunk, or return `nothing` for a model
without state. This method is generated by [`@dmg_storage`](@ref) and is the damage-model
analogue of `get_storage`.
"""
function get_dmg_storage end

get_dmg_storage(dmgmodel, solver, system) = nothing

"""
    init_damage_state(mat, solver, system)

$(internal_api_warning())

Return the value of the storage field that was declared with `dmg_state::DamageState`,
which by default is the state of the damage model of `mat`. This is the escape hatch for a
material family that assembles the state of its damage model itself.
"""
function init_damage_state(mat, solver, system)
    return get_dmg_storage(get_dmgmodel(mat), solver, system)
end

"""
    damage_state(storage)

$(extension_api_note())

Return the state of the damage model that is carried by a storage, i.e. the field declared
with `dmg_state::DamageState`, or `nothing` for a storage that does not declare one. This is
how [`calc_failure!`](@ref), [`calc_damage!`](@ref), [`kinematic_weight`](@ref) and
[`safe_degradation`](@ref) reach the per-bond variables of a stateful damage model.
"""
function damage_state end

@inline damage_state(::AbstractStorage) = nothing

"""
    has_damage_state(::Type{Storage})

$(internal_api_warning())

Return whether a storage type declares a field with the [`DamageState`](@ref) declaration
and can therefore carry the state of a damage model.
"""
function has_damage_state end

has_damage_state(::Type) = false

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

# --------------------------------------------------------------------------------------
# continuous degradation
#
# A damage model that softens a bond instead of deleting it enters the kinematics with
# `kinematic_weight` and the force with `degrade_bond_stress`; `safe_degradation` is the
# scalar factor both are built from. The defaults below are what every model that deletes
# bonds wants, and they cost nothing: the weight is the constant one and the stress passes
# through by dispatch.
# --------------------------------------------------------------------------------------

"""
    kinematic_weight(dmgmodel, storage, bond_id)

$(extension_api_note())

Return the factor in `[0, 1]` by which bond `bond_id` takes part in the kinematics of a
correspondence material, i.e. in the moment matrix and the gradient weights. It is the
*kinematic* counterpart of [`safe_degradation`](@ref): that one scales what a bond carries,
this one scales what it contributes to the deformation gradient.

The default is `1.0` for every damage model, i.e. a bond is either fully present or, once
it has failed, absent. A model that softens a bond instead defines a method reading its own
state, see [`@dmg_storage`](@ref) and [`damage_state`](@ref).

See also [`safe_degradation`](@ref), [`degrade_bond_stress`](@ref), [`calc_failure!`](@ref).
"""
function kinematic_weight end

@inline kinematic_weight(::AbstractDamageModel, ::AbstractStorage, ::Integer) = 1.0

"""
    safe_degradation(dmgmodel, storage, bond_id)

$(extension_api_note())

Return the factor in `[0, 1]` by which the stress carried by bond `bond_id` is degraded. It
is the *static* counterpart of [`kinematic_weight`](@ref): that one scales what a bond
contributes to the moment matrix and the gradient weights, this one scales what it carries.

The default is `1.0` for every damage model, i.e. no degradation, which is what a model that
deletes bonds outright wants. A model that softens a bond instead defines a method reading
its own state, e.g.

```julia
@inline function Peridynamics.safe_degradation(::MyDamage,
                                               storage::Peridynamics.AbstractStorage,
                                               bond_id)
    @inbounds d = Peridynamics.damage_state(storage).bond_damage[bond_id]
    return (1 - d)^2
end
```

Degrading a bond continuously rather than deleting it keeps the moment matrix of a
correspondence material a continuous function of the deformation, which is what makes
fragmentation stable; a deleted bond is a jump.

See also [`kinematic_weight`](@ref), [`calc_failure!`](@ref), [`@dmg_storage`](@ref).
"""
function safe_degradation end

@inline function safe_degradation(::AbstractDamageModel, ::AbstractStorage, bond_id)
    return 1.0
end

"""
    degrade_bond_stress(dmgmodel, storage, bond_id, P)

$(extension_api_note())

Return the first Piola-Kirchhoff stress that bond `bond_id` actually carries, given the
undegraded stress `P` its constitutive model produced. This is where a damage model that
softens a bond instead of deleting it enters the **force**, next to
[`kinematic_weight`](@ref), which is where it enters the **kinematics**.

The default returns `P` unchanged. Note that it does so by dispatch and not by multiplying
with `1.0`: a material whose damage model does not degrade must compile to exactly the code
it did before, and a floating-point multiplication by one is not something the compiler is
allowed to remove.

A model that degrades opts in with one line, reusing its own [`safe_degradation`](@ref):

```julia
@inline function Peridynamics.degrade_bond_stress(dmg::MyDamage,
                                                  storage::Peridynamics.AbstractStorage,
                                                  bond_id, P)
    return Peridynamics.safe_degradation(dmg, storage, bond_id) * P
end
```

See also [`safe_degradation`](@ref), [`kinematic_weight`](@ref), [`@dmg_storage`](@ref).
"""
function degrade_bond_stress end

@inline degrade_bond_stress(::AbstractDamageModel, ::AbstractStorage, bond_id, P) = P

# --------------------------------------------------------------------------------------
# logging
# --------------------------------------------------------------------------------------

"""
    log_dmgmodel(dmgmodel; indentation)

$(internal_api_warning())

Generate a logging message for a damage model. By default, logs only the damage model type.
A damage model with properties worth logging defines its own method.
"""
function log_dmgmodel(dmgmodel::AbstractDamageModel; indentation)
    return msg_qty("damage model type", typeof(dmgmodel); indentation)
end

function log_param_property(::Val{:Gc}, param; indentation)
    return msg_qty("critical energy release rate", param.Gc; indentation)
end

function log_param_property(::Val{:εc}, param; indentation)
    return msg_qty("critical stretch", param.εc; indentation)
end
