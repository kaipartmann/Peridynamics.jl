"""
    CriticalStretch

A damage model based on the stretch of the bond. The bond is considered to be broken
if the stretch exceeds a critical value.
The critical value can be defined via the fracture energy `Gc` or the critical stretch `εc`
using the [`material!`](@ref) function. The damage model is defined globally for the whole
body as part of the material.
"""
struct CriticalStretch <: AbstractDamageModel end


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
    critical_stretch(dmgmodel, mat, δ, K, Gc)

$(extension_api_note())

Return the critical stretch `εc` that belongs to the critical energy release rate `Gc`, for
the damage model `dmgmodel` of the material `mat` with the horizon `δ` and the bulk modulus
`K`.

# Default

`εc = sqrt(5 Gc / (9 K δ))`, the relation of the constant micro-modulus of bond-based
peridynamics. A material with another micro-modulus defines this method **and**
[`energy_release_rate`](@ref), always both, so that `Gc` and `epsilon_c` keep converting
into each other.

# Example

```julia
Peridynamics.critical_stretch(::CriticalStretch, ::MyMaterial, δ, K, Gc) =
    sqrt(2 * Gc / (3 * K * δ))
Peridynamics.energy_release_rate(::CriticalStretch, ::MyMaterial, δ, K, εc) =
    1.5 * K * δ * εc^2
```

See also [`energy_release_rate`](@ref), [`get_frac_params`](@ref).
"""
critical_stretch(dmgmodel, mat, δ, K, Gc) = sqrt(5.0 * Gc / (9.0 * K * δ))

"""
    energy_release_rate(dmgmodel, mat, δ, K, εc)

$(extension_api_note())

Return the critical energy release rate `Gc` that belongs to the critical stretch `εc`, the
inverse of [`critical_stretch`](@ref).

# Default

`Gc = 9/5 K δ εc^2`, the relation of the constant micro-modulus. Whoever defines one of the
two defines the other, so that `Gc` and `epsilon_c` stay consistent whichever one the user
gives.

See also [`critical_stretch`](@ref), [`get_frac_params`](@ref).
"""
energy_release_rate(dmgmodel, mat, δ, K, εc) = 9.0 / 5.0 * K * δ * εc^2

"""
    get_frac_params(dmgmodel, mat, δ, K; kwargs...)

$(extension_api_note())

Return the fracture parameters `Gc` and `εc` of a damage model as a `NamedTuple`, resolved
from the fracture keywords of [`material!`](@ref). This is the provider of the
[`FractureParameters`](@ref) block. A method of its own is only for a model that reads
*other* keywords, and a material with another micro-modulus defines the two conversion
hooks rather than this method.

# Arguments

- `dmgmodel`: The damage model.
- `mat`: The material, so that a conversion can depend on it.
- `δ`: The horizon.
- `K`: The bulk modulus.

# Keywords

- `Gc`: The critical energy release rate, or `nothing` if not given.
- `epsilon_c`: The critical stretch, or `nothing` if not given.

# Default

`Gc` and `epsilon_c` are converted into each other with [`critical_stretch`](@ref) and
[`energy_release_rate`](@ref). Giving both is an error, giving neither switches fracture
off with `Gc = εc = 0`. A damage model with the standard fracture keywords therefore
defines nothing here.

# Example

Every fracture keyword arrives as a keyword argument, and one the user did not give arrives
as `nothing`, which is how a method decides what it accepts.

```julia
function Peridynamics.get_frac_params(::MyDamage, mat, δ, K; tau_c=nothing, kwargs...)
    isnothing(tau_c) && return (; Gc=0.0, εc=0.0, τc=0.0)
    ...
end
```

See also [`critical_stretch`](@ref), [`energy_release_rate`](@ref).
"""
function get_frac_params(dmgmodel::AbstractDamageModel, mat, δ, K; Gc=nothing,
                         epsilon_c=nothing, kwargs...)
    if !isnothing(Gc) && !isnothing(epsilon_c)
        msg = "insufficient keywords for calculation of fracture parameters!\n"
        msg *= "Define either Gc or epsilon_c, not both!\n"
        throw(ArgumentError(msg))
    elseif !isnothing(Gc)
        return (; Gc=float(Gc), εc=float(critical_stretch(dmgmodel, mat, δ, K, Gc)))
    elseif !isnothing(epsilon_c)
        return (; Gc=float(energy_release_rate(dmgmodel, mat, δ, K, epsilon_c)),
                εc=float(epsilon_c))
    end
    return (; Gc=0.0, εc=0.0)
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
    calc_failure!(storage, system, mat, dmgmodel, paramsetup, t, Δt, i)

$(extension_api_note())

Decide which bonds of point `i` have failed and update the fracture bookkeeping
accordingly. This is the one method a damage model has to define. It is called once
per local point and per time step, right before [`force_density_point!`](@ref).

The bookkeeping lives in the state of the model itself, declared with
[`@dmg_storage`](@ref) by inheriting the block of the system family, e.g.
[`BondFracFields`](@ref), and read flat off the storage or with [`damage_state`](@ref).

Rules a method has to follow:

- Reset the count of the point with `storage.n_active_bonds[i] = 0`, set
  `storage.bond_active[bond_id] = false` for every bond that fails and add every bond that
  is still active to `storage.n_active_bonds[i]`, because that count is what
  [`calc_damage!`](@ref) turns into the damage of the point.
- A bond for which [`bond_may_fail`](@ref) is `false` must never fail, which is how
  [`no_failure!`](@ref) and the pre-cracks are honored.
- Read the stretch of a bond with [`bond_stretch`](@ref) and its current length with
  [`current_bond_length`](@ref), never by gathering the two positions and taking the norm.

# Arguments

- `storage`: The storage of the body chunk.
- `system`: The system of the body chunk.
- `mat`: The material.
- `dmgmodel`: The damage model, i.e. what a new model dispatches on.
- `paramsetup`: The parameters of the body chunk. Resolve them with [`get_params`](@ref).
- `t::Real`: The current simulation time.
- `Δt::Real`: The current time step.
- `i::Int`: The index of the local point that is evaluated.

# Default

An [`InterfaceError`](@ref) that asks for a method of the damage model.

# Example

The criterion of [`CriticalStretch`](@ref), which is what a model with another criterion
replaces:

```julia
function Peridynamics.calc_failure!(storage, system, mat, ::MyDamage, paramsetup, t, Δt, i)
    (; εc) = Peridynamics.get_params(paramsetup, i)
    storage.n_active_bonds[i] = 0
    for bond_id in Peridynamics.each_bond_idx(system, i)
        ε = Peridynamics.bond_stretch(storage, system, i, bond_id)
        if ε > εc && Peridynamics.bond_may_fail(system, bond_id)
            storage.bond_active[bond_id] = false
        end
        storage.n_active_bonds[i] += storage.bond_active[bond_id]
    end
    return nothing
end
```

See also [`calc_damage!`](@ref), [`bond_stretch`](@ref), [`@dmg_storage`](@ref).
"""
function calc_failure!(storage, system, mat, dmgmodel::AbstractDamageModel, paramsetup, t, Δt,
                       i)
    name = nameof(typeof(dmgmodel))
    hint = "Define the failure criterion of `$(name)`, e.g.\n"
    hint *= "        function Peridynamics.calc_failure!(storage, system, mat, ::$(name), "
    hint *= "paramsetup, t, Δt, i)\n"
    hint *= "            ...\n"
    hint *= "        end"
    return throw(InterfaceError(dmgmodel, "calc_failure!", hint))
end

"""
    calc_damage!(storage, system, mat, dmgmodel, paramsetup, i)

$(extension_api_note())

Reduce the bond-wise state of the damage model to the damage of point `i`, which is the
scalar written to `storage.damage` and exported as the `:damage` field. It is called once
per local point and per time step, directly after [`calc_failure!`](@ref), and
[`get_damage`](@ref) is the reading counterpart.

# Arguments

- `storage`: The storage of the body chunk. A stateful damage model reaches its own bond
    fields with [`damage_state`](@ref).
- `system`: The system of the body chunk.
- `mat`: The material.
- `dmgmodel`: The damage model, i.e. what a new model dispatches on.
- `paramsetup`: The parameters of the body chunk. Resolve them with [`get_params`](@ref).
- `i::Int`: The index of the local point that is evaluated.

# Default

The fraction of broken bonds, `1 - n_active_bonds[i] / n_neighbors[i]`, written for the
standard bookkeeping and a no-op when the storage does not carry it. A model that inherits
[`BondFracFields`](@ref) therefore needs no method of its own. A model that degrades a bond
continuously instead of deleting it defines one, so that a partially damaged bond is
counted with its degree of damage rather than as intact.

See also [`calc_failure!`](@ref), [`get_damage`](@ref), [`bond_integrity`](@ref),
[`kinematic_weight`](@ref), [`@dmg_storage`](@ref).
"""
function calc_damage! end

"""
    has_fracture(mat, params)
    has_fracture(dmgmodel, params)

$(extension_api_note())

Return whether the point parameters `params` enable fracture, which decides whether the
bonds of the points that [`material!`](@ref) assigns them to may fail.

# Default

Fracture is enabled when `params` carry `Gc` and `εc` and both are nonzero, which is the
case for every damage model that inherits [`FractureParameters`](@ref) as soon as the user
gives `Gc` or `epsilon_c`.

# Example

A model whose own parameters are the fracture parameters answers for itself.

```julia
Peridynamics.has_fracture(::MyDamage, params) = true
```

See also [`get_frac_params`](@ref), [`FractureParameters`](@ref).
"""
function has_fracture(mat::AbstractMaterial, params::AbstractPointParameters)
    return has_fracture(get_dmgmodel(mat), params)
end

function has_fracture(::AbstractDamageModel, params)
    (hasproperty(params, :Gc) && hasproperty(params, :εc)) || return false
    return !(isapprox(params.Gc, 0; atol=eps()) || isapprox(params.εc, 0; atol=eps()))
end

# a material without a damage model has nothing that could break a bond
function has_fracture(::Nothing, params)
    return false
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
function req_storage_fields(mat::AbstractMaterial, dmgmodel::AbstractDamageModel)
    return damage_storage_type(dmgmodel, system_type(mat)) === Nothing ? () : (:dmg_state,)
end
req_storage_fields(::AbstractMaterial, ::Nothing) = ()

# a material without a damage model has no failure to evaluate and no damage to reduce, so
# the standard force density loop runs with both steps as no-ops
function calc_failure!(storage, system, mat, dmgmodel::Nothing, paramsetup, t, Δt, i)
    return nothing
end

function calc_damage!(storage, system, mat, dmgmodel::Nothing, paramsetup, i)
    return nothing
end

# --------------------------------------------------------------------------------------
# state of a damage model
#
# The twin of the constitutive-model state in `core/constitutive_models.jl`: a damage
# model brings the per-bond variables it needs instead of every material having to
# allocate them for it. Note what is deliberately *absent*: nothing here touches
# `is_history_dependent`. A damage model integrates its state in `calc_failure!`, which
# runs exactly once per step under every solver that supports fracture, so a stateful
# damage model stays usable under solvers that declare
# `supports_history_dependence(solver) == false`.
# --------------------------------------------------------------------------------------

"""
    damage_storage_type(dmgmodel, ::Type{System})
    damage_storage_type(dmgmodel, ::Type{System}, ::Type{FT})

$(extension_api_note())

Return the type of the state a damage model carries on the system family `System`,
instantiated for the float type `FT` of the simulation. This method is generated by
[`@dmg_storage`](@ref), so a model declared with that macro needs nothing here. The system
type is dispatched on because the fracture bookkeeping is named by the system family, e.g.
[`CriticalStretch`](@ref) carries [`BondFracFields`](@ref) on a bond system and
[`InteractionFracFields`](@ref) on an [`InteractionSystem`](@ref).

# Default

`Nothing`, the answer for a damage model without state.

See also [`@dmg_storage`](@ref), [`damage_state`](@ref), [`storage_type`](@ref).
"""
function damage_storage_type end

function damage_storage_type(dmgmodel, ::Type{<:AbstractSystem},
                             ::Type=default_float_type())
    return Nothing
end

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
[`bond_integrity`](@ref) reach the per-bond variables of a stateful damage model.
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

"""
    has_storage_field(::Type{Storage}, ::Val{field})

$(internal_api_warning())

Return whether a storage type provides a field, flat or through one of its nested states.
The answer is computed when the generated function is expanded, so a guard like
`has_storage_field(typeof(storage), Val(:bond_active))` folds to a constant and the branch
that is not taken disappears from the compiled code. This is what makes the default methods
of [`bond_is_active`](@ref), [`get_damage`](@ref) and [`calc_damage!`](@ref) cost nothing
on a storage that carries the standard bookkeeping and stay correct on one that does not.
"""
@generated function has_storage_field(::Type{S}, ::Val{field}) where {S,field}
    in(field, fieldnames(S)) && return :(true)
    for name in fieldnames(S)
        FT = fieldtype(S, name)
        FT <: Union{AbstractConstitutiveState,AbstractDamageState} || continue
        in(field, fieldnames(FT)) && return :(true)
    end
    return :(false)
end

"""
    bond_is_active(storage, system, bond_id)

$(extension_api_note())

Return whether bond `bond_id` is intact. This is what a force density multiplies into the
influence function or branches on, and on an [`InteractionSystem`](@ref) `bond_id` is the
index of a one-neighbor interaction. Everything outside the damage model asks this function
instead of reading a storage field.

The call is resolved through the state of the damage model, see [`damage_state`](@ref). A
model with its own notion of a broken bond overrides
`bond_is_active(state::MyState, storage, system, bond_id)` on its state type.

# Default

`storage.bond_active` on a bond system and `storage.one_ni_active` on an interaction system
when the storage provides the field, flat or inside a state, and `true` otherwise. A model
that inherits [`BondFracFields`](@ref) therefore needs no method of its own.

# Example

```julia
for bond_id in Peridynamics.each_bond_idx(system, i)
    ω = Peridynamics.bond_is_active(storage, system, bond_id) *
        Peridynamics.surface_correction_factor(system, bond_id)
end
```

See also [`get_damage`](@ref), [`break_bond!`](@ref), [`damage_state`](@ref).
"""
@inline function bond_is_active(storage::AbstractStorage, system::AbstractSystem,
                                bond_id::Int)
    return bond_is_active(damage_state(storage), storage, system, bond_id)
end

"""
    get_damage(storage, i)

$(extension_api_note())

Return the damage of point `i`, the scalar that [`calc_damage!`](@ref) writes once per
point and time step. Everything outside the damage model reads the damage through this
function.

The call is resolved through the state of the damage model, see [`damage_state`](@ref). A
model that stores its damage differently overrides
`get_damage(state::MyState, storage, i)` on its state type.

# Default

`storage.damage` when the storage provides the field, flat or inside a state, and `0.0`
otherwise.

See also [`calc_damage!`](@ref), [`bond_is_active`](@ref), [`damage_state`](@ref).
"""
@inline function get_damage(storage::AbstractStorage, i::Int)
    return get_damage(damage_state(storage), storage, i)
end

"""
    break_bond!(storage, system, dmgmodel, i, bond_id)

$(extension_api_note())

Mark bond `bond_id` of point `i` as failed. This is how code outside the damage model
breaks a single bond, e.g. [`BACMaterial`](@ref) when a bond-associated family cannot carry
stress anymore. The count of active bonds is left alone, the next [`calc_failure!`](@ref)
recomputes it. A damage model with bookkeeping of its own defines a method dispatching on
its model type.

# Default

Writes the standard fracture bookkeeping when the storage provides it, and is a no-op
otherwise.

See also [`break_bonds!`](@ref), [`bond_is_active`](@ref), [`calc_failure!`](@ref).
"""
function break_bond! end

"""
    break_bonds!(storage, system, dmgmodel, i)

$(extension_api_note())

Mark every bond of point `i` as failed and zero its count of active bonds. This is how code
outside the damage model removes a whole point from the force calculation, e.g. the
`maxdmg` mechanism of [`CMaterial`](@ref). A damage model with bookkeeping of its own
defines a method dispatching on its model type.

# Default

Writes the standard fracture bookkeeping when the storage provides it, and is a no-op
otherwise.

See also [`break_bond!`](@ref), [`bond_is_active`](@ref), [`calc_failure!`](@ref).
"""
function break_bonds! end

# a material without a damage model has no bookkeeping to write
@inline break_bond!(storage, system, dmgmodel::Nothing, i, bond_id) = nothing
@inline break_bonds!(storage, system, dmgmodel::Nothing, i) = nothing

"""
    FractureParameters

$(extension_api_note())

Parameter block of the critical energy release rate `Gc` and the critical stretch `εc`,
resolved by [`get_frac_params`](@ref) of the damage model. The block belongs to the damage
model, so it is inherited inside a [`@dmg_params`](@ref) declaration, which is how
[`CriticalStretch`](@ref) declares its parameters. It reads the material `mat` and the
horizon `δ` and the bulk modulus `K` of the material parameters declared above the
`dmg_params::DamageParameters` marker. See [`@params_fields`](@ref).

$(block_table(FractureParameters))
"""
@params_fields FractureParameters begin
    @derived (; Gc, εc) = get_frac_params(model, mat, δ, K; Gc, epsilon_c)
    @log "critical energy release rate" Gc
    @log "critical stretch" εc
end

"""
    CriticalStretchParameters

$(internal_api_warning())

The point parameters of [`CriticalStretch`](@ref): the [`FractureParameters`](@ref) block,
declared with [`@dmg_params`](@ref). They occupy the `dmg_params::DamageParameters` marker
field of every material used with the standard damage model, and are read flat off the
point parameters, e.g. `params.Gc`.

$(block_table(CriticalStretchParameters))
"""
@dmg_params CriticalStretch struct CriticalStretchParameters
    @inherit FractureParameters
end

# --------------------------------------------------------------------------------------
# bond integrity and kinematic weight
#
# Damage acts on a bond in two conceptually different ways, so a damage model that softens
# bonds instead of deleting them answers two different questions. `bond_integrity` is
# constitutive: which fraction of its load does the bond still carry? `kinematic_weight`
# is kinematic: how much can the motion of the neighbor still be trusted when the
# deformation gradient is reconstructed? Both default to the constant one by dispatch,
# which is what every model that deletes bonds wants and costs nothing after inlining.
# --------------------------------------------------------------------------------------

"""
    bond_integrity(dmgmodel, storage, bond_id)

$(extension_api_note())

Return the integrity of bond `bond_id`: the fraction in `[0, 1]` of its undamaged
load-carrying capacity that the bond retains. An intact bond has integrity `1`, a bond that
stores and transmits nothing has integrity `0`. A damage model that softens bonds instead
of deleting them defines a method that reads its own state.

The integrity is honored by the materials that evaluate their constitutive model per bond,
i.e. [`RKCMaterial`](@ref) and [`RKCRMaterial`](@ref), where the stress and the strain
energy density of every bond are scaled by it. Combining a model that defines this method
with a material that ignores it fails once when the [`Job`](@ref) is created, see
[`supports_bond_integrity`](@ref).

# Default

`1.0` for every damage model. A failed bond is excluded through `bond_active` before the
integrity is ever asked, so a model that deletes bonds needs no method of its own.

# Example

```julia
@inline function Peridynamics.bond_integrity(::MyDamage,
                                             storage::Peridynamics.AbstractStorage,
                                             bond_id)
    @inbounds d = Peridynamics.damage_state(storage).bond_damage[bond_id]
    return (1 - d)^2
end
```

See also [`kinematic_weight`](@ref), [`supports_bond_integrity`](@ref),
[`damage_state`](@ref), [`@dmg_storage`](@ref).
"""
function bond_integrity end

# `bond_id` stays untyped so that a model's own method, which the docstring example
# leaves untyped as well, is strictly more specific instead of ambiguous
@inline bond_integrity(::AbstractDamageModel, ::AbstractStorage, bond_id) = 1.0

"""
    kinematic_weight(dmgmodel, storage, bond_id)

$(extension_api_note())

Return the weight in `[0, 1]` with which bond `bond_id` enters the reconstruction of the
deformation gradient, i.e. the moment matrix and the gradient weights of
[`RKCMaterial`](@ref) and [`RKCRMaterial`](@ref).

The kinematic weight is not the [`bond_integrity`](@ref). The integrity states how much
load a bond carries, the kinematic weight whether the motion of its neighbor can still be
trusted as data of the least-squares fit. A softened bond can remain perfectly valid data.
Combining a model that defines this method with a material that ignores it fails once when
the [`Job`](@ref) is created, see [`supports_kinematic_weight`](@ref).

!!! note
    The gradient weights are cached and recomputed only for points whose damage grew. A
    damage model whose kinematic weights evolve continuously has to set
    `storage.update_gradients[i] = true` for the affected points in its own
    [`calc_damage!`](@ref) method.

# Default

`1.0` for every damage model. A failed bond is excluded from the fit through `bond_active`,
so a model that deletes bonds needs no method of its own.

See also [`bond_integrity`](@ref), [`supports_kinematic_weight`](@ref),
[`damage_state`](@ref), [`@dmg_storage`](@ref).
"""
function kinematic_weight end

@inline kinematic_weight(::AbstractDamageModel, ::AbstractStorage, bond_id) = 1.0

"""
    supports_bond_integrity(mat)

$(extension_api_note())

Return whether the force path of a material scales the stress and the strain energy
density of every bond with its [`bond_integrity`](@ref). [`RKCMaterial`](@ref) and
[`RKCRMaterial`](@ref) declare `true`, and a custom material that applies the integrity in
its own force routines declares it the same way.

The declaration is checked once when a [`Job`](@ref) is created, see
[`check_damage_model`](@ref). A damage model that defines [`bond_integrity`](@ref)
combined with a material that ignores it fails there instead of silently not softening.

# Default

`false`.

# Example

```julia
Peridynamics.supports_bond_integrity(::MyMaterial) = true
```

See also [`supports_kinematic_weight`](@ref), [`bond_integrity`](@ref).
"""
function supports_bond_integrity end

supports_bond_integrity(::AbstractMaterial) = false

"""
    supports_kinematic_weight(mat)

$(extension_api_note())

Return whether a material weights the bonds with their [`kinematic_weight`](@ref) when it
reconstructs the deformation gradient. [`RKCMaterial`](@ref) and [`RKCRMaterial`](@ref)
declare `true`, and a custom material that applies the weight in its own gradient
reconstruction declares it the same way.

The declaration is checked once when a [`Job`](@ref) is created, see
[`check_damage_model`](@ref). A damage model that defines [`kinematic_weight`](@ref)
combined with a material that ignores it fails there instead of silently not softening.

# Default

`false`.

# Example

```julia
Peridynamics.supports_kinematic_weight(::MyMaterial) = true
```

See also [`supports_bond_integrity`](@ref), [`kinematic_weight`](@ref).
"""
function supports_kinematic_weight end

supports_kinematic_weight(::AbstractMaterial) = false

# whether the damage model brings its own method for a softening hook: the method that
# dispatch would pick for this model and storage is then not the default above
function overrides_softening_hook(hook, dmgmodel, ::Type{Storage}) where {Storage}
    default = which(hook, Tuple{AbstractDamageModel,AbstractStorage,Int})
    return which(hook, Tuple{typeof(dmgmodel),Storage,Int}) !== default
end

"""
    check_damage_model(spatial_setup)

$(internal_api_warning())

Check that the damage model and the material fit together. This is verified once when a
[`Job`](@ref) is created, next to [`check_constitutive_model`](@ref): a damage model that
defines [`bond_integrity`](@ref) or [`kinematic_weight`](@ref) must be combined with a
material whose force path calls the hooks, see [`supports_bond_integrity`](@ref) and
[`supports_kinematic_weight`](@ref). Throws a [`SofteningSupportError`](@ref) that names
the ignored hooks and how to fix the setup.

There is no structural check of the fracture bookkeeping. Every consumer outside the
damage model reads it through [`bond_is_active`](@ref) and [`get_damage`](@ref), which
behave neutrally when a model carries no bookkeeping, and a pre-crack that cannot be
applied throws its own error.
"""
function check_damage_model(mat::AbstractMaterial)
    dmgmodel = get_dmgmodel(mat)
    isnothing(dmgmodel) && return nothing
    Storage = storage_type(mat)
    ignored = String[]
    if !supports_bond_integrity(mat) &&
       overrides_softening_hook(bond_integrity, dmgmodel, Storage)
        push!(ignored, "bond_integrity")
    end
    if !supports_kinematic_weight(mat) &&
       overrides_softening_hook(kinematic_weight, dmgmodel, Storage)
        push!(ignored, "kinematic_weight")
    end
    isempty(ignored) && return nothing
    M, D = typeof(mat), typeof(dmgmodel)
    reason = "the force path of the material `$(nameof(M))` never calls "
    reason *= join(("`$(hook)`" for hook in ignored), " and ")
    reason *= ", so the softening that the damage model defines would be silently ignored"
    fix = "Use a material that supports the softening hooks, e.g. `RKCMaterial` or "
    fix *= "`RKCRMaterial`. A custom material that calls the hooks in its own force path "
    fix *= "declares that with `Peridynamics.supports_bond_integrity(::MyMaterial) = "
    fix *= "true` and `Peridynamics.supports_kinematic_weight(::MyMaterial) = true`."
    throw(SofteningSupportError(D, M, reason, fix))
end

check_damage_model(body::AbstractBody) = check_damage_model(body.mat)

function check_damage_model(ms::AbstractMultibodySetup)
    for body in each_body(ms)
        check_damage_model(body)
    end
    return nothing
end

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
