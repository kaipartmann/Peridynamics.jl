#=
Everything CriticalStretch does and everything the standard fracture bookkeeping means.

The damage model is a plug-in box: materials, systems and corrections never touch the
bookkeeping by field name, they go through the interface functions declared in
`physics/fracture.jl`. This file owns the blocks that make up the standard bookkeeping,
the states of `CriticalStretch` that carry them, the failure criterion, and the default
methods of the interface functions. The defaults interpret the standard bookkeeping and
are guarded by `has_storage_field`, so they serve every damage model that inherits the
blocks and behave neutrally for one that does not. The struct itself lives in
`physics/fracture.jl`, because the material constructors and the `@dmg_params` declaration
need it before the storage machinery is available.
=#

"""
    BondFracFields

$(extension_api_note())

The fracture bookkeeping of a bond system, see [`@storage_fields`](@ref). The block
belongs to the damage model: a model that deletes bonds inherits it in its
[`@dmg_storage`](@ref) declaration, which is how [`CriticalStretch`](@ref) carries it.
Every bond starts active, the count of active bonds starts at the neighbor count and the
damage starts at zero.

Inheriting this block also opts the model into the default methods of
[`bond_is_active`](@ref), [`get_damage`](@ref), [`break_bond!`](@ref),
[`break_bonds!`](@ref), [`calc_damage!`](@ref) and the pre-crack application, which all
read and write these fields.

$(block_table(BondFracFields))
"""
@storage_fields BondFracFields begin
    damage::PointScalar
    n_active_bonds::PointScalar{Int}
    bond_active::BondScalar{Bool} = true
end

"""
    InteractionFracFields

$(extension_api_note())

The fracture bookkeeping of an interaction system, see [`@storage_fields`](@ref). This is
the [`BondFracFields`](@ref) twin for the one-neighbor interactions: the block belongs to
the damage model and is inherited in its [`@dmg_storage`](@ref) declaration for the
[`InteractionSystem`](@ref), which is how [`CriticalStretch`](@ref) carries it. Every
interaction starts active, the count of active interactions starts at the interaction
count and the damage starts at zero. Like its twin, inheriting the block opts the model
into the default methods of the interface functions.

$(block_table(InteractionFracFields))
"""
@storage_fields InteractionFracFields begin
    damage::PointScalar
    n_active_one_nis::PointScalar{Int}
    one_ni_active::BondScalar{Bool} = true
end

"""
    BondFracState

$(extension_api_note())

The damage state of [`CriticalStretch`](@ref) on a bond system: the
[`BondFracFields`](@ref) block and nothing else. It is declared with
[`@dmg_storage`](@ref), so every material that carries a `dmg_state::DamageState` field
gets the fracture bookkeeping through the damage model, and the fields are read flat off
the storage, e.g. `storage.bond_active`.

$(block_table(BondFracState))
"""
@dmg_storage CriticalStretch AbstractBondSystem struct BondFracState
    @inherit BondFracFields
end

"""
    InteractionFracState

$(extension_api_note())

The damage state of [`CriticalStretch`](@ref) on an [`InteractionSystem`](@ref): the
[`InteractionFracFields`](@ref) block and nothing else, the twin of
[`BondFracState`](@ref) with the field names of the one-neighbor interactions.

$(block_table(InteractionFracState))
"""
@dmg_storage CriticalStretch InteractionSystem struct InteractionFracState
    @inherit InteractionFracFields
end

function calc_failure!(storage::AbstractStorage, system::AbstractBondSystem,
                       mat::AbstractMaterial, dmgmodel::CriticalStretch,
                       paramsetup::AbstractParameterSetup, t, Δt, i)
    (; εc) = get_params(paramsetup, i)
    (; n_active_bonds, bond_active) = storage
    n_active_bonds[i] = 0
    for bond_id in each_bond_idx(system, i)
        ε = bond_stretch(storage, system, i, bond_id)
        if ε > εc && bond_may_fail(system, bond_id)
            bond_active[bond_id] = false
        end
        n_active_bonds[i] += bond_active[bond_id]
    end
    return nothing
end

function calc_failure!(storage::AbstractStorage, system::InteractionSystem,
                       mat::AbstractInteractionSystemMaterial, dmgmodel::CriticalStretch,
                       paramsetup::AbstractParameterSetup, t, Δt, i)
    (; εc) = get_params(paramsetup, i)
    (; n_active_one_nis, one_ni_active) = storage
    n_active_one_nis[i] = 0
    for bond_id in each_one_ni_idx(system, i)
        ε = bond_stretch(storage, system, i, bond_id)
        if ε > εc && bond_may_fail(system, bond_id)
            one_ni_active[bond_id] = false
        end
        n_active_one_nis[i] += one_ni_active[bond_id]
    end
    return nothing
end

# The default damage is the fraction of broken bonds, guarded so that a model without the
# standard bookkeeping leaves the damage alone. A model that reduces its state differently
# defines its own method dispatching on its model type.
function calc_damage!(storage::AbstractStorage, system::AbstractBondSystem,
                      mat::AbstractMaterial, dmgmodel::AbstractDamageModel,
                      paramsetup::AbstractParameterSetup, i)
    S = typeof(storage)
    if has_storage_field(S, Val(:damage)) && has_storage_field(S, Val(:n_active_bonds))
        @inbounds storage.damage[i] = 1 - storage.n_active_bonds[i] / system.n_neighbors[i]
    end
    return nothing
end

function calc_damage!(storage::AbstractStorage, system::InteractionSystem,
                      mat::AbstractInteractionSystemMaterial,
                      dmgmodel::AbstractDamageModel, paramsetup::AbstractParameterSetup, i)
    S = typeof(storage)
    if has_storage_field(S, Val(:damage)) && has_storage_field(S, Val(:n_active_one_nis))
        @inbounds storage.damage[i] = 1 - storage.n_active_one_nis[i] / system.n_neighbors[i]
    end
    return nothing
end

# The default methods of the read hooks, resolved through the damage state. Any state that
# carries the standard bookkeeping gets the behavior of `CriticalStretch`, and the guards
# fold at compile time, see `has_storage_field`. The seam methods are deliberately untyped
# in everything but the state, so a custom model can override on its state type alone
# without creating an ambiguity, and the helpers dispatch on the system for the field name.
@inline function bond_is_active(::Union{Nothing,AbstractDamageState}, storage, system,
                                bond_id)
    return standard_bond_is_active(storage, system, bond_id)
end

# The flag read is `@inbounds`: a bond id always comes from `each_bond_idx` or one of
# its relatives, and the check would otherwise keep the length of the flag vector alive
# inside every force loop, which costs the cheap kernels measurably.
@inline function standard_bond_is_active(storage::AbstractStorage,
                                         system::AbstractBondSystem, bond_id::Int)
    has_storage_field(typeof(storage), Val(:bond_active)) || return true
    return @inbounds storage.bond_active[bond_id]
end

@inline function standard_bond_is_active(storage::AbstractStorage,
                                         system::InteractionSystem, one_ni_id::Int)
    has_storage_field(typeof(storage), Val(:one_ni_active)) || return true
    return @inbounds storage.one_ni_active[one_ni_id]
end

@inline function get_damage(::Union{Nothing,AbstractDamageState}, storage, i)
    has_storage_field(typeof(storage), Val(:damage)) || return 0.0
    return storage.damage[i]
end

# The default methods of the write hooks, dispatching on the damage model because the
# caller always has it at hand. Same guards, same neutrality, and again untyped in
# everything but the model so that an override on the model type alone dominates.
function break_bond!(storage, system, dmgmodel::AbstractDamageModel, i, bond_id)
    standard_break_bond!(storage, system, bond_id)
    return nothing
end

function standard_break_bond!(storage::AbstractStorage, system::AbstractBondSystem,
                              bond_id::Int)
    has_storage_field(typeof(storage), Val(:bond_active)) || return nothing
    storage.bond_active[bond_id] = false
    return nothing
end

function standard_break_bond!(storage::AbstractStorage, system::InteractionSystem,
                              one_ni_id::Int)
    has_storage_field(typeof(storage), Val(:one_ni_active)) || return nothing
    storage.one_ni_active[one_ni_id] = false
    return nothing
end

function break_bonds!(storage, system, dmgmodel::AbstractDamageModel, i)
    standard_break_bonds!(storage, system, i)
    return nothing
end

function standard_break_bonds!(storage::AbstractStorage, system::AbstractBondSystem, i)
    S = typeof(storage)
    if has_storage_field(S, Val(:bond_active))
        storage.bond_active[each_bond_idx(system, i)] .= false
    end
    if has_storage_field(S, Val(:n_active_bonds))
        storage.n_active_bonds[i] = 0
    end
    return nothing
end

function standard_break_bonds!(storage::AbstractStorage, system::InteractionSystem, i)
    S = typeof(storage)
    if has_storage_field(S, Val(:one_ni_active))
        storage.one_ni_active[each_one_ni_idx(system, i)] .= false
    end
    if has_storage_field(S, Val(:n_active_one_nis))
        storage.n_active_one_nis[i] = 0
    end
    return nothing
end

"""
    failure_by_sets!(storage, system, dmgmodel, set_a, set_b)

$(internal_api_warning())

Break every bond between a point of `set_a` and a point of `set_b`, which is how a
predefined crack is applied to a body chunk, and keep the fracture bookkeeping of the
storage in sync. The default writes the standard bookkeeping and throws an error for a
damage model that does not carry it, because silently ignoring a pre-crack would falsify
the simulation. A damage model with its own state defines a method that also writes the
crack into that state:

```julia
function Peridynamics.failure_by_sets!(storage, system::Peridynamics.AbstractBondSystem,
                                       dmg::MyDamage, set_a, set_b)
    Peridynamics.failure_by_sets!(storage, system, CriticalStretch(), set_a, set_b)
    # ... mark the broken bonds in `Peridynamics.damage_state(storage)`
    return nothing
end
```
"""
function failure_by_sets!(storage, system::AbstractBondSystem,
                          dmgmodel::AbstractDamageModel, set_a, set_b)
    check_precrack_bookkeeping(storage, dmgmodel, Val(:bond_active), Val(:n_active_bonds))
    (; n_active_bonds, bond_active) = storage
    n_active_bonds .= 0
    for i in each_point_idx(system)
        for bond_id in each_bond_idx(system, i)
            neighbor_id = get_neighbor(system, bond_id)
            point_in_a = in(i, set_a)
            point_in_b = in(i, set_b)
            neigh_in_a = in(neighbor_id, set_a)
            neigh_in_b = in(neighbor_id, set_b)
            if (point_in_a && neigh_in_b) || (point_in_b && neigh_in_a)
                bond_active[bond_id] = false
            end
            n_active_bonds[i] += bond_active[bond_id]
        end
    end
    return nothing
end

function failure_by_sets!(storage, system::InteractionSystem,
                          dmgmodel::AbstractDamageModel, set_a, set_b)
    check_precrack_bookkeeping(storage, dmgmodel, Val(:one_ni_active),
                               Val(:n_active_one_nis))
    storage.n_active_one_nis .= 0
    for point_id in each_point_idx(system)
        for bond_id in each_one_ni_idx(system, point_id)
            neighbor_id = get_neighbor(system, bond_id)
            point_in_a = in(point_id, set_a)
            point_in_b = in(point_id, set_b)
            neigh_in_a = in(neighbor_id, set_a)
            neigh_in_b = in(neighbor_id, set_b)
            if (point_in_a && neigh_in_b) || (point_in_b && neigh_in_a)
                storage.one_ni_active[bond_id] = false
            end
            storage.n_active_one_nis[point_id] += storage.one_ni_active[bond_id]
        end
    end
    return nothing
end

function check_precrack_bookkeeping(storage, dmgmodel, active_field::Val, count_field::Val)
    if has_storage_field(typeof(storage), active_field) &&
       has_storage_field(typeof(storage), count_field)
        return nothing
    end
    D = nameof(typeof(dmgmodel))
    msg = "cannot apply a pre-crack with the damage model `$(D)`!\n"
    msg *= "The model carries no standard fracture bookkeeping, so there is nothing to\n"
    msg *= "write the broken bonds into. Either create the pre-crack with\n"
    msg *= "`precrack!(body, :set_a, :set_b; update_dmg=false)` so that the bonds are\n"
    msg *= "filtered out during setup, inherit the bookkeeping block of the system family\n"
    msg *= "in the `@dmg_storage` declaration of the model, e.g.\n"
    msg *= "    Peridynamics.@dmg_storage $(D) struct MyState\n"
    msg *= "        @inherit BondFracFields\n"
    msg *= "    end\n"
    msg *= "or define `Peridynamics.failure_by_sets!` for the model."
    return error(msg)
end
