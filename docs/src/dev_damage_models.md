# Damage models

A damage model decides when a bond fails. The tutorial
[Writing your own damage model](@ref tutorial_custom_damage_model) writes one with a state
and a parameter of its own in full, this page is the lookup table.

## Contract

| what you define | signature | required | default or fallback |
|:---|:---|:---:|:---|
| the type | `struct MyDamage <: Peridynamics.AbstractDamageModel end` | yes | |
| the failure criterion | [`calc_failure!(storage, system, mat, ::MyDamage, paramsetup, t, Δt, i)`](@ref Peridynamics.calc_failure!) | yes | [`InterfaceError`](@ref Peridynamics.InterfaceError) |
| the parameters of the model | `Peridynamics.@dmg_params MyDamage struct MyDamageParameters ... end` | no | the model has no parameters |
| the state of the model | `Peridynamics.@dmg_storage MyDamage struct MyDamageState ... end`, or `Peridynamics.@dmg_storage MyDamage System struct MyDamageState ... end` for one system family | no | the model has no state, every bond is active and the damage is zero |
| the damage of a point | [`calc_damage!(storage, system, mat, dmgmodel, paramsetup, i)`](@ref Peridynamics.calc_damage!) | no | the fraction of broken bonds |
| another notion of an intact bond | [`bond_is_active(state, storage, system, bond_id)`](@ref Peridynamics.bond_is_active), on the **state** type | no | reads `bond_active`, or `one_ni_active` on an interaction system |
| another notion of the damage | [`get_damage(state, storage, i)`](@ref Peridynamics.get_damage), on the **state** type | no | reads `damage` |
| how an outsider breaks a bond | [`break_bond!(storage, system, dmgmodel, i, bond_id)`](@ref Peridynamics.break_bond!) and [`break_bonds!(storage, system, dmgmodel, i)`](@ref Peridynamics.break_bonds!), on the **model** type | no | writes the standard bookkeeping, does nothing without it |
| other fracture keywords | [`get_frac_params(dmgmodel, mat, δ, K; Gc=nothing, epsilon_c=nothing, kwargs...)`](@ref Peridynamics.get_frac_params) | no | converts `Gc` and `epsilon_c` into each other |
| whether fracture is enabled | [`has_fracture(dmgmodel, params)`](@ref Peridynamics.has_fracture) | no | reads `Gc` and `εc` and requires both to be nonzero |
| softening instead of deleting | [`bond_integrity(dmgmodel, storage, bond_id)`](@ref Peridynamics.bond_integrity), [`kinematic_weight(dmgmodel, storage, bond_id)`](@ref Peridynamics.kinematic_weight) | no | `1.0`, see the softening section |

Inheriting [`FractureParameters`](@ref Peridynamics.FractureParameters) inside the
`@dmg_params` body brings the standard pair `Gc` and `εc`, and inheriting
[`BondFracFields`](@ref Peridynamics.BondFracFields) or
[`InteractionFracFields`](@ref Peridynamics.InteractionFracFields) inside the `@dmg_storage`
body brings `bond_active`, `n_active_bonds` and `damage`, and with them every default above.

## Skeleton

```julia
struct MyDamage <: Peridynamics.AbstractDamageModel end

Peridynamics.@dmg_params MyDamage struct MyDamageParameters
    @inherit FractureParameters
    @log "failure delay" @kwarg tau τ
end

Peridynamics.@dmg_storage MyDamage struct MyDamageState
    @inherit BondFracFields
    bond_damage::BondScalar
end

function Peridynamics.calc_failure!(storage, system, mat, ::MyDamage, paramsetup, t, Δt, i)
    (; εc) = get_params(paramsetup, i)
    storage.n_active_bonds[i] = 0
    for bond_id in each_bond_idx(system, i)
        ε = bond_stretch(storage, system, i, bond_id)
        ...
        storage.n_active_bonds[i] += storage.bond_active[bond_id]
    end
    return nothing
end
```

## The plug-in box

A damage model is a box with walls. These are the walls, and each of them has a default as
soon as the state inherits `BondFracFields`.

| function | read or write | who calls it | default with `BondFracFields` |
|:---|:---|:---|:---|
| [`calc_failure!`](@ref Peridynamics.calc_failure!) | write | the force density loop, once per local point and time step | none, this is the one method a model defines |
| [`calc_damage!`](@ref Peridynamics.calc_damage!) | write | the force density loop, right after `calc_failure!` | the fraction of broken bonds |
| [`bond_is_active`](@ref Peridynamics.bond_is_active) | read | materials, corrections, every other consumer | the `bond_active` flag of the bond |
| [`get_damage`](@ref Peridynamics.get_damage) | read | export, the `maxdmg` mechanism, the gradient update decision | the `damage` of the point |
| [`break_bond!`](@ref Peridynamics.break_bond!) | write | a material that kills one bond, e.g. `BACMaterial` | deactivate the bond, the count is left to the next `calc_failure!` |
| [`break_bonds!`](@ref Peridynamics.break_bonds!) | write | a material that removes a whole point, e.g. `CMaterial` | deactivate every bond of the point and zero its count |

## Rules

- **Inside** the model read and write the state flat, e.g. `storage.bond_active`, or with
  [`damage_state`](@ref Peridynamics.damage_state). **Outside** the model everything goes
  through the interface functions above.
- A bond for which [`bond_may_fail(system, bond_id)`](@ref Peridynamics.bond_may_fail) is
  `false` must never break. That is how [`no_failure!`](@ref) and the pre-cracks are
  honored.
- Reset `storage.n_active_bonds[i]` at the top of `calc_failure!`, then add every bond that
  is still active, because that count is what the default `calc_damage!` turns into the
  damage of the point.
- Read the stretch with [`bond_stretch`](@ref Peridynamics.bond_stretch) and the length with
  [`current_bond_length`](@ref Peridynamics.current_bond_length), never by gathering the two
  positions and taking the norm.
- A damage state does **not** make anything history dependent. The state is advanced in
  `calc_failure!`, which every time solver calls exactly once per step, so a stateful damage
  model stays usable under [`NewtonKrylov`](@ref).
- A model without bookkeeping runs with every bond active and zero damage. Only a pre-crack
  it cannot apply is an error.
- A model with a state needs a material whose storage declares `dmg_state::DamageState`,
  which every storage of this package does. This is checked when a [`Job`](@ref) is created
  and throws a [`StorageContractError`](@ref Peridynamics.StorageContractError).

## Fracture parameters

| you want ... | define |
|:---|:---|
| the standard keywords `Gc` and `epsilon_c` | `@inherit FractureParameters` inside `@dmg_params`, and nothing else |
| another micro-modulus, so another relation between `Gc` and `εc` | [`critical_stretch`](@ref Peridynamics.critical_stretch) **and** [`energy_release_rate`](@ref Peridynamics.energy_release_rate) on the material, always both |
| other fracture keywords, e.g. a critical stress | [`get_frac_params`](@ref Peridynamics.get_frac_params) |
| to decide yourself whether the bonds of a point set may fail | [`has_fracture`](@ref Peridynamics.has_fracture) |

Giving both `Gc` and `epsilon_c` is an error, and giving neither switches fracture off.

## Softening a bond instead of deleting it

Deleting a bond is a jump in the moment matrix of a reproducing kernel material, and no
regularization of its inverse can absorb a jump in its input. A model can therefore let a
bond fade out instead, through two hooks that both default to `1.0`.

| hook | scales | meaning |
|:---|:---|:---|
| [`bond_integrity`](@ref Peridynamics.bond_integrity) | the stress and the strain energy the bond carries | the continuity `1 - d` of classical damage mechanics, evaluated per bond |
| [`kinematic_weight`](@ref Peridynamics.kinematic_weight) | what the bond contributes to the moment matrix and the gradient weights | takes a neighbor across a forming crack out of the least-squares fit smoothly, so the moment matrix stays a continuous function of the damage |

A material says with
[`supports_bond_integrity`](@ref Peridynamics.supports_bond_integrity) and
[`supports_kinematic_weight`](@ref Peridynamics.supports_kinematic_weight) whether its force
path calls the hooks. [`RKCMaterial`](@ref) and [`RKCRMaterial`](@ref) do. Combining a
softening model with a material that ignores the hooks throws a
[`SofteningSupportError`](@ref Peridynamics.SofteningSupportError) when the [`Job`](@ref) is
created, instead of silently not softening.

A model that softens also defines `calc_damage!`, because the default damage of a point is
the fraction of deleted bonds, which is not what a softening model means.

See also [Materials](@ref), [Storages](@ref), [Point parameters](@ref),
[Damage in peridynamics formulations](@ref expl_dmg), [Blocks you can inherit](@ref),
[Extension API](@ref).
