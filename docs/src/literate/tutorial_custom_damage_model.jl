# # [Writing your own damage model](@id tutorial_custom_damage_model)

# A damage model decides when a bond fails. The built-in [`CriticalStretch`](@ref) breaks a
# bond once its stretch exceeds ``\varepsilon_c``. Here a bond instead accumulates damage
# ``D``,
#
# ```math
# \dot{D} = \frac{1}{\tau} \left( \frac{\varepsilon}{\varepsilon_c} - 1 \right)_+ ,
# ```
#
# and breaks once ``D \geq 1``, with the delay ``\tau`` a new material parameter, written
# against the [Extension API](@ref).

using Peridynamics
using Peridynamics: BondSystem, each_bond_idx, get_params, get_n_loc_points, damage_state,
                    bond_stretch, bond_may_fail

# ## The type
#
# A damage model is a subtype of
# [`AbstractDamageModel`](@ref Peridynamics.AbstractDamageModel), with no fields of its own,
# its parameters living in the point parameters instead.

struct DelayedFailure <: Peridynamics.AbstractDamageModel end

# ## The parameters
#
# A damage model declares its parameters with [`@dmg_params`](@ref Peridynamics.@dmg_params).
# Inheriting [`FractureParameters`](@ref Peridynamics.FractureParameters) brings `Gc` and
# `εc`, and the delay is its own keyword, `tau`, stored as `τ`.

Peridynamics.@dmg_params DelayedFailure struct DelayedFailureParameters
    @inherit FractureParameters
    @log "failure delay" @kwarg tau τ
end

# Every material whose point parameters carry `dmg_params::DamageParameters` accepts these
# keywords now, read flat as `params.τ`, and fracture is enabled by default from `Gc`/`εc`.

# ## The state
#
# A damage model owns the fracture bookkeeping. Inheriting
# [`BondFracFields`](@ref Peridynamics.BondFracFields) brings `bond_active`, `n_active_bonds`
# and `damage`, plus the defaults of [`bond_is_active`](@ref Peridynamics.bond_is_active) and
# the rest of the interface. The accumulated damage per bond is a field of our own, declared
# with [`@storage`](@ref Peridynamics.@storage), never exchanged between chunks.

Peridynamics.@dmg_storage DelayedFailure struct DelayedFailureState
    @inherit BondFracFields
    bond_damage::BondScalar
end

# ## The criterion
#
# This is the one method a damage model has to define, run once per local point per time
# step, right before the force density: reset the point's bond count, deactivate failing
# bonds, count the active ones in `n_active_bonds`, and never break a bond for which
# [`bond_may_fail`](@ref Peridynamics.bond_may_fail) is `false`, honoring
# [`no_failure!`](@ref) and pre-cracks. The state is reached with
# [`damage_state`](@ref Peridynamics.damage_state), the stretch with
# [`bond_stretch`](@ref Peridynamics.bond_stretch).

function Peridynamics.calc_failure!(storage, system::BondSystem, mat, ::DelayedFailure,
                                    paramsetup, t, Δt, i)
    (; εc, τ) = get_params(paramsetup, i)
    (; bond_damage) = damage_state(storage)
    storage.n_active_bonds[i] = 0
    for bond_id in each_bond_idx(system, i)
        ε = bond_stretch(storage, system, i, bond_id)
        if storage.bond_active[bond_id] && bond_may_fail(system, bond_id) && ε > εc
            bond_damage[bond_id] += (ε / εc - 1) * Δt / τ
            if bond_damage[bond_id] ≥ 1
                storage.bond_active[bond_id] = false
            end
        end
        storage.n_active_bonds[i] += storage.bond_active[bond_id]
    end
    return nothing
end

# The package computes the damage of a point right after this, so a model that deletes bonds
# does not define [`calc_damage!`](@ref Peridynamics.calc_damage!). A model that softens
# bonds instead also defines [`bond_integrity`](@ref Peridynamics.bond_integrity) and
# [`kinematic_weight`](@ref Peridynamics.kinematic_weight).

# ## Exporting the accumulated damage
#
# `bond_damage` lives per bond, reduced to one value per point, the maximum over its bonds,
# and announced with [`custom_field`](@ref Peridynamics.custom_field) so a `Job` accepts it.

Peridynamics.custom_field(::Type{<:Peridynamics.AbstractStorage}, ::Val{:bond_damage}) = true

function Peridynamics.export_field(::Val{:bond_damage}, mat, system, storage, paramsetup, t)
    (; bond_damage) = damage_state(storage)
    out = zeros(get_n_loc_points(system))
    for i in eachindex(out)
        bond_ids = each_bond_idx(system, i)
        isempty(bond_ids) || (out[i] = maximum(@view bond_damage[bond_ids]))
    end
    return out
end

# ## Running it
#
# The model is finished. It is plugged into a material like the built-in one, here into the
# bond-based material of the package, and `material!` takes `tau` next to `Gc`:

l, Δx = 0.1, 0.002
pos, vol = uniform_box(l, 0.1l, 0.1l, Δx)
body = Body(BBMaterial(; dmgmodel=DelayedFailure()), pos, vol)
material!(body; horizon=3.015Δx, rho=2700, E=70e9, Gc=100, tau=2e-6)

point_set!(x -> x < -0.4l, body, :left)
point_set!(x -> x > 0.4l, body, :right)
velocity_bc!(t -> -10.0, body, :left, :x)
velocity_bc!(t -> 10.0, body, :right, :x)

job = Job(body, VelocityVerlet(steps=400);
          path=joinpath(tempdir(), "delayed_failure"),
          fields=(:displacement, :damage, :bond_damage))

#-
#md # ```julia
#md # submit(job)
#md # ```

# No bond breaks in the step it is first overstretched now, and the exported field records
# how far a point is from losing its most loaded bond.

# ## Where to go next
#
# - [Writing your own material](@ref tutorial_custom_material) for the material half.
# - [Writing your own constitutive model](@ref tutorial_custom_constitutive_model) for a
#   stress-strain relation with a history. A damage model can read the state of such a model
#   through [`constitutive_state`](@ref Peridynamics.constitutive_state), which is how a
#   ductile damage model reads the plastic strain.
# - [Materials](@ref) for the declaration language in full.
# - [Extension API](@ref) for every name used here.
