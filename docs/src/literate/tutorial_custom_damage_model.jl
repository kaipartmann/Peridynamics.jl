# # [Writing your own damage model](@id tutorial_custom_damage_model)

# A damage model decides when a bond fails. The built-in [`CriticalStretch`](@ref) breaks a
# bond the moment its stretch exceeds the critical value ``\varepsilon_c``. Many materials do
# not fail that abruptly: a bond that is overstretched for a short moment survives, one that
# stays overstretched breaks. That is a failure criterion with a delay, in the spirit of the
# time-integrated criteria used for spallation, and it needs two things the built-in model
# does not have: a state per bond, and the time step.
#
# We let a bond accumulate damage ``D`` while it is stretched beyond ``\varepsilon_c``,
#
# ```math
# \dot{D} = \frac{1}{\tau} \left( \frac{\varepsilon}{\varepsilon_c} - 1 \right)_+ ,
# ```
#
# and break it once ``D \geq 1``. The delay ``\tau`` is a new material parameter. With
# ``\tau \to 0`` the model turns into `CriticalStretch`. Everything here is written against
# the [Extension API](@ref), and the model works with every material of the package.

using Peridynamics
using Peridynamics: BondSystem, each_bond_idx, get_params, get_n_loc_points,
                    get_vector_diff, damage_state
using Peridynamics.LinearAlgebra: norm

# ## The type
#
# A damage model is a subtype of
# [`AbstractDamageModel`](@ref Peridynamics.AbstractDamageModel). It has no fields of its own
# here, because its parameters belong to the point parameters, where they can differ between
# point sets and appear in the simulation log.

struct DelayedFailure <: Peridynamics.AbstractDamageModel end

# ## The parameters
#
# A damage model owns its parameters and declares them with
# [`@dmg_params`](@ref Peridynamics.@dmg_params), in the same language as the point
# parameters of a material. Inheriting [`FractureParameters`](@ref Peridynamics.FractureParameters)
# brings the standard pair `Gc` and `εc`, including the conversion between them that the
# material provides, so `material!(...; Gc)` and `material!(...; epsilon_c)` both work. The
# delay is a keyword of its own, `tau`, stored as `τ`.

Peridynamics.@dmg_params DelayedFailure struct DelayedFailureParameters
    @inherit FractureParameters
    @log "failure delay" @kwarg tau τ
end

# Every material whose point parameters carry the marker `dmg_params::DamageParameters`
# accepts these keywords now, and all materials of the package do. The parameters are read
# flat off the point parameters, e.g. `params.τ`, next to the parameters of the material.
# Whether fracture is enabled follows from `Gc` and `εc` by default, so there is nothing to
# define for that either.

# ## The state
#
# The accumulated damage of every bond has to survive from one time step to the next.
# [`@dmg_storage`](@ref Peridynamics.@dmg_storage) declares it with the field shapes of
# [`@storage`](@ref Peridynamics.@storage). The state is allocated with the storage of
# whatever material the model is attached to, moves with it to another array backend and is
# never exchanged between chunks, because a bond belongs to one point. A model without a
# state simply skips this declaration.

Peridynamics.@dmg_storage DelayedFailure struct DelayedFailureState
    bond_damage::BondScalar
end

# ## The criterion
#
# This is the one method a damage model has to define. It is called once per local point and
# per time step, right before the force density, and it receives the time and the time step.
# Its contract is short: deactivate the bonds that fail, count the ones that are still active
# in `n_active_bonds`, and never break a bond whose `fail_permit` is `false`, because that is
# how [`no_failure!`](@ref) and the pre-cracks are honored.
#
# The state is reached with [`damage_state`](@ref Peridynamics.damage_state). The bonds are
# read exactly as in a force density: `each_bond_idx`, `system.bonds[bond_id]`, and the
# current positions from the storage.

function Peridynamics.calc_failure!(storage, system::BondSystem, mat, ::DelayedFailure,
                                    paramsetup, t, Δt, i)
    (; εc, τ) = get_params(paramsetup, i)
    (; bond_damage) = damage_state(storage)
    for bond_id in each_bond_idx(system, i)
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        Δxij = get_vector_diff(storage.position, i, j)
        ε = (norm(Δxij) - L) / L
        if storage.bond_active[bond_id] && bond.fail_permit && ε > εc
            bond_damage[bond_id] += (ε / εc - 1) * Δt / τ
            if bond_damage[bond_id] ≥ 1
                storage.bond_active[bond_id] = false
            end
        end
        storage.n_active_bonds[i] += storage.bond_active[bond_id]
    end
    return nothing
end

# The damage of a point, the fraction of its bonds that failed, is computed by the package
# right after this, so a model that deletes bonds does not define
# [`calc_damage!`](@ref Peridynamics.calc_damage!). A model that softens bonds instead of
# deleting them does, and it also defines
# [`bond_integrity`](@ref Peridynamics.bond_integrity) and
# [`kinematic_weight`](@ref Peridynamics.kinematic_weight), see the [Extension API](@ref).

# ## Exporting the accumulated damage
#
# `bond_damage` lives per bond. To look at it, it is reduced to one value per point, the
# maximum over the bonds of a point, and announced with
# [`custom_field`](@ref Peridynamics.custom_field), so that a `Job` accepts the name.

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
# how far every point is from losing its most loaded bond. The same file runs with
# `julia -t 6` and under `mpiexec` unchanged.

# ## Where to go next
#
# - [Writing your own material](@ref tutorial_custom_material) for the material half.
# - [Writing your own constitutive model](@ref tutorial_custom_constitutive_model) for a
#   stress-strain relation with a history. A damage model can read the state of such a model
#   through [`constitutive_state`](@ref Peridynamics.constitutive_state), which is how a
#   ductile damage model reads the plastic strain.
# - [Materials](@ref) for the declaration language in full.
# - [Extension API](@ref) for every name used here.
